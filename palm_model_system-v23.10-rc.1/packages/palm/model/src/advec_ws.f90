!> @file advec_ws.f90
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
! Authors:
! ! --------
! @author Matthias Suehring
!
!
! Description:
! ------------
!> Advection scheme for scalars and momentum using the flux formulation of Wicker and Skamarock 5th
!> order. Additionally the module contains of a routine using for initialisation and steering of the
!> statical evaluation.
!> The computation of turbulent fluxes takes place inside the advection routines.
!> Near non-cyclic boundaries the order of the applied advection scheme is degraded.
!> A divergence correction is applied. It is necessary for topography, since the divergence is not
!> sufficiently reduced, resulting in erroneous fluxes and could lead to numerical instabilities.
!>
!> @todo Implement monotonic flux limiter also for vector version.
!> @todo Move 3d arrays advc_flag, advc_flags_m from modules to advec_ws
!> @todo Move arrays flux_l_x from modules to advec_ws
!--------------------------------------------------------------------------------------------------!
 MODULE advec_ws

    USE arrays_3d,                                                                                 &
        ONLY:  ddzu,                                                                               &
               ddzw,                                                                               &
               diss_l_diss,                                                                        &
               diss_l_e,                                                                           &
               diss_l_pt,                                                                          &
               diss_l_q,                                                                           &
               diss_l_s,                                                                           &
               diss_l_u,                                                                           &
               diss_l_v,                                                                           &
               diss_l_w,                                                                           &
               diss_s_diss,                                                                        &
               diss_s_e,                                                                           &
               diss_s_pt,                                                                          &
               diss_s_q,                                                                           &
               diss_s_s,                                                                           &
               diss_s_u,                                                                           &
               diss_s_v,                                                                           &
               diss_s_w,                                                                           &
               drho_air,                                                                           &
               drho_air_zw,                                                                        &
               rho_air,                                                                            &
               rho_air_zw,                                                                         &
               flux_l_diss,                                                                        &
               flux_l_e,                                                                           &
               flux_l_pt,                                                                          &
               flux_l_q,                                                                           &
               flux_l_s,                                                                           &
               flux_l_u,                                                                           &
               flux_l_v,                                                                           &
               flux_l_w,                                                                           &
               flux_s_diss,                                                                        &
               flux_s_e,                                                                           &
               flux_s_pt,                                                                          &
               flux_s_q,                                                                           &
               flux_s_s,                                                                           &
               flux_s_u,                                                                           &
               flux_s_v,                                                                           &
               flux_s_w,                                                                           &
               tend,                                                                               &
               u,                                                                                  &
               u_stokes_zu,                                                                        &
               v,                                                                                  &
               v_stokes_zu,                                                                        &
               w

    USE control_parameters,                                                                        &
        ONLY:  bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               dt_3d,                                                                              &
               humidity,                                                                           &
               intermediate_timestep_count,                                                        &
               loop_optimization,                                                                  &
               passive_scalar,                                                                     &
               rans_tke_e,                                                                         &
               symmetry_flag,                                                                      &
               u_gtrans,                                                                           &
               v_gtrans,                                                                           &
               ws_scheme_mom,                                                                      &
               ws_scheme_sca

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz,                                                                     &
               exchange_horiz_int

    USE indices,                                                                                   &
        ONLY:  advc_flags_m,                                                                       &
               advc_flags_s,                                                                       &
               nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxlu,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nysv,                                                                               &
               nzb,                                                                                &
               nzb_max,                                                                            &
               nzt,                                                                                &
               topo_flags

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy

    USE kinds

    USE pegrid,                                                                                    &
        ONLY:  comm2d,                                                                             &
               threads_per_task

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               sums_salsa_ws_l,                                                                    &
               sums_us2_ws_l,                                                                      &
               sums_vs2_ws_l,                                                                      &
               sums_ws2_ws_l,                                                                      &
               sums_wschs_ws_l,                                                                    &
               sums_wsncs_ws_l,                                                                    &
               sums_wsnrs_ws_l,                                                                    &
               sums_wspts_ws_l,                                                                    &
               sums_wsqs_ws_l,                                                                     &
               sums_wsss_ws_l,                                                                     &
               sums_wsqcs_ws_l,                                                                    &
               sums_wsqrs_ws_l,                                                                    &
               sums_wsqis_ws_l,                                                                    &
               sums_wsnis_ws_l,                                                                    &
               sums_wsqgs_ws_l,                                                                    &
               sums_wsngs_ws_l,                                                                    &
               sums_wsqss_ws_l,                                                                    &
               sums_wsnss_ws_l,                                                                    &
               sums_wssas_ws_l,                                                                    &
               sums_wsus_ws_l,                                                                     &
               sums_wsvs_ws_l,                                                                     &
               weight_substep



    IMPLICIT NONE

    REAL(wp) ::  adv_mom_1            !< 1/4 - constant used in 5th-order advection scheme for momentum advection (1st-order part)
    REAL(wp) ::  adv_mom_3            !< 1/24 - constant used in 5th-order advection scheme for momentum advection (3rd-order part)
    REAL(wp) ::  adv_mom_5            !< 1/120 - constant used in 5th-order advection scheme for momentum advection (5th-order part)
    REAL(wp) ::  adv_sca_1            !< 1/2 - constant used in 5th-order advection scheme for scalar advection (1st-order part)
    REAL(wp) ::  adv_sca_3            !< 1/12 - constant used in 5th-order advection scheme for scalar advection (3rd-order part)
    REAL(wp) ::  adv_sca_5            !< 1/60 - constant used in 5th-order advection scheme for scalar advection (5th-order part)

    PRIVATE
    PUBLIC   advec_s_ws,                                                                           &
             advec_u_ws,                                                                           &
             advec_v_ws,                                                                           &
             advec_w_ws,                                                                           &
             ws_init,                                                                              &
             ws_init_flags_scalar,                                                                 &
             ws_statistics

    INTERFACE ws_init
       MODULE PROCEDURE ws_init
    END INTERFACE ws_init

    INTERFACE ws_init_flags_scalar
       MODULE PROCEDURE ws_init_flags_scalar
    END INTERFACE ws_init_flags_scalar

    INTERFACE ws_statistics
       MODULE PROCEDURE ws_statistics
    END INTERFACE ws_statistics

    INTERFACE advec_s_ws
       MODULE PROCEDURE advec_s_ws
       MODULE PROCEDURE advec_s_ws_div_bc
       MODULE PROCEDURE advec_s_ws_ij
    END INTERFACE advec_s_ws

    INTERFACE advec_u_ws
       MODULE PROCEDURE advec_u_ws
       MODULE PROCEDURE advec_u_ws_ij
    END INTERFACE advec_u_ws

    INTERFACE advec_v_ws
       MODULE PROCEDURE advec_v_ws
       MODULE PROCEDURE advec_v_ws_ij
    END INTERFACE advec_v_ws

    INTERFACE advec_w_ws
       MODULE PROCEDURE advec_w_ws
       MODULE PROCEDURE advec_w_ws_ij
    END INTERFACE advec_w_ws

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of WS-scheme
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ws_init

!
!-- Set factors for scalar and momentum advection.
    adv_sca_5 = 1.0_wp /  60.0_wp
    adv_sca_3 = 1.0_wp /  12.0_wp
    adv_sca_1 = 1.0_wp /   2.0_wp
    adv_mom_5 = 1.0_wp / 120.0_wp
    adv_mom_3 = 1.0_wp /  24.0_wp
    adv_mom_1 = 1.0_wp /   4.0_wp
!
!-- Arrays needed for statical evaluation of fluxes.
    IF ( ws_scheme_mom )  THEN

       ALLOCATE( sums_wsus_ws_l(nzb:nzt+1,0:threads_per_task-1),                                   &
                 sums_wsvs_ws_l(nzb:nzt+1,0:threads_per_task-1),                                   &
                 sums_us2_ws_l(nzb:nzt+1,0:threads_per_task-1),                                    &
                 sums_vs2_ws_l(nzb:nzt+1,0:threads_per_task-1),                                    &
                 sums_ws2_ws_l(nzb:nzt+1,0:threads_per_task-1) )

       sums_wsus_ws_l = 0.0_wp
       sums_wsvs_ws_l = 0.0_wp
       sums_us2_ws_l  = 0.0_wp
       sums_vs2_ws_l  = 0.0_wp
       sums_ws2_ws_l  = 0.0_wp

    ENDIF

    IF ( ws_scheme_sca )  THEN

       ALLOCATE( sums_wspts_ws_l(nzb:nzt+1,0:threads_per_task-1) )
       sums_wspts_ws_l = 0.0_wp

       IF ( humidity  )  THEN
          ALLOCATE( sums_wsqs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
          sums_wsqs_ws_l = 0.0_wp
       ENDIF

       IF ( passive_scalar )  THEN
          ALLOCATE( sums_wsss_ws_l(nzb:nzt+1,0:threads_per_task-1) )
          sums_wsss_ws_l = 0.0_wp
       ENDIF

    ENDIF

!
!-- Arrays needed for reasons of speed optimization
    IF ( ws_scheme_mom )  THEN

       ALLOCATE( flux_s_u(nzb+1:nzt,0:threads_per_task-1),                                         &
                 flux_s_v(nzb+1:nzt,0:threads_per_task-1),                                         &
                 flux_s_w(nzb+1:nzt,0:threads_per_task-1),                                         &
                 diss_s_u(nzb+1:nzt,0:threads_per_task-1),                                         &
                 diss_s_v(nzb+1:nzt,0:threads_per_task-1),                                         &
                 diss_s_w(nzb+1:nzt,0:threads_per_task-1) )
       ALLOCATE( flux_l_u(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                                 &
                 flux_l_v(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                                 &
                 flux_l_w(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                                 &
                 diss_l_u(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                                 &
                 diss_l_v(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                                 &
                 diss_l_w(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )

    ENDIF
!
!-- For the vector version the buffer arrays for scalars are not necessary, since internal arrays
!-- are used in the vector version.
    IF ( loop_optimization /= 'vector' )  THEN
       IF ( ws_scheme_sca )  THEN

          ALLOCATE( flux_s_pt(nzb+1:nzt,0:threads_per_task-1),                                     &
                    flux_s_e(nzb+1:nzt,0:threads_per_task-1),                                      &
                    diss_s_pt(nzb+1:nzt,0:threads_per_task-1),                                     &
                    diss_s_e(nzb+1:nzt,0:threads_per_task-1) )
          ALLOCATE( flux_l_pt(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                             &
                    flux_l_e(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                              &
                    diss_l_pt(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                             &
                    diss_l_e(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )

          IF ( rans_tke_e )  THEN
             ALLOCATE( flux_s_diss(nzb+1:nzt,0:threads_per_task-1),                                &
                       diss_s_diss(nzb+1:nzt,0:threads_per_task-1) )
             ALLOCATE( flux_l_diss(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                        &
                       diss_l_diss(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
          ENDIF

          IF ( humidity )  THEN
             ALLOCATE( flux_s_q(nzb+1:nzt,0:threads_per_task-1),                                   &
                       diss_s_q(nzb+1:nzt,0:threads_per_task-1) )
             ALLOCATE( flux_l_q(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                           &
                       diss_l_q(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
          ENDIF

          IF ( passive_scalar )  THEN
             ALLOCATE( flux_s_s(nzb+1:nzt,0:threads_per_task-1),                                   &
                       diss_s_s(nzb+1:nzt,0:threads_per_task-1) )
             ALLOCATE( flux_l_s(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                           &
                       diss_l_s(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
          ENDIF

       ENDIF
    ENDIF
!
!-- Initialize the flag arrays controlling degradation near walls, i.e. to decrease the numerical
!-- stencil appropriately. The order of the scheme is degraded near solid walls as well as near
!-- non-cyclic inflow and outflow boundaries. Do this separately for momentum and scalars.
    IF ( ws_scheme_mom )  CALL ws_init_flags_momentum

    IF ( ws_scheme_sca )  THEN
       ALLOCATE( advc_flags_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       advc_flags_s = 0
       CALL ws_init_flags_scalar( bc_dirichlet_l  .OR.  bc_radiation_l,                            &
                                  bc_dirichlet_n  .OR.  bc_radiation_n,                            &
                                  bc_dirichlet_r  .OR.  bc_radiation_r,                            &
                                  bc_dirichlet_s  .OR.  bc_radiation_s,                            &
                                  advc_flags_s )
    ENDIF

 END SUBROUTINE ws_init

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of flags to control the order of the advection scheme near solid walls and
!> non-cyclic inflow boundaries, where the order is sucessively degraded.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ws_init_flags_momentum


    INTEGER(iwp) ::  i     !< index variable along x
    INTEGER(iwp) ::  j     !< index variable along y
    INTEGER(iwp) ::  k     !< index variable along z
    INTEGER(iwp) ::  k_mm  !< dummy index along z
    INTEGER(iwp) ::  k_pp  !< dummy index along z
    INTEGER(iwp) ::  k_ppp !< dummy index along z

    LOGICAL      ::  flag_set !< steering variable for advection flags

    ALLOCATE( advc_flags_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    advc_flags_m = 0
!
!-- Set advc_flags_m to steer the degradation of the advection scheme in advec_ws near
!-- topography, inflow- and outflow boundaries as well as bottom and top of model domain.
!-- advc_flags_m remains zero for all non-prognostic grid points.
!-- u-component
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          At first, set flags to WS1.
!--          Since fluxes are swapped in advec_ws.f90, this is necessary to
!--          in order to handle the left/south flux.
!--          near vertical walls.
             advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 0 )
             advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 3 )
!
!--          u component - x-direction
!--          WS1 (0), WS3 (1), WS5 (2)
             IF ( .NOT. BTEST(topo_flags(k,j,i+1),1)  .OR.                                         &
                      ( ( bc_dirichlet_l .OR. bc_radiation_l ) .AND. i <= nxlu  )  .OR.            &
                      ( ( bc_dirichlet_r .OR. bc_radiation_r ) .AND. i == nxr   ) )                &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 0 )
             ELSEIF ( ( .NOT. BTEST(topo_flags(k,j,i+2),1)  .AND.                                  &
                              BTEST(topo_flags(k,j,i+1),1)  .OR.                                   &
                        .NOT. BTEST(topo_flags(k,j,i-1),1) )                                       &
                                              .OR.                                                 &
                      ( ( bc_dirichlet_r .OR. bc_radiation_r ) .AND. i == nxr-1 )  .OR.            &
                      ( ( bc_dirichlet_l .OR. bc_radiation_l ) .AND. i == nxlu+1) )                &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 1 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 0 )
             ELSEIF ( BTEST(topo_flags(k,j,i+1),1)  .AND.                                          &
                      BTEST(topo_flags(k,j,i+2),1)  .AND.                                          &
                      BTEST(topo_flags(k,j,i-1),1) )                                               &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 2 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 0 )
             ENDIF
!
!--          u component - y-direction
!--          WS1 (3), WS3 (4), WS5 (5)
             IF ( .NOT. BTEST(topo_flags(k,j+1,i),1)  .OR.                                         &
                      ( ( bc_dirichlet_s .OR. bc_radiation_s ) .AND. j == nys   )  .OR.            &
                      ( ( bc_dirichlet_n .OR. bc_radiation_n ) .AND. j == nyn   ) )                &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 3 )
             ELSEIF ( ( .NOT. BTEST(topo_flags(k,j+2,i),1)  .AND.                                  &
                              BTEST(topo_flags(k,j+1,i),1)  .OR.                                   &
                        .NOT. BTEST(topo_flags(k,j-1,i),1) )                                       &
                                              .OR.                                                 &
                      ( ( bc_dirichlet_s .OR. bc_radiation_s ) .AND. j == nysv  )  .OR.            &
                      ( ( bc_dirichlet_n .OR. bc_radiation_n ) .AND. j == nyn-1 ) )                &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 4 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 3 )
             ELSEIF ( BTEST(topo_flags(k,j+1,i),1)  .AND.                                          &
                      BTEST(topo_flags(k,j+2,i),1)  .AND.                                          &
                      BTEST(topo_flags(k,j-1,i),1) )                                               &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 5 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 3 )
             ENDIF
!
!--          u component - z-direction. Fluxes are calculated on w-grid level. Boundary u-values
!--          at/within walls aren't used.
!--          WS1 (6), WS3 (7), WS5 (8)
             IF ( k == nzb+1 )  THEN
                k_mm = nzb
             ELSE
                k_mm = k - 2
             ENDIF
             IF ( k > nzt-1 )  THEN
                k_pp = nzt+1
             ELSE
                k_pp = k + 2
             ENDIF
             IF ( k > nzt-2 )  THEN
                k_ppp = nzt+1
             ELSE
                k_ppp = k + 3
             ENDIF

             flag_set = .FALSE.
             IF ( ( .NOT. BTEST(topo_flags(k-1,j,i),1)                 .AND.                       &
                          BTEST(topo_flags(k,j,i),1)                   .AND.                       &
                          BTEST(topo_flags(k+1,j,i),1) )               .OR.                        &
                  ( .NOT. BTEST(topo_flags(k_pp,j,i),1)                .AND.                       &
                          BTEST(topo_flags(k+1,j,i),1)                 .AND.                       &
                          BTEST(topo_flags(k,j,i),1) )                 .OR.                        &
                  ( k == nzt .AND. symmetry_flag == 0 )                .OR.                        &
                  ( k == nzb+1  .AND.  BTEST(topo_flags(k,j,i),1)                                  &
                                .AND.  BTEST(topo_flags(k+1,j,i),1) ) )                            &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 6 )
                flag_set = .TRUE.
             ELSEIF ( ( .NOT. BTEST(topo_flags(k_mm,j,i),1)            .OR.                        &
                        .NOT. BTEST(topo_flags(k_ppp,j,i),1) )         .AND.                       &
                              BTEST(topo_flags(k-1,j,i),1)             .AND.                       &
                              BTEST(topo_flags(k,j,i),1)               .AND.                       &
                              BTEST(topo_flags(k+1,j,i),1)             .AND.                       &
                              BTEST(topo_flags(k_pp,j,i),1)            .AND.                       &
                        .NOT. flag_set                                 .OR.                        &
                      ( k == nzt-1  .AND.  symmetry_flag == 0 )        .OR.                        &
                      ( k == nzb+2  .AND.  BTEST(topo_flags(k,j,i),1)                              &
                                    .AND.  BTEST(topo_flags(k+1,j,i),1)                            &
                                    .AND.  BTEST(topo_flags(k_pp,j,i),1)                           &
                                    .AND.  BTEST(topo_flags(k-1,j,i),1) ) )                        &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 7 )
                flag_set = .TRUE.
             ELSEIF ( BTEST(topo_flags(k_mm,j,i),1)                    .AND.                       &
                      BTEST(topo_flags(k-1,j,i),1)                     .AND.                       &
                      BTEST(topo_flags(k,j,i),1)                       .AND.                       &
                      BTEST(topo_flags(k+1,j,i),1)                     .AND.                       &
                      BTEST(topo_flags(k_pp,j,i),1)                    .AND.                       &
                      BTEST(topo_flags(k_ppp,j,i),1)                   .AND.                       &
                      .NOT. flag_set )                                                             &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 8 )
             ENDIF

          ENDDO
!
!--       Special treatment for k = nzb in case of elevated nesting. If lower boundary is "open",
!--       set 1st order flag.
          k = nzb
          IF ( BTEST( topo_flags(k,j,i), 1 )  .AND.  BTEST( topo_flags(k+1,j,i), 1 ) )  THEN
             advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 6 )
          ENDIF

       ENDDO
    ENDDO
!
!-- v-component
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          At first, set flags to WS1.
!--          Since fluxes are swapped in advec_ws.f90, this is necessary to in order to handle the
!--          left/south flux.
             advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 9  )
             advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 12 )
!
!--          v component - x-direction
!--          WS1 (9), WS3 (10), WS5 (11)
             IF ( .NOT. BTEST(topo_flags(k,j,i+1),2)  .OR.                                         &
                      ( ( bc_dirichlet_l .OR. bc_radiation_l ) .AND. i == nxl   )  .OR.            &
                      ( ( bc_dirichlet_r .OR. bc_radiation_r ) .AND. i == nxr   ) )                &
            THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 9 )
!
!--          WS3
             ELSEIF ( ( .NOT. BTEST(topo_flags(k,j,i+2),2)    .AND.                                &
                              BTEST(topo_flags(k,j,i+1),2) )  .OR.                                 &
                        .NOT. BTEST(topo_flags(k,j,i-1),2)                                         &
                                              .OR.                                                 &
                      ( ( bc_dirichlet_r .OR. bc_radiation_r ) .AND. i == nxr-1 )  .OR.            &
                      ( ( bc_dirichlet_l .OR. bc_radiation_l ) .AND. i == nxlu  ) )                &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 10 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 9 )
             ELSEIF ( BTEST(topo_flags(k,j,i+1),2)  .AND.                                          &
                      BTEST(topo_flags(k,j,i+2),2)  .AND.                                          &
                      BTEST(topo_flags(k,j,i-1),2) )                                               &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 11 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 9 )
             ENDIF
!
!--          v component - y-direction
!--          WS1 (12), WS3 (13), WS5 (14)
             IF ( .NOT. BTEST(topo_flags(k,j+1,i),2)  .OR.                                         &
                      ( ( bc_dirichlet_s .OR. bc_radiation_s ) .AND. j <= nysv  )   .OR.           &
                      ( ( bc_dirichlet_n .OR. bc_radiation_n ) .AND. j == nyn   ) )                &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 12 )
             ELSEIF ( ( .NOT. BTEST(topo_flags(k,j+2,i),2)  .AND.                                  &
                              BTEST(topo_flags(k,j+1,i),2)  .OR.                                   &
                        .NOT. BTEST(topo_flags(k,j-1,i),2) )                                       &
                                              .OR.                                                 &
                      ( (  bc_dirichlet_s .OR. bc_radiation_s ) .AND. j == nysv+1)  .OR.           &
                      ( (  bc_dirichlet_n .OR. bc_radiation_n ) .AND. j == nyn-1 ) )               &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 13 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 12 )
             ELSEIF ( BTEST(topo_flags(k,j+1,i),2)  .AND.                                          &
                      BTEST(topo_flags(k,j+2,i),2)  .AND.                                          &
                      BTEST(topo_flags(k,j-1,i),2) )                                               &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 14 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 12 )
             ENDIF
!
!--          v component - z-direction. Fluxes are calculated on w-grid level. Boundary v-values
!--          at/within walls aren't used.
!--          WS1 (15), WS3 (16), WS5 (17)
             IF ( k == nzb+1 )  THEN
                k_mm = nzb
             ELSE
                k_mm = k - 2
             ENDIF
             IF ( k > nzt-1 )  THEN
                k_pp = nzt+1
             ELSE
                k_pp = k + 2
             ENDIF
             IF ( k > nzt-2 )  THEN
                k_ppp = nzt+1
             ELSE
                k_ppp = k + 3
             ENDIF

             flag_set = .FALSE.
             IF ( ( .NOT. BTEST(topo_flags(k-1,j,i),2)                 .AND.                       &
                          BTEST(topo_flags(k,j,i),2)                   .AND.                       &
                          BTEST(topo_flags(k+1,j,i),2) )               .OR.                        &
                  ( .NOT. BTEST(topo_flags(k_pp,j,i),2)                .AND.                       &
                          BTEST(topo_flags(k+1,j,i),2)                 .AND.                       &
                          BTEST(topo_flags(k,j,i),2) )                 .OR.                        &
                  ( k == nzt .AND. symmetry_flag == 0 )                .OR.                        &
                  ( k == nzb+1  .AND.  BTEST(topo_flags(k,j,i),2)                                  &
                                .AND.  BTEST(topo_flags(k+1,j,i),2) ) )                            &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 15 )
                flag_set = .TRUE.
             ELSEIF ( ( .NOT. BTEST(topo_flags(k_mm,j,i),2)            .OR.                        &
                        .NOT. BTEST(topo_flags(k_ppp,j,i),2) )         .AND.                       &
                              BTEST(topo_flags(k-1,j,i),2)             .AND.                       &
                              BTEST(topo_flags(k,j,i),2)               .AND.                       &
                              BTEST(topo_flags(k+1,j,i),2)             .AND.                       &
                              BTEST(topo_flags(k_pp,j,i),2)            .AND.                       &
                        .NOT. flag_set                                 .OR.                        &
                      ( k == nzt-1  .AND.  symmetry_flag == 0 )        .OR.                        &
                      ( k == nzb+2  .AND.  BTEST(topo_flags(k,j,i),2)                              &
                                    .AND.  BTEST(topo_flags(k+1,j,i),2)                            &
                                    .AND.  BTEST(topo_flags(k_pp,j,i),2)                           &
                                    .AND.  BTEST(topo_flags(k-1,j,i),2) ) )                        &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 16 )
                flag_set = .TRUE.
             ELSEIF ( BTEST(topo_flags(k_mm,j,i),2)                    .AND.                       &
                      BTEST(topo_flags(k-1,j,i),2)                     .AND.                       &
                      BTEST(topo_flags(k,j,i),2)                       .AND.                       &
                      BTEST(topo_flags(k+1,j,i),2)                     .AND.                       &
                      BTEST(topo_flags(k_pp,j,i),2)                    .AND.                       &
                      BTEST(topo_flags(k_ppp,j,i),2)                   .AND.                       &
                      .NOT. flag_set )                                                             &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 17 )
             ENDIF

          ENDDO
!
!--       Special treatment for k = nzb in case of elevated nesting. If lower boundary is "open",
!--       set 1st order flag.
          k = nzb
          IF ( BTEST( topo_flags(k,j,i), 2 )  .AND.  BTEST( topo_flags(k+1,j,i), 2 ) )  THEN
             advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 15 )
          ENDIF

       ENDDO
    ENDDO
!
!-- w - component
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          At first, set flags to WS1.
!--          Since fluxes are swapped in advec_ws.f90, this is necessary to in order to handle the
!--          left/south flux.
             advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 18 )
             advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 21 )
!
!--          w component - x-direction
!--          WS1 (18), WS3 (19), WS5 (20)
             IF ( .NOT. BTEST(topo_flags(k,j,i+1),3)  .OR.                                         &
                      ( (  bc_dirichlet_l .OR. bc_radiation_l ) .AND. i == nxl   )  .OR.           &
                      ( (  bc_dirichlet_r .OR. bc_radiation_r ) .AND. i == nxr   ) )               &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 18 )
             ELSEIF ( ( .NOT. BTEST(topo_flags(k,j,i+2),3)  .AND.                                  &
                              BTEST(topo_flags(k,j,i+1),3)  .OR.                                   &
                        .NOT. BTEST(topo_flags(k,j,i-1),3) )                                       &
                                              .OR.                                                 &
                      ( ( bc_dirichlet_r .OR. bc_radiation_r )  .AND. i == nxr-1 )  .OR.           &
                      ( ( bc_dirichlet_l .OR.  bc_radiation_l ) .AND. i == nxlu  ) )               &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 19 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 18 )
             ELSEIF ( BTEST(topo_flags(k,j,i+1),3)  .AND.                                          &
                      BTEST(topo_flags(k,j,i+2),3)  .AND.                                          &
                      BTEST(topo_flags(k,j,i-1),3) )                                               &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i),20 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 18 )
             ENDIF
!
!--          w component - y-direction
!--          WS1 (21), WS3 (22), WS5 (23)
             IF ( .NOT. BTEST(topo_flags(k,j+1,i),3)  .OR.                                         &
                      ( ( bc_dirichlet_s .OR. bc_radiation_s ) .AND. j == nys   )  .OR.            &
                      ( ( bc_dirichlet_n .OR. bc_radiation_n ) .AND. j == nyn   ) )                &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 21 )
             ELSEIF ( ( .NOT. BTEST(topo_flags(k,j+2,i),3)  .AND.                                  &
                              BTEST(topo_flags(k,j+1,i),3)  .OR.                                   &
                        .NOT. BTEST(topo_flags(k,j-1,i),3) )                                       &
                                              .OR.                                                 &
                      ( ( bc_dirichlet_s .OR. bc_radiation_s ) .AND. j == nysv  )  .OR.            &
                      ( ( bc_dirichlet_n .OR. bc_radiation_n ) .AND. j == nyn-1 ) )                &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 22 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 21 )
             ELSEIF ( BTEST(topo_flags(k,j+1,i),3)  .AND.                                          &
                      BTEST(topo_flags(k,j+2,i),3)  .AND.                                          &
                      BTEST(topo_flags(k,j-1,i),3) )                                               &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 23 )
!
!--             Clear flag for WS1
                advc_flags_m(k,j,i) = IBCLR( advc_flags_m(k,j,i), 21 )
             ENDIF
!
!--          w component - z-direction. Fluxes are calculated on scalar grid level. Boundary
!--          w-values at walls are used. Flux at k=i is defined at scalar position k=i+1 with i
!--          being an integer.
!--          WS1 (24), WS3 (25), WS5 (26)
             IF ( k == nzb+1 )  THEN
                k_mm = nzb
             ELSE
                k_mm = k - 2
             ENDIF
             IF ( k > nzt-1 )  THEN
                k_pp = nzt+1
             ELSE
                k_pp = k + 2
             ENDIF
             IF ( k > nzt-2 )  THEN
                k_ppp = nzt+1
             ELSE
                k_ppp = k + 3
             ENDIF

             flag_set = .FALSE.
             IF ( ( .NOT. BTEST(topo_flags(k,j,i),3)                  .AND.                        &
                          BTEST(topo_flags(k+1,j,i),3) )              .OR.                         &
                  ( .NOT. BTEST(topo_flags(k+1,j,i),3)                .AND.                        &
                          BTEST(topo_flags(k,j,i),3) )                .OR.                         &
                  ( k == nzt-1 ) )                                                                 &
             THEN
!
!--             At k = nzb (or the uppermost topography grid point) a flag is explicitly set,
!--             although this is not a prognostic level. However, contrary to the advection of
!--             u, v and s, this is necessary because flux_t(nzb) is used for the
!--             tendency at k = nzb+1.
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 24 )
                flag_set = .TRUE.
             ELSEIF ( ( .NOT. BTEST(topo_flags(k-1,j,i),3)            .AND.                        &
                              BTEST(topo_flags(k,j,i),3)              .AND.                        &
                              BTEST(topo_flags(k+1,j,i),3)            .AND.                        &
                              BTEST(topo_flags(k_pp,j,i),3) )         .OR.                         &
                      ( .NOT. BTEST(topo_flags(k_pp,j,i),3)           .AND.                        &
                              BTEST(topo_flags(k+1,j,i),3)            .AND.                        &
                              BTEST(topo_flags(k,j,i),3)              .AND.                        &
                              BTEST(topo_flags(k-1,j,i),3) )          .AND.                        &
                        .NOT. flag_set                                .OR.                         &
                      ( k == nzt-2 )                                  .OR.                         &
                      ( k == nzb+1  .AND.  BTEST(topo_flags(k,j,i),3)                              &
                                    .AND.  BTEST(topo_flags(k+1,j,i),3) ) )                        &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 25 )
                flag_set = .TRUE.
             ELSEIF ( BTEST(topo_flags(k-1,j,i),3)                    .AND.                        &
                      BTEST(topo_flags(k,j,i),3)                      .AND.                        &
                      BTEST(topo_flags(k+1,j,i),3)                    .AND.                        &
                      BTEST(topo_flags(k_pp,j,i),3)                   .AND.                        &
                      .NOT. flag_set )                                                             &
             THEN
                advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 26 )
             ENDIF

          ENDDO
!
!--       Special treatment for k = nzb in case of elevated nesting. If lower boundary is "open",
!--       set 1st order flag.
          k = nzb
          IF ( BTEST( topo_flags(k,j,i), 3 )  .AND.  BTEST(topo_flags(k+1,j,i), 3 ) )  THEN
             advc_flags_m(k,j,i) = IBSET( advc_flags_m(k,j,i), 24 )
          ENDIF

       ENDDO
    ENDDO
!
!-- Exchange ghost points for advection flags
    CALL exchange_horiz_int( advc_flags_m, nys, nyn, nxl, nxr, nzt, nbgp )
!
!-- Set boundary flags at inflow and outflow boundary in case of
!-- non-cyclic boundary conditions.
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       advc_flags_m(:,:,nxl-1) = advc_flags_m(:,:,nxl)
    ENDIF

    IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
       advc_flags_m(:,:,nxr+1) = advc_flags_m(:,:,nxr)
    ENDIF

    IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
       advc_flags_m(:,nyn+1,:) = advc_flags_m(:,nyn,:)
    ENDIF

    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       advc_flags_m(:,nys-1,:) = advc_flags_m(:,nys,:)
    ENDIF

 END SUBROUTINE ws_init_flags_momentum


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of flags to control the order of the advection scheme near solid walls and
!> non-cyclic inflow boundaries, where the order is sucessively degraded.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ws_init_flags_scalar( non_cyclic_l, non_cyclic_n, non_cyclic_r, non_cyclic_s,          &
                                  advc_flag, extensive_degrad, alternative_communicator )

    INTEGER(iwp), OPTIONAL ::  alternative_communicator  !< alternative MPI communicator to be used

    INTEGER(iwp) ::  i            !< index variable along x
    INTEGER(iwp) ::  j            !< index variable along y
    INTEGER(iwp) ::  k            !< index variable along z
    INTEGER(iwp) ::  k_mm         !< dummy index along z
    INTEGER(iwp) ::  k_pp         !< dummy index along z
    INTEGER(iwp) ::  k_ppp        !< dummy index along z

    INTEGER(iwp), INTENT(INOUT), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::                       &
                                           advc_flag !< flag array to control order of scalar advection

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  topo_flags_alt !< temporary alternative array for topography flags

    LOGICAL ::  flag_set     !< steering variable for advection flags
    LOGICAL ::  non_cyclic_l !< flag that indicates non-cyclic boundary on the left
    LOGICAL ::  non_cyclic_n !< flag that indicates non-cyclic boundary on the north
    LOGICAL ::  non_cyclic_r !< flag that indicates non-cyclic boundary on the right
    LOGICAL ::  non_cyclic_s !< flag that indicates non-cyclic boundary on the south

    LOGICAL, OPTIONAL ::  extensive_degrad !< flag indicating that extensive degradation is required, e.g. for
                                           !< passive scalars nearby topography along the horizontal directions,
                                           !< as no monotonic limiter can be applied there
!
!-- In case an alternative communicator is used, the topography flag on the ghost points needs
!-- to be set correspondingly. Therefore, set a temporary alternative topography array and exchange
!-- ghost points using the alternative communicator.
    ALLOCATE( topo_flags_alt(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    topo_flags_alt = topo_flags

    IF ( PRESENT( alternative_communicator ) )  THEN
       CALL exchange_horiz_int( topo_flags_alt, nys, nyn, nxl, nxr, nzt, nbgp,                     &
                                alternative_communicator = alternative_communicator )
!
!--    In case of non-cyclic boundaries set a Neumann condition for the topography.
       IF ( non_cyclic_l )  topo_flags_alt(:,:,nxl-1) = topo_flags_alt(:,:,nxl)
       IF ( non_cyclic_r )  topo_flags_alt(:,:,nxr+1) = topo_flags_alt(:,:,nxr)
       IF ( non_cyclic_n )  topo_flags_alt(:,nyn+1,:) = topo_flags_alt(:,nyn,:)
       IF ( non_cyclic_s )  topo_flags_alt(:,nys-1,:) = topo_flags_alt(:,nys,:)
    ENDIF
!
!-- Set flags to steer the degradation of the advection scheme in advec_ws near topography, inflow-
!-- and outflow boundaries as well as bottom and top of model domain. advc_flags_m remains zero for
!-- all non-prognostic grid points.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             IF ( .NOT.  BTEST(topo_flags_alt(k,j,i),0) )  CYCLE
!
!--          scalar - x-direction
!--          WS1 (0), WS3 (1), WS5 (2)
             IF ( ( .NOT. BTEST(topo_flags_alt(k,j,i+1),0)              .OR.                       &
                    .NOT. BTEST(topo_flags_alt(k,j,i+2),0)              .OR.                       &
                    .NOT. BTEST(topo_flags_alt(k,j,i-1),0) )            .OR.                       &
                    ( non_cyclic_l  .AND.  i == 0  )                .OR.                           &
                    ( non_cyclic_r  .AND.  i == nx ) )                                             &
             THEN
                advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 0 )
             ELSEIF ( ( .NOT. BTEST(topo_flags_alt(k,j,i+3),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j,i+1),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j,i+2),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j,i-1),0)                                     &
                      )                       .OR.                                                 &
                      ( .NOT. BTEST(topo_flags_alt(k,j,i-2),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j,i+1),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j,i+2),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j,i-1),0)                                     &
                      )                       .OR.                                                 &
                      ( non_cyclic_r  .AND.  i == nx-1 )            .OR.                           &
                      ( non_cyclic_l  .AND.  i == 1    ) )                                         &
             THEN
                advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 1 )
             ELSEIF ( BTEST(topo_flags_alt(k,j,i+1),0)                  .AND.                      &
                      BTEST(topo_flags_alt(k,j,i+2),0)                  .AND.                      &
                      BTEST(topo_flags_alt(k,j,i+3),0)                  .AND.                      &
                      BTEST(topo_flags_alt(k,j,i-1),0)                  .AND.                      &
                      BTEST(topo_flags_alt(k,j,i-2),0) )                                           &
             THEN
                advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 2 )
             ENDIF
!
!--          scalar - y-direction
!--          WS1 (3), WS3 (4), WS5 (5)
             IF ( ( .NOT. BTEST(topo_flags_alt(k,j+1,i),0)              .OR.                       &
                    .NOT. BTEST(topo_flags_alt(k,j+2,i),0)              .OR.                       &
                    .NOT. BTEST(topo_flags_alt(k,j-1,i),0))             .OR.                       &
                  ( non_cyclic_s  .AND.  j == 0  )                  .OR.                           &
                  ( non_cyclic_n  .AND.  j == ny ) )                                               &
             THEN
                advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 3 )
!
!--          WS3
             ELSEIF ( ( .NOT. BTEST(topo_flags_alt(k,j+3,i),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j+1,i),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j+2,i),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j-1,i),0)                                     &
                      )                       .OR.                                                 &
                      ( .NOT. BTEST(topo_flags_alt(k,j-2,i),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j+1,i),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j+2,i),0)          .AND.                      &
                              BTEST(topo_flags_alt(k,j-1,i),0)                                     &
                      )                       .OR.                                                 &
                      ( non_cyclic_s  .AND.  j == 1    )  .OR.                                     &
                      ( non_cyclic_n  .AND.  j == ny-1 ) )                                         &
             THEN
                advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 4 )
!
!--          WS5
             ELSEIF ( BTEST(topo_flags_alt(k,j+1,i),0)                  .AND.                      &
                      BTEST(topo_flags_alt(k,j+2,i),0)                  .AND.                      &
                      BTEST(topo_flags_alt(k,j+3,i),0)                  .AND.                      &
                      BTEST(topo_flags_alt(k,j-1,i),0)                  .AND.                      &
                      BTEST(topo_flags_alt(k,j-2,i),0) )                                           &
             THEN
                advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 5 )
             ENDIF
!
!--          Near topography, set horizontal advection scheme to 1st order for passive scalars, even
!--          if only one direction may be blocked by topography. These locations will be identified
!--          by topo_flags bit 31. Note, since several modules define advection flags but
!--          may apply different scalar boundary conditions, bit 31 is temporarily stored on
!--          advc_flags.
!--          Moreover, note that this extended degradtion for passive scalars is not required for
!--          the vertical direction as there the monotonic limiter can be applied.
             IF ( PRESENT( extensive_degrad ) )  THEN
                IF ( extensive_degrad )  THEN
!
!--                At all grid points that are within a three-grid point range to topography, set
!--                1st-order scheme.
                   IF( BTEST( advc_flag(k,j,i), 31 ) )  THEN
!
!--                   Clear flags that might indicate higher-order advection along x- and
!--                   y-direction.
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 1 )
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 2 )
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 4 )
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 5 )
!
!--                   Set flags that indicate 1st-order advection along x- and y-direction.
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 0 )
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 3 )
                   ENDIF
!
!--                Adjacent to this extended degradation zone, successively upgrade the order of the
!--                scheme if this grid point isn't flagged with bit 31 (indicating extended
!--                degradation zone).
                   IF ( .NOT. BTEST( advc_flag(k,j,i), 31 ) )  THEN
!
!--                   x-direction. First, clear all previous settings, then set flag for 3rd-order
!--                   scheme.
                      IF ( BTEST( advc_flag(k,j,i-1), 31 )  .AND.                                  &
                           BTEST( advc_flag(k,j,i+1), 31 ) )  THEN
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 0 )
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 1 )
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 2 )

                         advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 1 )
                      ENDIF
!
!--                   x-direction. First, clear all previous settings, then set flag for 5rd-order
!--                   scheme.
                      IF ( .NOT. BTEST( advc_flag(k,j,i-1), 31 )  .AND.                            &
                                 BTEST( advc_flag(k,j,i-2), 31 )  .AND.                            &
                           .NOT. BTEST( advc_flag(k,j,i+1), 31 )  .AND.                            &
                                 BTEST( advc_flag(k,j,i+2), 31 ) )  THEN
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 0 )
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 1 )
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 2 )

                         advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 2 )
                      ENDIF
!
!--                   y-direction. First, clear all previous settings, then set flag for 3rd-order
!--                   scheme.
                      IF ( BTEST( advc_flag(k,j-1,i), 31 )  .AND.                                  &
                           BTEST( advc_flag(k,j+1,i), 31 ) )  THEN
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 3 )
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 4 )
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 5 )

                         advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 4 )
                      ENDIF
!
!--                   y-direction. First, clear all previous settings, then set flag for 5rd-order
!--                   scheme.
                      IF ( .NOT. BTEST( advc_flag(k,j-1,i), 31 )  .AND.                            &
                                 BTEST( advc_flag(k,j-2,i), 31 )  .AND.                            &
                           .NOT. BTEST( advc_flag(k,j+1,i), 31 )  .AND.                            &
                                 BTEST( advc_flag(k,j+2,i), 31 ) )  THEN
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 3 )
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 4 )
                         advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 5 )

                         advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 5 )
                      ENDIF
                   ENDIF

                ENDIF

!
!--             Near lateral boundaries set the flags again. In order to avoid strong numerical
!--             oscillations near the boundaries, which may lead to scalar built-up, also employ
!--             extended degradation zones here.
!--             x-direction
                IF ( ( non_cyclic_l  .AND.  i <= 3  )  .OR.                                        &
                     ( non_cyclic_r  .AND.  i >= nx-3 ) )  THEN
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 0 )
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 1 )
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 2 )

                   advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 0 )
                ENDIF

                IF ( ( non_cyclic_l  .AND.  i == 4    )  .OR.                                      &
                     ( non_cyclic_r  .AND.  i == nx-4 ) )  THEN
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 0 )
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 1 )
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 2 )

                   advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 1 )
                ENDIF
!
!--             y-direction
                IF ( ( non_cyclic_n  .AND.  j <= 3  )  .OR.                                        &
                     ( non_cyclic_s  .AND.  j >= ny-3 ) )  THEN
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 3 )
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 4 )
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 5 )

                   advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 3 )
                ENDIF

                IF ( ( non_cyclic_n  .AND.  j == 4    )  .OR.                                      &
                     ( non_cyclic_s  .AND.  j == ny-4 ) )  THEN
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 3 )
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 4 )
                   advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 5 )

                   advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 4 )
                ENDIF

             ENDIF
!
!--          scalar - z-direction. Fluxes are calculated on w-grid level. Boundary values at/within
!--          walls aren't used.
!--          WS1 (6), WS3 (7), WS5 (8)
             IF ( k == nzb+1 )  THEN
                k_mm = nzb
             ELSE
                k_mm = k - 2
             ENDIF
             IF ( k > nzt-1 )  THEN
                k_pp = nzt+1
             ELSE
                k_pp = k + 2
             ENDIF
             IF ( k > nzt-2 )  THEN
                k_ppp = nzt+1
             ELSE
                k_ppp = k + 3
             ENDIF

             flag_set = .FALSE.
             IF ( ( .NOT. BTEST(topo_flags_alt(k-1,j,i),0)               .AND.                     &
                          BTEST(topo_flags_alt(k,j,i),0)                 .AND.                     &
                          BTEST(topo_flags_alt(k+1,j,i),0) )             .OR.                      &
                  ( .NOT. BTEST(topo_flags_alt(k_pp,j,i),0)              .AND.                     &
                          BTEST(topo_flags_alt(k+1,j,i),0)               .AND.                     &
                          BTEST(topo_flags_alt(k,j,i),0) )               .OR.                      &
                  ( k == nzt  .AND.  symmetry_flag == 0 )                .OR.                      &
                  ( k == nzb+1  .AND.  BTEST(topo_flags_alt(k,j,i),0)                              &
                                .AND.  BTEST(topo_flags_alt(k+1,j,i),0) ) )                        &
             THEN
                advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 6 )
                flag_set = .TRUE.
             ELSEIF ( ( .NOT. BTEST(topo_flags_alt(k_mm,j,i),0)          .OR.                      &
                        .NOT. BTEST(topo_flags_alt(k_ppp,j,i),0) )       .AND.                     &
                              BTEST(topo_flags_alt(k-1,j,i),0)           .AND.                     &
                              BTEST(topo_flags_alt(k,j,i),0)             .AND.                     &
                              BTEST(topo_flags_alt(k+1,j,i),0)           .AND.                     &
                              BTEST(topo_flags_alt(k_pp,j,i),0)          .AND.                     &
                        .NOT. flag_set                                   .OR.                      &
                      ( k == nzt-1  .AND.  symmetry_flag == 0 )          .OR.                      &
                      ( k == nzb+2  .AND.  BTEST(topo_flags_alt(k,j,i),0)                          &
                                    .AND.  BTEST(topo_flags_alt(k+1,j,i),0)                        &
                                    .AND.  BTEST(topo_flags_alt(k_pp,j,i),0)                       &
                                    .AND.  BTEST(topo_flags_alt(k-1,j,i),0) ) )                    &
             THEN
                advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 7 )
                flag_set = .TRUE.
             ELSEIF ( BTEST(topo_flags_alt(k_mm,j,i),0)                  .AND.                     &
                      BTEST(topo_flags_alt(k-1,j,i),0)                   .AND.                     &
                      BTEST(topo_flags_alt(k,j,i),0)                     .AND.                     &
                      BTEST(topo_flags_alt(k+1,j,i),0)                   .AND.                     &
                      BTEST(topo_flags_alt(k_pp,j,i),0)                  .AND.                     &
                      BTEST(topo_flags_alt(k_ppp,j,i),0)                 .AND.                     &
                     .NOT. flag_set )                                                              &
             THEN
                advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 8 )
             ENDIF

          ENDDO
!
!--       Special treatment for k = nzb in case of elevated nesting. If lower boundary is "open",
!--       set 1st order flag.
          k = nzb
          IF ( BTEST( topo_flags_alt(k,j,i), 0 )  .AND.  BTEST( topo_flags_alt(k+1,j,i), 0 ) )  THEN
             advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 6 )
          ENDIF

       ENDDO
    ENDDO
!
!-- Lateral boundary grid points need special treatment if an alternative communicator is used
!-- which indicates cyclic boundaries, though the boundaries are non-cyclic for the
!-- momentum, e.g., a nested domain with cyclic conditions for the chemistry. Consider the
!-- following situation: a flow from the north, no topography on the northern inflow boundary,
!-- but a hill on the southern outflow boundary boundary. The inflow of momentum is not disturb
!-- by any topography. However, due to the artificial cyclic conditions for the chemistry,
!-- topography is considered at the northern inflow boundary for these scalars. The mismatch
!-- between momentum and scalar boundary treatment can yield to advection of zero scalar values
!-- from the northern boundary since the flow components do not correspond to the artificial
!-- topography at the northern boundary considered for the scalar (there is topography but the
!-- v-component is not zero at the corresponding grid point). This situation is also not
!-- considered in the advc_flag array. To avoid such effects, set the chemistry specific
!-- advection flags to zero for these special situations.
    IF ( PRESENT( alternative_communicator ) )  THEN
!
!--    Treat left boundary.
       IF ( .NOT. non_cyclic_l  .AND.  bc_dirichlet_l )  THEN
          i = nxl
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                IF ( BTEST( topo_flags_alt(k,j,i), 0 ) )  THEN
!
!--                Topography at i-1 and first order along x - disable advection.
                   IF ( .NOT. BTEST( topo_flags_alt(k,j,i-1), 0 )  .AND.                           &
                        BTEST( advc_flag(k,j,i), 0 ) )                                             &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 0 )
!
!--                Topography at i-1 and first order along x - reduce to first order.
                   ELSEIF ( .NOT.  BTEST( topo_flags_alt(k,j,i-2), 0 )  .AND.                      &
                            BTEST( advc_flag(k,j,i), 1 ) )                                         &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 1 )
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 0 )
!
!--                Topography at i-1 and first order along x - reduce to third order.
                   ELSEIF ( .NOT.  BTEST( topo_flags_alt(k,j,i-3), 0 )  .AND.                      &
                            BTEST( advc_flag(k,j,i), 2 ) )                                         &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 2 )
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 1 )
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
!
!--    Treat right boundary.
       IF ( .NOT. non_cyclic_r  .AND.  bc_dirichlet_r )  THEN
          i = nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                IF ( BTEST( topo_flags_alt(k,j,i), 0 ) )  THEN
!
!--                Topography at i+1 and first order along x - disable advection.
                   IF ( .NOT. BTEST( topo_flags_alt(k,j,i+1), 0 )  .AND.                           &
                        BTEST( advc_flag(k,j,i), 0 ) )                                             &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 0 )
!
!--                Topography at i+1 and first order along x - reduce to first order.
                   ELSEIF ( .NOT.  BTEST( topo_flags_alt(k,j,i+2), 0 )  .AND.                      &
                            BTEST( advc_flag(k,j,i), 1 ) )                                         &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 1 )
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 0 )
!
!--                Topography at i+1 and first order along x - reduce to third order.
                   ELSEIF ( .NOT.  BTEST( topo_flags_alt(k,j,i+3), 0 )  .AND.                      &
                            BTEST( advc_flag(k,j,i), 2 ) )                                         &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 2 )
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 1 )
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
!
!--    Treat south boundary.
       IF ( .NOT. non_cyclic_s  .AND.  bc_dirichlet_s )  THEN
          j = nys
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                IF ( BTEST( topo_flags_alt(k,j,i), 0 ) )  THEN
!
!--                Topography at j-1 and first order along y - disable advection.
                   IF ( .NOT. BTEST( topo_flags_alt(k,j-1,i), 0 )  .AND.                           &
                        BTEST( advc_flag(k,j,i), 3 ) )                                             &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 3 )
!
!--                Topography at j-1 and third order along y - reduce to first order.
                   ELSEIF ( .NOT.  BTEST( topo_flags_alt(k,j-2,i), 0 )  .AND.                      &
                            BTEST( advc_flag(k,j,i), 4 ) )                                         &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 4 )
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 3 )
!
!--                Topography at j-1 and fifth order along y - reduce to third order.
                   ELSEIF ( .NOT.  BTEST( topo_flags_alt(k,j-3,i), 0 )  .AND.                      &
                            BTEST( advc_flag(k,j,i), 5 ) )                                         &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 5 )
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 4 )
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
!
!--    Treat north boundary.
       IF ( .NOT. non_cyclic_n  .AND.  bc_dirichlet_n )  THEN
          j = nyn
          DO  i = nxl, nxr
             DO  k = nzb+1, nzt
                IF ( BTEST( topo_flags_alt(k,j,i), 0 ) )  THEN
!
!--                Topography at j+1 and first order along y - disable advection.
                   IF ( .NOT. BTEST( topo_flags_alt(k,j+1,i), 0 )  .AND.                           &
                        BTEST( advc_flag(k,j,i), 3 ) )                                             &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 3 )
!
!--                Topography at j+1 and third order along y - reduce to first order.
                   ELSEIF ( .NOT.  BTEST( topo_flags_alt(k,j+2,i), 0 )  .AND.                      &
                            BTEST( advc_flag(k,j,i), 4 ) )                                         &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 4 )
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 3 )
!
!--                Topography at j+1 and fifth order along y - reduce to third order.
                   ELSEIF ( .NOT.  BTEST( topo_flags_alt(k,j+3,i), 0 )  .AND.                      &
                            BTEST( advc_flag(k,j,i), 5 ) )                                         &
                   THEN
                      advc_flag(k,j,i) = IBCLR( advc_flag(k,j,i), 5 )
                      advc_flag(k,j,i) = IBSET( advc_flag(k,j,i), 4 )
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF

    ENDIF
!
!-- Exchange ghost points for advection flags.
    IF ( PRESENT( alternative_communicator ) )  THEN
       CALL exchange_horiz_int( advc_flag, nys, nyn, nxl, nxr, nzt, nbgp,                          &
                                alternative_communicator = alternative_communicator )
    ELSE
       CALL exchange_horiz_int( advc_flag, nys, nyn, nxl, nxr, nzt, nbgp )
    ENDIF

!
!-- Set the flags at inflow and outflow boundary in case of non-cyclic boundary conditions.
    IF ( non_cyclic_l )  THEN
       advc_flag(:,:,nxl-1) = advc_flag(:,:,nxl)
    ENDIF

    IF ( non_cyclic_r )  THEN
      advc_flag(:,:,nxr+1) = advc_flag(:,:,nxr)
    ENDIF

    IF ( non_cyclic_n )  THEN
       advc_flag(:,nyn+1,:) = advc_flag(:,nyn,:)
    ENDIF

    IF ( non_cyclic_s )  THEN
       advc_flag(:,nys-1,:) = advc_flag(:,nys,:)
    ENDIF

    DEALLOCATE( topo_flags_alt )

 END SUBROUTINE ws_init_flags_scalar


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize variables used for storing statistic quantities (fluxes, variances)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ws_statistics


!
!-- The arrays needed for statistical evaluation are set to to 0 at the beginning of
!-- prognostic_equations.
    IF ( ws_scheme_mom )  THEN
       !$ACC KERNELS PRESENT(sums_wsus_ws_l, sums_wsvs_ws_l) &
       !$ACC PRESENT(sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l)
       sums_wsus_ws_l = 0.0_wp
       sums_wsvs_ws_l = 0.0_wp
       sums_us2_ws_l  = 0.0_wp
       sums_vs2_ws_l  = 0.0_wp
       sums_ws2_ws_l  = 0.0_wp
       !$ACC END KERNELS
    ENDIF

    IF ( ws_scheme_sca )  THEN
       !$ACC KERNELS PRESENT(sums_wspts_ws_l)
       sums_wspts_ws_l = 0.0_wp
       !$ACC END KERNELS
       IF ( humidity       )  sums_wsqs_ws_l = 0.0_wp
       IF ( passive_scalar )  sums_wsss_ws_l = 0.0_wp

    ENDIF

 END SUBROUTINE ws_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Scalar advection - Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_s_ws_ij( advc_flag, i, j, sk, sk_char, swap_flux_y_local, swap_diss_y_local,     &
                           swap_flux_x_local, swap_diss_x_local, i_omp, tn, non_cyclic_l,          &
                           non_cyclic_n, non_cyclic_r, non_cyclic_s, flux_limitation )


    CHARACTER (LEN = *), INTENT(IN) ::  sk_char !< string identifier, used for assign fluxes to the
                                                !<correct dimension in the analysis array

    INTEGER(iwp) ::  i         !< grid index along x-direction
    INTEGER(iwp) ::  i_omp     !< leftmost index on subdomain, or in case of OpenMP, on thread
    INTEGER(iwp) ::  j         !< grid index along y-direction
    INTEGER(iwp) ::  k         !< grid index along z-direction
    INTEGER(iwp) ::  k_mm      !< k-2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_mmm     !< k-3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_pp      !< k+2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_ppp     !< k+3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  nzb_max_l !< index indicating upper bound for order degradation of horizontal advection terms
    INTEGER(iwp) ::  tn        !< number of OpenMP thread

    INTEGER(iwp), INTENT(IN), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::                          &
                                        advc_flag !< flag array to control order of scalar advection

    LOGICAL           ::  limiter         !< control flag indicating the application of flux limitation
    LOGICAL           ::  non_cyclic_l    !< flag that indicates non-cyclic boundary on the left
    LOGICAL           ::  non_cyclic_n    !< flag that indicates non-cyclic boundary on the north
    LOGICAL           ::  non_cyclic_r    !< flag that indicates non-cyclic boundary on the right
    LOGICAL           ::  non_cyclic_s    !< flag that indicates non-cyclic boundary on the south
    LOGICAL, OPTIONAL ::  flux_limitation !< flag indicating flux limitation of the vertical advection

    REAL(wp) ::  diss_d        !< artificial dissipation term at grid box bottom
    REAL(wp) ::  div           !< velocity diverence on scalar grid
    REAL(wp) ::  div_in        !< vertical flux divergence of ingoing fluxes
    REAL(wp) ::  div_out       !< vertical flux divergence of outgoing fluxes
    REAL(wp) ::  f_corr_d      !< correction flux at grid-cell bottom, i.e. the difference between high and low-order flux
    REAL(wp) ::  f_corr_t      !< correction flux at grid-cell top, i.e. the difference between high and low-order flux
    REAL(wp) ::  f_corr_d_in   !< correction flux of ingoing flux part at grid-cell bottom
    REAL(wp) ::  f_corr_t_in   !< correction flux of ingoing flux part at grid-cell top
    REAL(wp) ::  f_corr_d_out  !< correction flux of outgoing flux part at grid-cell bottom
    REAL(wp) ::  f_corr_t_out  !< correction flux of outgoing flux part at grid-cell top
    REAL(wp) ::  fac_correction!< factor to limit the in- and outgoing fluxes
    REAL(wp) ::  flux_d        !< 6th-order flux at grid box bottom
    REAL(wp) ::  ibit0         !< flag indicating 1st-order scheme along x-direction
    REAL(wp) ::  ibit1         !< flag indicating 3rd-order scheme along x-direction
    REAL(wp) ::  ibit2         !< flag indicating 5th-order scheme along x-direction
    REAL(wp) ::  ibit3         !< flag indicating 1st-order scheme along y-direction
    REAL(wp) ::  ibit4         !< flag indicating 3rd-order scheme along y-direction
    REAL(wp) ::  ibit5         !< flag indicating 5th-order scheme along y-direction
    REAL(wp) ::  ibit6         !< flag indicating 1st-order scheme along z-direction
    REAL(wp) ::  ibit7         !< flag indicating 3rd-order scheme along z-direction
    REAL(wp) ::  ibit8         !< flag indicating 5th-order scheme along z-direction
    REAL(wp) ::  max_val       !< maximum value of the quanitity along the numerical stencil (in vertical direction)
    REAL(wp) ::  min_val       !< maximum value of the quanitity along the numerical stencil (in vertical direction)
    REAL(wp) ::  mon           !< monotone solution of the advection equation using 1st-order fluxes
    REAL(wp) ::  u_comp        !< advection velocity along x-direction
    REAL(wp) ::  v_comp        !< advection velocity along y-direction

    REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_n     !< discretized artificial dissipation at northward-side
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_r     !< discretized artificial dissipation at rightward-side
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_t     !< discretized artificial dissipation at top
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_n     !< discretized 6th-order flux at northward-side
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_r     !< discretized 6th-order flux at rightward-side
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_t     !< discretized 6th-order flux at top
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_t_1st !< discretized 1st-order flux at top

    REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1) ::  swap_diss_y_local !< discretized artificial dissipation at southward-side
    REAL(wp), DIMENSION(nzb+1:nzt,0:threads_per_task-1) ::  swap_flux_y_local !< discretized 6th-order flux at northward-side
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) ::  swap_diss_x_local !< discretized artificial dissipation at leftward-side
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,0:threads_per_task-1) ::  swap_flux_x_local !< discretized 6th-order flux at leftward-side

!
!-- sk is an array from parameter list. It should not be a pointer, because in that case the
!-- compiler can not assume a stride 1 and cannot perform a strided one vector load. Adding the
!-- CONTIGUOUS keyword makes things even worse, because the compiler cannot assume strided one in
!-- the caller side.
    REAL(wp), INTENT(IN),DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !< advected scalar

!
!-- Used local modified copy of nzb_max (used to degrade order of discretization) at non-cyclic
!-- boundaries. Modify only at relevant points instead of the entire subdomain. This should lead to
!-- betterload balance between boundary and non-boundary PEs.
    IF( non_cyclic_l  .AND.  i <= nxl + 2  .OR.                                                    &
        non_cyclic_r  .AND.  i >= nxr - 2  .OR.                                                    &
        non_cyclic_s  .AND.  j <= nys + 2  .OR.                                                    &
        non_cyclic_n  .AND.  j >= nyn - 2 )  THEN
       nzb_max_l = nzt
    ELSE
       nzb_max_l = nzb_max
    END IF
!
!-- Set control flag for flux limiter
    limiter = .FALSE.
    IF ( PRESENT( flux_limitation) )  limiter = flux_limitation
!
!-- Compute southside fluxes of the respective PE bounds.
    IF ( j == nys )  THEN
!
!--    Up to the top of the highest topography.
       DO  k = nzb+1, nzb_max_l

          ibit5 = REAL( IBITS(advc_flag(k,j-1,i),5,1), KIND = wp )
          ibit4 = REAL( IBITS(advc_flag(k,j-1,i),4,1), KIND = wp )
          ibit3 = REAL( IBITS(advc_flag(k,j-1,i),3,1), KIND = wp )
!
!--       Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
          v_comp                  = v(k,j,i) - v_gtrans + v_stokes_zu(k)                           &
                                    * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          swap_flux_y_local(k,tn) = v_comp *         (                                             &
                                            ( 37.0_wp * ibit5 * adv_sca_5                          &
                                         +     7.0_wp * ibit4 * adv_sca_3                          &
                                         +              ibit3 * adv_sca_1                          &
                                            ) * ( sk(k,j,i)   + sk(k,j-1,i) )                      &
                                         -  (  8.0_wp * ibit5 * adv_sca_5                          &
                                         +              ibit4 * adv_sca_3                          &
                                            ) * ( sk(k,j+1,i) + sk(k,j-2,i) )                      &
                                         +  (           ibit5 * adv_sca_5 )                        &
                                              * ( sk(k,j+2,i) + sk(k,j-3,i) )                      &
                                                     )

          swap_diss_y_local(k,tn) = -ABS( v_comp ) * (                                             &
                                            ( 10.0_wp * ibit5 * adv_sca_5                          &
                                         +     3.0_wp * ibit4 * adv_sca_3                          &
                                         +              ibit3 * adv_sca_1                          &
                                            ) * ( sk(k,j,i)   - sk(k,j-1,i) )                      &
                                         -  ( 5.0_wp  * ibit5 * adv_sca_5                          &
                                         +              ibit4 * adv_sca_3                          &
                                            ) * ( sk(k,j+1,i) - sk(k,j-2,i) )                      &
                                         +  (           ibit5 * adv_sca_5   )                      &
                                              * ( sk(k,j+2,i) - sk(k,j-3,i)  )                     &
                                                     )

       ENDDO
!
!--    Above to the top of the highest topography. No degradation necessary.
       DO  k = nzb_max_l+1, nzt

          v_comp                  = v(k,j,i) - v_gtrans + v_stokes_zu(k)
          swap_flux_y_local(k,tn) = v_comp * ( 37.0_wp * ( sk(k,j,i)   + sk(k,j-1,i) )             &
                                              - 8.0_wp * ( sk(k,j+1,i) + sk(k,j-2,i) )             &
                                              + ( sk(k,j+2,i) + sk(k,j-3,i) )                      &
                                             ) * adv_sca_5
          swap_diss_y_local(k,tn) = -ABS( v_comp ) * (                                             &
                                     10.0_wp * ( sk(k,j,i)   - sk(k,j-1,i) )                       &
                                    - 5.0_wp * ( sk(k,j+1,i) - sk(k,j-2,i) )                       &
                                    +            sk(k,j+2,i) - sk(k,j-3,i)                         &
                                                     ) * adv_sca_5

       ENDDO

    ENDIF
!
!-- Compute leftside fluxes of the respective PE bounds.
    IF ( i == i_omp )  THEN

       DO  k = nzb+1, nzb_max_l

          ibit2 = REAL( IBITS(advc_flag(k,j,i-1),2,1), KIND = wp )
          ibit1 = REAL( IBITS(advc_flag(k,j,i-1),1,1), KIND = wp )
          ibit0 = REAL( IBITS(advc_flag(k,j,i-1),0,1), KIND = wp )
!
!--       Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
          u_comp                    = u(k,j,i) - u_gtrans + u_stokes_zu(k)                         &
                                      * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          swap_flux_x_local(k,j,tn) = u_comp * (                                                   &
                                            ( 37.0_wp * ibit2 * adv_sca_5                          &
                                         +     7.0_wp * ibit1 * adv_sca_3                          &
                                         +              ibit0 * adv_sca_1                          &
                                            ) *   ( sk(k,j,i) + sk(k,j,i-1) )                      &
                                         -  (  8.0_wp * ibit2 * adv_sca_5                          &
                                         +              ibit1 * adv_sca_3                          &
                                            ) * ( sk(k,j,i+1) + sk(k,j,i-2) )                      &
                                         +  (           ibit2 * adv_sca_5                          &
                                            ) * ( sk(k,j,i+2) + sk(k,j,i-3) )                      &
                                               )

          swap_diss_x_local(k,j,tn) = -ABS( u_comp ) * (                                           &
                                            ( 10.0_wp * ibit2 * adv_sca_5                          &
                                         +     3.0_wp * ibit1 * adv_sca_3                          &
                                         +              ibit0 * adv_sca_1                          &
                                            ) *   ( sk(k,j,i) - sk(k,j,i-1) )                      &
                                         -  (  5.0_wp * ibit2 * adv_sca_5                          &
                                         +              ibit1 * adv_sca_3                          &
                                            ) * ( sk(k,j,i+1) - sk(k,j,i-2) )                      &
                                         +  (           ibit2 * adv_sca_5                          &
                                            ) * ( sk(k,j,i+2) - sk(k,j,i-3)    )                   &
                                                      )

       ENDDO

       DO  k = nzb_max_l+1, nzt

          u_comp                    = u(k,j,i) - u_gtrans + u_stokes_zu(k)
          swap_flux_x_local(k,j,tn) = u_comp  * (                                                  &
                                        37.0_wp * ( sk(k,j,i) + sk(k,j,i-1) )                      &
                                      -  8.0_wp * ( sk(k,j,i+1) + sk(k,j,i-2) )                    &
                                      +           ( sk(k,j,i+2) + sk(k,j,i-3) )                    &
                                                ) * adv_sca_5

          swap_diss_x_local(k,j,tn) = -ABS( u_comp ) * (                                           &
                                        10.0_wp * ( sk(k,j,i)   - sk(k,j,i-1) )                    &
                                       - 5.0_wp * ( sk(k,j,i+1) - sk(k,j,i-2) )                    &
                                       +          ( sk(k,j,i+2) - sk(k,j,i-3) )                    &
                                                       ) * adv_sca_5

       ENDDO

    ENDIF
!
!-- Now compute the fluxes for the horizontal termns up to the highest
!-- topography.
    DO  k = nzb+1, nzb_max_l

       ibit2 = REAL( IBITS(advc_flag(k,j,i),2,1), KIND = wp )
       ibit1 = REAL( IBITS(advc_flag(k,j,i),1,1), KIND = wp )
       ibit0 = REAL( IBITS(advc_flag(k,j,i),0,1), KIND = wp )
!
!--    Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
       u_comp    = u(k,j,i+1) - u_gtrans + u_stokes_zu(k)                                          &
                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 1 ) )
       flux_r(k) = u_comp * (                                                                      &
                   (  37.0_wp * ibit2 * adv_sca_5                                                  &
                   +   7.0_wp * ibit1 * adv_sca_3                                                  &
                   +            ibit0 * adv_sca_1 ) * ( sk(k,j,i+1) + sk(k,j,i)   )                &
                   - ( 8.0_wp * ibit2 * adv_sca_5                                                  &
                   +            ibit1 * adv_sca_3 ) * ( sk(k,j,i+2) + sk(k,j,i-1) )                &
                   + (          ibit2 * adv_sca_5 ) * ( sk(k,j,i+3) + sk(k,j,i-2) )                &
                            )

       diss_r(k) = -ABS( u_comp ) * (                                                              &
                   (  10.0_wp * ibit2 * adv_sca_5                                                  &
                   +   3.0_wp * ibit1 * adv_sca_3                                                  &
                   +            ibit0 * adv_sca_1 ) * ( sk(k,j,i+1) - sk(k,j,i) )                  &
                   - ( 5.0_wp * ibit2 * adv_sca_5                                                  &
                   +            ibit1 * adv_sca_3 ) * ( sk(k,j,i+2) - sk(k,j,i-1) )                &
                   + (          ibit2 * adv_sca_5 ) * ( sk(k,j,i+3) - sk(k,j,i-2) )                &
                                    )

       ibit5 = REAL( IBITS(advc_flag(k,j,i),5,1), KIND = wp )
       ibit4 = REAL( IBITS(advc_flag(k,j,i),4,1), KIND = wp )
       ibit3 = REAL( IBITS(advc_flag(k,j,i),3,1), KIND = wp )
!
!--    Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
       v_comp    = v(k,j+1,i) - v_gtrans + v_stokes_zu(k)                                          &
                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j+1,i), 2 ) )
       flux_n(k) = v_comp * (                                                                      &
                   (  37.0_wp * ibit5 * adv_sca_5                                                  &
                   +   7.0_wp * ibit4 * adv_sca_3                                                  &
                   +            ibit3 * adv_sca_1 ) * ( sk(k,j+1,i) + sk(k,j,i)   )                &
                   - ( 8.0_wp * ibit5 * adv_sca_5                                                  &
                   +            ibit4 * adv_sca_3 ) * ( sk(k,j+2,i) + sk(k,j-1,i) )                &
                   + (          ibit5 * adv_sca_5 ) * ( sk(k,j+3,i) + sk(k,j-2,i) )                &
                            )

       diss_n(k) = -ABS( v_comp ) * (                                                              &
                   (  10.0_wp * ibit5 * adv_sca_5                                                  &
                   +   3.0_wp * ibit4 * adv_sca_3                                                  &
                   +            ibit3 * adv_sca_1 ) * ( sk(k,j+1,i) - sk(k,j,i)   )                &
                   - ( 5.0_wp * ibit5 * adv_sca_5                                                  &
                   +            ibit4 * adv_sca_3 ) * ( sk(k,j+2,i) - sk(k,j-1,i) )                &
                   + (          ibit5 * adv_sca_5 ) * ( sk(k,j+3,i) - sk(k,j-2,i) )                &
                                    )
    ENDDO
!
!-- Now compute the fluxes for the horizontal terms above the topography
!-- where no degradation along the horizontal parts is necessary (except
!-- for the non-cyclic lateral boundaries treated by nzb_max_l).
    DO  k = nzb_max_l+1, nzt

       u_comp    = u(k,j,i+1) - u_gtrans + u_stokes_zu(k)
       flux_r(k) = u_comp * (                                                                      &
                     37.0_wp * ( sk(k,j,i+1) + sk(k,j,i)   )                                       &
                   -  8.0_wp * ( sk(k,j,i+2) + sk(k,j,i-1) )                                       &
                   +           ( sk(k,j,i+3) + sk(k,j,i-2) ) ) * adv_sca_5
       diss_r(k) = -ABS( u_comp ) * (                                                              &
                     10.0_wp * ( sk(k,j,i+1) - sk(k,j,i)   )                                       &
                   -  5.0_wp * ( sk(k,j,i+2) - sk(k,j,i-1) )                                       &
                   +           ( sk(k,j,i+3) - sk(k,j,i-2) ) ) * adv_sca_5

       v_comp    = v(k,j+1,i) - v_gtrans + v_stokes_zu(k)
       flux_n(k) = v_comp * (                                                                      &
                     37.0_wp * ( sk(k,j+1,i) + sk(k,j,i)   )                                       &
                   -  8.0_wp * ( sk(k,j+2,i) + sk(k,j-1,i) )                                       &
                   +           ( sk(k,j+3,i) + sk(k,j-2,i) ) ) * adv_sca_5
       diss_n(k) = -ABS( v_comp ) * (                                                              &
                     10.0_wp * ( sk(k,j+1,i) - sk(k,j,i)   )                                       &
                   -  5.0_wp * ( sk(k,j+2,i) - sk(k,j-1,i) )                                       &
                   +           ( sk(k,j+3,i) - sk(k,j-2,i) ) ) * adv_sca_5

    ENDDO
!
!-- Now, compute vertical fluxes. Split loop into a part treating the lowest grid points with
!-- indirect indexing, a main loop without indirect indexing, and a loop for the uppermost grid
!-- points with indirect indexing. This allows better vectorization for the main loop.
!-- First, compute the flux at model surface, which need has to be calculated explicetely for the
!-- tendency at the first w-level. For topography wall this is done implicitely by advc_flag.
    k = nzb
    ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )
    flux_t(nzb) = w(k,j,i)         * rho_air_zw(k) * ( ibit6 * adv_sca_1 )                         &
                                                   * ( sk(k+1,j,i) + sk(k,j,i) )
    diss_t(nzb) = -ABS( w(k,j,i) ) * rho_air_zw(k) * ( ibit6 * adv_sca_1 )                         &
                                                   * ( sk(k+1,j,i) - sk(k,j,i) )

    DO  k = nzb+1, nzb+1
       ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
       ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
       ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )
!
!--    k index has to be modified near bottom and top, else array subscripts will be exceeded.
       k_ppp = k + 3 * ibit8
       k_pp  = k + 2 * ( 1 - ibit6 )
       k_mm  = k - 2 * ibit8

       flux_t(k) = w(k,j,i)  * rho_air_zw(k) * (                                                   &
                   (  37.0_wp * ibit8 * adv_sca_5                                                  &
                   +   7.0_wp * ibit7 * adv_sca_3                                                  &
                   +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)  + sk(k,j,i)    )              &
                   - ( 8.0_wp * ibit8 * adv_sca_5                                                  &
                   +            ibit7 * adv_sca_3 ) * ( sk(k_pp,j,i) + sk(k-1,j,i)  )              &
                   + (          ibit8 * adv_sca_5 ) * ( sk(k_ppp,j,i)+ sk(k_mm,j,i) )              &
                                               )

       diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                                            &
                   (  10.0_wp * ibit8 * adv_sca_5                                                  &
                   +   3.0_wp * ibit7 * adv_sca_3                                                  &
                   +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)   - sk(k,j,i)    )             &
                   - ( 5.0_wp * ibit8 * adv_sca_5                                                  &
                   +            ibit7 * adv_sca_3 ) * ( sk(k_pp,j,i)  - sk(k-1,j,i)  )             &
                   + (          ibit8 * adv_sca_5 ) * ( sk(k_ppp,j,i) - sk(k_mm,j,i) )             &
                                                      )
    ENDDO

    DO  k = nzb+2, nzt-2
       ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
       ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
       ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )

       flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                                                    &
                   (  37.0_wp * ibit8 * adv_sca_5                                                  &
                   +   7.0_wp * ibit7 * adv_sca_3                                                  &
                   +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)  + sk(k,j,i)   )               &
                   - ( 8.0_wp * ibit8 * adv_sca_5                                                  &
                   +            ibit7 * adv_sca_3 ) * ( sk(k+2,j,i)  + sk(k-1,j,i) )               &
                   + (          ibit8 * adv_sca_5 ) * ( sk(k+3,j,i)  + sk(k-2,j,i) )               &
                                              )

       diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                                            &
                   (  10.0_wp * ibit8 * adv_sca_5                                                  &
                   +   3.0_wp * ibit7 * adv_sca_3                                                  &
                   +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i) - sk(k,j,i)   )                &
                   - ( 5.0_wp * ibit8 * adv_sca_5                                                  &
                   +            ibit7 * adv_sca_3 ) * ( sk(k+2,j,i) - sk(k-1,j,i) )                &
                   + (          ibit8 * adv_sca_5 ) * ( sk(k+3,j,i) - sk(k-2,j,i) )                &
                                                      )
    ENDDO

    DO  k = nzt-1, nzt-symmetry_flag
       ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
       ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
       ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )
!
!--    k index needs to be modified near bottom and top, else array subscripts will be exceeded.
!--    (k_ppp = k+3, k_mm = k) if 5th order scheme is employed, else (k_ppp = k, k_mm = k).
!--    Special treatment is required for k_pp. This must be limited at k = nzt, where either
!--    the first order scheme is employed (ibit6), or all bit flags are zero (in case there is a
!--    rigid lid defined at or near the domain top).
       k_ppp = k + 3 * ibit8
       k_pp  = k + 2 * ( 1 - ibit6  ) * ( ibit6 + ibit7 + ibit8 )
       k_mm  = k - 2 * ibit8


       flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                                                    &
                   (  37.0_wp * ibit8 * adv_sca_5                                                  &
                   +   7.0_wp * ibit7 * adv_sca_3                                                  &
                   +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)  + sk(k,j,i)    )              &
                   - ( 8.0_wp * ibit8 * adv_sca_5                                                  &
                   +            ibit7 * adv_sca_3 ) * ( sk(k_pp,j,i) + sk(k-1,j,i)  )              &
                   + (          ibit8 * adv_sca_5 ) * ( sk(k_ppp,j,i)+ sk(k_mm,j,i) )              &
                                              )

       diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                                            &
                   (  10.0_wp * ibit8 * adv_sca_5                                                  &
                   +   3.0_wp * ibit7 * adv_sca_3                                                  &
                   +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)   - sk(k,j,i)    )             &
                   - ( 5.0_wp * ibit8 * adv_sca_5                                                  &
                   +            ibit7 * adv_sca_3 ) * ( sk(k_pp,j,i)  - sk(k-1,j,i)  )             &
                   + (          ibit8 * adv_sca_5 ) * ( sk(k_ppp,j,i) - sk(k_mm,j,i) )             &
                                                      )
    ENDDO

!
!-- Set resolved/turbulent flux at model top to zero (w-level). In case that a symmetric behavior
!-- between bottom and top shall be guaranteed (closed channel flow), the flux at nzt is also set to
!-- zero.
    IF ( symmetry_flag == 1 ) THEN
       flux_t(nzt) = 0.0_wp
       diss_t(nzt) = 0.0_wp
    ENDIF
    flux_t(nzt+1) = 0.0_wp
    diss_t(nzt+1) = 0.0_wp


    IF ( limiter )  THEN
!
!--    Compute monotone first-order fluxes which are required for mononteflux limitation.
       flux_t_1st(nzb) = 0.0_wp
       DO  k = nzb+1, nzb_max_l
          flux_t_1st(k) = (      w(k,j,i)   * ( sk(k+1,j,i) + sk(k,j,i) )                          &
                           -ABS( w(k,j,i) ) * ( sk(k+1,j,i) - sk(k,j,i) ) )                        &
                           * rho_air_zw(k) * adv_sca_1
!
!--       In flux limitation the total flux will be corrected. For the sake of cleariness the
!--       higher-order advective and disspative fluxes will be merged onto flux_t.
          flux_t(k) = flux_t(k) + diss_t(k)
          diss_t(k) = 0.0_wp
       ENDDO
!
!--    Flux limitation of vertical fluxes according to Skamarock (2006).
!--    Please note, as flux limitation implies linear dependencies of fluxes, flux limitation is
!--    only made for the vertical advection term. Limitation of the horizontal terms cannot be
!--    parallelized.
!--    Due to the linear dependency, the following loop will not be vectorized.
!--    Further, note that the flux limiter is only applied within the urban layer, i.e up to the
!--    topography top.
       DO  k = nzb+1, nzb_max_l
!
!--       Compute one-dimensional divergence along the vertical direction, which is used to correct
!--       the advection discretization. This is necessary as in one-dimensional space the advection
!--       velocity should actually be constant.
          div = ( w(k,j,i)   * rho_air_zw(k)                                                       &
                - w(k-1,j,i) * rho_air_zw(k-1)                                                     &
                ) * drho_air(k) * ddzw(k)
!
!--       Compute monotone solution of the advection equation from 1st-order fluxes. Please note,
!--       the advection equation is corrected by the divergence term (in 1D the advective flow
!--       should be divergence free). Moreover, please note, as time-increment the full timestep is
!--       used, even though a Runge-Kutta scheme will be used. However, the length of the actual
!--       time increment is not important at all since it cancels out later when the fluxes are
!--       limited.
          mon = sk(k,j,i) + ( - ( flux_t_1st(k) - flux_t_1st(k-1) )                                &
                          * drho_air(k) * ddzw(k)                                                  &
                          + div * sk(k,j,i)                                                        &
                            ) * dt_3d
!
!--       Determine minimum and maximum values along the numerical stencil.
          k_mmm = MAX( k - 3, nzb + 1 )
          k_ppp = MIN( k + 3, nzt + 1 )

          min_val = MINVAL( sk(k_mmm:k_ppp,j,i) )
          max_val = MAXVAL( sk(k_mmm:k_ppp,j,i) )
!
!--       Compute difference between high- and low-order fluxes, which may act as correction fluxes
          f_corr_t = flux_t(k)   - flux_t_1st(k)
          f_corr_d = flux_t(k-1) - flux_t_1st(k-1)
!
!--       Determine outgoing fluxes, i.e. the part of the fluxes which can decrease the value within
!--       the grid box
          f_corr_t_out = MAX( 0.0_wp, f_corr_t )
          f_corr_d_out = MIN( 0.0_wp, f_corr_d )
!
!--       Determine ingoing fluxes, i.e. the part of the fluxes which can  increase the value within
!--       the grid box
          f_corr_t_in = MIN( 0.0_wp, f_corr_t)
          f_corr_d_in = MAX( 0.0_wp, f_corr_d)
!
!--       Compute divergence of outgoing correction fluxes
          div_out = - ( f_corr_t_out - f_corr_d_out ) * drho_air(k) * ddzw(k) * dt_3d
!
!--       Compute divergence of ingoing correction fluxes
          div_in = - ( f_corr_t_in - f_corr_d_in ) * drho_air(k) * ddzw(k) * dt_3d
!
!--       Check if outgoing fluxes can lead to undershoots, i.e. values smaller than the minimum
!--       value within the numerical stencil. If so, limit them.
          IF ( mon - min_val < - div_out  .AND.  ABS( div_out ) > 0.0_wp )  THEN
             fac_correction = ( mon - min_val ) / ( - div_out )
             f_corr_t_out = f_corr_t_out * fac_correction
             f_corr_d_out = f_corr_d_out * fac_correction
          ENDIF
!
!--       Check if ingoing fluxes can lead to overshoots, i.e. values larger than the maximum value
!--       within the numerical stencil. If so, limit them.
          IF ( mon - max_val > - div_in  .AND.  ABS( div_in ) > 0.0_wp )  THEN
             fac_correction = ( mon - max_val ) / ( - div_in )
             f_corr_t_in = f_corr_t_in * fac_correction
             f_corr_d_in = f_corr_d_in * fac_correction
          ENDIF
!
!--       Finally add the limited fluxes to the original ones. If no flux limitation was done, the
!--       fluxes equal the original ones.
          flux_t(k)   = flux_t_1st(k)   + f_corr_t_out + f_corr_t_in
          flux_t(k-1) = flux_t_1st(k-1) + f_corr_d_out + f_corr_d_in
       ENDDO
    ENDIF
!
!-- Now compute the tendency term including divergence correction.
    DO  k = nzb+1, nzb_max_l

       flux_d    = flux_t(k-1)
       diss_d    = diss_t(k-1)

       ibit2 = REAL( IBITS(advc_flag(k,j,i),2,1), KIND = wp )
       ibit1 = REAL( IBITS(advc_flag(k,j,i),1,1), KIND = wp )
       ibit0 = REAL( IBITS(advc_flag(k,j,i),0,1), KIND = wp )

       ibit5 = REAL( IBITS(advc_flag(k,j,i),5,1), KIND = wp )
       ibit4 = REAL( IBITS(advc_flag(k,j,i),4,1), KIND = wp )
       ibit3 = REAL( IBITS(advc_flag(k,j,i),3,1), KIND = wp )

       ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
       ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
       ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )
!
!--    Calculate the divergence of the velocity field. A respective correction is needed to overcome
!--    numerical instabilities introduced by an insufficient reduction of divergences near
!--    topography.
       div         =   ( u(k,j,i+1) * ( ibit0 + ibit1 + ibit2 )                                    &
                       - u(k,j,i)   * (                                                            &
                       REAL( IBITS(advc_flag(k,j,i-1),0,1), KIND = wp )                            &
                     + REAL( IBITS(advc_flag(k,j,i-1),1,1), KIND = wp )                            &
                     + REAL( IBITS(advc_flag(k,j,i-1),2,1), KIND = wp )                            &
                                      )                                                            &
                       ) * ddx                                                                     &
                     + ( v(k,j+1,i) * ( ibit3 + ibit4 + ibit5 )                                    &
                       - v(k,j,i)   * (                                                            &
                       REAL( IBITS(advc_flag(k,j-1,i),3,1), KIND = wp )                            &
                     + REAL( IBITS(advc_flag(k,j-1,i),4,1), KIND = wp )                            &
                     + REAL( IBITS(advc_flag(k,j-1,i),5,1), KIND = wp )                            &
                                      )                                                            &
                       ) * ddy                                                                     &
                     + ( w(k,j,i)   * rho_air_zw(k) * ( ibit6 + ibit7 + ibit8 )                    &
                       - w(k-1,j,i) * rho_air_zw(k-1) *                                            &
                                      (                                                            &
                       REAL( IBITS(advc_flag(k-1,j,i),6,1), KIND = wp )                            &
                     + REAL( IBITS(advc_flag(k-1,j,i),7,1), KIND = wp )                            &
                     + REAL( IBITS(advc_flag(k-1,j,i),8,1), KIND = wp )                            &
                                      )                                                            &
                       ) * drho_air(k) * ddzw(k)

       tend(k,j,i) = tend(k,j,i) - (                                                               &
                     ( flux_r(k) + diss_r(k) - swap_flux_x_local(k,j,tn)                           &
                     - swap_diss_x_local(k,j,tn) ) * ddx                                           &
                     + ( flux_n(k) + diss_n(k) - swap_flux_y_local(k,tn)                           &
                     - swap_diss_y_local(k,tn)   ) * ddy                                           &
                     + ( ( flux_t(k) + diss_t(k) ) - ( flux_d + diss_d )                           &
                       )             * drho_air(k) * ddzw(k)                                       &
                                   ) + sk(k,j,i) * div


       swap_flux_y_local(k,tn)   = flux_n(k)
       swap_diss_y_local(k,tn)   = diss_n(k)
       swap_flux_x_local(k,j,tn) = flux_r(k)
       swap_diss_x_local(k,j,tn) = diss_r(k)

    ENDDO

    DO  k = nzb_max_l+1, nzt

       flux_d    = flux_t(k-1)
       diss_d    = diss_t(k-1)
!
!--    Calculate the divergence of the velocity field. A respective correction is needed to overcome
!--    numerical instabilities introduced by an insufficient reduction of divergences near
!--    topography.
       div         =   ( u(k,j,i+1) - u(k,j,i) ) * ddx                                             &
                     + ( v(k,j+1,i) - v(k,j,i) ) * ddy                                             &
                     + ( w(k,j,i)   * rho_air_zw(k)                                                &
                       - w(k-1,j,i) * rho_air_zw(k-1)                                              &
                                               ) * drho_air(k) * ddzw(k)

       tend(k,j,i) = tend(k,j,i) - (                                                               &
                     ( flux_r(k) + diss_r(k) - swap_flux_x_local(k,j,tn)                           &
                     - swap_diss_x_local(k,j,tn) ) * ddx                                           &
                     + ( flux_n(k) + diss_n(k) - swap_flux_y_local(k,tn)                           &
                     - swap_diss_y_local(k,tn)   ) * ddy                                           &
                     + ( ( flux_t(k) + diss_t(k) ) - ( flux_d + diss_d )                           &
                                                 ) * drho_air(k) * ddzw(k)                         &
                                   ) + sk(k,j,i) * div


       swap_flux_y_local(k,tn)   = flux_n(k)
       swap_diss_y_local(k,tn)   = diss_n(k)
       swap_flux_x_local(k,j,tn) = flux_r(k)
       swap_diss_x_local(k,j,tn) = diss_r(k)

    ENDDO

!
!-- Evaluation of statistics.
    SELECT CASE ( sk_char )

       CASE ( 'pt' )

          DO  k = nzb+1, nzt
             sums_wspts_ws_l(k,tn) = sums_wspts_ws_l(k,tn) +                                       &
                                     ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )     &
                                                 * ( w(k,j,i) - hom(k,1,3,0)                 )     &
                                     + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )     &
                                                 *   ABS( w(k,j,i) - hom(k,1,3,0)            )     &
                                     ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'sa' )

          DO  k = nzb+1, nzt
             sums_wssas_ws_l(k,tn) = sums_wssas_ws_l(k,tn) +                                       &
                                     ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )     &
                                                 * ( w(k,j,i)    - hom(k,1,3,0)              )     &
                                     + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )     &
                                                 *   ABS( w(k,j,i) - hom(k,1,3,0)            )     &
                                     ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'q' )

          DO  k = nzb+1, nzt
             sums_wsqs_ws_l(k,tn)  = sums_wsqs_ws_l(k,tn) +                                        &
                                     ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )     &
                                                 * ( w(k,j,i) - hom(k,1,3,0)                 )     &
                                     + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )     &
                                                 *   ABS( w(k,j,i) - hom(k,1,3,0)            )     &
                                     ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'qc' )

          DO  k = nzb+1, nzt
             sums_wsqcs_ws_l(k,tn)  = sums_wsqcs_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'qg' )

          DO  k = nzb+1, nzt
             sums_wsqgs_ws_l(k,tn)  = sums_wsqgs_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'qi' )

          DO  k = nzb+1, nzt
             sums_wsqis_ws_l(k,tn)  = sums_wsqis_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'qr' )

          DO  k = nzb+1, nzt
             sums_wsqrs_ws_l(k,tn)  = sums_wsqrs_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'qs' )

          DO  k = nzb+1, nzt
             sums_wsqss_ws_l(k,tn)  = sums_wsqss_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'nc' )

          DO  k = nzb+1, nzt
             sums_wsncs_ws_l(k,tn)  = sums_wsncs_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'ng' )

          DO  k = nzb+1, nzt
             sums_wsngs_ws_l(k,tn)  = sums_wsngs_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO


       CASE ( 'ni' )

          DO  k = nzb+1, nzt
             sums_wsnis_ws_l(k,tn)  = sums_wsnis_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'nr' )

          DO  k = nzb+1, nzt
             sums_wsnrs_ws_l(k,tn)  = sums_wsnrs_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'ns' )

          DO  k = nzb+1, nzt
             sums_wsnss_ws_l(k,tn)  = sums_wsnss_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 's' )

          DO  k = nzb+1, nzt
             sums_wsss_ws_l(k,tn)  = sums_wsss_ws_l(k,tn) +                                        &
                                     ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )     &
                                                 * ( w(k,j,i) - hom(k,1,3,0)                 )     &
                                     + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )     &
                                                 *   ABS( w(k,j,i) - hom(k,1,3,0)            )     &
                                     ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'aerosol_mass', 'aerosol_number', 'salsa_gas' )

          DO  k = nzb+1, nzt
             sums_salsa_ws_l(k,tn)  = sums_salsa_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

       CASE ( 'kc' )
          DO  k = nzb+1, nzt
             sums_wschs_ws_l(k,tn)  = sums_wschs_ws_l(k,tn) +                                      &
                                      ( flux_t(k) / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )    &
                                                  * ( w(k,j,i) - hom(k,1,3,0)                 )    &
                                      + diss_t(k) / ( ABS(w(k,j,i)) + 1.0E-20_wp              )    &
                                                  *   ABS( w(k,j,i) - hom(k,1,3,0)            )    &
                                      ) * weight_substep(intermediate_timestep_count)
          ENDDO

    END SELECT

 END SUBROUTINE advec_s_ws_ij




!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of u-component - Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_u_ws_ij( i, j, i_omp, tn )


    INTEGER(iwp) ::  i         !< grid index along x-direction
    INTEGER(iwp) ::  i_omp     !< leftmost index on subdomain, or in case of OpenMP, on thread
    INTEGER(iwp) ::  j         !< grid index along y-direction
    INTEGER(iwp) ::  k         !< grid index along z-direction
    INTEGER(iwp) ::  k_mm      !< k-2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_pp      !< k+2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_ppp     !< k+3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  nzb_max_l !< index indicating upper bound for order degradation of horizontal advection terms
    INTEGER(iwp) ::  tn        !< number of OpenMP thread

    REAL(wp)    ::  diss_d   !< artificial dissipation term at grid box bottom
    REAL(wp)    ::  div      !< diverence on u-grid
    REAL(wp)    ::  flux_d   !< 6th-order flux at grid box bottom
    REAL(wp)    ::  gu       !< Galilei-transformation velocity along x
    REAL(wp)    ::  gv       !< Galilei-transformation velocity along y
    REAL(wp)    ::  ibit0   !< flag indicating 1st-order scheme along x-direction
    REAL(wp)    ::  ibit1   !< flag indicating 3rd-order scheme along x-direction
    REAL(wp)    ::  ibit2   !< flag indicating 5th-order scheme along x-direction
    REAL(wp)    ::  ibit3   !< flag indicating 1st-order scheme along y-direction
    REAL(wp)    ::  ibit4   !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)    ::  ibit5   !< flag indicating 5th-order scheme along y-direction
    REAL(wp)    ::  ibit6   !< flag indicating 1st-order scheme along z-direction
    REAL(wp)    ::  ibit7   !< flag indicating 3rd-order scheme along z-direction
    REAL(wp)    ::  ibit8   !< flag indicating 5th-order scheme along z-direction
    REAL(wp)    ::  u_comp_l !< advection velocity along x at leftmost grid point on subdomain

    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_n !< discretized artificial dissipation at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_r !< discretized artificial dissipation at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_t !< discretized artificial dissipation at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_n !< discretized 6th-order flux at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_r !< discretized 6th-order flux at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_t !< discretized 6th-order flux at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  u_comp !< advection velocity along x
    REAL(wp), DIMENSION(nzb:nzt+1) ::  v_comp !< advection velocity along y
    REAL(wp), DIMENSION(nzb:nzt+1) ::  w_comp !< advection velocity along z
!
!-- Used local modified copy of nzb_max (used to degrade order of discretization) at non-cyclic
!-- boundaries. Modify only at relevant points instead of the entire subdomain. This should lead to
!-- better load balance between boundary and non-boundary PEs.
    IF( ( bc_dirichlet_l  .OR.  bc_radiation_l )  .AND.  i <= nxl + 2  .OR.                        &
        ( bc_dirichlet_r  .OR.  bc_radiation_r )  .AND.  i >= nxr - 2  .OR.                        &
        ( bc_dirichlet_s  .OR.  bc_radiation_s )  .AND.  j <= nys + 2  .OR.                        &
        ( bc_dirichlet_n  .OR.  bc_radiation_n )  .AND.  j >= nyn - 2 )  THEN
       nzb_max_l = nzt
    ELSE
       nzb_max_l = nzb_max
    END IF

    gu = 2.0_wp * u_gtrans
    gv = 2.0_wp * v_gtrans
!
!-- Compute southside fluxes for the respective boundary of PE
    IF ( j == nys  )  THEN

       DO  k = nzb+1, nzb_max_l

          ibit5 = REAL( IBITS(advc_flags_m(k,j-1,i),5,1), KIND = wp )
          ibit4 = REAL( IBITS(advc_flags_m(k,j-1,i),4,1), KIND = wp )
          ibit3 = REAL( IBITS(advc_flags_m(k,j-1,i),3,1), KIND = wp )

          v_comp(k)      = v(k,j,i) + v(k,j,i-1) - gv
          flux_s_u(k,tn) = v_comp(k) * (                                                           &
                           (  37.0_wp * ibit5 * adv_mom_5                                          &
                           +   7.0_wp * ibit4 * adv_mom_3                                          &
                           +            ibit3 * adv_mom_1 ) * ( u(k,j,i)    + u(k,j-1,i) )         &
                           - ( 8.0_wp * ibit5 * adv_mom_5                                          &
                           +            ibit4 * adv_mom_3 ) * ( u(k,j+1,i)  + u(k,j-2,i) )         &
                           + (          ibit5 * adv_mom_5 ) *  ( u(k,j+2,i) + u(k,j-3,i) )         &
                                       )

          diss_s_u(k,tn) = -ABS ( v_comp(k) ) * (                                                  &
                           (  10.0_wp * ibit5 * adv_mom_5                                          &
                           +   3.0_wp * ibit4 * adv_mom_3                                          &
                           +            ibit3 * adv_mom_1 ) * ( u(k,j,i)   - u(k,j-1,i) )          &
                           - ( 5.0_wp * ibit5 * adv_mom_5                                          &
                           +            ibit4 * adv_mom_3 ) * ( u(k,j+1,i) - u(k,j-2,i) )          &
                           + (          ibit5 * adv_mom_5 ) * ( u(k,j+2,i) - u(k,j-3,i) )          &
                                                )

       ENDDO

       DO  k = nzb_max_l+1, nzt

          v_comp(k)      = v(k,j,i) + v(k,j,i-1) - gv
          flux_s_u(k,tn) = v_comp(k) * (                                                           &
                             37.0_wp * ( u(k,j,i)   + u(k,j-1,i) )                                 &
                           -  8.0_wp * ( u(k,j+1,i) + u(k,j-2,i) )                                 &
                           +           ( u(k,j+2,i) + u(k,j-3,i) ) ) * adv_mom_5
          diss_s_u(k,tn) = -ABS(v_comp(k)) * (                                                     &
                                    10.0_wp * ( u(k,j,i)   - u(k,j-1,i) )                          &
                                  -  5.0_wp * ( u(k,j+1,i) - u(k,j-2,i) )                          &
                                  +           ( u(k,j+2,i) - u(k,j-3,i) ) ) * adv_mom_5

       ENDDO

    ENDIF
!
!-- Compute leftside fluxes for the respective boundary of PE
    IF ( i == i_omp  .OR.  i == nxlu )  THEN

       DO  k = nzb+1, nzb_max_l

          ibit2 = REAL( IBITS(advc_flags_m(k,j,i-1),2,1), KIND = wp )
          ibit1 = REAL( IBITS(advc_flags_m(k,j,i-1),1,1), KIND = wp )
          ibit0 = REAL( IBITS(advc_flags_m(k,j,i-1),0,1), KIND = wp )

          u_comp_l         = u(k,j,i) + u(k,j,i-1) - gu
          flux_l_u(k,j,tn) = u_comp_l * (                                                          &
                             (  37.0_wp * ibit2 * adv_mom_5                                        &
                             +   7.0_wp * ibit1 * adv_mom_3                                        &
                             +            ibit0 * adv_mom_1 ) * ( u(k,j,i)   + u(k,j,i-1) )        &
                             - ( 8.0_wp * ibit2 * adv_mom_5                                        &
                             +            ibit1 * adv_mom_3 ) * ( u(k,j,i+1) + u(k,j,i-2) )        &
                             + (          ibit2 * adv_mom_5 ) * ( u(k,j,i+2) + u(k,j,i-3) )        &
                                        )

          diss_l_u(k,j,tn) = -ABS( u_comp_l ) * (                                                  &
                             (  10.0_wp * ibit2 * adv_mom_5                                        &
                             +   3.0_wp * ibit1 * adv_mom_3                                        &
                             +           ibit0  * adv_mom_1 ) * ( u(k,j,i)   - u(k,j,i-1) )        &
                             - ( 5.0_wp * ibit2 * adv_mom_5                                        &
                             +            ibit1 * adv_mom_3 ) * ( u(k,j,i+1) - u(k,j,i-2) )        &
                             + (          ibit2 * adv_mom_5 ) * ( u(k,j,i+2) - u(k,j,i-3) )        &
                                                )

       ENDDO

       DO  k = nzb_max_l+1, nzt

          u_comp_l         = u(k,j,i) + u(k,j,i-1) - gu
          flux_l_u(k,j,tn) = u_comp_l * (                                                          &
                               37.0_wp * ( u(k,j,i)   + u(k,j,i-1) )                               &
                             -  8.0_wp * ( u(k,j,i+1) + u(k,j,i-2) )                               &
                             +           ( u(k,j,i+2) + u(k,j,i-3) ) ) * adv_mom_5
          diss_l_u(k,j,tn) = -ABS(u_comp_l) * (                                                    &
                               10.0_wp * ( u(k,j,i)   - u(k,j,i-1) )                               &
                             -  5.0_wp * ( u(k,j,i+1) - u(k,j,i-2) )                               &
                             +           ( u(k,j,i+2) - u(k,j,i-3) ) ) * adv_mom_5

       ENDDO

    ENDIF
!
!-- Now compute the fluxes tendency terms for the horizontal and vertical parts.
    DO  k = nzb+1, nzb_max_l

       ibit2 = REAL( IBITS(advc_flags_m(k,j,i),2,1), KIND = wp )
       ibit1 = REAL( IBITS(advc_flags_m(k,j,i),1,1), KIND = wp )
       ibit0 = REAL( IBITS(advc_flags_m(k,j,i),0,1), KIND = wp )

       u_comp(k) = u(k,j,i+1) + u(k,j,i)
       flux_r(k) = ( u_comp(k) - gu ) * (                                                          &
                   (  37.0_wp * ibit2 * adv_mom_5                                                  &
                   +   7.0_wp * ibit1 * adv_mom_3                                                  &
                   +            ibit0 * adv_mom_1 ) * ( u(k,j,i+1) + u(k,j,i)   )                  &
                   - ( 8.0_wp * ibit2 * adv_mom_5                                                  &
                   +            ibit1 * adv_mom_3 ) * ( u(k,j,i+2) + u(k,j,i-1) )                  &
                   + (          ibit2 * adv_mom_5 ) * ( u(k,j,i+3) + u(k,j,i-2) )                  &
                                        )

       diss_r(k) = -ABS( u_comp(k) - gu ) * (                                                      &
                   (  10.0_wp * ibit2 * adv_mom_5                                                  &
                   +   3.0_wp * ibit1 * adv_mom_3                                                  &
                   +            ibit0 * adv_mom_1 ) * ( u(k,j,i+1) - u(k,j,i)   )                  &
                   - ( 5.0_wp * ibit2 * adv_mom_5                                                  &
                   +            ibit1 * adv_mom_3 ) * ( u(k,j,i+2) - u(k,j,i-1) )                  &
                   + (          ibit2 * adv_mom_5 ) * ( u(k,j,i+3) - u(k,j,i-2) )                  &
                                            )

       ibit5 = REAL( IBITS(advc_flags_m(k,j,i),5,1), KIND = wp )
       ibit4 = REAL( IBITS(advc_flags_m(k,j,i),4,1), KIND = wp )
       ibit3 = REAL( IBITS(advc_flags_m(k,j,i),3,1), KIND = wp )

       v_comp(k) = v(k,j+1,i) + v(k,j+1,i-1) - gv
       flux_n(k) = v_comp(k) * (                                                                   &
                   (  37.0_wp * ibit5 * adv_mom_5                                                  &
                   +   7.0_wp * ibit4 * adv_mom_3                                                  &
                   +            ibit3 * adv_mom_1 ) * ( u(k,j+1,i) + u(k,j,i)   )                  &
                   - ( 8.0_wp * ibit5 * adv_mom_5                                                  &
                   +            ibit4 * adv_mom_3 ) * ( u(k,j+2,i) + u(k,j-1,i) )                  &
                   + (          ibit5 * adv_mom_5 ) * ( u(k,j+3,i) + u(k,j-2,i) )                  &
                               )

       diss_n(k) = -ABS ( v_comp(k) ) * (                                                          &
                   (  10.0_wp * ibit5 * adv_mom_5                                                  &
                   +   3.0_wp * ibit4 * adv_mom_3                                                  &
                   +    ibit3 * adv_mom_1 ) * ( u(k,j+1,i) - u(k,j,i)   )                          &
                   - ( 5.0_wp * ibit5 * adv_mom_5                                                  &
                   +    ibit4 * adv_mom_3 ) * ( u(k,j+2,i) - u(k,j-1,i) )                          &
                   + (  ibit5 * adv_mom_5 ) * ( u(k,j+3,i) - u(k,j-2,i) )                          &
                                        )
    ENDDO

    DO  k = nzb_max_l+1, nzt

       u_comp(k) = u(k,j,i+1) + u(k,j,i)
       flux_r(k) = ( u_comp(k) - gu ) * (                                                          &
                      37.0_wp * ( u(k,j,i+1) + u(k,j,i)   )                                        &
                    -  8.0_wp * ( u(k,j,i+2) + u(k,j,i-1) )                                        &
                    +           ( u(k,j,i+3) + u(k,j,i-2) ) ) * adv_mom_5
       diss_r(k) = -ABS( u_comp(k) - gu ) * (                                                      &
                      10.0_wp * ( u(k,j,i+1) - u(k,j,i)   )                                        &
                    -  5.0_wp * ( u(k,j,i+2) - u(k,j,i-1) )                                        &
                    +           ( u(k,j,i+3) - u(k,j,i-2) ) ) * adv_mom_5

       v_comp(k) = v(k,j+1,i) + v(k,j+1,i-1) - gv
       flux_n(k) = v_comp(k) * (                                                                   &
                      37.0_wp * ( u(k,j+1,i) + u(k,j,i)   )                                        &
                    -  8.0_wp * ( u(k,j+2,i) + u(k,j-1,i) )                                        &
                    +           ( u(k,j+3,i) + u(k,j-2,i) ) ) * adv_mom_5
       diss_n(k) = -ABS( v_comp(k) ) * (                                                           &
                      10.0_wp * ( u(k,j+1,i) - u(k,j,i)   )                                        &
                    -  5.0_wp * ( u(k,j+2,i) - u(k,j-1,i) )                                        &
                    +           ( u(k,j+3,i) - u(k,j-2,i) ) ) * adv_mom_5

    ENDDO
!
!-- Now, compute vertical fluxes. Split loop into a part treating the lowest grid points with
!-- indirect indexing, a main loop without indirect indexing, and a loop for the uppermost grid
!-- points with indirect indexing. This allows better vectorization for the main loop.
!-- First, compute the flux at model surface, which need has to be calculated explicitly for the
!-- tendency at the first w-level. For topography wall this is done implicitely by advc_flags_m.
    k = nzb
    ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )
    w_comp(k) = w(k,j,i) + w(k,j,i-1)
    flux_t(nzb) = w_comp(k)         * rho_air_zw(k) * ( ibit6 * adv_mom_1 )                        &
                                                    * ( u(k+1,j,i) + u(k,j,i) )
    diss_t(nzb) = -ABS( w_comp(k) ) * rho_air_zw(k) * ( ibit6 * adv_mom_1 )                        &
                                                    * ( u(k+1,j,i) - u(k,j,i) )

    DO  k = nzb+1, nzb+1
!
!--    k index has to be modified near bottom and top, else array subscripts will be exceeded.
       ibit8 = REAL( IBITS(advc_flags_m(k,j,i),8,1), KIND = wp )
       ibit7 = REAL( IBITS(advc_flags_m(k,j,i),7,1), KIND = wp )
       ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )

       k_ppp = k + 3 * ibit8
       k_pp  = k + 2 * ( 1 - ibit6 )
       k_mm  = k - 2 * ibit8

       w_comp(k) = w(k,j,i) + w(k,j,i-1)
       flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                                   &
                   (  37.0_wp * ibit8 * adv_mom_5                                                  &
                   +   7.0_wp * ibit7 * adv_mom_3                                                  &
                   +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)   + u(k,j,i)    )               &
                   - ( 8.0_wp * ibit8 * adv_mom_5                                                  &
                   +            ibit7 * adv_mom_3 ) * ( u(k_pp,j,i)  + u(k-1,j,i)  )               &
                   + (          ibit8 * adv_mom_5 ) * ( u(k_ppp,j,i) + u(k_mm,j,i) )               &
                                               )

       diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                           &
                   (  10.0_wp * ibit8 * adv_mom_5                                                  &
                   +   3.0_wp * ibit7 * adv_mom_3                                                  &
                   +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)   - u(k,j,i)    )               &
                   - ( 5.0_wp * ibit8 * adv_mom_5                                                  &
                   +            ibit7 * adv_mom_3 ) * ( u(k_pp,j,i)  - u(k-1,j,i)  )               &
                   + (          ibit8 * adv_mom_5 ) * ( u(k_ppp,j,i) - u(k_mm,j,i) )               &
                                                       )
    ENDDO

    DO  k = nzb+2, nzt-2

       ibit8 = REAL( IBITS(advc_flags_m(k,j,i),8,1), KIND = wp )
       ibit7 = REAL( IBITS(advc_flags_m(k,j,i),7,1), KIND = wp )
       ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )

       w_comp(k) = w(k,j,i) + w(k,j,i-1)
       flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                                   &
                   (  37.0_wp * ibit8 * adv_mom_5                                                  &
                   +   7.0_wp * ibit7 * adv_mom_3                                                  &
                   +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)  + u(k,j,i)  )                  &
                   - ( 8.0_wp * ibit8 * adv_mom_5                                                  &
                   +            ibit7 * adv_mom_3 ) * ( u(k+2,j,i) + u(k-1,j,i) )                  &
                   + (          ibit8 * adv_mom_5 ) * ( u(k+3,j,i) + u(k-2,j,i) )                  &
                                               )

       diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                           &
                   (  10.0_wp * ibit8 * adv_mom_5                                                  &
                   +   3.0_wp * ibit7 * adv_mom_3                                                  &
                   +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)  - u(k,j,i)   )                 &
                   - ( 5.0_wp * ibit8 * adv_mom_5                                                  &
                   +            ibit7 * adv_mom_3 ) * ( u(k+2,j,i)  - u(k-1,j,i) )                 &
                   + (          ibit8 * adv_mom_5 ) * ( u(k+3,j,i) - u(k-2,j,i)  )                 &
                                                       )
    ENDDO

    DO  k = nzt-1, nzt-symmetry_flag
!
!--    k index has to be modified near bottom and top, else array subscripts will be exceeded.
       ibit8 = REAL( IBITS(advc_flags_m(k,j,i),8,1), KIND = wp )
       ibit7 = REAL( IBITS(advc_flags_m(k,j,i),7,1), KIND = wp )
       ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )
!
!--    k index needs to be modified near bottom and top, else array subscripts will be exceeded.
!--    (k_ppp = k+3, k_mm = k) if 5th order scheme is employed, else (k_ppp = k, k_mm = k).
!--    Special treatment is required for k_pp. This must be limited at k = nzt, where either
!--    the first order scheme is employed (ibit6), or all bit flags are zero (in case there is a
!--    rigid lid defined at or near the domain top).
       k_ppp = k + 3 * ibit8
       k_pp  = k + 2 * ( 1 - ibit6 ) * ( ibit6 + ibit7 + ibit8 )
       k_mm  = k - 2 * ibit8

       w_comp(k) = w(k,j,i) + w(k,j,i-1)
       flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                                   &
                   (  37.0_wp * ibit8 * adv_mom_5                                                  &
                   +   7.0_wp * ibit7 * adv_mom_3                                                  &
                   +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)   + u(k,j,i)    )               &
                   - ( 8.0_wp * ibit8 * adv_mom_5                                                  &
                   +            ibit7 * adv_mom_3 ) * ( u(k_pp,j,i)  + u(k-1,j,i)  )               &
                   + (          ibit8 * adv_mom_5 ) * ( u(k_ppp,j,i) + u(k_mm,j,i) )               &
                                               )

       diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                           &
                   (  10.0_wp * ibit8 * adv_mom_5                                                  &
                   +   3.0_wp * ibit7 * adv_mom_3                                                  &
                   +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)   - u(k,j,i)    )               &
                   - ( 5.0_wp * ibit8 * adv_mom_5                                                  &
                   +            ibit7 * adv_mom_3 ) * ( u(k_pp,j,i)  - u(k-1,j,i)  )               &
                   + (          ibit8 * adv_mom_5 ) * ( u(k_ppp,j,i) - u(k_mm,j,i) )               &
                                                       )
    ENDDO

!
!-- Set resolved/turbulent flux at model top to zero (w-level). In case that a symmetric behavior
!-- between bottom and top shall be guaranteed (closed channel flow), the flux at nzt is also set to
!-- zero.
    IF ( symmetry_flag == 1 )  THEN
       flux_t(nzt) = 0.0_wp
       diss_t(nzt) = 0.0_wp
       w_comp(nzt) = 0.0_wp
    ENDIF
    flux_t(nzt+1) = 0.0_wp
    diss_t(nzt+1) = 0.0_wp
    w_comp(nzt+1) = 0.0_wp

    DO  k = nzb+1, nzb_max_l

       flux_d    = flux_t(k-1)
       diss_d    = diss_t(k-1)

       ibit2 = REAL( IBITS(advc_flags_m(k,j,i),2,1), KIND = wp )
       ibit1 = REAL( IBITS(advc_flags_m(k,j,i),1,1), KIND = wp )
       ibit0 = REAL( IBITS(advc_flags_m(k,j,i),0,1), KIND = wp )

       ibit5 = REAL( IBITS(advc_flags_m(k,j,i),5,1), KIND = wp )
       ibit4 = REAL( IBITS(advc_flags_m(k,j,i),4,1), KIND = wp )
       ibit3 = REAL( IBITS(advc_flags_m(k,j,i),3,1), KIND = wp )

       ibit8 = REAL( IBITS(advc_flags_m(k,j,i),8,1), KIND = wp )
       ibit7 = REAL( IBITS(advc_flags_m(k,j,i),7,1), KIND = wp )
       ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )
!
!--    Calculate the divergence of the velocity field. A respective correction is needed to overcome
!--    numerical instabilities introduced by an insufficient reduction of divergences near
!--    topography.
       div = ( ( u_comp(k)       * ( ibit0 + ibit1 + ibit2 )                                       &
             - ( u(k,j,i)   + u(k,j,i-1)   )                                                       &
                                 * (                                                               &
                  REAL( IBITS(advc_flags_m(k,j,i-1),0,1), KIND = wp )                              &
                + REAL( IBITS(advc_flags_m(k,j,i-1),1,1), KIND = wp )                              &
                + REAL( IBITS(advc_flags_m(k,j,i-1),2,1), KIND = wp )                              &
                                   )                                                               &
               ) * ddx                                                                             &
             + ( ( v_comp(k) + gv ) * ( ibit3 + ibit4 + ibit5 )                                    &
               - ( v(k,j,i)   + v(k,j,i-1 )  )                                                     &
                                 * (                                                               &
                  REAL( IBITS(advc_flags_m(k,j-1,i),3,1), KIND = wp )                              &
                + REAL( IBITS(advc_flags_m(k,j-1,i),4,1), KIND = wp )                              &
                + REAL( IBITS(advc_flags_m(k,j-1,i),5,1), KIND = wp )                              &
                                   )                                                               &
               ) * ddy                                                                             &
             + ( w_comp(k)   * rho_air_zw(k) * ( ibit6 + ibit7 + ibit8 )                           &
             -   w_comp(k-1) * rho_air_zw(k-1)                                                     &
                                 * (                                                               &
                  REAL( IBITS(advc_flags_m(k-1,j,i),6,1), KIND = wp )                              &
                + REAL( IBITS(advc_flags_m(k-1,j,i),7,1), KIND = wp )                              &
                + REAL( IBITS(advc_flags_m(k-1,j,i),8,1), KIND = wp )                              &
                                   )                                                               &
               ) * drho_air(k) * ddzw(k)                                                           &
             ) * 0.5_wp

       tend(k,j,i) = tend(k,j,i) - (                                                               &
                         ( flux_r(k) + diss_r(k)                                                   &
                       -   flux_l_u(k,j,tn) - diss_l_u(k,j,tn) ) * ddx                             &
                       + ( flux_n(k) + diss_n(k)                                                   &
                       -   flux_s_u(k,tn) - diss_s_u(k,tn)     ) * ddy                             &
                       + ( ( flux_t(k) + diss_t(k) )                                               &
                       -   ( flux_d    + diss_d    )                                               &
                         )                         * drho_air(k) * ddzw(k)                         &
                                   ) + div * u(k,j,i)

       flux_l_u(k,j,tn) = flux_r(k)
       diss_l_u(k,j,tn) = diss_r(k)
       flux_s_u(k,tn)   = flux_n(k)
       diss_s_u(k,tn)   = diss_n(k)
!
!--    Statistical Evaluation of u'u'. The factor has to be applied for right evaluation when
!--    gallilei_trans = .T. .
       sums_us2_ws_l(k,tn) = sums_us2_ws_l(k,tn) +                                                 &
                             ( flux_r(k)                                                           &
                               * ( u_comp(k) - 2.0_wp * hom(k,1,1,0)                   )           &
                               / ( u_comp(k) - gu + SIGN( 1.0E-20_wp, u_comp(k) - gu ) )           &
                             + diss_r(k)                                                           &
                               *   ABS( u_comp(k) - 2.0_wp * hom(k,1,1,0)              )           &
                               / ( ABS( u_comp(k) - gu ) + 1.0E-20_wp                  )           &
                             ) *   weight_substep(intermediate_timestep_count)
!
!--    Statistical Evaluation of w'u'.
       sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn) +                                               &
                              ( flux_t(k)                                                          &
                                * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                   )          &
                                / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )           )          &
                              + diss_t(k)                                                          &
                                *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)              )          &
                                / ( ABS( w_comp(k) ) + 1.0E-20_wp                       )          &
                              ) *   weight_substep(intermediate_timestep_count)
    ENDDO

    DO  k = nzb_max_l+1, nzt

       flux_d    = flux_t(k-1)
       diss_d    = diss_t(k-1)
!
!--    Calculate the divergence of the velocity field. A respective correction is needed to overcome
!--    numerical instabilities introduced by an insufficient reduction of divergences near
!--    topography.
       div = ( ( u_comp(k)       - ( u(k,j,i)   + u(k,j,i-1) ) ) * ddx                             &
             + ( v_comp(k) + gv  - ( v(k,j,i)   + v(k,j,i-1) ) ) * ddy                             &
             + ( w_comp(k)   * rho_air_zw(k)                                                       &
             -   w_comp(k-1) * rho_air_zw(k-1)                                                     &
               )                                   * drho_air(k) * ddzw(k)                         &
             ) * 0.5_wp

       tend(k,j,i) = tend(k,j,i) - (                                                               &
                         ( flux_r(k) + diss_r(k)                                                   &
                       -   flux_l_u(k,j,tn) - diss_l_u(k,j,tn) ) * ddx                             &
                       + ( flux_n(k) + diss_n(k)                                                   &
                       -   flux_s_u(k,tn) - diss_s_u(k,tn)     ) * ddy                             &
                       + ( ( flux_t(k) + diss_t(k) )                                               &
                       -   ( flux_d    + diss_d    )                                               &
                                                 ) * drho_air(k) * ddzw(k)                         &
                                   ) + div * u(k,j,i)

       flux_l_u(k,j,tn) = flux_r(k)
       diss_l_u(k,j,tn) = diss_r(k)
       flux_s_u(k,tn)   = flux_n(k)
       diss_s_u(k,tn)   = diss_n(k)
!
!--    Statistical Evaluation of u'u'. The factor has to be applied for right evaluation when
!--    gallilei_trans = .T. .
       sums_us2_ws_l(k,tn) = sums_us2_ws_l(k,tn) +                                                 &
                             ( flux_r(k)                                                           &
                               * ( u_comp(k) - 2.0_wp * hom(k,1,1,0)                   )           &
                               / ( u_comp(k) - gu + SIGN( 1.0E-20_wp, u_comp(k) - gu ) )           &
                             + diss_r(k)                                                           &
                               *   ABS( u_comp(k) - 2.0_wp * hom(k,1,1,0)              )           &
                               / ( ABS( u_comp(k) - gu ) + 1.0E-20_wp                  )           &
                             ) *   weight_substep(intermediate_timestep_count)
!
!--    Statistical Evaluation of w'u'.
       sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn) +                                               &
                              ( flux_t(k)                                                          &
                                * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                   )          &
                                / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )           )          &
                              + diss_t(k)                                                          &
                                *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)              )          &
                                / ( ABS( w_comp(k) ) + 1.0E-20_wp                       )          &
                              ) *   weight_substep(intermediate_timestep_count)
    ENDDO



 END SUBROUTINE advec_u_ws_ij



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of v-component - Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_v_ws_ij( i, j, i_omp, tn )


    INTEGER(iwp)  ::  i         !< grid index along x-direction
    INTEGER(iwp)  ::  i_omp     !< leftmost index on subdomain, or in case of OpenMP, on thread
    INTEGER(iwp)  ::  j         !< grid index along y-direction
    INTEGER(iwp)  ::  k         !< grid index along z-direction
    INTEGER(iwp)  ::  k_mm      !< k-2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp)  ::  k_pp      !< k+2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp)  ::  k_ppp     !< k+3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp)  ::  nzb_max_l !< index indicating upper bound for order degradation of horizontal advection terms
    INTEGER(iwp)  ::  tn        !< number of OpenMP thread

    REAL(wp)      ::  diss_d   !< artificial dissipation term at grid box bottom
    REAL(wp)      ::  div      !< divergence on v-grid
    REAL(wp)      ::  flux_d   !< 6th-order flux at grid box bottom
    REAL(wp)      ::  gu       !< Galilei-transformation velocity along x
    REAL(wp)      ::  gv       !< Galilei-transformation velocity along y
    REAL(wp)      ::  ibit9    !< flag indicating 1st-order scheme along x-direction
    REAL(wp)      ::  ibit10   !< flag indicating 3rd-order scheme along x-direction
    REAL(wp)      ::  ibit11   !< flag indicating 5th-order scheme along x-direction
    REAL(wp)      ::  ibit12   !< flag indicating 1st-order scheme along y-direction
    REAL(wp)      ::  ibit13   !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)      ::  ibit14   !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)      ::  ibit15   !< flag indicating 1st-order scheme along z-direction
    REAL(wp)      ::  ibit16   !< flag indicating 3rd-order scheme along z-direction
    REAL(wp)      ::  ibit17   !< flag indicating 3rd-order scheme along z-direction
    REAL(wp)      ::  v_comp_l !< advection velocity along y on leftmost grid point on subdomain

    REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_n !< discretized artificial dissipation at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_r !< discretized artificial dissipation at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_t !< discretized artificial dissipation at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_n !< discretized 6th-order flux at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_r !< discretized 6th-order flux at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_t !< discretized 6th-order flux at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  u_comp !< advection velocity along x
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  v_comp !< advection velocity along y
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  w_comp !< advection velocity along z
!
!-- Used local modified copy of nzb_max (used to degrade order of discretization) at non-cyclic
!-- boundaries. Modify only at relevant points instead of the entire subdomain. This should lead to
!-- better load balance between boundary and non-boundary PEs.
    IF( ( bc_dirichlet_l  .OR.  bc_radiation_l )  .AND.  i <= nxl + 2  .OR.                        &
        ( bc_dirichlet_r  .OR.  bc_radiation_r )  .AND.  i >= nxr - 2  .OR.                        &
        ( bc_dirichlet_s  .OR.  bc_radiation_s )  .AND.  j <= nys + 2  .OR.                        &
        ( bc_dirichlet_n  .OR.  bc_radiation_n )  .AND.  j >= nyn - 2 )  THEN
       nzb_max_l = nzt
    ELSE
       nzb_max_l = nzb_max
    END IF

    gu = 2.0_wp * u_gtrans
    gv = 2.0_wp * v_gtrans

!
!-- Compute leftside fluxes for the respective boundary.
    IF ( i == i_omp )  THEN

       DO  k = nzb+1, nzb_max_l

          ibit11 = REAL( IBITS(advc_flags_m(k,j,i-1),11,1), KIND = wp )
          ibit10 = REAL( IBITS(advc_flags_m(k,j,i-1),10,1), KIND = wp )
          ibit9  = REAL( IBITS(advc_flags_m(k,j,i-1),9,1),  KIND = wp )

          u_comp(k)        = u(k,j-1,i) + u(k,j,i) - gu
          flux_l_v(k,j,tn) = u_comp(k) * (                                                         &
                             (  37.0_wp * ibit11 * adv_mom_5                                       &
                             +   7.0_wp * ibit10 * adv_mom_3                                       &
                             +            ibit9  * adv_mom_1 ) * ( v(k,j,i)   + v(k,j,i-1) )       &
                             - ( 8.0_wp * ibit11 * adv_mom_5                                       &
                             +            ibit10 * adv_mom_3 ) * ( v(k,j,i+1) + v(k,j,i-2) )       &
                             + (          ibit11 * adv_mom_5 ) * ( v(k,j,i+2) + v(k,j,i-3) )       &
                                         )

          diss_l_v(k,j,tn) = -ABS( u_comp(k) ) * (                                                 &
                             (  10.0_wp * ibit11 * adv_mom_5                                       &
                             +   3.0_wp * ibit10 * adv_mom_3                                       &
                             +            ibit9  * adv_mom_1 ) * ( v(k,j,i)   - v(k,j,i-1) )       &
                             - ( 5.0_wp * ibit11 * adv_mom_5                                       &
                             +            ibit10 * adv_mom_3 ) * ( v(k,j,i+1) - v(k,j,i-2) )       &
                             + (          ibit11 * adv_mom_5 ) * ( v(k,j,i+2) - v(k,j,i-3) )       &
                                                 )

       ENDDO

       DO  k = nzb_max_l+1, nzt

          u_comp(k)        = u(k,j-1,i) + u(k,j,i) - gu
          flux_l_v(k,j,tn) = u_comp(k) * (                                                         &
                               37.0_wp * ( v(k,j,i)   + v(k,j,i-1) )                               &
                             -  8.0_wp * ( v(k,j,i+1) + v(k,j,i-2) )                               &
                             +           ( v(k,j,i+2) + v(k,j,i-3) ) ) * adv_mom_5
          diss_l_v(k,j,tn) = -ABS( u_comp(k) ) * (                                                 &
                               10.0_wp * ( v(k,j,i)   - v(k,j,i-1) )                               &
                             -  5.0_wp * ( v(k,j,i+1) - v(k,j,i-2) )                               &
                             +           ( v(k,j,i+2) - v(k,j,i-3) ) ) * adv_mom_5

       ENDDO

    ENDIF
!
!-- Compute southside fluxes for the respective boundary.
    IF ( j == nysv )  THEN

       DO  k = nzb+1, nzb_max_l

          ibit14 = REAL( IBITS(advc_flags_m(k,j-1,i),14,1), KIND = wp )
          ibit13 = REAL( IBITS(advc_flags_m(k,j-1,i),13,1), KIND = wp )
          ibit12 = REAL( IBITS(advc_flags_m(k,j-1,i),12,1), KIND = wp )

          v_comp_l       = v(k,j,i) + v(k,j-1,i) - gv
          flux_s_v(k,tn) = v_comp_l * (                                                            &
                         (  37.0_wp * ibit14 * adv_mom_5                                           &
                         +   7.0_wp * ibit13 * adv_mom_3                                           &
                         +            ibit12 * adv_mom_1 ) * ( v(k,j,i)   + v(k,j-1,i) )           &
                         - ( 8.0_wp * ibit14 * adv_mom_5                                           &
                         +            ibit13 * adv_mom_3 ) * ( v(k,j+1,i) + v(k,j-2,i) )           &
                         + (          ibit14 * adv_mom_5 ) * ( v(k,j+2,i) + v(k,j-3,i) )           &
                                      )

          diss_s_v(k,tn) = -ABS( v_comp_l ) * (                                                    &
                           (  10.0_wp * ibit14 * adv_mom_5                                         &
                           +   3.0_wp * ibit13 * adv_mom_3                                         &
                           +            ibit12 * adv_mom_1 ) * ( v(k,j,i)   - v(k,j-1,i) )         &
                           - ( 5.0_wp * ibit14 * adv_mom_5                                         &
                           +            ibit13 * adv_mom_3 ) * ( v(k,j+1,i) - v(k,j-2,i) )         &
                           + (          ibit14 * adv_mom_5 ) * ( v(k,j+2,i) - v(k,j-3,i) )         &
                                              )

       ENDDO

       DO  k = nzb_max_l+1, nzt

          v_comp_l       = v(k,j,i) + v(k,j-1,i) - gv
          flux_s_v(k,tn) = v_comp_l * (                                                            &
                        37.0_wp * ( v(k,j,i)   + v(k,j-1,i)   )                                    &
                      -  8.0_wp * ( v(k,j+1,i) + v(k,j-2,i) )                                      &
                      +           ( v(k,j+2,i) + v(k,j-3,i) ) ) * adv_mom_5
          diss_s_v(k,tn) = -ABS( v_comp_l ) * (                                                    &
                        10.0_wp * ( v(k,j,i)   - v(k,j-1,i)   )                                    &
                      -  5.0_wp * ( v(k,j+1,i) - v(k,j-2,i) )                                      &
                      +           ( v(k,j+2,i) - v(k,j-3,i) ) ) * adv_mom_5

       ENDDO

    ENDIF
!
!--    Now compute the fluxes and tendency terms for the horizontal and
!--    verical parts.
    DO  k = nzb+1, nzb_max_l

       ibit11 = REAL( IBITS(advc_flags_m(k,j,i),11,1), KIND = wp )
       ibit10 = REAL( IBITS(advc_flags_m(k,j,i),10,1), KIND = wp )
       ibit9  = REAL( IBITS(advc_flags_m(k,j,i),9,1),  KIND = wp )

       u_comp(k) = u(k,j-1,i+1) + u(k,j,i+1) - gu
       flux_r(k) = u_comp(k) * (                                                                   &
                  (  37.0_wp * ibit11 * adv_mom_5                                                  &
                  +   7.0_wp * ibit10 * adv_mom_3                                                  &
                  +            ibit9  * adv_mom_1 ) * ( v(k,j,i+1) + v(k,j,i)   )                  &
                  - ( 8.0_wp * ibit11 * adv_mom_5                                                  &
                  +            ibit10 * adv_mom_3 ) * ( v(k,j,i+2) + v(k,j,i-1) )                  &
                  + (          ibit11 * adv_mom_5 ) * ( v(k,j,i+3) + v(k,j,i-2) )                  &
                               )

       diss_r(k) = -ABS( u_comp(k) ) * (                                                           &
                  (  10.0_wp * ibit11 * adv_mom_5                                                  &
                  +   3.0_wp * ibit10 * adv_mom_3                                                  &
                  +            ibit9  * adv_mom_1 ) * ( v(k,j,i+1) - v(k,j,i)  )                   &
                  - ( 5.0_wp * ibit11 * adv_mom_5                                                  &
                  +            ibit10 * adv_mom_3 ) * ( v(k,j,i+2) - v(k,j,i-1) )                  &
                  + (          ibit11 * adv_mom_5 ) * ( v(k,j,i+3) - v(k,j,i-2) )                  &
                                       )

       ibit14 = REAL( IBITS(advc_flags_m(k,j,i),14,1), KIND = wp )
       ibit13 = REAL( IBITS(advc_flags_m(k,j,i),13,1), KIND = wp )
       ibit12 = REAL( IBITS(advc_flags_m(k,j,i),12,1), KIND = wp )


       v_comp(k) = v(k,j+1,i) + v(k,j,i)
       flux_n(k) = ( v_comp(k) - gv ) * (                                                          &
                   (  37.0_wp * ibit14 * adv_mom_5                                                 &
                   +   7.0_wp * ibit13 * adv_mom_3                                                 &
                   +            ibit12 * adv_mom_1 ) * ( v(k,j+1,i) + v(k,j,i)   )                 &
                   - ( 8.0_wp * ibit14 * adv_mom_5                                                 &
                   +            ibit13 * adv_mom_3 ) * ( v(k,j+2,i) + v(k,j-1,i) )                 &
                   + (          ibit14 * adv_mom_5 ) * ( v(k,j+3,i) + v(k,j-2,i) )                 &
                                        )

       diss_n(k) = -ABS( v_comp(k) - gv ) * (                                                      &
                   (  10.0_wp * ibit14 * adv_mom_5                                                 &
                   +   3.0_wp * ibit13 * adv_mom_3                                                 &
                   +            ibit12 * adv_mom_1 ) * ( v(k,j+1,i) - v(k,j,i)   )                 &
                   - ( 5.0_wp * ibit14 * adv_mom_5                                                 &
                   +            ibit13 * adv_mom_3 ) * ( v(k,j+2,i) - v(k,j-1,i) )                 &
                   + (          ibit14 * adv_mom_5 ) * ( v(k,j+3,i) - v(k,j-2,i) )                 &
                                            )
    ENDDO

    DO  k = nzb_max_l+1, nzt

       u_comp(k) = u(k,j-1,i+1) + u(k,j,i+1) - gu
       flux_r(k) = u_comp(k) * (                                                                   &
                     37.0_wp * ( v(k,j,i+1) + v(k,j,i)   )                                         &
                   -  8.0_wp * ( v(k,j,i+2) + v(k,j,i-1) )                                         &
                   +           ( v(k,j,i+3) + v(k,j,i-2) ) ) * adv_mom_5

       diss_r(k) = -ABS( u_comp(k) ) * (                                                           &
                     10.0_wp * ( v(k,j,i+1) - v(k,j,i) )                                           &
                   -  5.0_wp * ( v(k,j,i+2) - v(k,j,i-1) )                                         &
                   +           ( v(k,j,i+3) - v(k,j,i-2) ) ) * adv_mom_5


       v_comp(k) = v(k,j+1,i) + v(k,j,i)
       flux_n(k) = ( v_comp(k) - gv ) * (                                                          &
                     37.0_wp * ( v(k,j+1,i) + v(k,j,i)   )                                         &
                   -  8.0_wp * ( v(k,j+2,i) + v(k,j-1,i) )                                         &
                   +           ( v(k,j+3,i) + v(k,j-2,i) ) ) * adv_mom_5

       diss_n(k) = -ABS( v_comp(k) - gv ) * (                                                      &
                     10.0_wp * ( v(k,j+1,i) - v(k,j,i)   )                                         &
                   -  5.0_wp * ( v(k,j+2,i) - v(k,j-1,i) )                                         &
                   +           ( v(k,j+3,i) - v(k,j-2,i) ) ) * adv_mom_5
    ENDDO
!
!-- Now, compute vertical fluxes. Split loop into a part treating the lowest grid points with
!-- indirect indexing, a main loop without indirect indexing, and a loop for the uppermost grid
!-- points with indirect indexing. This allows better vectorization for the main loop.
!-- First, compute the flux at model surface, which need has to be calculated explicitly for the
!-- tendency at the first w-level. For topography wall this is done implicitely by advc_flags_m.
    k = nzb
    ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )
    w_comp(k) = w(k,j-1,i) + w(k,j,i)
    flux_t(nzb) = w_comp(k)         * rho_air_zw(k) * ( ibit15 * adv_mom_1 )                       &
                                                    * ( v(k+1,j,i) + v(k,j,i) )
    diss_t(nzb) = -ABS( w_comp(k) ) * rho_air_zw(k) * ( ibit15 * adv_mom_1 )                       &
                                                    * ( v(k+1,j,i) - v(k,j,i) )

    DO  k = nzb+1, nzb+1
!
!--    k index has to be modified near bottom and top, else array
!--    subscripts will be exceeded.
       ibit17 = REAL( IBITS(advc_flags_m(k,j,i),17,1), KIND = wp )
       ibit16 = REAL( IBITS(advc_flags_m(k,j,i),16,1), KIND = wp )
       ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )

       k_ppp = k + 3 * ibit17
       k_pp  = k + 2 * ( 1 - ibit15  )
       k_mm  = k - 2 * ibit17

       w_comp(k) = w(k,j-1,i) + w(k,j,i)
       flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                                   &
                   (  37.0_wp * ibit17 * adv_mom_5                                                 &
                   +   7.0_wp * ibit16 * adv_mom_3                                                 &
                   +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i)   + v(k,j,i)    )              &
                   - ( 8.0_wp * ibit17 * adv_mom_5                                                 &
                   +            ibit16 * adv_mom_3 ) * ( v(k_pp,j,i)  + v(k-1,j,i)  )              &
                   + (          ibit17 * adv_mom_5 ) * ( v(k_ppp,j,i) + v(k_mm,j,i) )              &
                                               )

       diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                           &
                   (  10.0_wp * ibit17 * adv_mom_5                                                 &
                   +   3.0_wp * ibit16 * adv_mom_3                                                 &
                   +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i)   - v(k,j,i)    )              &
                   - ( 5.0_wp * ibit17 * adv_mom_5                                                 &
                   +            ibit16 * adv_mom_3 ) * ( v(k_pp,j,i)  - v(k-1,j,i)  )              &
                   + (          ibit17 * adv_mom_5 ) * ( v(k_ppp,j,i) - v(k_mm,j,i) )              &
                                                       )
    ENDDO

    DO  k = nzb+2, nzt-2

       ibit17 = REAL( IBITS(advc_flags_m(k,j,i),17,1), KIND = wp )
       ibit16 = REAL( IBITS(advc_flags_m(k,j,i),16,1), KIND = wp )
       ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )

       w_comp(k) = w(k,j-1,i) + w(k,j,i)
       flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                                   &
                  (  37.0_wp * ibit17 * adv_mom_5                                                  &
                  +   7.0_wp * ibit16 * adv_mom_3                                                  &
                  +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i) + v(k,j,i)   )                  &
                  - ( 8.0_wp * ibit17 * adv_mom_5                                                  &
                  +            ibit16 * adv_mom_3 ) * ( v(k+2,j,i) + v(k-1,j,i) )                  &
                  + (          ibit17 * adv_mom_5 ) * ( v(k+3,j,i) + v(k-2,j,i) )                  &
                                               )

       diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                           &
                   (  10.0_wp * ibit17 * adv_mom_5                                                 &
                   +   3.0_wp * ibit16 * adv_mom_3                                                 &
                   +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i) - v(k,j,i)   )                 &
                   - ( 5.0_wp * ibit17 * adv_mom_5                                                 &
                   +            ibit16 * adv_mom_3 ) * ( v(k+2,j,i) - v(k-1,j,i) )                 &
                   + (          ibit17 * adv_mom_5 ) * ( v(k+3,j,i) - v(k-2,j,i) )                 &
                                                       )
    ENDDO

    DO  k = nzt-1, nzt-symmetry_flag
!
!--    k index has to be modified near bottom and top, else array subscripts will be exceeded.
       ibit17 = REAL( IBITS(advc_flags_m(k,j,i),17,1), KIND = wp )
       ibit16 = REAL( IBITS(advc_flags_m(k,j,i),16,1), KIND = wp )
       ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )
!
!--    k index needs to be modified near bottom and top, else array subscripts will be exceeded.
!--    (k_ppp = k+3, k_mm = k) if 5th order scheme is employed, else (k_ppp = k, k_mm = k).
!--    Special treatment is required for k_pp. This must be limited at k = nzt, where either
!--    the first order scheme is employed (ibit15), or all bit flags are zero (in case there is a
!--    rigid lid defined at or near the domain top).
       k_ppp = k + 3 * ibit17
       k_pp  = k + 2 * ( 1 - ibit15  ) * ( ibit15 + ibit16 + ibit17 )
       k_mm  = k - 2 * ibit17

       w_comp(k) = w(k,j-1,i) + w(k,j,i)
       flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                                   &
                   (  37.0_wp * ibit17 * adv_mom_5                                                 &
                   +   7.0_wp * ibit16 * adv_mom_3                                                 &
                   +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i)   + v(k,j,i)    )              &
                   - ( 8.0_wp * ibit17 * adv_mom_5                                                 &
                   +            ibit16 * adv_mom_3 ) * ( v(k_pp,j,i)  + v(k-1,j,i)  )              &
                   + (          ibit17 * adv_mom_5 ) * ( v(k_ppp,j,i) + v(k_mm,j,i) )              &
                                               )

       diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                           &
                   (  10.0_wp * ibit17 * adv_mom_5                                                 &
                   +   3.0_wp * ibit16 * adv_mom_3                                                 &
                   +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i)   - v(k,j,i)    )              &
                   - ( 5.0_wp * ibit17 * adv_mom_5                                                 &
                   +            ibit16 * adv_mom_3 ) * ( v(k_pp,j,i)  - v(k-1,j,i)  )              &
                   + (          ibit17 * adv_mom_5 ) * ( v(k_ppp,j,i) - v(k_mm,j,i) )              &
                                                       )
    ENDDO

!
!-- Set resolved/turbulent flux at model top to zero (w-level). In case that a symmetric behavior
!-- between bottom and top shall be guaranteed (closed channel flow), the flux at nzt is also set to
!-- zero.
    IF ( symmetry_flag == 1 ) THEN
       flux_t(nzt) = 0.0_wp
       diss_t(nzt) = 0.0_wp
       w_comp(nzt) = 0.0_wp
    ENDIF
    flux_t(nzt+1) = 0.0_wp
    diss_t(nzt+1) = 0.0_wp
    w_comp(nzt+1) = 0.0_wp

    DO  k = nzb+1, nzb_max_l

       flux_d    = flux_t(k-1)
       diss_d    = diss_t(k-1)

       ibit11 = REAL( IBITS(advc_flags_m(k,j,i),11,1), KIND = wp )
       ibit10 = REAL( IBITS(advc_flags_m(k,j,i),10,1), KIND = wp )
       ibit9  = REAL( IBITS(advc_flags_m(k,j,i),9,1),  KIND = wp )

       ibit14 = REAL( IBITS(advc_flags_m(k,j,i),14,1), KIND = wp )
       ibit13 = REAL( IBITS(advc_flags_m(k,j,i),13,1), KIND = wp )
       ibit12 = REAL( IBITS(advc_flags_m(k,j,i),12,1), KIND = wp )

       ibit17 = REAL( IBITS(advc_flags_m(k,j,i),17,1), KIND = wp )
       ibit16 = REAL( IBITS(advc_flags_m(k,j,i),16,1), KIND = wp )
       ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )
!
!--    Calculate the divergence of the velocity field. A respective correction is needed to overcome
!--    numerical instabilities introduced by an insufficient reduction of divergences near
!--    topography.
       div = ( ( ( u_comp(k)     + gu )                                                            &
                                    * ( ibit9 + ibit10 + ibit11 )                                  &
               - ( u(k,j-1,i) + u(k,j,i) )                                                         &
                                    * (                                                            &
                     REAL( IBITS(advc_flags_m(k,j,i-1),9,1),  KIND = wp )                          &
                   + REAL( IBITS(advc_flags_m(k,j,i-1),10,1), KIND = wp )                          &
                   + REAL( IBITS(advc_flags_m(k,j,i-1),11,1), KIND = wp )                          &
                                      )                                                            &
               ) * ddx                                                                             &
             + ( v_comp(k)                                                                         &
                                    * ( ibit12 + ibit13 + ibit14 )                                 &
             - ( v(k,j,i)     + v(k,j-1,i) )                                                       &
                                    * (                                                            &
                     REAL( IBITS(advc_flags_m(k,j-1,i),12,1), KIND = wp )                          &
                   + REAL( IBITS(advc_flags_m(k,j-1,i),13,1), KIND = wp )                          &
                   + REAL( IBITS(advc_flags_m(k,j-1,i),14,1), KIND = wp )                          &
                                      )                                                            &
               ) * ddy                                                                             &
             + ( w_comp(k)   * rho_air_zw(k) * ( ibit15 + ibit16 + ibit17 )                        &
             -   w_comp(k-1) * rho_air_zw(k-1)                                                     &
                                    * (                                                            &
                     REAL( IBITS(advc_flags_m(k-1,j,i),15,1), KIND = wp )                          &
                   + REAL( IBITS(advc_flags_m(k-1,j,i),16,1), KIND = wp )                          &
                   + REAL( IBITS(advc_flags_m(k-1,j,i),17,1), KIND = wp )                          &
                                      )                                                            &
                ) * drho_air(k) * ddzw(k)                                                          &
             ) * 0.5_wp

       tend(k,j,i) = tend(k,j,i) - (                                                               &
                      ( flux_r(k) + diss_r(k)                                                      &
                    -   flux_l_v(k,j,tn) - diss_l_v(k,j,tn)   ) * ddx                              &
                    + ( flux_n(k) + diss_n(k)                                                      &
                    -   flux_s_v(k,tn) - diss_s_v(k,tn)       ) * ddy                              &
                    + ( ( flux_t(k) + diss_t(k) )                                                  &
                    -   ( flux_d    + diss_d    )                                                  &
                                                ) * drho_air(k) * ddzw(k)                          &
                                   ) + v(k,j,i) * div

       flux_l_v(k,j,tn) = flux_r(k)
       diss_l_v(k,j,tn) = diss_r(k)
       flux_s_v(k,tn)   = flux_n(k)
       diss_s_v(k,tn)   = diss_n(k)
!
!--    Statistical Evaluation of v'v'. The factor has to be applied for right evaluation when
!--    gallilei_trans = .T. .
       sums_vs2_ws_l(k,tn) = sums_vs2_ws_l(k,tn) +                                                 &
                             ( flux_n(k)                                                           &
                               * ( v_comp(k) - 2.0_wp * hom(k,1,2,0)                   )           &
                               / ( v_comp(k) - gv + SIGN( 1.0E-20_wp, v_comp(k) - gv ) )           &
                             + diss_n(k)                                                           &
                               *   ABS( v_comp(k) - 2.0_wp * hom(k,1,2,0)              )           &
                               / ( ABS( v_comp(k) - gv ) + 1.0E-20_wp                  )           &
                             ) *   weight_substep(intermediate_timestep_count)
!
!--    Statistical Evaluation of w'u'.
       sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn) +                                               &
                              ( flux_t(k)                                                          &
                                * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                   )          &
                                / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )           )          &
                              + diss_t(k)                                                          &
                                *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)              )          &
                                / ( ABS( w_comp(k) ) + 1.0E-20_wp                       )          &
                              ) *   weight_substep(intermediate_timestep_count)

    ENDDO

    DO  k = nzb_max_l+1, nzt

       flux_d    = flux_t(k-1)
       diss_d    = diss_t(k-1)
!
!--    Calculate the divergence of the velocity field. A respective correction is needed to overcome
!--    numerical instabilities introduced by an insufficient reduction of divergences near
!--    topography.
       div = ( ( u_comp(k) + gu - ( u(k,j-1,i) + u(k,j,i)   ) ) * ddx                              &
             + ( v_comp(k)      - ( v(k,j,i)   + v(k,j-1,i) ) ) * ddy                              &
             + ( w_comp(k)   * rho_air_zw(k)                                                       &
             -   w_comp(k-1) * rho_air_zw(k-1)                                                     &
               ) * drho_air(k) * ddzw(k)                                                           &
             ) * 0.5_wp

       tend(k,j,i) = tend(k,j,i) - (                                                               &
                      ( flux_r(k) + diss_r(k)                                                      &
                    -   flux_l_v(k,j,tn) - diss_l_v(k,j,tn)   ) * ddx                              &
                    + ( flux_n(k) + diss_n(k)                                                      &
                    -   flux_s_v(k,tn) - diss_s_v(k,tn)       ) * ddy                              &
                    + ( ( flux_t(k) + diss_t(k) )                                                  &
                    -   ( flux_d    + diss_d    )                                                  &
                                                ) * drho_air(k) * ddzw(k)                          &
                                   ) + v(k,j,i) * div

       flux_l_v(k,j,tn) = flux_r(k)
       diss_l_v(k,j,tn) = diss_r(k)
       flux_s_v(k,tn)   = flux_n(k)
       diss_s_v(k,tn)   = diss_n(k)
!
!--    Statistical Evaluation of v'v'. The factor has to be applied for right evaluation when
!--    gallilei_trans = .T. .
       sums_vs2_ws_l(k,tn) = sums_vs2_ws_l(k,tn) +                                                 &
                             ( flux_n(k)                                                           &
                               * ( v_comp(k) - 2.0_wp * hom(k,1,2,0)                   )           &
                               / ( v_comp(k) - gv + SIGN( 1.0E-20_wp, v_comp(k) - gv ) )           &
                             + diss_n(k)                                                           &
                               *   ABS( v_comp(k) - 2.0_wp * hom(k,1,2,0)              )           &
                               / ( ABS( v_comp(k) - gv ) + 1.0E-20_wp                  )           &
                             ) *   weight_substep(intermediate_timestep_count)
!
!--    Statistical Evaluation of w'u'.
       sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn) +                                               &
                              ( flux_t(k)                                                          &
                                * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                   )          &
                                / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )           )          &
                              + diss_t(k)                                                          &
                                *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)              )          &
                                / ( ABS( w_comp(k) ) + 1.0E-20_wp                       )          &
                              ) *   weight_substep(intermediate_timestep_count)

    ENDDO


 END SUBROUTINE advec_v_ws_ij



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of w-component - Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_w_ws_ij( i, j, i_omp, tn )


    INTEGER(iwp) ::  i         !< grid index along x-direction
    INTEGER(iwp) ::  i_omp     !< leftmost index on subdomain, or in case of OpenMP, on thread
    INTEGER(iwp) ::  j         !< grid index along y-direction
    INTEGER(iwp) ::  k         !< grid index along z-direction
    INTEGER(iwp) ::  k_mm      !< k-2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_pp      !< k+2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_ppp     !< k+3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  nzb_max_l !< index indicating upper bound for order degradation of horizontal advection terms
    INTEGER(iwp) ::  tn        !< number of OpenMP thread

    REAL(wp)    ::  diss_d  !< discretized artificial dissipation at top of the grid box
    REAL(wp)    ::  div     !< divergence on w-grid
    REAL(wp)    ::  flux_d  !< discretized 6th-order flux at top of the grid box
    REAL(wp)    ::  gu      !< Galilei-transformation velocity along x
    REAL(wp)    ::  gv      !< Galilei-transformation velocity along y
    REAL(wp)    ::  ibit18  !< flag indicating 1st-order scheme along x-direction
    REAL(wp)    ::  ibit19  !< flag indicating 3rd-order scheme along x-direction
    REAL(wp)    ::  ibit20  !< flag indicating 5th-order scheme along x-direction
    REAL(wp)    ::  ibit21  !< flag indicating 1st-order scheme along y-direction
    REAL(wp)    ::  ibit22  !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)    ::  ibit23  !< flag indicating 5th-order scheme along y-direction
    REAL(wp)    ::  ibit24  !< flag indicating 1st-order scheme along z-direction
    REAL(wp)    ::  ibit25  !< flag indicating 3rd-order scheme along z-direction
    REAL(wp)    ::  ibit26  !< flag indicating 5th-order scheme along z-direction

    REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_n !< discretized artificial dissipation at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_r !< discretized artificial dissipation at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  diss_t !< discretized artificial dissipation at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_n !< discretized 6th-order flux at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_r !< discretized 6th-order flux at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  flux_t !< discretized 6th-order flux at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  u_comp !< advection velocity along x
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  v_comp !< advection velocity along y
    REAL(wp), DIMENSION(nzb:nzt+1)  ::  w_comp !< advection velocity along z
!
!-- Used local modified copy of nzb_max (used to degrade order of discretization) at non-cyclic
!-- boundaries. Modify only at relevant points instead of the entire subdomain. This should lead to
!-- better load balance between boundary and non-boundary PEs.
    IF( ( bc_dirichlet_l  .OR.  bc_radiation_l )  .AND.  i <= nxl + 2  .OR.                        &
        ( bc_dirichlet_r  .OR.  bc_radiation_r )  .AND.  i >= nxr - 2  .OR.                        &
        ( bc_dirichlet_s  .OR.  bc_radiation_s )  .AND.  j <= nys + 2  .OR.                        &
        ( bc_dirichlet_n  .OR.  bc_radiation_n )  .AND.  j >= nyn - 2 )  THEN
       nzb_max_l = nzt - 1
    ELSE
       nzb_max_l = nzb_max
    END IF

    gu = 2.0_wp * u_gtrans
    gv = 2.0_wp * v_gtrans
!
!-- Compute southside fluxes for the respective boundary.
    IF ( j == nys )  THEN

       DO  k = nzb+1, nzb_max_l
          ibit23 = REAL( IBITS(advc_flags_m(k,j-1,i),23,1), KIND = wp )
          ibit22 = REAL( IBITS(advc_flags_m(k,j-1,i),22,1), KIND = wp )
          ibit21 = REAL( IBITS(advc_flags_m(k,j-1,i),21,1), KIND = wp )

          v_comp(k)      = v(k+1,j,i) + v(k,j,i) - gv
          flux_s_w(k,tn) = v_comp(k) * (                                                           &
                           (  37.0_wp * ibit23 * adv_mom_5                                         &
                           +   7.0_wp * ibit22 * adv_mom_3                                         &
                           +            ibit21 * adv_mom_1 ) * ( w(k,j,i)   + w(k,j-1,i) )         &
                           - ( 8.0_wp * ibit23 * adv_mom_5                                         &
                           +            ibit22 * adv_mom_3 ) * ( w(k,j+1,i) + w(k,j-2,i) )         &
                           + (          ibit23 * adv_mom_5 ) * ( w(k,j+2,i) + w(k,j-3,i) )         &
                                       )

          diss_s_w(k,tn) = -ABS( v_comp(k) ) * (                                                   &
                           (  10.0_wp * ibit23 * adv_mom_5                                         &
                           +   3.0_wp * ibit22 * adv_mom_3                                         &
                           +            ibit21 * adv_mom_1 ) * ( w(k,j,i)   - w(k,j-1,i) )         &
                           - ( 5.0_wp * ibit23 * adv_mom_5                                         &
                           +            ibit22 * adv_mom_3 ) * ( w(k,j+1,i) - w(k,j-2,i) )         &
                           + (          ibit23 * adv_mom_5 ) * ( w(k,j+2,i) - w(k,j-3,i) )         &
                                               )

       ENDDO

       DO  k = nzb_max_l+1, nzt-1

          v_comp(k)      = v(k+1,j,i) + v(k,j,i) - gv
          flux_s_w(k,tn) = v_comp(k) * (                                                           &
                             37.0_wp * ( w(k,j,i)   + w(k,j-1,i) )                                 &
                           -  8.0_wp * ( w(k,j+1,i) + w(k,j-2,i) )                                 &
                           +           ( w(k,j+2,i) + w(k,j-3,i) ) ) * adv_mom_5
          diss_s_w(k,tn) = -ABS( v_comp(k) ) * (                                                   &
                             10.0_wp * ( w(k,j,i)   - w(k,j-1,i) )                                 &
                           -  5.0_wp * ( w(k,j+1,i) - w(k,j-2,i) )                                 &
                           +           ( w(k,j+2,i) - w(k,j-3,i) ) ) * adv_mom_5

       ENDDO

    ENDIF
!
!-- Compute leftside fluxes for the respective boundary.
    IF ( i == i_omp ) THEN

       DO  k = nzb+1, nzb_max_l

          ibit20 = REAL( IBITS(advc_flags_m(k,j,i-1),20,1), KIND = wp )
          ibit19 = REAL( IBITS(advc_flags_m(k,j,i-1),19,1), KIND = wp )
          ibit18 = REAL( IBITS(advc_flags_m(k,j,i-1),18,1), KIND = wp )

          u_comp(k)        = u(k+1,j,i) + u(k,j,i) - gu
          flux_l_w(k,j,tn) = u_comp(k) * (                                                         &
                             (  37.0_wp * ibit20 * adv_mom_5                                       &
                             +   7.0_wp * ibit19 * adv_mom_3                                       &
                             +            ibit18 * adv_mom_1 ) * ( w(k,j,i)   + w(k,j,i-1) )       &
                             - ( 8.0_wp * ibit20 * adv_mom_5                                       &
                             +            ibit19 * adv_mom_3 ) * ( w(k,j,i+1) + w(k,j,i-2) )       &
                             + (          ibit20 * adv_mom_5 ) * ( w(k,j,i+2) + w(k,j,i-3) )       &
                                         )

          diss_l_w(k,j,tn) = -ABS( u_comp(k) ) * (                                                 &
                             (  10.0_wp * ibit20 * adv_mom_5                                       &
                             +   3.0_wp * ibit19 * adv_mom_3                                       &
                             +            ibit18 * adv_mom_1 ) * ( w(k,j,i)   - w(k,j,i-1) )       &
                             - ( 5.0_wp * ibit20 * adv_mom_5                                       &
                             +            ibit19 * adv_mom_3 ) * ( w(k,j,i+1) - w(k,j,i-2) )       &
                             + (          ibit20 * adv_mom_5 ) * ( w(k,j,i+2) - w(k,j,i-3) )       &
                                                 )

       ENDDO

       DO  k = nzb_max_l+1, nzt-1

          u_comp(k)        = u(k+1,j,i) + u(k,j,i) - gu
          flux_l_w(k,j,tn) = u_comp(k) * (                                                         &
                               37.0_wp * ( w(k,j,i)   + w(k,j,i-1) )                               &
                             -  8.0_wp * ( w(k,j,i+1) + w(k,j,i-2) )                               &
                             +           ( w(k,j,i+2) + w(k,j,i-3) ) ) * adv_mom_5
          diss_l_w(k,j,tn) = -ABS( u_comp(k) ) * (                                                 &
                               10.0_wp * ( w(k,j,i)   - w(k,j,i-1) )                               &
                             -  5.0_wp * ( w(k,j,i+1) - w(k,j,i-2) )                               &
                             +           ( w(k,j,i+2) - w(k,j,i-3) ) ) * adv_mom_5

       ENDDO

    ENDIF
!
!-- Now compute the fluxes and tendency terms for the horizontal and vertical parts.
    DO  k = nzb+1, nzb_max_l

       ibit20 = REAL( IBITS(advc_flags_m(k,j,i),20,1), KIND = wp )
       ibit19 = REAL( IBITS(advc_flags_m(k,j,i),19,1), KIND = wp )
       ibit18 = REAL( IBITS(advc_flags_m(k,j,i),18,1), KIND = wp )

       u_comp(k) = u(k+1,j,i+1) + u(k,j,i+1) - gu
       flux_r(k) = u_comp(k) * (                                                                   &
                   (  37.0_wp * ibit20 * adv_mom_5                                                 &
                   +   7.0_wp * ibit19 * adv_mom_3                                                 &
                   +            ibit18 * adv_mom_1 ) * ( w(k,j,i+1) + w(k,j,i)   )                 &
                   - ( 8.0_wp * ibit20 * adv_mom_5                                                 &
                   +            ibit19 * adv_mom_3 ) * ( w(k,j,i+2) + w(k,j,i-1) )                 &
                   + (          ibit20 * adv_mom_5 ) * ( w(k,j,i+3) + w(k,j,i-2) )                 &
                               )

       diss_r(k) = -ABS( u_comp(k) ) * (                                                           &
                   (  10.0_wp * ibit20 * adv_mom_5                                                 &
                   +   3.0_wp * ibit19 * adv_mom_3                                                 &
                   +            ibit18 * adv_mom_1 ) * ( w(k,j,i+1) - w(k,j,i)   )                 &
                   - ( 5.0_wp * ibit20 * adv_mom_5                                                 &
                   +            ibit19 * adv_mom_3 ) * ( w(k,j,i+2) - w(k,j,i-1) )                 &
                   + (          ibit20 * adv_mom_5 ) * ( w(k,j,i+3) - w(k,j,i-2) )                 &
                                       )

       ibit23 = REAL( IBITS(advc_flags_m(k,j,i),23,1), KIND = wp )
       ibit22 = REAL( IBITS(advc_flags_m(k,j,i),22,1), KIND = wp )
       ibit21 = REAL( IBITS(advc_flags_m(k,j,i),21,1), KIND = wp )

       v_comp(k) = v(k+1,j+1,i) + v(k,j+1,i) - gv
       flux_n(k) = v_comp(k) * (                                                                   &
                  (  37.0_wp * ibit23 * adv_mom_5                                                  &
                  +   7.0_wp * ibit22 * adv_mom_3                                                  &
                  +            ibit21 * adv_mom_1 ) * ( w(k,j+1,i) + w(k,j,i)   )                  &
                  - ( 8.0_wp * ibit23 * adv_mom_5                                                  &
                  +            ibit22 * adv_mom_3 ) * ( w(k,j+2,i) + w(k,j-1,i) )                  &
                  + (          ibit23 * adv_mom_5)  * ( w(k,j+3,i) + w(k,j-2,i) )                  &
                               )

       diss_n(k) = -ABS( v_comp(k) ) * (                                                           &
                   (  10.0_wp * ibit23 * adv_mom_5                                                 &
                   +   3.0_wp * ibit22 * adv_mom_3                                                 &
                   +            ibit21 * adv_mom_1 ) * ( w(k,j+1,i) - w(k,j,i)   )                 &
                   - ( 5.0_wp * ibit23 * adv_mom_5                                                 &
                   +            ibit22 * adv_mom_3 ) * ( w(k,j+2,i)  - w(k,j-1,i) )                &
                   + (          ibit23 * adv_mom_5 ) * ( w(k,j+3,i)  - w(k,j-2,i) )                &
                                       )
    ENDDO

    DO  k = nzb_max_l+1, nzt-1

       u_comp(k) = u(k+1,j,i+1) + u(k,j,i+1) - gu
       flux_r(k) = u_comp(k) * (                                                                   &
                     37.0_wp * ( w(k,j,i+1) + w(k,j,i)   )                                         &
                   -  8.0_wp * ( w(k,j,i+2) + w(k,j,i-1) )                                         &
                   +           ( w(k,j,i+3) + w(k,j,i-2) ) ) * adv_mom_5

       diss_r(k) = -ABS( u_comp(k) ) * (                                                           &
                     10.0_wp * ( w(k,j,i+1) - w(k,j,i)   )                                         &
                   -  5.0_wp * ( w(k,j,i+2) - w(k,j,i-1) )                                         &
                   +           ( w(k,j,i+3) - w(k,j,i-2) ) ) * adv_mom_5

       v_comp(k) = v(k+1,j+1,i) + v(k,j+1,i) - gv
       flux_n(k) = v_comp(k) * (                                                                   &
                     37.0_wp * ( w(k,j+1,i) + w(k,j,i)   )                                         &
                   -  8.0_wp * ( w(k,j+2,i) + w(k,j-1,i) )                                         &
                   +           ( w(k,j+3,i) + w(k,j-2,i) ) ) * adv_mom_5

       diss_n(k) = -ABS( v_comp(k) ) * (                                                           &
                     10.0_wp * ( w(k,j+1,i) - w(k,j,i)   )                                         &
                   -  5.0_wp * ( w(k,j+2,i) - w(k,j-1,i) )                                         &
                   +           ( w(k,j+3,i) - w(k,j-2,i) ) ) * adv_mom_5
    ENDDO

!
!-- Now, compute vertical fluxes. Split loop into a part treating the lowest grid points with
!-- indirect indexing, a main loop without indirect indexing, and a loop for the uppermost grid
!-- points with indirect indexing. This allows better vectorization for the main loop.
!-- First, compute the flux at model surface, which need has to be calculated explicitly for the
!-- tendency at the first w-level. For topography wall this is done implicitely by advc_flags_m.
!-- First, compute flux at lowest level, located at z=dz/2.
    k         = nzb + 1
    w_comp(k) = w(k,j,i) + w(k-1,j,i)
    flux_t(0) =  w_comp(k)        * rho_air(k) * ( w(k,j,i) + w(k-1,j,i) ) * adv_mom_1
    diss_t(0) = -ABS( w_comp(k) ) * rho_air(k) * ( w(k,j,i) - w(k-1,j,i) ) * adv_mom_1

    DO  k = nzb+1, nzb+1
!
!--    k index has to be modified near bottom and top, else array subscripts will be exceeded.
       ibit26 = REAL( IBITS(advc_flags_m(k,j,i),26,1), KIND = wp )
       ibit25 = REAL( IBITS(advc_flags_m(k,j,i),25,1), KIND = wp )
       ibit24 = REAL( IBITS(advc_flags_m(k,j,i),24,1), KIND = wp )

       k_ppp = k + 3 * ibit26
       k_pp  = k + 2 * ( 1 - ibit24  )
       k_mm  = k - 2 * ibit26

       w_comp(k) = w(k+1,j,i) + w(k,j,i)
       flux_t(k) = w_comp(k) * rho_air(k+1) * (                                                    &
                   (  37.0_wp * ibit26 * adv_mom_5                                                 &
                   +   7.0_wp * ibit25 * adv_mom_3                                                 &
                   +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)   + w(k,j,i)    )              &
                   - ( 8.0_wp * ibit26 * adv_mom_5                                                 &
                   +            ibit25 * adv_mom_3 ) * ( w(k_pp,j,i)  + w(k-1,j,i)  )              &
                   + (          ibit26 * adv_mom_5 ) * ( w(k_ppp,j,i) + w(k_mm,j,i) )              &
                                              )

       diss_t(k) = -ABS( w_comp(k) ) * rho_air(k+1) * (                                            &
                   (  10.0_wp * ibit26 * adv_mom_5                                                 &
                   +   3.0_wp * ibit25 * adv_mom_3                                                 &
                   +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)   - w(k,j,i)    )              &
                   - ( 5.0_wp * ibit26 * adv_mom_5                                                 &
                   +            ibit25 * adv_mom_3 ) * ( w(k_pp,j,i)  - w(k-1,j,i)  )              &
                   + (          ibit26 * adv_mom_5 ) * ( w(k_ppp,j,i) - w(k_mm,j,i) )              &
                                                      )
    ENDDO

    DO  k = nzb+2, nzt-2

       ibit26 = REAL( IBITS(advc_flags_m(k,j,i),26,1), KIND = wp )
       ibit25 = REAL( IBITS(advc_flags_m(k,j,i),25,1), KIND = wp )
       ibit24 = REAL( IBITS(advc_flags_m(k,j,i),24,1), KIND = wp )

       w_comp(k) = w(k+1,j,i) + w(k,j,i)
       flux_t(k) = w_comp(k) * rho_air(k+1) * (                                                    &
                   (  37.0_wp * ibit26 * adv_mom_5                                                 &
                   +   7.0_wp * ibit25 * adv_mom_3                                                 &
                   +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)  + w(k,j,i)   )                &
                   - ( 8.0_wp * ibit26 * adv_mom_5                                                 &
                   +            ibit25 * adv_mom_3 ) * ( w(k+2,j,i)  + w(k-1,j,i) )                &
                   + (          ibit26 * adv_mom_5 ) * ( w(k+3,j,i)  + w(k-2,j,i) )                &
                                              )

       diss_t(k) = -ABS( w_comp(k) ) * rho_air(k+1) * (                                            &
                   (  10.0_wp * ibit26 * adv_mom_5                                                 &
                   +   3.0_wp * ibit25 * adv_mom_3                                                 &
                   +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i) - w(k,j,i)    )                &
                   - ( 5.0_wp * ibit26 * adv_mom_5                                                 &
                   +            ibit25 * adv_mom_3 ) * ( w(k+2,j,i) - w(k-1,j,i)  )                &
                   + (          ibit26 * adv_mom_5 ) * ( w(k+3,j,i) - w(k-2,j,i)  )                &
                                                      )
    ENDDO

    DO  k = nzt-1, nzt-1
!
!--    k index has to be modified near bottom and top, else array subscripts will be exceeded.
       ibit26 = REAL( IBITS(advc_flags_m(k,j,i),26,1), KIND = wp )
       ibit25 = REAL( IBITS(advc_flags_m(k,j,i),25,1), KIND = wp )
       ibit24 = REAL( IBITS(advc_flags_m(k,j,i),24,1), KIND = wp )
!
!--    k index needs to be modified near bottom and top, else array subscripts will be exceeded.
!--    (k_ppp = k+3, k_mm = k) if 5th order scheme is employed, else (k_ppp = k, k_mm = k).
!--    Special treatment is required for k_pp. This must be limited at k = nzt, where either
!--    the first order scheme is employed (ibit24), or all bit flags are zero (in case there is a
!--    rigid lid defined at or near the domain top).
       k_ppp = k + 3 * ibit26
       k_pp  = k + 2 * ( 1 - ibit24  ) * ( ibit24 + ibit25 + ibit26 )
       k_mm  = k - 2 * ibit26

       w_comp(k) = w(k+1,j,i) + w(k,j,i)
       flux_t(k) = w_comp(k) * rho_air(k+1) * (                                                    &
                   (  37.0_wp * ibit26 * adv_mom_5                                                 &
                   +   7.0_wp * ibit25 * adv_mom_3                                                 &
                   +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)   + w(k,j,i)    )              &
                   - ( 8.0_wp * ibit26 * adv_mom_5                                                 &
                   +            ibit25 * adv_mom_3 ) * ( w(k_pp,j,i)  + w(k-1,j,i)  )              &
                   + (          ibit26 * adv_mom_5 ) * ( w(k_ppp,j,i) + w(k_mm,j,i) )              &
                                              )

       diss_t(k) = -ABS( w_comp(k) ) * rho_air(k+1) * (                                            &
                   (  10.0_wp * ibit26 * adv_mom_5                                                 &
                   +   3.0_wp * ibit25 * adv_mom_3                                                 &
                   +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)   - w(k,j,i)    )              &
                   - ( 5.0_wp * ibit26 * adv_mom_5                                                 &
                   +            ibit25 * adv_mom_3 ) * ( w(k_pp,j,i)  - w(k-1,j,i)  )              &
                   + (          ibit26 * adv_mom_5 ) * ( w(k_ppp,j,i) - w(k_mm,j,i) )              &
                                                      )
    ENDDO

!
!-- Set resolved/turbulent flux at model top to zero (w-level). Hint: The flux at nzt is defined at
!-- the scalar grid point nzt+1. Therefore, the flux at nzt+1 is already outside of the model
!-- domain.
    flux_t(nzt) = 0.0_wp
    diss_t(nzt) = 0.0_wp
    w_comp(nzt) = 0.0_wp

    flux_t(nzt+1) = 0.0_wp
    diss_t(nzt+1) = 0.0_wp
    w_comp(nzt+1) = 0.0_wp

    DO  k = nzb+1, nzb_max_l

       flux_d    = flux_t(k-1)
       diss_d    = diss_t(k-1)

       ibit20 = REAL( IBITS(advc_flags_m(k,j,i),20,1), KIND = wp )
       ibit19 = REAL( IBITS(advc_flags_m(k,j,i),19,1), KIND = wp )
       ibit18 = REAL( IBITS(advc_flags_m(k,j,i),18,1), KIND = wp )

       ibit23 = REAL( IBITS(advc_flags_m(k,j,i),23,1), KIND = wp )
       ibit22 = REAL( IBITS(advc_flags_m(k,j,i),22,1), KIND = wp )
       ibit21 = REAL( IBITS(advc_flags_m(k,j,i),21,1), KIND = wp )

       ibit26 = REAL( IBITS(advc_flags_m(k,j,i),26,1), KIND = wp )
       ibit25 = REAL( IBITS(advc_flags_m(k,j,i),25,1), KIND = wp )
       ibit24 = REAL( IBITS(advc_flags_m(k,j,i),24,1), KIND = wp )
!
!--    Calculate the divergence of the velocity field. A respective correction is needed to overcome
!--    numerical instabilities introduced by an insufficient reduction of divergences near
!--    topography.
       div = ( ( ( u_comp(k) + gu ) * ( ibit18 + ibit19 + ibit20 )                                 &
               - ( u(k+1,j,i) + u(k,j,i)   )                                                       &
                                 * (                                                               &
                  REAL( IBITS(advc_flags_m(k,j,i-1),18,1), KIND = wp )                             &
                + REAL( IBITS(advc_flags_m(k,j,i-1),19,1), KIND = wp )                             &
                + REAL( IBITS(advc_flags_m(k,j,i-1),20,1), KIND = wp )                             &
                                   )                                                               &
               ) * ddx                                                                             &
             + ( ( v_comp(k) + gv ) * ( ibit21 + ibit22 + ibit23 )                                 &
               - ( v(k+1,j,i) + v(k,j,i)   )                                                       &
                                 * (                                                               &
                  REAL( IBITS(advc_flags_m(k,j-1,i),21,1), KIND = wp )                             &
                + REAL( IBITS(advc_flags_m(k,j-1,i),22,1), KIND = wp )                             &
                + REAL( IBITS(advc_flags_m(k,j-1,i),23,1), KIND = wp )                             &
                                   )                                                               &
               ) * ddy                                                                             &
             + ( w_comp(k)               * rho_air(k+1)                                            &
                                         * ( ibit24 + ibit25 + ibit26 )                            &
             - ( w(k,j,i) + w(k-1,j,i) ) * rho_air(k)                                              &
                                 * (                                                               &
                  REAL( IBITS(advc_flags_m(k-1,j,i),24,1), KIND = wp )                             &
                + REAL( IBITS(advc_flags_m(k-1,j,i),25,1), KIND = wp )                             &
                + REAL( IBITS(advc_flags_m(k-1,j,i),26,1), KIND = wp )                             &
                                   )                                                               &
               ) * drho_air_zw(k) * ddzu(k+1)                                                      &
             ) * 0.5_wp

       tend(k,j,i) = tend(k,j,i) - (                                                               &
                     ( flux_r(k) + diss_r(k)                                                       &
                   -   flux_l_w(k,j,tn) - diss_l_w(k,j,tn)   ) * ddx                               &
                   + ( flux_n(k) + diss_n(k)                                                       &
                   -   flux_s_w(k,tn) - diss_s_w(k,tn)       ) * ddy                               &
                   + ( ( flux_t(k) + diss_t(k) )                                                   &
                   -   ( flux_d    + diss_d    )                                                   &
                                           ) * drho_air_zw(k) * ddzu(k+1)                          &
                                   ) + div * w(k,j,i)

       flux_l_w(k,j,tn) = flux_r(k)
       diss_l_w(k,j,tn) = diss_r(k)
       flux_s_w(k,tn)   = flux_n(k)
       diss_s_w(k,tn)   = diss_n(k)
!
!--    Statistical Evaluation of w'w'.
       sums_ws2_ws_l(k,tn)  = sums_ws2_ws_l(k,tn) +                                                &
                              ( flux_t(k)                                                          &
                              * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                )               &
                              / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )        )               &
                               + diss_t(k)                                                         &
                              *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)           )               &
                              / ( ABS( w_comp(k) ) + 1.0E-20_wp                    )               &
                              ) *   weight_substep(intermediate_timestep_count)

    ENDDO

    DO  k = nzb_max_l+1, nzt-1

       flux_d    = flux_t(k-1)
       diss_d    = diss_t(k-1)
!
!--    Calculate the divergence of the velocity field. A respective correction is needed to overcome
!--    numerical instabilities introduced by an insufficient reduction of divergences near
!--    topography.
       div = ( ( u_comp(k) + gu - ( u(k+1,j,i) + u(k,j,i)   ) ) * ddx                              &
             + ( v_comp(k) + gv - ( v(k+1,j,i) + v(k,j,i)   ) ) * ddy                              &
             + ( w_comp(k)               * rho_air(k+1)                                            &
             - ( w(k,j,i) + w(k-1,j,i) ) * rho_air(k)                                              &
               ) * drho_air_zw(k) * ddzu(k+1)                                                      &
             ) * 0.5_wp

       tend(k,j,i) = tend(k,j,i) - (                                                               &
                       ( flux_r(k) + diss_r(k)                                                     &
                     -   flux_l_w(k,j,tn) - diss_l_w(k,j,tn)   ) * ddx                             &
                     + ( flux_n(k) + diss_n(k)                                                     &
                     -   flux_s_w(k,tn) - diss_s_w(k,tn)       ) * ddy                             &
                     + ( ( flux_t(k) + diss_t(k) )                                                 &
                     -   ( flux_d    + diss_d    )                                                 &
                                              ) * drho_air_zw(k) * ddzu(k+1)                       &
                                   ) + div * w(k,j,i)

       flux_l_w(k,j,tn) = flux_r(k)
       diss_l_w(k,j,tn) = diss_r(k)
       flux_s_w(k,tn)   = flux_n(k)
       diss_s_w(k,tn)   = diss_n(k)
!
!--    Statistical Evaluation of w'w'.
       sums_ws2_ws_l(k,tn)  = sums_ws2_ws_l(k,tn) +                                                &
                              ( flux_t(k)                                                          &
                                * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                )             &
                                / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )        )             &
                              + diss_t(k)                                                          &
                                *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)           )             &
                                / ( ABS( w_comp(k) ) + 1.0E-20_wp                    )             &
                              ) *   weight_substep(intermediate_timestep_count)

    ENDDO

 END SUBROUTINE advec_w_ws_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Scalar advection - Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_s_ws( advc_flag, sk, sk_char,                                                    &
                        non_cyclic_l, non_cyclic_n,                                                &
                        non_cyclic_r, non_cyclic_s )


    CHARACTER (LEN = *), INTENT(IN) ::  sk_char !< string identifier, used for assign fluxes
                                                !< to the correct dimension in the analysis array

    INTEGER(iwp) ::  i         !< grid index along x-direction
    INTEGER(iwp) ::  j         !< grid index along y-direction
    INTEGER(iwp) ::  k         !< grid index along z-direction
    INTEGER(iwp) ::  k_mm      !< k-2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_pp      !< k+2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_ppp     !< k+3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  nzb_max_l !< index indicating upper bound for order degradation of horizontal advection terms
    INTEGER(iwp) ::  sk_num    !< integer identifier, used for assign fluxes to the correct dimension in the analysis array
    INTEGER(iwp) ::  tn = 0    !< number of OpenMP thread (is always zero here)

    INTEGER(iwp), INTENT(IN), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::                          &
                                               advc_flag !< flag array to control order of scalar advection

    LOGICAL  ::  non_cyclic_l    !< flag that indicates non-cyclic boundary on the left
    LOGICAL  ::  non_cyclic_n    !< flag that indicates non-cyclic boundary on the north
    LOGICAL  ::  non_cyclic_r    !< flag that indicates non-cyclic boundary on the right
    LOGICAL  ::  non_cyclic_s    !< flag that indicates non-cyclic boundary on the south

    REAL(wp) ::  diss_d        !< artificial dissipation term at grid box bottom
    REAL(wp) ::  div           !< velocity diverence on scalar grid
    REAL(wp) ::  flux_d        !< 6th-order flux at grid box bottom
    REAL(wp) ::  ibit0         !< flag indicating 1st-order scheme along x-direction
    REAL(wp) ::  ibit1         !< flag indicating 3rd-order scheme along x-direction
    REAL(wp) ::  ibit2         !< flag indicating 5th-order scheme along x-direction
    REAL(wp) ::  ibit3         !< flag indicating 1st-order scheme along y-direction
    REAL(wp) ::  ibit4         !< flag indicating 3rd-order scheme along y-direction
    REAL(wp) ::  ibit5         !< flag indicating 5th-order scheme along y-direction
    REAL(wp) ::  ibit6         !< flag indicating 1st-order scheme along z-direction
    REAL(wp) ::  ibit7         !< flag indicating 3rd-order scheme along z-direction
    REAL(wp) ::  ibit8         !< flag indicating 5th-order scheme along z-direction
#ifdef _OPENACC
    REAL(wp) ::  ibit0_l  !< flag indicating 1st-order scheme along x-direction
    REAL(wp) ::  ibit1_l  !< flag indicating 3rd-order scheme along x-direction
    REAL(wp) ::  ibit2_l  !< flag indicating 5th-order scheme along x-direction
    REAL(wp) ::  ibit3_s  !< flag indicating 1st-order scheme along y-direction
    REAL(wp) ::  ibit4_s  !< flag indicating 3rd-order scheme along y-direction
    REAL(wp) ::  ibit5_s  !< flag indicating 5th-order scheme along y-direction
#endif
    REAL(wp) ::  u_comp        !< advection velocity along x-direction
    REAL(wp) ::  v_comp        !< advection velocity along y-direction
#ifdef _OPENACC
    REAL(wp) ::  u_comp_l !< advection velocity along x-direction
    REAL(wp) ::  v_comp_s !< advection velocity along y-direction
#endif
!
!--    sk is an array from parameter list. It should not be a pointer, because in that case the
!--    compiler can not assume a stride 1 and cannot perform a strided one vector load. Adding the
!--    CONTIGUOUS keyword makes things even worse, because the compiler cannot assume strided one in
!--    the caller side.


    REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_n     !< discretized artificial dissipation at northward-side
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_r     !< discretized artificial dissipation at rightward-side
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_t     !< discretized artificial dissipation at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_n     !< discretized 6th-order flux at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_r     !< discretized 6th-order flux at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_t     !< discretized 6th-order flux at top of the grid box

    REAL(wp), DIMENSION(nzb+1:nzt)         ::  diss_l            !< discretized artificial dissipation at leftward-side
    REAL(wp), DIMENSION(nzb+1:nzt)         ::  diss_s            !< discretized artificial dissipation at southward-side
    REAL(wp), DIMENSION(nzb+1:nzt)         ::  flux_l            !< discretized 6th-order flux at leftward-side
    REAL(wp), DIMENSION(nzb+1:nzt)         ::  flux_s            !< discretized 6th-order flux at southward-side
#ifndef _OPENACC
    REAL(wp), DIMENSION(nzb+1:nzt)         ::  swap_diss_y_local !< discretized artificial dissipation at southward-side
    REAL(wp), DIMENSION(nzb+1:nzt)         ::  swap_flux_y_local !< discretized 6th-order flux at northward-side

    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_diss_x_local !< discretized artificial dissipation at leftward-side
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_flux_x_local !< discretized 6th-order flux at leftward-side
#endif

    REAL(wp), INTENT(IN),DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<  advected scalar

    CALL cpu_log( log_point_s(49), 'advec_s_ws', 'start' )

    SELECT CASE ( sk_char )

        CASE ( 'pt' )
           sk_num = 1
        CASE ( 'sa' )
           sk_num = 2
        CASE ( 'q' )
           sk_num = 3
        CASE ( 'qc' )
           sk_num = 4
        CASE ( 'qr' )
           sk_num = 5
        CASE ( 'nc' )
           sk_num = 6
        CASE ( 'nr' )
           sk_num = 7
        CASE ( 's' )
           sk_num = 8
        CASE ( 'aerosol_mass', 'aerosol_number', 'salsa_gas' )
           sk_num = 9
        CASE ( 'ni' )
           sk_num = 10
        CASE ( 'qi' )
           sk_num = 11
        CASE ( 'ng' )
           sk_num = 12
        CASE ( 'qg' )
           sk_num = 13
        CASE ( 'ns' )
           sk_num = 14
        CASE ( 'qs' )
           sk_num = 15
        CASE DEFAULT
           sk_num = 0

    END SELECT

    !$ACC PARALLEL LOOP COLLAPSE(2) FIRSTPRIVATE(tn, sk_num) &
    !$ACC PRIVATE(i, j, k, k_mm, k_pp, k_ppp) &
    !$ACC PRIVATE(ibit0, ibit1, ibit2, ibit3, ibit4, ibit5) &
    !$ACC PRIVATE(ibit0_l, ibit1_l, ibit2_l) &
    !$ACC PRIVATE(ibit3_s, ibit4_s, ibit5_s) &
    !$ACC PRIVATE(ibit6, ibit7, ibit8) &
    !$ACC PRIVATE(nzb_max_l) &
    !$ACC PRIVATE(diss_l, diss_r, flux_l, flux_r) &
    !$ACC PRIVATE(diss_n, diss_s, flux_s, flux_n) &
    !$ACC PRIVATE(flux_t, diss_t, flux_d, diss_d) &
    !$ACC PRIVATE(div, u_comp, u_comp_l, v_comp, v_comp_s) &
    !$ACC PRESENT(advc_flag) &
    !$ACC PRESENT(sk, u, v, w, u_stokes_zu, v_stokes_zu) &
    !$ACC PRESENT(drho_air, rho_air_zw, ddzw) &
    !$ACC PRESENT(tend) &
    !$ACC PRESENT(hom(:,1,1:3,0)) &
    !$ACC PRESENT(weight_substep(intermediate_timestep_count)) &
    !$ACC PRESENT(sums_wspts_ws_l, sums_wssas_ws_l) &
    !$ACC PRESENT(sums_wsqs_ws_l, sums_wsqcs_ws_l) &
    !$ACC PRESENT(sums_wsqrs_ws_l, sums_wsncs_ws_l) &
    !$ACC PRESENT(sums_wsnrs_ws_l, sums_wsss_ws_l) &
    !$ACC PRESENT(sums_wsnis_ws_l, sums_wsqis_ws_l) &
    !$ACC PRESENT(sums_wsngs_ws_l, sums_wsqgs_ws_l) &
    !$ACC PRESENT(sums_wsnss_ws_l, sums_wsqss_ws_l) &
    !$ACC PRESENT(sums_salsa_ws_l)
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Used local modified copy of nzb_max (used to degrade order of discretization) at
!--       non-cyclic boundaries. Modify only at relevant points instead of the entire subdomain.
!--       This should lead to better load balance between boundary and non-boundary PEs.
          IF( non_cyclic_l  .AND.  i <= nxl + 2  .OR.                                              &
              non_cyclic_r  .AND.  i >= nxr - 2  .OR.                                              &
              non_cyclic_s  .AND.  j <= nys + 2  .OR.                                              &
              non_cyclic_n  .AND.  j >= nyn - 2 )  THEN
             nzb_max_l = nzt
          ELSE
             nzb_max_l = nzb_max
          END IF
#ifndef _OPENACC
!
!--       Compute leftside fluxes of the respective PE bounds.
          IF ( i == nxl )  THEN

             DO  k = nzb+1, nzb_max_l

                ibit2 = REAL( IBITS(advc_flag(k,j,i-1),2,1), KIND = wp )
                ibit1 = REAL( IBITS(advc_flag(k,j,i-1),1,1), KIND = wp )
                ibit0 = REAL( IBITS(advc_flag(k,j,i-1),0,1), KIND = wp )
!
!--             Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
                u_comp                    = u(k,j,i) - u_gtrans + u_stokes_zu(k)                   &
                                            * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
                swap_flux_x_local(k,j)    = u_comp * (                                             &
                                              ( 37.0_wp * ibit2 * adv_sca_5                        &
                                            +    7.0_wp * ibit1 * adv_sca_3                        &
                                            +             ibit0 * adv_sca_1                        &
                                              ) *                                                  &
                                              ( sk(k,j,i)   + sk(k,j,i-1)    )                     &
                                            - (  8.0_wp * ibit2 * adv_sca_5                        &
                                            +             ibit1 * adv_sca_3                        &
                                              ) *                                                  &
                                              ( sk(k,j,i+1) + sk(k,j,i-2)    )                     &
                                            + (           ibit2 * adv_sca_5                        &
                                              ) *                                                  &
                                              ( sk(k,j,i+2) + sk(k,j,i-3)    )                     &
                                                     )

                swap_diss_x_local(k,j)    = -ABS( u_comp ) * (                                     &
                                              ( 10.0_wp * ibit2 * adv_sca_5                        &
                                            +    3.0_wp * ibit1 * adv_sca_3                        &
                                            +             ibit0 * adv_sca_1                        &
                                              ) *                                                  &
                                              ( sk(k,j,i)   - sk(k,j,i-1)    )                     &
                                            - (  5.0_wp * ibit2 * adv_sca_5                        &
                                            +             ibit1 * adv_sca_3                        &
                                              ) *                                                  &
                                              ( sk(k,j,i+1) - sk(k,j,i-2)    )                     &
                                            + (           ibit2 * adv_sca_5                        &
                                              ) *                                                  &
                                              ( sk(k,j,i+2) - sk(k,j,i-3)    )                     &
                                                             )

             ENDDO

             DO  k = nzb_max_l+1, nzt

                u_comp                    = u(k,j,i) - u_gtrans + u_stokes_zu(k)
                swap_flux_x_local(k,j)    = u_comp * (                                             &
                                              37.0_wp * ( sk(k,j,i)   + sk(k,j,i-1) )              &
                                            -  8.0_wp * ( sk(k,j,i+1) + sk(k,j,i-2) )              &
                                            +           ( sk(k,j,i+2) + sk(k,j,i-3) )              &
                                                     ) * adv_sca_5

                swap_diss_x_local(k,j)    = -ABS( u_comp ) * (                                     &
                                              10.0_wp * ( sk(k,j,i)   - sk(k,j,i-1) )              &
                                            -  5.0_wp * ( sk(k,j,i+1) - sk(k,j,i-2) )              &
                                            +           ( sk(k,j,i+2) - sk(k,j,i-3) )              &
                                                             ) * adv_sca_5

             ENDDO

          ENDIF
!
!--       Compute southside fluxes of the respective PE bounds.
          IF ( j == nys )  THEN
!
!--          Up to the top of the highest topography.
             DO  k = nzb+1, nzb_max_l

                ibit5 = REAL( IBITS(advc_flag(k,j-1,i),5,1), KIND = wp )
                ibit4 = REAL( IBITS(advc_flag(k,j-1,i),4,1), KIND = wp )
                ibit3 = REAL( IBITS(advc_flag(k,j-1,i),3,1), KIND = wp )
!
!--             Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
                v_comp                  = v(k,j,i) - v_gtrans + v_stokes_zu(k)                     &
                                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                swap_flux_y_local(k)    = v_comp *         (                                       &
                                             ( 37.0_wp * ibit5 * adv_sca_5                         &
                                          +     7.0_wp * ibit4 * adv_sca_3                         &
                                          +              ibit3 * adv_sca_1                         &
                                             ) *                                                   &
                                             ( sk(k,j,i)   + sk(k,j-1,i)    )                      &
                                          -  (  8.0_wp * ibit5 * adv_sca_5                         &
                                          +              ibit4 * adv_sca_3                         &
                                             ) *                                                   &
                                             ( sk(k,j+1,i) + sk(k,j-2,i)    )                      &
                                          +  (           ibit5 * adv_sca_5                         &
                                             ) *                                                   &
                                             ( sk(k,j+2,i) + sk(k,j-3,i)    )                      &
                                                           )

                swap_diss_y_local(k)    = -ABS( v_comp ) * (                                       &
                                            ( 10.0_wp * ibit5 * adv_sca_5                          &
                                          +    3.0_wp * ibit4 * adv_sca_3                          &
                                          +             ibit3 * adv_sca_1                          &
                                            ) *                                                    &
                                            ( sk(k,j,i)   - sk(k,j-1,i)  )                         &
                                          - (  5.0_wp * ibit5 * adv_sca_5                          &
                                          +             ibit4 * adv_sca_3                          &
                                            ) *                                                    &
                                            ( sk(k,j+1,i) - sk(k,j-2,i)  )                         &
                                          + (           ibit5 * adv_sca_5                          &
                                            ) *                                                    &
                                            ( sk(k,j+2,i) - sk(k,j-3,i)  )                         &
                                                           )

             ENDDO
!
!--          Above to the top of the highest topography. No degradation necessary.
             DO  k = nzb_max_l+1, nzt

                v_comp                  = v(k,j,i) - v_gtrans + v_stokes_zu(k)
                swap_flux_y_local(k)    = v_comp * (                                               &
                                            37.0_wp * ( sk(k,j,i)   + sk(k,j-1,i) )                &
                                          -  8.0_wp * ( sk(k,j+1,i) + sk(k,j-2,i) )                &
                                          +           ( sk(k,j+2,i) + sk(k,j-3,i) )                &
                                                   ) * adv_sca_5
                swap_diss_y_local(k)    = -ABS( v_comp ) * (                                       &
                                            10.0_wp * ( sk(k,j,i)   - sk(k,j-1,i) )                &
                                          -  5.0_wp * ( sk(k,j+1,i) - sk(k,j-2,i) )                &
                                          +             sk(k,j+2,i) - sk(k,j-3,i)                  &
                                                           ) * adv_sca_5

             ENDDO

          ENDIF
#endif
!
!--       Now compute the fluxes and tendency terms for the horizontal and vertical parts up to the
!--       top of the highest topography.
          DO  k = nzb+1, nzb_max_l
!
!--          Note: It is faster to conduct all multiplications explicitly, e.g. * adv_sca_5 ... than
!--          to determine a factor and multiplicate the flux at the end.
             ibit2 = REAL( IBITS(advc_flag(k,j,i),2,1), KIND = wp )
             ibit1 = REAL( IBITS(advc_flag(k,j,i),1,1), KIND = wp )
             ibit0 = REAL( IBITS(advc_flag(k,j,i),0,1), KIND = wp )
!
!--          Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
             u_comp    = u(k,j,i+1) - u_gtrans + u_stokes_zu(k)                                    &
                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 1 ) )
             flux_r(k) = u_comp * (                                                                &
                         (  37.0_wp * ibit2 * adv_sca_5                                            &
                         +   7.0_wp * ibit1 * adv_sca_3                                            &
                         +            ibit0 * adv_sca_1 ) * ( sk(k,j,i+1) + sk(k,j,i)   )          &
                         - ( 8.0_wp * ibit2 * adv_sca_5                                            &
                         +            ibit1 * adv_sca_3 ) * ( sk(k,j,i+2) + sk(k,j,i-1) )          &
                         + (          ibit2 * adv_sca_5 ) * ( sk(k,j,i+3) + sk(k,j,i-2) )          &
                                  )

             diss_r(k) = -ABS( u_comp ) * (                                                        &
                           (  10.0_wp * ibit2 * adv_sca_5                                          &
                           +   3.0_wp * ibit1 * adv_sca_3                                          &
                           +            ibit0 * adv_sca_1 ) * ( sk(k,j,i+1) - sk(k,j,i)  )         &
                           - ( 5.0_wp * ibit2 * adv_sca_5                                          &
                           +            ibit1 * adv_sca_3 ) * ( sk(k,j,i+2) - sk(k,j,i-1) )        &
                           + (          ibit2 * adv_sca_5 ) * ( sk(k,j,i+3) - sk(k,j,i-2) )        &
                                          )
#ifdef _OPENACC
!
!--          Recompute the left fluxes.
             ibit2_l = REAL( IBITS(advc_flag(k,j,i-1),2,1), KIND = wp )
             ibit1_l = REAL( IBITS(advc_flag(k,j,i-1),1,1), KIND = wp )
             ibit0_l = REAL( IBITS(advc_flag(k,j,i-1),0,1), KIND = wp )
!
!--          Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
             u_comp_l                  = u(k,j,i) - u_gtrans + u_stokes_zu(k)                      &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
             flux_l(k)                 = u_comp_l * (                                              &
                                           ( 37.0_wp * ibit2_l * adv_sca_5                         &
                                         +    7.0_wp * ibit1_l * adv_sca_3                         &
                                         +             ibit0_l * adv_sca_1                         &
                                           ) *                                                     &
                                           ( sk(k,j,i)   + sk(k,j,i-1)    )                        &
                                         - (  8.0_wp * ibit2_l * adv_sca_5                         &
                                         +              ibit1_l * adv_sca_3                        &
                                           ) *                                                     &
                                           ( sk(k,j,i+1) + sk(k,j,i-2)    )                        &
                                         + (           ibit2_l * adv_sca_5                         &
                                           ) *                                                     &
                                           ( sk(k,j,i+2) + sk(k,j,i-3)    )                        &
                                                    )

             diss_l(k)                 = -ABS( u_comp_l ) * (                                      &
                                           ( 10.0_wp * ibit2_l * adv_sca_5                         &
                                         +    3.0_wp * ibit1_l * adv_sca_3                         &
                                         +             ibit0_l * adv_sca_1                         &
                                           ) *                                                     &
                                           ( sk(k,j,i)   - sk(k,j,i-1) )                           &
                                         - (  5.0_wp * ibit2_l * adv_sca_5                         &
                                         +              ibit1_l * adv_sca_3                        &
                                           ) *                                                     &
                                           ( sk(k,j,i+1) - sk(k,j,i-2) )                           &
                                         + (           ibit2_l * adv_sca_5                         &
                                           ) *                                                     &
                                           ( sk(k,j,i+2) - sk(k,j,i-3) )                           &
                                                            )
#else
             flux_l(k) = swap_flux_x_local(k,j)
             diss_l(k) = swap_diss_x_local(k,j)
#endif
             ibit5 = REAL( IBITS(advc_flag(k,j,i),5,1), KIND = wp )
             ibit4 = REAL( IBITS(advc_flag(k,j,i),4,1), KIND = wp )
             ibit3 = REAL( IBITS(advc_flag(k,j,i),3,1), KIND = wp )
!
!--          Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
             v_comp    = v(k,j+1,i) - v_gtrans + v_stokes_zu(k)                                    &
                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j+1,i), 2 ) )
             flux_n(k) = v_comp * (                                                                &
                         (  37.0_wp * ibit5 * adv_sca_5                                            &
                         +   7.0_wp * ibit4 * adv_sca_3                                            &
                         +            ibit3 * adv_sca_1 ) * ( sk(k,j+1,i) + sk(k,j,i)   )          &
                         - ( 8.0_wp * ibit5 * adv_sca_5                                            &
                         +            ibit4 * adv_sca_3 ) * ( sk(k,j+2,i) + sk(k,j-1,i) )          &
                         + (          ibit5 * adv_sca_5 ) * ( sk(k,j+3,i) + sk(k,j-2,i) )          &
                                  )

             diss_n(k) = -ABS( v_comp ) * (                                                        &
                         (  10.0_wp * ibit5 * adv_sca_5                                            &
                         +   3.0_wp * ibit4 * adv_sca_3                                            &
                         +            ibit3 * adv_sca_1 ) * ( sk(k,j+1,i) - sk(k,j,i)   )          &
                         - ( 5.0_wp * ibit5 * adv_sca_5                                            &
                         +            ibit4 * adv_sca_3 ) * ( sk(k,j+2,i) - sk(k,j-1,i) )          &
                         + (          ibit5 * adv_sca_5 ) * ( sk(k,j+3,i) - sk(k,j-2,i) )          &
                                          )
#ifdef _OPENACC
!
!--          Recompute the south fluxes.
             ibit5_s = REAL( IBITS(advc_flag(k,j-1,i),5,1), KIND = wp )
             ibit4_s = REAL( IBITS(advc_flag(k,j-1,i),4,1), KIND = wp )
             ibit3_s = REAL( IBITS(advc_flag(k,j-1,i),3,1), KIND = wp )
!
!--          Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
             v_comp_s                = v(k,j,i) - v_gtrans + v_stokes_zu(k)                        &
                                       * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
             flux_s(k)               = v_comp_s * (                                                &
                                         ( 37.0_wp * ibit5_s * adv_sca_5                           &
                                       +    7.0_wp * ibit4_s * adv_sca_3                           &
                                       +             ibit3_s * adv_sca_1                           &
                                         ) *                                                       &
                                       ( sk(k,j,i)   + sk(k,j-1,i)     )                           &
                                       - (  8.0_wp * ibit5_s * adv_sca_5                           &
                                       +             ibit4_s * adv_sca_3                           &
                                         ) *                                                       &
                                       ( sk(k,j+1,i) + sk(k,j-2,i)    )                            &
                                       + (           ibit5_s * adv_sca_5                           &
                                         ) *                                                       &
                                       ( sk(k,j+2,i) + sk(k,j-3,i)    )                            &
                                                  )

             diss_s(k)               = -ABS( v_comp_s ) * (                                        &
                                         ( 10.0_wp * ibit5_s * adv_sca_5                           &
                                       +    3.0_wp * ibit4_s * adv_sca_3                           &
                                       +             ibit3_s * adv_sca_1                           &
                                         ) *                                                       &
                                       ( sk(k,j,i)   - sk(k,j-1,i)    )                            &
                                       - (  5.0_wp * ibit5_s * adv_sca_5                           &
                                       +             ibit4_s * adv_sca_3                           &
                                         ) *                                                       &
                                       ( sk(k,j+1,i) - sk(k,j-2,i)    )                            &
                                       + (           ibit5_s * adv_sca_5                           &
                                          ) *                                                      &
                                       ( sk(k,j+2,i) - sk(k,j-3,i)    )                            &
                                                         )
#else
             flux_s(k) = swap_flux_y_local(k)
             diss_s(k) = swap_diss_y_local(k)
#endif
          ENDDO
!
!--       Now compute the fluxes and tendency terms for the horizontal and vertical parts above the
!--       top of the highest topography. No degradation for the horizontal parts, but for the
!--       vertical it is still needed.
          DO  k = nzb_max_l+1, nzt

             u_comp    = u(k,j,i+1) - u_gtrans + u_stokes_zu(k)
             flux_r(k) = u_comp * (                                                                &
                           37.0_wp * ( sk(k,j,i+1) + sk(k,j,i)   )                                 &
                         -  8.0_wp * ( sk(k,j,i+2) + sk(k,j,i-1) )                                 &
                         +           ( sk(k,j,i+3) + sk(k,j,i-2) ) ) * adv_sca_5
             diss_r(k) = -ABS( u_comp ) * (                                                        &
                           10.0_wp * ( sk(k,j,i+1) - sk(k,j,i)   )                                 &
                         -  5.0_wp * ( sk(k,j,i+2) - sk(k,j,i-1) )                                 &
                         +           ( sk(k,j,i+3) - sk(k,j,i-2) ) ) * adv_sca_5
#ifdef _OPENACC
!
!--          Recompute the left fluxes.
             u_comp_l  = u(k,j,i) - u_gtrans + u_stokes_zu(k)
             flux_l(k) = u_comp_l * (                                                              &
                           37.0_wp * ( sk(k,j,i)   + sk(k,j,i-1) )                                 &
                         -  8.0_wp * ( sk(k,j,i+1) + sk(k,j,i-2) )                                 &
                         +           ( sk(k,j,i+2) + sk(k,j,i-3) ) ) * adv_sca_5

             diss_l(k) = -ABS( u_comp_l ) * (                                                      &
                            10.0_wp * ( sk(k,j,i)   - sk(k,j,i-1) )                                &
                         -   5.0_wp * ( sk(k,j,i+1) - sk(k,j,i-2) )                                &
                         +            ( sk(k,j,i+2) - sk(k,j,i-3) ) ) * adv_sca_5
#else
             flux_l(k) = swap_flux_x_local(k,j)
             diss_l(k) = swap_diss_x_local(k,j)

#endif

             v_comp    = v(k,j+1,i) - v_gtrans + v_stokes_zu(k)
             flux_n(k) = v_comp * (                                                                &
                           37.0_wp * ( sk(k,j+1,i) + sk(k,j,i)   )                                 &
                         -  8.0_wp * ( sk(k,j+2,i) + sk(k,j-1,i) )                                 &
                         +           ( sk(k,j+3,i) + sk(k,j-2,i) ) ) * adv_sca_5
             diss_n(k) = -ABS( v_comp ) * (                                                        &
                           10.0_wp * ( sk(k,j+1,i) - sk(k,j,i)   )                                 &
                         -  5.0_wp * ( sk(k,j+2,i) - sk(k,j-1,i) )                                 &
                         +           ( sk(k,j+3,i) - sk(k,j-2,i) ) ) * adv_sca_5
#ifdef _OPENACC
!
!--          Recompute the south fluxes.
             v_comp_s  = v(k,j,i) - v_gtrans + v_stokes_zu(k)
             flux_s(k) = v_comp_s  * (                                                             &
                           37.0_wp * ( sk(k,j,i)   + sk(k,j-1,i) )                                 &
                         -  8.0_wp * ( sk(k,j+1,i) + sk(k,j-2,i) )                                 &
                         +           ( sk(k,j+2,i) + sk(k,j-3,i) ) ) * adv_sca_5
             diss_s(k) = -ABS( v_comp_s ) * (                                                      &
                           10.0_wp * ( sk(k,j,i)   - sk(k,j-1,i) )                                 &
                         -  5.0_wp * ( sk(k,j+1,i) - sk(k,j-2,i) )                                 &
                         +           ( sk(k,j+2,i) - sk(k,j-3,i) ) ) * adv_sca_5
#else
             flux_s(k) = swap_flux_y_local(k)
             diss_s(k) = swap_diss_y_local(k)
#endif
          ENDDO
!
!--       Now, compute vertical fluxes. Split loop into a part treating the lowest grid points with
!--       indirect indexing, a main loop without indirect indexing, and a loop for the uppermost
!--       grid points with indirect indexing. This allows better vectorization for the main loop.
!--       First, compute the flux at model surface, which need has to be calculated explicetely for
!--       the tendency at the first w-level. For topography wall this is done implicitely by
!--       advc_flag.
          k = nzb
          ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )
          flux_t(nzb) = w(k,j,i)         * rho_air_zw(k) * ( ibit6 * adv_sca_1 )                   &
                                                         * ( sk(k+1,j,i) + sk(k,j,i) )
          diss_t(nzb) = -ABS( w(k,j,i) ) * rho_air_zw(k) * ( ibit6 * adv_sca_1 )                   &
                                                         * ( sk(k+1,j,i) - sk(k,j,i) )

          DO  k = nzb+1, nzb+1
             ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )

             k_ppp = k + 3 * ibit8
             k_pp  = k + 2 * ( 1 - ibit6  )
             k_mm  = k - 2 * ibit8

             flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                                              &
                         (  37.0_wp * ibit8 * adv_sca_5                                            &
                         +   7.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)   + sk(k,j,i)    )       &
                         - ( 8.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( sk(k_pp,j,i)  + sk(k-1,j,i)  )       &
                         + (          ibit8 * adv_sca_5 ) * ( sk(k_ppp,j,i) + sk(k_mm,j,i) )       &
                                                    )

             diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                                      &
                         (  10.0_wp * ibit8 * adv_sca_5                                            &
                         +   3.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)   - sk(k,j,i)    )       &
                         - ( 5.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( sk(k_pp,j,i)  - sk(k-1,j,i)  )       &
                         + (          ibit8 * adv_sca_5 ) * ( sk(k_ppp,j,i) - sk(k_mm,j,i) )       &
                                                            )
          ENDDO

          DO  k = nzb+2, nzt-2
             ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )

             flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                                              &
                         (  37.0_wp * ibit8 * adv_sca_5                                            &
                         +   7.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i) + sk(k,j,i)   )          &
                         - ( 8.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( sk(k+2,j,i) + sk(k-1,j,i) )          &
                         + (          ibit8 * adv_sca_5 ) * ( sk(k+3,j,i) + sk(k-2,j,i) )          &
                                                    )

             diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                                      &
                         (  10.0_wp * ibit8 * adv_sca_5                                            &
                         +   3.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)   - sk(k,j,i)    )       &
                         - ( 5.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( sk(k+2,j,i)  - sk(k-1,j,i)  )        &
                         + (          ibit8 * adv_sca_5 ) * ( sk(k+3,j,i) - sk(k-2,j,i) )          &
                                                            )
          ENDDO

          DO  k = nzt-1, nzt-symmetry_flag
             ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )
!
!--          k index needs to be modified near bottom and top, else array subscripts will be
!--          exceeded. (k_ppp = k+3, k_mm = k) if 5th order scheme is employed, else
!--          (k_ppp = k, k_mm = k). Special treatment is required for k_pp. This must be limited
!--          at k = nzt, where either the first order scheme is employed (ibit6), or all bit flags
!--          are zero (in case there is a rigid lid defined at or near the domain top).
             k_ppp = k + 3 * ibit8
             k_pp  = k + 2 * ( 1 - ibit6  ) * ( ibit6 + ibit7 + ibit8 )
             k_mm  = k - 2 * ibit8

             flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                                              &
                         (  37.0_wp * ibit8 * adv_sca_5                                            &
                         +   7.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)  + sk(k,j,i)    )        &
                         - ( 8.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( sk(k_pp,j,i) + sk(k-1,j,i)  )        &
                         + (          ibit8 * adv_sca_5 ) * ( sk(k_ppp,j,i)+ sk(k_mm,j,i) )        &
                                                    )

             diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                                      &
                         (  10.0_wp * ibit8 * adv_sca_5                                            &
                         +   3.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( sk(k+1,j,i)   - sk(k,j,i)    )       &
                         - ( 5.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( sk(k_pp,j,i)  - sk(k-1,j,i)  )       &
                         + (          ibit8 * adv_sca_5 ) * ( sk(k_ppp,j,i) - sk(k_mm,j,i) )       &
                                                            )
          ENDDO

!
!--       Set resolved/turbulent flux at model top to zero (w-level). In case that a symmetric
!--       behavior between bottom and top shall be guaranteed (closed channel flow), the flux at nzt
!--       is also set to zero.
          IF ( symmetry_flag == 1 ) THEN
             flux_t(nzt) = 0.0_wp
             diss_t(nzt) = 0.0_wp
          ENDIF
          flux_t(nzt+1) = 0.0_wp
          diss_t(nzt+1) = 0.0_wp

          DO  k = nzb+1, nzb_max_l

             flux_d = flux_t(k-1)
             diss_d = diss_t(k-1)

             ibit2 = REAL( IBITS(advc_flag(k,j,i),2,1), KIND = wp )
             ibit1 = REAL( IBITS(advc_flag(k,j,i),1,1), KIND = wp )
             ibit0 = REAL( IBITS(advc_flag(k,j,i),0,1), KIND = wp )

             ibit5 = REAL( IBITS(advc_flag(k,j,i),5,1), KIND = wp )
             ibit4 = REAL( IBITS(advc_flag(k,j,i),4,1), KIND = wp )
             ibit3 = REAL( IBITS(advc_flag(k,j,i),3,1), KIND = wp )

             ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )
!
!--          Calculate the divergence of the velocity field. A respective correction is needed to
!--          overcome numerical instabilities introduced by an insufficient reduction of divergences
!--          near topography.
             div         =   ( u(k,j,i+1) * ( ibit0 + ibit1 + ibit2 )                              &
                             - u(k,j,i)   * (                                                      &
                             REAL( IBITS(advc_flag(k,j,i-1),0,1), KIND = wp )                      &
                           + REAL( IBITS(advc_flag(k,j,i-1),1,1), KIND = wp )                      &
                           + REAL( IBITS(advc_flag(k,j,i-1),2,1), KIND = wp )                      &
                                            )                                                      &
                             ) * ddx                                                               &
                           + ( v(k,j+1,i) * ( ibit3 + ibit4 + ibit5 )                              &
                             - v(k,j,i)   * (                                                      &
                             REAL( IBITS(advc_flag(k,j-1,i),3,1), KIND = wp )                      &
                           + REAL( IBITS(advc_flag(k,j-1,i),4,1), KIND = wp )                      &
                           + REAL( IBITS(advc_flag(k,j-1,i),5,1), KIND = wp )                      &
                                            )                                                      &
                             ) * ddy                                                               &
                           + ( w(k,j,i)   * rho_air_zw(k)   *                                      &
                                            ( ibit6 + ibit7 + ibit8 )                              &
                             - w(k-1,j,i) * rho_air_zw(k-1) *                                      &
                                            (                                                      &
                             REAL( IBITS(advc_flag(k-1,j,i),6,1), KIND = wp )                      &
                           + REAL( IBITS(advc_flag(k-1,j,i),7,1), KIND = wp )                      &
                           + REAL( IBITS(advc_flag(k-1,j,i),8,1), KIND = wp )                      &
                                            )                                                      &
                             ) * drho_air(k) * ddzw(k)

             tend(k,j,i) = tend(k,j,i) - (                                                         &
                             ( flux_r(k) + diss_r(k) -                                             &
                               flux_l(k) - diss_l(k) ) * ddx                                       &
                           + ( flux_n(k) + diss_n(k) -                                             &
                               flux_s(k) - diss_s(k) ) * ddy                                       &
                           + ( flux_t(k) + diss_t(k) -                                             &
                               flux_d    - diss_d    ) * drho_air(k) * ddzw(k)                     &
                                         ) + sk(k,j,i) * div

#ifndef _OPENACC
             swap_flux_y_local(k)   = flux_n(k)
             swap_diss_y_local(k)   = diss_n(k)
             swap_flux_x_local(k,j) = flux_r(k)
             swap_diss_x_local(k,j) = diss_r(k)
#endif

          ENDDO

          DO  k = nzb_max_l+1, nzt

             flux_d = flux_t(k-1)
             diss_d = diss_t(k-1)
!
!--          Calculate the divergence of the velocity field. A respective correction is needed to
!--          overcome numerical instabilities introduced by an insufficient reduction of divergences
!--          near topography.
             div         =   ( u(k,j,i+1) - u(k,j,i) ) * ddx                                       &
                           + ( v(k,j+1,i) - v(k,j,i) ) * ddy                                       &
                           + ( w(k,j,i)   * rho_air_zw(k)                                          &
                             - w(k-1,j,i) * rho_air_zw(k-1)                                        &
                                                     ) * drho_air(k) * ddzw(k)

             tend(k,j,i) = tend(k,j,i) - (                                                         &
                             ( flux_r(k) + diss_r(k) -                                             &
                               flux_l(k) - diss_l(k) ) * ddx                                       &
                           + ( flux_n(k) + diss_n(k) -                                             &
                               flux_s(k) - diss_s(k) ) * ddy                                       &
                           + ( flux_t(k) + diss_t(k) -                                             &
                               flux_d    - diss_d    ) * drho_air(k) * ddzw(k)                     &
                                         ) + sk(k,j,i) * div

#ifndef _OPENACC
             swap_flux_y_local(k)   = flux_n(k)
             swap_diss_y_local(k)   = diss_n(k)
             swap_flux_x_local(k,j) = flux_r(k)
             swap_diss_x_local(k,j) = diss_r(k)
#endif
          ENDDO

!
!--       Evaluation of statistics.
          DO  k = nzb+1, nzt
             SELECT CASE ( sk_num )

                CASE ( 1 )
                   !$ACC ATOMIC
                   sums_wspts_ws_l(k,tn) = sums_wspts_ws_l(k,tn) +                                 &
                                           ( flux_t(k)                                             &
                                               / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )       &
                                               * ( w(k,j,i) - hom(k,1,3,0)                 )       &
                                           + diss_t(k)                                             &
                                               / ( ABS(w(k,j,i)) + 1.0E-20_wp              )       &
                                               *   ABS(w(k,j,i)  - hom(k,1,3,0)            )       &
                                           ) * weight_substep(intermediate_timestep_count)
                CASE ( 2 )
                   !$ACC ATOMIC
                   sums_wssas_ws_l(k,tn) = sums_wssas_ws_l(k,tn) +                                 &
                                           ( flux_t(k)                                             &
                                               / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )       &
                                               * ( w(k,j,i) - hom(k,1,3,0)                 )       &
                                           + diss_t(k)                                             &
                                               / ( ABS(w(k,j,i)) + 1.0E-20_wp              )       &
                                               *   ABS(w(k,j,i)  - hom(k,1,3,0)            )       &
                                           ) * weight_substep(intermediate_timestep_count)
                CASE ( 3 )
                   !$ACC ATOMIC
                   sums_wsqs_ws_l(k,tn)  = sums_wsqs_ws_l(k,tn) +                                  &
                                           ( flux_t(k)                                             &
                                               / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )       &
                                               * ( w(k,j,i) - hom(k,1,3,0)                 )       &
                                           + diss_t(k)                                             &
                                               / ( ABS(w(k,j,i)) + 1.0E-20_wp              )       &
                                               *   ABS(w(k,j,i)  - hom(k,1,3,0)            )       &
                                           ) * weight_substep(intermediate_timestep_count)
                CASE ( 4 )
                   !$ACC ATOMIC
                   sums_wsqcs_ws_l(k,tn)  = sums_wsqcs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 5 )
                   !$ACC ATOMIC
                   sums_wsqrs_ws_l(k,tn)  = sums_wsqrs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 6 )
                   !$ACC ATOMIC
                   sums_wsncs_ws_l(k,tn)  = sums_wsncs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 7 )
                   !$ACC ATOMIC
                   sums_wsnrs_ws_l(k,tn)  = sums_wsnrs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 8 )
                   !$ACC ATOMIC
                   sums_wsss_ws_l(k,tn)  = sums_wsss_ws_l(k,tn) +                                  &
                                           ( flux_t(k)                                             &
                                               / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )       &
                                               * ( w(k,j,i) - hom(k,1,3,0)                 )       &
                                           + diss_t(k)                                             &
                                               / ( ABS(w(k,j,i)) + 1.0E-20_wp              )       &
                                               *   ABS(w(k,j,i) - hom(k,1,3,0)             )       &
                                           ) * weight_substep(intermediate_timestep_count)
                CASE ( 9 )
                   !$ACC ATOMIC
                   sums_salsa_ws_l(k,tn)  = sums_salsa_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 10 )
                   !$ACC ATOMIC
                   sums_wsnis_ws_l(k,tn)  = sums_wsnis_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 11 )
                   !$ACC ATOMIC
                   sums_wsqis_ws_l(k,tn)  = sums_wsqis_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 12 )
                   !$ACC ATOMIC
                   sums_wsngs_ws_l(k,tn)  = sums_wsngs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 13 )
                   !$ACC ATOMIC
                   sums_wsqgs_ws_l(k,tn)  = sums_wsqgs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)

                CASE ( 14 )
                   !$ACC ATOMIC
                   sums_wsnss_ws_l(k,tn)  = sums_wsnss_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 15 )
                   !$ACC ATOMIC
                   sums_wsqss_ws_l(k,tn)  = sums_wsqss_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)


             END SELECT

          ENDDO
       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(49), 'advec_s_ws', 'stop' )

 END SUBROUTINE advec_s_ws


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Scalar advection with advanced divergence correction. Call for all grid points.
!> The fifth-order scheme of Wicker and Skamarock is only conditionally convervative, meaning
!> scalar is only conserved as long as the flow divergence is negligible. In case of topography,
!> however, the remaining flow divergence can be still significant near walls, leading to numerical
!> errors in the advection discretization affecting the scalar conservation. To maintain scalar
!> conservation, an advanced divergence correction has been implemented according to the method
!> implemented in the Bott-Clond scheme. With this method, scalar is almost fully conserved.
!> Please note, this method has not been extensively tested yet. Further, please note, this method
!> alters the stability properties of the advection discretization. Tests showed that scalar
!> advection can become unstable for Courant numbers >= 0.5.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_s_ws_div_bc( advc_flag, sk, sk_char,                                             &
                               non_cyclic_l, non_cyclic_n, non_cyclic_r, non_cyclic_s,             &
                               advanced_divergence_correction )


    CHARACTER (LEN = *), INTENT(IN) ::  sk_char !< string identifier, used for assign fluxes
                                                !< to the correct dimension in the analysis array

    INTEGER(iwp) ::  i         !< grid index along x-direction
    INTEGER(iwp) ::  j         !< grid index along y-direction
    INTEGER(iwp) ::  k         !< grid index along z-direction
    INTEGER(iwp) ::  k_mm      !< k-2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_pp      !< k+2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_ppp     !< k+3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  nzb_max_l !< index indicating upper bound for order degradation of horizontal advection terms
    INTEGER(iwp) ::  sk_num    !< integer identifier, used for assign fluxes to the correct dimension in the analysis array
    INTEGER(iwp) ::  tn = 0    !< number of OpenMP thread (is always zero here)

    INTEGER(iwp), INTENT(IN), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::                          &
                                               advc_flag !< flag array to control order of scalar advection

    LOGICAL  ::  advanced_divergence_correction !< flag to trigger advanced flow-divergence correction
    LOGICAL  ::  non_cyclic_l                   !< flag that indicates non-cyclic boundary on the left
    LOGICAL  ::  non_cyclic_n                   !< flag that indicates non-cyclic boundary on the north
    LOGICAL  ::  non_cyclic_r                   !< flag that indicates non-cyclic boundary on the right
    LOGICAL  ::  non_cyclic_s                   !< flag that indicates non-cyclic boundary on the south

    REAL(wp) ::  d_new         !< component-wise divergence
    REAL(wp) ::  ibit0         !< flag indicating 1st-order scheme along x-direction
    REAL(wp) ::  ibit1         !< flag indicating 3rd-order scheme along x-direction
    REAL(wp) ::  ibit2         !< flag indicating 5th-order scheme along x-direction
    REAL(wp) ::  ibit3         !< flag indicating 1st-order scheme along y-direction
    REAL(wp) ::  ibit4         !< flag indicating 3rd-order scheme along y-direction
    REAL(wp) ::  ibit5         !< flag indicating 5th-order scheme along y-direction
    REAL(wp) ::  ibit6         !< flag indicating 1st-order scheme along z-direction
    REAL(wp) ::  ibit7         !< flag indicating 3rd-order scheme along z-direction
    REAL(wp) ::  ibit8         !< flag indicating 5th-order scheme along z-direction
    REAL(wp) ::  tendcy        !< tendency resulting from component-wise flux divergence
    REAL(wp) ::  u_comp        !< advection velocity along x-direction
    REAL(wp) ::  v_comp        !< advection velocity along y-direction


    REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_n     !< discretized artificial dissipation at northward-side
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_r     !< discretized artificial dissipation at rightward-side
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  diss_t     !< discretized artificial dissipation at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_n     !< discretized 6th-order flux at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_r     !< discretized 6th-order flux at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1)         ::  flux_t     !< discretized 6th-order flux at top of the grid box

    REAL(wp), DIMENSION(nzb+1:nzt)         ::  swap_diss_y_local !< discretized artificial dissipation at southward-side
    REAL(wp), DIMENSION(nzb+1:nzt)         ::  swap_flux_y_local !< discretized 6th-order flux at northward-side

    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_diss_x_local !< discretized artificial dissipation at leftward-side
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn) ::  swap_flux_x_local !< discretized 6th-order flux at leftward-side

    REAL(wp), DIMENSION(:,:,:), POINTER    ::  s1                !< pointer used to swap old and updated scalar values
    REAL(wp), DIMENSION(:,:,:), POINTER    ::  s2                !< pointer used to swap old and updated scalar values
!
!-- sk is an array from parameter list. It should not be a pointer, because in that case the
!-- compiler can not assume a stride 1 and cannot perform a strided one vector load. Adding the
!-- CONTIGUOUS keyword makes things even worse, because the compiler cannot assume strided one in
!-- the caller side.
    REAL(wp), DIMENSION(nzb+1:nzt,nysg:nyng,nxlg:nxrg)             ::  d      !< component-wise divergence
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  sk     !< advected scalar
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), TARGET     ::  sk_t1  !< target array used to swap old and updated scalar values
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), TARGET     ::  sk_t2  !< target array used to swap old and updated scalar values

    CALL cpu_log( log_point_s(49), 'advec_s_ws_div_bc', 'start' )

    SELECT CASE ( sk_char )

       CASE ( 'pt' )
          sk_num = 1
       CASE ( 'sa' )
          sk_num = 2
       CASE ( 'q' )
          sk_num = 3
       CASE ( 'qc' )
          sk_num = 4
       CASE ( 'qr' )
          sk_num = 5
       CASE ( 'nc' )
          sk_num = 6
       CASE ( 'nr' )
          sk_num = 7
       CASE ( 's' )
          sk_num = 8
       CASE ( 'aerosol_mass', 'aerosol_number', 'salsa_gas' )
          sk_num = 9
       CASE ( 'ni' )
          sk_num = 10
       CASE ( 'qi' )
          sk_num = 11
       CASE ( 'ng' )
          sk_num = 12
       CASE ( 'qg' )
          sk_num = 13
       CASE ( 'ns' )
          sk_num = 14
       CASE ( 'qs' )
          sk_num = 15
       CASE DEFAULT
          sk_num = 0

    END SELECT
!
!-- Store originally passed scalar array.
    s1 => sk_t1
    s2 => sk_t2
    s1 = sk
!
!-- To avoid compiler complaints.
    advanced_divergence_correction = advanced_divergence_correction

    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Used local modified copy of nzb_max (used to degrade order of discretization) at
!--       non-cyclic boundaries. Modify only at relevant points instead of the entire subdomain.
!--       This should lead to better load balance between boundary and non-boundary PEs.
          IF( non_cyclic_l  .AND.  i <= nxl + 2  .OR.                                              &
              non_cyclic_r  .AND.  i >= nxr - 2  .OR.                                              &
              non_cyclic_s  .AND.  j <= nys + 2  .OR.                                              &
              non_cyclic_n  .AND.  j >= nyn - 2 )  THEN
             nzb_max_l = nzt
          ELSE
             nzb_max_l = nzb_max
          END IF
!
!--       Compute leftside fluxes of the respective PE bounds.
          IF ( i == nxl )  THEN

             DO  k = nzb+1, nzb_max_l

                ibit2 = REAL( IBITS(advc_flag(k,j,i-1),2,1), KIND = wp )
                ibit1 = REAL( IBITS(advc_flag(k,j,i-1),1,1), KIND = wp )
                ibit0 = REAL( IBITS(advc_flag(k,j,i-1),0,1), KIND = wp )
!
!--             Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
                u_comp                    = u(k,j,i) - u_gtrans + u_stokes_zu(k)                   &
                                            * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
                swap_flux_x_local(k,j)    = u_comp * (                                             &
                                              ( 37.0_wp * ibit2 * adv_sca_5                        &
                                            +    7.0_wp * ibit1 * adv_sca_3                        &
                                            +             ibit0 * adv_sca_1                        &
                                              ) *                                                  &
                                              ( s1(k,j,i)   + s1(k,j,i-1)    )                     &
                                            - (  8.0_wp * ibit2 * adv_sca_5                        &
                                            +             ibit1 * adv_sca_3                        &
                                              ) *                                                  &
                                              ( s1(k,j,i+1) + s1(k,j,i-2)    )                     &
                                            + (           ibit2 * adv_sca_5                        &
                                              ) *                                                  &
                                              ( s1(k,j,i+2) + s1(k,j,i-3)    )                     &
                                                     )

                swap_diss_x_local(k,j)    = -ABS( u_comp ) * (                                     &
                                              ( 10.0_wp * ibit2 * adv_sca_5                        &
                                            +    3.0_wp * ibit1 * adv_sca_3                        &
                                            +             ibit0 * adv_sca_1                        &
                                              ) *                                                  &
                                              ( s1(k,j,i)   - s1(k,j,i-1)    )                     &
                                            - (  5.0_wp * ibit2 * adv_sca_5                        &
                                            +             ibit1 * adv_sca_3                        &
                                              ) *                                                  &
                                              ( s1(k,j,i+1) - s1(k,j,i-2)    )                     &
                                            + (           ibit2 * adv_sca_5                        &
                                              ) *                                                  &
                                              ( s1(k,j,i+2) - s1(k,j,i-3)    )                     &
                                                             )

             ENDDO

             DO  k = nzb_max_l+1, nzt

                u_comp                    = u(k,j,i) - u_gtrans + u_stokes_zu(k)
                swap_flux_x_local(k,j)    = u_comp * (                                             &
                                              37.0_wp * ( s1(k,j,i)   + s1(k,j,i-1) )              &
                                            -  8.0_wp * ( s1(k,j,i+1) + s1(k,j,i-2) )              &
                                            +           ( s1(k,j,i+2) + s1(k,j,i-3) )              &
                                                     ) * adv_sca_5

                swap_diss_x_local(k,j)    = -ABS( u_comp ) * (                                     &
                                              10.0_wp * ( s1(k,j,i)   - s1(k,j,i-1) )              &
                                            -  5.0_wp * ( s1(k,j,i+1) - s1(k,j,i-2) )              &
                                            +           ( s1(k,j,i+2) - s1(k,j,i-3) )              &
                                                             ) * adv_sca_5

             ENDDO

          ENDIF

!
!--       Now compute the fluxes and tendency terms for the horizontal and vertical parts up to the
!--       top of the highest topography.
          DO  k = nzb+1, nzb_max_l
!
!--          Note: It is faster to conduct all multiplications explicitly, e.g. * adv_sca_5 ... than
!--          to determine a factor and multiplicate the flux at the end.
             ibit2 = REAL( IBITS(advc_flag(k,j,i),2,1), KIND = wp )
             ibit1 = REAL( IBITS(advc_flag(k,j,i),1,1), KIND = wp )
             ibit0 = REAL( IBITS(advc_flag(k,j,i),0,1), KIND = wp )
!
!--          Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
             u_comp    = u(k,j,i+1) - u_gtrans + u_stokes_zu(k)                                    &
                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 1 ) )
             flux_r(k) = u_comp * (                                                                &
                         (  37.0_wp * ibit2 * adv_sca_5                                            &
                         +   7.0_wp * ibit1 * adv_sca_3                                            &
                         +            ibit0 * adv_sca_1 ) * ( s1(k,j,i+1) + s1(k,j,i)   )          &
                         - ( 8.0_wp * ibit2 * adv_sca_5                                            &
                         +            ibit1 * adv_sca_3 ) * ( s1(k,j,i+2) + s1(k,j,i-1) )          &
                         + (          ibit2 * adv_sca_5 ) * ( s1(k,j,i+3) + s1(k,j,i-2) )          &
                                  )

             diss_r(k) = -ABS( u_comp ) * (                                                        &
                           (  10.0_wp * ibit2 * adv_sca_5                                          &
                           +   3.0_wp * ibit1 * adv_sca_3                                          &
                           +            ibit0 * adv_sca_1 ) * ( s1(k,j,i+1) - s1(k,j,i)  )         &
                           - ( 5.0_wp * ibit2 * adv_sca_5                                          &
                           +            ibit1 * adv_sca_3 ) * ( s1(k,j,i+2) - s1(k,j,i-1) )        &
                           + (          ibit2 * adv_sca_5 ) * ( s1(k,j,i+3) - s1(k,j,i-2) )        &
                                          )
!
!--          Component-wise divergence correction. Therefore, compute temporay new solution using
!--          an Euler scheme with only advection along x (attention: conditionally unstable
!--          with WS-scheme!).
             d_new = d(k,j,i) - ( u(k,j,i+1) - u(k,j,i) ) * dt_3d * ddx

             tendcy = - ( flux_r(k)             + diss_r(k) - (                                    &
                          swap_flux_x_local(k,j) + swap_diss_x_local(k,j) ) ) * ddx

             s2(k,j,i) = ( ( 1.0_wp + d(k,j,i) ) * s1(k,j,i) + tendcy * dt_3d ) / ( 1.0_wp + d_new )

             d(k,j,i) = d_new

             swap_flux_x_local(k,j) = flux_r(k)
             swap_diss_x_local(k,j) = diss_r(k)

          ENDDO
!
!--       Now compute the fluxes and tendency terms for the horizontal and vertical parts above the
!--       top of the highest topography. No degradation for the horizontal parts, but for the
!--       vertical it is still needed.
          DO  k = nzb_max_l+1, nzt

             u_comp    = u(k,j,i+1) - u_gtrans + u_stokes_zu(k)
             flux_r(k) = u_comp * (                                                                &
                           37.0_wp * ( s1(k,j,i+1) + s1(k,j,i)   )                                 &
                         -  8.0_wp * ( s1(k,j,i+2) + s1(k,j,i-1) )                                 &
                         +           ( s1(k,j,i+3) + s1(k,j,i-2) ) ) * adv_sca_5
             diss_r(k) = -ABS( u_comp ) * (                                                        &
                           10.0_wp * ( s1(k,j,i+1) - s1(k,j,i)   )                                 &
                         -  5.0_wp * ( s1(k,j,i+2) - s1(k,j,i-1) )                                 &
                         +           ( s1(k,j,i+3) - s1(k,j,i-2) ) ) * adv_sca_5
!
!--          Component-wise divergence correction. Therefore, compute temporay new solution using
!--          an Euler scheme with only advection along x (attention: conditionally unstable
!--          with WS-scheme!).
             d_new = d(k,j,i) - ( u(k,j,i+1) - u(k,j,i) ) * dt_3d * ddx

             tendcy = - ( flux_r(k)             + diss_r(k) - (                                    &
                          swap_flux_x_local(k,j) + swap_diss_x_local(k,j) ) ) * ddx

             s2(k,j,i) = ( ( 1.0_wp + d(k,j,i) ) * s1(k,j,i) + tendcy * dt_3d ) / ( 1.0_wp + d_new )

             d(k,j,i) = d_new

             swap_flux_x_local(k,j) = flux_r(k)
             swap_diss_x_local(k,j) = diss_r(k)

          ENDDO
       ENDDO
    ENDDO

    CALL exchange_horiz( s2, nbgp )

    s1 => sk_t2
    s2 => sk_t1

    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Used local modified copy of nzb_max (used to degrade order of discretization) at
!--       non-cyclic boundaries. Modify only at relevant points instead of the entire subdomain.
!--       This should lead to better load balance between boundary and non-boundary PEs.
          IF( non_cyclic_l  .AND.  i <= nxl + 2  .OR.                                              &
              non_cyclic_r  .AND.  i >= nxr - 2  .OR.                                              &
              non_cyclic_s  .AND.  j <= nys + 2  .OR.                                              &
              non_cyclic_n  .AND.  j >= nyn - 2 )  THEN
             nzb_max_l = nzt
          ELSE
             nzb_max_l = nzb_max
          END IF
!
!--       Compute southside fluxes of the respective PE bounds.
          IF ( j == nys )  THEN
!
!--          Up to the top of the highest topography.
             DO  k = nzb+1, nzb_max_l

                ibit5 = REAL( IBITS(advc_flag(k,j-1,i),5,1), KIND = wp )
                ibit4 = REAL( IBITS(advc_flag(k,j-1,i),4,1), KIND = wp )
                ibit3 = REAL( IBITS(advc_flag(k,j-1,i),3,1), KIND = wp )
!
!--             Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
                v_comp                  = v(k,j,i) - v_gtrans + v_stokes_zu(k)                     &
                                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                swap_flux_y_local(k)    = v_comp *         (                                       &
                                             ( 37.0_wp * ibit5 * adv_sca_5                         &
                                          +     7.0_wp * ibit4 * adv_sca_3                         &
                                          +              ibit3 * adv_sca_1                         &
                                             ) *                                                   &
                                             ( s1(k,j,i)   + s1(k,j-1,i)    )                      &
                                          -  (  8.0_wp * ibit5 * adv_sca_5                         &
                                          +              ibit4 * adv_sca_3                         &
                                             ) *                                                   &
                                             ( s1(k,j+1,i) + s1(k,j-2,i)    )                      &
                                          +  (           ibit5 * adv_sca_5                         &
                                             ) *                                                   &
                                             ( s1(k,j+2,i) + s1(k,j-3,i)    )                      &
                                                           )

                swap_diss_y_local(k)    = -ABS( v_comp ) * (                                       &
                                            ( 10.0_wp * ibit5 * adv_sca_5                          &
                                          +    3.0_wp * ibit4 * adv_sca_3                          &
                                          +             ibit3 * adv_sca_1                          &
                                            ) *                                                    &
                                            ( s1(k,j,i)   - s1(k,j-1,i)  )                         &
                                          - (  5.0_wp * ibit5 * adv_sca_5                          &
                                          +             ibit4 * adv_sca_3                          &
                                            ) *                                                    &
                                            ( s1(k,j+1,i) - s1(k,j-2,i)  )                         &
                                          + (           ibit5 * adv_sca_5                          &
                                            ) *                                                    &
                                            ( s1(k,j+2,i) - s1(k,j-3,i)  )                         &
                                                           )

             ENDDO
!
!--          Above to the top of the highest topography. No degradation necessary.
             DO  k = nzb_max_l+1, nzt

                v_comp                  = v(k,j,i) - v_gtrans + v_stokes_zu(k)
                swap_flux_y_local(k)    = v_comp * (                                               &
                                            37.0_wp * ( s1(k,j,i)   + s1(k,j-1,i) )                &
                                          -  8.0_wp * ( s1(k,j+1,i) + s1(k,j-2,i) )                &
                                          +           ( s1(k,j+2,i) + s1(k,j-3,i) )                &
                                                   ) * adv_sca_5
                swap_diss_y_local(k)    = -ABS( v_comp ) * (                                       &
                                            10.0_wp * ( s1(k,j,i)   - s1(k,j-1,i) )                &
                                          -  5.0_wp * ( s1(k,j+1,i) - s1(k,j-2,i) )                &
                                          +             s1(k,j+2,i) - s1(k,j-3,i)                  &
                                                           ) * adv_sca_5

             ENDDO

          ENDIF
!
!--       Now compute the fluxes and tendency terms for the horizontal and vertical parts up to the
!--       top of the highest topography.
          DO  k = nzb+1, nzb_max_l
             ibit5 = REAL( IBITS(advc_flag(k,j,i),5,1), KIND = wp )
             ibit4 = REAL( IBITS(advc_flag(k,j,i),4,1), KIND = wp )
             ibit3 = REAL( IBITS(advc_flag(k,j,i),3,1), KIND = wp )
!
!--          Multiply Stokes velocity with topography flags to avoid possible fluxes at walls.
             v_comp    = v(k,j+1,i) - v_gtrans + v_stokes_zu(k)                                    &
                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j+1,i), 2 ) )
             flux_n(k) = v_comp * (                                                                &
                         (  37.0_wp * ibit5 * adv_sca_5                                            &
                         +   7.0_wp * ibit4 * adv_sca_3                                            &
                         +            ibit3 * adv_sca_1 ) * ( s1(k,j+1,i) + s1(k,j,i)   )          &
                         - ( 8.0_wp * ibit5 * adv_sca_5                                            &
                         +            ibit4 * adv_sca_3 ) * ( s1(k,j+2,i) + s1(k,j-1,i) )          &
                         + (          ibit5 * adv_sca_5 ) * ( s1(k,j+3,i) + s1(k,j-2,i) )          &
                                  )

             diss_n(k) = -ABS( v_comp ) * (                                                        &
                         (  10.0_wp * ibit5 * adv_sca_5                                            &
                         +   3.0_wp * ibit4 * adv_sca_3                                            &
                         +            ibit3 * adv_sca_1 ) * ( s1(k,j+1,i) - s1(k,j,i)   )          &
                         - ( 5.0_wp * ibit5 * adv_sca_5                                            &
                         +            ibit4 * adv_sca_3 ) * ( s1(k,j+2,i) - s1(k,j-1,i) )          &
                         + (          ibit5 * adv_sca_5 ) * ( s1(k,j+3,i) - s1(k,j-2,i) )          &
                                          )
!
!--          Component-wise divergence correction. Therefore, compute temporay new solution using
!--          an Euler scheme with only advection along y (attention: conditionally unstable
!--          with WS-scheme!).
             d_new = d(k,j,i) - ( v(k,j+1,i) - v(k,j,i) ) * dt_3d * ddy

             tendcy = - ( flux_n(k)            + diss_n(k) - (                                     &
                          swap_flux_y_local(k) + swap_diss_y_local(k) ) ) * ddy

             s2(k,j,i) = ( ( 1.0_wp + d(k,j,i) ) * s1(k,j,i) + tendcy * dt_3d ) / ( 1.0_wp + d_new )

             d(k,j,i) = d_new

             swap_flux_y_local(k)   = flux_n(k)
             swap_diss_y_local(k)   = diss_n(k)
          ENDDO

          DO  k = nzb_max_l+1, nzt

             v_comp    = v(k,j+1,i) - v_gtrans + v_stokes_zu(k)
             flux_n(k) = v_comp * (                                                                &
                           37.0_wp * ( s1(k,j+1,i) + s1(k,j,i)   )                                 &
                         -  8.0_wp * ( s1(k,j+2,i) + s1(k,j-1,i) )                                 &
                         +           ( s1(k,j+3,i) + s1(k,j-2,i) ) ) * adv_sca_5
             diss_n(k) = -ABS( v_comp ) * (                                                        &
                           10.0_wp * ( s1(k,j+1,i) - s1(k,j,i)   )                                 &
                         -  5.0_wp * ( s1(k,j+2,i) - s1(k,j-1,i) )                                 &
                         +           ( s1(k,j+3,i) - s1(k,j-2,i) ) ) * adv_sca_5
!
!--          Component-wise divergence correction. Therefore, compute temporay new solution using
!--          an Euler scheme with only advection along y (attention: conditionally unstable
!--          with WS-scheme!).
             d_new = d(k,j,i) - ( v(k,j+1,i) - v(k,j,i) ) * dt_3d * ddy

             tendcy = - ( flux_n(k)            + diss_n(k) - (                                     &
                          swap_flux_y_local(k) + swap_diss_y_local(k) ) ) * ddy

             s2(k,j,i) = ( ( 1.0_wp + d(k,j,i) ) * s1(k,j,i) + tendcy * dt_3d ) / ( 1.0_wp + d_new )

             d(k,j,i) = d_new

             swap_flux_y_local(k)   = flux_n(k)
             swap_diss_y_local(k)   = diss_n(k)
          ENDDO

       ENDDO
    ENDDO

    CALL exchange_horiz( s2, nbgp )

    s1 => sk_t1
    s2 => sk_t2

    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Now, compute vertical fluxes. Split loop into a part treating the lowest grid points with
!--       indirect indexing, a main loop without indirect indexing, and a loop for the uppermost
!--       grid points with indirect indexing. This allows better vectorization for the main loop.
!--       First, compute the flux at model surface, which need has to be calculated explicetely for
!--       the tendency at the first w-level. For topography wall this is done implicitely by
!--       advc_flag.
          k = nzb
          ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )
          flux_t(nzb) = w(k,j,i)         * rho_air_zw(k) * ( ibit6 * adv_sca_1 )                   &
                                                         * ( s1(k+1,j,i) + s1(k,j,i) )
          diss_t(nzb) = -ABS( w(k,j,i) ) * rho_air_zw(k) * ( ibit6 * adv_sca_1 )                   &
                                                         * ( s1(k+1,j,i) - s1(k,j,i) )

          DO  k = nzb+1, nzb+1
             ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )

             k_ppp = k + 3 * ibit8
             k_pp  = k + 2 * ( 1 - ibit6  )
             k_mm  = k - 2 * ibit8

             flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                                              &
                         (  37.0_wp * ibit8 * adv_sca_5                                            &
                         +   7.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( s1(k+1,j,i)   + s1(k,j,i)    )       &
                         - ( 8.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( s1(k_pp,j,i)  + s1(k-1,j,i)  )       &
                         + (          ibit8 * adv_sca_5 ) * ( s1(k_ppp,j,i) + s1(k_mm,j,i) )       &
                                                    )

             diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                                      &
                         (  10.0_wp * ibit8 * adv_sca_5                                            &
                         +   3.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( s1(k+1,j,i)   - s1(k,j,i)    )       &
                         - ( 5.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( s1(k_pp,j,i)  - s1(k-1,j,i)  )       &
                         + (          ibit8 * adv_sca_5 ) * ( s1(k_ppp,j,i) - s1(k_mm,j,i) )       &
                                                            )
          ENDDO

          DO  k = nzb+2, nzt-2
             ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )

             flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                                              &
                         (  37.0_wp * ibit8 * adv_sca_5                                            &
                         +   7.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( s1(k+1,j,i) + s1(k,j,i)   )          &
                         - ( 8.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( s1(k+2,j,i) + s1(k-1,j,i) )          &
                         + (          ibit8 * adv_sca_5 ) * ( s1(k+3,j,i) + s1(k-2,j,i) )          &
                                                    )

             diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                                      &
                         (  10.0_wp * ibit8 * adv_sca_5                                            &
                         +   3.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( s1(k+1,j,i) - s1(k,j,i)    )         &
                         - ( 5.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( s1(k+2,j,i) - s1(k-1,j,i)  )         &
                         + (          ibit8 * adv_sca_5 ) * ( s1(k+3,j,i) - s1(k-2,j,i) )          &
                                                            )
          ENDDO

          DO  k = nzt-1, nzt-symmetry_flag
             ibit8 = REAL( IBITS(advc_flag(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flag(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flag(k,j,i),6,1), KIND = wp )
!
!--          k index needs to be modified near bottom and top, else array subscripts will be
!--          exceeded. (k_ppp = k+3, k_mm = k) if 5th order scheme is employed, else
!--          (k_ppp = k, k_mm = k). Special treatment is required for k_pp. This must be limited
!--          at k = nzt, where either the first order scheme is employed (ibit6), or all bit flags
!--          are zero (in case there is a rigid lid defined at or near the domain top).
             k_ppp = k + 3 * ibit8
             k_pp  = k + 2 * ( 1 - ibit6  ) * ( ibit6 + ibit7 + ibit8 )
             k_mm  = k - 2 * ibit8

             flux_t(k) = w(k,j,i) * rho_air_zw(k) * (                                              &
                         (  37.0_wp * ibit8 * adv_sca_5                                            &
                         +   7.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( s1(k+1,j,i)  + s1(k,j,i)    )        &
                         - ( 8.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( s1(k_pp,j,i) + s1(k-1,j,i)  )        &
                         + (          ibit8 * adv_sca_5 ) * ( s1(k_ppp,j,i)+ s1(k_mm,j,i) )        &
                                                    )

             diss_t(k) = -ABS( w(k,j,i) ) * rho_air_zw(k) * (                                      &
                         (  10.0_wp * ibit8 * adv_sca_5                                            &
                         +   3.0_wp * ibit7 * adv_sca_3                                            &
                         +            ibit6 * adv_sca_1 ) * ( s1(k+1,j,i)   - s1(k,j,i)    )       &
                         - ( 5.0_wp * ibit8 * adv_sca_5                                            &
                         +            ibit7 * adv_sca_3 ) * ( s1(k_pp,j,i)  - s1(k-1,j,i)  )       &
                         + (          ibit8 * adv_sca_5 ) * ( s1(k_ppp,j,i) - s1(k_mm,j,i) )       &
                                                            )
          ENDDO
!
!--       Set resolved/turbulent flux at model top to zero (w-level). In case that a symmetric
!--       behavior between bottom and top shall be guaranteed (closed channel flow), the flux at nzt
!--       is also set to zero.
          IF ( symmetry_flag == 1 ) THEN
             flux_t(nzt) = 0.0_wp
             diss_t(nzt) = 0.0_wp
          ENDIF
          flux_t(nzt+1) = 0.0_wp
          diss_t(nzt+1) = 0.0_wp
!
!--       Component-wise divergence correction. Therefore, compute temporay new solution using
!--       an Euler scheme with only advection along y (attention: conditionally unstable
!--       with WS-scheme!).
          DO  k = nzb+1, nzt
             d_new = d(k,j,i) - ( w(k,j,i) * rho_air_zw(k) - w(k-1,j,i) * rho_air_zw(k-1) )         &
                              * drho_air(k) * ddzw(k) * dt_3d

             tendcy = -( flux_t(k)   + diss_t(k) - (                                                &
                         flux_t(k-1) + diss_t(k-1) ) ) * ddzw(k) * drho_air(k)

             s2(k,j,i) = ( ( 1.0_wp + d(k,j,i) ) * s1(k,j,i) + tendcy * dt_3d ) / ( 1.0_wp + d_new )
             d(k,j,i) = d_new
          ENDDO
!
!--       Evaluation of statistics.
          DO  k = nzb+1, nzt
             SELECT CASE ( sk_num )

                CASE ( 1 )
                   sums_wspts_ws_l(k,tn) = sums_wspts_ws_l(k,tn) +                                 &
                                           ( flux_t(k)                                             &
                                               / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )       &
                                               * ( w(k,j,i) - hom(k,1,3,0)                 )       &
                                           + diss_t(k)                                             &
                                               / ( ABS(w(k,j,i)) + 1.0E-20_wp              )       &
                                               *   ABS(w(k,j,i)  - hom(k,1,3,0)            )       &
                                           ) * weight_substep(intermediate_timestep_count)
                CASE ( 2 )
                   sums_wssas_ws_l(k,tn) = sums_wssas_ws_l(k,tn) +                                 &
                                           ( flux_t(k)                                             &
                                               / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )       &
                                               * ( w(k,j,i) - hom(k,1,3,0)                 )       &
                                           + diss_t(k)                                             &
                                               / ( ABS(w(k,j,i)) + 1.0E-20_wp              )       &
                                               *   ABS(w(k,j,i)  - hom(k,1,3,0)            )       &
                                           ) * weight_substep(intermediate_timestep_count)
                CASE ( 3 )
                   sums_wsqs_ws_l(k,tn)  = sums_wsqs_ws_l(k,tn) +                                  &
                                           ( flux_t(k)                                             &
                                               / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )       &
                                               * ( w(k,j,i) - hom(k,1,3,0)                 )       &
                                           + diss_t(k)                                             &
                                               / ( ABS(w(k,j,i)) + 1.0E-20_wp              )       &
                                               *   ABS(w(k,j,i)  - hom(k,1,3,0)            )       &
                                           ) * weight_substep(intermediate_timestep_count)
                CASE ( 4 )
                   sums_wsqcs_ws_l(k,tn)  = sums_wsqcs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 5 )
                   sums_wsqrs_ws_l(k,tn)  = sums_wsqrs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 6 )
                   sums_wsncs_ws_l(k,tn)  = sums_wsncs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 7 )
                   sums_wsnrs_ws_l(k,tn)  = sums_wsnrs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 8 )
                   sums_wsss_ws_l(k,tn)  = sums_wsss_ws_l(k,tn) +                                  &
                                           ( flux_t(k)                                             &
                                               / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )       &
                                               * ( w(k,j,i) - hom(k,1,3,0)                 )       &
                                           + diss_t(k)                                             &
                                               / ( ABS(w(k,j,i)) + 1.0E-20_wp              )       &
                                               *   ABS(w(k,j,i) - hom(k,1,3,0)             )       &
                                           ) * weight_substep(intermediate_timestep_count)
                CASE ( 9 )
                   sums_salsa_ws_l(k,tn)  = sums_salsa_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i)  - hom(k,1,3,0)            )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 10 )
                   sums_wsnis_ws_l(k,tn)  = sums_wsnis_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 11 )
                   sums_wsqis_ws_l(k,tn)  = sums_wsqis_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 12 )
                   sums_wsngs_ws_l(k,tn)  = sums_wsngs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 13 )
                   sums_wsqgs_ws_l(k,tn)  = sums_wsqgs_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)

                CASE ( 14 )
                   sums_wsnss_ws_l(k,tn)  = sums_wsnss_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)
                CASE ( 15 )
                   sums_wsqss_ws_l(k,tn)  = sums_wsqss_ws_l(k,tn) +                                &
                                            ( flux_t(k)                                            &
                                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )      &
                                                * ( w(k,j,i) - hom(k,1,3,0)                 )      &
                                            + diss_t(k)                                            &
                                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )      &
                                                *   ABS(w(k,j,i) - hom(k,1,3,0)             )      &
                                            ) * weight_substep(intermediate_timestep_count)

             END SELECT

          ENDDO

       ENDDO
    ENDDO
!
!-- Re-calculate tendency term which is further used within the RK time-integration.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
              tend(k,j,i) = ( s2(k,j,i) - sk(k,j,i) ) / dt_3d
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(49), 'advec_s_ws_div_bc', 'stop' )

 END SUBROUTINE advec_s_ws_div_bc

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of u - Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_u_ws


    INTEGER(iwp) ::  i         !< grid index along x-direction
    INTEGER(iwp) ::  j         !< grid index along y-direction
    INTEGER(iwp) ::  k         !< grid index along z-direction
    INTEGER(iwp) ::  k_mm      !< k-2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_pp      !< k+2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_ppp     !< k+3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  nzb_max_l !< index indicating upper bound for order degradation of horizontal advection terms
    INTEGER(iwp) ::  tn = 0    !< number of OpenMP thread (is always zero here)

    REAL(wp)    ::  ibit0 !< flag indicating 1st-order scheme along x-direction
    REAL(wp)    ::  ibit1 !< flag indicating 3rd-order scheme along x-direction
    REAL(wp)    ::  ibit2 !< flag indicating 5th-order scheme along x-direction
#ifdef _OPENACC
    REAL(wp)    ::  ibit0_l !< flag indicating 1st-order scheme along x-direction
    REAL(wp)    ::  ibit1_l !< flag indicating 3rd-order scheme along x-direction
    REAL(wp)    ::  ibit2_l !< flag indicating 5th-order scheme along x-direction
#endif
    REAL(wp)    ::  ibit3 !< flag indicating 1st-order scheme along y-direction
    REAL(wp)    ::  ibit4 !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)    ::  ibit5 !< flag indicating 5th-order scheme along y-direction
#ifdef _OPENACC
    REAL(wp)    ::  ibit3_s !< flag indicating 1st-order scheme along y-direction
    REAL(wp)    ::  ibit4_s !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)    ::  ibit5_s !< flag indicating 5th-order scheme along y-direction
#endif
    REAL(wp)    ::  diss_d !< artificial dissipation term at grid box bottom
    REAL(wp)    ::  div    !< diverence on u-grid
    REAL(wp)    ::  flux_d !< 6th-order flux at grid box bottom
    REAL(wp)    ::  gu     !< Galilei-transformation velocity along x
    REAL(wp)    ::  gv     !< Galilei-transformation velocity along y
    REAL(wp)    ::  ibit6  !< flag indicating 1st-order scheme along z-direction
    REAL(wp)    ::  ibit7  !< flag indicating 3rd-order scheme along z-direction
    REAL(wp)    ::  ibit8  !< flag indicating 5th-order scheme along z-direction
    REAL(wp)    ::  u_comp_l !<
    REAL(wp)    ::  v_comp_s !< advection velocity along y

    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_n !< discretized artificial dissipation at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_r !< discretized artificial dissipation at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_t !< discretized artificial dissipation at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_n !< discretized 6th-order flux at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_r !< discretized 6th-order flux at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_t !< discretized 6th-order flux at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  u_comp !< advection velocity along x
    REAL(wp), DIMENSION(nzb:nzt+1) ::  v_comp !< advection velocity along y
    REAL(wp), DIMENSION(nzb:nzt+1) ::  w_comp !< advection velocity along z

    CALL cpu_log( log_point_s(68), 'advec_u_ws', 'start' )

    gu = 2.0_wp * u_gtrans
    gv = 2.0_wp * v_gtrans

    !$ACC PARALLEL LOOP COLLAPSE(2) FIRSTPRIVATE(tn, gu, gv) &
    !$ACC PRIVATE(i, j, k, k_mm, k_pp, k_ppp) &
    !$ACC PRIVATE(ibit0, ibit1, ibit2, ibit3, ibit4, ibit5) &
    !$ACC PRIVATE(ibit0_l, ibit1_l, ibit2_l) &
    !$ACC PRIVATE(ibit3_s, ibit4_s, ibit5_s) &
    !$ACC PRIVATE(nzb_max_l) &
    !$ACC PRIVATE(ibit6, ibit7, ibit8) &
    !$ACC PRIVATE(flux_r, diss_r) &
    !$ACC PRIVATE(flux_n, diss_n) &
    !$ACC PRIVATE(flux_t, diss_t, flux_d, diss_d) &
    !$ACC PRIVATE(flux_l_u, diss_l_u, flux_s_u, diss_s_u) &
    !$ACC PRIVATE(div, u_comp, u_comp_l, v_comp, v_comp_s, w_comp) &
    !$ACC PRESENT(advc_flags_m) &
    !$ACC PRESENT(u, v, w) &
    !$ACC PRESENT(drho_air, rho_air_zw, ddzw) &
    !$ACC PRESENT(tend) &
    !$ACC PRESENT(hom(:,1,1:3,0)) &
    !$ACC PRESENT(weight_substep(intermediate_timestep_count)) &
    !$ACC PRESENT(sums_us2_ws_l, sums_wsus_ws_l)
    DO  i = nxlu, nxr

       DO  j = nys, nyn
!
!--       Used local modified copy of nzb_max (used to degrade order of discretization) at
!--       non-cyclic boundaries. Modify only at relevant points instead of the entire subdomain.
!--       This should lead to better load balance between boundary and non-boundary PEs.
          IF( ( bc_dirichlet_l  .OR.  bc_radiation_l )  .AND.  i <= nxl + 2  .OR.                  &
              ( bc_dirichlet_r  .OR.  bc_radiation_r )  .AND.  i >= nxr - 2  .OR.                  &
              ( bc_dirichlet_s  .OR.  bc_radiation_s )  .AND.  j <= nys + 2  .OR.                  &
              ( bc_dirichlet_n  .OR.  bc_radiation_n )  .AND.  j >= nyn - 2 )  THEN
             nzb_max_l = nzt
          ELSE
             nzb_max_l = nzb_max
          END IF
#ifndef _OPENACC
!
!--       Compute southside fluxes for the respective boundary of PE
          IF ( j == nys  )  THEN

             DO  k = nzb+1, nzb_max_l

                ibit5 = REAL( IBITS(advc_flags_m(k,j-1,i),5,1), KIND = wp )
                ibit4 = REAL( IBITS(advc_flags_m(k,j-1,i),4,1), KIND = wp )
                ibit3 = REAL( IBITS(advc_flags_m(k,j-1,i),3,1), KIND = wp )

                v_comp_s       = v(k,j,i) + v(k,j,i-1) - gv
                flux_s_u(k,tn) = v_comp_s * (                                                      &
                               (  37.0_wp * ibit5 * adv_mom_5                                      &
                               +   7.0_wp * ibit4 * adv_mom_3                                      &
                               +            ibit3 * adv_mom_1 ) * ( u(k,j,i)   + u(k,j-1,i) )      &
                               - ( 8.0_wp * ibit5 * adv_mom_5                                      &
                               +            ibit4 * adv_mom_3 ) * ( u(k,j+1,i) + u(k,j-2,i) )      &
                               + (          ibit5 * adv_mom_5 ) * ( u(k,j+2,i) + u(k,j-3,i) )      &
                                            )

                diss_s_u(k,tn) = -ABS ( v_comp_s ) * (                                             &
                               (  10.0_wp * ibit5 * adv_mom_5                                      &
                               +   3.0_wp * ibit4 * adv_mom_3                                      &
                               +            ibit3 * adv_mom_1 ) * ( u(k,j,i)   - u(k,j-1,i) )      &
                               - ( 5.0_wp * ibit5 * adv_mom_5                                      &
                               +            ibit4 * adv_mom_3 ) * ( u(k,j+1,i) - u(k,j-2,i) )      &
                               + (          ibit5 * adv_mom_5 ) * ( u(k,j+2,i) - u(k,j-3,i) )      &
                                                     )

             ENDDO

             DO  k = nzb_max_l+1, nzt

                v_comp_s       = v(k,j,i) + v(k,j,i-1) - gv
                flux_s_u(k,tn) = v_comp_s * (                                                      &
                                 37.0_wp * ( u(k,j,i)   + u(k,j-1,i) )                             &
                               -  8.0_wp * ( u(k,j+1,i) + u(k,j-2,i) )                             &
                               +           ( u(k,j+2,i) + u(k,j-3,i) ) ) * adv_mom_5
                diss_s_u(k,tn) = -ABS( v_comp_s ) * (                                              &
                                 10.0_wp * ( u(k,j,i)   - u(k,j-1,i) )                             &
                               -  5.0_wp * ( u(k,j+1,i) - u(k,j-2,i) )                             &
                               +           ( u(k,j+2,i) - u(k,j-3,i) ) ) * adv_mom_5

             ENDDO

          ENDIF
!
!--       Compute leftside fluxes for the respective boundary of PE
          IF ( i == nxlu )  THEN

             DO  k = nzb+1, nzb_max_l

                ibit2 = REAL( IBITS(advc_flags_m(k,j,i-1),2,1), KIND = wp )
                ibit1 = REAL( IBITS(advc_flags_m(k,j,i-1),1,1), KIND = wp )
                ibit0 = REAL( IBITS(advc_flags_m(k,j,i-1),0,1), KIND = wp )

                u_comp_l         = u(k,j,i) + u(k,j,i-1) - gu
                flux_l_u(k,j,tn) = u_comp_l * (                                                    &
                                 (  37.0_wp * ibit2 * adv_mom_5                                    &
                                 +   7.0_wp * ibit1 * adv_mom_3                                    &
                                 +            ibit0 * adv_mom_1 ) * ( u(k,j,i)   + u(k,j,i-1) )    &
                                 - ( 8.0_wp * ibit2 * adv_mom_5                                    &
                                 +            ibit1 * adv_mom_3 ) * ( u(k,j,i+1) + u(k,j,i-2) )    &
                                 + (          ibit2 * adv_mom_5 ) * ( u(k,j,i+2) + u(k,j,i-3) )    &
                                              )

                diss_l_u(k,j,tn) = -ABS( u_comp_l ) * (                                            &
                                 (  10.0_wp * ibit2 * adv_mom_5                                    &
                                 +   3.0_wp * ibit1 * adv_mom_3                                    &
                                 +           ibit0  * adv_mom_1 ) * ( u(k,j,i)   - u(k,j,i-1) )    &
                                 - ( 5.0_wp * ibit2 * adv_mom_5                                    &
                                 +            ibit1 * adv_mom_3 ) * ( u(k,j,i+1) - u(k,j,i-2) )    &
                                 + (          ibit2 * adv_mom_5 ) * ( u(k,j,i+2) - u(k,j,i-3) )    &
                                                      )

             ENDDO

             DO  k = nzb_max_l+1, nzt

                u_comp_l         = u(k,j,i) + u(k,j,i-1) - gu
                flux_l_u(k,j,tn) = u_comp_l * (                                                    &
                                     37.0_wp * ( u(k,j,i)   + u(k,j,i-1) )                         &
                                   -  8.0_wp * ( u(k,j,i+1) + u(k,j,i-2) )                         &
                                   +           ( u(k,j,i+2) + u(k,j,i-3) ) ) * adv_mom_5
                diss_l_u(k,j,tn) = -ABS(u_comp_l) * (                                              &
                                     10.0_wp * ( u(k,j,i)   - u(k,j,i-1) )                         &
                                   -  5.0_wp * ( u(k,j,i+1) - u(k,j,i-2) )                         &
                                   +           ( u(k,j,i+2) - u(k,j,i-3) ) ) * adv_mom_5

             ENDDO

          ENDIF
#endif

!
!--       Now compute the fluxes for the horizontal and parts.
          DO  k = nzb+1, nzb_max_l

             ibit2 = REAL( IBITS(advc_flags_m(k,j,i),2,1), KIND = wp )
             ibit1 = REAL( IBITS(advc_flags_m(k,j,i),1,1), KIND = wp )
             ibit0 = REAL( IBITS(advc_flags_m(k,j,i),0,1), KIND = wp )

             u_comp(k) = u(k,j,i+1) + u(k,j,i)
             flux_r(k) = ( u_comp(k) - gu ) * (                                                    &
                         (  37.0_wp * ibit2 * adv_mom_5                                            &
                         +   7.0_wp * ibit1 * adv_mom_3                                            &
                         +            ibit0 * adv_mom_1 ) * ( u(k,j,i+1) + u(k,j,i)   )            &
                         - ( 8.0_wp * ibit2 * adv_mom_5                                            &
                         +            ibit1 * adv_mom_3 ) * ( u(k,j,i+2) + u(k,j,i-1) )            &
                         + (          ibit2 * adv_mom_5 ) * ( u(k,j,i+3) + u(k,j,i-2) )            &
                                              )

             diss_r(k) = -ABS( u_comp(k) - gu ) * (                                                &
                         (  10.0_wp * ibit2 * adv_mom_5                                            &
                         +   3.0_wp * ibit1 * adv_mom_3                                            &
                         +            ibit0 * adv_mom_1 ) * ( u(k,j,i+1) - u(k,j,i)   )            &
                         - ( 5.0_wp * ibit2 * adv_mom_5                                            &
                         +            ibit1 * adv_mom_3 ) * ( u(k,j,i+2) - u(k,j,i-1) )            &
                         + (          ibit2 * adv_mom_5 ) * ( u(k,j,i+3) - u(k,j,i-2) )            &
                                                  )

#ifdef _OPENACC
!
!--          Recompute the left fluxes.
             ibit2_l = REAL( IBITS(advc_flags_m(k,j,i-1),2,1), KIND = wp )
             ibit1_l = REAL( IBITS(advc_flags_m(k,j,i-1),1,1), KIND = wp )
             ibit0_l = REAL( IBITS(advc_flags_m(k,j,i-1),0,1), KIND = wp )

             u_comp_l          = u(k,j,i) + u(k,j,i-1) - gu
             flux_l_u(k,j,tn)  = u_comp_l * (                                                      &
                                 (  37.0_wp * ibit2_l * adv_mom_5                                  &
                                 +   7.0_wp * ibit1_l * adv_mom_3                                  &
                                 +            ibit0_l * adv_mom_1 ) * ( u(k,j,i)   + u(k,j,i-1) )  &
                                 - ( 8.0_wp * ibit2_l * adv_mom_5                                  &
                                 +            ibit1_l * adv_mom_3 ) * ( u(k,j,i+1) + u(k,j,i-2) )  &
                                 + (          ibit2_l * adv_mom_5 ) * ( u(k,j,i+2) + u(k,j,i-3) )  &
                                            )

             diss_l_u(k,j,tn)  = -ABS( u_comp_l ) * (                                              &
                                 (  10.0_wp * ibit2_l * adv_mom_5                                  &
                                 +   3.0_wp * ibit1_l * adv_mom_3                                  &
                                 +           ibit0_l  * adv_mom_1 ) * ( u(k,j,i)   - u(k,j,i-1) )  &
                                 - ( 5.0_wp * ibit2_l * adv_mom_5                                  &
                                 +            ibit1_l * adv_mom_3 ) * ( u(k,j,i+1) - u(k,j,i-2) )  &
                                 + (          ibit2_l * adv_mom_5 ) * ( u(k,j,i+2) - u(k,j,i-3) )  &
                                                    )
#endif


             ibit5 = REAL( IBITS(advc_flags_m(k,j,i),5,1), KIND = wp )
             ibit4 = REAL( IBITS(advc_flags_m(k,j,i),4,1), KIND = wp )
             ibit3 = REAL( IBITS(advc_flags_m(k,j,i),3,1), KIND = wp )

             v_comp(k) = v(k,j+1,i) + v(k,j+1,i-1) - gv
             flux_n(k) = v_comp(k) * (                                                             &
                         (  37.0_wp * ibit5 * adv_mom_5                                            &
                         +   7.0_wp * ibit4 * adv_mom_3                                            &
                         +            ibit3 * adv_mom_1 ) * ( u(k,j+1,i) + u(k,j,i)   )            &
                         - ( 8.0_wp * ibit5 * adv_mom_5                                            &
                         +            ibit4 * adv_mom_3 ) * ( u(k,j+2,i) + u(k,j-1,i) )            &
                         + (          ibit5 * adv_mom_5 ) * ( u(k,j+3,i) + u(k,j-2,i) )            &
                                     )

             diss_n(k) = -ABS ( v_comp(k) ) * (                                                    &
                         (  10.0_wp * ibit5 * adv_mom_5                                            &
                         +   3.0_wp * ibit4 * adv_mom_3                                            &
                         +            ibit3 * adv_mom_1 ) * ( u(k,j+1,i) - u(k,j,i)   )            &
                         - ( 5.0_wp * ibit5 * adv_mom_5                                            &
                         +            ibit4 * adv_mom_3 ) * ( u(k,j+2,i) - u(k,j-1,i) )            &
                         + (          ibit5 * adv_mom_5 ) * ( u(k,j+3,i) - u(k,j-2,i) )            &
                                              )

#ifdef _OPENACC
!
!--          Recompute the south fluxes.
             ibit5_s = REAL( IBITS(advc_flags_m(k,j-1,i),5,1), KIND = wp )
             ibit4_s = REAL( IBITS(advc_flags_m(k,j-1,i),4,1), KIND = wp )
             ibit3_s = REAL( IBITS(advc_flags_m(k,j-1,i),3,1), KIND = wp )

             v_comp_s         = v(k,j,i) + v(k,j,i-1) - gv
             flux_s_u(k,tn)   = v_comp_s * (                                                       &
                                (  37.0_wp * ibit5_s * adv_mom_5                                   &
                                +   7.0_wp * ibit4_s * adv_mom_3                                   &
                                +            ibit3_s * adv_mom_1 ) * ( u(k,j,i)   + u(k,j-1,i) )   &
                                - ( 8.0_wp * ibit5_s * adv_mom_5                                   &
                                +            ibit4_s * adv_mom_3 ) * ( u(k,j+1,i) + u(k,j-2,i) )   &
                                + (          ibit5_s * adv_mom_5 ) * ( u(k,j+2,i) + u(k,j-3,i) )   &
                                           )

             diss_s_u(k,tn)   = -ABS ( v_comp_s ) * (                                              &
                                (  10.0_wp * ibit5_s * adv_mom_5                                   &
                                +   3.0_wp * ibit4_s * adv_mom_3                                   &
                                +            ibit3_s * adv_mom_1 ) * ( u(k,j,i)   - u(k,j-1,i) )   &
                                - ( 5.0_wp * ibit5_s * adv_mom_5                                   &
                                +            ibit4_s * adv_mom_3 ) * ( u(k,j+1,i) - u(k,j-2,i) )   &
                                + (          ibit5_s * adv_mom_5 ) * ( u(k,j+2,i) - u(k,j-3,i) )   &
                                                    )
#endif
          ENDDO

          DO  k = nzb_max_l+1, nzt

             u_comp(k) = u(k,j,i+1) + u(k,j,i)
             flux_r(k) = ( u_comp(k) - gu ) * (                                                    &
                            37.0_wp * ( u(k,j,i+1) + u(k,j,i)   )                                  &
                         -   8.0_wp * ( u(k,j,i+2) + u(k,j,i-1) )                                  &
                         +           ( u(k,j,i+3) + u(k,j,i-2) ) ) * adv_mom_5
             diss_r(k) = -ABS( u_comp(k) - gu ) * (                                                &
                            10.0_wp * ( u(k,j,i+1) - u(k,j,i)   )                                  &
                          -  5.0_wp * ( u(k,j,i+2) - u(k,j,i-1) )                                  &
                          +           ( u(k,j,i+3) - u(k,j,i-2) ) ) * adv_mom_5
#ifdef _OPENACC
!
!--          Recompute the left fluxes.
             u_comp_l           = u(k,j,i) + u(k,j,i-1) - gu
             flux_l_u(k,j,tn)   = u_comp_l * (                                                     &
                                    37.0_wp * ( u(k,j,i) + u(k,j,i-1)   )                          &
                                  -  8.0_wp * ( u(k,j,i+1) + u(k,j,i-2) )                          &
                                  +           ( u(k,j,i+2) + u(k,j,i-3) ) ) * adv_mom_5
             diss_l_u(k,j,tn)   = -ABS(u_comp_l) * (                                               &
                                    10.0_wp * ( u(k,j,i) - u(k,j,i-1)   )                          &
                                  -  5.0_wp * ( u(k,j,i+1) - u(k,j,i-2) )                          &
                                  +           ( u(k,j,i+2) - u(k,j,i-3) ) ) * adv_mom_5

#endif
             v_comp(k) = v(k,j+1,i) + v(k,j+1,i-1) - gv
             flux_n(k) = v_comp(k) * (                                                             &
                            37.0_wp * ( u(k,j+1,i) + u(k,j,i)   )                                  &
                         -   8.0_wp * ( u(k,j+2,i) + u(k,j-1,i) )                                  &
                         +            ( u(k,j+3,i) + u(k,j-2,i) ) ) * adv_mom_5
             diss_n(k) = -ABS( v_comp(k) ) * (                                                     &
                            10.0_wp * ( u(k,j+1,i) - u(k,j,i)   )                                  &
                         -   5.0_wp * ( u(k,j+2,i) - u(k,j-1,i) )                                  &
                         +            ( u(k,j+3,i) - u(k,j-2,i) ) ) * adv_mom_5
#ifdef _OPENACC
!
!--          Recompute the south fluxes.
             v_comp_s         = v(k,j,i) + v(k,j,i-1) - gv
             flux_s_u(k,tn)   = v_comp_s * (                                                       &
                                  37.0_wp * ( u(k,j,i) + u(k,j-1,i)   )                            &
                                -  8.0_wp * ( u(k,j+1,i) + u(k,j-2,i) )                            &
                                +           ( u(k,j+2,i) + u(k,j-3,i) ) ) * adv_mom_5
             diss_s_u(k,tn)   = -ABS( v_comp_s ) * (                                               &
                                   10.0_wp * ( u(k,j,i) - u(k,j-1,i)   )                           &
                                -   5.0_wp * ( u(k,j+1,i) - u(k,j-2,i) )                           &
                                +            ( u(k,j+2,i) - u(k,j-3,i) ) ) * adv_mom_5
#endif
          ENDDO
!
!--       Now, compute vertical fluxes. Split loop into a part treating the lowest 2 grid points
!--       with indirect indexing, a main loop without indirect indexing, and a loop for the
!--       uppermost 2 grid points with indirect indexing. This allows better vectorization for the
!--       main loop. First, compute the flux at model surface, which need has to be calculated
!--       explicetely for the tendency at the first w-level. For topography wall this is done
!--       implicitely by advc_flags_m.
          k = nzb
          ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )
          w_comp(k) = w(k,j,i) + w(k,j,i-1)
          flux_t(nzb) = w_comp(k)         * rho_air_zw(k) * ( ibit6 * adv_mom_1 )                  &
                                                          * ( u(k+1,j,i) + u(k,j,i) )
          diss_t(nzb) = -ABS( w_comp(k) ) * rho_air_zw(k) * ( ibit6 * adv_mom_1 )                  &
                                                          * ( u(k+1,j,i) - u(k,j,i) )
          DO  k = nzb+1, nzb+1
!
!--          k index has to be modified near bottom and top, else array subscripts will be exceeded.
             ibit8 = REAL( IBITS(advc_flags_m(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flags_m(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )

             k_ppp = k + 3 * ibit8
             k_pp  = k + 2 * ( 1 - ibit6 )
             k_mm  = k - 2 * ibit8

             w_comp(k) = w(k,j,i) + w(k,j,i-1)
             flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                             &
                        (  37.0_wp * ibit8 * adv_mom_5                                             &
                        +   7.0_wp * ibit7 * adv_mom_3                                             &
                        +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)   + u(k,j,i)    )          &
                        - ( 8.0_wp * ibit8 * adv_mom_5                                             &
                        +            ibit7 * adv_mom_3 ) * ( u(k_pp,j,i)  + u(k-1,j,i)  )          &
                        + (          ibit8 * adv_mom_5 ) * ( u(k_ppp,j,i) + u(k_mm,j,i) )          &
                                                     )

             diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                     &
                        (  10.0_wp * ibit8 * adv_mom_5                                             &
                        +   3.0_wp * ibit7 * adv_mom_3                                             &
                        +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)   - u(k,j,i)    )          &
                        - ( 5.0_wp * ibit8 * adv_mom_5                                             &
                        +            ibit7 * adv_mom_3 ) * ( u(k_pp,j,i)  - u(k-1,j,i)  )          &
                        + (          ibit8 * adv_mom_5 ) * ( u(k_ppp,j,i) - u(k_mm,j,i) )          &
                                                             )
          ENDDO

          DO  k = nzb+2, nzt-2

             ibit8 = REAL( IBITS(advc_flags_m(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flags_m(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )

             w_comp(k) = w(k,j,i) + w(k,j,i-1)
             flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                             &
                        (  37.0_wp * ibit8 * adv_mom_5                                             &
                        +   7.0_wp * ibit7 * adv_mom_3                                             &
                        +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)  + u(k,j,i)   )            &
                        - ( 8.0_wp * ibit8 * adv_mom_5                                             &
                        +            ibit7 * adv_mom_3 ) * ( u(k+2,j,i)  + u(k-1,j,i) )            &
                        + (          ibit8 * adv_mom_5 ) * ( u(k+3,j,i)  + u(k-2,j,i) )            &
                                                     )

             diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                     &
                         (  10.0_wp * ibit8 * adv_mom_5                                            &
                         +   3.0_wp * ibit7 * adv_mom_3                                            &
                         +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)  - u(k,j,i)    )          &
                         - ( 5.0_wp * ibit8 * adv_mom_5                                            &
                         +            ibit7 * adv_mom_3 ) * ( u(k+2,j,i)  - u(k-1,j,i)  )          &
                         + (          ibit8 * adv_mom_5 ) * ( u(k+3,j,i)  - u(k-2,j,i)  )          &
                                                             )
          ENDDO

          DO  k = nzt-1, nzt-symmetry_flag
!
!--          k index has to be modified near bottom and top, else array subscripts will be exceeded.
             ibit8 = REAL( IBITS(advc_flags_m(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flags_m(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )
!
!--          k index needs to be modified near bottom and top, else array subscripts will be
!--          exceeded. (k_ppp = k+3, k_mm = k) if 5th order scheme is employed, else
!--          (k_ppp = k, k_mm = k). Special treatment is required for k_pp. This must be limited
!--          at k = nzt, where either the first order scheme is employed (ibit6), or all bit flags
!--          are zero (in case there is a rigid lid defined at or near the domain top).
             k_ppp = k + 3 * ibit8
             k_pp  = k + 2 * ( 1 - ibit6 ) * ( ibit6 + ibit7 + ibit8 )
             k_mm  = k - 2 * ibit8

             w_comp(k) = w(k,j,i) + w(k,j,i-1)
             flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                             &
                         (  37.0_wp * ibit8 * adv_mom_5                                            &
                         +   7.0_wp * ibit7 * adv_mom_3                                            &
                         +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)   + u(k,j,i)    )         &
                         - ( 8.0_wp * ibit8 * adv_mom_5                                            &
                         +            ibit7 * adv_mom_3 ) * ( u(k_pp,j,i)  + u(k-1,j,i)  )         &
                         + (          ibit8 * adv_mom_5 ) * ( u(k_ppp,j,i) + u(k_mm,j,i) )         &
                                                     )

             diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                     &
                         (  10.0_wp * ibit8 * adv_mom_5                                            &
                         +   3.0_wp * ibit7 * adv_mom_3                                            &
                         +            ibit6 * adv_mom_1 ) * ( u(k+1,j,i)   - u(k,j,i)    )         &
                         - ( 5.0_wp * ibit8 * adv_mom_5                                            &
                         +            ibit7 * adv_mom_3 ) * ( u(k_pp,j,i)  - u(k-1,j,i)  )         &
                         + (          ibit8 * adv_mom_5 ) * ( u(k_ppp,j,i) - u(k_mm,j,i) )         &
                                                             )
          ENDDO
!
!--       Set resolved/turbulent flux at model top to zero (w-level). In case that a symmetric
!--       behavior between bottom and top shell be guaranteed (closed channel flow), the flux at nzt
!--       is also set to zero.
          IF ( symmetry_flag == 1 )  THEN
             flux_t(nzt) = 0.0_wp
             diss_t(nzt) = 0.0_wp
             w_comp(nzt) = 0.0_wp
          ENDIF
          flux_t(nzt+1) = 0.0_wp
          diss_t(nzt+1) = 0.0_wp
          w_comp(nzt+1) = 0.0_wp

          DO  k = nzb+1, nzb_max_l

             flux_d    = flux_t(k-1)
             diss_d    = diss_t(k-1)

             ibit2 = REAL( IBITS(advc_flags_m(k,j,i),2,1), KIND = wp )
             ibit1 = REAL( IBITS(advc_flags_m(k,j,i),1,1), KIND = wp )
             ibit0 = REAL( IBITS(advc_flags_m(k,j,i),0,1), KIND = wp )

             ibit5 = REAL( IBITS(advc_flags_m(k,j,i),5,1), KIND = wp )
             ibit4 = REAL( IBITS(advc_flags_m(k,j,i),4,1), KIND = wp )
             ibit3 = REAL( IBITS(advc_flags_m(k,j,i),3,1), KIND = wp )

             ibit8 = REAL( IBITS(advc_flags_m(k,j,i),8,1), KIND = wp )
             ibit7 = REAL( IBITS(advc_flags_m(k,j,i),7,1), KIND = wp )
             ibit6 = REAL( IBITS(advc_flags_m(k,j,i),6,1), KIND = wp )
!
!--          Calculate the divergence of the velocity field. A respective correction is needed to
!--          overcome numerical instabilities introduced by an insufficient reduction of divergences
!--          near topography.
             div =   ( ( u_comp(k)        * ( ibit0 + ibit1 + ibit2 )                              &
                     - ( u(k,j,i)   + u(k,j,i-1)   )                                               &
                                       * (                                                         &
                        REAL( IBITS(advc_flags_m(k,j,i-1),0,1), KIND = wp )                        &
                      + REAL( IBITS(advc_flags_m(k,j,i-1),1,1), KIND = wp )                        &
                      + REAL( IBITS(advc_flags_m(k,j,i-1),2,1), KIND = wp )                        &
                                         )                                                         &
                     ) * ddx                                                                       &
                   + ( ( v_comp(k) + gv ) * ( ibit3 + ibit4 + ibit5 )                              &
                     - ( v(k,j,i)   + v(k,j,i-1 )  )                                               &
                                       * (                                                         &
                        REAL( IBITS(advc_flags_m(k,j-1,i),3,1), KIND = wp )                        &
                      + REAL( IBITS(advc_flags_m(k,j-1,i),4,1), KIND = wp )                        &
                      + REAL( IBITS(advc_flags_m(k,j-1,i),5,1), KIND = wp )                        &
                                         )                                                         &
                     ) * ddy                                                                       &
                   + ( w_comp(k)   * rho_air_zw(k) * ( ibit6 + ibit7 + ibit8 )                     &
                   -   w_comp(k-1) * rho_air_zw(k-1)                                               &
                                       * (                                                         &
                        REAL( IBITS(advc_flags_m(k-1,j,i),6,1), KIND = wp )                        &
                      + REAL( IBITS(advc_flags_m(k-1,j,i),7,1), KIND = wp )                        &
                      + REAL( IBITS(advc_flags_m(k-1,j,i),8,1), KIND = wp )                        &
                                         )                                                         &
                     ) * drho_air(k) * ddzw(k)                                                     &
                   ) * 0.5_wp

             tend(k,j,i) = tend(k,j,i) - (                                                         &
                               ( flux_r(k) + diss_r(k)                                             &
                             -   flux_l_u(k,j,tn) - diss_l_u(k,j,tn) ) * ddx                       &
                             + ( flux_n(k) + diss_n(k)                                             &
                             -   flux_s_u(k,tn) - diss_s_u(k,tn)     ) * ddy                       &
                             + ( ( flux_t(k) + diss_t(k) )                                         &
                             -   ( flux_d    + diss_d    )           ) * drho_air(k) * ddzw(k)     &
                                         ) + div * u(k,j,i)
#ifndef _OPENACC
!
!--          Swap fluxes. Note, in the OPENACC case these are computed again.
             flux_l_u(k,j,tn) = flux_r(k)
             diss_l_u(k,j,tn) = diss_r(k)
             flux_s_u(k,tn)   = flux_n(k)
             diss_s_u(k,tn)   = diss_n(k)
#endif
!
!--          Statistical Evaluation of u'u'. The factor has to be applied for right evaluation when
!--          gallilei_trans = .T. .
             sums_us2_ws_l(k,tn) = sums_us2_ws_l(k,tn) +                                           &
                                   ( flux_r(k)                                                     &
                                     * ( u_comp(k) - 2.0_wp * hom(k,1,1,0)                   )     &
                                     / ( u_comp(k) - gu + SIGN( 1.0E-20_wp, u_comp(k) - gu ) )     &
                                   + diss_r(k)                                                     &
                                     *   ABS( u_comp(k) - 2.0_wp * hom(k,1,1,0)              )     &
                                     / ( ABS( u_comp(k) - gu ) + 1.0E-20_wp                  )     &
                                   ) *   weight_substep(intermediate_timestep_count)
!
!--          Statistical Evaluation of w'u'.
             sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn) +                                         &
                                    ( flux_t(k)                                                    &
                                      * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                   )    &
                                      / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )           )    &
                                    + diss_t(k)                                                    &
                                      *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)              )    &
                                      / ( ABS( w_comp(k) ) + 1.0E-20_wp                       )    &
                                    ) *   weight_substep(intermediate_timestep_count)
          ENDDO

          DO  k = nzb_max_l+1, nzt

             flux_d    = flux_t(k-1)
             diss_d    = diss_t(k-1)
!
!--          Calculate the divergence of the velocity field. A respective correction is needed to
!--          overcome numerical instabilities introduced by an insufficient reduction of divergences
!--          near topography.
             div = ( ( u_comp(k)       - ( u(k,j,i)   + u(k,j,i-1) ) ) * ddx                       &
                   + ( v_comp(k) + gv  - ( v(k,j,i)   + v(k,j,i-1) ) ) * ddy                       &
                   + ( w_comp(k)   * rho_air_zw(k)                                                 &
                   -   w_comp(k-1) * rho_air_zw(k-1)                                               &
                     ) * drho_air(k) * ddzw(k)                                                     &
                   ) * 0.5_wp

             tend(k,j,i) = tend(k,j,i) - (                                                         &
                               ( flux_r(k) + diss_r(k)                                             &
                             -   flux_l_u(k,j,tn) - diss_l_u(k,j,tn) )     * ddx                   &
                             + ( flux_n(k) + diss_n(k)                                             &
                             -   flux_s_u(k,tn) - diss_s_u(k,tn)     )     * ddy                   &
                             + ( ( flux_t(k) + diss_t(k) )                                         &
                             -   ( flux_d    + diss_d    ) ) * drho_air(k) * ddzw(k)               &
                                         ) + div * u(k,j,i)
#ifndef _OPENACC
             flux_l_u(k,j,tn) = flux_r(k)
             diss_l_u(k,j,tn) = diss_r(k)
             flux_s_u(k,tn)   = flux_n(k)
             diss_s_u(k,tn)   = diss_n(k)
#endif
!
!--          Statistical Evaluation of u'u'. The factor has to be applied for right evaluation when
!--          gallilei_trans = .T. .
             sums_us2_ws_l(k,tn) = sums_us2_ws_l(k,tn) +                                           &
                                   ( flux_r(k)                                                     &
                                     * ( u_comp(k) - 2.0_wp * hom(k,1,1,0)                   )     &
                                     / ( u_comp(k) - gu + SIGN( 1.0E-20_wp, u_comp(k) - gu ) )     &
                                   + diss_r(k)                                                     &
                                     *   ABS( u_comp(k) - 2.0_wp * hom(k,1,1,0)              )     &
                                     / ( ABS( u_comp(k) - gu ) + 1.0E-20_wp                  )     &
                                   ) *   weight_substep(intermediate_timestep_count)
!
!--          Statistical Evaluation of w'u'.
             sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn) +                                         &
                                    ( flux_t(k)                                                    &
                                      * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                   )    &
                                      / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )           )    &
                                    + diss_t(k)                                                    &
                                      *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)              )    &
                                      / ( ABS( w_comp(k) ) + 1.0E-20_wp                       )    &
                                    ) *   weight_substep(intermediate_timestep_count)
          ENDDO
       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(68), 'advec_u_ws', 'stop' )

 END SUBROUTINE advec_u_ws


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of v - Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_v_ws


    INTEGER(iwp) ::  i         !< grid index along x-direction
    INTEGER(iwp) ::  j         !< grid index along y-direction
    INTEGER(iwp) ::  k         !< grid index along z-direction
    INTEGER(iwp) ::  k_mm      !< k-2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_pp      !< k+2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_ppp     !< k+3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  nzb_max_l !< index indicating upper bound for order degradation of horizontal advection terms
    INTEGER(iwp) ::  tn = 0    !< number of OpenMP thread

    REAL(wp)    ::  diss_d   !< artificial dissipation term at grid box bottom
    REAL(wp)    ::  div      !< diverence on v-grid
    REAL(wp)    ::  flux_d   !< artificial 6th-order flux at grid box bottom
    REAL(wp)    ::  gu       !< Galilei-transformation velocity along x
    REAL(wp)    ::  gv       !< Galilei-transformation velocity along y
    REAL(wp)    ::  ibit9 !< flag indicating 1st-order scheme along x-direction
    REAL(wp)    ::  ibit10 !< flag indicating 3rd-order scheme along x-direction
    REAL(wp)    ::  ibit11 !< flag indicating 5th-order scheme along x-direction
#ifdef _OPENACC
    REAL(wp)    ::  ibit9_l !< flag indicating 1st-order scheme along x-direction
    REAL(wp)    ::  ibit10_l !< flag indicating 3rd-order scheme along x-direction
    REAL(wp)    ::  ibit11_l !< flag indicating 5th-order scheme along x-direction
#endif
    REAL(wp)    ::  ibit12 !< flag indicating 1st-order scheme along y-direction
    REAL(wp)    ::  ibit13 !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)    ::  ibit14 !< flag indicating 5th-order scheme along y-direction
#ifdef _OPENACC
    REAL(wp)    ::  ibit12_s !< flag indicating 1st-order scheme along y-direction
    REAL(wp)    ::  ibit13_s !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)    ::  ibit14_s !< flag indicating 5th-order scheme along y-direction
#endif
    REAL(wp)    ::  ibit15   !< flag indicating 1st-order scheme along z-direction
    REAL(wp)    ::  ibit16   !< flag indicating 3rd-order scheme along z-direction
    REAL(wp)    ::  ibit17   !< flag indicating 5th-order scheme along z-direction
#ifdef _OPENACC
    REAL(wp)    ::  u_comp_l !< advection velocity along x at leftward side
    REAL(wp)    ::  v_comp_s !< advection velocity along y at southward side
#endif

    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_n !< discretized artificial dissipation at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_r !< discretized artificial dissipation at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_t !< discretized artificial dissipation at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_n !< discretized 6th-order flux at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_r !< discretized 6th-order flux at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_t !< discretized 6th-order flux at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  u_comp !< advection velocity along x
    REAL(wp), DIMENSION(nzb:nzt+1) ::  v_comp !< advection velocity along y
    REAL(wp), DIMENSION(nzb:nzt+1) ::  w_comp !< advection velocity along z

    CALL cpu_log( log_point_s(69), 'advec_v_ws', 'start' )

    gu = 2.0_wp * u_gtrans
    gv = 2.0_wp * v_gtrans

    !$ACC PARALLEL LOOP COLLAPSE(2) FIRSTPRIVATE(tn, gu, gv) &
    !$ACC PRIVATE(i, j, k, k_mm, k_pp, k_ppp) &
    !$ACC PRIVATE(ibit9, ibit10, ibit11, ibit12, ibit13, ibit14) &
    !$ACC PRIVATE(ibit15, ibit16, ibit17) &
    !$ACC PRIVATE(ibit9_l, ibit10_l, ibit11_l) &
    !$ACC PRIVATE(ibit12_s, ibit13_s, ibit14_s) &
    !$ACC PRIVATE(nzb_max_l) &
    !$ACC PRIVATE(flux_r, diss_r) &
    !$ACC PRIVATE(flux_n, diss_n) &
    !$ACC PRIVATE(flux_t, diss_t, flux_d, diss_d) &
    !$ACC PRIVATE(flux_l_v, diss_l_v, flux_s_v, diss_s_v) &
    !$ACC PRIVATE(div, u_comp, u_comp_l, v_comp, v_comp_s, w_comp) &
    !$ACC PRESENT(advc_flags_m) &
    !$ACC PRESENT(u, v, w) &
    !$ACC PRESENT(drho_air, rho_air_zw, ddzw) &
    !$ACC PRESENT(tend) &
    !$ACC PRESENT(hom(:,1,1:3,0)) &
    !$ACC PRESENT(weight_substep(intermediate_timestep_count)) &
    !$ACC PRESENT(sums_vs2_ws_l, sums_wsvs_ws_l)
    DO  i = nxl, nxr
       DO  j = nysv, nyn
!
!--       Used local modified copy of nzb_max (used to degrade order of discretization) at
!--       non-cyclic boundaries. Modify only at relevant points instead of the entire subdomain.
!--       This should lead to better load balance between boundary and non-boundary PEs.
          IF( ( bc_dirichlet_l  .OR.  bc_radiation_l )  .AND.  i <= nxl + 2  .OR. &
              ( bc_dirichlet_r  .OR.  bc_radiation_r )  .AND.  i >= nxr - 2  .OR. &
              ( bc_dirichlet_s  .OR.  bc_radiation_s )  .AND.  j <= nys + 2  .OR. &
              ( bc_dirichlet_n  .OR.  bc_radiation_n )  .AND.  j >= nyn - 2 )  THEN
             nzb_max_l = nzt
          ELSE
             nzb_max_l = nzb_max
          END IF

#ifndef _OPENACC
          IF ( i == nxl )  THEN
             DO  k = nzb+1, nzb_max_l

                ibit11 = REAL( IBITS(advc_flags_m(k,j,i-1),11,1), KIND = wp )
                ibit10 = REAL( IBITS(advc_flags_m(k,j,i-1),10,1), KIND = wp )
                ibit9  = REAL( IBITS(advc_flags_m(k,j,i-1),9,1),  KIND = wp )

                u_comp(k)        = u(k,j-1,i) + u(k,j,i) - gu
                flux_l_v(k,j,tn) = u_comp(k) * (                                                   &
                                   (  37.0_wp * ibit11 * adv_mom_5                                 &
                                   +   7.0_wp * ibit10 * adv_mom_3                                 &
                                   +             ibit9 * adv_mom_1 ) * ( v(k,j,i)   + v(k,j,i-1) ) &
                                   - ( 8.0_wp * ibit11 * adv_mom_5                                 &
                                   +            ibit10 * adv_mom_3 ) * ( v(k,j,i+1) + v(k,j,i-2) ) &
                                   + (          ibit11 * adv_mom_5 ) * ( v(k,j,i+2) + v(k,j,i-3) ) &
                                               )

                diss_l_v(k,j,tn) = -ABS( u_comp(k) ) * (                                           &
                                   (  10.0_wp * ibit11 * adv_mom_5                                 &
                                   +   3.0_wp * ibit10 * adv_mom_3                                 &
                                   +             ibit9 * adv_mom_1 ) * ( v(k,j,i)   - v(k,j,i-1) ) &
                                   - ( 5.0_wp * ibit11 * adv_mom_5                                 &
                                   +            ibit10 * adv_mom_3 ) * ( v(k,j,i+1) - v(k,j,i-2) ) &
                                   + (          ibit11 * adv_mom_5 ) * ( v(k,j,i+2) - v(k,j,i-3) ) &
                                                       )

             ENDDO

             DO  k = nzb_max_l+1, nzt

                u_comp(k)        = u(k,j-1,i) + u(k,j,i) - gu
                flux_l_v(k,j,tn) = u_comp(k) * (                                                   &
                                     37.0_wp * ( v(k,j,i) + v(k,j,i-1)   )                         &
                                   -  8.0_wp * ( v(k,j,i+1) + v(k,j,i-2) )                         &
                                   +           ( v(k,j,i+2) + v(k,j,i-3) ) ) * adv_mom_5
                diss_l_v(k,j,tn) = -ABS( u_comp(k) ) * (                                           &
                                     10.0_wp * ( v(k,j,i) - v(k,j,i-1)   )                         &
                                   -  5.0_wp * ( v(k,j,i+1) - v(k,j,i-2) )                         &
                                   +           ( v(k,j,i+2) - v(k,j,i-3) ) ) * adv_mom_5

             ENDDO
          ENDIF

          IF ( j == nysv )  THEN
             DO  k = nzb+1, nzb_max_l

                ibit14 = REAL( IBITS(advc_flags_m(k,j-1,i),14,1), KIND = wp )
                ibit13 = REAL( IBITS(advc_flags_m(k,j-1,i),13,1), KIND = wp )
                ibit12 = REAL( IBITS(advc_flags_m(k,j-1,i),12,1), KIND = wp )

                v_comp(k)      = v(k,j,i) + v(k,j-1,i) - gv
                flux_s_v(k,tn) = v_comp(k) * (                                                     &
                                 (  37.0_wp * ibit14 * adv_mom_5                                   &
                                 +   7.0_wp * ibit13 * adv_mom_3                                   &
                                 +            ibit12 * adv_mom_1 ) * ( v(k,j,i)   + v(k,j-1,i) )   &
                                 - ( 8.0_wp * ibit14 * adv_mom_5                                   &
                                 +            ibit13 * adv_mom_3 ) * ( v(k,j+1,i) + v(k,j-2,i) )   &
                                 + (          ibit14 * adv_mom_5 ) * ( v(k,j+2,i) + v(k,j-3,i) )   &
                                             )

                diss_s_v(k,tn) = -ABS( v_comp(k) ) * (                                             &
                                   (  10.0_wp * ibit14 * adv_mom_5                                 &
                                   +   3.0_wp * ibit13 * adv_mom_3                                 &
                                   +            ibit12 * adv_mom_1 ) * ( v(k,j,i)   - v(k,j-1,i) ) &
                                   - ( 5.0_wp * ibit14 * adv_mom_5                                 &
                                   +            ibit13 * adv_mom_3 ) * ( v(k,j+1,i) - v(k,j-2,i) ) &
                                   + (          ibit14 * adv_mom_5 ) * ( v(k,j+2,i) - v(k,j-3,i) ) &
                                                     )

             ENDDO

             DO  k = nzb_max_l+1, nzt

                v_comp(k)      = v(k,j,i) + v(k,j-1,i) - gv
                flux_s_v(k,tn) = v_comp(k) * (                                                     &
                                   37.0_wp * ( v(k,j,i) + v(k,j-1,i)   )                           &
                                 -  8.0_wp * ( v(k,j+1,i) + v(k,j-2,i) )                           &
                                 +           ( v(k,j+2,i) + v(k,j-3,i) ) ) * adv_mom_5
                diss_s_v(k,tn) = -ABS( v_comp(k) ) * (                                             &
                                   10.0_wp * ( v(k,j,i) - v(k,j-1,i)   )                           &
                                 -  5.0_wp * ( v(k,j+1,i) - v(k,j-2,i) )                           &
                                 +           ( v(k,j+2,i) - v(k,j-3,i) ) ) * adv_mom_5

             ENDDO
          ENDIF
#endif

          DO  k = nzb+1, nzb_max_l

             ibit11 = REAL( IBITS(advc_flags_m(k,j,i),11,1), KIND = wp )
             ibit10 = REAL( IBITS(advc_flags_m(k,j,i),10,1), KIND = wp )
             ibit9  = REAL( IBITS(advc_flags_m(k,j,i),9,1),  KIND = wp )

             u_comp(k) = u(k,j-1,i+1) + u(k,j,i+1) - gu
             flux_r(k) = u_comp(k) * (                                                             &
                         (  37.0_wp * ibit11 * adv_mom_5                                           &
                         +   7.0_wp * ibit10 * adv_mom_3                                           &
                         +             ibit9 * adv_mom_1 ) * ( v(k,j,i+1) + v(k,j,i)   )           &
                         - ( 8.0_wp * ibit11 * adv_mom_5                                           &
                         +            ibit10 * adv_mom_3 ) * ( v(k,j,i+2) + v(k,j,i-1) )           &
                         + (          ibit11 * adv_mom_5 ) * ( v(k,j,i+3) + v(k,j,i-2) )           &
                                     )

             diss_r(k) = -ABS( u_comp(k) ) * (                                                     &
                         (  10.0_wp * ibit11 * adv_mom_5                                           &
                         +   3.0_wp * ibit10 * adv_mom_3                                           &
                         +             ibit9 * adv_mom_1 ) * ( v(k,j,i+1) - v(k,j,i)  )            &
                         - ( 5.0_wp * ibit11 * adv_mom_5                                           &
                         +            ibit10 * adv_mom_3 ) * ( v(k,j,i+2) - v(k,j,i-1) )           &
                         + (          ibit11 * adv_mom_5 ) * ( v(k,j,i+3) - v(k,j,i-2) )           &
                                             )

#ifdef _OPENACC
!
!--          Recompute the left fluxes.
             ibit11_l = REAL( IBITS(advc_flags_m(k,j,i-1),11,1), KIND = wp )
             ibit10_l = REAL( IBITS(advc_flags_m(k,j,i-1),10,1), KIND = wp )
             ibit9_l  = REAL( IBITS(advc_flags_m(k,j,i-1),9,1),  KIND = wp )

             u_comp_l          = u(k,j-1,i) + u(k,j,i) - gu
             flux_l_v(k,j,tn)  = u_comp_l * (                                                      &
                                 (  37.0_wp * ibit11_l * adv_mom_5                                 &
                                 +   7.0_wp * ibit10_l * adv_mom_3                                 &
                                 +             ibit9_l * adv_mom_1 ) * ( v(k,j,i)   + v(k,j,i-1) ) &
                                 - ( 8.0_wp * ibit11_l * adv_mom_5                                 &
                                 +            ibit10_l * adv_mom_3 ) * ( v(k,j,i+1) + v(k,j,i-2) ) &
                                 + (          ibit11_l * adv_mom_5 ) * ( v(k,j,i+2) + v(k,j,i-3) ) &
                                            )

              diss_l_v(k,j,tn) = -ABS( u_comp_l ) * (                                              &
                                 (  10.0_wp * ibit11_l * adv_mom_5                                 &
                                 +   3.0_wp * ibit10_l * adv_mom_3                                 &
                                 +             ibit9_l * adv_mom_1 ) * ( v(k,j,i)   - v(k,j,i-1) ) &
                                 - ( 5.0_wp * ibit11_l * adv_mom_5                                 &
                                 +            ibit10_l * adv_mom_3 ) * ( v(k,j,i+1) - v(k,j,i-2) ) &
                                 + (          ibit11_l * adv_mom_5 ) * ( v(k,j,i+2) - v(k,j,i-3) ) &
                                                    )
#endif

             ibit14 = REAL( IBITS(advc_flags_m(k,j,i),14,1), KIND = wp )
             ibit13 = REAL( IBITS(advc_flags_m(k,j,i),13,1), KIND = wp )
             ibit12 = REAL( IBITS(advc_flags_m(k,j,i),12,1), KIND = wp )

             v_comp(k) = v(k,j+1,i) + v(k,j,i)
             flux_n(k) = ( v_comp(k) - gv ) * (                                                    &
                         (  37.0_wp * ibit14 * adv_mom_5                                           &
                         +   7.0_wp * ibit13 * adv_mom_3                                           &
                         +            ibit12 * adv_mom_1 ) * ( v(k,j+1,i) + v(k,j,i)   )           &
                         - ( 8.0_wp * ibit14 * adv_mom_5                                           &
                         +            ibit13 * adv_mom_3 ) * ( v(k,j+2,i) + v(k,j-1,i) )           &
                         + (          ibit14 * adv_mom_5 ) * ( v(k,j+3,i) + v(k,j-2,i) )           &
                                              )

             diss_n(k) = -ABS( v_comp(k) - gv ) * (                                                &
                         (  10.0_wp * ibit14 * adv_mom_5                                           &
                         +   3.0_wp * ibit13 * adv_mom_3                                           &
                         +            ibit12 * adv_mom_1 ) * ( v(k,j+1,i) - v(k,j,i)  )            &
                         - ( 5.0_wp * ibit14 * adv_mom_5                                           &
                         +            ibit13 * adv_mom_3 ) * ( v(k,j+2,i) - v(k,j-1,i) )           &
                         + (          ibit14 * adv_mom_5 ) * ( v(k,j+3,i) - v(k,j-2,i) )           &
                                                  )

#ifdef _OPENACC
!
!--          Recompute the south fluxes.
             ibit14_s = REAL( IBITS(advc_flags_m(k,j-1,i),14,1), KIND = wp )
             ibit13_s = REAL( IBITS(advc_flags_m(k,j-1,i),13,1), KIND = wp )
             ibit12_s = REAL( IBITS(advc_flags_m(k,j-1,i),12,1), KIND = wp )

             v_comp_s         = v(k,j,i) + v(k,j-1,i) - gv
             flux_s_v(k,tn)   = v_comp_s * (                                                       &
                                (  37.0_wp * ibit14_s * adv_mom_5                                  &
                                +   7.0_wp * ibit13_s * adv_mom_3                                  &
                                +            ibit12_s * adv_mom_1 ) * ( v(k,j,i)   + v(k,j-1,i) )  &
                                - ( 8.0_wp * ibit14_s * adv_mom_5                                  &
                                +            ibit13_s * adv_mom_3 ) * ( v(k,j+1,i) + v(k,j-2,i) )  &
                                + (          ibit14_s * adv_mom_5 ) * ( v(k,j+2,i) + v(k,j-3,i) )  &
                                           )

            diss_s_v(k,tn)   = -ABS( v_comp_s ) * (                                                &
                               (  10.0_wp * ibit14_s * adv_mom_5                                   &
                               +   3.0_wp * ibit13_s * adv_mom_3                                   &
                               +            ibit12_s * adv_mom_1 ) * ( v(k,j,i)   - v(k,j-1,i) )   &
                               - ( 5.0_wp * ibit14_s * adv_mom_5                                   &
                               +            ibit13_s * adv_mom_3 ) * ( v(k,j+1,i) - v(k,j-2,i) )   &
                               + (          ibit14_s * adv_mom_5 ) * ( v(k,j+2,i) - v(k,j-3,i) )   &
                                                  )
#endif
          ENDDO

          DO  k = nzb_max_l+1, nzt

             u_comp(k) = u(k,j-1,i+1) + u(k,j,i+1) - gu
             flux_r(k) = u_comp(k) * (                                                             &
                           37.0_wp * ( v(k,j,i+1) + v(k,j,i)   )                                   &
                         -  8.0_wp * ( v(k,j,i+2) + v(k,j,i-1) )                                   &
                         +           ( v(k,j,i+3) + v(k,j,i-2) ) ) * adv_mom_5

             diss_r(k) = -ABS( u_comp(k) ) * (                                                     &
                           10.0_wp * ( v(k,j,i+1) - v(k,j,i) )                                     &
                         -  5.0_wp * ( v(k,j,i+2) - v(k,j,i-1) )                                   &
                         +           ( v(k,j,i+3) - v(k,j,i-2) ) ) * adv_mom_5

#ifdef _OPENACC
!
!--          Recompute the left fluxes.
             u_comp_l           = u(k,j-1,i) + u(k,j,i) - gu
             flux_l_v(k,j,tn)   = u_comp_l * (                                                     &
                                    37.0_wp * ( v(k,j,i) + v(k,j,i-1)   )                          &
                                  -  8.0_wp * ( v(k,j,i+1) + v(k,j,i-2) )                          &
                                  +           ( v(k,j,i+2) + v(k,j,i-3) ) ) * adv_mom_5
             diss_l_v(k,j,tn)   = -ABS( u_comp_l ) * (                                             &
                                    10.0_wp * ( v(k,j,i) - v(k,j,i-1)   )                          &
                                  -  5.0_wp * ( v(k,j,i+1) - v(k,j,i-2) )                          &
                                  +           ( v(k,j,i+2) - v(k,j,i-3) ) ) * adv_mom_5
#endif

             v_comp(k) = v(k,j+1,i) + v(k,j,i)
             flux_n(k) = ( v_comp(k) - gv ) * (                                                    &
                           37.0_wp * ( v(k,j+1,i) + v(k,j,i)   )                                   &
                         -  8.0_wp * ( v(k,j+2,i) + v(k,j-1,i) )                                   &
                         +           ( v(k,j+3,i) + v(k,j-2,i) ) ) * adv_mom_5

             diss_n(k) = -ABS( v_comp(k) - gv ) * (                                                &
                           10.0_wp * ( v(k,j+1,i) - v(k,j,i)   )                                   &
                         -  5.0_wp * ( v(k,j+2,i) - v(k,j-1,i) )                                   &
                         +           ( v(k,j+3,i) - v(k,j-2,i) ) ) * adv_mom_5

#ifdef _OPENACC
!
!--          Recompute the south fluxes.
             v_comp_s = v(k,j,i) + v(k,j-1,i) - gv
             flux_s_v(k,tn)   = v_comp_s  * (                                                      &
                                  37.0_wp * ( v(k,j,i) + v(k,j-1,i)   )                            &
                                -  8.0_wp * ( v(k,j+1,i) + v(k,j-2,i) )                            &
                                +           ( v(k,j+2,i) + v(k,j-3,i) ) ) * adv_mom_5
             diss_s_v(k,tn)   = -ABS( v_comp_s ) * (                                               &
                                  10.0_wp * ( v(k,j,i) - v(k,j-1,i)   )                            &
                                -  5.0_wp * ( v(k,j+1,i) - v(k,j-2,i) )                            &
                                +           ( v(k,j+2,i) - v(k,j-3,i) ) ) * adv_mom_5
#endif
          ENDDO
!
!--       Now, compute vertical fluxes. Split loop into a part treating the lowest 2 grid points
!--       with indirect indexing, a main loop without indirect indexing, and a loop for the
!--       uppermost 2 grid points with indirect indexing. This allows better vectorization for the
!--       main loop.
!--       First, compute the flux at model surface, which has to be calculated explicetely for the
!--       tendency at the first w-level. For topography wall this is done implicitely by
!--       advc_flags_m.
          k = nzb
          ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )
          w_comp(k) = w(k,j-1,i) + w(k,j,i)
          flux_t(nzb) = w_comp(k)         * rho_air_zw(k) * ( ibit15 * adv_mom_1 )                 &
                                                          * ( v(k+1,j,i) + v(k,j,i) )
          diss_t(nzb) = -ABS( w_comp(k) ) * rho_air_zw(k) * ( ibit15 * adv_mom_1 )                 &
                                                          * ( v(k+1,j,i) - v(k,j,i) )
          DO  k = nzb+1, nzb+1
!
!--          k index has to be modified near bottom and top, else array subscripts will be exceeded.
             ibit17 = REAL( IBITS(advc_flags_m(k,j,i),17,1), KIND = wp )
             ibit16 = REAL( IBITS(advc_flags_m(k,j,i),16,1), KIND = wp )
             ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )

             k_ppp = k + 3 * ibit17
             k_pp  = k + 2 * ( 1 - ibit15  )
             k_mm  = k - 2 * ibit17

             w_comp(k) = w(k,j-1,i) + w(k,j,i)
             flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                             &
                         (  37.0_wp * ibit17 * adv_mom_5                                           &
                         +   7.0_wp * ibit16 * adv_mom_3                                           &
                         +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i)   + v(k,j,i)    )        &
                         - ( 8.0_wp * ibit17 * adv_mom_5                                           &
                         +            ibit16 * adv_mom_3 ) * ( v(k_pp,j,i)  + v(k-1,j,i)  )        &
                         + (          ibit17 * adv_mom_5 ) * ( v(k_ppp,j,i) + v(k_mm,j,i) )        &
                                                     )

             diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                     &
                         (  10.0_wp * ibit17 * adv_mom_5                                           &
                         +   3.0_wp * ibit16 * adv_mom_3                                           &
                         +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i)   - v(k,j,i)    )        &
                         - ( 5.0_wp * ibit17 * adv_mom_5                                           &
                         +            ibit16 * adv_mom_3 ) * ( v(k_pp,j,i)  - v(k-1,j,i)  )        &
                         + (          ibit17 * adv_mom_5 ) * ( v(k_ppp,j,i) - v(k_mm,j,i) )        &
                                                             )
          ENDDO

          DO  k = nzb+2, nzt-2
             ibit17 = REAL( IBITS(advc_flags_m(k,j,i),17,1), KIND = wp )
             ibit16 = REAL( IBITS(advc_flags_m(k,j,i),16,1), KIND = wp )
             ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )

             w_comp(k) = w(k,j-1,i) + w(k,j,i)
             flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                             &
                       (  37.0_wp * ibit17 * adv_mom_5                                             &
                       +   7.0_wp * ibit16 * adv_mom_3                                             &
                       +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i)  + v(k,j,i)    )           &
                       - ( 8.0_wp * ibit17 * adv_mom_5                                             &
                       +            ibit16 * adv_mom_3 ) * ( v(k+2,j,i) + v(k-1,j,i)   )           &
                       + (          ibit17 * adv_mom_5 ) * ( v(k+3,j,i) + v(k-2,j,i)   )           &
                                                     )

             diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                     &
                         (  10.0_wp * ibit17 * adv_mom_5                                           &
                         +   3.0_wp * ibit16 * adv_mom_3                                           &
                         +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i) - v(k,j,i)     )         &
                         - ( 5.0_wp * ibit17 * adv_mom_5                                           &
                         +            ibit16 * adv_mom_3 ) * ( v(k+2,j,i) - v(k-1,j,i)   )         &
                         + (          ibit17 * adv_mom_5 ) * ( v(k+3,j,i) - v(k-2,j,i)   )         &
                                                             )
          ENDDO

          DO  k = nzt-1, nzt-symmetry_flag
!
!--          k index has to be modified near bottom and top, else array subscripts will be exceeded.
             ibit17 = REAL( IBITS(advc_flags_m(k,j,i),17,1), KIND = wp )
             ibit16 = REAL( IBITS(advc_flags_m(k,j,i),16,1), KIND = wp )
             ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )
!
!--          k index needs to be modified near bottom and top, else array subscripts will be
!--          exceeded. (k_ppp = k+3, k_mm = k) if 5th order scheme is employed, else
!--          (k_ppp = k, k_mm = k). Special treatment is required for k_pp. This must be limited
!--          at k = nzt, where either the first order scheme is employed (ibit15), or all bit flags
!--          are zero (in case there is a rigid lid defined at or near the domain top).
             k_ppp = k + 3 * ibit17
             k_pp  = k + 2 * ( 1 - ibit15  ) * ( ibit15 + ibit16 + ibit17 )
             k_mm  = k - 2 * ibit17

             w_comp(k) = w(k,j-1,i) + w(k,j,i)
             flux_t(k) = w_comp(k) * rho_air_zw(k) * (                                             &
                         (  37.0_wp * ibit17 * adv_mom_5                                           &
                         +   7.0_wp * ibit16 * adv_mom_3                                           &
                         +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i)   + v(k,j,i)    )        &
                         - ( 8.0_wp * ibit17 * adv_mom_5                                           &
                         +            ibit16 * adv_mom_3 ) * ( v(k_pp,j,i)  + v(k-1,j,i)  )        &
                         + (          ibit17 * adv_mom_5 ) * ( v(k_ppp,j,i) + v(k_mm,j,i) )        &
                                                     )

             diss_t(k) = -ABS( w_comp(k) ) * rho_air_zw(k) * (                                     &
                         (  10.0_wp * ibit17 * adv_mom_5                                           &
                         +   3.0_wp * ibit16 * adv_mom_3                                           &
                         +            ibit15 * adv_mom_1 ) * ( v(k+1,j,i)   - v(k,j,i)    )        &
                         - ( 5.0_wp * ibit17 * adv_mom_5                                           &
                         +            ibit16 * adv_mom_3 ) * ( v(k_pp,j,i)  - v(k-1,j,i)  )        &
                         + (          ibit17 * adv_mom_5 ) * ( v(k_ppp,j,i) - v(k_mm,j,i) )        &
                                                             )
          ENDDO

!
!--       Set resolved/turbulent flux at model top to zero (w-level). In case that a symmetric
!--       behavior between bottom and top shall be guaranteed (closed channel flow), the flux at nzt
!--       is also set to zero.
          IF ( symmetry_flag == 1 ) THEN
             flux_t(nzt) = 0.0_wp
             diss_t(nzt) = 0.0_wp
             w_comp(nzt) = 0.0_wp
          ENDIF
          flux_t(nzt+1) = 0.0_wp
          diss_t(nzt+1) = 0.0_wp
          w_comp(nzt+1) = 0.0_wp

          DO k = nzb+1, nzb_max_l

             flux_d    = flux_t(k-1)
             diss_d    = diss_t(k-1)

             ibit11 = REAL( IBITS(advc_flags_m(k,j,i),11,1), KIND = wp )
             ibit10 = REAL( IBITS(advc_flags_m(k,j,i),10,1), KIND = wp )
             ibit9  = REAL( IBITS(advc_flags_m(k,j,i),9,1),  KIND = wp )

             ibit14 = REAL( IBITS(advc_flags_m(k,j,i),14,1), KIND = wp )
             ibit13 = REAL( IBITS(advc_flags_m(k,j,i),13,1), KIND = wp )
             ibit12 = REAL( IBITS(advc_flags_m(k,j,i),12,1), KIND = wp )

             ibit17 = REAL( IBITS(advc_flags_m(k,j,i),17,1), KIND = wp )
             ibit16 = REAL( IBITS(advc_flags_m(k,j,i),16,1), KIND = wp )
             ibit15 = REAL( IBITS(advc_flags_m(k,j,i),15,1), KIND = wp )
!
!--          Calculate the divergence of the velocity field. A respective correction is needed to
!--          overcome numerical instabilities caused by an insufficient reduction of divergences
!--          near topography.
             div = ( ( ( u_comp(k)     + gu )                                                      &
                                          * ( ibit9 + ibit10 + ibit11 )                            &
                   - ( u(k,j-1,i)   + u(k,j,i) )                                                   &
                                          * (                                                      &
                           REAL( IBITS(advc_flags_m(k,j,i-1),9,1), KIND = wp )                     &
                         + REAL( IBITS(advc_flags_m(k,j,i-1),10,1), KIND = wp )                    &
                         + REAL( IBITS(advc_flags_m(k,j,i-1),11,1), KIND = wp )                    &
                                            )                                                      &
                     ) * ddx                                                                       &
                   + ( v_comp(k)                                                                   &
                                          * ( ibit12 + ibit13 + ibit14 )                           &
                   - ( v(k,j,i)     + v(k,j-1,i) )                                                 &
                                          * (                                                      &
                           REAL( IBITS(advc_flags_m(k,j-1,i),12,1), KIND = wp )                    &
                         + REAL( IBITS(advc_flags_m(k,j-1,i),13,1), KIND = wp )                    &
                         + REAL( IBITS(advc_flags_m(k,j-1,i),14,1), KIND = wp )                    &
                                            )                                                      &
                     ) * ddy                                                                       &
                   + ( w_comp(k) * rho_air_zw(k)                                                   &
                                          * ( ibit15 + ibit16 + ibit17 )                           &
                   - ( w(k-1,j-1,i) + w(k-1,j,i) ) * rho_air_zw(k-1)                               &
                                          * (                                                      &
                           REAL( IBITS(advc_flags_m(k-1,j,i),15,1), KIND = wp )                    &
                         + REAL( IBITS(advc_flags_m(k-1,j,i),16,1), KIND = wp )                    &
                         + REAL( IBITS(advc_flags_m(k-1,j,i),17,1), KIND = wp )                    &
                                            )                                                      &
                      ) * drho_air(k) * ddzw(k)                                                    &
                   ) * 0.5_wp


             tend(k,j,i) = tend(k,j,i) - (                                                         &
                             ( ( flux_r(k) + diss_r(k) )                                           &
                           -   ( flux_l_v(k,j,tn) + diss_l_v(k,j,tn) ) ) * ddx                     &
                           + ( ( flux_n(k) + diss_n(k) )                                           &
                           -   ( flux_s_v(k,tn) + diss_s_v(k,tn) ) ) * ddy                         &
                           + ( ( flux_t(k) + diss_t(k) )                                           &
                           -   ( flux_d + diss_d ) ) * drho_air(k) * ddzw(k)                       &
                                                  )  + v(k,j,i) * div

#ifndef _OPENACC
!
!--          Swap fluxes. Note, in the OPENACC case these are computed again.
             flux_l_v(k,j,tn) = flux_r(k)
             diss_l_v(k,j,tn) = diss_r(k)
             flux_s_v(k,tn)   = flux_n(k)
             diss_s_v(k,tn)   = diss_n(k)
#endif

!
!--          Statistical Evaluation of v'v'. The factor has to be applied for right evaluation when
!--          gallilei_trans = .T. .
             !$ACC ATOMIC
             sums_vs2_ws_l(k,tn) = sums_vs2_ws_l(k,tn) +                                           &
                                   ( flux_n(k)                                                     &
                                     * ( v_comp(k) - 2.0_wp * hom(k,1,2,0)                   )     &
                                     / ( v_comp(k) - gv + SIGN( 1.0E-20_wp, v_comp(k) - gv ) )     &
                                   + diss_n(k)                                                     &
                                     *   ABS( v_comp(k) - 2.0_wp * hom(k,1,2,0)              )     &
                                     / ( ABS( v_comp(k) - gv ) + 1.0E-20_wp                  )     &
                                   ) *   weight_substep(intermediate_timestep_count)
!
!--          Statistical Evaluation of w'u'.
             !$ACC ATOMIC
             sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn) +                                         &
                                    ( flux_t(k)                                                    &
                                      * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                   )    &
                                      / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )           )    &
                                    + diss_t(k)                                                    &
                                      *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)              )    &
                                      / ( ABS( w_comp(k) ) + 1.0E-20_wp                       )    &
                                    ) *   weight_substep(intermediate_timestep_count)

          ENDDO

          DO  k = nzb_max_l+1, nzt

             flux_d    = flux_t(k-1)
             diss_d    = diss_t(k-1)
!
!--          Calculate the divergence of the velocity field. A respective correction is needed to
!--          overcome numerical instabilities caused by an insufficient reduction of divergences
!--          near topography.
             div = (  ( u_comp(k) + gu - ( u(k,j-1,i)  + u(k,j,i)   ) ) * ddx                      &
                   +  ( v_comp(k)      - ( v(k,j,i)    + v(k,j-1,i) ) ) * ddy                      &
                   +  ( w_comp(k)                     * rho_air_zw(k) -                            &
                        ( w(k-1,j-1,i) + w(k-1,j,i) ) * rho_air_zw(k-1)                            &
                      ) * drho_air(k) * ddzw(k)                                                    &
                   ) * 0.5_wp

             tend(k,j,i) = tend(k,j,i) - (                                                         &
                             ( ( flux_r(k) + diss_r(k) )                                           &
                           -   ( flux_l_v(k,j,tn) + diss_l_v(k,j,tn) ) ) * ddx                     &
                           + ( ( flux_n(k) + diss_n(k) )                                           &
                           -   ( flux_s_v(k,tn)   + diss_s_v(k,tn)   ) ) * ddy                     &
                           + ( ( flux_t(k) + diss_t(k) )                                           &
                           -   ( flux_d + diss_d ) ) * drho_air(k) * ddzw(k)                       &
                                         )  + v(k,j,i) * div

#ifndef _OPENACC
!
!--          Swap fluxes. Note, in the OPENACC case these are computed again.
             flux_l_v(k,j,tn) = flux_r(k)
             diss_l_v(k,j,tn) = diss_r(k)
             flux_s_v(k,tn)   = flux_n(k)
             diss_s_v(k,tn)   = diss_n(k)
#endif

!
!--          Statistical Evaluation of v'v'. The factor has to be applied for right evaluation when
!--          gallilei_trans = .T. .
             !$ACC ATOMIC
             sums_vs2_ws_l(k,tn) = sums_vs2_ws_l(k,tn) +                                           &
                                   ( flux_n(k)                                                     &
                                     * ( v_comp(k) - 2.0_wp * hom(k,1,2,0)                   )     &
                                     / ( v_comp(k) - gv + SIGN( 1.0E-20_wp, v_comp(k) - gv ) )     &
                                   +   diss_n(k)                                                   &
                                     *   ABS( v_comp(k) - 2.0_wp * hom(k,1,2,0)              )     &
                                     / ( ABS( v_comp(k) - gv ) + 1.0E-20_wp                  )     &
                                   ) *   weight_substep(intermediate_timestep_count)
!
!--          Statistical Evaluation of w'u'.
             !$ACC ATOMIC
             sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn) +                                         &
                                    ( flux_t(k)                                                    &
                                      * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)                   )    &
                                      / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )           )    &
                                    + diss_t(k)                                                    &
                                      *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)              )    &
                                      / ( ABS( w_comp(k) ) + 1.0E-20_wp                       )    &
                                    ) *  weight_substep(intermediate_timestep_count)

          ENDDO

       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(69), 'advec_v_ws', 'stop' )

 END SUBROUTINE advec_v_ws


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of w - Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_w_ws


    INTEGER(iwp) ::  i         !< grid index along x-direction
    INTEGER(iwp) ::  j         !< grid index along y-direction
    INTEGER(iwp) ::  k         !< grid index along z-direction
    INTEGER(iwp) ::  k_mm      !< k-2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_pp      !< k+2 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  k_ppp     !< k+3 index in disretization, can be modified to avoid segmentation faults
    INTEGER(iwp) ::  nzb_max_l !< index indicating upper bound for order degradation of horizontal advection terms
    INTEGER(iwp) ::  tn = 0    !< number of OpenMP thread

    REAL(wp)    ::  diss_d   !< artificial dissipation term at grid box bottom
    REAL(wp)    ::  div      !< divergence on w-grid
    REAL(wp)    ::  flux_d   !< 6th-order flux at grid box bottom
    REAL(wp)    ::  gu       !< Galilei-transformation velocity along x
    REAL(wp)    ::  gv       !< Galilei-transformation velocity along y
    REAL(wp)    ::  ibit18 !< flag indicating 1st-order scheme along x-direction
    REAL(wp)    ::  ibit19 !< flag indicating 3rd-order scheme along x-direction
    REAL(wp)    ::  ibit20 !< flag indicating 5th-order scheme along x-direction
#ifdef _OPENACC
    REAL(wp)    ::  ibit18_l !< flag indicating 1st-order scheme along x-direction
    REAL(wp)    ::  ibit19_l !< flag indicating 3rd-order scheme along x-direction
    REAL(wp)    ::  ibit20_l !< flag indicating 5th-order scheme along x-direction
#endif
    REAL(wp)    ::  ibit21 !< flag indicating 1st-order scheme along y-direction
    REAL(wp)    ::  ibit22 !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)    ::  ibit23 !< flag indicating 5th-order scheme along y-direction
#ifdef _OPENACC
    REAL(wp)    ::  ibit21_s !< flag indicating 1st-order scheme along y-direction
    REAL(wp)    ::  ibit22_s !< flag indicating 3rd-order scheme along y-direction
    REAL(wp)    ::  ibit23_s !< flag indicating 5th-order scheme along y-direction
#endif
    REAL(wp)    ::  ibit24   !< flag indicating 1st-order scheme along z-direction
    REAL(wp)    ::  ibit25   !< flag indicating 3rd-order scheme along z-direction
    REAL(wp)    ::  ibit26   !< flag indicating 5th-order scheme along z-direction
#ifdef _OPENACC
    REAL(wp)    ::  u_comp_l !< advection velocity along x
    REAL(wp)    ::  v_comp_s !< advection velocity along y
#endif

    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_n !< discretized artificial dissipation at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_r !< discretized artificial dissipation at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  diss_t !< discretized artificial dissipation at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_n !< discretized 6th-order flux at northward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_r !< discretized 6th-order flux at rightward-side of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  flux_t !< discretized 6th-order flux at top of the grid box
    REAL(wp), DIMENSION(nzb:nzt+1) ::  u_comp !< advection velocity along x
    REAL(wp), DIMENSION(nzb:nzt+1) ::  v_comp !< advection velocity along y
    REAL(wp), DIMENSION(nzb:nzt+1) ::  w_comp !< advection velocity along z


    CALL cpu_log( log_point_s(87), 'advec_w_ws', 'start' )

    gu = 2.0_wp * u_gtrans
    gv = 2.0_wp * v_gtrans

    !$ACC PARALLEL LOOP COLLAPSE(2) FIRSTPRIVATE(tn, gu, gv) &
    !$ACC PRIVATE(i, j, k, k_mm, k_pp, k_ppp) &
    !$ACC PRIVATE(ibit18, ibit19, ibit20, ibit21, ibit22, ibit23) &
    !$ACC PRIVATE(ibit24, ibit25, ibit26) &
    !$ACC PRIVATE(ibit18_l, ibit19_l, ibit20_l) &
    !$ACC PRIVATE(ibit21_s, ibit22_s, ibit23_s) &
    !$ACC PRIVATE(nzb_max_l) &
    !$ACC PRIVATE(flux_r, diss_r) &
    !$ACC PRIVATE(flux_n, diss_n) &
    !$ACC PRIVATE(flux_t, diss_t, flux_d, diss_d) &
    !$ACC PRIVATE(flux_l_w, diss_l_w, flux_s_w, diss_s_w) &
    !$ACC PRIVATE(div, u_comp, u_comp_l, v_comp, v_comp_s, w_comp) &
    !$ACC PRESENT(advc_flags_m) &
    !$ACC PRESENT(u, v, w) &
    !$ACC PRESENT(drho_air, rho_air_zw, ddzw) &
    !$ACC PRESENT(tend) &
    !$ACC PRESENT(hom(:,1,1:3,0)) &
    !$ACC PRESENT(weight_substep(intermediate_timestep_count)) &
    !$ACC PRESENT(sums_ws2_ws_l)
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Used local modified copy of nzb_max (used to degrade order of discretization) at
!--       non-cyclic boundaries. Modify only at relevant points instead of the entire subdomain.
!--       This should lead to better load balance between boundary and non-boundary PEs.
          IF( ( bc_dirichlet_l  .OR.  bc_radiation_l )  .AND.  i <= nxl + 2  .OR.                  &
              ( bc_dirichlet_r  .OR.  bc_radiation_r )  .AND.  i >= nxr - 2  .OR.                  &
              ( bc_dirichlet_s  .OR.  bc_radiation_s )  .AND.  j <= nys + 2  .OR.                  &
              ( bc_dirichlet_n  .OR.  bc_radiation_n )  .AND.  j >= nyn - 2 )  THEN
             nzb_max_l = nzt - 1
          ELSE
             nzb_max_l = nzb_max
          END IF

#ifndef _OPENACC
          IF ( i == nxl )  THEN
             DO  k = nzb+1, nzb_max_l
                ibit20 = REAL( IBITS(advc_flags_m(k,j,i-1),20,1), KIND = wp )
                ibit19 = REAL( IBITS(advc_flags_m(k,j,i-1),19,1), KIND = wp )
                ibit18 = REAL( IBITS(advc_flags_m(k,j,i-1),18,1), KIND = wp )

                u_comp(k)        = u(k+1,j,i) + u(k,j,i) - gu
                flux_l_w(k,j,tn) = u_comp(k) * (                                                   &
                                   (  37.0_wp * ibit20 * adv_mom_5                                 &
                                   +   7.0_wp * ibit19 * adv_mom_3                                 &
                                   +            ibit18 * adv_mom_1 ) * ( w(k,j,i)   + w(k,j,i-1) ) &
                                   - ( 8.0_wp * ibit20 * adv_mom_5                                 &
                                   +            ibit19 * adv_mom_3 ) * ( w(k,j,i+1) + w(k,j,i-2) ) &
                                   + (          ibit20 * adv_mom_5 ) * ( w(k,j,i+2) + w(k,j,i-3) ) &
                                               )

                diss_l_w(k,j,tn) = -ABS( u_comp(k) ) * (                                           &
                                   (  10.0_wp * ibit20 * adv_mom_5                                 &
                                   +   3.0_wp * ibit19 * adv_mom_3                                 &
                                   +            ibit18 * adv_mom_1 ) * ( w(k,j,i)   - w(k,j,i-1) ) &
                                   - ( 5.0_wp * ibit20 * adv_mom_5                                 &
                                   +            ibit19 * adv_mom_3 ) * ( w(k,j,i+1) - w(k,j,i-2) ) &
                                   + (          ibit20 * adv_mom_5 ) * ( w(k,j,i+2) - w(k,j,i-3) ) &
                                                       )

             ENDDO

             DO  k = nzb_max_l+1, nzt-1

                u_comp(k)        = u(k+1,j,i) + u(k,j,i) - gu
                flux_l_w(k,j,tn) = u_comp(k) * (                                                   &
                                     37.0_wp * ( w(k,j,i)   + w(k,j,i-1)   )                       &
                                   -  8.0_wp * ( w(k,j,i+1) + w(k,j,i-2) )                         &
                                   +           ( w(k,j,i+2) + w(k,j,i-3) ) ) * adv_mom_5
                diss_l_w(k,j,tn) = -ABS( u_comp(k) ) * (                                           &
                                     10.0_wp * ( w(k,j,i)   - w(k,j,i-1)   )                       &
                                   -  5.0_wp * ( w(k,j,i+1) - w(k,j,i-2) )                         &
                                   +           ( w(k,j,i+2) - w(k,j,i-3) ) ) * adv_mom_5

             ENDDO

          ENDIF

          IF ( j == nys )  THEN
             DO  k = nzb+1, nzb_max_l

                ibit23 = REAL( IBITS(advc_flags_m(k,j-1,i),23,1), KIND = wp )
                ibit22 = REAL( IBITS(advc_flags_m(k,j-1,i),22,1), KIND = wp )
                ibit21 = REAL( IBITS(advc_flags_m(k,j-1,i),21,1), KIND = wp )

                v_comp(k)      = v(k+1,j,i) + v(k,j,i) - gv
                flux_s_w(k,tn) = v_comp(k) * (                                                     &
                                 (  37.0_wp * ibit23 * adv_mom_5                                   &
                                 +   7.0_wp * ibit22 * adv_mom_3                                   &
                                 +            ibit21 * adv_mom_1 ) * ( w(k,j,i)   + w(k,j-1,i) )   &
                                 - ( 8.0_wp * ibit23 * adv_mom_5                                   &
                                 +            ibit22 * adv_mom_3 ) * ( w(k,j+1,i) + w(k,j-2,i) )   &
                                 + (          ibit23 * adv_mom_5 ) * ( w(k,j+2,i) + w(k,j-3,i) )   &
                                             )

                diss_s_w(k,tn) = -ABS( v_comp(k) ) * (                                             &
                                 (  10.0_wp * ibit23 * adv_mom_5                                   &
                                 +   3.0_wp * ibit22 * adv_mom_3                                   &
                                 +            ibit21 * adv_mom_1 ) * ( w(k,j,i)   - w(k,j-1,i) )   &
                                 - ( 5.0_wp * ibit23 * adv_mom_5                                   &
                                 +            ibit22 * adv_mom_3 ) * ( w(k,j+1,i) - w(k,j-2,i) )   &
                                 + (          ibit23 * adv_mom_5 ) * ( w(k,j+2,i) - w(k,j-3,i) )   &
                                                     )

             ENDDO

             DO  k = nzb_max_l+1, nzt-1

                v_comp(k)      = v(k+1,j,i) + v(k,j,i) - gv
                flux_s_w(k,tn) = v_comp(k) * (                                                     &
                                   37.0_wp * ( w(k,j,i) + w(k,j-1,i)   )                           &
                                 -  8.0_wp * ( w(k,j+1,i) +w(k,j-2,i)  )                           &
                                 +           ( w(k,j+2,i) + w(k,j-3,i) ) ) * adv_mom_5
                diss_s_w(k,tn) = -ABS( v_comp(k) ) * (                                             &
                                   10.0_wp * ( w(k,j,i) - w(k,j-1,i)   )                           &
                                 -  5.0_wp * ( w(k,j+1,i) - w(k,j-2,i) )                           &
                                 +           ( w(k,j+2,i) - w(k,j-3,i) ) ) * adv_mom_5

             ENDDO
          ENDIF
#endif
          DO  k = nzb+1, nzb_max_l

             ibit20 = REAL( IBITS(advc_flags_m(k,j,i),20,1), KIND = wp )
             ibit19 = REAL( IBITS(advc_flags_m(k,j,i),19,1), KIND = wp )
             ibit18 = REAL( IBITS(advc_flags_m(k,j,i),18,1), KIND = wp )

             u_comp(k) = u(k+1,j,i+1) + u(k,j,i+1) - gu
             flux_r(k) = u_comp(k) * (                                                             &
                         (  37.0_wp * ibit20 * adv_mom_5                                           &
                         +   7.0_wp * ibit19 * adv_mom_3                                           &
                         +            ibit18 * adv_mom_1 ) * ( w(k,j,i+1) + w(k,j,i)   )           &
                         - ( 8.0_wp * ibit20 * adv_mom_5                                           &
                         +            ibit19 * adv_mom_3 ) * ( w(k,j,i+2) + w(k,j,i-1) )           &
                         + (          ibit20 * adv_mom_5 ) * ( w(k,j,i+3) + w(k,j,i-2) )           &
                                     )

             diss_r(k) = -ABS( u_comp(k) ) * (                                                     &
                         (  10.0_wp * ibit20 * adv_mom_5                                           &
                         +   3.0_wp * ibit19 * adv_mom_3                                           &
                         +            ibit18 * adv_mom_1 ) * ( w(k,j,i+1) - w(k,j,i)  )            &
                         - ( 5.0_wp * ibit20 * adv_mom_5                                           &
                         +            ibit19 * adv_mom_3 ) * ( w(k,j,i+2) - w(k,j,i-1) )           &
                         + (          ibit20 * adv_mom_5 ) * ( w(k,j,i+3) - w(k,j,i-2) )           &
                                             )

#ifdef _OPENACC
!
!--          Recompute the left fluxes.
             ibit20_l = REAL( IBITS(advc_flags_m(k,j,i-1),20,1), KIND = wp )
             ibit19_l = REAL( IBITS(advc_flags_m(k,j,i-1),19,1), KIND = wp )
             ibit18_l = REAL( IBITS(advc_flags_m(k,j,i-1),18,1), KIND = wp )

             u_comp_l           = u(k+1,j,i) + u(k,j,i) - gu
             flux_l_w(k,j,tn)   = u_comp_l * (                                                     &
                                  (  37.0_wp * ibit20_l * adv_mom_5                                &
                                  +   7.0_wp * ibit19_l * adv_mom_3                                &
                                  +            ibit18_l * adv_mom_1 ) * ( w(k,j,i)   + w(k,j,i-1) )&
                                  - ( 8.0_wp * ibit20_l * adv_mom_5                                &
                                  +            ibit19_l * adv_mom_3 ) * ( w(k,j,i+1) + w(k,j,i-2) )&
                                  + (          ibit20_l * adv_mom_5 ) * ( w(k,j,i+2) + w(k,j,i-3) )&
                                             )

             diss_l_w(k,j,tn)   = -ABS( u_comp_l ) * (                                             &
                                  (  10.0_wp * ibit20_l * adv_mom_5                                &
                                  +   3.0_wp * ibit19_l * adv_mom_3                                &
                                  +            ibit18_l * adv_mom_1 ) * ( w(k,j,i)   - w(k,j,i-1) )&
                                  - ( 5.0_wp * ibit20_l * adv_mom_5                                &
                                  +            ibit19_l * adv_mom_3 ) * ( w(k,j,i+1) - w(k,j,i-2) )&
                                  + (          ibit20_l * adv_mom_5 ) * ( w(k,j,i+2) - w(k,j,i-3) )&
                                                     )
#endif


             ibit23 = REAL( IBITS(advc_flags_m(k,j,i),23,1), KIND = wp )
             ibit22 = REAL( IBITS(advc_flags_m(k,j,i),22,1), KIND = wp )
             ibit21 = REAL( IBITS(advc_flags_m(k,j,i),21,1), KIND = wp )

             v_comp(k) = v(k+1,j+1,i) + v(k,j+1,i) - gv
             flux_n(k) = v_comp(k) * (                                                             &
                         (  37.0_wp * ibit23 * adv_mom_5                                           &
                         +   7.0_wp * ibit22 * adv_mom_3                                           &
                         +            ibit21 * adv_mom_1 ) * ( w(k,j+1,i) + w(k,j,i)   )           &
                         - ( 8.0_wp * ibit23 * adv_mom_5                                           &
                         +            ibit22 * adv_mom_3 ) * ( w(k,j+2,i) + w(k,j-1,i) )           &
                         + (          ibit23 * adv_mom_5 ) * ( w(k,j+3,i) + w(k,j-2,i) )           &
                                     )

             diss_n(k) = -ABS( v_comp(k) ) * (                                                     &
                       (  10.0_wp * ibit23 * adv_mom_5                                             &
                       +   3.0_wp * ibit22 * adv_mom_3                                             &
                       +            ibit21 * adv_mom_1 ) * ( w(k,j+1,i) - w(k,j,i)  )              &
                       - ( 5.0_wp * ibit23 * adv_mom_5                                             &
                       +            ibit22 * adv_mom_3 ) * ( w(k,j+2,i) - w(k,j-1,i) )             &
                       + (          ibit23 * adv_mom_5 ) * ( w(k,j+3,i) - w(k,j-2,i) )             &
                                             )

#ifdef _OPENACC
!
!--          Recompute the south fluxes.
             ibit23_s = REAL( IBITS(advc_flags_m(k,j-1,i),23,1), KIND = wp )
             ibit22_s = REAL( IBITS(advc_flags_m(k,j-1,i),22,1), KIND = wp )
             ibit21_s = REAL( IBITS(advc_flags_m(k,j-1,i),21,1), KIND = wp )

             v_comp_s         = v(k+1,j,i) + v(k,j,i) - gv
             flux_s_w(k,tn)   = v_comp_s * (                                                       &
                                (  37.0_wp * ibit23_s * adv_mom_5                                  &
                                +   7.0_wp * ibit22_s * adv_mom_3                                  &
                                +            ibit21_s * adv_mom_1 ) * ( w(k,j,i)   + w(k,j-1,i) )  &
                                - ( 8.0_wp * ibit23_s * adv_mom_5                                  &
                                +            ibit22_s * adv_mom_3 ) * ( w(k,j+1,i) + w(k,j-2,i) )  &
                                + (          ibit23_s * adv_mom_5 ) * ( w(k,j+2,i) + w(k,j-3,i) )  &
                                           )

             diss_s_w(k,tn)   = -ABS( v_comp_s ) * (                                               &
                                (  10.0_wp * ibit23_s * adv_mom_5                                  &
                                +   3.0_wp * ibit22_s * adv_mom_3                                  &
                                +            ibit21_s * adv_mom_1 ) * ( w(k,j,i)   - w(k,j-1,i) )  &
                                - ( 5.0_wp * ibit23_s * adv_mom_5                                  &
                                +            ibit22_s * adv_mom_3 ) * ( w(k,j+1,i) - w(k,j-2,i) )  &
                                + (          ibit23_s * adv_mom_5 ) * ( w(k,j+2,i) - w(k,j-3,i) )  &
                                                   )
#endif
          ENDDO

          DO  k = nzb_max_l+1, nzt-1

             u_comp(k) = u(k+1,j,i+1) + u(k,j,i+1) - gu
             flux_r(k) = u_comp(k) * (                                                             &
                           37.0_wp * ( w(k,j,i+1) + w(k,j,i)   )                                   &
                         -  8.0_wp * ( w(k,j,i+2) + w(k,j,i-1) )                                   &
                         +           ( w(k,j,i+3) + w(k,j,i-2) ) ) * adv_mom_5

             diss_r(k) = -ABS( u_comp(k) ) * (                                                     &
                           10.0_wp * ( w(k,j,i+1) - w(k,j,i)   )                                   &
                         -  5.0_wp * ( w(k,j,i+2) - w(k,j,i-1) )                                   &
                         +           ( w(k,j,i+3) - w(k,j,i-2) ) ) * adv_mom_5

#ifdef _OPENACC
!
!--          Recompute the left fluxes.
             u_comp_l           = u(k+1,j,i) + u(k,j,i) - gu
             flux_l_w(k,j,tn)   = u_comp_l * (                                                     &
                                    37.0_wp * ( w(k,j,i)   + w(k,j,i-1)   )                        &
                                  -  8.0_wp * ( w(k,j,i+1) + w(k,j,i-2) )                          &
                                  +           ( w(k,j,i+2) + w(k,j,i-3) ) ) * adv_mom_5
             diss_l_w(k,j,tn)   = -ABS( u_comp_l ) * (                                             &
                                    10.0_wp * ( w(k,j,i)   - w(k,j,i-1)   )                        &
                                  -  5.0_wp * ( w(k,j,i+1) - w(k,j,i-2) )                          &
                                  +           ( w(k,j,i+2) - w(k,j,i-3) ) ) * adv_mom_5
#endif

             v_comp(k) = v(k+1,j+1,i) + v(k,j+1,i) - gv
             flux_n(k) = v_comp(k) * (                                                             &
                           37.0_wp * ( w(k,j+1,i) + w(k,j,i)   )                                   &
                         -  8.0_wp * ( w(k,j+2,i) + w(k,j-1,i) )                                   &
                         +           ( w(k,j+3,i) + w(k,j-2,i) ) ) * adv_mom_5

             diss_n(k) = -ABS( v_comp(k) ) * (                                                     &
                           10.0_wp * ( w(k,j+1,i) - w(k,j,i)   )                                   &
                         -  5.0_wp * ( w(k,j+2,i) - w(k,j-1,i) )                                   &
                         +           ( w(k,j+3,i) - w(k,j-2,i) ) ) * adv_mom_5

#ifdef _OPENACC
!
!--          Recompute the south fluxes.
             v_comp_s         = v(k+1,j,i) + v(k,j,i) - gv
             flux_s_w(k,tn)   = v_comp_s * (                                                       &
                                  37.0_wp * ( w(k,j,i)   + w(k,j-1,i) )                            &
                                -  8.0_wp * ( w(k,j+1,i) + w(k,j-2,i)  )                           &
                                +           ( w(k,j+2,i) + w(k,j-3,i) ) ) * adv_mom_5
             diss_s_w(k,tn)   = -ABS( v_comp_s ) * (                                               &
                                  10.0_wp * ( w(k,j,i)   - w(k,j-1,i) )                            &
                                -  5.0_wp * ( w(k,j+1,i) - w(k,j-2,i) )                            &
                                +           ( w(k,j+2,i) - w(k,j-3,i) ) ) * adv_mom_5
#endif
          ENDDO
!
!--       Now, compute vertical fluxes. Split loop into a part treating the lowest grid points with
!--       indirect indexing, a main loop without indirect indexing, and a loop for the uppermost
!--       grid points with indirect indexing. This allows better vectorization for the main loop.
!--       First, compute the flux at model surface, which need has to be calculated explicitly for
!--       the tendency at the first w-level. For topography wall this is done implicitely by
!--       advc_flags_m.
          k         = nzb + 1
          w_comp(k) = w(k,j,i) + w(k-1,j,i)
          flux_t(0) =  w_comp(k)        * rho_air(k) * ( w(k,j,i) + w(k-1,j,i) ) * adv_mom_1
          diss_t(0) = -ABS( w_comp(k) ) * rho_air(k) * ( w(k,j,i) - w(k-1,j,i) ) * adv_mom_1

          DO  k = nzb+1, nzb+1
!
!--          k index has to be modified near bottom and top, else array subscripts will be exceeded.
             ibit26 = REAL( IBITS(advc_flags_m(k,j,i),26,1), KIND = wp )
             ibit25 = REAL( IBITS(advc_flags_m(k,j,i),25,1), KIND = wp )
             ibit24 = REAL( IBITS(advc_flags_m(k,j,i),24,1), KIND = wp )

             k_ppp = k + 3 * ibit26
             k_pp  = k + 2 * ( 1 - ibit24  )
             k_mm  = k - 2 * ibit26

             w_comp(k) = w(k+1,j,i) + w(k,j,i)
             flux_t(k) = w_comp(k) * rho_air(k+1) * (                                              &
                        (  37.0_wp * ibit26 * adv_mom_5                                            &
                        +   7.0_wp * ibit25 * adv_mom_3                                            &
                        +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)   + w(k,j,i)    )         &
                        - ( 8.0_wp * ibit26 * adv_mom_5                                            &
                        +            ibit25 * adv_mom_3 ) * ( w(k_pp,j,i)  + w(k-1,j,i)  )         &
                        + (          ibit26 * adv_mom_5 ) * ( w(k_ppp,j,i) + w(k_mm,j,i) )         &
                                                    )

             diss_t(k) = -ABS( w_comp(k) ) * rho_air(k+1) * (                                      &
                         (  10.0_wp * ibit26 * adv_mom_5                                           &
                         +   3.0_wp * ibit25 * adv_mom_3                                           &
                         +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)   - w(k,j,i)    )        &
                         - ( 5.0_wp * ibit26 * adv_mom_5                                           &
                         +            ibit25 * adv_mom_3 ) * ( w(k_pp,j,i)  - w(k-1,j,i)  )        &
                         + (          ibit26 * adv_mom_5 ) * ( w(k_ppp,j,i) - w(k_mm,j,i) )        &
                                                            )
          ENDDO

          DO  k = nzb+2, nzt-2

             ibit26 = REAL( IBITS(advc_flags_m(k,j,i),26,1), KIND = wp )
             ibit25 = REAL( IBITS(advc_flags_m(k,j,i),25,1), KIND = wp )
             ibit24 = REAL( IBITS(advc_flags_m(k,j,i),24,1), KIND = wp )

             w_comp(k) = w(k+1,j,i) + w(k,j,i)
             flux_t(k) = w_comp(k) * rho_air(k+1) * (                                              &
                        (  37.0_wp * ibit26 * adv_mom_5                                            &
                        +   7.0_wp * ibit25 * adv_mom_3                                            &
                        +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)  + w(k,j,i)   )           &
                        - ( 8.0_wp * ibit26 * adv_mom_5                                            &
                        +            ibit25 * adv_mom_3 ) * ( w(k+2,j,i)  + w(k-1,j,i) )           &
                        + (          ibit26 * adv_mom_5 ) * ( w(k+3,j,i)  + w(k-2,j,i) )           &
                                                    )

             diss_t(k) = -ABS( w_comp(k) ) * rho_air(k+1) * (                                      &
                         (  10.0_wp * ibit26 * adv_mom_5                                           &
                         +   3.0_wp * ibit25 * adv_mom_3                                           &
                         +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i) - w(k,j,i)    )          &
                         - ( 5.0_wp * ibit26 * adv_mom_5                                           &
                         +            ibit25 * adv_mom_3 ) * ( w(k+2,j,i) - w(k-1,j,i)  )          &
                         + (          ibit26 * adv_mom_5 ) * ( w(k+3,j,i) - w(k-2,j,i)  )          &
                                                            )
          ENDDO

          DO  k = nzt-1, nzt-1
!
!--          k index has to be modified near bottom and top, else array subscripts will be exceeded.
             ibit26 = REAL( IBITS(advc_flags_m(k,j,i),26,1), KIND = wp )
             ibit25 = REAL( IBITS(advc_flags_m(k,j,i),25,1), KIND = wp )
             ibit24 = REAL( IBITS(advc_flags_m(k,j,i),24,1), KIND = wp )
!
!--          k index needs to be modified near bottom and top, else array subscripts will be
!--          exceeded. (k_ppp = k+3, k_mm = k) if 5th order scheme is employed, else
!--          (k_ppp = k, k_mm = k). Special treatment is required for k_pp. This must be limited
!--          at k = nzt, where either the first order scheme is employed (ibit24), or all bit flags
!--          are zero (in case there is a rigid lid defined at or near the domain top).
             k_ppp = k + 3 * ibit26
             k_pp  = k + 2 * ( 1 - ibit24  ) * ( ibit24 + ibit25 + ibit26 )
             k_mm  = k - 2 * ibit26

             w_comp(k) = w(k+1,j,i) + w(k,j,i)
             flux_t(k) = w_comp(k) * rho_air(k+1) * (                                              &
                         (  37.0_wp * ibit26 * adv_mom_5                                           &
                         +   7.0_wp * ibit25 * adv_mom_3                                           &
                         +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)   + w(k,j,i)    )        &
                         - ( 8.0_wp * ibit26 * adv_mom_5                                           &
                         +            ibit25 * adv_mom_3 ) * ( w(k_pp,j,i)  + w(k-1,j,i)  )        &
                         + (          ibit26 * adv_mom_5 ) * ( w(k_ppp,j,i) + w(k_mm,j,i) )        &
                                                    )

             diss_t(k) = -ABS( w_comp(k) ) * rho_air(k+1) * (                                      &
                         (  10.0_wp * ibit26 * adv_mom_5                                           &
                         +   3.0_wp * ibit25 * adv_mom_3                                           &
                         +            ibit24 * adv_mom_1 ) * ( w(k+1,j,i)   - w(k,j,i)    )        &
                         - ( 5.0_wp * ibit26 * adv_mom_5                                           &
                         +            ibit25 * adv_mom_3 ) * ( w(k_pp,j,i)  - w(k-1,j,i)  )        &
                         + (          ibit26 * adv_mom_5 ) * ( w(k_ppp,j,i) - w(k_mm,j,i) )        &
                                                            )
          ENDDO

!
!--       Set resolved/turbulent flux at model top to zero (w-level). Hint: The flux at nzt is
!--       defined at the scalar grid point nzt+1. Therefore, the flux at nzt+1 is already outside of
!--       the model domain
          flux_t(nzt) = 0.0_wp
          diss_t(nzt) = 0.0_wp
          w_comp(nzt) = 0.0_wp

          flux_t(nzt+1) = 0.0_wp
          diss_t(nzt+1) = 0.0_wp
          w_comp(nzt+1) = 0.0_wp

          DO  k = nzb+1, nzb_max_l

             flux_d    = flux_t(k-1)
             diss_d    = diss_t(k-1)

             ibit20 = REAL( IBITS(advc_flags_m(k,j,i),20,1), KIND = wp )
             ibit19 = REAL( IBITS(advc_flags_m(k,j,i),19,1), KIND = wp )
             ibit18 = REAL( IBITS(advc_flags_m(k,j,i),18,1), KIND = wp )

             ibit23 = REAL( IBITS(advc_flags_m(k,j,i),23,1), KIND = wp )
             ibit22 = REAL( IBITS(advc_flags_m(k,j,i),22,1), KIND = wp )
             ibit21 = REAL( IBITS(advc_flags_m(k,j,i),21,1), KIND = wp )

             ibit26 = REAL( IBITS(advc_flags_m(k,j,i),26,1), KIND = wp )
             ibit25 = REAL( IBITS(advc_flags_m(k,j,i),25,1), KIND = wp )
             ibit24 = REAL( IBITS(advc_flags_m(k,j,i),24,1), KIND = wp )
!
!--          Calculate the divergence of the velocity field. A respective correction is needed to
!--          overcome numerical instabilities introduced by an insufficient reduction of divergences
!--          near topography.
             div = ( ( ( u_comp(k) + gu ) * ( ibit18 + ibit19 + ibit20 )                           &
                     - ( u(k+1,j,i) + u(k,j,i)   )                                                 &
                                       * (                                                         &
                        REAL( IBITS(advc_flags_m(k,j,i-1),18,1), KIND = wp )                       &
                      + REAL( IBITS(advc_flags_m(k,j,i-1),19,1), KIND = wp )                       &
                      + REAL( IBITS(advc_flags_m(k,j,i-1),20,1), KIND = wp )                       &
                                         )                                                         &
                     ) * ddx                                                                       &
                   + ( ( v_comp(k) + gv ) * ( ibit21 + ibit22 + ibit23 )                           &
                     - ( v(k+1,j,i) + v(k,j,i)   )                                                 &
                                       * (                                                         &
                        REAL( IBITS(advc_flags_m(k,j-1,i),21,1), KIND = wp )                       &
                      + REAL( IBITS(advc_flags_m(k,j-1,i),22,1), KIND = wp )                       &
                      + REAL( IBITS(advc_flags_m(k,j-1,i),23,1), KIND = wp )                       &
                                         )                                                         &
                     ) * ddy                                                                       &
                   + ( w_comp(k)               * rho_air(k+1)                                      &
                                               * ( ibit24 + ibit25 + ibit26 )                      &
                   - ( w(k,j,i) + w(k-1,j,i) ) * rho_air(k)                                        &
                                       * (                                                         &
                        REAL( IBITS(advc_flags_m(k-1,j,i),24,1), KIND = wp )                       &
                      + REAL( IBITS(advc_flags_m(k-1,j,i),25,1), KIND = wp )                       &
                      + REAL( IBITS(advc_flags_m(k-1,j,i),26,1), KIND = wp )                       &
                                         )                                                         &
                     ) * drho_air_zw(k) * ddzu(k+1)                                                &
                   ) * 0.5_wp

             tend(k,j,i) = tend(k,j,i) - (                                                         &
                              ( flux_r(k) + diss_r(k)                                              &
                            -   flux_l_w(k,j,tn) - diss_l_w(k,j,tn)   ) * ddx                      &
                            + ( flux_n(k) + diss_n(k)                                              &
                            -   flux_s_w(k,tn) - diss_s_w(k,tn)       ) * ddy                      &
                            + ( ( flux_t(k) + diss_t(k) )                                          &
                            -   ( flux_d    + diss_d    ) ) * drho_air_zw(k) * ddzu(k+1)           &
                                         ) + div * w(k,j,i)
#ifndef _OPENACC
             flux_l_w(k,j,tn) = flux_r(k)
             diss_l_w(k,j,tn) = diss_r(k)
             flux_s_w(k,tn)   = flux_n(k)
             diss_s_w(k,tn)   = diss_n(k)
#endif
!
!--          Statistical Evaluation of w'w'.
             sums_ws2_ws_l(k,tn)  = sums_ws2_ws_l(k,tn) +                                          &
                                    ( flux_t(k)                                                    &
                                      * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)           )            &
                                      / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )   )            &
                                    + diss_t(k)                                                    &
                                      *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)      )            &
                                      / ( ABS( w_comp(k) ) + 1.0E-20_wp               )            &
                                    ) *   weight_substep(intermediate_timestep_count)

          ENDDO

          DO  k = nzb_max_l+1, nzt-1

             flux_d    = flux_t(k-1)
             diss_d    = diss_t(k-1)
!
!--          Calculate the divergence of the velocity field. A respective correction is needed to
!--          overcome numerical instabilities introduced by an insufficient reduction of divergences
!--          near topography.
             div = ( ( u_comp(k) + gu - ( u(k+1,j,i) + u(k,j,i)   ) ) * ddx                        &
                   + ( v_comp(k) + gv - ( v(k+1,j,i) + v(k,j,i)   ) ) * ddy                        &
                   + ( w_comp(k)               * rho_air(k+1)                                      &
                   - ( w(k,j,i) + w(k-1,j,i) ) * rho_air(k)                                        &
                     ) * drho_air_zw(k) * ddzu(k+1)                                                &
                   ) * 0.5_wp

             tend(k,j,i) = tend(k,j,i) - (                                                         &
                             ( flux_r(k) + diss_r(k)                                               &
                           -   flux_l_w(k,j,tn) - diss_l_w(k,j,tn)   ) * ddx                       &
                           + ( flux_n(k) + diss_n(k)                                               &
                           -   flux_s_w(k,tn) - diss_s_w(k,tn)       ) * ddy                       &
                           + ( ( flux_t(k) + diss_t(k) )                                           &
                           -   ( flux_d    + diss_d    ) ) * drho_air_zw(k) * ddzu(k+1)            &
                                         ) + div * w(k,j,i)
#ifndef _OPENACC
             flux_l_w(k,j,tn) = flux_r(k)
             diss_l_w(k,j,tn) = diss_r(k)
             flux_s_w(k,tn)   = flux_n(k)
             diss_s_w(k,tn)   = diss_n(k)
#endif
!
!--          Statistical Evaluation of w'w'.
             sums_ws2_ws_l(k,tn)  = sums_ws2_ws_l(k,tn) +                                          &
                                    ( flux_t(k)                                                    &
                                      * ( w_comp(k) - 2.0_wp * hom(k,1,3,0)          )             &
                                      / ( w_comp(k) + SIGN( 1.0E-20_wp, w_comp(k) )  )             &
                                    + diss_t(k)                                                    &
                                      *   ABS( w_comp(k) - 2.0_wp * hom(k,1,3,0)     )             &
                                      / ( ABS( w_comp(k) ) + 1.0E-20_wp              )             &
                                    ) *   weight_substep(intermediate_timestep_count)

          ENDDO

       ENDDO
    ENDDO

    CALL cpu_log( log_point_s(87), 'advec_w_ws', 'stop' )

 END SUBROUTINE advec_w_ws

 END MODULE advec_ws
