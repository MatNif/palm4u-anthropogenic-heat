!> @file diffusion_v.f90
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
! Description:
! ------------
!> Diffusion term of the v-component
!--------------------------------------------------------------------------------------------------!
 MODULE diffusion_v_mod

    USE arrays_3d,                                                                                 &
        ONLY:  ddzu,                                                                               &
               ddzw,                                                                               &
               drho_air,                                                                           &
               dzw,                                                                                &
               km,                                                                                 &
               tend,                                                                               &
               u,                                                                                  &
               v,                                                                                  &
               w,                                                                                  &
               rho_air_zw

    USE control_parameters,                                                                        &
        ONLY:  constant_top_momentumflux,                                                          &
               use_surface_fluxes,                                                                 &
               use_top_fluxes

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nysv,                                                                               &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_flags

    USE kinds

    USE surface_mod,                                                                               &
        ONLY:  surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_top,                                                                           &
               surf_usm

    PRIVATE
    PUBLIC diffusion_v

    INTERFACE diffusion_v
       MODULE PROCEDURE diffusion_v
       MODULE PROCEDURE diffusion_v_ij
    END INTERFACE diffusion_v

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_v

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  kmxm          !< diffusion coefficient on leftward side of the v-gridbox - interpolated onto yv-yv grid
       REAL(wp) ::  kmxp          !< diffusion coefficient on rightward side of the v-gridbox - interpolated onto yv-yv grid
       REAL(wp) ::  kmzm          !< diffusion coefficient on bottom of the gridbox - interpolated onto yv-zw grid
       REAL(wp) ::  kmzp          !< diffusion coefficient on top of the gridbox - interpolated onto yv-zw grid
       REAL(wp) ::  mask_east     !< flag to mask vertical surface east of the grid point
       REAL(wp) ::  mask_west     !< flag to mask vertical surface west of the grid point

       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_d   !< flux at k-1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_t   !< flux at k+1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_l   !< flux at i-1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_r   !< flux at i+1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_n   !< flux at j+1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_s   !< flux at j-1/2 (from v-grid point)

       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
       !$ACC PRIVATE(surf_e, surf_s, flag, kmxm, kmxp, kmzm, kmzp) &
       !$ACC PRIVATE(flux_r(nzb+1:nzt), flux_l(nzb+1:nzt), flux_n(nzb+1:nzt), flux_s(nzb+1:nzt), flux_t(nzb+1:nzt), flux_d(nzb+1:nzt)) &
       !$ACC PRIVATE(mask_east, mask_west) &
       !$ACC PRESENT(topo_flags, km) &
       !$ACC PRESENT(u, v, w) &
       !$ACC PRESENT(ddzu, ddzw, drho_air, rho_air_zw) &
       !$ACC PRESENT(surf_def, surf_lsm, surf_top, surf_usm) &
       !$ACC PRESENT(tend)
       DO  i = nxl, nxr
          DO  j = nysv, nyn
!
!--          Compute diffusion fluxes at the grid-cell boundaries using deformation approach.
!--          Note, at wall-bounded grid points the respective fluxes will be overwritten later.
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt
                mask_east = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 2 ) )
                mask_west = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 2 ) )
!
!--             Interpolate eddy diffusivities onto staggered grid
                kmxp = 0.25_wp * ( km(k,j,i) + km(k,j,i+1) + km(k,j-1,i) + km(k,j-1,i+1) )
                kmxm = 0.25_wp * ( km(k,j,i) + km(k,j,i-1) + km(k,j-1,i) + km(k,j-1,i-1) )
                kmzp = 0.25_wp * ( km(k,j,i) + km(k+1,j,i) + km(k,j-1,i) + km(k+1,j-1,i) )
                kmzm = 0.25_wp * ( km(k,j,i) + km(k-1,j,i) + km(k,j-1,i) + km(k-1,j-1,i) )
!
!--             Mask diffusion at grid points bounded by east- or west facing walls. For the moment
!--             this is necessary especially at corner grid points (at eastward facing walls). There,
!--             the v-grid point is actually still bounded by a wall (topo_flags, bit 2 is zero),
!--             but no surface is defined at this point (v(k,j,i) is topography, but s(k,j,i) is
!--             atmosphere). As a result, no surface flux is added, meaning that the momentum sink
!--             is missing at these corner gird points. In order to add at least any momentum sink,
!--             the masking of these points should be avoided. Please note, this does not have any
!--             effect on the diffusion at east- or west-facing walls which are no corner grid
!--             points. This is because the diffusive fluxes will be replaced by surface fluxes
!--             further below.
                flux_r(k) = - kmxp * ( ( v(k,j,i+1) - v(k,j,i)     ) * ddx                         &
                                     + ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy                         &
                                     ) * mask_east
                flux_l(k) = - kmxm * ( ( v(k,j,i)   - v(k,j,i-1)   ) * ddx                         &
                                     + ( u(k,j,i)   - u(k,j-1,i)   ) * ddy                         &
                                     ) * mask_west
                flux_n(k) = - 2.0_wp * km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i)   ) * ddy
                flux_s(k) = - 2.0_wp * km(k,j-1,i) * ( v(k,j,i)   - v(k,j-1,i) ) * ddy
                flux_t(k) = - kmzp * ( ( v(k+1,j,i) - v(k,j,i)     ) * ddzu(k+1)                   &
                                     + ( w(k,j,i)   - w(k,j-1,i)   ) * ddy                         &
                                     ) * rho_air_zw(k)
                flux_d(k) = - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)                     &
                                     + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy                         &
                                     ) * rho_air_zw(k-1)

             ENDDO
!
!--          If surface fluxes should be used, overwrite the fluxes. Therefore, first determine
!--          start and end indices for the surface elements stored for grid point (j,i). If there
!--          are no wall grid points at given (j,i), loops won't be entered. Further, the surface
!--          fluxes are also multiplied with the relevant normal vector component, which includes
!--          the sign of the flux as well as weights the flux according to its surface orientation
!--          in case of slanted walls. For the v-component tendency, v'u' is added at walls
!--          that are bounded along the east-/ west direction. No fluxes are added along the
!--          north-/south direction because these wouldn't result from any shear stress but from
!--          momentum of the wall itself (which is 0).
             IF ( use_surface_fluxes )  THEN
!
!--             Default-type surfaces
                surf_s = surf_def%start_index(j,i)
                surf_e = surf_def%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k           = surf_def%k(m)
                   flux_r(k)   = MERGE( surf_def%vsus(m), flux_r(k), surf_def%westward(m) )
                   flux_l(k)   = MERGE( surf_def%vsus(m), flux_l(k), surf_def%eastward(m) )
                   flux_t(k)   = MERGE( surf_def%vsws(m), flux_t(k), surf_def%downward(m) )
                   flux_d(k)   = MERGE( surf_def%vsws(m), flux_d(k), surf_def%upward(m)   )
                ENDDO
!
!--             Natural-type surfaces
                surf_s = surf_lsm%start_index(j,i)
                surf_e = surf_lsm%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k           = surf_lsm%k(m)
                   flux_r(k)   = MERGE( surf_lsm%vsus(m), flux_r(k), surf_lsm%westward(m) )
                   flux_l(k)   = MERGE( surf_lsm%vsus(m), flux_l(k), surf_lsm%eastward(m) )
                   flux_t(k)   = MERGE( surf_lsm%vsws(m), flux_t(k), surf_lsm%downward(m) )
                   flux_d(k)   = MERGE( surf_lsm%vsws(m), flux_d(k), surf_lsm%upward(m)   )
                ENDDO
!
!--             Urban-type surfaces
                surf_s = surf_usm%start_index(j,i)
                surf_e = surf_usm%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k           = surf_usm%k(m)
                   flux_r(k)   = MERGE( surf_usm%vsus(m), flux_r(k), surf_usm%westward(m) )
                   flux_l(k)   = MERGE( surf_usm%vsus(m), flux_l(k), surf_usm%eastward(m) )
                   flux_t(k)   = MERGE( surf_usm%vsws(m), flux_t(k), surf_usm%downward(m) )
                   flux_d(k)   = MERGE( surf_usm%vsws(m), flux_d(k), surf_usm%upward(m)   )
                ENDDO
             ENDIF
!
!--          Add momentum flux at model top
             IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
                surf_s = surf_top%start_index(j,i)
                surf_e = surf_top%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k         = surf_top%k(m)
                   flux_t(k) = surf_top%vsws(m)
                ENDDO
             ENDIF
!
!--          Compute tendency. Flag 2 is used to mask topography on v-grid.
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt
                flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                tend(k,j,i) = tend(k,j,i) - (   ( flux_r(k) - flux_l(k) ) * ddx                    &
                                              + ( flux_n(k) - flux_s(k) ) * ddy                    &
                                              + ( flux_t(k) - flux_d(k) ) * ddzw(k) * drho_air(k)  &
                                            ) * flag
             ENDDO

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_v


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_v_ij( i, j )

       IMPLICIT NONE


       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  kmxm          !< diffusion coefficient on leftward side of the v-gridbox - interpolated onto yv-yv grid
       REAL(wp) ::  kmxp          !< diffusion coefficient on rightward side of the v-gridbox - interpolated onto yv-yv grid
       REAL(wp) ::  kmzm          !< diffusion coefficient on bottom of the gridbox - interpolated onto yv-zw grid
       REAL(wp) ::  kmzp          !< diffusion coefficient on top of the gridbox - interpolated onto yv-zw grid
       REAL(wp) ::  mask_east     !< flag to mask vertical surface east of the grid point
       REAL(wp) ::  mask_west     !< flag to mask vertical surface west of the grid point

       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_d   !< flux at k-1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_t   !< flux at k+1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_l   !< flux at i-1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_r   !< flux at i+1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_n   !< flux at j+1/2 (from v-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_s   !< flux at j-1/2 (from v-grid point)
!
!--    Compute diffusion fluxes at the grid-cell boundaries using deformation approach.
!--    Note, at wall-bounded grid points the respective fluxes will be overwritten later.
       DO  k = nzb+1, nzt

          mask_east = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 2 ) )
          mask_west = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 2 ) )
!
!--       Interpolate eddy diffusivities onto staggered grid
          kmxp = 0.25_wp * ( km(k,j,i) + km(k,j,i+1) + km(k,j-1,i) + km(k,j-1,i+1) )
          kmxm = 0.25_wp * ( km(k,j,i) + km(k,j,i-1) + km(k,j-1,i) + km(k,j-1,i-1) )
          kmzp = 0.25_wp * ( km(k,j,i) + km(k+1,j,i) + km(k,j-1,i) + km(k+1,j-1,i) )
          kmzm = 0.25_wp * ( km(k,j,i) + km(k-1,j,i) + km(k,j-1,i) + km(k-1,j-1,i) )
!
!--       Mask diffusion at grid points bounded by east- or west facing walls. For the moment
!--       this is necessary especially at corner grid points (at eastward facing walls). There,
!--       the v-grid point is actually still bounded by a wall (topo_flags, bit 2 is zero), but
!--       no surface is defined at this point (v(k,j,i) is topography, but s(k,j,i) is atmosphere).
!--       As a result, no surface flux is added, meaning that the momentum sink is missing
!--       at these corner gird points. In order to add at least any momentum sink, the masking
!--       of these points should be avoided. Please note, this does not have any effect on the
!--       diffusion at east- or west-facing walls which are no corner grid points. This is
!--       because the diffusive fluxes will be replaced by surface fluxes further below.
          flux_r(k) = - kmxp * ( ( v(k,j,i+1) - v(k,j,i)     ) * ddx                               &
                               + ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy                               &
                               ) * mask_east
          flux_l(k) = - kmxm * ( ( v(k,j,i)   - v(k,j,i-1)   ) * ddx                               &
                               + ( u(k,j,i)   - u(k,j-1,i)   ) * ddy                               &
                               ) * mask_west
          flux_n(k) = - 2.0_wp * km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i)   ) * ddy
          flux_s(k) = - 2.0_wp * km(k,j-1,i) * ( v(k,j,i)   - v(k,j-1,i) ) * ddy
          flux_t(k) = - kmzp * ( ( v(k+1,j,i) - v(k,j,i)     ) * ddzu(k+1)                         &
                               + ( w(k,j,i)   - w(k,j-1,i)   ) * ddy                               &
                               ) * rho_air_zw(k)
          flux_d(k) = - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)                           &
                               + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy                               &
                               ) * rho_air_zw(k-1)
       ENDDO

!
!--    If surface fluxes should be used, overwrite the fluxes. Therefore, first determine start and
!--    end indices for the surface elements stored for grid point (j,i). If there are no wall grid
!--    points at given (j,i), loops won't be entered. Further, the surface fluxes are also
!--    multiplied with the relevant normal vector component, which includes the sign of the flux.
!--    For the v-component tendency, v'u' is added at walls that are bounded along the east-/
!--    west direction. No fluxes are added along the north-/south direction because these wouldn't
!--    result from any shear stress but from momentum of the wall itself (which is 0).
       IF ( use_surface_fluxes )  THEN
!
!--       Default-type surfaces
          surf_s = surf_def%start_index(j,i)
          surf_e = surf_def%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def%k(m)
             flux_r(k)   = MERGE( surf_def%vsus(m), flux_r(k), surf_def%westward(m) )
             flux_l(k)   = MERGE( surf_def%vsus(m), flux_l(k), surf_def%eastward(m) )
             flux_t(k)   = MERGE( surf_def%vsws(m), flux_t(k), surf_def%downward(m) )
             flux_d(k)   = MERGE( surf_def%vsws(m), flux_d(k), surf_def%upward(m)   )
          ENDDO
!
!--       Natural-type surfaces
          surf_s = surf_lsm%start_index(j,i)
          surf_e = surf_lsm%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_lsm%k(m)
             flux_r(k)   = MERGE( surf_lsm%vsus(m), flux_r(k), surf_lsm%westward(m) )
             flux_l(k)   = MERGE( surf_lsm%vsus(m), flux_l(k), surf_lsm%eastward(m) )
             flux_t(k)   = MERGE( surf_lsm%vsws(m), flux_t(k), surf_lsm%downward(m) )
             flux_d(k)   = MERGE( surf_lsm%vsws(m), flux_d(k), surf_lsm%upward(m)   )
          ENDDO
!
!--       Urban-type surfaces
          surf_s = surf_usm%start_index(j,i)
          surf_e = surf_usm%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_usm%k(m)
             flux_r(k)   = MERGE( surf_usm%vsus(m), flux_r(k), surf_usm%westward(m) )
             flux_l(k)   = MERGE( surf_usm%vsus(m), flux_l(k), surf_usm%eastward(m) )
             flux_t(k)   = MERGE( surf_usm%vsws(m), flux_t(k), surf_usm%downward(m) )
             flux_d(k)   = MERGE( surf_usm%vsws(m), flux_d(k), surf_usm%upward(m)   )
          ENDDO
       ENDIF
!
!--    Add momentum flux at model top
       IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
          surf_s = surf_top%start_index(j,i)
          surf_e = surf_top%end_index(j,i)
          DO  m = surf_s, surf_e
             k         = surf_top%k(m)
             flux_t(k) = surf_top%vsws(m)
          ENDDO
       ENDIF
!
!--    Compute tendency. Flag 2 is used to mask topography on v-grid.
       DO  k = nzb+1, nzt
          flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          tend(k,j,i) = tend(k,j,i) - (   ( flux_r(k) - flux_l(k) ) * ddx                          &
                                        + ( flux_n(k) - flux_s(k) ) * ddy                          &
                                        + ( flux_t(k) - flux_d(k) ) * ddzw(k) * drho_air(k)        &
                                      ) * flag
       ENDDO

    END SUBROUTINE diffusion_v_ij

 END MODULE diffusion_v_mod
