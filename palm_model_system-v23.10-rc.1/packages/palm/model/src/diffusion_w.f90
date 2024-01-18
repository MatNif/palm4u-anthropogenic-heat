!> @file diffusion_w.f90
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
!> Diffusion term of the w-component
!--------------------------------------------------------------------------------------------------!
 MODULE diffusion_w_mod

    USE arrays_3d,                                                                                 &
        ONLY:  ddzu,                                                                               &
               ddzw,                                                                               &
               drho_air_zw,                                                                        &
               dzu,                                                                                &
               km,                                                                                 &
               rho_air,                                                                            &
               tend,                                                                               &
               u,                                                                                  &
               v,                                                                                  &
               w

    USE control_parameters,                                                                        &
        ONLY:  use_surface_fluxes

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_flags

    USE kinds

    USE surface_mod,                                                                               &
        ONLY:  surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_usm

    PRIVATE
    PUBLIC diffusion_w

    INTERFACE diffusion_w
       MODULE PROCEDURE diffusion_w
       MODULE PROCEDURE diffusion_w_ij
    END INTERFACE diffusion_w

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_w

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  kmxm              !< diffusion coefficient on leftward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmxp              !<diffusion coefficient on rightward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmym              !< diffusion coefficient on southward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  kmyp              !< diffusion coefficient on northward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  mask_east         !< flag to mask vertical wall east of the grid point
       REAL(wp) ::  mask_north        !< flag to mask vertical wall north of the grid point
       REAL(wp) ::  mask_south        !< flag to mask vertical wall south of the grid point
       REAL(wp) ::  mask_west         !< flag to mask vertical wall west of the grid point

       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_d   !< flux at k-1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_t   !< flux at k+1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_l   !< flux at i-1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_r   !< flux at i+1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_n   !< flux at j+1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_s   !< flux at j-1/2 (from w-grid point)

       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
       !$ACC PRIVATE(surf_e, surf_s, flag, kmxm, kmxp, kmym, kmyp) &
       !$ACC PRIVATE(flux_r(nzb+1:nzt), flux_l(nzb+1:nzt), flux_n(nzb+1:nzt), flux_s(nzb+1:nzt), flux_t(nzb+1:nzt), flux_d(nzb+1:nzt)) &
       !$ACC PRIVATE(mask_west, mask_east, mask_south, mask_north) &
       !$ACC PRESENT(topo_flags, km) &
       !$ACC PRESENT(u, v, w) &
       !$ACC PRESENT(ddzu, ddzw, rho_air, drho_air_zw) &
       !$ACC PRESENT(surf_def, surf_lsm, surf_usm) &
       !$ACC PRESENT(tend)
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Compute diffusion fluxes at the grid-cell boundaries using deformation approach.
!--          Note, at wall-bounded grid points the respective fluxes will be overwritten later.
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt-1
!
!--             Predetermine flags mask wall-bounded grid points. For the
!--             moment this is necessary especially at the upper corner grid points.
!--             There, the w-grid point is actually still bounded by a wall
!--             (topo_flags, bit 3 is zero), but no surface is defined at this point (w(k,j,i)
!--             is topography, but s(k,j,i) is atmosphere).
                mask_east  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 3 ) )
                mask_west  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 3 ) )
                mask_south = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j-1,i), 3 ) )
                mask_north = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j+1,i), 3 ) )
!
!--             Interpolate eddy diffusivities onto staggered grid
                kmxp = 0.25_wp * ( km(k,j,i) + km(k,j,i+1) + km(k+1,j,i) + km(k+1,j,i+1) )
                kmxm = 0.25_wp * ( km(k,j,i) + km(k,j,i-1) + km(k+1,j,i) + km(k+1,j,i-1) )
                kmyp = 0.25_wp * ( km(k,j,i) + km(k+1,j,i) + km(k,j+1,i) + km(k+1,j+1,i) )
                kmym = 0.25_wp * ( km(k,j,i) + km(k+1,j,i) + km(k,j-1,i) + km(k+1,j-1,i) )

                flux_r(k) = - kmxp * ( ( w(k,j,i+1)   - w(k,j,i)   ) * ddx                         &
                                     + ( u(k+1,j,i+1) - u(k,j,i+1) ) * ddzu(k+1)                   &
                                     ) * mask_east
                flux_l(k) = - kmxm * ( ( w(k,j,i)     - w(k,j,i-1) ) * ddx                         &
                                     + ( u(k+1,j,i)   - u(k,j,i)   ) * ddzu(k+1)                   &
                                     ) * mask_west
                flux_n(k) = - kmyp * ( ( w(k,j+1,i)   - w(k,j,i)   ) * ddy                         &
                                     + ( v(k+1,j+1,i) - v(k,j+1,i) ) * ddzu(k+1)                   &
                                     ) * mask_north
                flux_s(k) = - kmym * ( ( w(k,j,i)     - w(k,j-1,i) ) * ddy                         &
                                     + ( v(k+1,j,i)   - v(k,j,i)   ) * ddzu(k+1)                   &
                                     ) * mask_south
                flux_t(k) = - 2.0_wp * km(k+1,j,i) * ( w(k+1,j,i) - w(k,j,i)   ) * ddzw(k+1)       &
                                                                                 * rho_air(k+1)
                flux_d(k) = - 2.0_wp * km(k,j,i)   * ( w(k,j,i)   - w(k-1,j,i) ) * ddzw(k)         &
                                                                                 * rho_air(k)
             ENDDO
!
!--          If surface fluxes should be used, overwrite the fluxes. Therefore, first determine
!--          start and end indices for the surface elements stored for grid point (j,i). If there
!--          are no wall grid points at given (j,i), loops won't be entered. Further, the surface
!--          fluxes are also multiplied with the relevant normal vector component, which includes
!--          the sign of the flux. For the w-component tendency, u'w' and v'w' is
!--          added at walls that are bounded along the east-/ west- and north/south-direction,
!--          respectively. No fluxes are added along the vertical because these wouldn't result
!--          from any shear stress but from momentum of the wall itself (which is 0).
             IF ( use_surface_fluxes )  THEN
!.
!--             Default-type surfaces
                surf_s = surf_def%start_index(j,i)
                surf_e = surf_def%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k         = surf_def%k(m)
                   flux_r(k) = MERGE( surf_def%wsus_wsvs(m), flux_r(k), surf_def%westward(m)  )
                   flux_l(k) = MERGE( surf_def%wsus_wsvs(m), flux_l(k), surf_def%eastward(m)  )
                   flux_n(k) = MERGE( surf_def%wsus_wsvs(m), flux_n(k), surf_def%southward(m) )
                   flux_s(k) = MERGE( surf_def%wsus_wsvs(m), flux_s(k), surf_def%northward(m) )
                ENDDO
!
!--             Natural-type surfaces
                surf_s = surf_lsm%start_index(j,i)
                surf_e = surf_lsm%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k         = surf_lsm%k(m)
                   flux_r(k) = MERGE( surf_lsm%wsus_wsvs(m), flux_r(k), surf_lsm%westward(m)  )
                   flux_l(k) = MERGE( surf_lsm%wsus_wsvs(m), flux_l(k), surf_lsm%eastward(m)  )
                   flux_n(k) = MERGE( surf_lsm%wsus_wsvs(m), flux_n(k), surf_lsm%southward(m) )
                   flux_s(k) = MERGE( surf_lsm%wsus_wsvs(m), flux_s(k), surf_lsm%northward(m) )
                ENDDO
!
!--             Urban-type surfaces
                surf_s = surf_usm%start_index(j,i)
                surf_e = surf_usm%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k         = surf_usm%k(m)
                   flux_r(k) = MERGE( surf_usm%wsus_wsvs(m), flux_r(k), surf_usm%westward(m)  )
                   flux_l(k) = MERGE( surf_usm%wsus_wsvs(m), flux_l(k), surf_usm%eastward(m)  )
                   flux_n(k) = MERGE( surf_usm%wsus_wsvs(m), flux_n(k), surf_usm%southward(m) )
                   flux_s(k) = MERGE( surf_usm%wsus_wsvs(m), flux_s(k), surf_usm%northward(m) )
                ENDDO
             ENDIF
!
!--          Compute tendency. Flag 3 is used to mask topography on w-grid.
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt-1
                flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
                tend(k,j,i) = tend(k,j,i) -                                                        &
                                        (   ( flux_r(k) - flux_l(k) ) * ddx                        &
                                          + ( flux_n(k) - flux_s(k) ) * ddy                        &
                                          + ( flux_t(k) - flux_d(k) ) * ddzu(k+1) * drho_air_zw(k) &
                                        ) * flag
             ENDDO

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_w


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_w_ij( i, j )

       IMPLICIT NONE


       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  kmxm          !< diffusion coefficient on leftward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmxp          !< diffusion coefficient on rightward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmym          !< diffusion coefficient on southward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  kmyp          !< diffusion coefficient on northward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  mask_east     !< flag to mask vertical wall east of the grid point
       REAL(wp) ::  mask_north    !< flag to mask vertical wall north of the grid point
       REAL(wp) ::  mask_south    !< flag to mask vertical wall south of the grid point
       REAL(wp) ::  mask_west     !< flag to mask vertical wall west of the grid point

       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_d   !< flux at k-1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_t   !< flux at k+1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_l   !< flux at i-1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_r   !< flux at i+1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_n   !< flux at j+1/2 (from w-grid point)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_s   !< flux at j-1/2 (from w-grid point)

!
!--    Compute diffusion fluxes at the grid-cell boundaries using deformation approach.
!--    Note, at wall-bounded grid points the respective fluxes will be overwritten later.
       DO  k = nzb+1, nzt-1
!
!--       Predetermine flags mask wall-bounded grid points.
          mask_east  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 3 ) )
          mask_west  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 3 ) )
          mask_south = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j-1,i), 3 ) )
          mask_north = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j+1,i), 3 ) )
!
!--       Interpolate eddy diffusivities onto staggered grid
          kmxp = 0.25_wp * ( km(k,j,i) + km(k,j,i+1) + km(k+1,j,i) + km(k+1,j,i+1) )
          kmxm = 0.25_wp * ( km(k,j,i) + km(k,j,i-1) + km(k+1,j,i) + km(k+1,j,i-1) )
          kmyp = 0.25_wp * ( km(k,j,i) + km(k+1,j,i) + km(k,j+1,i) + km(k+1,j+1,i) )
          kmym = 0.25_wp * ( km(k,j,i) + km(k+1,j,i) + km(k,j-1,i) + km(k+1,j-1,i) )

          flux_r(k) = - kmxp * ( ( w(k,j,i+1)   - w(k,j,i)   ) * ddx                               &
                               + ( u(k+1,j,i+1) - u(k,j,i+1) ) * ddzu(k+1)                         &
                               ) * mask_east
          flux_l(k) = - kmxm * ( ( w(k,j,i)     - w(k,j,i-1) ) * ddx                               &
                               + ( u(k+1,j,i)   - u(k,j,i)   ) * ddzu(k+1)                         &
                               ) * mask_west
          flux_n(k) = - kmyp * ( ( w(k,j+1,i)   - w(k,j,i)   ) * ddy                               &
                               + ( v(k+1,j+1,i) - v(k,j+1,i) ) * ddzu(k+1)                         &
                               ) * mask_north
          flux_s(k) = - kmym * ( ( w(k,j,i)     - w(k,j-1,i) ) * ddy                               &
                               + ( v(k+1,j,i)   - v(k,j,i)   ) * ddzu(k+1)                         &
                               ) * mask_south
          flux_t(k) = - 2.0_wp * km(k+1,j,i) * ( w(k+1,j,i) - w(k,j,i)   ) * ddzw(k+1) * rho_air(k+1)
          flux_d(k) = - 2.0_wp * km(k,j,i)   * ( w(k,j,i)   - w(k-1,j,i) ) * ddzw(k)   * rho_air(k)
       ENDDO
!
!--    If surface fluxes should be used, overwrite the fluxes. Therefore, first determine start and
!--    end indices for the surface elements stored for grid point (j,i). If there are no wall grid
!--    points at given (j,i), loops won't be entered. Further, the surface fluxes are also
!--    multiplied with the relevant normal vector component, which includes the sign of the flux.
!--    For the w-component tendency, u'w' and v'w' is added at walls that are bounded along the east-/
!--    west- and north/south-direction, respectively. No fluxes are added along the vertical direction
!--    because these wouldn't result from any shear stress but from momentum of the wall itself
!--    (which is 0).
       IF ( use_surface_fluxes )  THEN
!.
!--       Default-type surfaces
          surf_s = surf_def%start_index(j,i)
          surf_e = surf_def%end_index(j,i)
          DO  m = surf_s, surf_e
             k         = surf_def%k(m)
             flux_r(k) = MERGE( surf_def%wsus_wsvs(m), flux_r(k), surf_def%westward(m)  )
             flux_l(k) = MERGE( surf_def%wsus_wsvs(m), flux_l(k), surf_def%eastward(m)  )
             flux_n(k) = MERGE( surf_def%wsus_wsvs(m), flux_n(k), surf_def%southward(m) )
             flux_s(k) = MERGE( surf_def%wsus_wsvs(m), flux_s(k), surf_def%northward(m) )
          ENDDO
!
!--       Natural-type surfaces
          surf_s = surf_lsm%start_index(j,i)
          surf_e = surf_lsm%end_index(j,i)
          DO  m = surf_s, surf_e
             k         = surf_lsm%k(m)
             flux_r(k) = MERGE( surf_lsm%wsus_wsvs(m), flux_r(k), surf_lsm%westward(m)  )
             flux_l(k) = MERGE( surf_lsm%wsus_wsvs(m), flux_l(k), surf_lsm%eastward(m)  )
             flux_n(k) = MERGE( surf_lsm%wsus_wsvs(m), flux_n(k), surf_lsm%southward(m) )
             flux_s(k) = MERGE( surf_lsm%wsus_wsvs(m), flux_s(k), surf_lsm%northward(m) )
          ENDDO
!
!--       Urban-type surfaces
          surf_s = surf_usm%start_index(j,i)
          surf_e = surf_usm%end_index(j,i)
          DO  m = surf_s, surf_e
             k         = surf_usm%k(m)
             flux_r(k) = MERGE( surf_usm%wsus_wsvs(m), flux_r(k), surf_usm%westward(m)  )
             flux_l(k) = MERGE( surf_usm%wsus_wsvs(m), flux_l(k), surf_usm%eastward(m)  )
             flux_n(k) = MERGE( surf_usm%wsus_wsvs(m), flux_n(k), surf_usm%southward(m) )
             flux_s(k) = MERGE( surf_usm%wsus_wsvs(m), flux_s(k), surf_usm%northward(m) )
          ENDDO
       ENDIF

!
!--    Compute tendency. Flag 3 is used to mask topography on w-grid.
       DO  k = nzb+1, nzt-1
          flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
          tend(k,j,i) = tend(k,j,i) - (   ( flux_r(k) - flux_l(k) ) * ddx                          &
                                        + ( flux_n(k) - flux_s(k) ) * ddy                          &
                                        + ( flux_t(k) - flux_d(k) ) * ddzu(k+1) * drho_air_zw(k)   &
                                      ) * flag
       ENDDO

    END SUBROUTINE diffusion_w_ij

 END MODULE diffusion_w_mod
