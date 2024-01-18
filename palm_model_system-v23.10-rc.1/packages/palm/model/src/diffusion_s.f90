!> @file diffusion_s.f90
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
!> Diffusion term of scalar quantities (temperature and water content)
!--------------------------------------------------------------------------------------------------!
 MODULE diffusion_s_mod

    USE arrays_3d,                                                                                 &
        ONLY:  ddzu,                                                                               &
               ddzw,                                                                               &
               drho_air,                                                                           &
               dzw,                                                                                &
               kh,                                                                                 &
               tend,                                                                               &
               rho_air_zw

    USE control_parameters,                                                                        &
        ONLY: use_surface_fluxes,                                                                  &
              use_top_fluxes

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

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
               nzt,                                                                                &
               topo_flags

    USE kinds

    USE surface_mod,                                                                               &
        ONLY:  surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_top,                                                                           &
               surf_usm

    PRIVATE
    PUBLIC diffusion_s

    INTERFACE diffusion_s
       MODULE PROCEDURE diffusion_s
       MODULE PROCEDURE diffusion_s_ij
    END INTERFACE diffusion_s

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_s( s, s_flux_t, s_flux_def, s_flux_lsm, s_flux_usm )



       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  flux_tmp          !< temporary variable used to map surface flux onto the correct grid-cell flux

       REAL(wp), DIMENSION(1:surf_def%ns) ::  s_flux_def  !< fluxes defined at default-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm%ns) ::  s_flux_lsm  !< fluxes defined at natural-type surfaces
       REAL(wp), DIMENSION(1:surf_usm%ns) ::  s_flux_usm  !< fluxes defined at urban-type surfaces
       REAL(wp), DIMENSION(1:surf_top%ns) ::  s_flux_t    !< fluxes at model top

       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_d   !< flux defined at k-1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_t   !< flux defined at k+1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_l   !< flux defined at i-1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_r   !< flux defined at i+1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_n   !< flux defined at j+1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_s   !< flux defined at j-1/2 (from grid-cell center)

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s  !< treated scalar


       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j) &
       !$ACC PRIVATE(surf_e, surf_s, flag, flux_tmp) &
       !$ACC PRIVATE(flux_r(nzb+1:nzt), flux_l(nzb+1:nzt), flux_n(nzb+1:nzt), flux_s(nzb+1:nzt), flux_t(nzb+1:nzt), flux_d(nzb+1:nzt)) &
       !$ACC PRESENT(topo_flags, kh) &
       !$ACC PRESENT(s) &
       !$ACC PRESENT(ddzu, ddzw, drho_air, dzw, rho_air_zw) &
       !$ACC PRESENT(surf_def, surf_lsm, surf_top, surf_usm) &
       !$ACC PRESENT(s_flux_def, s_flux_t, s_flux_lsm, s_flux_usm) &
       !$ACC PRESENT(tend)
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Compute horizontal and vertical diffusive fluxes. Note, no masking of topography or
!--          wall-bounded grid cells is done here. This will done further below.
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt
                flux_r(k) = -0.5_wp * ( kh(k,j,i) + kh(k,j,i+1) )                                  &
                                    * ( s(k,j,i+1) - s(k,j,i)   ) * ddx
                flux_l(k) = -0.5_wp * ( kh(k,j,i) + kh(k,j,i-1) )                                  &
                                    * ( s(k,j,i)   - s(k,j,i-1) ) * ddx
                flux_n(k) = -0.5_wp * ( kh(k,j,i) + kh(k,j+1,i) )                                  &
                                    * ( s(k,j+1,i) - s(k,j,i)   ) * ddy
                flux_s(k) = -0.5_wp * ( kh(k,j,i) + kh(k,j-1,i) )                                  &
                                    * ( s(k,j,i)   - s(k,j-1,i) ) * ddy
                flux_t(k) = -0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )                                  &
                                    * ( s(k+1,j,i) - s(k,j,i)   )                                  &
                                    * ddzu(k+1) * rho_air_zw(k)
                flux_d(k) = -0.5_wp * ( kh(k,j,i) + kh(k-1,j,i) )                                  &
                                    * ( s(k,j,i)   - s(k-1,j,i) )                                  &
                                    * ddzu(k) * rho_air_zw(k-1)
             ENDDO
!
!--          If surface fluxes should be used, overwrite the fluxes. Therefore, first determine
!--          start and end indices for the surface elements stored for grid point (j,i). If there
!--          are no wall grid points at given (j,i), loops won't be entered. Further, the surface
!--          fluxes are also multiplied with the relevant normal vector component, which includes
!--          the sign of the flux.
             IF ( use_surface_fluxes )  THEN
!
!--             First, for default-type surfaces - all orientations
                surf_s = surf_def%start_index(j,i)
                surf_e = surf_def%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k         = surf_def%k(m)
                   flux_tmp  = s_flux_def(m) * surf_def%n_eff(m)
                   flux_r(k) = MERGE( flux_tmp, flux_r(k), surf_def%westward(m)  )
                   flux_l(k) = MERGE( flux_tmp, flux_l(k), surf_def%eastward(m)  )
                   flux_n(k) = MERGE( flux_tmp, flux_n(k), surf_def%southward(m) )
                   flux_s(k) = MERGE( flux_tmp, flux_s(k), surf_def%northward(m) )
                   flux_t(k) = MERGE( flux_tmp, flux_t(k), surf_def%downward(m)  )
                   flux_d(k) = MERGE( flux_tmp, flux_d(k), surf_def%upward(m)    )
                ENDDO
!
!--             Now, natural-type surfaces - all orientations
                surf_s = surf_lsm%start_index(j,i)
                surf_e = surf_lsm%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k         = surf_lsm%k(m)
                   flux_tmp  = s_flux_lsm(m) * surf_lsm%n_eff(m)
                   flux_r(k) = MERGE( flux_tmp, flux_r(k), surf_lsm%westward(m)  )
                   flux_l(k) = MERGE( flux_tmp, flux_l(k), surf_lsm%eastward(m)  )
                   flux_n(k) = MERGE( flux_tmp, flux_n(k), surf_lsm%southward(m) )
                   flux_s(k) = MERGE( flux_tmp, flux_s(k), surf_lsm%northward(m) )
                   flux_t(k) = MERGE( flux_tmp, flux_t(k), surf_lsm%downward(m)  )
                   flux_d(k) = MERGE( flux_tmp, flux_d(k), surf_lsm%upward(m)    )
                ENDDO
!
!--             Now, for urban-type surfaces - all orientations
                surf_s = surf_usm%start_index(j,i)
                surf_e = surf_usm%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k         = surf_usm%k(m)
                   flux_tmp  = s_flux_usm(m) * surf_usm%n_eff(m)
                   flux_r(k) = MERGE( flux_tmp, flux_r(k), surf_usm%westward(m)  )
                   flux_l(k) = MERGE( flux_tmp, flux_l(k), surf_usm%eastward(m)  )
                   flux_n(k) = MERGE( flux_tmp, flux_n(k), surf_usm%southward(m) )
                   flux_s(k) = MERGE( flux_tmp, flux_s(k), surf_usm%northward(m) )
                   flux_t(k) = MERGE( flux_tmp, flux_t(k), surf_usm%downward(m)  )
                   flux_d(k) = MERGE( flux_tmp, flux_d(k), surf_usm%upward(m)    )
                ENDDO
             ENDIF
!
!--          Vertical diffusion along z-direction at the top boundary
             IF ( use_top_fluxes )  THEN
                surf_s = surf_top%start_index(j,i)
                surf_e = surf_top%end_index(j,i)
                !$ACC LOOP PRIVATE(k, m)
                DO  m = surf_s, surf_e
                   k         = surf_top%k(m)
                   flux_t(k) = s_flux_t(m)
                ENDDO
             ENDIF
!
!--          Compute tendency. Flag is used to mask topography.
             !$ACC LOOP PRIVATE(k)
             DO  k = nzb+1, nzt
                flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                tend(k,j,i) = tend(k,j,i) - (   ( flux_r(k) - flux_l(k) ) * ddx                    &
                                              + ( flux_n(k) - flux_s(k) ) * ddy                    &
                                              + ( flux_t(k) - flux_d(k) ) * ddzw(k) * drho_air(k)  &
                                            ) * flag
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE diffusion_s

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_s_ij( i, j, s, s_flux_t, s_flux_def, s_flux_lsm, s_flux_usm )

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !< running index x direction
       INTEGER(iwp) ::  j      !< running index y direction
       INTEGER(iwp) ::  k      !< running index z direction
       INTEGER(iwp) ::  m      !< running index surface elements
       INTEGER(iwp) ::  surf_e !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag       !< flag to mask topography grid points
       REAL(wp) ::  flux_tmp   !< temporary variable used to map surface flux onto the correct grid-cell flux

       REAL(wp), DIMENSION(1:surf_def%ns) ::  s_flux_def  !< fluxes defined at default-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm%ns) ::  s_flux_lsm  !< fluxes defined at natural-type surfaces
       REAL(wp), DIMENSION(1:surf_usm%ns) ::  s_flux_usm  !< fluxes defined at urban-type surfaces
       REAL(wp), DIMENSION(1:surf_top%ns) ::  s_flux_t    !< fluxes at model top

       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_d   !< flux defined at k-1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_t   !< flux defined at k+1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_l   !< flux defined at i-1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_r   !< flux defined at i+1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_n   !< flux defined at j+1/2 (from grid-cell center)
       REAL(wp), DIMENSION(nzb+1:nzt) ::  flux_s   !< flux defined at j-1/2 (from grid-cell center)

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s  !< treated scalar
!
!--    Compute horizontal and vertical diffusive fluxes. Note, no masking of topography or
!--    wall-bounded grid cells is done here. This will done further below.
       DO  k = nzb+1, nzt
          flux_r(k) = - 0.5_wp * ( kh(k,j,i) + kh(k,j,i+1) ) * ( s(k,j,i+1) - s(k,j,i)   ) * ddx
          flux_l(k) = - 0.5_wp * ( kh(k,j,i) + kh(k,j,i-1) ) * ( s(k,j,i)   - s(k,j,i-1) ) * ddx
          flux_n(k) = - 0.5_wp * ( kh(k,j,i) + kh(k,j+1,i) ) * ( s(k,j+1,i) - s(k,j,i)   ) * ddy
          flux_s(k) = - 0.5_wp * ( kh(k,j,i) + kh(k,j-1,i) ) * ( s(k,j,i)   - s(k,j-1,i) ) * ddy
          flux_t(k) = - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) ) * ( s(k+1,j,i) - s(k,j,i)   )         &
                               * ddzu(k+1) * rho_air_zw(k)
          flux_d(k) = - 0.5_wp * ( kh(k,j,i) + kh(k-1,j,i) ) * ( s(k,j,i)   - s(k-1,j,i) )         &
                               * ddzu(k) * rho_air_zw(k-1)
       ENDDO
!
!--    If surface fluxes should be used, overwrite the fluxes. Therefore, first determine start and
!--    end indices for the surface elements stored for grid point (j,i). If there are no wall grid
!--    points at given (j,i), loops won't be entered. Further, the surface fluxes are also
!--    multiplied with the relevant normal vector component, which includes the sign of the flux.
       IF ( use_surface_fluxes )  THEN
!
!--       First, for default-type surfaces - all orientations
          surf_s = surf_def%start_index(j,i)
          surf_e = surf_def%end_index(j,i)
          DO  m = surf_s, surf_e
             k         = surf_def%k(m)
             flux_tmp  = s_flux_def(m) * surf_def%n_eff(m)
             flux_r(k) = MERGE( flux_tmp, flux_r(k), surf_def%westward(m)  )
             flux_l(k) = MERGE( flux_tmp, flux_l(k), surf_def%eastward(m)  )
             flux_n(k) = MERGE( flux_tmp, flux_n(k), surf_def%southward(m) )
             flux_s(k) = MERGE( flux_tmp, flux_s(k), surf_def%northward(m) )
             flux_t(k) = MERGE( flux_tmp, flux_t(k), surf_def%downward(m)  )
             flux_d(k) = MERGE( flux_tmp, flux_d(k), surf_def%upward(m)    )
          ENDDO
!
!--       Now, natural-type surfaces - all orientations
          surf_s = surf_lsm%start_index(j,i)
          surf_e = surf_lsm%end_index(j,i)
          DO  m = surf_s, surf_e
             k         = surf_lsm%k(m)
             flux_tmp  = s_flux_lsm(m) * surf_lsm%n_eff(m)
             flux_r(k) = MERGE( flux_tmp, flux_r(k), surf_lsm%westward(m)  )
             flux_l(k) = MERGE( flux_tmp, flux_l(k), surf_lsm%eastward(m)  )
             flux_n(k) = MERGE( flux_tmp, flux_n(k), surf_lsm%southward(m) )
             flux_s(k) = MERGE( flux_tmp, flux_s(k), surf_lsm%northward(m) )
             flux_t(k) = MERGE( flux_tmp, flux_t(k), surf_lsm%downward(m)  )
             flux_d(k) = MERGE( flux_tmp, flux_d(k), surf_lsm%upward(m)    )
          ENDDO
!
!--       Now, for urban-type surfaces - all orientations
          surf_s = surf_usm%start_index(j,i)
          surf_e = surf_usm%end_index(j,i)
          DO  m = surf_s, surf_e
             k         = surf_usm%k(m)
             flux_tmp  = s_flux_usm(m) * surf_usm%n_eff(m)
             flux_r(k) = MERGE( flux_tmp, flux_r(k), surf_usm%westward(m)  )
             flux_l(k) = MERGE( flux_tmp, flux_l(k), surf_usm%eastward(m)  )
             flux_n(k) = MERGE( flux_tmp, flux_n(k), surf_usm%southward(m) )
             flux_s(k) = MERGE( flux_tmp, flux_s(k), surf_usm%northward(m) )
             flux_t(k) = MERGE( flux_tmp, flux_t(k), surf_usm%downward(m)  )
             flux_d(k) = MERGE( flux_tmp, flux_d(k), surf_usm%upward(m)    )
          ENDDO
       ENDIF
!
!--    Vertical diffusion along z-direction at the top boundary
       IF ( use_top_fluxes )  THEN
          surf_s = surf_top%start_index(j,i)
          surf_e = surf_top%end_index(j,i)
          DO  m = surf_s, surf_e
             k         = surf_top%k(m)
             flux_t(k) = s_flux_t(m)
          ENDDO
       ENDIF
!
!--    Compute tendency. Flag is used to mask topography.
       DO  k = nzb+1, nzt
          flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
          tend(k,j,i) = tend(k,j,i) - (   ( flux_r(k) - flux_l(k) ) * ddx                          &
                                        + ( flux_n(k) - flux_s(k) ) * ddy                          &
                                        + ( flux_t(k) - flux_d(k) ) * ddzw(k) * drho_air(k)        &
                                      ) * flag
       ENDDO


    END SUBROUTINE diffusion_s_ij

 END MODULE diffusion_s_mod
