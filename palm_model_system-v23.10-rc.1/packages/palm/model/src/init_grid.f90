!> @file init_grid.f90
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
! -------------------------------------------------------------------------------------------------!
!> Creating grid depending constants
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_grid

    USE arrays_3d,                                                                                 &
        ONLY:  dd2zu,                                                                              &
               ddzu,                                                                               &
               ddzu_pres,                                                                          &
               ddzw,                                                                               &
               dzu,                                                                                &
               dzw,                                                                                &
               x,                                                                                  &
               xu,                                                                                 &
               y,                                                                                  &
               yv,                                                                                 &
               zu,                                                                                 &
               zw

    USE chem_modules,                                                                              &
        ONLY:  bc_cs_b,                                                                            &
               bc_cs_t

    USE control_parameters,                                                                        &
        ONLY:  bc_p_b,                                                                             &
               bc_p_t,                                                                             &
               bc_pt_b,                                                                            &
               bc_pt_t,                                                                            &
               bc_q_b,                                                                             &
               bc_q_t,                                                                             &
               bc_s_b,                                                                             &
               bc_s_t,                                                                             &
               bc_uv_b,                                                                            &
               bc_uv_t,                                                                            &
               child_domain,                                                                       &
               constant_flux_layer,                                                                &
               dz,                                                                                 &
               dz_max,                                                                             &
               dz_stretch_factor,                                                                  &
               dz_stretch_factor_array,                                                            &
               dz_stretch_level,                                                                   &
               dz_stretch_level_end,                                                               &
               dz_stretch_level_end_index,                                                         &
               dz_stretch_level_start_index,                                                       &
               dz_stretch_level_start,                                                             &
               message_string,                                                                     &
               nested_run,                                                                         &
               number_dz,                                                                          &
               number_stretch_level_end,                                                           &
               number_stretch_level_start,                                                         &
               ocean_mode,                                                                         &
               psolver,                                                                            &
               symmetry_flag,                                                                      &
               topography,                                                                         &
               use_surface_fluxes

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddx2,                                                                               &
               ddy,                                                                                &
               ddy2,                                                                               &
               dx,                                                                                 &
               dx2,                                                                                &
               dy,                                                                                 &
               dy2

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nzb,                                                                                &
               nzb_diff,                                                                           &
               nzt

    USE kinds

    USE ocean_mod,                                                                                 &
        ONLY:  bc_sa_b

    USE pegrid

    USE pmc_interface,                                                                             &
        ONLY:  nest_shift_z

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !< index variable along x
    INTEGER(iwp) ::  j             !< index variable along y
    INTEGER(iwp) ::  k             !< index variable along z
    INTEGER(iwp) ::  n             !< loop variable for stretching

    REAL(wp) ::  dz_level_end  !< distance between calculated height level for u/v-grid and user-specified end level for stretching
    REAL(wp) ::  dz_stretched  !< stretched vertical grid spacing

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  min_dz_stretch_level_end !< Array that contains all minimum heights where the stretching
                                                                     !< can end

!
!-- Calculation of horizontal array bounds including ghost layers
    nxlg = nxl - nbgp
    nxrg = nxr + nbgp
    nysg = nys - nbgp
    nyng = nyn + nbgp

!
!-- Allocate grid arrays
    ALLOCATE( x(-nbgp:nx+nbgp) )
    ALLOCATE( xu(-nbgp:nx+nbgp) )

    DO i = -nbgp, nx+nbgp
       xu(i) = i * dx
       x(i)  = i * dx + 0.5_wp * dx
    ENDDO

    ALLOCATE( y(-nbgp:ny+nbgp) )
    ALLOCATE( yv(-nbgp:ny+nbgp) )

    DO j = -nbgp, ny+nbgp
       yv(j) = j * dy
       y(j)  = j * dy + 0.5_wp * dy
    ENDDO

    ALLOCATE( ddzu(1:nzt+1) )
    ALLOCATE( ddzw(1:nzt+1) )
    ALLOCATE( dd2zu(1:nzt) )
    ALLOCATE( dzu(1:nzt+1) )
    ALLOCATE( dzw(1:nzt+1) )
    ALLOCATE( zu(nzb:nzt+1) )
    ALLOCATE( zw(nzb:nzt+1) )

!
!-- For constructing an appropriate grid, the vertical grid spacing dz has to be specified with a
!-- non-negative value in the parameter file.
    IF ( dz(1) == -1.0_wp )  THEN
       message_string = 'missing dz'
       CALL message( 'init_grid', 'PAC0210', 1, 2, 0, 6, 0 )
    ELSEIF ( dz(1) <= 0.0_wp )  THEN
       WRITE( message_string, * ) 'dz=',dz(1),' <= 0.0'
       CALL message( 'init_grid', 'PAC0211', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize dz_stretch_level_start with the value of dz_stretch_level if it was set by the user.
    IF ( dz_stretch_level /= -9999999.9_wp ) THEN
       dz_stretch_level_start(1) = dz_stretch_level
    ENDIF

!
!-- Determine number of dz values and stretching levels specified by the user to allow right
!-- controlling of the stretching mechanism and to perform error checks. The additional requirement
!-- that dz /= dz_max for counting number of user-specified dz values is necessary. Otherwise
!-- restarts would abort if the old stretching mechanism with dz_stretch_level is used (Attention:
!-- The user is not allowed to specify a dz value equal to the default of dz_max = 999.0).
    number_dz = COUNT( dz /= -1.0_wp  .AND.  dz /= dz_max)
    number_stretch_level_start = COUNT( dz_stretch_level_start /= -9999999.9_wp )
    number_stretch_level_end = COUNT( dz_stretch_level_end /= 9999999.9_wp )

!
!-- The number of specified end levels +1 has to be the same as the number
!-- of specified dz values
    IF ( number_dz /= number_stretch_level_end + 1 ) THEN
       WRITE( message_string, * )  'number of values for dz = ', number_dz,                        &
                                   'has to be the same as& ', 'the number of values for ',         &
                                   'dz_stretch_level_end + 1 = ', number_stretch_level_end+1
       CALL message( 'init_grid', 'PAC0212', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- The number of specified start levels has to be the same or one less than the number of specified
!-- dz values
    IF ( number_dz /= number_stretch_level_start + 1  .AND.                                        &
         number_dz /= number_stretch_level_start )  THEN
       WRITE( message_string, * )  'number of values for dz = ', number_dz,                        &
                                   'has to be the same as or one ',                                &
                                   'more than& the number of values for ',                         &
                                   'dz_stretch_level_start = ', number_stretch_level_start
       CALL message( 'init_grid', 'PAC0213', 1, 2, 0, 6, 0 )
    ENDIF

!-- The number of specified start levels has to be the same or one more than the number of specified
!-- end levels
    IF ( number_stretch_level_start /= number_stretch_level_end + 1  .AND.                         &
         number_stretch_level_start /= number_stretch_level_end ) THEN
       WRITE( message_string, * )  'The number of values for ',                                    &
                                   'dz_stretch_level_start = ', dz_stretch_level_start,            &
                                   'has to be the ', 'same or one more than& the number of ',      &
                                   'values for dz_stretch_level_end = ', number_stretch_level_end
       CALL message( 'init_grid', 'PAC0214', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize dz for the free atmosphere with the value of dz_max
    IF ( dz(number_stretch_level_start+1) == -1.0_wp  .AND.  number_stretch_level_start /= 0 ) THEN
       dz(number_stretch_level_start+1) = dz_max
    ENDIF

!
!-- Initialize the stretching factor if (infinitely) stretching in the free atmosphere is desired
!-- (dz_stretch_level_end was not specified for the free atmosphere)
    IF ( number_stretch_level_start == number_stretch_level_end + 1 )  THEN
       dz_stretch_factor_array(number_stretch_level_start) = dz_stretch_factor
    ENDIF

!
!-- Allocation of arrays for stretching
    ALLOCATE( min_dz_stretch_level_end(number_stretch_level_start) )

!
!-- Define the vertical grid levels. Start with atmosphere branch
    IF ( .NOT. ocean_mode )  THEN

!
!--    The stretching region has to be large enough to allow for a smooth transition between two
!--    different grid spacings. The number 4 is an empirical value.
       DO n = 1, number_stretch_level_start
          min_dz_stretch_level_end(n) = dz_stretch_level_start(n) + 4 * MAX( dz(n),dz(n+1) )
       ENDDO

       IF ( ANY( min_dz_stretch_level_end(1:number_stretch_level_start) >                          &
                 dz_stretch_level_end(1:number_stretch_level_start) ) )  THEN
          message_string= 'each dz_stretch_level_end has to be larger ' //                         &
                          'than its corresponding value for &' //                                  &
                          'dz_stretch_level_start + 4*MAX(dz(n),dz(n+1)) '//                       &
                          'to allow for smooth grid stretching'
          CALL message( 'init_grid', 'PAC0215', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Stretching must not be applied within the surface layer (first two grid points). For the
!--    default case dz_stretch_level_start is negative. Therefore the absolut value is checked here.
       IF ( ANY( ABS( dz_stretch_level_start ) <= dz(1) * 1.5_wp ) ) THEN
          WRITE( message_string, * )  'at least one value of ABS( dz_stretch_level_start ) ',      &
                                      'is smaller or equal ', dz(1) * 1.5, '(dz(1) * 1.5)'
          CALL message( 'init_grid', 'PAC0216', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    The stretching has to start and end on a grid level. Therefore user-specified values are
!--    mapped to the next lowest level. The calculation of the first level is realized differently
!--    just because of historical reasons (the advanced/new stretching mechanism was realized in a
!--    way that results don't change if the old parameters dz_stretch_level, dz_stretch_factor and
!--    dz_max are used).
       IF ( number_stretch_level_start /= 0 )  THEN
          dz_stretch_level_start(1) = INT( (dz_stretch_level_start(1) - dz(1)/2.0) / dz(1) )       &
                                      * dz(1) + dz(1)/2.0
       ENDIF

       IF ( number_stretch_level_start > 1 ) THEN
          DO n = 2, number_stretch_level_start
             dz_stretch_level_start(n) = INT( dz_stretch_level_start(n) / dz(n) ) * dz(n)
          ENDDO
       ENDIF

       IF ( number_stretch_level_end /= 0 ) THEN
          DO n = 1, number_stretch_level_end
             dz_stretch_level_end(n) = INT( dz_stretch_level_end(n) / dz(n+1) ) * dz(n+1)
          ENDDO
       ENDIF

!
!--    Determine stretching factor if necessary
       IF ( number_stretch_level_end >= 1 ) THEN
          CALL calculate_stretching_factor( number_stretch_level_end )
       ENDIF

!
!--    Grid for atmosphere with surface at z=0. This corresponds to k=0 on the w-grid AND
!--    the u-grid. Vertical shift of elevated childs is considered via nest_shift_z, which is
!--    zero by default.
!--    First compute the u- and v-levels.
!--    The second u-level (k=1) corresponds to the top of the surface layer.
       IF ( nest_shift_z == 0.0_wp )  THEN
!
!--       Bottom boundary is a real surface and zu(0) is defined on that one. In such a case,
!--       zu(0) is used only for output purposes and does not affect the calculations.
          zu(0) = nest_shift_z
          zu(1) = zu(0) + dz(1) * 0.5_wp
       ELSE
!
!--       For elevated nests the bottom grid point positions are required for correct calculation
!--       of vertical advection and diffusion of u, v, and scalars at k=1. There is no surface, so
!--       the points have a distance of one grid spacing.
          zu(0) = nest_shift_z - dz(1) * 0.5_wp
          zu(1) = nest_shift_z + dz(1) * 0.5_wp
       ENDIF

!
!--    Determine u and v height levels considering the possibility of grid stretching in several
!--    heights.
       n = 1
       dz_stretch_level_start_index = nzt+1
       dz_stretch_level_end_index = nzt+1
       dz_stretched = dz(1)

!--    The default value of dz_stretch_level_start is negative, thus the first condition is true
!--    even if no stretching shall be applied. Hence, the second condition is also necessary.
       DO  k = 2, nzt+1-symmetry_flag
          IF ( dz_stretch_level_start(n) <= zu(k-1)  .AND.                                         &
               dz_stretch_level_start(n) /= -9999999.9_wp )  THEN
             dz_stretched = dz_stretched * dz_stretch_factor_array(n)

             IF ( dz(n) > dz(n+1) )  THEN
                dz_stretched = MAX( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (higher) dz
             ELSE
                dz_stretched = MIN( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (lower) dz
             ENDIF

             IF ( dz_stretch_level_start_index(n) == nzt+1 )  dz_stretch_level_start_index(n) = k-1

          ENDIF

          zu(k) = zu(k-1) + dz_stretched

!
!--       Make sure that the stretching ends exactly at dz_stretch_level_end
          dz_level_end = ABS( zu(k) - dz_stretch_level_end(n) )

          IF ( dz_level_end  <  dz(n+1) / 3.0 )  THEN
             zu(k) = dz_stretch_level_end(n)
             dz_stretched = dz(n+1)
             dz_stretch_level_end_index(n) = k
             n = n + 1
          ENDIF
       ENDDO

!
!--    If a closed channel flow is simulated, make sure that grid structure is the same for both
!--    bottom and top boundary. (Hint: Using a different dz at the bottom and at the top makes no
!--    sense due to symmetric boundaries where dz should be equal. Therefore, different dz at the
!--    bottom and top causes an abort (see check_parameters).)
       IF ( topography == 'closed_channel' )  THEN
          zu(nzt+1) = zu(nzt) + dz(1) * 0.5_wp
       ENDIF

!
!--    Compute the w-levels. They are always staggered half-way between the corresponding u-levels.
!--    In case of dirichlet bc for u and v at the ground the first u- and w-level (k=0) are defined
!--    at same height (z=0). Vertical shift of elevated childs is considered via nest_shift_z,
!--    which is zero by default.
!--    Per default, the top w-level is extrapolated linearly. In case of a closed channel flow,
!--    zu(nzt+1) and zw(nzt) must be set explicitely.
!--    (Hint: Using a different dz at the bottom and at the top makes no sense due to symmetric
!--    boundaries where dz should be equal. Therefore, different dz at the bottom and top causes an
!--    abort (see check_parameters).)
       zw(0) = nest_shift_z
       DO  k = 1, nzt-symmetry_flag
          zw(k) = ( zu(k) + zu(k+1) ) * 0.5_wp
       ENDDO
       IF ( topography == 'closed_channel' )  THEN
          zw(nzt)   = zw(nzt-1) + dz(1)
          zw(nzt+1) = zw(nzt) + dz(1)
       ELSE
          zw(nzt+1) = zw(nzt) + 2.0_wp * ( zu(nzt+1) - zw(nzt) )
       ENDIF

    ELSE !ocean branch

!
!--    The stretching region has to be large enough to allow for a smooth transition between two
!--    different grid spacings. The number 4 is an empirical value
       DO n = 1, number_stretch_level_start
          min_dz_stretch_level_end(n) = dz_stretch_level_start(n) - 4 * MAX( dz(n),dz(n+1) )
       ENDDO

       IF ( ANY( min_dz_stretch_level_end (1:number_stretch_level_start) <                         &
                 dz_stretch_level_end(1:number_stretch_level_start) ) )  THEN
             message_string= 'Each dz_stretch_level_end has to be less ' //                        &
                             'than its corresponding value for &' //                               &
                             'dz_stretch_level_start - 4*MAX(dz(n),dz(n+1)) '//                    &
                             'to allow for smooth grid stretching'
             CALL message( 'init_grid', 'PAC0215', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Stretching must not be applied close to the surface (last two grid points). For the default
!--    case dz_stretch_level_start is negative.
       IF ( ANY( dz_stretch_level_start >= -dz(1) * 1.5_wp ) )  THEN
          WRITE( message_string, * )  'at least one value of dz_stretch_level_start ',             &
                                      'is greater or equal ', -dz(1) * 1.5, '(-dz(1) * 1.5)'
             CALL message( 'init_grid', 'PAC0216', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    The stretching has to start and end on a grid level. Therefore user-specified values are
!--    mapped to the next highest level. The calculation of the first level is realized differently
!--    just because of historical reasons (the advanced/new stretching mechanism was realized in a
!--    way that results don't change if the old parameters dz_stretch_level, dz_stretch_factor and
!--    dz_max are used)
       IF ( number_stretch_level_start /= 0 )  THEN
          dz_stretch_level_start(1) = INT( (dz_stretch_level_start(1) + dz(1)/2.0) / dz(1) )       &
                                      * dz(1) - dz(1)/2.0
       ENDIF

       IF ( number_stretch_level_start > 1 )  THEN
          DO n = 2, number_stretch_level_start
             dz_stretch_level_start(n) = INT( dz_stretch_level_start(n) / dz(n) ) * dz(n)
          ENDDO
       ENDIF

       IF ( number_stretch_level_end /= 0 )  THEN
          DO n = 1, number_stretch_level_end
             dz_stretch_level_end(n) = INT( dz_stretch_level_end(n) / dz(n+1) ) * dz(n+1)
          ENDDO
       ENDIF

!
!--    Determine stretching factor if necessary
       IF ( number_stretch_level_end >= 1 ) THEN
          CALL calculate_stretching_factor( number_stretch_level_end )
       ENDIF

!
!--    Grid for ocean with free water surface is at k=nzt (w-grid).
!--    In case of neumann bc at the ground the first first u-level (k=0) lies below the first
!--    w-level (k=0). In case of dirichlet bc the first u- and w-level are defined at same height,
!--    but staggered from the second level. Vertical shift of elevated childs is considered via
!--    nest_shift_z, which is zero by default.
!--    The second u-level (k=nzt) corresponds to the top of the surface layer.
!--    z values are negative starting from z=0 (surface)
       zu(nzt+1) =   dz(1) * 0.5_wp + nest_shift_z
       zu(nzt)   = - dz(1) * 0.5_wp + nest_shift_z

!
!--    Determine u and v height levels considering the possibility of grid stretching in several
!--    heights.
       n = 1
       dz_stretch_level_start_index = 0
       dz_stretch_level_end_index = 0
       dz_stretched = dz(1)

       DO  k = nzt-1, 0, -1

          IF ( dz_stretch_level_start(n) >= zu(k+1) )  THEN
             dz_stretched = dz_stretched * dz_stretch_factor_array(n)

             IF ( dz(n) > dz(n+1) )  THEN
                dz_stretched = MAX( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (higher) dz
             ELSE
                dz_stretched = MIN( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (lower) dz
             ENDIF

             IF ( dz_stretch_level_start_index(n) == 0 )  dz_stretch_level_start_index(n) = k+1

          ENDIF

          zu(k) = zu(k+1) - dz_stretched

!
!--       Make sure that the stretching ends exactly at dz_stretch_level_end
          dz_level_end = ABS( zu(k) - dz_stretch_level_end(n) )

          IF ( dz_level_end  < dz(n+1)/3.0 )  THEN
             zu(k) = dz_stretch_level_end(n)
             dz_stretched = dz(n+1)
             dz_stretch_level_end_index(n) = k
             n = n + 1
          ENDIF
       ENDDO

!
!--    Compute the w-levels. They are always staggered half-way between the corresponding u-levels,
!--    except for u and v at the ground. In this case the first u- and
!--    w-level are defined at same height. The top w-level (nzt+1) is not used but set for
!--    consistency, since w and all scalar variables are defined up to nzt+1.
!--    Vertical shift of elevated childs is considered via nest_shift_z, which is zero by default.
       zw(nzt+1) = dz(1)  + nest_shift_z
       zw(nzt)   = 0.0_wp + nest_shift_z
       DO  k = 0, nzt
          zw(k) = ( zu(k) + zu(k+1) ) * 0.5_wp
       ENDDO

!
!--    Like for the atmosphere, the first u-grid and w-grid-level are defined at same height.
       zu(0) = zw(0)

    ENDIF !End of defining the vertical grid levels

!
!-- Compute grid lengths.
    DO  k = 1, nzt+1
       dzu(k)  = zu(k) - zu(k-1)
       ddzu(k) = 1.0_wp / dzu(k)
       dzw(k)  = zw(k) - zw(k-1)
       ddzw(k) = 1.0_wp / dzw(k)
    ENDDO

    DO  k = 1, nzt
       dd2zu(k) = 1.0_wp / ( dzu(k) + dzu(k+1) )
    ENDDO

!
!-- The FFT- SOR-pressure solvers assume grid spacings of a staggered grid everywhere. For the
!-- actual grid, the grid spacing at the lowest level is only dz/2, but should be dz. Therefore, an
!-- additional array containing with appropriate grid information is created for these solvers.
    IF ( psolver(1:9) /= 'multigrid' )  THEN
       ALLOCATE( ddzu_pres(1:nzt+1) )
       ddzu_pres = ddzu
       ddzu_pres(1) = ddzu_pres(2)  ! change for lowest level
    ENDIF

!
!-- Compute the reciprocal values of the horizontal grid lengths.
    ddx = 1.0_wp / dx
    ddy = 1.0_wp / dy
    dx2 = dx * dx
    dy2 = dy * dy
    ddx2 = 1.0_wp / dx2
    ddy2 = 1.0_wp / dy2

!
!-- Define vertical gridpoint from (or to) which on the usual finite difference form (which does not
!-- use surface fluxes) is applied. Note, this is only used in model_1d_mod.
    IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
       nzb_diff = nzb + 2
    ELSE
       nzb_diff = nzb + 1
    ENDIF

!
!-- Set bottom and top boundary conditions for childs depending on its vertical position.
    IF ( nested_run  .AND.  child_domain )  THEN

       IF ( .NOT. ocean_mode )  THEN
!
!--       Bottom condition for elevated (shifted) childs. Then the bottom boundary is "open".
          IF ( nest_shift_z /= 0.0_wp )  THEN
             bc_uv_b = 'nested'
             bc_pt_b = 'nested'
             bc_q_b  = 'nested'
             bc_s_b  = 'nested'
             bc_cs_b = 'nested'
             bc_p_b  = 'neumann'
          ENDIF
!
!--       Child top is always "open".
!> TODO:  Should be adjusted for the case if the top of the child matches the top of the root.
          bc_uv_t  = 'nested'
          bc_pt_t  = 'nested'
          bc_q_t   = 'nested'
          bc_s_t   = 'nested'
          bc_cs_t  = 'nested'
          bc_p_t   = 'neumann'

       ELSE
!
!--       The child top must match the sea surface, so the boundary conditions keep unchanged there.
!--       Child bottom is always "open".
!> TODO:  Should be adjusted for the case if the bottom of the child matches the bottom of the root.
          bc_uv_b = 'nested'
          bc_pt_b = 'nested'
          bc_q_b  = 'nested'
          bc_s_b  = 'nested'
          bc_sa_b = 'nested'
          bc_cs_b = 'nested'
          bc_p_b  = 'neumann'

       ENDIF

    ENDIF

 END SUBROUTINE init_grid


! Description:
! -------------------------------------------------------------------------------------------------!
!> Calculation of the stretching factor through an iterative method. Ideas were taken from the paper
!> "Regional stretched grid generation and its application to the NCAR RegCM (1999)". Normally, no
!> analytic solution exists because the system of equations has two variables (r,l) but four
!> requirements  (l=integer, r=[0,88;1,2], Eq(6), Eq(5) starting from index j=1) which results into
!> an overdetermined system.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calculate_stretching_factor( number_end )

    USE control_parameters,                                                                        &
        ONLY:  dz,                                                                                 &
               dz_stretch_factor_array,                                                            &
               dz_stretch_level_end,                                                               &
               dz_stretch_level_start,                                                             &
               message_string

    USE kinds

    IMPLICIT NONE

    REAL(wp), PARAMETER ::  stretch_factor_interval = 1.0E-06_wp  !< interval for sampling possible stretching factors
    REAL(wp), PARAMETER ::  stretch_factor_lower_limit = 0.88_wp  !< lowest possible stretching factor
    REAL(wp), PARAMETER ::  stretch_factor_upper_limit = 1.12_wp  !< highest possible stretching factor

    INTEGER(iwp) ::  iterations  !< number of iterations until stretch_factor_lower/upper_limit is reached
    INTEGER(iwp) ::  l_rounded   !< after l_rounded grid levels dz(n) is strechted to dz(n+1) with stretch_factor_2
    INTEGER(iwp) ::  n           !< loop variable for stretching

    INTEGER(iwp), INTENT(IN) ::  number_end !< number of user-specified end levels for stretching

    REAL(wp) ::  delta_l               !< absolute difference between l and l_rounded
    REAL(wp) ::  delta_stretch_factor  !< absolute difference between stretch_factor_1 and stretch_factor_2
    REAL(wp) ::  delta_total_new       !< sum of delta_l and delta_stretch_factor for the next iteration (should be as small as
                                       !< possible)
    REAL(wp) ::  delta_total_old       !< sum of delta_l and delta_stretch_factor for the last iteration
    REAL(wp) ::  distance              !< distance between dz_stretch_level_start and dz_stretch_level_end (stretching region)
    REAL(wp) ::  l                     !< value that fulfil Eq. (5) in the paper mentioned above together with stretch_factor_1
                                       !< exactly
    REAL(wp) ::  numerator             !< numerator of the quotient
    REAL(wp) ::  stretch_factor_1      !< stretching factor that fulfil Eq. (5) togehter with l exactly
    REAL(wp) ::  stretch_factor_2      !< stretching factor that fulfil Eq. (6) togehter with l_rounded exactly

    REAL(wp) ::  dz_stretch_factor_array_2(9) = 1.08_wp  !< Array that contains all stretch_factor_2 that belongs to
                                                         !< stretch_factor_1


    l = 0
    DO  n = 1, number_end

       iterations = 1
       stretch_factor_1 = 1.0_wp
       stretch_factor_2 = 1.0_wp
       delta_total_old = 1.0_wp

!
!--    First branch for stretching from rough to fine
       IF ( dz(n) > dz(n+1) )  THEN
          DO WHILE ( stretch_factor_1 >= stretch_factor_lower_limit )

             stretch_factor_1 = 1.0_wp - iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) - dz_stretch_level_start(n) )
             numerator = distance * stretch_factor_1 / dz(n) + stretch_factor_1 - distance / dz(n)

             IF ( numerator > 0.0_wp )  THEN
                l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0_wp
                l_rounded = NINT( l )
                delta_l = ABS( l_rounded - l ) / l
             ENDIF

             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )

             delta_stretch_factor = ABS( stretch_factor_1 - stretch_factor_2 ) / stretch_factor_2

             delta_total_new = delta_l + delta_stretch_factor

!
!--          stretch_factor_1 is taken to guarantee that the stretching procedure ends as close as
!--          possible to dz_stretch_level_end.
!--          stretch_factor_2 would guarantee that the stretched dz(n) is equal to dz(n+1) after
!--          l_rounded grid levels.
             IF (delta_total_new < delta_total_old)  THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF

             iterations = iterations + 1

          ENDDO

!
!--    Second branch for stretching from fine to rough
       ELSEIF ( dz(n) < dz(n+1) )  THEN

          DO WHILE ( stretch_factor_1 <= stretch_factor_upper_limit )

             stretch_factor_1 = 1.0_wp + iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) - dz_stretch_level_start(n) )
             numerator = distance * stretch_factor_1 / dz(n) + stretch_factor_1 - distance / dz(n)

             l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0_wp
             l_rounded = NINT( l )
             delta_l = ABS( l_rounded - l ) / l

             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )

             delta_stretch_factor = ABS( stretch_factor_1 - stretch_factor_2 ) / stretch_factor_2

             delta_total_new = delta_l + delta_stretch_factor

!
!--          stretch_factor_1 is taken to guarantee that the stretching procedure ends as close as
!--          possible to dz_stretch_level_end.
!--          stretch_factor_2 would guarantee that the stretched dz(n) is equal to dz(n+1) after
!--          l_rounded grid levels.
             IF (delta_total_new < delta_total_old)  THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF

             iterations = iterations + 1
          ENDDO

       ELSE

          message_string= 'adjacent values of dz have equal values'
          CALL message( 'init_grid', 'PAC0017', 1, 2, 0, 6, 0 )

       ENDIF

!
!--    Check if also the second stretching factor fits into the allowed interval. If not, print a
!--    warning for the user.
       IF ( dz_stretch_factor_array_2(n) < stretch_factor_lower_limit  .OR.                        &
            dz_stretch_factor_array_2(n) > stretch_factor_upper_limit )  THEN
          WRITE( message_string, * ) 'stretch_factor_2 = ', dz_stretch_factor_array_2(n),          &
                                     ' which is responsible for exactly reaching& dz =',           &
                                      dz(n+1), 'after a specific amount of',                       &
                                     ' grid levels& exceeds the upper',                            &
                                     ' limit =', stretch_factor_upper_limit,                       &
                                     ' &or lower limit = ', stretch_factor_lower_limit
          CALL message( 'init_grid', 'PAC0218', 0, 1, 0, 6, 0 )

       ENDIF

    ENDDO

 END SUBROUTINE calculate_stretching_factor
