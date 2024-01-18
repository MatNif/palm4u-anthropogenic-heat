!> @file src/inifor_transform.f90
!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2017-2021 Leibniz Universitaet Hannover
! Copyright 2017-2021 Deutscher Wetterdienst Offenbach
!------------------------------------------------------------------------------!
!
! Authors:
! --------
!> @author Eckhard Kadasch (Deutscher Wetterdienst, Offenbach)
!
! Description:
! ------------
!> The transform module provides INIFOR's low-level transformation and
!> interpolation routines. The rotated-pole transformation routines phirot2phi,
!> phi2phirot, rlarot2rla, rla2rlarot, uv2uvrot, and uvrot2uv are adapted from
!> int2lm's utility routines.
!------------------------------------------------------------------------------!
#if defined ( __netcdf )
 MODULE inifor_transform

    USE inifor_control
    USE inifor_defs,                                                           &
        ONLY: BETA, G, P_SL, PI, RD, T_SL, TO_DEGREES, TO_RADIANS, WATER_ID,   &
              iwp, wp
    USE inifor_types
    USE inifor_util,                                                           &
        ONLY: get_basic_state, real_to_str, str

    IMPLICIT NONE

 CONTAINS


 SUBROUTINE interpolate_1d(in_arr, out_arr, outgrid)
    TYPE(grid_definition), INTENT(IN) ::  outgrid
    REAL(wp), INTENT(IN)              ::  in_arr(:)
    REAL(wp), INTENT(OUT)             ::  out_arr(:)

    INTEGER(iwp) :: k, l, nz

    nz = UBOUND(out_arr, 1)

    DO  k = nz, LBOUND(out_arr, 1), -1

!
!--    TODO: Remove IF clause and extrapolate based on a critical vertical
!--    TODO: index marking the lower bound of COSMO-DE data coverage.
!--    Check for negative interpolation weights indicating grid points 
!--    below COSMO-DE domain and extrapolate from the top in such cells.
       IF (outgrid%w(1,k,1) < -1.0_wp .AND. k < nz)  THEN
          out_arr(k) = out_arr(k+1)
       ELSE
          out_arr(k) = 0.0_wp
          DO  l = 1, 2
             out_arr(k) = out_arr(k) +                                      &
                 outgrid%w(1,k,l) * in_arr(outgrid%kkk(1,k,l) )
          ENDDO
       ENDIF
    ENDDO

 END SUBROUTINE interpolate_1d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates linearly in the vertical direction in very column (i,j) of the
!> output array outvar(i,j,:) using values of the source array invar. In cells
!> that are outside the COSMO-DE domain, indicated by negative interpolation
!> weights, extrapolate constantly from the cell above.
!> 
!> Input parameters:
!> -----------------
!> invar : Array of source data
!> 
!> outgrid%kk : Array of vertical neighbour indices. kk(i,j,k,:) contain the
!>     indices of the two vertical neighbors of PALM-4U point (i,j,k) on the
!>     input grid corresponding to the source data invar.
!> 
!> outgrid%w_verti : Array of weights for vertical linear interpolation
!>     corresponding to neighbour points indexed by kk.
!>
!> Output papameters:
!> ------------------
!> outvar : Array of interpolated data
!------------------------------------------------------------------------------!
 SUBROUTINE interpolate_1d_arr(in_arr, out_arr, outgrid)
    TYPE(grid_definition), INTENT(IN) ::  outgrid
    REAL(wp), INTENT(IN)              ::  in_arr(0:,0:,0:)
    REAL(wp), INTENT(OUT)             ::  out_arr(0:,0:,:)

    INTEGER(iwp) :: i, j, k, l, nz

    nz = UBOUND(out_arr, 3)

    DO  j = LBOUND(out_arr, 2), UBOUND(out_arr, 2)
    DO  i = LBOUND(out_arr, 1), UBOUND(out_arr, 1)
    DO  k = nz, LBOUND(out_arr, 3), -1

!
!--    TODO: Remove IF clause and extrapolate based on a critical vertical
!--    TODO: index marking the lower bound of COSMO-DE data coverage.
!--    Check for negative interpolation weights indicating grid points 
!--    below COSMO-DE domain and extrapolate from the top in such cells.
       IF (outgrid%w_verti(i,j,k,1) < -1.0_wp .AND. k < nz)  THEN
          out_arr(i,j,k) = out_arr(i,j,k+1)
       ELSE
          out_arr(i,j,k) = 0.0_wp
          DO  l = 1, 2
             out_arr(i,j,k) = out_arr(i,j,k) +                              &
                 outgrid%w_verti(i,j,k,l) *                               &
                 in_arr(i,j,outgrid%kk(i,j,k, l) )
          ENDDO
       ENDIF
    ENDDO
    ENDDO
    ENDDO
 END SUBROUTINE interpolate_1d_arr


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates bi-linearly in horizontal planes on every k level of the output
!> array outvar(:,:,k) using values of the source array invar(:,:,:). The source
!> (invar) and interpolation array (outvar) need to have matching dimensions.
!>
!> Input parameters:
!> -----------------
!> invar : Array of source data
!> 
!> outgrid%ii,%jj : Array of neighbour indices in x and y direction.
!>     ii(i,j,k,:), and jj(i,j,k,:) contain the four horizontal neighbour points
!>     of PALM-4U point (i,j,k) on the input grid corresponding to the source
!>     data invar. (The outgrid carries the relationship with the ingrid in the
!      form of the interpolation weights.)
!> 
!> outgrid%w_horiz: Array of weights for horizontal bi-linear interpolation
!>     corresponding to neighbour points indexed by ii and jj.
!>
!> Output papameters:
!> ------------------
!> outvar : Array of interpolated data
!------------------------------------------------------------------------------!
 SUBROUTINE interpolate_2d(invar, outvar, outgrid, ncvar)
!
!-- I index 0-based for the indices of the outvar to be consistent with the
!-- outgrid indices and interpolation weights.
    TYPE(grid_definition), INTENT(IN)  ::  outgrid
    REAL(wp), INTENT(IN)               ::  invar(0:,0:,0:)
    REAL(wp), INTENT(OUT)              ::  outvar(0:,0:,0:)
    TYPE(nc_var), INTENT(IN), OPTIONAL ::  ncvar

    INTEGER(iwp) ::  i, j, k, l

!
!-- TODO: check if input dimensions are consistent, i.e. ranges are correct
    IF ( UBOUND(outvar, 3) .GT. UBOUND(invar, 3) )  THEN
        message = "Output array for '" // TRIM(ncvar%name) // "' has ' more levels (" // &
           TRIM(str(UBOUND(outvar, 3, kind=iwp))) // ") than input variable ("//&
           TRIM(str(UBOUND(invar, 3, kind=iwp))) // ")."
        CALL inifor_abort('interpolate_2d', message)
    ENDIF

    DO  k = 0, UBOUND(outvar, 3)
    DO  j = 0, UBOUND(outvar, 2)
    DO  i = 0, UBOUND(outvar, 1)
       outvar(i,j,k) = 0.0_wp
       DO  l = 1, 4
          
          outvar(i,j,k) = outvar(i,j,k) +                                      &
             outgrid%w_horiz(i,j,l) * invar( outgrid%ii(i,j,l),                &
                                               outgrid%jj(i,j,l),              &
                                               k )
       ENDDO
    ENDDO
    ENDDO
    ENDDO
        
 END SUBROUTINE interpolate_2d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the horizontal average of the in_arr(:,:,:) and store it in
!> out_arr(:)
!------------------------------------------------------------------------------!
 SUBROUTINE average_2d(in_arr, out_arr, ii, jj)
    REAL(wp), INTENT(IN)                   ::  in_arr(0:,0:,0:)
    REAL(wp), INTENT(OUT)                  ::  out_arr(0:)
    INTEGER(iwp), INTENT(IN), DIMENSION(:) ::  ii, jj

    INTEGER(iwp) ::  i, j, k, l
    REAL(wp)     ::  ni

    IF (SIZE(ii) /= SIZE(jj))  THEN
       message = "Length of 'ii' and 'jj' index lists do not match." //     &
          NEW_LINE(' ') // "ii has " // str(SIZE(ii, kind=iwp)) // " elements, " //   &
          NEW_LINE(' ') // "jj has " // str(SIZE(jj, kind=iwp)) // "."
       CALL inifor_abort('average_2d', message)
    ENDIF

    IF (SIZE(ii) == 0)  THEN
       message = "No columns to average over; " //                          &
                 "size of index lists 'ii' and 'jj' is zero."
       CALL inifor_abort('average_2d', message)
    ENDIF

    DO  k = 0, UBOUND(out_arr, 1)

       out_arr(k) = 0.0_wp
       DO  l = 1, UBOUND(ii, 1)
          i = ii(l)
          j = jj(l)
          out_arr(k) = out_arr(k) + in_arr(i, j, k)
       ENDDO

    ENDDO

    ni = 1.0_wp / SIZE(ii)
    out_arr(:) = out_arr(:) * ni

 END SUBROUTINE average_2d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Three-dimensional interpolation driver. Interpolates from the source_array to
!> the given palm_grid.
!>
!> The routine separates horizontal and vertical
!> interpolation. In the horizontal interpolation step, the source_array data is
!> interpolated along COSMO grid levels onto the intermediate grid (vertically
!> as coarse as COSMO, horizontally as fine as PALM).
!------------------------------------------------------------------------------!
 SUBROUTINE interpolate_3d(source_array, palm_array, palm_intermediate, palm_grid)
    TYPE(grid_definition), INTENT(IN) ::  palm_intermediate, palm_grid
    REAL(wp), DIMENSION(:,:,:), INTENT(IN)  ::  source_array
    REAL(wp), DIMENSION(:,:,:), INTENT(OUT) ::  palm_array
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  intermediate_array
    INTEGER(iwp) ::  nx, ny, nlev

    nx = palm_intermediate%nx
    ny = palm_intermediate%ny
    nlev = palm_intermediate%nz

!
!-- Interpolate from COSMO to intermediate grid. Allocating with one
!-- less point in the vertical, since scalars like T have 50 instead of 51
!-- points in COSMO.
    ALLOCATE(intermediate_array(0:nx, 0:ny, 0:nlev-1)) !

    CALL interpolate_2d(source_array, intermediate_array, palm_intermediate)

!
!-- Interpolate from intermediate grid to palm_grid grid, includes
!-- extrapolation for cells below COSMO domain.
    CALL interpolate_1d_arr(intermediate_array, palm_array, palm_grid)

    DEALLOCATE(intermediate_array)

 END SUBROUTINE interpolate_3d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Average data horizontally from the source_array over the region given by the
!> averaging grid 'avg_grid' and store the result in 'profile_array'.
!------------------------------------------------------------------------------!
 SUBROUTINE interp_average_profile(source_array, profile_array, avg_grid)
    TYPE(grid_definition), INTENT(IN)      ::  avg_grid
    REAL(wp), DIMENSION(:,:,:), INTENT(IN) ::  source_array(0:,0:,:)
    REAL(wp), DIMENSION(:), INTENT(OUT)    ::  profile_array

    INTEGER(iwp) ::  i_source, j_source, k_profile, k_source, l, m

    REAL ::  ni_columns

    profile_array(:) = 0.0_wp

    DO  l = 1, avg_grid%n_columns
       i_source = avg_grid%iii(l)
       j_source = avg_grid%jjj(l)

!
!--    Loop over PALM levels
       DO  k_profile = avg_grid%k_min, UBOUND(profile_array, 1)

!
!--       Loop over vertical interpolation neighbours
          DO  m = 1, 2

             k_source = avg_grid%kkk(l, k_profile, m)

             profile_array(k_profile) = profile_array(k_profile)            &
                + avg_grid%w(l, k_profile, m)                             &
                * source_array(i_source, j_source, k_source)
!
!--       Loop over vertical interpolation neighbours m
          ENDDO

!
!--    Loop over PALM levels k_profile
       ENDDO

!
!-- Loop over horizontal neighbours l
    ENDDO

    ni_columns = 1.0_wp / avg_grid%n_columns
    profile_array(:) = profile_array(:) * ni_columns

!
!-- Constant extrapolation to the bottom
    profile_array(1:avg_grid%k_min-1) = profile_array(avg_grid%k_min)

 END SUBROUTINE interp_average_profile


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Average data horizontally from the source_array over the region given by the
!> averaging grid 'avg_grid' and store the result in 'profile_array'.
!------------------------------------------------------------------------------!
 SUBROUTINE average_profile( source_array, profile_array, avg_grid )

    TYPE(grid_definition), INTENT(IN) ::  avg_grid
    REAL(wp), INTENT(IN)              ::  source_array(0:,0:,:)
    REAL(wp), INTENT(OUT)             ::  profile_array(:)

    INTEGER(iwp) ::  i_source, j_source, l, nz, nlev

    REAL(wp) ::  ni_columns

    nlev = SIZE( source_array, 3 )
    nz   = SIZE( profile_array, 1 )

    IF ( nlev /= nz )  THEN
       message = "Lengths of input and output profiles do not match: " //   &
                 "cosmo_pressure(" // TRIM( str( nlev ) ) //                &
                 "), profile_array(" // TRIM( str( nz ) )  // ")."
       CALL inifor_abort( 'average_profile', message )
    ENDIF
    
    profile_array(:) = 0.0_wp

    DO  l = 1, avg_grid%n_columns

       i_source = avg_grid%iii(l)
       j_source = avg_grid%jjj(l)

       profile_array(:) = profile_array(:)                                  &
                        + source_array(i_source, j_source, :)

    ENDDO

    ni_columns = 1.0_wp / avg_grid%n_columns
    profile_array(:) = profile_array(:) * ni_columns

 END SUBROUTINE average_profile


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This is a sister routine to average_profile() and differes from it in that
!> it removes the COSMO basic state pressure from the input array before
!> averaging.
!------------------------------------------------------------------------------!
 SUBROUTINE average_pressure_perturbation( cosmo_pressure, profile_array,   &
                                           cosmo_grid, avg_grid )

    TYPE(grid_definition), INTENT(IN) ::  cosmo_grid, avg_grid
    REAL(wp), INTENT(IN)              ::  cosmo_pressure(0:,0:,:)
    REAL(wp), INTENT(OUT)             ::  profile_array(:)

    INTEGER(iwp)          ::  i_source, j_source, l, nz, nlev
    REAL(wp)              ::  ni_columns
    REAL(wp), ALLOCATABLE ::  basic_state_pressure(:)

    nlev = SIZE( cosmo_pressure, 3 )
    nz   = SIZE( profile_array, 1 )

    IF ( nlev /= nz )  THEN
       message = "Lengths of input and output profiles do not match: " //   &
                 "cosmo_pressure(" // TRIM( str( nlev ) ) //                &
                 "), profile_array(" // TRIM( str( nz ) )  // ")."
       CALL inifor_abort('average_pressure_perturbation', message)
    ENDIF

    ALLOCATE( basic_state_pressure(nz) )
    profile_array(:) = 0.0_wp

    DO  l = 1, avg_grid%n_columns
       i_source = avg_grid%iii(l)
       j_source = avg_grid%jjj(l)

!
!--    Compute pressure perturbation by removing COSMO basic state pressure
       CALL get_basic_state( cosmo_grid%hfl(i_source,j_source,:), BETA,   &
                             P_SL, T_SL, RD, G, basic_state_pressure )

       profile_array(:) = profile_array(:)                                  &
                        + cosmo_pressure(i_source, j_source, :)             &
                        - basic_state_pressure(:)

!
!-- Loop over horizontal neighbours l
    ENDDO

    DEALLOCATE( basic_state_pressure )

    ni_columns = 1.0_wp / avg_grid%n_columns
    profile_array(:) = profile_array(:) * ni_columns

 END SUBROUTINE average_pressure_perturbation


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes average soil values in COSMO-DE water cells from neighbouring
!> non-water cells. This is done as a preprocessing step for the COSMO-DE
!> soil input arrays, which contain unphysical values for water cells.
!>
!> This routine computes the average of up to all nine neighbouring cells
!> or keeps the original value, if not at least one non-water neightbour
!> is available.
!>
!> By repeatedly applying this step, soil data can be extrapolated into
!> 'water' regions occupying multiple cells, one cell per iteration.
!> 
!> Input parameters:
!> -----------------
!> soiltyp : 2d map of COSMO-DE soil types
!> nz : number of layers in the COSMO-DE soil
!> niter : number iterations
!> 
!> Output parameters:
!> ------------------
!> array : the soil array (i.e. water content or temperature)
!------------------------------------------------------------------------------!
 SUBROUTINE fill_water_cells(soiltyp, array, nz, niter)
    INTEGER(iwp), DIMENSION(:,:,:), INTENT(IN) :: soiltyp
    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT)  :: array
    INTEGER(iwp), INTENT(IN)                   :: nz, niter

    REAL(wp), DIMENSION(nz)                    :: column
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE  :: old_soiltyp, new_soiltyp
    INTEGER(iwp)                               :: l, i, j, nx, ny, n_cells, ii, jj, iter
    INTEGER(iwp), DIMENSION(8)                 :: di, dj

    nx = SIZE(array, 1)
    ny = SIZE(array, 2)
    di = (/ -1, -1, -1, 0,  0,  1, 1, 1 /)
    dj = (/ -1,  0,  1, -1, 1, -1, 0, 1 /)

    ALLOCATE(old_soiltyp(SIZE(soiltyp,1), &
                         SIZE(soiltyp,2) ))

    ALLOCATE(new_soiltyp(SIZE(soiltyp,1), &
                         SIZE(soiltyp,2) ))

    old_soiltyp(:,:) = soiltyp(:,:,1)
    new_soiltyp(:,:) = soiltyp(:,:,1)

    DO  iter = 1, niter

       DO  j = 1, ny
       DO  i = 1, nx
       
          IF (old_soiltyp(i,j) == WATER_ID)  THEN 

             n_cells = 0
             column(:) = 0.0_wp
             DO  l = 1, SIZE(di)

                ii = MIN(nx, MAX(1_iwp, i + di(l)))
                jj = MIN(ny, MAX(1_iwp, j + dj(l)))

                IF (old_soiltyp(ii,jj) .NE. WATER_ID)  THEN
                   n_cells = n_cells + 1
                   column(:) = column(:) + array(ii,jj,:)
                ENDIF

             ENDDO

!
!--          Overwrite if at least one non-water neighbour cell is available
             IF (n_cells > 0)  THEN
                array(i,j,:) = column(:) / n_cells
                new_soiltyp(i,j) = 0
             ENDIF

          ENDIF

       ENDDO
       ENDDO

       old_soiltyp(:,:) = new_soiltyp(:,:)

    ENDDO

    DEALLOCATE(old_soiltyp, new_soiltyp)

 END SUBROUTINE fill_water_cells


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Takes the averaged pressure profile p_palm and sets the lowest entry to the
!> extrapolated pressure at the surface, assuming a linear density profile.
!------------------------------------------------------------------------------!
 SUBROUTINE get_surface_pressure(p_palm, rho_cosmo, avg_grid)
    REAL(wp), DIMENSION(:), INTENT(IN)    ::  rho_cosmo
    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  p_palm
    TYPE(grid_definition), INTENT(IN)     ::  avg_grid

    REAL(wp)     ::  drhodz_surface_cosmo
    INTEGER(iwp) ::  k_min_palm

    drhodz_surface_cosmo =                                                     &
       ( rho_cosmo(2) - rho_cosmo(1) ) /                                       &
       ( avg_grid%intermediate_h(1,1,2) - avg_grid%intermediate_h(1,1,1) )

    k_min_palm = avg_grid%k_min

    p_palm(1) = linear_density_pressure(                                       &
       p0 = p_palm(k_min_palm),                                                &
       z0 = avg_grid%z(k_min_palm),                                            &
       rho00 = rho_cosmo(1),                                                   &
       z00 = avg_grid%intermediate_h(1,1,1),                                   &
       drhodz = drhodz_surface_cosmo,                                          &
       g = G,                                                                  &
       z = 0.0_wp                                                              &
    )

 END SUBROUTINE get_surface_pressure


!------------------------------------------------------------------------------!
! Description:
! -----------
!> Computes the hydrostatic pressure p at height z given the pressure p0 at
!> height z0. The pressure is computed based on the solution of the hydrostatic
!> equation assuming a linear density profile with density rho00 at z00 and the
!> vertical density gradient drhodz. 
!------------------------------------------------------------------------------!
 FUNCTION linear_density_pressure(p0, z0, rho00, z00, drhodz, g, z)  RESULT(p)

    REAL(wp), INTENT(IN)  ::  p0, z0, rho00, z00, drhodz, g, z
    REAL(wp) ::  p

    p = p0 - ( z - z0 ) * g * (                                                &
           rho00 + 0.5_wp * drhodz * ( z + z0 - 2.0_wp * z00 )                 &
        )

 END FUNCTION linear_density_pressure

!------------------------------------------------------------------------------!
! Description:
! -----------
!> This routine computes profiles of the two geostrophic wind components ug and
!> vg.
!------------------------------------------------------------------------------!
 SUBROUTINE geostrophic_winds(p_north, p_south, p_east, p_west, rho, f3,    &
                              Lx, Ly, phi_n, lam_n, phi_g, lam_g, ug, vg)

    REAL(wp), DIMENSION(:), INTENT(IN)  ::  p_north, p_south, p_east, p_west,  &
                                            rho
    REAL(wp), INTENT(IN)                ::  f3, Lx, Ly, phi_n, lam_n, phi_g, lam_g
    REAL(wp), DIMENSION(:), INTENT(OUT) ::  ug, vg
    REAL(wp)                            ::  facx, facy

    facx = 1.0_wp / (Lx * f3)
    facy = 1.0_wp / (Ly * f3)
    ug(:) = - facy / rho(:) * (p_north(:) - p_south(:))
    vg(:) =   facx / rho(:) * (p_east(:) - p_west(:))

    CALL rotate_vector_field(                                                  &
       ug, vg, angle = meridian_convergence_rotated(phi_n, lam_n, phi_g, lam_g)&
    )

 END SUBROUTINE geostrophic_winds


!------------------------------------------------------------------------------!
! Description:
! -----------
!> This routine computes the inverse Plate Carree projection, i.e. in projects
!> Cartesian coordinates (x,y) onto a sphere. It returns the latitude and
!> lngitude of a geographical system centered at x0 and y0.
!------------------------------------------------------------------------------!
 SUBROUTINE inv_plate_carree(x, y, x0, y0, r, lat, lon)
    REAL(wp), INTENT(IN)  ::  x(:), y(:), x0, y0, r
    REAL(wp), INTENT(OUT) ::  lat(:), lon(:)
    
    REAL(wp) :: ri

!
!-- TODO check dimensions of lat/lon and y/x match

    ri = 1.0_wp / r
    
    lat(:) = (y(:) - y0) * ri
    lon(:) = (x(:) - x0) * ri
 END SUBROUTINE 


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the reverse Plate-Carree projection of a x or y position on a
!> Cartesian grid.
!>
!> Input parameters:
!> -----------------
!> xy : x or y coordinate of the Cartasian grid point [m].
!>
!> xy0 : x or y coordinate that coincides with the origin of the underlying
!>     sperical system (crossing point of the equator and prime meridian) [m].
!>
!> r : Radius of the of the underlying sphere, e.g. EARTH_RADIUS [m].
!> 
!> Returns:
!> --------
!> project : Longitude (in case xy = x) or latitude (xy = y) of the given input
!>     coordinate xy.
!------------------------------------------------------------------------------!
 ELEMENTAL REAL(wp) FUNCTION project(xy, xy0, r)
    REAL(wp), INTENT(IN)  ::  xy, xy0, r
    REAL(wp) :: ri

!
!-- If this elemental function is called with a large array as xy, it is
!-- computationally more efficient to precompute the inverse radius and
!-- then muliply.
    ri = 1.0_wp / r

    project = (xy - xy0) * ri

 END FUNCTION project


!------------------------------------------------------------------------------!
! Description:
! ------------
!> For a rotated-pole system with the origin at geographical latitude 'phi_c',
!> compute the geographical latitude of its rotated north pole.
!------------------------------------------------------------------------------!
 REAL(wp) FUNCTION phic_to_phin(phi_c)
     REAL(wp), INTENT(IN) ::  phi_c

     phic_to_phin = 0.5_wp * PI - ABS(phi_c)

 END FUNCTION phic_to_phin

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> For a rotated-pole system with the origin at geographical latitude 'phi_c'
!> and longitude 'lam_c', compute the geographical longitude of its rotated
!> north pole.
!------------------------------------------------------------------------------!
 REAL(wp) FUNCTION lamc_to_lamn(phi_c, lam_c)
    REAL(wp), INTENT(IN) ::  phi_c, lam_c
     
    lamc_to_lamn = lam_c
    IF (phi_c > 0.0_wp)  THEN
       lamc_to_lamn = lam_c - SIGN(PI, lam_c)
    ENDIF

 END FUNCTION lamc_to_lamn


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set gamma according to whether PALM domain is in the northern or southern
!> hemisphere of the COSMO rotated-pole system. Gamma assumes either the
!> value 0 or PI and is needed to work around around a bug in the
!> rotated-pole coordinate transformations.
!------------------------------------------------------------------------------!
 REAL(wp) FUNCTION gamma_from_hemisphere(phi_cg, phi_ref)
    REAL(wp), INTENT(IN) ::  phi_cg
    REAL(wp), INTENT(IN) ::  phi_ref

    LOGICAL ::  palm_origin_is_south_of_cosmo_origin
    
    palm_origin_is_south_of_cosmo_origin = (phi_cg < phi_ref)

    IF (palm_origin_is_south_of_cosmo_origin)  THEN
        gamma_from_hemisphere = PI
    ELSE
        gamma_from_hemisphere = 0.0_wp
    ENDIF
 END FUNCTION gamma_from_hemisphere


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes the geographical coordinates corresponding to the given rotated-pole
!> coordinates.
!>
!> In INIFOR, this routine is used to convert coordinates between two
!> rotated-pole systems: COSMO-DE's rotated-pole system, and one centred at the
!> PALM-4U domain centre. In this case, the PALM-4U system is thought of as the
!> rotated-pole system and the routine is used to rotate back to COSMO-DE's
!> system which is thought of as the geographical one.
!> 
!> Input parameters:
!> -----------------
!> phir(:), lamr(: ): latitudes and longitudes of the rotated-pole grid
!> 
!> phip, lamp: latitude and longitude of the rotated north pole
!>
!> gam: "angle between the north poles. If [gam] is not present, the other
!>       system is the real geographical system." (original phiro2rot
!>       description)
!> 
!> Output parameters:
!> ------------------
!> phi(:,:), lam(:,:): geographical latitudes and logitudes
!------------------------------------------------------------------------------!
 SUBROUTINE rotate_to_cosmo(phir, lamr, phip, lamp, phi, lam, gam)
    REAL(wp), INTENT(IN)  ::  phir(0:), lamr(0:), phip, lamp, gam
    REAL(wp), INTENT(OUT) ::  phi(0:,0:), lam(0:,0:)

    INTEGER(iwp) ::  i, j
    
    IF ( SIZE(phi, 1) .NE. SIZE(lam, 1) .OR. &
         SIZE(phi, 2) .NE. SIZE(lam, 2) )  THEN
       PRINT *, "inifor: rotate_to_cosmo: Dimensions of phi and lambda do not match. Dimensions are:"
       PRINT *, "inifor: rotate_to_cosmo: phi: ", SIZE(phi, 1), SIZE(phi, 2)
       PRINT *, "inifor: rotate_to_cosmo: lam: ", SIZE(lam, 1), SIZE(lam, 2)
       STOP
    ENDIF

    IF ( SIZE(phir) .NE. SIZE(phi, 2) .OR. &
         SIZE(lamr) .NE. SIZE(phi, 1) )  THEN
       PRINT *, "inifor: rotate_to_cosmo: Dimensions of phir and lamr do not match. Dimensions are:"
       PRINT *, "inifor: rotate_to_cosmo: phir: ", SIZE(phir), SIZE(phi, 2)
       PRINT *, "inifor: rotate_to_cosmo: lamr: ", SIZE(lamr), SIZE(phi, 1)
       STOP
    ENDIF
    
    DO  j = 0, UBOUND(phir, 1)
       DO  i = 0, UBOUND(lamr, 1)

          phi(i,j) = phirot2phi(phir(j) * TO_DEGREES,                       &
                                lamr(i) * TO_DEGREES,                       &
                                phip * TO_DEGREES,                          &
                                gam  * TO_DEGREES) * TO_RADIANS

          lam(i,j) = rlarot2rla(phir(j) * TO_DEGREES,                       &
                                lamr(i) * TO_DEGREES,                       &
                                phip * TO_DEGREES,                          &
                                lamp * TO_DEGREES,                          &
                                gam  * TO_DEGREES) * TO_RADIANS

       ENDDO
    ENDDO

 END SUBROUTINE rotate_to_cosmo
       

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Rotate the given vector field (x(:), y(:)) by the given 'angle'.
!------------------------------------------------------------------------------!
 SUBROUTINE rotate_vector_field(x, y, angle)
    REAL(wp), DIMENSION(:), INTENT(INOUT) :: x, y  !< x and y coodrinate in arbitrary units
    REAL(wp), INTENT(IN)                  :: angle !< rotation angle [deg]

    INTEGER(iwp) :: i
    REAL(wp)     :: sine, cosine, v_rot(2), rotation(2,2)

    sine = SIN(angle * TO_RADIANS)
    cosine = COS(angle * TO_RADIANS)
!
!-- RESAHPE() fills columns first, so the rotation matrix becomes
!-- 
!-- rotation = [ cosine   -sine  ]
!--            [  sine    cosine ]
    rotation = RESHAPE( (/cosine, sine, -sine, cosine/), (/2, 2/) )

    DO  i = LBOUND(x, 1), UBOUND(x, 1)

       v_rot(:) = MATMUL(rotation, (/x(i), y(i)/))

       x(i) = v_rot(1)
       y(i) = v_rot(2)

    ENDDO

 END SUBROUTINE rotate_vector_field



!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine computes the local meridian convergence between a rotated-pole
!> and a geographical system using the Eq. (6) given in the DWD manual
!> 
!>    Baldauf et al. (2018), Beschreibung des operationelle Kürzestfrist-
!>    vorhersagemodells COSMO-D2 und COSMO-D2-EPS und seiner Ausgabe in die
!>    Datenbanken des DWD.
!>    https://www.dwd.de/SharedDocs/downloads/DE/modelldokumentationen/nwv/cosmo_d2/cosmo_d2_dbbeschr_aktuell.pdf?__blob=publicationFile&v=2
!------------------------------------------------------------------------------!
 FUNCTION meridian_convergence_rotated(phi_n, lam_n, phi_g, lam_g)          &
    RESULT(delta)

    REAL(wp), INTENT(IN) ::  phi_n, lam_n, phi_g, lam_g
    REAL(wp)             ::  delta

    delta = atan2( COS(phi_n) * SIN(lam_n - lam_g),                         &
                   COS(phi_g) * SIN(phi_n) -                                &
                   SIN(phi_g) * COS(phi_n) * COS(lam_n - lam_g) )

 END FUNCTION meridian_convergence_rotated

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute indices of PALM-4U grid point neighbours in the target 
!> system (COSMO-DE) by rounding up and down. (i,j) are the indices of
!> the PALM-4U grid and (ii(i,j,1-4), jj(i,j,1-4)) contain the indices
!> of the its four neigbouring points in the COSMO-DE grid.
!>
!>
!>                     COSMO-DE grid
!>                     -------------
!>           jj, lat
!>              ^        j
!>              |         \          i
!>  jj(i,j,2/3) + ... 2 ---\--------/------ 3
!>              |     | ^   \      /        |
!>              |     | |wp  \    /         |
!>              |     | v     \  /          |
!>       latpos + ............ o/ (i,j)     |
!>              |     |        :            |
!>              |     |        :<----wl---->|
!>  jj(i,j,1/4) + ... 1 -------:----------- 4
!>              |     :        :            :
!>              |     :        :            :
!>              |     :      lonpos         :
!>              L-----+--------+------------+------> ii, lon 
!>               ii(i,j,1/2)        ii(i,j,3/4)
!>
!> 
!> Input parameters:
!> -----------------
!> source_lat, source_lon : (rotated-pole) coordinates of the source grid (e.g.
!>    COSMO-DE)
!>
!> source_dxi, source_dyi : inverse grid spacings of the source grid.
!>
!> target_lat, target_lon : (rotated-pole) coordinates of the target grid (e.g.
!>    COSMO-DE)
!> 
!> Output parameters:
!> ------------------
!> palm_ii, palm_jj : x and y index arrays of horizontal neighbour columns
!> 
!------------------------------------------------------------------------------!
 SUBROUTINE find_horizontal_neighbours(cosmo_lat, cosmo_lon,                &
                                       palm_clat, palm_clon,                &
                                       palm_ii, palm_jj)

    REAL(wp), DIMENSION(0:), INTENT(IN)            ::  cosmo_lat, cosmo_lon
    REAL(wp), DIMENSION(0:,0:), INTENT(IN)         ::  palm_clat, palm_clon
    REAL(wp)                                       ::  cosmo_dxi, cosmo_dyi
    INTEGER(iwp), DIMENSION(0:,0:,1:), INTENT(OUT) ::  palm_ii, palm_jj

    REAL(wp)     ::  lonpos, latpos, lon0, lat0
    INTEGER(iwp) ::  i, j

    lon0 = cosmo_lon(0)
    lat0 = cosmo_lat(0)
    cosmo_dxi = 1.0_wp / (cosmo_lon(1) - cosmo_lon(0))
    cosmo_dyi = 1.0_wp / (cosmo_lat(1) - cosmo_lat(0))

    DO  j = 0, UBOUND(palm_clon, 2)!palm_grid%ny
    DO  i = 0, UBOUND(palm_clon, 1)!palm_grid%nx
!
!--    Compute the floating point index corrseponding to PALM-4U grid point
!--    location along target grid (COSMO-DE) axes.
       lonpos = (palm_clon(i,j) - lon0) * cosmo_dxi
       latpos = (palm_clat(i,j) - lat0) * cosmo_dyi

       IF (lonpos < 0.0_wp .OR. latpos < 0.0_wp)  THEN
          message = "lonpos or latpos out of bounds " //                    &
             "while finding interpolation neighbours!" // NEW_LINE(' ') //  &
             "          (i,j) = (" //                                       &
             TRIM(str(i)) // ", " // TRIM(str(j)) // ")" // NEW_LINE(' ') //&
             "          lonpos " // TRIM(real_to_str(lonpos*TO_DEGREES)) // &
             ", latpos " // TRIM(real_to_str(latpos*TO_DEGREES)) // NEW_LINE(' ') // &
             "          lon0 " // TRIM(real_to_str(lon0  *TO_DEGREES)) //   &
             ", lat0   " // TRIM(real_to_str(lat0*TO_DEGREES)) // NEW_LINE(' ') // &
             "          PALM lon " // TRIM(real_to_str(palm_clon(i,j)*TO_DEGREES)) // &
             ", PALM lat " // TRIM(real_to_str(palm_clat(i,j)*TO_DEGREES))
          CALL inifor_abort('find_horizontal_neighbours', message)
       ENDIF

       palm_ii(i,j,1) = FLOOR(lonpos)
       palm_ii(i,j,2) = FLOOR(lonpos)
       palm_ii(i,j,3) = CEILING(lonpos)
       palm_ii(i,j,4) = CEILING(lonpos)

       palm_jj(i,j,1) = FLOOR(latpos)
       palm_jj(i,j,2) = CEILING(latpos)
       palm_jj(i,j,3) = CEILING(latpos)
       palm_jj(i,j,4) = FLOOR(latpos)
    ENDDO
    ENDDO

 END SUBROUTINE find_horizontal_neighbours

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes linear vertical interpolation neighbour indices and weights for each
!> column of the given palm grid.
!------------------------------------------------------------------------------!
 SUBROUTINE find_vertical_neighbours_and_weights_interp( palm_grid,            &
                                                         palm_intermediate )
    TYPE(grid_definition), INTENT(INOUT) ::  palm_grid
    TYPE(grid_definition), INTENT(IN)    ::  palm_intermediate

    INTEGER(iwp) ::  i, j, k, nx, ny, nz, nlev, k_intermediate
    LOGICAL      ::  point_is_below_grid, point_is_above_grid,                 &
                     point_is_in_current_cell
    REAL(wp)     ::  current_height, column_base, column_top, h_top, h_bottom, &
                     weight

    nx   = palm_grid%nx
    ny   = palm_grid%ny
    nz   = palm_grid%nz
    nlev = palm_intermediate%nz

!
!-- in each column of the fine grid, find vertical neighbours of every cell
    DO  j = 0, ny
    DO  i = 0, nx

       k_intermediate = 0

       column_base = palm_intermediate%intermediate_h(i,j,0)
       column_top  = palm_intermediate%intermediate_h(i,j,nlev)

!
!--    scan through palm_grid column and set neighbour indices in
!--    case current_height is either below column_base, in the current
!--    cell, or above column_top. Keep increasing current cell index until
!--    the current cell overlaps with the current_height.
       DO  k = 1, nz

!
!--       Memorize the top and bottom boundaries of the coarse cell and the
!--       current height within it
          current_height = palm_grid%z(k) + palm_grid%z0
          h_top    = palm_intermediate%intermediate_h(i,j,k_intermediate+1)
          h_bottom = palm_intermediate%intermediate_h(i,j,k_intermediate)

          point_is_above_grid = (current_height > column_top) !22000m, very unlikely
          point_is_below_grid = (current_height < column_base)

          point_is_in_current_cell = (                                      &
             current_height >= h_bottom .AND.                               &
             current_height <  h_top                                        &
          )

!
!--       set default weights
          palm_grid%w_verti(i,j,k,1:2) = 0.0_wp

          IF (point_is_above_grid)  THEN

             palm_grid%kk(i,j,k,1:2) = nlev
             palm_grid%w_verti(i,j,k,1:2) = - 2.0_wp

             message =                                                         &
                "PALM-4U grid extends above COSMO model top. " //              &
                "Height of COSMO model top: " //                               &
                TRIM( real_to_str( column_top, format = '(F11.3)') ) //        &
                " m, height of current PALM point: " //                        &
                TRIM( real_to_str( current_height , format = '(F11.3)') ) //   &
                " m, height of PALM model top: " //                            &
                TRIM( real_to_str( MAXVAL(palm_grid%z) + palm_grid%z0, '(F11.3)') ) // " m."
             CALL inifor_abort('find_vertical_neighbours_and_weights_interp', message)

          ELSE IF (point_is_below_grid)  THEN

             palm_grid%kk(i,j,k,1:2) = 0
             palm_grid%w_verti(i,j,k,1:2) = - 2.0_wp

          ELSE
!
!--          cycle through intermediate levels until current
!--          intermediate-grid cell overlaps with current_height
             DO  WHILE (.NOT. point_is_in_current_cell .AND. k_intermediate <= nlev-1)
                k_intermediate = k_intermediate + 1

                h_top    = palm_intermediate%intermediate_h(i,j,k_intermediate+1)
                h_bottom = palm_intermediate%intermediate_h(i,j,k_intermediate)
                point_is_in_current_cell = (                                &
                   current_height >= h_bottom .AND.                         &
                   current_height <  h_top                                  &
                )
             ENDDO

             IF (k_intermediate > nlev-1)  THEN
                message = "Index " // TRIM(str(k_intermediate)) //          &
                          " is above intermediate grid range."
                CALL inifor_abort('find_vertical_neighbours_interp', message)
             ENDIF
   
             palm_grid%kk(i,j,k,1) = k_intermediate
             palm_grid%kk(i,j,k,2) = k_intermediate + 1

!
!--          compute vertical weights
             weight = (h_top - current_height) / (h_top - h_bottom)
             palm_grid%w_verti(i,j,k,1) = weight
             palm_grid%w_verti(i,j,k,2) = 1.0_wp - weight
          ENDIF

       ENDDO

    ENDDO
    ENDDO

 END SUBROUTINE find_vertical_neighbours_and_weights_interp


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computes linear vertical interpolation neighbour indices and weights for each
!> column of the given averaging grid.
!>
!> The difference to the _interp variant of this routine lies in how columns
!> are adressed. While the _interp variant loops over all PALM grid columns
!> given by combinations of all index indices (i,j), this routine loops over a
!> subset of COSMO columns, which are sequentially stored in the index lists
!> iii(:) and jjj(:).
!------------------------------------------------------------------------------!
 SUBROUTINE find_vertical_neighbours_and_weights_average(                   &
    avg_grid, level_based_averaging                                         &
 )

    TYPE(grid_definition), INTENT(INOUT), TARGET ::  avg_grid
    LOGICAL                                      ::  level_based_averaging

    INTEGER(iwp)      ::  i, j, k_palm, k_intermediate, l, nlev
    LOGICAL           ::  point_is_below_grid, point_is_above_grid,         &
                          point_is_in_current_cell
    REAL(wp)          ::  current_height, column_base, column_top, h_top,   &
                          h_bottom, weight
    REAL(wp), POINTER ::  cosmo_h(:,:,:)


    avg_grid%k_min = LBOUND(avg_grid%z, 1)

    nlev = SIZE(avg_grid%cosmo_h, 3)

!
!-- For level-based averaging, use the profile of averaged vertical mesoscale
!-- levels computed in init_averaging_grid().
    IF (level_based_averaging)  THEN
       cosmo_h => avg_grid%intermediate_h
    ELSE
       cosmo_h => avg_grid%cosmo_h
    ENDIF

!
!-- in each column of the fine grid, find vertical neighbours of every cell
    DO  l = 1, avg_grid%n_columns

!--    The profile of averaged vertical mesoscale levels stored in
!--    intermediate_h only contains one column. By using the same column -- and
!--    consequently the same vertical interpolation neighbours and weights --
!--    
       IF (level_based_averaging)  THEN
          i = 1
          j = 1
       ELSE
          i = avg_grid%iii(l)
          j = avg_grid%jjj(l)
       ENDIF

       column_base = cosmo_h(i,j,1)
       column_top  = cosmo_h(i,j,nlev)

!
!--    Scan through avg_grid column until and set neighbour indices in
!--    case current_height is either below column_base, in the current
!--    cell, or above column_top. Keep increasing current cell index until
!--    the current cell overlaps with the current_height.
       k_intermediate = 1 !avg_grid%cosmo_h is indexed 1-based.
       DO  k_palm = 1, avg_grid%nz

!
!--       Memorize the top and bottom boundaries of the coarse cell and the
!--       current height within it
          current_height = avg_grid%z(k_palm) + avg_grid%z0
          h_top    = cosmo_h(i,j,k_intermediate+1)
          h_bottom = cosmo_h(i,j,k_intermediate)

!
!--       COSMO column top is located at 22000m, point_is_above_grid is very
!--       unlikely.
          point_is_above_grid = (current_height > column_top)
          point_is_below_grid = (current_height < column_base)

          point_is_in_current_cell = (                                      &
             current_height >= h_bottom .AND.                               &
             current_height <  h_top                                        &
          )

!
!--       set default weights
          avg_grid%w(l,k_palm,1:2) = 0.0_wp

          IF (point_is_above_grid)  THEN

             avg_grid%kkk(l,k_palm,1:2) = nlev
             avg_grid%w(l,k_palm,1:2) = - 2.0_wp

             message =                                                         &
                "PALM-4U grid extends above COSMO model top. " //              &
                "Height of COSMO model top: " //                               &
                TRIM( real_to_str( column_top, format = '(F11.3)') ) //        &
                " m, height of current PALM point: " //                        &
                TRIM( real_to_str( current_height , format = '(F11.3)') ) //   &
                " m, height of PALM model top: " //                            &
                TRIM( real_to_str( MAXVAL(avg_grid%z) + avg_grid%z0, '(F11.3)') ) // " m."
             CALL inifor_abort('find_vertical_neighbours_and_weights_average', message)

          ELSE IF (point_is_below_grid)  THEN

             avg_grid%kkk(l,k_palm,1:2) = 0
             avg_grid%w(l,k_palm,1:2) = - 2.0_wp
             avg_grid%k_min = MAX(k_palm + 1, avg_grid%k_min)
          ELSE
!
!--          cycle through intermediate levels until current
!--          intermediate-grid cell overlaps with current_height
             DO  WHILE (.NOT. point_is_in_current_cell .AND. k_intermediate <= nlev-1)
                k_intermediate = k_intermediate + 1

                h_top    = cosmo_h(i,j,k_intermediate+1)
                h_bottom = cosmo_h(i,j,k_intermediate)
                point_is_in_current_cell = (                                &
                   current_height >= h_bottom .AND.                         &
                   current_height <  h_top                                  &
                )
             ENDDO

!
!--          k_intermediate = 48 indicates the last section (indices 48 and 49), i.e.
!--          k_intermediate = 49 is not the beginning of a valid cell.
             IF (k_intermediate > nlev-1)  THEN
                message = "Index " // TRIM(str(k_intermediate)) //          &
                          " is above intermediate grid range."
                CALL inifor_abort('find_vertical_neighbours_average', message)
             ENDIF
   
             avg_grid%kkk(l,k_palm,1) = k_intermediate
             avg_grid%kkk(l,k_palm,2) = k_intermediate + 1

!
!--          compute vertical weights
             weight = (h_top - current_height) / (h_top - h_bottom)
             avg_grid%w(l,k_palm,1) = weight
             avg_grid%w(l,k_palm,2) = 1.0_wp - weight
          ENDIF

!
!--    Loop over PALM levels k
       ENDDO

!
!--    Loop over averaging columns l
    ENDDO
 
 END SUBROUTINE find_vertical_neighbours_and_weights_average


 SUBROUTINE compute_sigma( z, zs, zt, sigma )
     REAL(wp), INTENT(IN)  ::  z(:,:,:)     ! vertical grid levels
     REAL(wp), INTENT(IN)  ::  zs(:,:,:)    ! terrain elevation / z "surface"
     REAL(wp), INTENT(IN)  ::  zt           ! model top
     REAL(wp), INTENT(OUT) ::  sigma(:,:,:) ! sigma-z coordinate
     INTEGER ::  i,j,k

     DO k = 1, UBOUND(z,3)
     DO j = 1, UBOUND(z,2)
     DO i = 1, UBOUND(z,1)
        sigma(i,j,k) = sigma_from_z( z(i,j,k), zs(i,j,1), zt )
     END DO
     END DO
     END DO
 END SUBROUTINE compute_sigma



 REAL(wp) FUNCTION sigma_from_z( z, zs, zt )
     REAL(wp), INTENT(IN)  ::  z            ! vertical grid levels
     REAL(wp), INTENT(IN)  ::  zs           ! terrain elevation / z "surface"
     REAL(wp), INTENT(IN)  ::  zt           ! model top

     sigma_from_z = ( z - zt ) / ( zs - zt )

 END FUNCTION sigma_from_z


!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Shift horizontally interpolated levels of the intermediate grid
!> (intermediate_grid%intermediate_hh) so that they follow the microscale
!> terrain.
!------------------------------------------------------------------------------!
 SUBROUTINE map_terrain_driver( intermediate_grid, zs_prime,                   &
                                blending_dz, blending_z_ubound )

    TYPE(grid_definition), TARGET, INTENT(INOUT) ::  intermediate_grid
    REAL(wp), INTENT(IN)                         ::  zs_prime(:,:,:)
    REAL(wp), INTENT(IN)                         ::  blending_dz
    REAL(wp), INTENT(IN)                         ::  blending_z_ubound
    REAL(wp), POINTER                            ::  zs(:,:,:), z(:,:,:), sigma(:,:,:), zt

    sigma => intermediate_grid%sigma
    zs => intermediate_grid%zs
    zt => intermediate_grid%zt
    z => intermediate_grid%intermediate_h

    CALL compute_blending_function( z, intermediate_grid%z0, zs, zt, sigma,    &
                                    blending_dz, blending_z_ubound )

    CALL map_terrain( z, sigma, zs(:,:,1), zs_prime(:,:,1), zt )

    IF ( .NOT. is_monotone( z ) )  THEN
       message = 'Grid levels are not monoteone after terrain mapping.'
       CALL inifor_abort( 'map_terrain_driver', message )
    ENDIF

 END SUBROUTINE map_terrain_driver


 LOGICAL FUNCTION is_monotone( z )
    REAL(wp), INTENT(IN) ::  z(0:,0:,0:)
    INTEGER              ::  nz

    nz = UBOUND( z, 3 )
    is_monotone = ALL( (z(:,:,1:nz) - z(:,:,0:nz-1)) .GT. 0.0_wp )

 END FUNCTION is_monotone

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fills the height array z with values of an S-shaped blending function that is
!> used to blend between the mesoscale and microscale vertical grid spacing
!> scale (zs - zt) and (zs_prime - zt), respectively.
!------------------------------------------------------------------------------!
 SUBROUTINE compute_blending_function( z, z0, zs, zt, sigma, blending_dz,      &
                                       blending_z_ubound )
     REAL(wp), INTENT(INOUT) ::  z(0:,0:,0:)
     REAL(wp), INTENT(IN)    ::  z0
     REAL(wp), INTENT(IN)    ::  zs(0:,0:,:), zt
     REAL(wp), INTENT(IN)    ::  sigma(0:,0:,0:)
     REAL(wp), INTENT(IN)    ::  blending_dz
     REAL(wp), INTENT(IN)    ::  blending_z_ubound
     REAL(wp)                ::  blending_z_lbound
     REAL(wp)                ::  sigma_blending_ubound, sigma_blending_lbound
     INTEGER                 ::  i, j, k

     DO k = 0, UBOUND(sigma, 3)
     DO j = 0, UBOUND(sigma, 2)
     DO i = 0, UBOUND(sigma, 1)
        blending_z_lbound = zs(i,j,1) + blending_dz
        sigma_blending_ubound = sigma_from_z( blending_z_ubound, zs(i,j,1), zt )
        sigma_blending_lbound = sigma_from_z( blending_z_lbound, zs(i,j,1), zt )
        z(i,j,k) = dfunc_cosine(                                               &
           sigma(i,j,k), sigma_blending_ubound, sigma_blending_lbound          &
        )
     ENDDO
     ENDDO
     ENDDO
 END SUBROUTINE compute_blending_function


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine computes the vertically mapped intermediate grid with levels
!> shifted such that the lowest half layer matches the PALM terrain. The shift
!> of each level is governed by the value of a blending function at the
!> corresponding sigma-z height.
!>
!> This routine expects the blending function values to be stored in the z array
!> and overwrites it with the corresponding shifted level heights. The blending
!> funtion's values are precomputed using the compute_blending_funtion()
!> routine.
!------------------------------------------------------------------------------!
 SUBROUTINE map_terrain( z, sigma, zs, zs_prime, zt )

     REAL(wp), INTENT(INOUT) ::  z(0:,0:,0:)
     REAL(wp), INTENT(IN)    ::  sigma(0:,0:,0:), zs(0:,0:), zs_prime(0:,0:), zt
     INTEGER :: i, j, k

     DO k = 0, UBOUND(sigma, 3)
     DO j = 0, UBOUND(sigma, 2)
     DO i = 0, UBOUND(sigma, 1)

        z(i,j,k) = zt + sigma(i,j,k) * ( (1.0_wp - z(i,j,k)) * ( zs(i,j)       - zt ) &
                                                  + z(i,j,k) * ( zs_prime(i,j) - zt ) )

     ENDDO
     ENDDO
     ENDDO

 END SUBROUTINE map_terrain


 REAL(wp) FUNCTION dfunc_cosine( sigma, sigma_0, sigma_1 )
    REAL(wp), INTENT(IN) :: sigma
    REAL(wp), INTENT(IN) :: sigma_0, sigma_1

    ! upper part, no grid shift
    IF (sigma < sigma_0)  THEN
       dfunc_cosine = 0.0_wp
    ! surface
    ELSEIF (sigma > sigma_1)  THEN
       dfunc_cosine = 1.0_wp
    ! blending zone: sigma_0 >= sigma >= sigma_1
    ELSE
       dfunc_cosine = 0.5_wp * ( 1.0_wp - COS(                                 &
          PI * ( sigma - sigma_0 ) / ( sigma_1 - sigma_0 )                     &
       ) )
    ENDIF

 END FUNCTION dfunc_cosine


 REAL(wp) FUNCTION dfunc( sigma )
    REAL(wp), INTENT(IN) :: sigma

    dfunc = 1.0_wp

 END FUNCTION dfunc


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the four weights for horizontal bilinear interpolation given the
!> coordinates clon(i,j) clat(i,j) of the PALM-4U grid in the COSMO-DE
!> rotated-pole grid and the neightbour indices ii(i,j,1-4) and jj(i,j,1-4).
!>
!> Input parameters:
!> -----------------
!> palm_grid%clon : longitudes of PALM-4U scalars (cell centres) in COSMO-DE's rotated-pole grid [rad]
!>
!> palm_grid%clat : latitudes of PALM-4U cell centres in COSMO-DE's rotated-pole grid [rad]
!>
!> cosmo_grid%lon : rotated-pole longitudes of scalars (cell centres) of the COSMO-DE grid [rad] 
!>
!> cosmo_grid%lat : rotated-pole latitudes of scalars (cell centers) of the COSMO-DE grid [rad]
!>
!> cosmo_grid%dxi : inverse grid spacing in the first dimension [m^-1]
!>
!> cosmo_grid%dyi : inverse grid spacing in the second dimension [m^-1]
!>
!> Output parameters:
!> ------------------
!> palm_grid%w_horiz(:,:,1-4) : weights for bilinear horizontal interpolation
!
!                               COSMO-DE grid
!                               -------------
!                     jj, lat
!                        ^        j
!                        |         \          i
!            jj(i,j,2/3) + ... 2 ---\--------/------ 3
!                        |     | ^   \      /        |
!                        |     | |wp  \    /         |
!                        |     | v     \  /          |
!                 latpos + ............ o/ (i,j)     |
!                        |     |        :            |
!                        |     |        :<----wl---->|
!            jj(i,j,1/4) + ... 1 -------:----------- 4
!                        |     :        :            :
!                        |     :        :            :
!                        |     :      lonpos         :
!                        L-----+--------+------------+------> ii, lon 
!                         ii(i,j,1/2)        ii(i,j,3/4)
!          
!------------------------------------------------------------------------------!
 SUBROUTINE compute_horizontal_interp_weights(cosmo_lat, cosmo_lon,         &
    palm_clat, palm_clon, palm_ii, palm_jj, palm_w_horiz)
       
    REAL(wp), DIMENSION(0:), INTENT(IN)           ::  cosmo_lat, cosmo_lon
    REAL(wp)                                      ::  cosmo_dxi, cosmo_dyi
    REAL(wp), DIMENSION(0:,0:), INTENT(IN)        ::  palm_clat, palm_clon
    INTEGER(iwp), DIMENSION(0:,0:,1:), INTENT(IN) ::  palm_ii, palm_jj

    REAL(wp), DIMENSION(0:,0:,1:), INTENT(OUT) ::  palm_w_horiz

    REAL(wp)     ::  wlambda, wphi
    INTEGER(iwp) ::  i, j

    cosmo_dxi = 1.0_wp / (cosmo_lon(1) - cosmo_lon(0))
    cosmo_dyi = 1.0_wp / (cosmo_lat(1) - cosmo_lat(0))

    DO  j = 0, UBOUND(palm_clon, 2)
    DO  i = 0, UBOUND(palm_clon, 1)
    
!
!--    weight in lambda direction
       wlambda = ( cosmo_lon(palm_ii(i,j,4)) - palm_clon(i,j) ) * cosmo_dxi

!
!--    weight in phi direction
       wphi = ( cosmo_lat(palm_jj(i,j,2)) - palm_clat(i,j) ) * cosmo_dyi

       IF (wlambda > 1.0_wp .OR. wlambda < 0.0_wp)  THEN
           message = "Horizontal weight wlambda = " // TRIM(real_to_str(wlambda)) //   &
                     " is out bounds."
           CALL inifor_abort('compute_horizontal_interp_weights', message)
       ENDIF
       IF (wphi > 1.0_wp .OR. wphi < 0.0_wp)  THEN
           message = "Horizontal weight wphi = " // TRIM(real_to_str(wphi)) //   &
                     " is out bounds."
           CALL inifor_abort('compute_horizontal_interp_weights', message)
       ENDIF

       palm_w_horiz(i,j,1) = wlambda * wphi
       palm_w_horiz(i,j,2) = wlambda * (1.0_wp - wphi)
       palm_w_horiz(i,j,3) = (1.0_wp - wlambda) * (1.0_wp - wphi)
       palm_w_horiz(i,j,4) = 1.0_wp - SUM( palm_w_horiz(i,j,1:3) )

    ENDDO
    ENDDO
       
 END SUBROUTINE compute_horizontal_interp_weights


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolates u and v components of velocities located at cell faces to the
!> cell centres by averaging neighbouring values.
!>
!> This routine is designed to be used with COSMO-DE arrays where there are the
!> same number of grid points for scalars (centres) and velocities (faces). In
!> COSMO-DE the velocity points are staggared one half grid spaceing up-grid
!> which means the first centre point has to be omitted and is set to zero.
!------------------------------------------------------------------------------!
 SUBROUTINE centre_velocities(u_face, v_face, u_centre, v_centre)
    REAL(wp), DIMENSION(0:,0:,0:), INTENT(IN)  ::  u_face, v_face
    REAL(wp), DIMENSION(0:,0:,0:), INTENT(OUT) ::  u_centre, v_centre
    INTEGER(iwp) ::  nx, ny

    nx = UBOUND(u_face, 1)
    ny = UBOUND(u_face, 2)

    u_centre(0,:,:)  = 0.0_wp
    u_centre(1:,:,:) = 0.5_wp * ( u_face(0:nx-1,:,:) + u_face(1:,:,:) )

    v_centre(:,0,:)  = 0.0_wp
    v_centre(:,1:,:) = 0.5_wp * ( v_face(:,0:ny-1,:) + v_face(:,1:,:) )
 END SUBROUTINE centre_velocities


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the geographical latitude of a point given in rotated-pole cordinates
!------------------------------------------------------------------------------!
 FUNCTION phirot2phi (phirot, rlarot, polphi, polgam)
 
    REAL(wp), INTENT (IN) ::  polphi      !< latitude of the rotated north pole
    REAL(wp), INTENT (IN) ::  phirot      !< latitude in the rotated system
    REAL(wp), INTENT (IN) ::  rlarot      !< longitude in the rotated system
    REAL(wp), INTENT (IN) ::  polgam      !< angle between the north poles of the systems

    REAL(wp)              ::  phirot2phi  !< latitude in the geographical system
    
    REAL(wp)              ::  zsinpol, zcospol, zphis, zrlas, zarg, zgam
 
    zsinpol = SIN(polphi * TO_RADIANS)
    zcospol = COS(polphi * TO_RADIANS)
    zphis   = phirot * TO_RADIANS

    IF (rlarot > 180.0_wp)  THEN
       zrlas = rlarot - 360.0_wp
    ELSE
       zrlas = rlarot
    ENDIF
    zrlas = zrlas * TO_RADIANS
  
    IF (polgam /= 0.0_wp)  THEN
       zgam = polgam * TO_RADIANS
       zarg = zsinpol * SIN (zphis) +                                       &
              zcospol * COS(zphis) * ( COS(zrlas) * COS(zgam) -             &
                                       SIN(zgam)  * SIN(zrlas) )
    ELSE
       zarg = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
    ENDIF
   
    phirot2phi = ASIN (zarg) * TO_DEGREES
 
 END FUNCTION phirot2phi


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the geographical latitude of a point given in rotated-pole cordinates
!------------------------------------------------------------------------------!
 FUNCTION phi2phirot (phi, rla, polphi, pollam)
 
    REAL(wp), INTENT (IN) ::  polphi !< latitude of the rotated north pole
    REAL(wp), INTENT (IN) ::  pollam !< longitude of the rotated north pole
    REAL(wp), INTENT (IN) ::  phi    !< latitude in the geographical system
    REAL(wp), INTENT (IN) ::  rla    !< longitude in the geographical system
    
    REAL(wp) ::  phi2phirot          !< longitude in the rotated system
    
    REAL(wp) ::  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1
    
    zsinpol = SIN(polphi * TO_RADIANS)
    zcospol = COS(polphi * TO_RADIANS)
    zlampol = pollam * TO_RADIANS
    zphi    = phi * TO_RADIANS

    IF (rla > 180.0_wp)  THEN
       zrla1 = rla - 360.0_wp
    ELSE
       zrla1 = rla
    ENDIF
    zrla = zrla1 * TO_RADIANS
    
    zarg1 = SIN(zphi) * zsinpol
    zarg2 = COS(zphi) * zcospol * COS(zrla - zlampol)
    
    phi2phirot = ASIN(zarg1 + zarg2) * TO_DEGREES
 
 END FUNCTION phi2phirot


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the geographical longitude of a point given in rotated-pole cordinates
!------------------------------------------------------------------------------!
 FUNCTION rlarot2rla(phirot, rlarot, polphi, pollam, polgam)
 
    REAL(wp), INTENT (IN) ::  polphi !< latitude of the rotated north pole
    REAL(wp), INTENT (IN) ::  pollam !< longitude of the rotated north pole
    REAL(wp), INTENT (IN) ::  phirot !< latitude in the rotated system
    REAL(wp), INTENT (IN) ::  rlarot !< longitude in the rotated system
    REAL(wp), INTENT (IN) ::  polgam !< angle between the north poles of the systems
    
    REAL(wp) ::  rlarot2rla          !< latitude in the geographical system
    
    REAL(wp) ::  zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam
    
    zsinpol = SIN(TO_RADIANS * polphi)
    zcospol = COS(TO_RADIANS * polphi)
    zlampol = TO_RADIANS * pollam
    zphis   = TO_RADIANS * phirot

    IF (rlarot > 180.0_wp)  THEN
       zrlas = rlarot - 360.0_wp
    ELSE
       zrlas = rlarot
    ENDIF
    zrlas   = TO_RADIANS * zrlas
   
    IF (polgam /= 0.0_wp)  THEN
       zgam  = TO_RADIANS * polgam
       zarg1 = SIN(zlampol) * (zcospol * SIN(zphis) - zsinpol*COS(zphis) *  &
               (COS(zrlas) * COS(zgam) - SIN(zrlas) * SIN(zgam)) ) -        &
               COS(zlampol) * COS(zphis) * ( SIN(zrlas) * COS(zgam) +       &
                                             COS(zrlas) * SIN(zgam) )
    
       zarg2 = COS (zlampol) * (zcospol * SIN(zphis) - zsinpol*COS(zphis) * &
               (COS(zrlas) * COS(zgam) - SIN(zrlas) * SIN(zgam)) ) +        &
               SIN(zlampol) * COS(zphis) * ( SIN(zrlas) * COS(zgam) +       &
                                             COS(zrlas) * SIN(zgam) )
    ELSE
       zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +     &
                                   zcospol *              SIN(zphis)) -     &
                 COS (zlampol) *             SIN(zrlas) * COS(zphis)
       zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +     &
                                   zcospol *              SIN(zphis)) +     &
                 SIN (zlampol) *             SIN(zrlas) * COS(zphis)
    ENDIF
   
    IF (zarg2 == 0.0_wp)  zarg2 = 1.0E-20_wp
   
    rlarot2rla = ATAN2(zarg1,zarg2) * TO_DEGREES
     
 END FUNCTION rlarot2rla


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the rotated-pole longitude of a point given in geographical cordinates
!------------------------------------------------------------------------------!
 FUNCTION rla2rlarot ( phi, rla, polphi, pollam, polgam )

    REAL(wp), INTENT (IN) ::  polphi !< latitude of the rotated north pole
    REAL(wp), INTENT (IN) ::  pollam !< longitude of the rotated north pole
    REAL(wp), INTENT (IN) ::  phi    !< latitude in geographical system
    REAL(wp), INTENT (IN) ::  rla    !< longitude in geographical system
    REAL(wp), INTENT (IN) ::  polgam !< angle between the north poles of the systems
    
    REAL(wp) ::  rla2rlarot    !< latitude in the the rotated system
    
    REAL(wp) ::  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1
    
    zsinpol = SIN(polphi * TO_RADIANS)
    zcospol = COS(polphi * TO_RADIANS)
    zlampol = pollam * TO_RADIANS
    zphi    = phi * TO_RADIANS

    IF (rla > 180.0_wp)  THEN
       zrla1 = rla - 360.0_wp
    ELSE
       zrla1 = rla
    ENDIF
    zrla = zrla1 * TO_RADIANS
    
    zarg1 = - SIN (zrla-zlampol) * COS(zphi)
    zarg2 = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)
    
    IF (zarg2 == 0.0_wp)  zarg2 = 1.0E-20_wp
    
    rla2rlarot = ATAN2 (zarg1,zarg2) * TO_DEGREES
    
    IF (polgam /= 0.0_wp )  THEN
       rla2rlarot = polgam + rla2rlarot
       IF (rla2rlarot > 180._wp)  rla2rlarot = rla2rlarot - 360.0_wp
    ENDIF
    
 END FUNCTION rla2rlarot


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Rotate the given velocity vector (u,v) from the geographical to the
!> rotated-pole system
!------------------------------------------------------------------------------!
 SUBROUTINE uv2uvrot(u, v, rlat, rlon, pollat, pollon, urot, vrot)
 
    REAL(wp), INTENT (IN)  ::  u, v           !< wind components in the true geographical system
    REAL(wp), INTENT (IN)  ::  rlat, rlon     !< coordinates in the true geographical system
    REAL(wp), INTENT (IN)  ::  pollat, pollon !< latitude and longitude of the north pole of the rotated grid
    
    REAL(wp), INTENT (OUT) ::  urot, vrot     !< wind components in the rotated grid             
    
    REAL (wp) ::  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm
    
    zsinpol = SIN(pollat * TO_RADIANS)
    zcospol = COS(pollat * TO_RADIANS)
    zlonp   = (pollon-rlon) * TO_RADIANS
    zlat    = rlat * TO_RADIANS
    
    zarg1 = zcospol * SIN(zlonp)
    zarg2 = zsinpol * COS(zlat) - zcospol * SIN(zlat) * COS(zlonp)
    znorm = 1.0_wp / SQRT(zarg1*zarg1 + zarg2*zarg2)
    
    urot = u * zarg2 * znorm - v * zarg1 * znorm
    vrot = u * zarg1 * znorm + v * zarg2 * znorm
 
 END SUBROUTINE uv2uvrot


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Rotate the given velocity vector (urot, vrot) from the rotated-pole to the
!> geographical system
!------------------------------------------------------------------------------!
 SUBROUTINE uvrot2uv (urot, vrot, rlat, rlon, pollat, pollon, u, v)
 
    REAL(wp), INTENT(IN) ::  urot, vrot     !< wind components in the rotated grid
    REAL(wp), INTENT(IN) ::  rlat, rlon     !< latitude and longitude in the true geographical system
    REAL(wp), INTENT(IN) ::  pollat, pollon !< latitude and longitude of the north pole of the rotated grid
    
    REAL(wp), INTENT(OUT) ::  u, v          !< wind components in the true geographical system
    
    REAL(wp) ::  zsinpol, zcospol, zlonp, zlat, zarg1, zarg2, znorm
  
    zsinpol = SIN(pollat * TO_RADIANS)
    zcospol = COS(pollat * TO_RADIANS)
    zlonp   = (pollon-rlon) * TO_RADIANS
    zlat    = rlat * TO_RADIANS
  
    zarg1 = zcospol * SIN(zlonp)
    zarg2 = zsinpol * COS(zlat) - zcospol * SIN(zlat) * COS(zlonp)
    znorm = 1.0_wp / SQRT(zarg1*zarg1 + zarg2*zarg2)
  
    u =   urot * zarg2 * znorm + vrot * zarg1 * znorm
    v = - urot * zarg1 * znorm + vrot * zarg2 * znorm
 
 END SUBROUTINE uvrot2uv

 END MODULE inifor_transform
#endif

