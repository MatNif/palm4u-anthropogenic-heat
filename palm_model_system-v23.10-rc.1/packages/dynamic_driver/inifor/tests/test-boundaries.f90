!> @file tests/test-boundaries.f90
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
! @author Eckhard Kadasch
!
! Description:
! ------------
!> This program tests the boundary grid mode of INIFOR's init_grid_definition()
!> routine.
!------------------------------------------------------------------------------!
 PROGRAM test_boundaries

    USE inifor_defs,                                                           &
        ONLY :  iwp, wp
    USE inifor_grid,                                                           &
        ONLY :  init_grid_definition
    USE inifor_types,                                                          &
        ONLY :  grid_definition
    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER     ::  title = 'boundary initialization'
    CHARACTER(LEN=20), DIMENSION(2) ::  kind_list = (/ 'boundary', 'boundary' /)
    LOGICAL                         ::  res

    INTEGER(iwp)                    ::  i, nx, ny, nz
    TYPE(grid_definition)           ::  boundary_grid

    REAL(wp) ::  dx, dy, dz, lx, ly, lz, x(2), y(10)
    REAL(wp), TARGET :: z(10)

    CALL begin_test(title, res)

    ! Arange
    dx = 1e-3_wp
    dy = 1.0_wp
    dz = 10.0_wp
    nx = 9_iwp
    ny = 9_iwp
    nz = 9_iwp
    lx = 1.0_wp
    ly = 1e1_wp
    lz = 1e2_wp
    x =   (/ -0.5_wp*dx, lx + 0.5_wp*dx /)
    y = ( (/0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp/) + 0.5_wp )
    z = ( (/0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp/) + 0.5_wp ) * 10

    DO i = 1, SIZE(kind_list)
   
       ! Act
       CALL init_grid_definition(                                              &
          kind = kind_list(i), grid = boundary_grid,                           &
          xmin = x(i), xmax = x(i),                                            &
          ymin =  0.5_wp * dy, ymax = ly - 0.5_wp * dy,                        &
          x0 = 0.0_wp, y0 = 0.0_wp, z0 = 0.0_wp,                               &
          nx = 0_iwp, ny = ny, nz = nz, z = z                                  &
       )
          
   
       ! Assert
       ! asserting that grid % x has exactly two entries and that they match
       ! expected coordinates
       res = res .AND. assert_equal(boundary_grid % x, (/ x(i) /), 'x coordinates')

       ! asserting that grid % y and % z have expected ranges and coordinates
       res = res .AND. assert_equal( boundary_grid % y, y, 'y coordinates')
       res = res .AND. assert_equal( boundary_grid % z, z, 'z coordinates')
   
       CALL fini_grid_definition(boundary_grid)
    ENDDO

    CALL end_test(title, res)

 CONTAINS

 SUBROUTINE fini_grid_definition(grid)
    TYPE(grid_definition), INTENT(INOUT) ::  grid

    DEALLOCATE( grid % x, grid % y )
    DEALLOCATE( grid % kk )
    DEALLOCATE( grid % w_verti )

 END SUBROUTINE fini_grid_definition

 END PROGRAM test_boundaries
