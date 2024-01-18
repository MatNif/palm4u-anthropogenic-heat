!> @file tests/test-grid.f90
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
!> This program tests the PALM grid mode of INIFOR's init_grid_definition()
!> routine.
!------------------------------------------------------------------------------!
 PROGRAM test_grid

    USE inifor_defs,                                                           &
        ONLY : wp, iwp

    USE inifor_grid,                                                           &
        ONLY :  grid_definition, init_grid_definition, dx, dy, dz

    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER ::  title = "grid initialization"
    LOGICAL                     ::  res

    TYPE(grid_definition)       ::  mygrid
    INTEGER(iwp)                ::  i
    INTEGER(iwp), PARAMETER     ::  nx = 9,   ny = 19,   nz = 29
    REAL(wp), PARAMETER         ::  lx = 100., ly = 200., lz = 300.
    REAL(wp), DIMENSION(0:nx)   ::  x, xu
    REAL(wp), DIMENSION(0:ny)   ::  y, yv
    REAL(wp), DIMENSION(1:nz)   ::  z
    REAL(wp), DIMENSION(1:nz-1) ::  zw

    CALL begin_test(title, res)

    ! Arange
    dx = lx / (nx + 1)
    DO i = 0, nx
       xu(i) = REAL(i, KIND=wp) / (nx+1) * lx
       x(i)  = 0.5*dx + xu(i)
    ENDDO

    dy = ly / (ny + 1)
    DO i = 0, ny
       yv(i) = REAL(i, KIND=wp) / (ny+1) * ly
       y(i)  = 0.5*dy + yv(i)
    ENDDO

    dz(:) = lz / (nz + 1)
    DO i = 1, nz
       z(i) = REAL(i, KIND=wp) / (nz+1) * lz - 0.5_wp * dz(1)
       IF (i < nz) zw(i) = REAL(i, KIND=wp) / (nz+1) * lz
    ENDDO

    ! Act
    CALL init_grid_definition('palm', grid = mygrid,                           &
                              xmin = 0.0_wp, xmax = lx,                        &
                              ymin = 0.0_wp, ymax = ly,                        &
                              x0 = 0.0_wp, y0 = 0.0_wp, z0 = 0.0_wp,           &
                              nx = nx, ny = ny, nz = nz,                       &
                              z = z, zw = zw)

    ! Assert coordinates match
    res = res .AND. assert_equal(x,      mygrid%x,  "x" )
    res = res .AND. assert_equal(xu(1:), mygrid%xu, "xu")
    res = res .AND. assert_equal(y,      mygrid%y,  "y" )
    res = res .AND. assert_equal(yv(1:), mygrid%yv, "yv")
    res = res .AND. assert_equal(z,      mygrid%z,  "z" )
    res = res .AND. assert_equal(zw(1:), mygrid%zw, "zw")

    CALL end_test(title, res)

 END PROGRAM test_grid
