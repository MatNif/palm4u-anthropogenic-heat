!> @file tests/test-stretching.f90
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
!> This program tests INIFOR's implementation of PALM's grid stretching.
!> 
!------------------------------------------------------------------------------!
 PROGRAM test_stretching

    USE inifor_defs,                                                           &
        ONLY :  wp

    USE inifor_grid,                                                           &
        ONLY :  stretched_z

    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER ::  title = "stretched grid"
    LOGICAL                     ::  res

    INTEGER, PARAMETER ::  nz = 9
    INTEGER            ::  k  = nz - 2

    REAL(wp) ::  z(1:nz)
    REAL(wp) ::  dz(10)            = -1.0_wp
    REAL(wp) ::  dz_max            = 1000.0_wp
    REAL(wp) ::  dz_stretch_factor = 1.08_wp
    REAL(wp) ::  dz_stretch_level  = 2.0_wp
    REAL(wp) ::  dz_stretch_level_start(9) = -9999999.9_wp
    REAL(wp) ::  dz_stretch_level_end(9) = 9999999.9_wp
    REAL(wp) ::  dz_stretch_factor_array(9) = 1.08_wp

    CALL begin_test(title, res)

    ! Arange
    z(:)   = 0.0_wp
    dz(1)  = 1.0_wp

    ! Act
    CALL stretched_z(z, dz, dz_max=dz_max, &
                     dz_stretch_factor=dz_stretch_factor,                   &
                     dz_stretch_level=dz_stretch_level,                     &
                     dz_stretch_level_start=dz_stretch_level_start,         &
                     dz_stretch_level_end=dz_stretch_level_end,             &
                     dz_stretch_factor_array=dz_stretch_factor_array)

    ! Assert, that the total distance covered by the stretched region 
    ! matches the therotetial distance, i.e. the sum over the finite
    ! exponential series
    !
    !          Sum_{i=0}^n( dz[i] ) = dz[0] * Sum_{i=0}^n( f^i )
    !                               = dz[0] * (1 - f^(i+1)) / (1-f)
    !
    ! with f being the stretch factor.
    res = res .AND. &
          assert_equal( (/ z(UBOUND(z, 1)) - z(1)              /),             &
                        (/ (1.0_wp - dz_stretch_factor**(k+1)) /               &
                           (1.0_wp - dz_stretch_factor)        /),             &
                        'length of stretched grid' )

    CALL end_test(title, res)

 END PROGRAM test_stretching
