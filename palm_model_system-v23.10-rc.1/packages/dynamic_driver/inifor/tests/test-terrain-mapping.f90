!> @file tests/test-terrain-mapping.f90
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
! Current revisions:
! -----------------
! Initial revision
! 
! 
! Former revisions:
! -----------------
! 
!
! Authors:
! --------
! @author Eckhard Kadasch
!
! Description:
! ------------
!> Tests for INIFOR's mesoscale-microscale terrain mapping
!------------------------------------------------------------------------------!
 PROGRAM test_terrain_mapping

    USE inifor_defs,                                                           &
        ONLY : iwp, wp
    USE inifor_transform,                                                      &
        ONLY : is_monotone
    USE test_utils

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER ::  title = "monotonicity check"
    LOGICAL                     ::  res

    REAL(wp)     ::  z(3,3,10)
    INTEGER(iwp) ::  i,j,k

    CALL begin_test(title, res)

    DO k = 1, 10
    DO j = 1, 3
    DO i = 1, 3
       z(i,j,k) = k
    ENDDO
    ENDDO
    ENDDO

    res = res .AND. is_monotone( z )
    PRINT *, z(1,1,:)
    PRINT *, is_monotone( z )


    z(1,1,5) = z(1,1,6) + 0.1_wp
    res = res .AND. .NOT. is_monotone( z )
    PRINT *, z(1,1,:)
    PRINT *, is_monotone( z )

    CALL end_test(title, res)

 END PROGRAM test_terrain_mapping
