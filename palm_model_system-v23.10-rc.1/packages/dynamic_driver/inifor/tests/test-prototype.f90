!> @file tests/test-prototype.f90
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
!> This prototype test does nothing, but serves as a template for programming
!> new INIFOR tests.
!------------------------------------------------------------------------------!
 PROGRAM test_prototype

    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER ::  title = "Prototype"
    LOGICAL                     ::  res

    CALL begin_test(title, res)

    ! Arange
    !define parameters and reference values

    ! Act
    !compute result

    ! Assert
    !res = res .AND. assert_equal(<result_array>, <reference_array>, 'description')

    CALL end_test(title, res)

 END PROGRAM test_prototype
