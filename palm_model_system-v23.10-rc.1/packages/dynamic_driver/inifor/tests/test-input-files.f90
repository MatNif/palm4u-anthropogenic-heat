!> @file tests/test-input-files.f90
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
!> This program tests INIFOR's timestamping used for generating input file
!> names.
!------------------------------------------------------------------------------!
 PROGRAM test_input_files

    USE inifor_defs,                                                           &
        ONLY :  PATH, wp, iwp
    USE inifor_io,                                                             &
        ONLY :  get_datetime_file_list
    USE test_utils
    
    IMPLICIT NONE

    CHARACTER(LEN=60)                              ::  title
    CHARACTER(LEN=PATH), ALLOCATABLE, DIMENSION(:) ::  file_list, ref_list
    LOGICAL                                        ::  res
    INTEGER(iwp)                                   ::  i

    title = "input files - daylight saving to standard time"
    CALL begin_test(title, res)

    ! Arange
    ! ...a date range that inlcudes a shift from daylight saving time to
    ! standard time (29.10.2017). Since all time stamps in COSMO-DE input files
    ! are in UTC, this should not the naming cadence.
    ALLOCATE( ref_list(6) )
    ref_list(1)  = './laf2017102823-test.nc'
    ref_list(2)  = './laf2017102900-test.nc'
    ref_list(3)  = './laf2017102901-test.nc'
    ref_list(4)  = './laf2017102902-test.nc'
    ref_list(5)  = './laf2017102903-test.nc'
    ref_list(6)  = './laf2017102904-test.nc'

    ! Act
    CALL get_datetime_file_list(                                               &
       start_date_string='2017102823',                                         &
       start_hour=0_iwp, end_hour=5_iwp, step_hour=1_iwp,                      &
       input_path='./', prefix="laf", suffix='-test',                          &
       file_list=file_list                                                     &
    )

    ! Assert
    DO i = 1, 6
       res = res .AND. (TRIM(ref_list(i)) .EQ. TRIM(file_list(i)))
    ENDDO

    DEALLOCATE( ref_list, file_list )
    CALL end_test(title, res)


    title = "input files - leap day"
    CALL begin_test(title, res)

    ! Arange
    ! ...a date range that inlcudes a leap day (29. Feb. 2016) which should be
    ! inlcuded in UTC time stamps.
    ALLOCATE( ref_list(2) )
    ref_list(1)  = './laf2016022823-test.nc'
    ref_list(2)  = './laf2016022900-test.nc'

    ! Act
    CALL get_datetime_file_list(                                               &
       start_date_string='2016022823',                                         &
       start_hour=0_iwp, end_hour=1_iwp, step_hour=1_iwp,                      &
       input_path='./', prefix="laf", suffix='-test',                          &
       file_list=file_list                                                     &
    )

    ! Assert
    DO i = 1, 2
       res = res .AND. (TRIM(ref_list(i)) .EQ. TRIM(file_list(i)))
    ENDDO

    DEALLOCATE( ref_list, file_list )
    CALL end_test(title, res)



    title = "input files - negative start_hour and step_hour > 1 hour"
    CALL begin_test(title, res)

    ! Arange
    ! ...a date range that inlcudes a leap day (29. Feb. 2016) which should be
    ! inlcuded in UTC time stamps.
    ALLOCATE( ref_list(4) )
    ref_list(1)  = './laf2017102823-test.nc'
    ref_list(2)  = './laf2017102901-test.nc'
    ref_list(3)  = './laf2017102903-test.nc'
    ref_list(4)  = './laf2017102904-test.nc'

    ! Act
    CALL get_datetime_file_list(                                               &
       start_date_string='2017102901',                                         &
       start_hour=-2_iwp, end_hour=3_iwp, step_hour=2_iwp,                     &
       input_path='./', prefix="laf", suffix='-test',                          &
       file_list=file_list                                                     &
    )

    PRINT *, file_list

    ! Assert
    DO i = 1, 2
       res = res .AND. (TRIM(ref_list(i)) .EQ. TRIM(file_list(i)))
    ENDDO

    DEALLOCATE( ref_list, file_list )
    CALL end_test(title, res)

 END PROGRAM test_input_files
