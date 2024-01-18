!> @file tests/util.f90
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
!> This module provides utiliy functions used in all test programs.
!------------------------------------------------------------------------------!
 MODULE test_utils
 
    USE inifor_defs,                                                           &
        ONLY :  iwp, wp

    IMPLICIT NONE

 CONTAINS

    SUBROUTINE begin_test(title, res)
       CHARACTER(LEN=*), INTENT(IN) ::  title
       LOGICAL, INTENT(OUT)         ::  res

       res = .TRUE.

       PRINT '(/A)',       "******************************************************"
       PRINT '(A, A, A)',  " [  ]  Test '", TRIM(title), "' started."
       PRINT '(A/)',       "******************************************************"
    END SUBROUTINE begin_test

    SUBROUTINE end_test(title, res)
       CHARACTER(LEN=*), INTENT(IN) ::  title
       CHARACTER(LEN=30)            ::  msg, label
       LOGICAL, INTENT(IN)          ::  res


       IF (res .EQV. .TRUE.)  THEN
          msg = 'completed successfully.'
          label = ' [OK]'
       ELSE
          msg = 'failed.'
          label = ' [XX]'
       ENDIF

       PRINT '(/A, A, A, A)', TRIM(label) // "  Test '", TRIM(title), "' ", TRIM(msg)
       PRINT '(A/)',       "******************************************************"

    END SUBROUTINE end_test

    LOGICAL FUNCTION assert_equal(a, b, msg, ratio)
       REAL(wp), OPTIONAL, INTENT(IN)     ::  ratio
       REAL(wp), DIMENSION(:), INTENT(IN) ::  a, b
       CHARACTER(LEN=*), INTENT(IN)   ::  msg

       IF ( PRESENT(ratio) )  THEN
           assert_equal = assert(a, b, 'eq', ratio)
       ELSE
           assert_equal = assert(a, b, 'eq')
       ENDIF

       IF (assert_equal .eqv. .TRUE.)  THEN
           PRINT *, "Equality assertion for ", msg, " was successful."
       ELSE 
           PRINT *, "Equality assertion for ", msg, " failed. Maximum error is ",               &
              MAXVAL( ABS( a - b))
       ENDIF

    END FUNCTION assert_equal

    LOGICAL FUNCTION assert(a, b, mode, ratio)

       REAL(wp), DIMENSION(:), INTENT(IN) ::  a, b
       REAL(wp), OPTIONAL, INTENT(IN)     ::  ratio
       CHARACTER(LEN=*), INTENT(IN)   ::  mode

       REAL(wp)     ::  diff, mag, max_rel_diff
       INTEGER(iwp) ::  i

       max_rel_diff = 10 * EPSILON(1.0)
       IF (PRESENT(ratio)) max_rel_diff = ratio

       SELECT CASE( TRIM(mode) )

       ! This case is inspired by
       ! https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
       CASE('eq')
          IF ( ALL(a(:) == b(:)) )  THEN
             PRINT *, "Checking for exact equality"
             assert = .TRUE.
          ELSE
             assert = .TRUE.
             PRINT *, "Checking for near equality"
             DO i = 1, SIZE(a)
                diff   = ABS(a(i) - b(i)) 
                mag    = MAX( ABS(a(i)), ABS(b(i)) )
                assert = assert .AND. (diff < mag * max_rel_diff )
             ENDDO
          ENDIF 

       CASE DEFAULT
          PRINT *, " Error: Assert mode ", mode, " not implemented. Stopping."
          STOP

       END SELECT
    
    END FUNCTION assert

 END MODULE test_utils

