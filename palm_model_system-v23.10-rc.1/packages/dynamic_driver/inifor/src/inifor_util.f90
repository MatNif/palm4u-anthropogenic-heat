!> @file src/inifor_util.f90
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
!> The util module provides miscellaneous utility routines for INIFOR.
!------------------------------------------------------------------------------!
 MODULE inifor_util

    USE inifor_defs,                                                           &
        ONLY :  PI, DATE, SNAME, iwp, wp
    USE inifor_types,                                                          &
        ONLY :  grid_definition
    USE, INTRINSIC :: ISO_C_BINDING,                                           &
        ONLY :  C_CHAR, C_INT, C_PTR, C_SIZE_T

    IMPLICIT NONE

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fortran implementation of C's struct tm for representing points in time
!------------------------------------------------------------------------------!
    TYPE, BIND(c) :: tm_struct
       INTEGER(C_INT) :: tm_sec     !< seconds after the minute [0, 61]
       INTEGER(C_INT) :: tm_min     !< minutes after the hour [0, 59]
       INTEGER(C_INT) :: tm_hour    !< hours since midnight [0, 23]
       INTEGER(C_INT) :: tm_mday    !< day of the month [1, 31]
       INTEGER(C_INT) :: tm_mon     !< month since January [0, 11]
       INTEGER(C_INT) :: tm_year    !< years since 1900
       INTEGER(C_INT) :: tm_wday    !< days since Sunday [0, 6]
       INTEGER(C_INT) :: tm_yday    !< days since January 1st [0, 356]
       INTEGER(C_INT) :: tm_isdst   !< Daylight Saving Time flag
    END TYPE

 INTERFACE

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interface to C's strptime function, which converts a string to a tm
!> structure.
!------------------------------------------------------------------------------!
    FUNCTION strptime(string, format, timeinfo) BIND(c, NAME='strptime')
       IMPORT :: C_CHAR, C_SIZE_T, tm_struct

       IMPLICIT NONE

       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) ::  string, format
       TYPE(tm_struct), INTENT(OUT)                     ::  timeinfo

       INTEGER(C_SIZE_T)                                ::  strptime
    END FUNCTION


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interface to C's strftime function, which converts the given 'timeinfo'
!> structure to a string in the given 'format'.
!------------------------------------------------------------------------------!
    FUNCTION strftime(string, string_len, format, timeinfo) BIND(c, NAME='strftime')
       IMPORT :: C_CHAR, C_SIZE_T, tm_struct

       IMPLICIT NONE

       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) ::  string
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN)  ::  format
       INTEGER(C_SIZE_T), INTENT(IN)                     ::  string_len
       TYPE(tm_struct), INTENT(IN)                       ::  timeinfo

       INTEGER(C_SIZE_T)                                 ::  strftime
    END FUNCTION


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interface to C's mktime function, which converts the given tm structure to
!> Unix time (number of seconds since 0 UTC, 1st January 1970). INIFOR uses the
!> side effect that a call to mktime noramlizes the given 'timeinfo' structure,
!> e.g. increments the date if hours overfow 24.
!------------------------------------------------------------------------------!
    FUNCTION mktime(timeinfo) BIND(c, NAME='mktime')
       IMPORT :: C_PTR, tm_struct

       IMPLICIT NONE

       TYPE(tm_struct), INTENT(IN) ::  timeinfo

       TYPE(C_PTR)                 ::  mktime
    END FUNCTION

 END INTERFACE

 
 INTERFACE str
    MODULE PROCEDURE str_iwp
    MODULE PROCEDURE str_4
 END INTERFACE


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Takes a string of the form YYYYMMDDHH, adds the given number of hours to its
!> tm structure representation, and returns the corresponding string in the same
!> format
!------------------------------------------------------------------------------!
 CHARACTER(LEN=DATE) FUNCTION add_hours_to(date_string, hours)
    CHARACTER(LEN=DATE), INTENT(IN)          ::  date_string
    INTEGER(iwp), INTENT(IN)                 ::  hours

    CHARACTER(KIND=C_CHAR, LEN=*), PARAMETER ::  format_string = "%Y%m%d%H"
    CHARACTER(KIND=C_CHAR, LEN=DATE)         ::  c_date_string
    TYPE(C_PTR)                              ::  c_pointer
    TYPE(tm_struct)                          ::  time_info
    INTEGER(iwp)                             ::  err

    c_date_string = date_string

    ! Convert C string to C tm struct
    CALL init_tm(time_info)
    err = strptime(c_date_string, format_string, time_info)
 
    ! Manipulate and normalize C tm struct
    time_info%tm_hour = time_info%tm_hour + hours
    c_pointer = mktime(time_info)

    ! Convert back to C string
    err = strftime(c_date_string, INT(DATE, KIND=C_SIZE_T),                 &
                   format_string, time_info)

    add_hours_to = c_date_string
 END FUNCTION


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Print all members of the given tm structure
!------------------------------------------------------------------------------!
 SUBROUTINE print_tm(timeinfo)
    TYPE(tm_struct), INTENT(IN) :: timeinfo

    PRINT *, "sec: ", timeinfo%tm_sec,  &  !< seconds after the minute [0, 61]
             "min: ", timeinfo%tm_min,  &  !< minutes after the hour [0, 59]
             "hr:  ", timeinfo%tm_hour, &  !< hours since midnight [0, 23]
             "day: ", timeinfo%tm_mday, &  !< day of the month [1, 31]
             "mon: ", timeinfo%tm_mon,  &  !< month since January [0, 11]
             "yr:  ", timeinfo%tm_year, &  !< years since 1900
             "wday:", timeinfo%tm_wday, &  !< days since Sunday [0, 6]
             "yday:", timeinfo%tm_yday, &  !< days since January 1st [0, 356]
             "dst: ", timeinfo%tm_isdst    !< Daylight Saving time flag
 END SUBROUTINE print_tm

    
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize the given tm structure with zero values
!------------------------------------------------------------------------------!
 SUBROUTINE init_tm(timeinfo)
    TYPE(tm_struct), INTENT(INOUT) :: timeinfo

    timeinfo%tm_sec   = 0
    timeinfo%tm_min   = 0
    timeinfo%tm_hour  = 0
    timeinfo%tm_mday  = 0
    timeinfo%tm_mon   = 0
    timeinfo%tm_year  = 0
    timeinfo%tm_wday  = 0
    timeinfo%tm_yday  = 0

    ! We use UTC times, so marking Daylight Saving Time (DST) 'not available'
    ! (< 0). If this is set to 0, mktime will convert the timeinfo to DST and
    ! add one hour.
    timeinfo%tm_isdst = -1
 END SUBROUTINE init_tm


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fill the given array with values equally spaced between and including start
!> and stop
!------------------------------------------------------------------------------!
 SUBROUTINE linspace(start, stop, array)

    REAL(wp), INTENT(IN)    ::  start, stop
    REAL(wp), INTENT(INOUT) ::  array(0:)
    INTEGER(iwp)            ::  i, n

    n = UBOUND(array, 1)

    IF (n .EQ. 0)  THEN

       array(0) = start

    ELSE

       DO i = 0, n
          array(i) = start + REAL(i, wp) / n * (stop - start)
       ENDDO

    ENDIF
    
 END SUBROUTINE linspace


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reverse the order of the third (vertical) array dimension from top-down
!> (COSMO) to bottom-up (PALM)
!------------------------------------------------------------------------------!
 SUBROUTINE reverse(input_arr)

    INTEGER(iwp) ::  idx, opposite_idx, half_idx
    INTEGER(iwp) ::  size_1st_dimension
    INTEGER(iwp) ::  size_2nd_dimension
    INTEGER(iwp) ::  size_3rd_dimension
    INTEGER(iwp) ::  lbound_3rd_dimension
    INTEGER(iwp) ::  ubound_3rd_dimension

    REAL(wp), INTENT(INOUT) ::  input_arr(:,:,:)
    REAL(wp), ALLOCATABLE  :: buffer_arr(:,:)

    lbound_3rd_dimension = LBOUND(input_arr, 3)
    ubound_3rd_dimension = UBOUND(input_arr, 3)
    size_1st_dimension = SIZE(input_arr, 1)
    size_2nd_dimension = SIZE(input_arr, 2)
    size_3rd_dimension = SIZE(input_arr, 3)
    half_idx = lbound_3rd_dimension + size_3rd_dimension / 2 - 1

    ALLOCATE( buffer_arr(size_1st_dimension, size_2nd_dimension) )

    DO  idx = lbound_3rd_dimension, half_idx
       opposite_idx = ubound_3rd_dimension - ( idx - lbound_3rd_dimension )
       buffer_arr(:,:) = input_arr(:,:,idx)
       input_arr(:,:,idx) = input_arr(:,:,opposite_idx)
       input_arr(:,:,opposite_idx) = buffer_arr(:,:)
    ENDDO
    
    DEALLOCATE( buffer_arr )

 END SUBROUTINE reverse


!------------------------------------------------------------------------------!
! Description:
! ------------
!> 
!------------------------------------------------------------------------------!
 SUBROUTINE deaverage(avg_1, t1, avg_2, t2, avg_3, t3)

    REAL(wp), DIMENSION(:,:,:), INTENT(IN)  ::  avg_1, avg_2
    REAL(wp), INTENT(IN)                    ::  t1, t2, t3
    REAL(wp), DIMENSION(:,:,:), INTENT(OUT) ::  avg_3

    REAL(wp)                                ::  ti
 
    ti = 1.0_wp / t3

    avg_3(:,:,:) = ti * ( t2 * avg_2(:,:,:) - t1 * avg_1(:,:,:) )

 END SUBROUTINE deaverage


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the COSMO-DE/-D2 basic state pressure profile
!------------------------------------------------------------------------------!
 SUBROUTINE get_basic_state(z, beta, p_sl, t_sl, rd, g, p0)

    REAL(wp), INTENT(IN)  ::  z(1:)  !< height [m]
    REAL(wp), INTENT(IN)  ::  beta   !< logarithmic lapse rate, dT / d ln(p) [K]
    REAL(wp), INTENT(IN)  ::  p_sl   !< reference pressure [Pa]
    REAL(wp), INTENT(IN)  ::  t_sl   !< reference tempereature [K]
    REAL(wp), INTENT(IN)  ::  rd     !< ideal gas constant of dry air [J/kg/K]
    REAL(wp), INTENT(IN)  ::  g      !< acceleration of Earth's gravity [m/s^2]
    REAL(wp), INTENT(OUT) ::  p0(1:) !< COSMO basic state pressure [Pa]
    REAL(wp) ::  root_frac, factor   !< precomputed factors

    factor = - t_sl / beta
    root_frac = (2.0_wp * beta * g) / (rd * t_sl*t_sl)

    p0(:) = p_sl * EXP(                                                     &
               factor * ( 1.0_wp - SQRT( 1.0_wp - root_frac * z(:) ) )      &
    )

 END SUBROUTINE get_basic_state


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Converts the absolute temperature to the potential temperature in place using
!> the identity a^b = e^(b ln(a)).
!>
!>     theta = T * (p_ref/p)^(R/c_p) = T * e^( R/c_p * ln(p_ref/p) )
!------------------------------------------------------------------------------!
 SUBROUTINE potential_temperature(t, p, p_ref, r, cp)
    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::  t
    REAL(wp), DIMENSION(:,:,:), INTENT(IN)    ::  p
    REAL(wp), INTENT(IN)                      ::  p_ref, r, cp
    REAL(wp)                                  ::  rcp

    rcp = r/cp
    t(:,:,:) =  t(:,:,:) * EXP( rcp * LOG(p_ref / p(:,:,:)) )

 END SUBROUTINE potential_temperature


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the density in place of the given temperature (t_rho).
!------------------------------------------------------------------------------!
 SUBROUTINE moist_density(t_rho, p, qv, rd, rv)
    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) ::  t_rho
    REAL(wp), DIMENSION(:,:,:), INTENT(IN)    ::  p, qv
    REAL(wp), INTENT(IN)                      ::  rd, rv

    t_rho(:,:,:) = p(:,:,:) / (                                             &
       (rv * qv(:,:,:) + rd * (1.0_wp - qv(:,:,:))) * t_rho(:,:,:)          &
    )

 END SUBROUTINE moist_density

!------------------------------------------------------------------------------!
! Description:
! ------------
! Convert a real number to a string in scientific notation showing four 
! significant digits.
!------------------------------------------------------------------------------!
 CHARACTER(LEN=SNAME) FUNCTION real_to_str(val, format)

    REAL(wp), INTENT(IN)                   ::  val
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) ::  format

    IF (PRESENT( format ) )  THEN
       WRITE( real_to_str, format ) val
    ELSE
       WRITE( real_to_str, '(E11.4)' ) val
    ENDIF
    real_to_str = ADJUSTL( real_to_str )

 END FUNCTION real_to_str


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Converts the given real value to a string
!------------------------------------------------------------------------------!
 CHARACTER(LEN=16) FUNCTION real_to_str_f(val)

     REAL(wp), INTENT(IN) ::  val

     WRITE(real_to_str_f, '(F16.8)') val
     real_to_str_f = ADJUSTL(real_to_str_f)

 END FUNCTION real_to_str_f


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Converts the given integer value to a string
!------------------------------------------------------------------------------!
 CHARACTER(LEN=10) FUNCTION str_iwp(val)

     INTEGER(iwp), INTENT(IN) ::  val

     WRITE(str_iwp, '(i10)') val
     str_iwp= ADJUSTL(str_iwp)

 END FUNCTION str_iwp


 CHARACTER(LEN=10) FUNCTION str_4(val)

     INTEGER(4), INTENT(IN) ::  val

     WRITE(str_4, '(i10)') val
     str_4 = ADJUSTL(str_4)

 END FUNCTION str_4
!------------------------------------------------------------------------------!
! Description:
! ------------
!> If the given path is not conlcuded by a slash, add one.
!------------------------------------------------------------------------------!
 SUBROUTINE normalize_path(path)
     
     CHARACTER(LEN=*), INTENT(INOUT) ::  path
     INTEGER(iwp) ::  n

     n = LEN_TRIM(path)

     IF (path(n:n) .NE. '/')  THEN
        path = TRIM(path) // '/'
     ENDIF

 END SUBROUTINE


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check for exact or near equality of two floating point numbers. Inspired by
!> https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
!------------------------------------------------------------------------------!
 LOGICAL FUNCTION nearly_equal(a, b, max_rel_diff)

    REAL(wp), INTENT(IN) ::  a, b, max_rel_diff
    REAL(wp)             ::  diff, mag
 
    diff = ABS( a - b ) 
    mag = MAX( ABS(a), ABS(b) )
    nearly_equal = ( diff .LE. mag * max_rel_diff )

 END FUNCTION nearly_equal

 END MODULE inifor_util
