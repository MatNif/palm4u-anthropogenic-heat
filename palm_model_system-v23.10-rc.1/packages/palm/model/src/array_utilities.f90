!> @file array_utilities.f90
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
! Copyright 2021-2022 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
!> @author Helge Knoop (design and initial implementation)
!
! Description:
! ------------
!> A collection of handy array utilities.
!--------------------------------------------------------------------------------------------------!
 MODULE array_utilities

    IMPLICIT NONE

    PUBLIC                                                                                         &
       quicksort

    INTERFACE quicksort
       MODULE PROCEDURE quicksort_int8
       MODULE PROCEDURE quicksort_int16
       MODULE PROCEDURE quicksort_int32
       MODULE PROCEDURE quicksort_int32_2dim_e1
       MODULE PROCEDURE quicksort_real32
       MODULE PROCEDURE quicksort_real64
    END INTERFACE quicksort

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform the quicksort algorithm on a 1D int8 input array.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_int8( a )

    INTEGER(KIND=1), DIMENSION(:) ::  a  !< input array

    INTEGER(KIND=1) ::  t  !< variable for temporary value cache
    INTEGER(KIND=1) ::  x  !< variable for temporary value cache

    INTEGER ::  first = 1  !< index variable to mark a sorting pivot point
    INTEGER ::  last       !< index variable to mark a sorting pivot point
    INTEGER ::  i          !< index variable
    INTEGER ::  j          !< index variable


    last = SIZE( a, 1 )
    x = a((first+last)/2)
    i = first
    j = last

    DO
       DO WHILE ( a(i) < x )
          i = i + 1
       ENDDO
       DO WHILE ( x < a(j) )
          j = j - 1
       ENDDO
       IF ( i >= j )  EXIT
       t    = a(i)
       a(i) = a(j)
       a(j) = t
       i = i + 1
       j = j - 1
    ENDDO

    IF ( first < i - 1 )  CALL quicksort_int8( a(first:i-1) )
    IF ( j + 1 < last  )  CALL quicksort_int8( a(j+1:last) )

 END SUBROUTINE quicksort_int8


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform the quicksort algorithm on a 1D int16 input array.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_int16( a )

    INTEGER(KIND=2), DIMENSION(:) ::  a  !< input array

    INTEGER(KIND=2) ::  t  !< variable for temporary value cache
    INTEGER(KIND=2) ::  x  !< variable for temporary value cache

    INTEGER ::  first = 1  !< index variable to mark a sorting pivot point
    INTEGER ::  last       !< index variable to mark a sorting pivot point
    INTEGER ::  i          !< index variable
    INTEGER ::  j          !< index variable


    last = SIZE( a, 1 )
    x = a((first+last)/2)
    i = first
    j = last

    DO
       DO WHILE ( a(i) < x )
          i = i + 1
       ENDDO
       DO WHILE ( x < a(j) )
          j = j - 1
       ENDDO
       IF ( i >= j )  EXIT
       t    = a(i)
       a(i) = a(j)
       a(j) = t
       i = i + 1
       j = j - 1
    ENDDO

    IF ( first < i - 1 )  CALL quicksort_int16( a(first:i-1) )
    IF ( j + 1 < last  )  CALL quicksort_int16( a(j+1:last) )

 END SUBROUTINE quicksort_int16


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform the quicksort algorithm on a 1D int32 input array.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_int32( a )

    INTEGER(KIND=4), DIMENSION(:) ::  a  !< input array

    INTEGER(KIND=4) ::  t  !< variable for temporary value cache
    INTEGER(KIND=4) ::  x  !< variable for temporary value cache

    INTEGER ::  first = 1  !< index variable to mark a sorting pivot point
    INTEGER ::  last  !< index variable to mark a sorting pivot point
    INTEGER ::  i  !< index variable
    INTEGER ::  j  !< index variable


    last = SIZE( a, 1 )
    x = a((first+last)/2)
    i = first
    j = last

    DO
       DO WHILE ( a(i) < x )
          i = i + 1
       ENDDO
       DO WHILE ( x < a(j) )
          j = j - 1
       ENDDO
       IF ( i >= j ) EXIT
       t    = a(i)
       a(i) = a(j)
       a(j) = t
       i = i + 1
       j = j - 1
    ENDDO

    IF ( first < i - 1 )  CALL quicksort_int32( a(first:i-1) )
    IF ( j + 1 < last  )  CALL quicksort_int32( a(j+1:last) )

 END SUBROUTINE quicksort_int32


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform the quicksort algorithm on a 2D int32 input array.
!> The values a(1,...) are base of the sort. All elements of dim1 are sorted relative to a(1,...)
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_int32_2dim_e1( a )

    INTEGER(KIND=4), DIMENSION(:,:) ::  a  !< input array

    INTEGER(KIND=4), DIMENSION(size(a,1)) ::  t  !< variable for temporary value cache
    INTEGER(KIND=4), DIMENSION(size(a,1)) ::  x  !< variable for temporary value cache

    INTEGER ::  first = 1  !< index variable to mark a sorting pivot point
    INTEGER ::  last       !< index variable to mark a sorting pivot point
    INTEGER ::  i          !< index variable
    INTEGER ::  j          !< index variable


    IF ( SIZE( a, 2 ) == 0 )  RETURN

    last = SIZE( a, 2 )
    x = a(:,(first+last)/2)
    i = first
    j = last

    DO
       DO WHILE ( a(1,i) < x(1) )
          i = i + 1
       ENDDO
       DO WHILE ( x(1) < a(1,j) )
          j = j - 1
       ENDDO
       IF ( i >= j ) EXIT
       t    = a(:,i)
       a(:,i) = a(:,j)
       a(:,j) = t
       i = i + 1
       j = j - 1
    ENDDO

    IF ( first < i - 1 )  CALL quicksort( a(:,first:i-1) )
    IF ( j + 1 < last  )  CALL quicksort( a(:,j+1:last) )

 END SUBROUTINE quicksort_int32_2dim_e1


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform the quicksort algorithm on a 1D real32 input array.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_real32( a )

    REAL(KIND=4), DIMENSION(:) ::  a  !< input array

    INTEGER ::  first = 1  !< index variable to mark a sorting pivot point
    INTEGER ::  last  !< index variable to mark a sorting pivot point
    INTEGER ::  i  !< index variable
    INTEGER ::  j  !< index variable

    REAL(KIND=4) ::  t  !< variable for temporary value cache
    REAL(KIND=4) ::  x  !< variable for temporary value cache


    last = SIZE( a, 1 )
    x = a((first+last)/2)
    i = first
    j = last

    DO
       DO WHILE ( a(i) < x )
          i = i + 1
       ENDDO
       DO WHILE ( x < a(j) )
          j = j - 1
       ENDDO
       IF ( i >= j ) EXIT
       t    = a(i)
       a(i) = a(j)
       a(j) = t
       i = i + 1
       j = j - 1
    ENDDO

    IF ( first < i - 1 )  CALL quicksort_real32( a(first:i-1) )
    IF ( j + 1 < last  )  CALL quicksort_real32( a(j+1:last) )

 END SUBROUTINE quicksort_real32


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform the quicksort algorithm on a 1D real64 input array.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_real64( a )

    REAL(KIND=8), DIMENSION(:) ::  a  !< input array

    INTEGER ::  first = 1  !< index variable to mark a sorting pivot point
    INTEGER ::  last  !< index variable to mark a sorting pivot point
    INTEGER ::  i  !< index variable
    INTEGER ::  j  !< index variable

    REAL(KIND=8) ::  t  !< variable for temporary value cache
    REAL(KIND=8) ::  x  !< variable for temporary value cache


    last = SIZE( a, 1 )
    x = a((first+last)/2)
    i = first
    j = last

    DO
       DO WHILE ( a(i) < x )
          i = i + 1
       ENDDO
       DO WHILE ( x < a(j) )
          j = j - 1
       ENDDO
       IF ( i >= j ) EXIT
       t = a(i);  a(i) = a(j);  a(j) = t
       i = i + 1
       j = j - 1
    ENDDO

    IF ( first < i - 1 )  CALL quicksort_real64( a(first:i-1) )
    IF ( j + 1 < last  )  CALL quicksort_real64( a(j+1:last) )

 END SUBROUTINE quicksort_real64

 END MODULE array_utilities
