!> @file posix_interface_mod.f90
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
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Interface to some POSIX system calls, mainly used for read/write of restart files in non-parallel
!> mode in MPI-IO compatible format.
!--------------------------------------------------------------------------------------------------!
 MODULE posix_interface

    USE ISO_C_BINDING

    USE kinds

    IMPLICIT NONE

    PRIVATE

    SAVE

!
!-- Definitions copied from C include file fcntl.h
    INTEGER, PARAMETER ::  o_creat  = 64  !< 0100 octal
    INTEGER, PARAMETER ::  o_rdonly =  0  !<
    INTEGER, PARAMETER ::  o_rdwr   =  2  !<
    INTEGER, PARAMETER ::  o_wronly =  1  !<
    INTEGER, PARAMETER ::  seek_set =  0  !<


!
!-- Interfaces for POSIX calls
    INTERFACE
       INTEGER(C_INT)  FUNCTION C_OPEN( pathname, flags, mode )  BIND( C, NAME = 'open' )
          USE ISO_C_BINDING
          IMPLICIT NONE
          CHARACTER(KIND=C_CHAR), DIMENSION(128) ::  pathname  !<
          INTEGER(KIND=C_INT), VALUE             ::  flags     !<
          INTEGER(KIND=C_INT), VALUE             ::  mode      !<
       END FUNCTION C_OPEN
    END INTERFACE

   INTERFACE
      INTEGER(C_SIZE_T)  FUNCTION C_LSEEK( fd, offset, whence )  BIND( C, NAME = 'lseek' )
         USE ISO_C_BINDING
         IMPLICIT NONE
         INTEGER(KIND=C_INT), VALUE    ::  fd      !<
         INTEGER(KIND=C_SIZE_T), VALUE ::  offset  !<
         INTEGER(KIND=C_INT), VALUE    ::  whence  !<
      END FUNCTION C_LSEEK
   END INTERFACE

!
!-- The read system call uses values of type off_t. There is no Fortran C_OFF_T, therefore  C_SIZE_T
!-- has been used here, assuming both are 8 byte integers.
    INTERFACE
       INTEGER(C_SIZE_T)  FUNCTION C_READ( fd, buf, nr_byte )  BIND(C, NAME = 'read' )
          USE ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE    ::  fd       !<
          INTEGER(KIND=C_SIZE_T), VALUE ::  nr_byte  !<
          TYPE(C_PTR), VALUE            ::  buf      !<
          END FUNCTION C_READ
    END INTERFACE

    INTERFACE
       INTEGER(C_SIZE_T)  FUNCTION C_WRITE( fd, buf, nr_byte )  BIND( C, NAME = 'write' )
          USE ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE    ::  fd       !<
          INTEGER(KIND=C_SIZE_T), VALUE ::  nr_byte  !<
          TYPE(C_PTR), VALUE            ::  buf      !<
       END FUNCTION C_WRITE
    END INTERFACE

    INTERFACE
       INTEGER(C_INT)  FUNCTION C_CLOSE( fd )  BIND( C, NAME = 'close' )
          USE ISO_C_BINDING
          IMPLICIT NONE
          INTEGER(KIND=C_INT), VALUE ::  fd  !<
       END FUNCTION C_CLOSE
    END INTERFACE

!
!-- Sleep function from C library
    INTERFACE
       FUNCTION fsleep( seconds )  BIND( C, NAME = 'sleep' )
          IMPORT
          INTEGER(C_INT)                    ::  fsleep   !<
          INTEGER(C_INT), INTENT(IN), VALUE ::  seconds  !<
       END FUNCTION fsleep
    END INTERFACE

!
!-- PALM interfaces
    INTERFACE fortran_sleep
       MODULE PROCEDURE fortran_sleep
    END INTERFACE fortran_sleep

    INTERFACE posix_close
       MODULE PROCEDURE posix_close
    END INTERFACE posix_close

    INTERFACE posix_lseek
       MODULE PROCEDURE posix_lseek
    END INTERFACE posix_lseek

    INTERFACE posix_open
       MODULE PROCEDURE posix_open
    END INTERFACE posix_open

    INTERFACE posix_read
       MODULE PROCEDURE posix_read
       MODULE PROCEDURE posix_read_char_array
       MODULE PROCEDURE posix_read_int_1d
       MODULE PROCEDURE posix_read_int_2d_i4
       MODULE PROCEDURE posix_read_int_2d_i8
       MODULE PROCEDURE posix_read_i4_3d
       MODULE PROCEDURE posix_read_i8_3d
       MODULE PROCEDURE posix_read_offset_1d
       MODULE PROCEDURE posix_read_real_1d
       MODULE PROCEDURE posix_read_real_2d
       MODULE PROCEDURE posix_read_real_3d
    END INTERFACE posix_read

    INTERFACE posix_write
       MODULE PROCEDURE posix_write
       MODULE PROCEDURE posix_write_char_array
       MODULE PROCEDURE posix_write_int_1d
       MODULE PROCEDURE posix_write_int_2d_i4
       MODULE PROCEDURE posix_write_int_2d_i8
       MODULE PROCEDURE posix_write_i4_3d
       MODULE PROCEDURE posix_write_i8_3d
       MODULE PROCEDURE posix_write_offset_1d
       MODULE PROCEDURE posix_write_real_1d
       MODULE PROCEDURE posix_write_real_2d
       MODULE PROCEDURE posix_write_real_3d
    END INTERFACE posix_write

    PUBLIC fortran_sleep,                                                                          &
           posix_close,                                                                            &
           posix_lseek,                                                                            &
           posix_open,                                                                             &
           posix_read,                                                                             &
           posix_write

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Wait a specified amount of seconds
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fortran_sleep( seconds )

    INTEGER, INTENT(IN) ::  seconds             !< seconds to wait

    INTEGER(C_INT)      ::  seconds_in_c        !< same as seconds
    INTEGER(C_INT)      ::  sleep_return_value  !< returned value to sleep

    seconds_in_c = seconds

    sleep_return_value = fsleep( seconds_in_c )

 END SUBROUTINE fortran_sleep


 INTEGER FUNCTION posix_open( file_name, rd_flag )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)            ::  file_name     !<
    CHARACTER(LEN=1), DIMENSION(:), POINTER ::  f_string      !<
    CHARACTER(LEN=128), TARGET              ::  lo_file_name  !<

    INTEGER(C_INT)                          ::  fd            !<
    INTEGER(C_INT)                          ::  flags         !<
    INTEGER(C_INT)                          ::  name_len      !<
    INTEGER(C_INT)                          ::  mode          !<
    INTEGER, DIMENSION(1)                   ::  bufshape      !<

    LOGICAL, INTENT(IN)                     ::  rd_flag       !<

    TYPE(C_PTR)                             ::  ptr           !<


!
!-- Note: There might be better solutions to convert FORTRAN string to C string but this works on
!-- different FORTRAN compiler
    name_len     = LEN( TRIM( file_name ) ) + 1
    lo_file_name = TRIM( file_name ) // CHAR( 0 )
    ptr          = C_LOC( lo_file_name(1:1) )
    bufshape(1)  = name_len
    CALL C_F_POINTER( ptr, f_string, bufshape )

    mode = 420  ! Mode 644

    IF ( rd_flag )  THEN
       flags = o_rdonly
       fd    = C_OPEN( f_string, flags, mode )  ! Open for reading
    ELSE
       flags = o_wronly + o_creat
       fd    = C_OPEN( f_string, flags, mode )  ! Open for writing
    ENDIF

    posix_open = fd

 END FUNCTION posix_open



 SUBROUTINE posix_lseek( fid, offset )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                ::  fid     !<
    INTEGER(KIND=C_INT)                ::  my_fid  !<
    INTEGER(KIND=C_SIZE_T), INTENT(IN) ::  offset  !<
    INTEGER(KIND=C_SIZE_T)             ::  retval  !<
    INTEGER(KIND=C_INT)                ::  whence  !<


    my_fid = fid
    whence = seek_set

    retval = C_LSEEK( my_fid, offset, whence )

 END SUBROUTINE posix_lseek



 SUBROUTINE posix_read_int_1d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                             ::  nr_byte  !<
    INTEGER, INTENT(IN)                                 ::  fid      !<
    INTEGER, INTENT(IN)                                 ::  nw       !<
    INTEGER(KIND=iwp), INTENT(IN), TARGET, DIMENSION(:) ::  data     !<
    TYPE(C_PTR)                                         ::  buf      !<


    nr_byte = nw*iwp
    buf     = C_LOC( data )

    CALL posix_read( fid, buf, nr_byte )

 END SUBROUTINE posix_read_int_1d



 SUBROUTINE posix_read_int_2d_i4( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                               ::  nr_byte  !<
    INTEGER, INTENT(IN)                                   ::  fid      !<
    INTEGER, INTENT(IN)                                   ::  nw       !<
    INTEGER(KIND=isp), INTENT(IN), TARGET, DIMENSION(:,:) ::  data     !<
    TYPE(C_PTR)                                           ::  buf      !<


    nr_byte = nw * isp
    buf     = C_LOC( data )

    CALL posix_read( fid, buf, nr_byte )

 END SUBROUTINE posix_read_int_2d_i4



 SUBROUTINE posix_read_int_2d_i8 ( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                               ::  nr_byte  !<
    INTEGER, INTENT(IN)                                   ::  fid      !<
    INTEGER, INTENT(IN)                                   ::  nw       !<
    INTEGER(KIND=idp), INTENT(IN), TARGET, DIMENSION(:,:) ::  data     !<
    TYPE(C_PTR)                                           ::  buf      !<


    nr_byte = nw * idp
    buf     = C_LOC( data )

    CALL posix_read( fid, buf, nr_byte )

 END SUBROUTINE posix_read_int_2d_i8



 SUBROUTINE posix_read_i4_3d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                                 ::  nr_byte  !<
    INTEGER, INTENT(IN)                                     ::  fid      !<
    INTEGER, INTENT(IN)                                     ::  nw       !<
    INTEGER(KIND=isp), INTENT(IN), TARGET, DIMENSION(:,:,:) ::  data     !<
    TYPE(C_PTR)                                             ::  buf      !<


    nr_byte = nw * isp
    buf     = C_LOC( data )

    CALL posix_read( fid, buf, nr_byte )

 END SUBROUTINE posix_read_i4_3d



 SUBROUTINE posix_read_i8_3d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                                 ::  nr_byte  !<
    INTEGER, INTENT(IN)                                     ::  fid      !<
    INTEGER, INTENT(IN)                                     ::  nw       !<
    INTEGER(KIND=idp), INTENT(IN), TARGET, DIMENSION(:,:,:) ::  data     !<
    TYPE(C_PTR)                                             ::  buf      !<


    nr_byte = nw * idp
    buf     = C_LOC( data )

    CALL posix_read( fid, buf, nr_byte )

 END SUBROUTINE posix_read_i8_3d



 SUBROUTINE posix_read_offset_1d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                                  ::  nr_byte  !<
    INTEGER, INTENT(IN)                                      ::  fid      !<
    INTEGER, INTENT(IN)                                      ::  nw       !<
    INTEGER(KIND=C_SIZE_T), INTENT(IN), TARGET, DIMENSION(:) ::  data     !<
    TYPE(C_PTR)                                              ::  buf      !<


    nr_byte = nw * C_SIZE_T
    buf     = C_LOC( data )

    CALL posix_read( fid, buf, nr_byte )

 END SUBROUTINE posix_read_offset_1d



 SUBROUTINE posix_read_real_1d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                         ::  nr_byte  !<
    INTEGER, INTENT(IN)                             ::  fid      !<
    INTEGER, INTENT(IN)                             ::  nw       !<

    REAL(KIND=wp), INTENT(IN), TARGET, DIMENSION(:) ::  data     !<

    TYPE(C_PTR)                                     ::  buf      !<


    nr_byte = nw * wp
    buf     = C_LOC( data )

    CALL posix_read( fid, buf, nr_byte )

 END SUBROUTINE posix_read_real_1d



 SUBROUTINE posix_read_real_2d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                              ::  nr_byte  !<
    INTEGER, INTENT(IN)                                  ::  fid      !<
    INTEGER, INTENT(IN)                                  ::  nw       !<

    REAL(KIND=wp), INTENT(INOUT), TARGET, DIMENSION(:,:) ::  data     !<

    TYPE(C_PTR)                                          ::  buf      !<


    nr_byte = nw * wp
    buf     = C_LOC( data )

    CALL posix_read( fid, buf, nr_byte )

 END SUBROUTINE posix_read_real_2d



 SUBROUTINE posix_read_real_3d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                                ::  nr_byte  !<
    INTEGER, INTENT(IN)                                    ::  fid      !<
    INTEGER, INTENT(IN)                                    ::  nw       !<

    REAL(KIND=wp), INTENT(INOUT), TARGET, DIMENSION(:,:,:) ::  data     !<

    TYPE(C_PTR)                                            ::  buf      !<


    nr_byte = nw * wp
    buf     = C_LOC( data )

    CALL posix_read( fid, buf, nr_byte )

 END SUBROUTINE posix_read_real_3d



 SUBROUTINE posix_read( fid, buf, nb )

    IMPLICIT NONE

    INTEGER, INTENT(IN)    ::  fid      !<
    INTEGER , INTENT(IN)   ::  nb       !<
    INTEGER(KIND=C_INT)    ::  my_fid   !<
    INTEGER(KIND=C_SIZE_T) ::  nr_byte  !<
    INTEGER(KIND=C_SIZE_T) ::  retval   !<

    TYPE(C_PTR)            ::  buf      !<


    my_fid  = fid
    nr_byte = nb

    retval = C_READ( my_fid, buf, nr_byte )

!
!-- The posix standard says that it is not guaranteed that all bytes are read in one read system
!-- call. If retval is not equal to nr_byte, another system call has to be issued.
!-- However, in all Unix distributions it is commonly accepted, that all bytes are read in one call
!-- during during disk-IO. Therefore, here is only an error query and no reading in a while loop.
    IF ( retval /= nr_byte )  THEN
        WRITE( 6, * ) 'Number of bytes read does not match the number of requested bytes'
        CALL abort
    ENDIF

 END SUBROUTINE posix_read



 SUBROUTINE posix_read_char_array( fid, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), DIMENSION(:)   ::  data      !<
    CHARACTER(LEN=LEN(data)), TARGET ::  data_buf  !<

    INTEGER                          ::  i         !<
    INTEGER, INTENT(IN)              ::  fid       !<
    INTEGER(KIND=C_INT)              ::  my_fid    !<
    INTEGER(KIND=C_SIZE_T)           ::  name_len  !<
    INTEGER(KIND=C_SIZE_T)           ::  retval    !<

    TYPE(C_PTR)                      ::  ptr       !<


    my_fid  = fid

    DO  i = 1, SIZE( data )
       data_buf = data(i)
       name_len = LEN( data(i) )
       ptr      = C_LOC( data_buf(1:1) )
       retval   = C_READ( my_fid, ptr, name_len )
       data(i)  = data_buf
    ENDDO

 END SUBROUTINE posix_read_char_array



 SUBROUTINE posix_write_int_1d( fid, data, nw )

    IMPLICIT NONE

    INTEGER, INTENT(IN)                                 ::  fid      !<
    INTEGER ,INTENT(IN)                                 ::  nw       !<
    INTEGER(KIND=C_INT)                                 ::  my_fid   !<
    INTEGER(KIND=C_SIZE_T)                              ::  nr_byte  !<
    INTEGER(KIND=C_SIZE_T)                              ::  retval   !<

    INTEGER(KIND=iwp), INTENT(IN), TARGET, DIMENSION(:) ::  data     !<

    TYPE(C_PTR)                                         ::  buf      !<


    my_fid  = fid
    nr_byte = nw * iwp
    buf     = C_LOC( data )

    retval = C_WRITE( my_fid, buf, nr_byte )

 END SUBROUTINE posix_write_int_1d



 SUBROUTINE posix_write_int_2d_i4( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                               ::  nr_byte  !<
    INTEGER, INTENT(IN)                                   ::  fid      !<
    INTEGER, INTENT(IN)                                   ::  nw       !<

    INTEGER(KIND=isp), INTENT(IN), TARGET, DIMENSION(:,:) ::  data     !<

    TYPE(C_PTR)                                           :: buf       !<


    nr_byte = nw * isp
    buf     = C_LOC( data )

    CALL posix_write( fid, buf, nr_byte )

 END SUBROUTINE posix_write_int_2d_i4

 SUBROUTINE posix_write_int_2d_i8( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                               ::  nr_byte  !<
    INTEGER, INTENT(IN)                                   ::  fid      !<
    INTEGER, INTENT(IN)                                   ::  nw       !<

    INTEGER(KIND=idp), INTENT(IN), TARGET, DIMENSION(:,:) ::  data     !<

    TYPE(C_PTR)                                           :: buf       !<


    nr_byte = nw * idp
    buf     = C_LOC( data )

    CALL posix_write( fid, buf, nr_byte )

 END SUBROUTINE posix_write_int_2d_i8


 SUBROUTINE posix_write_i4_3d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                                 ::  nr_byte  !<
    INTEGER, INTENT(IN)                                     ::  fid      !<
    INTEGER, INTENT(IN)                                     ::  nw       !<

    INTEGER(KIND=isp), INTENT(IN), TARGET, DIMENSION(:,:,:) ::  data     !<

    TYPE(C_PTR)                                             ::  buf      !<


    nr_byte = nw * isp
    buf     = C_LOC( data )

    CALL posix_write( fid, buf, nr_byte )

 END SUBROUTINE posix_write_i4_3d



 SUBROUTINE posix_write_i8_3d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                                 ::  nr_byte  !<
    INTEGER, INTENT(IN)                                     ::  fid      !<
    INTEGER, INTENT(IN)                                     ::  nw       !<

    INTEGER(KIND=idp), INTENT(IN), TARGET, DIMENSION(:,:,:) ::  data     !<

    TYPE(C_PTR)                                             ::  buf      !<


    nr_byte = nw * idp
    buf     = C_LOC( data )

    CALL posix_write( fid, buf, nr_byte )

 END SUBROUTINE posix_write_i8_3d



 SUBROUTINE posix_write_offset_1d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                                  ::  nr_byte  !<
    INTEGER, INTENT(IN)                                      ::  fid      !<
    INTEGER, INTENT(IN)                                      ::  nw       !<

    INTEGER(KIND=C_SIZE_T), INTENT(IN), TARGET, DIMENSION(:) ::  data     !<

    TYPE(C_PTR)                                              ::  buf      !<


    nr_byte = nw * STORAGE_SIZE( data(1) ) / 8
    buf     = C_LOC( data )

    CALL posix_write(fid, buf, nr_byte )

 END SUBROUTINE posix_write_offset_1d



 SUBROUTINE posix_write_real_1d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                         ::  nr_byte   !<
    INTEGER, INTENT(IN)                             ::  fid       !<
    INTEGER, INTENT(IN)                             ::  nw        !<

    REAL(KIND=wp), INTENT(IN), TARGET, DIMENSION(:) ::  data      !<

    TYPE(C_PTR)                                     ::  buf       !<


    nr_byte = nw * wp
    buf     = C_LOC( data )

    CALL posix_write( fid, buf, nr_byte )

 END SUBROUTINE posix_write_real_1d



 SUBROUTINE posix_write_real_2d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                              ::  nr_byte  !<
    INTEGER, INTENT(IN)                                  ::  fid      !<
    INTEGER, INTENT(IN)                                  ::  nw       !<

    REAL(KIND=wp), INTENT(INOUT), TARGET, DIMENSION(:,:) ::  data     !<

    TYPE(C_PTR)                                          ::  buf      !<


    nr_byte = nw * wp
    buf     = C_LOC( data )

    CALL posix_write( fid, buf, nr_byte )

 END SUBROUTINE posix_write_real_2d



 SUBROUTINE posix_write_real_3d( fid, data, nw )

    IMPLICIT NONE

    INTEGER                                                ::  nr_byte  !<
    INTEGER, INTENT(IN)                                    ::  fid      !<
    INTEGER, INTENT(IN)                                    ::  nw       !<

    REAL(KIND=wp), INTENT(INOUT), TARGET, DIMENSION(:,:,:) ::  data     !<

    TYPE(C_PTR)                                            ::  buf      !<


    nr_byte = nw * wp
    buf     = C_LOC( data )

    CALL posix_write( fid, buf, nr_byte )

 END SUBROUTINE posix_write_real_3d



 SUBROUTINE posix_write( fid, buf, nb )

    IMPLICIT NONE

    INTEGER, INTENT(IN)    ::  fid      !<
    INTEGER , INTENT(IN)   ::  nb       !<
    INTEGER(KIND=C_INT)    ::  my_fid   !<
    INTEGER(KIND=C_SIZE_T) ::  nr_byte  !<
    INTEGER(KIND=C_SIZE_T) ::  retval   !<

    TYPE(C_PTR)            ::  buf      !<


    my_fid  = fid
    nr_byte = nb

    retval = C_WRITE( my_fid, buf, nr_byte )

    IF ( retval /= nr_byte )  THEN
       WRITE( 6, * ) 'Number of bytes to write does not match the number of requested bytes'
       CALL abort
    ENDIF

 END SUBROUTINE posix_write



 SUBROUTINE posix_write_char_array( fid, data )

    IMPLICIT NONE

    CHARACTER(LEN=*), DIMENSION(:)     ::  data      !<
    CHARACTER(LEN=LEN(data)+1), TARGET ::  data_buf  !<

    INTEGER                            ::  i         !<
    INTEGER, INTENT(IN)                ::  fid       !<
    INTEGER(KIND=C_INT)                ::  my_fid    !<
    INTEGER(KIND=C_SIZE_T)             ::  name_len  !<
    INTEGER(KIND=C_SIZE_T)             ::  retval    !<

    TYPE(C_PTR)                        ::  ptr       !<


    my_fid  = fid

    DO  i = 1, SIZE( data )
       data_buf = data(i) // CHAR( 0 )
       name_len = LEN( data(i) )
       ptr      = C_LOC( data_buf(1:1) )
       retval   = C_WRITE( my_fid, ptr, name_len )
    ENDDO

 END SUBROUTINE posix_write_char_array



 SUBROUTINE posix_close( fid )

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  fid     !<
    INTEGER(KIND=C_INT) ::  my_fid  !<
    INTEGER(KIND=C_INT) ::  retval  !<


    my_fid = fid

    retval = C_CLOSE( my_fid )

 END SUBROUTINE posix_close


 END MODULE posix_interface
