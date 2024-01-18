!> @file posix_interface_mod.f90
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
! Copyright 1997-2021  Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Matthias Suehring
! @author Klaus Ketelsen
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interface to posix systemcalls. The module allows to call the following
!> C system calls from FORTRAN: ftell, lseek.
!------------------------------------------------------------------------------!
 MODULE posix_interface

    INTERFACE
       INTEGER (C_SIZE_T) FUNCTION C_LSEEK (fd, offset, whence)                &
                                   BIND(C, NAME='lseek')

          USE ISO_C_BINDING

          IMPLICIT NONE

          INTEGER(KIND=C_INT),    VALUE ::  fd
          INTEGER(KIND=C_SIZE_T), VALUE ::  offset
          INTEGER(KIND=C_INT),    VALUE ::  whence

       END FUNCTION C_LSEEK

    END INTERFACE

    INTERFACE
       INTEGER (C_SIZE_T) FUNCTION C_FTELL ( fd )                              &
                                   BIND(C, NAME='ftell')

          USE ISO_C_BINDING

          IMPLICIT NONE

          INTEGER(KIND=C_INT),    VALUE ::  fd

       END FUNCTION C_FTELL

    END INTERFACE

    PUBLIC posix_lseek, posix_ftell

    CONTAINS
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interface for the C-routine lseek.
!------------------------------------------------------------------------------!
    SUBROUTINE posix_lseek( fid, offset )

       USE ISO_C_BINDING

       IMPLICIT NONE

       INTEGER,INTENT(IN)                     ::  fid    !< file unit
       INTEGER(KIND=C_SIZE_T),INTENT(IN)      ::  offset !< file offset from the beginning

       INTEGER(KIND=C_INT)                    ::  my_fid !< file unit
       INTEGER(KIND=C_SIZE_T)                 ::  retval !< return value
       INTEGER(KIND=C_INT)                    ::  whence !< mode, here start from the beginning

       my_fid = fid
       whence = 0

       retval = C_LSEEK( fid, offset, whence )

    END SUBROUTINE posix_lseek
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interface for the C-routine ftell.
!------------------------------------------------------------------------------!
    SUBROUTINE posix_ftell( fid, filepos )

       USE ISO_C_BINDING

       IMPLICIT NONE

       INTEGER,INTENT(IN)                      ::  fid     !< file unit
       INTEGER(KIND=C_SIZE_T), INTENT(INOUT)   ::  filepos !< file position

       filepos = C_FTELL( fid )

    END SUBROUTINE posix_ftell

 END MODULE posix_interface
