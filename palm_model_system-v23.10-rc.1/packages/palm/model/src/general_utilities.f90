!> @file general_utilities.f90
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
! Copyright 2021-2022 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> A collection of handy utilities.
!--------------------------------------------------------------------------------------------------!
 MODULE general_utilities

    USE kinds

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               ny,                                                                                 &
               nz

    IMPLICIT NONE

!
!-- Public functions
    PUBLIC                                                                                         &
       gridpoint_id,                                                                               &
       interpolate_linear

    INTERFACE gridpoint_id
       MODULE PROCEDURE gridpoint_id_2d
       MODULE PROCEDURE gridpoint_id_3d
    END INTERFACE gridpoint_id


    INTERFACE interpolate_linear
       MODULE PROCEDURE interpolate_linear_0d_wp
    END INTERFACE interpolate_linear

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> This functions computes an unique ID for each given grid point (j,i).
!--------------------------------------------------------------------------------------------------!
 FUNCTION gridpoint_id_2d( j, i )

    INTEGER(idp) ::  gridpoint_id_2d !< grid point ID, unique for each (j,i)
    INTEGER(iwp) ::  i               !< grid index in x-direction
    INTEGER(iwp) ::  j               !< grid index in y-direction


    gridpoint_id_2d = i + ( nx + 1 ) * j

 END FUNCTION gridpoint_id_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> This functions computes an unique ID for each given grid point (k,j,i).
!--------------------------------------------------------------------------------------------------!
 FUNCTION gridpoint_id_3d( k, j, i )

    INTEGER(idp) ::  gridpoint_id_3d !< grid point ID, unique for each (k,j,i)
    INTEGER(iwp) ::  i               !< grid index in x-direction
    INTEGER(iwp) ::  j               !< grid index in y-direction
    INTEGER(iwp) ::  k               !< grid index in k-direction


    gridpoint_id_3d = i * ( ny + 1 ) * ( nz + 1 ) + j * ( nz + 1 ) + k + 1

 END FUNCTION gridpoint_id_3d


!--------------------------------------------------------------------------------------------------!
! Description:
!--------------------------------------------------------------------------------------------------!
!> Interpolation function, used to interpolate between two real-type scalar values, e.g. in space
!> or time.
!--------------------------------------------------------------------------------------------------!
 FUNCTION interpolate_linear_0d_wp( var_x1, var_x2, fac  )

    REAL(wp)            :: fac                       !< interpolation factor
    REAL(wp)            :: interpolate_linear_0d_wp  !< interpolated value
    REAL(wp)            :: var_x1                    !< value at x1
    REAL(wp)            :: var_x2                    !< value at x2


    interpolate_linear_0d_wp = ( 1.0_wp - fac ) * var_x1 + fac * var_x2

 END FUNCTION interpolate_linear_0d_wp

 END MODULE general_utilities
