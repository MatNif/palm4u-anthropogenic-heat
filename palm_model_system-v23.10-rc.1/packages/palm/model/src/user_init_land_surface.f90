!> @file user_init_land_surface.f90
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
!
! Description:
! ------------
!> Execution of user-defined actions to initiate the land surface model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_init_land_surface

    USE kinds

    USE land_surface_model_mod

    USE surface_mod

    IMPLICIT NONE

!    INTEGER(iwp) ::  i  !< grid index
!    INTEGER(iwp) ::  j  !< grid index
!    INTEGER(iwp) ::  m  !< running index on 1D wall-type grid

!
!-- Here the user-defined land surface initialization actions follow.
!-- Example: set roughness length at natural land-surface
!     DO  m = 1, surf_lsm%ns
!        surf_lsm%z0(m) = 0.1_wp
!     ENDDO


 END SUBROUTINE user_init_land_surface

