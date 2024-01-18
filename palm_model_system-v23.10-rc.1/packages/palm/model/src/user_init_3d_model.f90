!> @file user_init_3d_model.f90
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
!> Allows the complete initialization of the 3d model.
!>
!> @attention The user is responsible to set at least all those quantities which
!>            are normally set within init_3d_model!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_init_3d_model

    USE arrays_3d

    USE control_parameters

    USE indices

    USE kinds

    USE surface_mod

    USE user

    IMPLICIT NONE

!    INTEGER(iwp) ::  l  !< running index surface orientation
!    INTEGER(iwp) ::  m  !< running index surface elements

!
!-- Initialization of surface-related quantities.
!-- The following example shows required initialization of surface quantitites at default-type
!-- surfaces. The following loop includes all surface orientations. If only surfaces with
!-- specific surface orientation shall be treated, the user can distinguish surfaces via
!-- the attribute %koff(m), or alternatively by the logical attributions upward, downward,
!-- northward, southward, eastward, and westward.
!   DO  m = 1, surf_def%ns
!      surf_def%ol(m)   = ...    ! Obukhov length
!      surf_def%us(m  ) = ...    ! friction velocity
!      surf_def%usws(m) = ...    ! vertical momentum flux, u-component
!      surf_def%vsws(m) = ...    ! vertical momentum flux, v-component
!      surf_def%z0(m)   = ...    ! roughness length for momentum
!      IF ( .NOT. neutral )  THEN
!         surf_def%ts(m)   = ... ! scaling parameter
!         surf_def%shf(m)  = ... ! surface sensible heat flux
!         surf_def%z0h(m)  = ... ! roughness length for heat
!      ENDIF
!      IF ( humditiy )  THEN
!         surf_def%qs(m)   = ... ! scaling parameter
!         surf_def%qsws(m) = ... ! surface latent heat flux
!         surf_def%z0q(m)  = ... ! roughness length for moisture
!      ENDIF
!      IF ( passive_scalar )  THEN
!         surf_def%ss(m)   = ... ! scaling parameter
!         surf_def%ssws(m) = ... ! surface latent heat flux
!      ENDIF
!   ENDDO
!
!-- Same for natural and urban type surfaces
!   DO  m = 1, surf_lsm%ns
!      ...
!   ENDDO
!   DO  m = 1, surf_usm%ns
!      ...
!   ENDDO
!
!
!-- In the following, initialize 3D quantities, e.g. u, v, w, pt, etc..

 END SUBROUTINE user_init_3d_model

