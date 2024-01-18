!> @file boundary_settings_mod.f90
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
!> This module contains some general settings of specific boundary conditions used by various
!> modules.
!--------------------------------------------------------------------------------------------------!
 MODULE boundary_settings_mod

    USE control_parameters,                                                                        &
        ONLY:  bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg

    USE kinds

    IMPLICIT NONE


    INTERFACE set_lateral_neumann_bc
       MODULE PROCEDURE set_lateral_neumann_bc_int1
       MODULE PROCEDURE set_lateral_neumann_bc_int4
       MODULE PROCEDURE set_lateral_neumann_bc_real
    END INTERFACE set_lateral_neumann_bc

!
!-- Public routines
    PUBLIC set_lateral_neumann_bc

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set Neumann boundary conditions at lateral boundaries for 1-byte 2D integer arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE set_lateral_neumann_bc_int1( ar )

    INTEGER(iwp) ::  i                        !< running index over ghost points

    INTEGER(ibp) ::  ar(nysg:nyng,nxlg:nxrg)  !< treated array


!
!-- Neumann-conditions at inflow/outflow/nested boundaries.
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       DO  i = nbgp, 1, -1
          ar(:,nxl-i) = ar(:,nxl)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
       DO  i = 1, nbgp
          ar(:,nxr+i) = ar(:,nxr)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       DO  i = nbgp, 1, -1
          ar(nys-i,:) = ar(nys,:)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
       DO  i = 1, nbgp
          ar(nyn+i,:) = ar(nyn,:)
       ENDDO
    ENDIF

 END SUBROUTINE set_lateral_neumann_bc_int1


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set Neumann boundary conditions at lateral boundaries for 4-byte 2D integer arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE set_lateral_neumann_bc_int4( ar )

    INTEGER(iwp) ::  i                        !< running index over ghost points

    INTEGER(iwp) ::  ar(nysg:nyng,nxlg:nxrg)  !< treated array


!
!-- Neumann-conditions at inflow/outflow/nested boundaries.
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       DO  i = nbgp, 1, -1
          ar(:,nxl-i) = ar(:,nxl)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
       DO  i = 1, nbgp
          ar(:,nxr+i) = ar(:,nxr)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       DO  i = nbgp, 1, -1
          ar(nys-i,:) = ar(nys,:)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
       DO  i = 1, nbgp
          ar(nyn+i,:) = ar(nyn,:)
       ENDDO
    ENDIF

 END SUBROUTINE set_lateral_neumann_bc_int4


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set Neumann boundary conditions at lateral boundaries for real-type 2D arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE set_lateral_neumann_bc_real( ar )

    INTEGER(iwp) ::  i                    !< running index over ghost points

    REAL(wp) ::  ar(nysg:nyng,nxlg:nxrg)  !< treated array

!
!-- Neumann-conditions at inflow/outflow/nested boundaries.
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       DO  i = nbgp, 1, -1
          ar(:,nxl-i) = ar(:,nxl)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
       DO  i = 1, nbgp
          ar(:,nxr+i) = ar(:,nxr)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       DO  i = nbgp, 1, -1
          ar(nys-i,:) = ar(nys,:)
       ENDDO
    ENDIF
    IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
       DO  i = 1, nbgp
          ar(nyn+i,:) = ar(nyn,:)
       ENDDO
    ENDIF

 END SUBROUTINE set_lateral_neumann_bc_real


 END MODULE boundary_settings_mod
