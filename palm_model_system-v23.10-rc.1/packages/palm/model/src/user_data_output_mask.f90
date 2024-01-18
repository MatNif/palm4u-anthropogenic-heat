!> @file user_data_output_mask.f90
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
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with
!> indices (i,j,k) for masked data output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_mask( av, variable, found, local_pf, mid )

    USE control_parameters

    USE indices

    USE kinds

    USE user

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  av              !<
    INTEGER(iwp) ::  mid             !< masked output running index
!    INTEGER(iwp) ::  i               !<
!    INTEGER(iwp) ::  j               !<
!    INTEGER(iwp) ::  k               !<
!    INTEGER(iwp) ::  topo_top_index  !< k index of highest horizontal surface

    LOGICAL ::  found  !<

    REAL(wp), DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  local_pf  !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av == 0  .OR.                                                                             &
         local_pf(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) == 0.0_wp )  CONTINUE


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!--    Uncomment and extend the following lines, if necessary.
!--    The arrays for storing the user defined quantities (here u2 and u2_av) have to be declared
!--    and defined by the user!
!--    Sample for user-defined output:
!       CASE ( 'u2' )
!          IF ( av == 0 )  THEN
!             IF ( .NOT. mask_surface(mid) )  THEN
!!
!!--             Default masked output
!                DO  i = 1, mask_size_l(mid,1)
!                   DO  j = 1, mask_size_l(mid,2)
!                      DO  k = 1, mask_size_l(mid,3)
!                         local_pf(i,j,k) = u2(mask_k(mid,k),                                       &
!                                              mask_j(mid,j),                                       &
!                                              mask_i(mid,i))
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ELSE
!!
!!--             Terrain-following masked output
!                DO  i = 1, mask_size_l(mid,1)
!                   DO  j = 1, mask_size_l(mid,2)
!!
!!--                   Get k index of highest horizontal surface
!                      topo_top_index = topo_top_ind( mask_j(mid,j), mask_i(mid,i), 1 )
!!
!!--                   Save output array
!                      DO  k = 1, mask_size_l(mid,3)
!                         local_pf(i,j,k) = u2(MIN( topo_top_index + mask_k(mid,k), nzt+1 ),        &
!                                              mask_j(mid,j), mask_i(mid,i) )
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ENDIF
!          ELSE
!             IF ( .NOT. mask_surface(mid) )  THEN
!!
!!--             Default masked output
!                DO  i = 1, mask_size_l(mid,1)
!                   DO  j = 1, mask_size_l(mid,2)
!                      DO  k = 1, mask_size_l(mid,3)
!                          local_pf(i,j,k) = u2_av(mask_k(mid,k), mask_j(mid,j), mask_i(mid,i) )
!                       ENDDO
!                    ENDDO
!                 ENDDO
!             ELSE
!!
!!--             Terrain-following masked output
!                DO  i = 1, mask_size_l(mid,1)
!                   DO  j = 1, mask_size_l(mid,2)
!!
!!--                   Get k index of highest horizontal surface
!                      topo_top_index = topo_top_ind( mask_j(mid,j), mask_i(mid,i), 1 )
!!
!!--                   Save output array
!                      DO  k = 1, mask_size_l(mid,3)
!                         local_pf(i,j,k) = u2_av( MIN( topo_top_index+mask_k(mid,k), nzt+1 ),      &
!                                                  mask_j(mid,j), mask_i(mid,i) )
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ENDIF
!          ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE user_data_output_mask
