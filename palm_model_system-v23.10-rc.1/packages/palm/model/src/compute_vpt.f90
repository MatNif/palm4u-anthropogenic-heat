!> @file compute_vpt.f90
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
! -------------
!> Computation of the virtual potential temperature 
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE compute_vpt
 
    USE arrays_3d,                                                                                 &
        ONLY:  d_exner,                                                                            &
               pt,                                                                                 &
               q,                                                                                  &
               qf,                                                                                 &
               ql,                                                                                 &
               vpt

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  ls_d_cp,                                                                            &
               lv_d_cp

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model,                                                                   &
               microphysics_ice_phase

    USE control_parameters,                                                                        &
        ONLY:  cloud_droplets

    USE indices,                                                                                   &
        ONLY:  nzb,                                                                                &
               nzt

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) :: k   !< 

    IF ( .NOT. bulk_cloud_model  .AND.  .NOT. cloud_droplets )  THEN
       vpt = pt * ( 1.0_wp + 0.61_wp * q )
    ELSEIF ( bulk_cloud_model  .AND.  .NOT. microphysics_ice_phase )  THEN
       DO  k = nzb, nzt+1
              vpt(k,:,:) = ( pt(k,:,:) + d_exner(k) * lv_d_cp * ql(k,:,:) ) *                      &
                           ( 1.0_wp + 0.61_wp * q(k,:,:) - 1.61_wp *  ql(k,:,:)  )
       ENDDO
    ELSEIF ( bulk_cloud_model  .AND.  microphysics_ice_phase )  THEN
       DO  k = nzb, nzt+1
          vpt(k,:,:) = ( pt(k,:,:) + d_exner(k) * lv_d_cp * ql(k,:,:)   +                          &
                                     d_exner(k) * ls_d_cp * qf(k,:,:) ) *                          &
                       ( 1.0_wp + 0.61_wp * q(k,:,:) - 1.61_wp * ( ql(k,:,:) + qf(k,:,:) ) )
       ENDDO
    ELSE
       vpt = pt * ( 1.0_wp + 0.61_wp * q - ql ) 
    ENDIF

 END SUBROUTINE compute_vpt
