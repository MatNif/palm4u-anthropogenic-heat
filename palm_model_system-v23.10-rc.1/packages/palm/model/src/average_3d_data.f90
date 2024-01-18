!> @file average_3d_data.f90
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
!> Time-averaging of 3d-data-arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE average_3d_data
 

    USE averaging

    USE control_parameters,                                                                        &
        ONLY:  average_count_3d,                                                                   &
               doav,                                                                               &
               doav_n,                                                                             &
               varnamelength

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz_2d

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nzb,                                                                                &
               nzt

    USE kinds

    USE module_interface,                                                                          &
        ONLY:  module_interface_3d_data_averaging




    IMPLICIT NONE

    INTEGER(iwp) ::  i   !< loop index
    INTEGER(iwp) ::  ii  !< loop index
    INTEGER(iwp) ::  j   !< loop index
    INTEGER(iwp) ::  k   !< loop index

    CHARACTER (LEN=varnamelength) ::  trimvar  !< TRIM of output-variable string


    CALL cpu_log (log_point(35),'average_3d_data','start')

!
!-- Check, if averaging is necessary
    IF ( average_count_3d <= 1 )  RETURN

!
!-- Loop of all variables to be averaged.
    DO  ii = 1, doav_n

       trimvar = TRIM( doav(ii) )

!
!--    Store the array chosen on the temporary array.
       SELECT CASE ( trimvar )

          CASE ( 'e' )
             IF ( ALLOCATED( e_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         e_av(k,j,i) = e_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ghf*' )
             IF ( ALLOCATED( ghf_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      ghf_av(j,i) = ghf_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( ghf_av )
             ENDIF

          CASE ( 'qsws*' )
             IF ( ALLOCATED( qsws_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      qsws_av(j,i) = qsws_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( qsws_av )
             ENDIF

          CASE ( 'thetal' )
             IF ( ALLOCATED( lpt_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         lpt_av(k,j,i) = lpt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'lwp*' )
             IF ( ALLOCATED( lwp_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      lwp_av(j,i) = lwp_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

         CASE ( 'ol*' )
             IF ( ALLOCATED( ol_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      ol_av(j,i) = ol_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( ol_av )
             ENDIF

          CASE ( 'p' )
             IF ( ALLOCATED( p_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         p_av(k,j,i) = p_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pc' )
             IF ( ALLOCATED( pc_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         pc_av(k,j,i) = pc_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pr' )
             IF ( ALLOCATED( pr_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         pr_av(k,j,i) = pr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pres_drag_x*' )
             IF ( ALLOCATED( pres_drag_x_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      pres_drag_x_av(j,i) = pres_drag_x_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pres_drag_y*' )
             IF ( ALLOCATED( pres_drag_y_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      pres_drag_y_av(j,i) = pres_drag_y_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'theta' )
             IF ( ALLOCATED( pt_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         pt_av(k,j,i) = pt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'q' )
             IF ( ALLOCATED( q_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         q_av(k,j,i) = q_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql' )
             IF ( ALLOCATED( ql_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_av(k,j,i) = ql_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_c' )
             IF ( ALLOCATED( ql_c_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_c_av(k,j,i) = ql_c_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_v' )
             IF ( ALLOCATED( ql_v_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_v_av(k,j,i) = ql_v_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_vp' )
             IF ( ALLOCATED( ql_vp_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_vp_av(k,j,i) = ql_vp_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qv' )
             IF ( ALLOCATED( qv_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         qv_av(k,j,i) = qv_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

         CASE ( 'qsurf*' )
             IF ( ALLOCATED( qsurf_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      qsurf_av(j,i) = qsurf_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( qsurf_av )
             ENDIF

         CASE ( 'r_a*' )
             IF ( ALLOCATED( r_a_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      r_a_av(j,i) = r_a_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( r_a_av )
             ENDIF

          CASE ( 's' )
             IF ( ALLOCATED( s_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         s_av(k,j,i) = s_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

         CASE ( 'shf*' )
             IF ( ALLOCATED( shf_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      shf_av(j,i) = shf_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( shf_av )
             ENDIF

          CASE ( 'ssurf*' )
             IF ( ALLOCATED( ssurf_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      ssurf_av(j,i) = ssurf_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( ssurf_av )
             ENDIF

          CASE ( 'ssws*' )
             IF ( ALLOCATED( ssws_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      ssws_av(j,i) = ssws_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( ssws_av )
             ENDIF

          CASE ( 't*' )
             IF ( ALLOCATED( ts_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      ts_av(j,i) = ts_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( ts_av )
             ENDIF

         CASE ( 'tsurf*' )
             IF ( ALLOCATED( tsurf_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      tsurf_av(j,i) = tsurf_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( tsurf_av )
             ENDIF

          CASE ( 'u' )
             IF ( ALLOCATED( u_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         u_av(k,j,i) = u_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'us*' )
             IF ( ALLOCATED( us_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      us_av(j,i) = us_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( us_av )
             ENDIF

          CASE ( 'v' )
             IF ( ALLOCATED( v_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         v_av(k,j,i) = v_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'thetav' )
             IF ( ALLOCATED( vpt_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vpt_av(k,j,i) = vpt_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'w' )
             IF ( ALLOCATED( w_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         w_av(k,j,i) = w_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'z0*' )
             IF ( ALLOCATED( z0_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      z0_av(j,i) = z0_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( z0_av )
             ENDIF

          CASE ( 'z0h*' )
             IF ( ALLOCATED( z0h_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      z0h_av(j,i) = z0h_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( z0h_av )
             ENDIF

          CASE ( 'z0q*' )
             IF ( ALLOCATED( z0q_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      z0q_av(j,i) = z0q_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( z0q_av )
             ENDIF

          CASE DEFAULT
!
!--          Averaging of data from all other modules
             CALL module_interface_3d_data_averaging( 'average', trimvar )

       END SELECT

    ENDDO

!
!-- Reset the counter
    average_count_3d = 0_iwp

    CALL cpu_log( log_point(35), 'average_3d_data', 'stop' )


 END SUBROUTINE average_3d_data
