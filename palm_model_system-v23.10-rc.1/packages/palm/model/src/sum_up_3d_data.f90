!> @file sum_up_3d_data.f90
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
! Copyright 2022-2022 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> Sum-up the values of 3d-arrays. The real averaging is later done in routine average_3d_data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sum_up_3d_data

    USE arrays_3d,                                                                                 &
        ONLY:  dzw,                                                                                &
               d_exner,                                                                            &
               e,                                                                                  &
               heatflux_output_conversion,                                                         &
               p,                                                                                  &
               pt,                                                                                 &
               q,                                                                                  &
               ql,                                                                                 &
               ql_c,                                                                               &
               ql_v,                                                                               &
               s,                                                                                  &
               scalarflux_output_conversion,                                                       &
               u,                                                                                  &
               v,                                                                                  &
               vpt,                                                                                &
               w,                                                                                  &
               waterflux_output_conversion

    USE averaging

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  lv_d_cp

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model

    USE control_parameters,                                                                        &
        ONLY:  average_count_3d,                                                                   &
               doav,                                                                               &
               doav_n,                                                                             &
               dz,                                                                                 &
               interpolate_to_grid_center,                                                         &
               rho_surface,                                                                        &
               urban_surface,                                                                      &
               varnamelength

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

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
               nzt,                                                                                &
               topo_top_ind

    USE kinds

    USE module_interface,                                                                          &
        ONLY:  module_interface_3d_data_averaging

    USE particle_attributes,                                                                       &
        ONLY:  grid_particles,                                                                     &
               number_of_particles,                                                                &
               particles,                                                                          &
               prt_count

    USE surface_mod,                                                                               &
        ONLY:  ind_pav_green,                                                                      &
               ind_veg_wall,                                                                       &
               ind_wat_win,                                                                        &
               surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_usm

    USE urban_surface_mod,                                                                         &
        ONLY:  usm_3d_data_averaging


    IMPLICIT NONE

    CHARACTER(LEN=varnamelength) ::  trimvar  !< TRIM of output-variable string

    INTEGER(iwp) ::  i       !< grid index x direction
    INTEGER(iwp) ::  ii      !< running index
    INTEGER(iwp) ::  j       !< grid index y direction
    INTEGER(iwp) ::  k       !< grid index x direction
    INTEGER(iwp) ::  kl      !< vertical index used to limit lower interpolation index
    INTEGER(iwp) ::  m       !< running index over surfacle elements
    INTEGER(iwp) ::  n       !< running index over number of particles per grid box

    REAL(wp) ::  mean_r   !< mean-particle radius witin grid box
    REAL(wp) ::  s_r2     !< mean-particle radius witin grid box to the power of two
    REAL(wp) ::  s_r3     !< mean-particle radius witin grid box to the power of three



    CALL cpu_log (log_point(34),'sum_up_3d_data','start')

!
!-- Allocate and initialize the summation arrays if called for the very first time or the first time
!-- after average_3d_data has been called (some or all of the arrays may have been already allocated
!-- in rrd_local)
    IF ( average_count_3d == 0 )  THEN

       DO  ii = 1, doav_n

          trimvar = TRIM( doav(ii) )

          SELECT CASE ( trimvar )

             CASE ( 'ghf*' )
                IF ( .NOT. ALLOCATED( ghf_av ) )  THEN
                   ALLOCATE( ghf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ghf_av = 0.0_wp

             CASE ( 'e' )
                IF ( .NOT. ALLOCATED( e_av ) )  THEN
                   ALLOCATE( e_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                e_av = 0.0_wp

             CASE ( 'lwp*' )
                IF ( .NOT. ALLOCATED( lwp_av ) )  THEN
                   ALLOCATE( lwp_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                lwp_av = 0.0_wp

             CASE ( 'ol*' )
                IF ( .NOT. ALLOCATED( ol_av ) )  THEN
                   ALLOCATE( ol_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ol_av = 0.0_wp

             CASE ( 'p' )
                IF ( .NOT. ALLOCATED( p_av ) )  THEN
                   ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                p_av = 0.0_wp

             CASE ( 'pc' )
                IF ( .NOT. ALLOCATED( pc_av ) )  THEN
                   ALLOCATE( pc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pc_av = 0.0_wp

             CASE ( 'pr' )
                IF ( .NOT. ALLOCATED( pr_av ) )  THEN
                   ALLOCATE( pr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pr_av = 0.0_wp

             CASE ( 'pres_drag_x*' )
                IF ( .NOT. ALLOCATED( pres_drag_x_av ) )  THEN
                   ALLOCATE( pres_drag_x_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                pres_drag_x_av = 0.0_wp

             CASE ( 'pres_drag_y*' )
                IF ( .NOT. ALLOCATED( pres_drag_y_av ) )  THEN
                   ALLOCATE( pres_drag_y_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                pres_drag_y_av = 0.0_wp

             CASE ( 'q' )
                IF ( .NOT. ALLOCATED( q_av ) )  THEN
                   ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                q_av = 0.0_wp

             CASE ( 'ql' )
                IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                   ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_av = 0.0_wp

             CASE ( 'ql_c' )
                IF ( .NOT. ALLOCATED( ql_c_av ) )  THEN
                   ALLOCATE( ql_c_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_c_av = 0.0_wp

             CASE ( 'ql_v' )
                IF ( .NOT. ALLOCATED( ql_v_av ) )  THEN
                   ALLOCATE( ql_v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_v_av = 0.0_wp

             CASE ( 'ql_vp' )
                IF ( .NOT. ALLOCATED( ql_vp_av ) )  THEN
                   ALLOCATE( ql_vp_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_vp_av = 0.0_wp

             CASE ( 'qsurf*' )
                IF ( .NOT. ALLOCATED( qsurf_av ) )  THEN
                   ALLOCATE( qsurf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsurf_av = 0.0_wp

             CASE ( 'qsws*' )
                IF ( .NOT. ALLOCATED( qsws_av ) )  THEN
                   ALLOCATE( qsws_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                qsws_av = 0.0_wp

             CASE ( 'qv' )
                IF ( .NOT. ALLOCATED( qv_av ) )  THEN
                   ALLOCATE( qv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qv_av = 0.0_wp

             CASE ( 'r_a*' )
                IF ( .NOT. ALLOCATED( r_a_av ) )  THEN
                   ALLOCATE( r_a_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                r_a_av = 0.0_wp

             CASE ( 's' )
                IF ( .NOT. ALLOCATED( s_av ) )  THEN
                   ALLOCATE( s_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                s_av = 0.0_wp

             CASE ( 'shf*' )
                IF ( .NOT. ALLOCATED( shf_av ) )  THEN
                   ALLOCATE( shf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                shf_av = 0.0_wp

             CASE ( 'ssurf*' )
                IF ( .NOT. ALLOCATED( ssurf_av ) )  THEN
                   ALLOCATE( ssurf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ssurf_av = 0.0_wp

             CASE ( 'ssws*' )
                IF ( .NOT. ALLOCATED( ssws_av ) )  THEN
                   ALLOCATE( ssws_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ssws_av = 0.0_wp

             CASE ( 't*' )
                IF ( .NOT. ALLOCATED( ts_av ) )  THEN
                   ALLOCATE( ts_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                ts_av = 0.0_wp

             CASE ( 'theta' )
                IF ( .NOT. ALLOCATED( pt_av ) )  THEN
                   ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pt_av = 0.0_wp

             CASE ( 'thetal' )
                IF ( .NOT. ALLOCATED( lpt_av ) )  THEN
                   ALLOCATE( lpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                lpt_av = 0.0_wp

             CASE ( 'tsurf*' )
                IF ( .NOT. ALLOCATED( tsurf_av ) )  THEN
                   ALLOCATE( tsurf_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                tsurf_av = 0.0_wp

             CASE ( 'u' )
                IF ( .NOT. ALLOCATED( u_av ) )  THEN
                   ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                u_av = 0.0_wp

             CASE ( 'us*' )
                IF ( .NOT. ALLOCATED( us_av ) )  THEN
                   ALLOCATE( us_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                us_av = 0.0_wp

             CASE ( 'v' )
                IF ( .NOT. ALLOCATED( v_av ) )  THEN
                   ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                v_av = 0.0_wp

             CASE ( 'thetav' )
                IF ( .NOT. ALLOCATED( vpt_av ) )  THEN
                   ALLOCATE( vpt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                vpt_av = 0.0_wp

             CASE ( 'w' )
                IF ( .NOT. ALLOCATED( w_av ) )  THEN
                   ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                w_av = 0.0_wp

             CASE ( 'z0*' )
                IF ( .NOT. ALLOCATED( z0_av ) )  THEN
                   ALLOCATE( z0_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                z0_av = 0.0_wp

             CASE ( 'z0h*' )
                IF ( .NOT. ALLOCATED( z0h_av ) )  THEN
                   ALLOCATE( z0h_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                z0h_av = 0.0_wp

             CASE ( 'z0q*' )
                IF ( .NOT. ALLOCATED( z0q_av ) )  THEN
                   ALLOCATE( z0q_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                z0q_av = 0.0_wp


             CASE DEFAULT

!
!--             Allocating and initializing data arrays for all other modules
                CALL module_interface_3d_data_averaging( 'allocate', trimvar )


          END SELECT

       ENDDO

    ENDIF

!
!-- Loop of all variables to be averaged.
    DO  ii = 1, doav_n

       trimvar = TRIM( doav(ii) )
!
!--    Store the array chosen on the temporary array.
       SELECT CASE ( trimvar )

          CASE ( 'ghf*' )
             IF ( ALLOCATED( ghf_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
!
!--                   Take ground-heat flux from the uppermost upward-facing surface, which is
!--                   is either an LSM or an USM surface.
                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  ghf_av(j,i) = ghf_av(j,i) + surf_lsm%ghf(m)
                      ENDDO
!
!--                   For urban-type surfaces, aggregate resistance from tile approach.
                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            ghf_av(j,i) = ghf_av(j,i) +                                            &
                                     surf_usm%frac(m,ind_veg_wall)  * surf_usm%wghf_eb(m)       +  &
                                     surf_usm%frac(m,ind_pav_green) * surf_usm%wghf_eb_green(m) +  &
                                     surf_usm%frac(m,ind_wat_win)   * surf_usm%wghf_eb_window(m)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'e' )
             IF ( ALLOCATED( e_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         e_av(k,j,i) = e_av(k,j,i) + e(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'thetal' )
             IF ( ALLOCATED( lpt_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         lpt_av(k,j,i) = lpt_av(k,j,i) + pt(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'lwp*' )
             IF ( ALLOCATED( lwp_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      lwp_av(j,i) = lwp_av(j,i) + SUM( ql(nzb:nzt,j,i) * dzw(1:nzt+1) )            &
                                                                       * rho_surface
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ol*' )
             IF ( ALLOCATED( ol_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  ol_av(j,i) = ol_av(j,i) + surf_def%ol(m)
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  ol_av(j,i) = ol_av(j,i) + surf_lsm%ol(m)
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  ol_av(j,i) = ol_av(j,i) + surf_usm%ol(m)
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'p' )
             IF ( ALLOCATED( p_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         p_av(k,j,i) = p_av(k,j,i) + p(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pc' )
             IF ( ALLOCATED( pc_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         pc_av(k,j,i) = pc_av(k,j,i) + prt_count(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pr' )
             IF ( ALLOCATED( pr_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         number_of_particles = prt_count(k,j,i)
                         IF ( number_of_particles <= 0 )  CYCLE
                         particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                         s_r2 = 0.0_wp
                         s_r3 = 0.0_wp

                         DO  n = 1, number_of_particles
                            IF ( particles(n)%particle_mask )  THEN
                               s_r2 = s_r2 + particles(n)%radius**2 * particles(n)%weight_factor
                               s_r3 = s_r3 + particles(n)%radius**3 * particles(n)%weight_factor
                            ENDIF
                         ENDDO

                         IF ( s_r2 > 0.0_wp )  THEN
                            mean_r = s_r3 / s_r2
                         ELSE
                            mean_r = 0.0_wp
                         ENDIF
                         pr_av(k,j,i) = pr_av(k,j,i) + mean_r
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pres_drag_x*' )
             IF ( ALLOCATED( pres_drag_x_av ) ) THEN

                DO  m = 1, surf_def%ns
                   i = surf_def%i(m)
                   j = surf_def%j(m)
                   k = surf_def%k(m)
!
!--                Pressure drag on right (east) faces
                   pres_drag_x_av(j,i) = pres_drag_x_av(j,i)                                       &
                                   + MERGE( -p(k,j,i) * dy * dz(1), 0.0_wp, surf_def%eastward(m) )
!
!--                Pressure drag on left (west) faces
                   pres_drag_x_av(j,i) = pres_drag_x_av(j,i)                                       &
                                   + MERGE(  p(k,j,i) * dy * dz(1), 0.0_wp, surf_def%westward(m) )
                ENDDO

                DO  m = 1, surf_lsm%ns
                   i = surf_lsm%i(m)
                   j = surf_lsm%j(m)
                   k = surf_lsm%k(m)
!
!--                Pressure drag on right (east) faces
                   pres_drag_x_av(j,i) = pres_drag_x_av(j,i)                                       &
                                   + MERGE( -p(k,j,i) * dy * dz(1), 0.0_wp, surf_lsm%eastward(m) )
!
!--                Pressure drag on left (west) faces
                   pres_drag_x_av(j,i) = pres_drag_x_av(j,i)                                       &
                                   + MERGE(  p(k,j,i) * dy * dz(1), 0.0_wp, surf_lsm%westward(m) )
                ENDDO
                DO  m = 1, surf_usm%ns
                   i = surf_usm%i(m)
                   j = surf_usm%j(m)
                   k = surf_usm%k(m)
!
!--                Pressure drag on right (east) faces
                   pres_drag_x_av(j,i) = pres_drag_x_av(j,i)                                       &
                                   + MERGE( -p(k,j,i) * dy * dz(1), 0.0_wp, surf_usm%eastward(m) )
!
!--                Pressure drag on left (west) faces
                   pres_drag_x_av(j,i) = pres_drag_x_av(j,i)                                       &
                                   + MERGE(  p(k,j,i) * dy * dz(1), 0.0_wp, surf_usm%westward(m) )
                ENDDO

             ENDIF

          CASE ( 'pres_drag_y*' )
             IF ( ALLOCATED( pres_drag_y_av ) ) THEN

               DO  m = 1, surf_def%ns
                   i = surf_def%i(m)
                   j = surf_def%j(m)
                   k = surf_def%k(m)
!
!--                Pressure drag on north faces
                   pres_drag_y_av(j,i) = pres_drag_y_av(j,i)                                       &
                                   + MERGE( -p(k,j,i) * dx * dz(1), 0.0_wp, surf_def%northward(m) )
!
!--                Pressure drag on south faces
                   pres_drag_y_av(j,i) = pres_drag_y_av(j,i)                                       &
                                   + MERGE(  p(k,j,i) * dx * dz(1), 0.0_wp, surf_def%southward(m) )
                ENDDO
                DO  m = 1, surf_lsm%ns
                   i = surf_lsm%i(m)
                   j = surf_lsm%j(m)
                   k = surf_lsm%k(m)
!
!--                Pressure drag on north faces
                   pres_drag_y_av(j,i) = pres_drag_y_av(j,i)                                       &
                                    + MERGE( -p(k,j,i) * dx * dz(1), 0.0_wp,surf_lsm%northward(m) )
!
!--                Pressure drag on south faces
                   pres_drag_y_av(j,i) = pres_drag_y_av(j,i)                                       &
                                    + MERGE(  p(k,j,i) * dx * dz(1), 0.0_wp, surf_lsm%southward(m) )
                ENDDO
                DO  m = 1, surf_usm%ns
                   i = surf_usm%i(m)
                   j = surf_usm%j(m)
                   k = surf_usm%k(m)
!
!--                Pressure drag on north faces
                   pres_drag_y_av(j,i) = pres_drag_y_av(j,i)                                       &
                                    + MERGE( -p(k,j,i) * dx * dz(1), 0.0_wp, surf_usm%northward(m) )
!
!--                Pressure drag on south faces
                   pres_drag_y_av(j,i) = pres_drag_y_av(j,i)                                       &
                                    + MERGE(  p(k,j,i) * dx * dz(1), 0.0_wp, surf_usm%southward(m) )
                ENDDO

             ENDIF

          CASE ( 'q' )
             IF ( ALLOCATED( q_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         q_av(k,j,i) = q_av(k,j,i) + q(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql' )
             IF ( ALLOCATED( ql_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_av(k,j,i) = ql_av(k,j,i) + ql(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_c' )
             IF ( ALLOCATED( ql_c_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_c_av(k,j,i) = ql_c_av(k,j,i) + ql_c(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_v' )
             IF ( ALLOCATED( ql_v_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         ql_v_av(k,j,i) = ql_v_av(k,j,i) + ql_v(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ql_vp' )
             IF ( ALLOCATED( ql_vp_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         number_of_particles = prt_count(k,j,i)
                         IF ( number_of_particles <= 0 )  CYCLE
                         particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                         DO  n = 1, number_of_particles
                            IF ( particles(n)%particle_mask )  THEN
                               ql_vp_av(k,j,i) = ql_vp_av(k,j,i) +                                 &
                                                 particles(n)%weight_factor / number_of_particles
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qsurf*' )
             IF ( ALLOCATED( qsurf_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  THEN
                            qsurf_av(j,i) = qsurf_av(j,i) + surf_def%q_surface(m)
                         ENDIF
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            qsurf_av(j,i) = qsurf_av(j,i) + surf_lsm%q_surface(m)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            qsurf_av(j,i) = qsurf_av(j,i) + surf_usm%q_surface(m)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qsws*' )
!
!--          In case of default surfaces, clean-up flux by density.
!--          In case of land- and urban-surfaces, convert fluxes into dynamic units.
!--          Question (maronga): are the .NOT. statements really required?
             IF ( ALLOCATED( qsws_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  THEN
                            k = surf_def%k(m)
                            qsws_av(j,i) = qsws_av(j,i) + surf_def%qsws(m) *                       &
                                                          waterflux_output_conversion(k)
                         ENDIF
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            k = surf_lsm%k(m)
                            qsws_av(j,i) = qsws_av(j,i) + surf_lsm%qsws(m) *                       &
                                                          waterflux_output_conversion(k)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            k = surf_usm%k(m)
                            qsws_av(j,i) = qsws_av(j,i) + surf_usm%qsws(m) *                       &
                                                          waterflux_output_conversion(k)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qv' )
             IF ( ALLOCATED( qv_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         qv_av(k,j,i) = qv_av(k,j,i) + q(k,j,i) - ql(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'r_a*' )
             IF ( ALLOCATED( r_a_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  r_a_av(j,i) = r_a_av(j,i) + surf_lsm%r_a(m)
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            r_a_av(j,i) = r_a_av(j,i) +                                            &
                                       surf_usm%frac(m,ind_veg_wall)  * surf_usm%r_a(m)       +    &
                                       surf_usm%frac(m,ind_pav_green) * surf_usm%r_a_green(m) +    &
                                       surf_usm%frac(m,ind_wat_win)   * surf_usm%r_a_window(m)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 's' )
             IF ( ALLOCATED( s_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         s_av(k,j,i) = s_av(k,j,i) + s(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'shf*' )
!
!--          In case of default surfaces, clean-up flux by density.
!--          In case of land- and urban-surfaces, convert fluxes into dynamic units.
             IF ( ALLOCATED( shf_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  THEN
                            k = surf_def%k(m)
                            shf_av(j,i) = shf_av(j,i) + surf_def%shf(m) *                          &
                                                        heatflux_output_conversion(k)
                         ENDIF
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            k = surf_lsm%k(m)
                            shf_av(j,i) = shf_av(j,i) + surf_lsm%shf(m) *                          &
                                                        heatflux_output_conversion(k)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            k = surf_usm%k(m)
                            shf_av(j,i) = shf_av(j,i) + surf_usm%shf(m) *                          &
                                                        heatflux_output_conversion(k)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ssurf*' )
             IF ( ALLOCATED( ssurf_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      k = topo_top_ind(j,i,0)
                      ssurf_av(j,i) = ssurf_av(j,i) + s(k,j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ssws*' )
             IF ( ALLOCATED( ssws_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  THEN
                            k = surf_def%k(m)
                            ssws_av(j,i) = ssws_av(j,i) + surf_def%ssws(m) *                       &
                                                          scalarflux_output_conversion(k)
                         ENDIF
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            k =  surf_lsm%k(m)
                            ssws_av(j,i) = ssws_av(j,i) + surf_lsm%ssws(m) *                       &
                                                          scalarflux_output_conversion(k)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            k = surf_usm%k(m)
                            ssws_av(j,i) = ssws_av(j,i) + surf_usm%ssws(m) *                       &
                                                          scalarflux_output_conversion(k)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 't*' )
             IF ( ALLOCATED( ts_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  ts_av(j,i) = ts_av(j,i) + surf_def%ts(m)
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  ts_av(j,i) = ts_av(j,i) + surf_lsm%ts(m)
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  ts_av(j,i) = ts_av(j,i) + surf_usm%ts(m)
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'theta' )
             IF ( ALLOCATED( pt_av ) ) THEN
                IF ( .NOT. bulk_cloud_model ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                            pt_av(k,j,i) = pt_av(k,j,i) + pt(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                            pt_av(k,j,i) = pt_av(k,j,i) + pt(k,j,i) + lv_d_cp * d_exner(k)         &
                                                                              * ql(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF

          CASE ( 'tsurf*' )
             IF ( ALLOCATED( tsurf_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  THEN
                            tsurf_av(j,i) = tsurf_av(j,i) + surf_def%pt_surface(m)
                         ENDIF
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            tsurf_av(j,i) = tsurf_av(j,i) + surf_lsm%pt_surface(m)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            tsurf_av(j,i) = tsurf_av(j,i) + surf_usm%pt_surface(m)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'u' )
             IF ( ALLOCATED( u_av ) ) THEN
                IF ( interpolate_to_grid_center )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nzt+1
                            u_av(k,j,i) = u_av(k,j,i) + 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            u_av(k,j,i) = u_av(k,j,i) + u(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF

          CASE ( 'us*' )
             IF ( ALLOCATED( us_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  us_av(j,i) = us_av(j,i) + surf_def%us(m)
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  us_av(j,i) = us_av(j,i) + surf_lsm%us(m)
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  us_av(j,i) = us_av(j,i) + surf_usm%us(m)
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'v' )
             IF ( ALLOCATED( v_av ) ) THEN
                IF ( interpolate_to_grid_center )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         DO  k = nzb, nzt+1
                            v_av(k,j,i) = v_av(k,j,i) + 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            v_av(k,j,i) = v_av(k,j,i) + v(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF

          CASE ( 'thetav' )
             IF ( ALLOCATED( vpt_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         vpt_av(k,j,i) = vpt_av(k,j,i) + vpt(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'w' )
             IF ( ALLOCATED( w_av ) ) THEN
                IF ( interpolate_to_grid_center )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         w_av(nzb,j,i) = w_av(nzb,j,i) + w(nzb,j,i)
                         DO  k = nzb, nzt+1
                            kl = MAX( k-1, nzb )
                            w_av(k,j,i) = w_av(k,j,i) + 0.5_wp * ( w(k,j,i) + w(kl,j,i) )
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            w_av(k,j,i) = w_av(k,j,i) + w(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF

          CASE ( 'z0*' )
             IF ( ALLOCATED( z0_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  z0_av(j,i) = z0_av(j,i) + surf_def%z0(m)
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  z0_av(j,i) = z0_av(j,i) + surf_lsm%z0(m)
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  z0_av(j,i) = z0_av(j,i) + surf_usm%z0(m)
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'z0h*' )
             IF ( ALLOCATED( z0h_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  z0h_av(j,i) = z0h_av(j,i) + surf_def%z0h(m)
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  z0h_av(j,i) = z0h_av(j,i) + surf_lsm%z0h(m)
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  z0h_av(j,i) = z0h_av(j,i) + surf_usm%z0h(m)
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'z0q*' )
             IF ( ALLOCATED( z0q_av ) ) THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn


                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward_top(m) )  z0q_av(j,i) = z0q_av(j,i) + surf_def%z0q(m)
                      ENDDO

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  z0q_av(j,i) = z0q_av(j,i) + surf_lsm%z0q(m)
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  z0q_av(j,i) = z0q_av(j,i) + surf_usm%z0q(m)
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT

!--          In case of urban surface variables it should be always checked if respective arrays are
!--          allocated, at least in case of a restart run, as averaged usm arrays are not read from
!--          file at the moment.
             IF ( urban_surface )  THEN
                CALL usm_3d_data_averaging( 'allocate', trimvar )
             ENDIF

!
!--          Summing up data from all other modules
             CALL module_interface_3d_data_averaging( 'sum', trimvar )


       END SELECT

    ENDDO

    CALL cpu_log( log_point(34), 'sum_up_3d_data', 'stop' )


 END SUBROUTINE sum_up_3d_data
