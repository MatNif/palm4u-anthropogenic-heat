!> @file diagnostic_output_quantities_mod.f90
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
! Authors:
! --------
! @author Farah Kanani-Suehring
!
!
! Description:
! ------------
!> Calculation of diagnostic quantities like wind speed and direction, absolute air temperature,
!> relative humidity, etc.
!--------------------------------------------------------------------------------------------------!
 MODULE diagnostic_output_quantities_mod

    USE arrays_3d,                                                                                 &
        ONLY:  d_exner,                                                                            &
               ddzu,                                                                               &
               ddzw,                                                                               &
               dd2zu,                                                                              &
               dzw,                                                                                &
               exner,                                                                              &
               hyp,                                                                                &
               km,                                                                                 &
               pt,                                                                                 &
               q,                                                                                  &
               ql,                                                                                 &
               s,                                                                                  &
               u,                                                                                  &
               v,                                                                                  &
               w,                                                                                  &
               zu,                                                                                 &
               zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  barometric_formula,                                                                 &
               degc_to_k,                                                                          &
               exner_function,                                                                     &
               kappa,                                                                              &
               pi,                                                                                 &
               magnus,                                                                             &
               rd_d_rv

    USE boundary_settings_mod,                                                                     &
        ONLY:  set_lateral_neumann_bc

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model

    USE control_parameters,                                                                        &
        ONLY:  average_count_3d,                                                                   &
               current_timestep_number,                                                            &
               cyclic_fill_initialization,                                                         &
               data_output,                                                                        &
               data_output_pr,                                                                     &
               do2d,                                                                               &
               do3d,                                                                               &
               domask,                                                                             &
               humidity,                                                                           &
               interpolate_to_grid_center,                                                         &
               length,                                                                             &
               masks,                                                                              &
               mask_i,                                                                             &
               mask_j,                                                                             &
               mask_k,                                                                             &
               mask_size_l,                                                                        &
               mask_surface,                                                                       &
               message_string,                                                                     &
               neutral,                                                                            &
               passive_scalar,                                                                     &
               plant_canopy,                                                                       &
               pt_surface,                                                                         &
               surface_pressure,                                                                   &
               restart_data_format_output,                                                         &
               restart_string,                                                                     &
               varnamelength

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz_2d

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nzb,                                                                                &
               nzb_max,                                                                            &
               nzt,                                                                                &
               topo_flags,                                                                         &
               topo_top_ind

    USE kinds

    USE palm_date_time_mod,                                                                        &
        ONLY:  seconds_per_hour

    USE pegrid

    USE plant_canopy_model_mod,                                                                    &
        ONLY: pch_index_ji,                                                                        &
              pcm_latentrate,                                                                      &
              pcm_sensiblerate

    USE profil_parameter,                                                                          &
        ONLY:  dopr_index

    USE radiation_model_mod,                                                                       &
        ONLY:  rad_lw_hr,                                                                          &
               radiation_scheme,                                                                   &
               rad_sw_hr

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array,                                                              &
               rrd_mpi_io,                                                                         &
               wrd_mpi_io

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               statistic_regions

    USE surface_layer_fluxes_mod,                                                                  &
        ONLY:  psi_m,                                                                              &
               psi_h

    USE surface_mod,                                                                               &
        ONLY:  surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_type,                                                                          &
               surf_usm

    IMPLICIT NONE

    CHARACTER(LEN=varnamelength), DIMENSION(500) ::  do_all = ' '

    INTEGER(iwp) ::  timestep_number_at_prev_calc = 0  !< ...at previous diagnostic output calculation

    LOGICAL ::  initialized_diagnostic_output_quantities = .FALSE. !< flag indicating whether output is initialized
    LOGICAL ::  prepared_diagnostic_output_quantities = .FALSE.    !< flag indicating whether output is p

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  count_kfd  !< average counter for katabatic-flow depth

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  kfd         !< katabatic-flow depth
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  kfd_av      !< time-averaged katabatic-flow depth
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  slope_angle !< terrain slope angle
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pt_2m       !< 2-m air potential temperature
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pt_2m_av    !< averaged 2-m air potential temperature
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  qv_2m       !< 2-m water vapor mixing ratio
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  qv_2m_av    !< averaged 2-m water vapor mixing ratio
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ta_2m       !< 2-m air temperature
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ta_2m_av    !< averaged 2-m air temperature
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_10m      !< horizontal wind speed at 10m
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  uv_10m_av   !< averaged horizontal wind speed at 10m
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_25m      !< volume flux integrated up to 25m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_25m_av   !< time-averaged volume flux integrated up to 25m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_50m      !< volume flux integrated up to 50m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_50m_av   !< time-averaged volume flux integrated up to 50m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_75m      !< volume flux integrated up to 75m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_75m_av   !< time-averaged volume flux integrated up to 75m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_100m     !< volume flux integrated up to 100m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_100m_av  !< time-averaged volume flux integrated up to 100m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_xxm      !< volume flux integrated up to katabatic-flow depth (variable height)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vf_xxm_av   !< time-averaged volume flux integrated up to katabatic-flow depth (variable height)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_25m     !< volume flux density integrated up to 25m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_25m_av  !< time-averaged volume-flux density integrated up to 25m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_50m     !< volume flux density integrated up to 50m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_50m_av  !< time-averaged volume-flux density  integrated up to 50m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_75m     !< volume flux density integrated up to 75m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_75m_av  !< time-averaged volume-flux density integrated up to 75m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_100m    !< volume flux density integrated up to 100m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_100m_av !< time-averaged volume-flux density integrated up to 100m above ground
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_xxm     !< volume flux density integrated up to katabatic-flow depth (variable height)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vfd_xxm_av  !< time-averaged volume-flux densitiy integrated up to katabatic-flow depth

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  div_new      !< divergence after calling pressure solver
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  div_new_av   !< time-averaged divergence after calling pressure solver
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  div_old      !< divergence before calling pressure solver
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  div_old_av   !< time-averaged divergence before calling pressure solver
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  hr           !< heating rate
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  hr_av        !< time-averaged heating rate
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  rh           !< relative humidity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  rh_av        !< avg. relative humidity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ta           !< air temperature
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ta_av        !< avg. air temperature
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ti           !< rotation(u,v,w) aka turbulence intensity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ti_av        !< avg. rotation(u,v,w) aka turbulence intensity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u_center     !< u at center of grid box
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u_center_av  !< mean of u_center
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uu           !< uu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uu_av        !< mean of uu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uv           !< uv
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uv_av        !< mean of uv
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uw           !< uw
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uw_av        !< mean of uw
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v_center     !< v at center of grid box
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v_center_av  !< mean of v_center
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vu           !< vu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vu_av        !< mean of vu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vv           !< vv
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vv_av        !< mean of vv
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vw           !< vw
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vw_av        !< mean of vw
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wdir         !< wind direction
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wdir_av      !< mean wind direction
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wq           !< wq
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wq_av        !< mean of wq
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ws           !< ws
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ws_av        !< mean of ws
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wspeed       !< horizontal wind speed
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wspeed_av    !< mean of horizotal wind speed
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wtheta       !< wtheta
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wtheta_av    !< mean of wtheta
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ww           !< ww
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ww_av        !< mean of ww
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wu           !< wu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wu_av        !< mean of wu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wv           !< wv
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wv_av        !< mean of wv

    SAVE

    PRIVATE

!
!-- Public variables
    PUBLIC div_new,                                                                                &
           div_old,                                                                                &
           initialized_diagnostic_output_quantities,                                               &
           prepared_diagnostic_output_quantities,                                                  &
           timestep_number_at_prev_calc
!
!-- Public routines
    PUBLIC doq_3d_data_averaging,                                                                  &
           doq_actions,                                                                            &
           doq_calculate,                                                                          &
           doq_check_data_output,                                                                  &
           doq_check_data_output_pr,                                                               &
           doq_define_netcdf_grid,                                                                 &
           doq_init,                                                                               &
           doq_output_2d,                                                                          &
           doq_output_3d,                                                                          &
           doq_output_mask,                                                                        &
           doq_statistics,                                                                         &
           doq_rrd_local,                                                                          &
           doq_wrd_local


    INTERFACE doq_3d_data_averaging
       MODULE PROCEDURE doq_3d_data_averaging
    END INTERFACE doq_3d_data_averaging

    INTERFACE doq_actions
       MODULE PROCEDURE doq_actions
       MODULE PROCEDURE doq_actions_ij
    END INTERFACE doq_actions

    INTERFACE doq_calculate
       MODULE PROCEDURE doq_calculate
    END INTERFACE doq_calculate

    INTERFACE doq_check_data_output
       MODULE PROCEDURE doq_check_data_output
    END INTERFACE doq_check_data_output

    INTERFACE doq_check_data_output_pr
       MODULE PROCEDURE doq_check_data_output_pr
    END INTERFACE doq_check_data_output_pr

    INTERFACE doq_define_netcdf_grid
       MODULE PROCEDURE doq_define_netcdf_grid
    END INTERFACE doq_define_netcdf_grid

    INTERFACE doq_output_2d
       MODULE PROCEDURE doq_output_2d
    END INTERFACE doq_output_2d

    INTERFACE doq_output_3d
       MODULE PROCEDURE doq_output_3d
    END INTERFACE doq_output_3d

    INTERFACE doq_output_mask
       MODULE PROCEDURE doq_output_mask
    END INTERFACE doq_output_mask

    INTERFACE doq_init
       MODULE PROCEDURE doq_init
    END INTERFACE doq_init

    INTERFACE doq_statistics
       MODULE PROCEDURE doq_statistics
    END INTERFACE doq_statistics

    INTERFACE doq_prepare
       MODULE PROCEDURE doq_prepare
    END INTERFACE doq_prepare

   INTERFACE doq_rrd_local
       MODULE PROCEDURE doq_rrd_local_ftn
       MODULE PROCEDURE doq_rrd_local_mpi
    END INTERFACE doq_rrd_local

    INTERFACE doq_wrd_local
       MODULE PROCEDURE doq_wrd_local
    END INTERFACE doq_wrd_local

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average diagnostic output quantities as well as allocate the array necessary for
!> storing the average.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_3d_data_averaging( mode, variable )

    CHARACTER (LEN=*) ::  mode     !<
    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  i !<
    INTEGER(iwp) ::  j !<
    INTEGER(iwp) ::  k !<


    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'div_new' )
             IF ( .NOT. ALLOCATED( div_new_av ) )  THEN
                ALLOCATE( div_new_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             div_new_av = 0.0_wp

          CASE ( 'div_old' )
             IF ( .NOT. ALLOCATED( div_old_av ) )  THEN
                ALLOCATE( div_old_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             div_old_av = 0.0_wp

          CASE ( 'kfd*' )
             IF ( .NOT. ALLOCATED( kfd_av ) )  THEN
                ALLOCATE( kfd_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
!
!--          Allocate cound_kfd without ghost layers since this has not been implemented in the
!--          restart IO yet.
             IF ( .NOT. ALLOCATED( count_kfd ) )  THEN
                ALLOCATE( count_kfd(nys:nyn,nxl:nxr) )
             ENDIF
             kfd_av    = 0.0_wp
             count_kfd = 0

          CASE ( 'hr' )
             IF ( .NOT. ALLOCATED( hr_av ) )  THEN
                ALLOCATE( hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             hr_av = 0.0_wp

          CASE ( 'vf25m*' )
             IF ( .NOT. ALLOCATED( vf_25m_av ) )  THEN
                ALLOCATE( vf_25m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vf_25m_av = 0.0_wp

          CASE ( 'vf50m*' )
             IF ( .NOT. ALLOCATED( vf_50m_av ) )  THEN
                ALLOCATE( vf_50m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vf_50m_av = 0.0_wp

          CASE ( 'vf75m*' )
             IF ( .NOT. ALLOCATED( vf_75m_av ) )  THEN
                ALLOCATE( vf_75m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vf_75m_av = 0.0_wp

          CASE ( 'vf100m*' )
             IF ( .NOT. ALLOCATED( vf_100m_av ) )  THEN
                ALLOCATE( vf_100m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vf_100m_av = 0.0_wp

          CASE ( 'vfxxm*' )
             IF ( .NOT. ALLOCATED( vf_xxm_av  ) )  ALLOCATE( vf_xxm_av(nysg:nyng,nxlg:nxrg)  )
             vf_xxm_av  = 0.0_wp

          CASE ( 'vfd25m*' )
             IF ( .NOT. ALLOCATED( vfd_25m_av ) )  THEN
                ALLOCATE( vfd_25m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vfd_25m_av = 0.0_wp

          CASE ( 'vfd50m*' )
             IF ( .NOT. ALLOCATED( vfd_50m_av ) )  THEN
                ALLOCATE( vfd_50m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vfd_50m_av = 0.0_wp

          CASE ( 'vfd75m*' )
             IF ( .NOT. ALLOCATED( vfd_75m_av ) )  THEN
                ALLOCATE( vfd_75m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vfd_75m_av = 0.0_wp

          CASE ( 'vfd100m*' )
             IF ( .NOT. ALLOCATED( vfd_100m_av ) )  THEN
                ALLOCATE( vfd_100m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vfd_100m_av = 0.0_wp

          CASE ( 'vfdxxm*' )
             IF ( .NOT. ALLOCATED( vfd_xxm_av  ) )  THEN
                ALLOCATE( vfd_xxm_av(nysg:nyng,nxlg:nxrg)  )
             ENDIF
             vfd_xxm_av = 0.0_wp

          CASE ( 'rh' )
             IF ( .NOT. ALLOCATED( rh_av ) )  THEN
                ALLOCATE( rh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             rh_av = 0.0_wp

          CASE ( 'ta' )
             IF ( .NOT. ALLOCATED( ta_av ) )  THEN
                ALLOCATE( ta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             ta_av = 0.0_wp

          CASE ( 'ti' )
             IF ( .NOT. ALLOCATED( ti_av ) )  THEN
                ALLOCATE( ti_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             ti_av = 0.0_wp

          CASE ( 'uu_product' )
             IF ( .NOT. ALLOCATED( uu_av ) )  THEN
                ALLOCATE( uu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             uu_av = 0.0_wp

          CASE ( 'uv_product' )
             IF ( .NOT. ALLOCATED( uv_av ) )  THEN
                ALLOCATE( uv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             uv_av = 0.0_wp

          CASE ( 'uw_product' )
             IF ( .NOT. ALLOCATED( uw_av ) )  THEN
                ALLOCATE( uw_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             uw_av = 0.0_wp

          CASE ( 'vu_product' )
             IF ( .NOT. ALLOCATED( vu_av ) )  THEN
                ALLOCATE( vu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             vu_av = 0.0_wp

          CASE ( 'vv_product' )
             IF ( .NOT. ALLOCATED( vv_av ) )  THEN
                ALLOCATE( vv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             vv_av = 0.0_wp

          CASE ( 'vw_product' )
             IF ( .NOT. ALLOCATED( vw_av ) )  THEN
                ALLOCATE( vw_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             vw_av = 0.0_wp

          CASE ( 'wu_product' )
             IF ( .NOT. ALLOCATED( wu_av ) )  THEN
                ALLOCATE( wu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wu_av = 0.0_wp

           CASE ( 'wv_product' )
             IF ( .NOT. ALLOCATED( wv_av ) )  THEN
                ALLOCATE( wv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wv_av = 0.0_wp

          CASE ( 'ww_product' )
             IF ( .NOT. ALLOCATED( ww_av ) )  THEN
                ALLOCATE( ww_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             ww_av = 0.0_wp

           CASE ( 'wtheta_product' )
             IF ( .NOT. ALLOCATED( wtheta_av ) )  THEN
                ALLOCATE( wtheta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wtheta_av = 0.0_wp

           CASE ( 'wq_product' )
             IF ( .NOT. ALLOCATED( wq_av ) )  THEN
                ALLOCATE( wq_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wq_av = 0.0_wp

          CASE ( 'ws_product' )
             IF ( .NOT. ALLOCATED( ws_av ) )  THEN
                ALLOCATE( ws_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             ws_av = 0.0_wp

          CASE ( 'theta_2m*' )
             IF ( .NOT. ALLOCATED( pt_2m_av ) )  THEN
                ALLOCATE( pt_2m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             pt_2m_av = 0.0_wp

          CASE ( 'qv_2m*' )
             IF ( .NOT. ALLOCATED( qv_2m_av ) )  THEN
                ALLOCATE( qv_2m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             qv_2m_av = 0.0_wp

          CASE ( 'ta_2m*' )
             IF ( .NOT. ALLOCATED( ta_2m_av ) )  THEN
                ALLOCATE( ta_2m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             ta_2m_av = 0.0_wp

          CASE ( 'wspeed_10m*' )
             IF ( .NOT. ALLOCATED( uv_10m_av ) )  THEN
                ALLOCATE( uv_10m_av(nysg:nyng,nxlg:nxrg) )
             ENDIF
             uv_10m_av = 0.0_wp

          CASE ( 'wspeed' )
             IF ( .NOT. ALLOCATED( wspeed_av ) )  THEN
                ALLOCATE( wspeed_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             wspeed_av = 0.0_wp

          CASE ( 'wdir' )
             IF ( .NOT. ALLOCATED( u_center_av ) )  THEN
                ALLOCATE( u_center_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( .NOT. ALLOCATED( v_center_av ) )  THEN
                ALLOCATE( v_center_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             u_center_av = 0.0_wp
             v_center_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'div_new' )
             IF ( ALLOCATED( div_new_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         div_new_av(k,j,i) = div_new_av(k,j,i) + div_new(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'div_old' )
             IF ( ALLOCATED( div_old_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         div_old_av(k,j,i) = div_old_av(k,j,i) + div_old(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'kfd*' )
             IF ( ALLOCATED( kfd_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      kfd_av(j,i) = kfd_av(j,i) + kfd(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'hr' )
             IF ( ALLOCATED( hr_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         hr_av(k,j,i) = hr_av(k,j,i) + hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vf25m*' )
             IF ( ALLOCATED( vf_25m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_25m_av(j,i) = vf_25m_av(j,i) + vf_25m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vf50m*' )
             IF ( ALLOCATED( vf_50m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_50m_av(j,i) = vf_50m_av(j,i) + vf_50m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vf75m*' )
             IF ( ALLOCATED( vf_75m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_75m_av(j,i) = vf_75m_av(j,i) + vf_75m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vf100m*' )
             IF ( ALLOCATED( vf_100m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_100m_av(j,i) = vf_100m_av(j,i) + vf_100m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfxxm*' )
             IF ( ALLOCATED( vf_xxm_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_xxm_av(j,i) = vf_xxm_av(j,i) + vf_xxm(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfd25m*' )
             IF ( ALLOCATED( vfd_25m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_25m_av(j,i) = vfd_25m_av(j,i) + vfd_25m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfd50m*' )
             IF ( ALLOCATED( vfd_50m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_50m_av(j,i) = vfd_50m_av(j,i) + vfd_50m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfd75m*' )
             IF ( ALLOCATED( vfd_75m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_75m_av(j,i) = vfd_75m_av(j,i) + vfd_75m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfd100m*' )
             IF ( ALLOCATED( vfd_100m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_100m_av(j,i) = vfd_100m_av(j,i) + vfd_100m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfdxxm*' )
             IF ( ALLOCATED( vfd_xxm_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_xxm_av(j,i) = vfd_xxm_av(j,i) + vfd_xxm(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rh' )
             IF ( ALLOCATED( rh_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         rh_av(k,j,i) = rh_av(k,j,i) + rh(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ta' )
             IF ( ALLOCATED( ta_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ta_av(k,j,i) = ta_av(k,j,i) + ta(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ti' )
             IF ( ALLOCATED( ti_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ti_av(k,j,i) = ti_av(k,j,i) + ti(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uu_product' )
             IF ( ALLOCATED( uu_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         uu_av(k,j,i) = uu_av(k,j,i) + uu(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uv_product' )
             IF ( ALLOCATED( uv_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         uv_av(k,j,i) = uv_av(k,j,i) + uv(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uw_product' )
             IF ( ALLOCATED( uw_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         uw_av(k,j,i) = uw_av(k,j,i) + uw(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vu_product' )
             IF ( ALLOCATED( vu_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         vu_av(k,j,i) = vu_av(k,j,i) + vu(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vv_product' )
             IF ( ALLOCATED( vv_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         vv_av(k,j,i) = vv_av(k,j,i) + vv(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vw_product' )
             IF ( ALLOCATED( vw_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         vw_av(k,j,i) = vw_av(k,j,i) + vw(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wu_product' )
             IF ( ALLOCATED( wu_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wu_av(k,j,i) = wu_av(k,j,i) + wu(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wv_product' )
             IF ( ALLOCATED( wv_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wv_av(k,j,i) = wv_av(k,j,i) + wv(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ww_product' )
             IF ( ALLOCATED( ww_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ww_av(k,j,i) = ww_av(k,j,i) + ww(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wtheta_product' )
             IF ( ALLOCATED( wtheta_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wtheta_av(k,j,i) = wtheta_av(k,j,i) + wtheta(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wq_product' )
             IF ( ALLOCATED( wq_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wq_av(k,j,i) = wq_av(k,j,i) + wq(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ws_product' )
             IF ( ALLOCATED( ws_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ws_av(k,j,i) = ws_av(k,j,i) + ws(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'theta_2m*' )
             IF ( ALLOCATED( pt_2m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      pt_2m_av(j,i) = pt_2m_av(j,i) + pt_2m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'qv_2m*' )
             IF ( ALLOCATED( qv_2m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      qv_2m_av(j,i) = qv_2m_av(j,i) + qv_2m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ta_2m*' )
             IF ( ALLOCATED( ta_2m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      ta_2m_av(j,i) = ta_2m_av(j,i) + ta_2m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wspeed_10m*' )
             IF ( ALLOCATED( uv_10m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      uv_10m_av(j,i) = uv_10m_av(j,i) + uv_10m(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wspeed' )
            IF ( ALLOCATED( wspeed_av ) )  THEN
               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     DO  k = nzb, nzt+1
                         wspeed_av(k,j,i) = wspeed_av(k,j,i) + wspeed(k,j,i)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

          CASE ( 'wdir' )
             IF ( ALLOCATED( u_center_av )  .AND.  ALLOCATED( v_center_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                        u_center_av(k,j,i) = u_center_av(k,j,i) + u_center(k,j,i)
                        v_center_av(k,j,i) = v_center_av(k,j,i) + v_center(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'div_new' )
             IF ( ALLOCATED( div_new_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         div_new_av(k,j,i) = div_new_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'div_old' )
             IF ( ALLOCATED( div_old_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         div_old_av(k,j,i) = div_old_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'kfd*' )
             IF ( ALLOCATED( kfd_av ) )  THEN
!
!--             Average katabatic flow depth. As katabatic flows occur sometimes only occasionally,
!--             divide by the number of detected occasions that enter the average. The output will
!--             be a depth that describe the mean characteristics of the katabatic flow. If
!--             every timestep would be used but the katabatic flow might occur maybe only
!--             half the time, the mean flow depth would indicate deceptive flow characteristics.
!--             In order to prevent divisions by zero, limit the counter variable to 1.
                count_kfd = MERGE( count_kfd, 1, count_kfd /= 0 )
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      kfd_av(j,i) = kfd_av(j,i) / REAL( count_kfd(j,i), KIND=wp )
                   ENDDO
                ENDDO
             ENDIF
             count_kfd = 0

          CASE ( 'hr' )
             IF ( ALLOCATED( hr_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         hr_av(k,j,i) = hr_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vf25m*' )
             IF ( ALLOCATED( vf_25m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_25m_av(j,i) = vf_25m_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vf50m*' )
             IF ( ALLOCATED( vf_50m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_50m_av(j,i) = vf_50m_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vf75m*' )
             IF ( ALLOCATED( vf_75m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_75m_av(j,i) = vf_75m_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vf100m*' )
             IF ( ALLOCATED( vf_100m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_100m_av(j,i) = vf_100m_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF
!
!--       Volume flux integrated up to a detected depth. Note, in contrast to the katabatic flow
!--       depth, which is a characterist of the katabatic flow, this variables gives an estimate
!--       of the mean volume flux rate to estimate the fresh air supply in an area. If there is
!--       no katabatic jet detected, a zero value will enter the average. In contrast to the
!--       katabatic flow depth this is no pure katabatic flow characteristic but a derived
!--       measure so that summation and averaging goes as usual.
          CASE ( 'vfxxm*' )
             IF ( ALLOCATED( vf_xxm_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vf_xxm_av(j,i) = vf_xxm_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfd25m*' )
             IF ( ALLOCATED( vfd_25m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_25m_av(j,i) = vfd_25m_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfd50m*' )
             IF ( ALLOCATED( vfd_50m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_50m_av(j,i) = vfd_50m_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfd75m*' )
             IF ( ALLOCATED( vfd_75m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_75m_av(j,i) = vfd_75m_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vfd100m*' )
             IF ( ALLOCATED( vfd_100m_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_100m_av(j,i) = vfd_100m_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF
!
!--       Volume flux integrated up to a detected depth. Note, in contrast to the katabatic flow
!--       depth, which is a characterist of the katabatic flow, this variables gives an estimate
!--       of the mean volume flux rate to estimate the fresh air supply in an area. If there is
!--       no katabatic jet detected, a zero value will enter the average. In contrast to the
!--       katabatic flow depth this is no pure katabatic flow characteristic but a derived
!--       measure so that summation and averaging goes as usual.
          CASE ( 'vfdxxm*' )
             IF ( ALLOCATED( vfd_xxm_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      vfd_xxm_av(j,i) = vfd_xxm_av(j,i) / REAL( average_count_3d, KIND=wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rh' )
             IF ( ALLOCATED( rh_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         rh_av(k,j,i) = rh_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ta' )
             IF ( ALLOCATED( ta_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ta_av(k,j,i) = ta_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ti' )
             IF ( ALLOCATED( ti_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ti_av(k,j,i) = ti_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uu_product' )
             IF ( ALLOCATED( uu_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         uu_av(k,j,i) = uu_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uv_product' )
             IF ( ALLOCATED( uv_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         uv_av(k,j,i) = uv_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'uw_product' )
             IF ( ALLOCATED( uw_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         uw_av(k,j,i) = uw_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vu_product' )
             IF ( ALLOCATED( vu_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         vu_av(k,j,i) = vu_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vv_product' )
             IF ( ALLOCATED( vv_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         vv_av(k,j,i) = vv_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'vw_product' )
             IF ( ALLOCATED( vw_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         vw_av(k,j,i) = vw_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wu_product' )
             IF ( ALLOCATED( wu_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wu_av(k,j,i) = wu_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wv_product' )
             IF ( ALLOCATED( wv_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wv_av(k,j,i) = wv_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ww_product' )
             IF ( ALLOCATED( ww_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ww_av(k,j,i) = ww_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wtheta_product' )
             IF ( ALLOCATED( wtheta_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wtheta_av(k,j,i) = wtheta_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wq_product' )
             IF ( ALLOCATED( wq_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wq_av(k,j,i) = wq_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'ws_product' )
             IF ( ALLOCATED( ws_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         ws_av(k,j,i) = ws_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

         CASE ( 'theta_2m*' )
            IF ( ALLOCATED( pt_2m_av ) )  THEN
               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     pt_2m_av(j,i) = pt_2m_av(j,i) / REAL( average_count_3d, KIND=wp )
                  ENDDO
               ENDDO
            ENDIF

         CASE ( 'qv_2m*' )
            IF ( ALLOCATED( qv_2m_av ) )  THEN
               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     qv_2m_av(j,i) = qv_2m_av(j,i) / REAL( average_count_3d, KIND=wp )
                  ENDDO
               ENDDO
            ENDIF

         CASE ( 'ta_2m*' )
            IF ( ALLOCATED( ta_2m_av ) )  THEN
               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     ta_2m_av(j,i) = ta_2m_av(j,i) / REAL( average_count_3d, KIND=wp )
                  ENDDO
               ENDDO
            ENDIF

         CASE ( 'wspeed_10m*' )
            IF ( ALLOCATED( uv_10m_av ) )  THEN
               DO  i = nxlg, nxrg
                  DO  j = nysg, nyng
                     uv_10m_av(j,i) = uv_10m_av(j,i) / REAL( average_count_3d, KIND=wp )
                  ENDDO
               ENDDO
            ENDIF

         CASE ( 'wspeed' )
             IF ( ALLOCATED( wspeed_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         wspeed_av(k,j,i) = wspeed_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wdir' )
             IF ( ALLOCATED( u_center_av )  .AND.  ALLOCATED( v_center_av ) )  THEN

                IF ( .NOT. ALLOCATED( wdir_av ) )  THEN
                   ALLOCATE( wdir_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                wdir_av = 0.0_wp

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         u_center_av(k,j,i) = u_center_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         v_center_av(k,j,i) = v_center_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                         wdir_av(k,j,i) = ATAN2( u_center_av(k,j,i), v_center_av(k,j,i) )          &
                                          / pi * 180.0_wp + 180.0_wp
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

       END SELECT

    ENDIF

 END SUBROUTINE doq_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Specific actions before, during, and after time-integration.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_actions( location )

    CHARACTER(LEN=*), INTENT(IN) ::  location  !< call location string

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

       CASE ( 'after_integration' )

       CASE DEFAULT

          CONTINUE

    END SELECT

 END SUBROUTINE doq_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Specific actions during time-integration - call for each i,j.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_actions_ij( i, j, location )

    CHARACTER(LEN=*), INTENT(IN) ::  location  !< call location string

    INTEGER(iwp)             ::  dummy     !< dummy value to avoid 'un-used variable' compiler error
    INTEGER(iwp), INTENT(IN) ::  i         !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  j         !< grid index in y-direction

    dummy = i + j

    SELECT CASE ( location )

       CASE DEFAULT

          CONTINUE

    END SELECT

 END SUBROUTINE doq_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for diagnostic output
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_check_data_output( var, unit, i, ilen, k )

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit  !<
    CHARACTER (LEN=*) ::  var   !<

    INTEGER(iwp), OPTIONAL, INTENT(IN) ::  i     !< Current element of data_output
    INTEGER(iwp), OPTIONAL, INTENT(IN) ::  ilen  !< Length of current entry in data_output
    INTEGER(iwp), OPTIONAL, INTENT(IN) ::  k     !< Output is xy mode? 0 = no, 1 = yes


    SELECT CASE ( TRIM( var ) )

       CASE ( 'div_new' )
          unit = '1/s'

       CASE ( 'div_old' )
          unit = '1/s'

       CASE ( 'hr' )
          IF ( neutral )  THEN
             message_string = 'data_output = ' // TRIM( var ) //                                   &
                              ' is not implemented for neutral = .TRUE.'
             CALL message( 'diagnostic_output', 'PAC0194', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'K/h'

       CASE ( 'rh' )
          unit = '%'

       CASE ( 'ta' )
          unit = 'degree_C'

       CASE ( 'ti' )
          unit = '1/s'

       CASE ( 'uu_product', 'uv_product', 'uw_product', 'vu_product', 'vv_product', 'vw_product',  &
              'wu_product', 'wv_product', 'ww_product' )
          unit = 'm2/s2'

       CASE ( 'wtheta_product' )
          unit = 'Km/s'

       CASE ( 'wq_product' )
          unit = 'm/s'

       CASE ( 'ws_product' )
          IF ( .NOT.  passive_scalar )  THEN
             message_string = 'data_output = ' // TRIM( var ) //                                   &
                              ' is not implemented for passive_scalar = .FALSE.'
             CALL message( 'diagnostic_output', 'PAC0195', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'm/s'

       CASE ( 'wspeed' )
          unit = 'm/s'

       CASE ( 'wdir' )
          unit = 'degree'
!
!--    Treat horizotal cross-section output quantities.
       CASE ( 'kfd*', 'qv_2m*', 'theta_2m*', 'ta_2m*', 'vf25m*', 'vf50m*', 'vf75m*', 'vf100m*',    &
              'vfxxm*', 'vfd25m*', 'vfd50m*', 'vfd75m*', 'vfd100m*', 'vfdxxm*', 'wspeed_10m*' )
!
          unit = 'm3/m/s'
!--       Check if output quantity is _xy only.
          IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
             message_string = 'illegal value for data_output: "' //                                &
                              TRIM( var ) // '" & only 2d-horizontal ' //                          &
                              'cross sections are allowed for this value'
             CALL message( 'diagnostic_output', 'PAC0127', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( TRIM( var ) == 'kfd*'        )  unit = 'm'
          IF ( TRIM( var ) == 'qv_2m*'      )  unit = 'kg/kg'
          IF ( TRIM( var ) == 'theta_2m*'   )  unit = 'K'
          IF ( TRIM( var ) == 'ta_2m*'      )  unit = 'degree_C'
          IF ( TRIM( var ) == 'vf25m*'      )  unit = 'm3/s'
          IF ( TRIM( var ) == 'vf50m*'      )  unit = 'm3/s'
          IF ( TRIM( var ) == 'vf75m*'      )  unit = 'm3/s'
          IF ( TRIM( var ) == 'vf100m*'     )  unit = 'm3/s'
          IF ( TRIM( var ) == 'vfxxm*'      )  unit = 'm3/s'
          IF ( TRIM( var ) == 'vfd25m*'     )  unit = 'm3/m/s'
          IF ( TRIM( var ) == 'vfd50m*'     )  unit = 'm3/m/s'
          IF ( TRIM( var ) == 'vfd75m*'     )  unit = 'm3/m/s'
          IF ( TRIM( var ) == 'vfd100m*'    )  unit = 'm3/m/s'
          IF ( TRIM( var ) == 'vfdxxm*'     )  unit = 'm3/m/s'
          IF ( TRIM( var ) == 'wspeed_10m*' )  unit = 'm/s'


       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE doq_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of user defined profile output quantities. For those variables not recognized by the
!> user, the parameter unit is set to "illegal", which tells the calling routine that the
!> output variable is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_check_data_output_pr( variable, var_count, unit, dopr_unit )

    CHARACTER (LEN=*) ::  unit      !<
    CHARACTER (LEN=*) ::  variable  !<
    CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit

    INTEGER(iwp) ::  pr_index       !<
    INTEGER(iwp) ::  var_count      !<

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'rh' )
          IF ( .NOT.  humidity )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                              ' requires humidity'
             CALL message( 'doq_check_data_output_pr', 'PAC0196', 1, 2, 0, 6, 0 )
          ENDIF
          pr_index = 130
          dopr_index(var_count) = pr_index
          dopr_unit     = '%'
          unit = dopr_unit
          hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE doq_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Compute slope angle of terrain. This is used as critirion to detect katabatic flows. Small-scale
!> terrain features are filtered by computing spatially running averages over topography height
!> until an averaging radius is reached.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_compute_slope_angle

    INTEGER(iwp) ::  i      !< grid index in x-direction
    INTEGER(iwp) ::  j      !< grid index in y-direction
    INTEGER(iwp) ::  k      !< grid index in z-direction
    INTEGER(iwp) ::  n1     !< running index in x
    INTEGER(iwp) ::  n2     !< running index in y
    INTEGER(iwp) ::  np = 1 !< grid point stencil used to compute moving average in an iteration

    REAL(wp) ::  filter_radius  !< filter radius for topography height
    REAL(wp) ::  view_distance  !< current view distance in the filtering procedure
    REAL(wp) ::  norm           !< absolute value of the normal vector
    REAL(wp) ::  normal_x       !< normal vector component in x
    REAL(wp) ::  normal_y       !< normal vector component in y
    REAL(wp) ::  normal_z       !< normal vector component in z

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  topo_runave     !< filtered topography height
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  topo_runave_dum !< filtered topography height
!
!-- If slope angle has not been allocated, allocate the array and compute its values.
    IF ( .NOT. ALLOCATED( slope_angle ) )  THEN
       ALLOCATE( slope_angle(nysg:nyng,nxlg:nxrg) )
!
!-- If slope angle has been already allocated, everything is done and no further action is required.
    ELSE
       RETURN
    ENDIF

    IF ( .NOT. ALLOCATED( topo_runave     ) )  ALLOCATE( topo_runave(nysg:nyng,nxlg:nxrg)     )
    IF ( .NOT. ALLOCATED( topo_runave_dum ) )  ALLOCATE( topo_runave_dum(nysg:nyng,nxlg:nxrg) )
!
!-- Obtain the highest terrain-grid point at (i,j)-index.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt+1
             IF ( BTEST( topo_flags(k,j,i), 5 ) )  THEN
                topo_runave(j,i) = zw(k)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL exchange_horiz_2d( topo_runave )
    CALL set_lateral_neumann_bc( topo_runave )
!
!-- Smooth topography by a weighted moving average algorithm. A moving average is employed in order
!-- to smooth-out small-scale terrain features. A view distance of 50m seems to be a reasonable
!-- trade-off between filtering small-scale features and maintaining terrain features. To compute
!-- the moving average, adjoining grid points need to be accessed. At the subdomain boundaries,
!-- however, this means that ghost points need to be accessed. Depending on the employed advection
!-- scheme, either 3 (WS-scheme) or only 1 ghost point (PW-scheme, BC-scheme) are defined. For this
!-- reason, compute the moving average in an iterative approach. During the first iteration the
!-- terrain information within a defined stencil enters the moving average. After each iteration
!-- a ghost point exchange is carried out for the smoothed terrain height. In a second iteration,
!-- smoothed terrain information within the defined stencil enters the moving average. Iterations
!-- are repeated until a view distance has been reached (number of accessed grid points in a
!-- direction x horizontal grid spacing x number of iterations). This iterative approach creates
!-- a weighted moving average, putting the largest weight on the local terrain height but smoothing
!-- small-scale features. A further note: Since the number of ghost points depend on the advection
!-- scheme (3 for WS-scheme and 1 for PW-scheme), and thus also the possible grid-point stencil
!-- would depend on the advection scheme, the weighted moving average would become a function of
!-- the number of available ghost points. In order to avoid this, only define a 9-point stencil for
!-- each iteration (only view 1 grid point in each direction), independent on the maximum available
!-- number of ghost points.
    filter_radius = 50.0_wp

    view_distance = 0.0_wp
    DO WHILE ( view_distance < filter_radius )
       topo_runave_dum = topo_runave
       topo_runave = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  n1 = -np, np
                DO  n2 = -np, np
                   topo_runave(j,i) = topo_runave(j,i) + topo_runave_dum(j+n2,i+n1)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       topo_runave = topo_runave / ( ( 2.0_wp * np + 1 ) + ( 2.0_wp * np + 1 ) )
       CALL exchange_horiz_2d( topo_runave )
       CALL set_lateral_neumann_bc( topo_runave )

       view_distance = view_distance + SQRT( ( np * dx )**2 + ( np * dy )**2 )
    ENDDO

    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Based on the filtered terrain, now compute the normal vector components from topography
!--       gradients. For the vertical component presume that there is no downward-facing terrain.
          normal_x = ( topo_runave(j,i+1) - topo_runave(j,i-1) ) * ddx * 0.5_wp
          normal_y = ( topo_runave(j+1,i) - topo_runave(j-1,i) ) * ddy * 0.5_wp
          normal_z = 1.0_wp
!
!--       Normalize normal vector
          norm = SQRT( normal_x**2 + normal_y**2 + normal_z**2 )
          normal_x = normal_x / norm
          normal_y = normal_y / norm
          normal_z = normal_z / norm
!
!--       From noralized normal vector now compute the slope angle. Note, slope angle is defined
!--       in the same coordinate system as the wind-direction, in order to make direct comparisons.
          slope_angle(j,i) = ATAN2( -normal_x, -normal_y ) / pi * 180.0_wp + 180.0_wp
       ENDDO
    ENDDO
    CALL exchange_horiz_2d( slope_angle )
    CALL set_lateral_neumann_bc( slope_angle )

    DEALLOCATE( topo_runave     )
    DEALLOCATE( topo_runave_dum )

 END SUBROUTINE doq_compute_slope_angle


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)  ::  variable    !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !<

    LOGICAL, INTENT(OUT)           ::  found       !<


    found  = .TRUE.

    SELECT CASE ( TRIM( variable ) )
!
!--    s grid
       CASE ( 'div_new', 'div_new_xy', 'div_new_xz', 'div_new_yz',                                 &
              'div_old', 'div_old_xy', 'div_old_xz', 'div_old_yz',                                 &
              'hr', 'hr_xy', 'hr_xz', 'hr_yz',                                                     &
              'rh', 'rh_xy', 'rh_xz', 'rh_yz',                                                     &
              'ta', 'ta_xy', 'ta_xz', 'ta_yz',                                                     &
              'ti', 'ti_xy', 'ti_xz', 'ti_yz',                                                     &
              'wspeed', 'wspeed_xy', 'wspeed_xz', 'wspeed_yz',                                     &
              'wdir', 'wdir_xy', 'wdir_xz', 'wdir_yz' )

          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'
!
!--    s grid surface variables
       CASE ( 'kfd*_xy', 'qv_2m*_xy', 'theta_2m*_xy', 'ta_2m*_xy', 'wspeed_10m*_xy', 'vf25m*_xy',  &
              'vf50m*_xy', 'vf75m*_xy', 'vf100m*_xy', 'vfxxm*_xy', 'vfd25m*_xy', 'vfd50m*_xy',     &
              'vfd75m*_xy', 'vfd100m*_xy', 'vfdxxm*_xy' )

          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'
!
!--    u grid
       CASE ( 'uu_product', 'uu_product_xy', 'uu_product_xz', 'uu_product_yz',                     &
              'uv_product', 'uv_product_xy', 'uv_product_xz', 'uv_product_yz',                     &
              'uw_product', 'uw_product_xy', 'uw_product_xz', 'uw_product_yz' )
!
!--       s grid
          IF ( interpolate_to_grid_center )  THEN
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'
!
!--       v grid
          ELSE
             grid_x = 'xu'
             grid_y = 'y'
             grid_z = 'zu'
          ENDIF
!
!--    v grid
       CASE ( 'vu_product', 'vu_product_xy', 'vu_product_xz', 'vu_product_yz',                     &
              'vv_product', 'vv_product_xy', 'vv_product_xz', 'vv_product_yz',                     &
              'vw_product', 'vw_product_xy', 'vw_product_xz', 'vw_product_yz' )
!
!--       s grid
          IF ( interpolate_to_grid_center )  THEN
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'
!
!--       v grid
          ELSE
             grid_x = 'x'
             grid_y = 'yv'
             grid_z = 'zu'
          ENDIF

       CASE ( 'wq_product', 'wq_product_xy', 'wq_product_xz', 'wq_product_yz',                     &
              'ws_product', 'ws_product_xy', 'ws_product_xz', 'ws_product_yz',                     &
              'wtheta_product', 'wtheta_product_xy', 'wtheta_product_xz', 'wtheta_product_yz',     &
              'wu_product', 'wu_product_xy', 'wu_product_xz', 'wu_product_yz',                     &
              'wv_product', 'wv_product_xy', 'wv_product_xz', 'wv_product_yz',                     &
              'ww_product', 'ww_product_xy', 'ww_product_xz', 'ww_product_yz' )
!
!--       s grid
          IF ( interpolate_to_grid_center )  THEN
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'
!
!--       w grid
          ELSE
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zw'
          ENDIF


       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT

 END SUBROUTINE doq_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 2D output variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  grid     !<
    CHARACTER (LEN=*) ::  mode     !<
    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av       !< value indicating averaged or non-averaged output
    INTEGER(iwp) ::  flag_nr  !< number of the topography flag (0: scalar, 1: u, 2: v, 3: w)
    INTEGER(iwp) ::  i        !< grid index x-direction
    INTEGER(iwp) ::  j        !< grid index y-direction
    INTEGER(iwp) ::  k        !< grid index z-direction
    INTEGER(iwp) ::  nzb_do   !<
    INTEGER(iwp) ::  nzt_do   !<

    LOGICAL ::  found             !< true if variable is in list
    LOGICAL ::  resorted          !< true if array is resorted
    LOGICAL ::  two_d             !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !<
    REAL(wp), DIMENSION(:,:,:), POINTER ::                 to_be_resorted  !< points to array which needs to be resorted for output


    flag_nr  = 0
    found    = .TRUE.
    resorted = .FALSE.
    two_d    = .FALSE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'div_new_xy', 'div_new_xz', 'div_new_yz' )
           IF ( av == 0 )  THEN
              to_be_resorted => div_new
           ELSE
              IF ( .NOT. ALLOCATED( div_new_av ) )  THEN
                 ALLOCATE( div_new_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                 div_new_av = 0.0_wp
              ENDIF
              to_be_resorted => div_new_av
           ENDIF
           flag_nr = 0

           IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'div_old_xy', 'div_old_xz', 'div_old_yz' )
           IF ( av == 0 )  THEN
              to_be_resorted => div_old
           ELSE
              IF ( .NOT. ALLOCATED( div_old_av ) )  THEN
                 ALLOCATE( div_old_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                 div_old_av = 0.0_wp
              ENDIF
              to_be_resorted => div_old_av
           ENDIF
           flag_nr = 0

           IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'kfd*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = kfd(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( kfd_av ) )  THEN
                ALLOCATE( kfd_av(nysg:nyng,nxlg:nxrg) )
                kfd_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = kfd_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'hr_xy', 'hr_xz', 'hr_yz' )
           IF ( av == 0 )  THEN
              to_be_resorted => hr
           ELSE
              IF ( .NOT. ALLOCATED( hr_av ) )  THEN
                 ALLOCATE( hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                 hr_av = 0.0_wp
              ENDIF
              to_be_resorted => hr_av
           ENDIF
           flag_nr = 0

           IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'vf25m*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_25m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vf_25m_av ) )  THEN
                ALLOCATE( vf_25m_av(nysg:nyng,nxlg:nxrg) )
                vf_25m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_25m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'vf50m*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_50m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vf_50m_av ) )  THEN
                ALLOCATE( vf_50m_av(nysg:nyng,nxlg:nxrg) )
                vf_50m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_50m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'vf75m*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_75m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vf_75m_av ) )  THEN
                ALLOCATE( vf_75m_av(nysg:nyng,nxlg:nxrg) )
                vf_75m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_75m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'vf100m*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_100m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vf_100m_av ) )  THEN
                ALLOCATE( vf_100m_av(nysg:nyng,nxlg:nxrg) )
                vf_100m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_100m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'vfxxm*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_xxm(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vf_xxm_av ) )  THEN
                ALLOCATE( vf_xxm_av(nysg:nyng,nxlg:nxrg) )
                vf_xxm_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vf_xxm_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'vfd25m*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vfd_25m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vfd_25m_av ) )  THEN
                ALLOCATE( vfd_25m_av(nysg:nyng,nxlg:nxrg) )
                vfd_25m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vfd_25m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'vfd50m*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vfd_50m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vfd_50m_av ) )  THEN
                ALLOCATE( vfd_50m_av(nysg:nyng,nxlg:nxrg) )
                vfd_50m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vfd_50m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'vfd75m*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                    &
                      local_pf(i,j,nzb+1) = vfd_75m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vfd_75m_av ) )  THEN
                ALLOCATE( vfd_75m_av(nysg:nyng,nxlg:nxrg) )
                vfd_75m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vfd_75m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'vfd100m*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vfd_100m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vfd_100m_av ) )  THEN
                ALLOCATE( vfd_100m_av(nysg:nyng,nxlg:nxrg) )
                vfd_100m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vfd_100m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'vfdxxm*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vfd_xxm(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( vfd_xxm_av ) )  THEN
                ALLOCATE( vfd_xxm_av(nysg:nyng,nxlg:nxrg) )
                vfd_xxm_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( surf_def%start_index(j,i) <= surf_def%end_index(j,i)  .OR.                 &
                        surf_lsm%start_index(j,i) <= surf_lsm%end_index(j,i) )                     &
                      local_pf(i,j,nzb+1) = vfd_xxm_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'rh_xy', 'rh_xz', 'rh_yz' )
           IF ( av == 0 )  THEN
              to_be_resorted => rh
           ELSE
              IF ( .NOT. ALLOCATED( rh_av ) )  THEN
                 ALLOCATE( rh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                 rh_av = 0.0_wp
              ENDIF
              to_be_resorted => rh_av
           ENDIF
           flag_nr = 0

           IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'ta_xy', 'ta_xz', 'ta_yz' )
           IF ( av == 0 )  THEN
              to_be_resorted => ta
           ELSE
              IF ( .NOT. ALLOCATED( ta_av ) )  THEN
                 ALLOCATE( ta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                 ta_av = 0.0_wp
              ENDIF
              to_be_resorted => ta_av
           ENDIF
           flag_nr = 0

           IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'ti_xy', 'ti_xz', 'ti_yz' )
           IF ( av == 0 )  THEN
              to_be_resorted => ti
           ELSE
              IF ( .NOT. ALLOCATED( ti_av ) )  THEN
                 ALLOCATE( ti_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                 ti_av = 0.0_wp
              ENDIF
              to_be_resorted => ti_av
           ENDIF
           flag_nr = 0

           IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'uu_xy', 'uu_xz', 'uu_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => uu
          ELSE
             IF ( .NOT. ALLOCATED( uu_av ) )  THEN
                ALLOCATE( uu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uu_av = 0.0_wp
             ENDIF
             to_be_resorted => uu_av
          ENDIF
          flag_nr = 1

          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'uv_xy', 'uv_xz', 'uv_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => uv
          ELSE
             IF ( .NOT. ALLOCATED( uv_av ) )  THEN
                ALLOCATE( uv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uv_av = 0.0_wp
             ENDIF
             to_be_resorted => uv_av
          ENDIF
          flag_nr = 1

          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'uw_xy', 'uw_xz', 'uw_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => uw
          ELSE
             IF ( .NOT. ALLOCATED( uw_av ) )  THEN
                ALLOCATE( uw_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uw_av = 0.0_wp
             ENDIF
             to_be_resorted => uw_av
          ENDIF
          flag_nr = 1

          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'vu_xy', 'vu_xz', 'vu_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => vu
          ELSE
             IF ( .NOT. ALLOCATED( vu_av ) )  THEN
                ALLOCATE( vu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vu_av = 0.0_wp
             ENDIF
             to_be_resorted => vu_av
          ENDIF
          flag_nr = 2

          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'vv_xy', 'vv_xz', 'vv_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => vv
          ELSE
             IF ( .NOT. ALLOCATED( vv_av ) )  THEN
                ALLOCATE( vv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vv_av = 0.0_wp
             ENDIF
             to_be_resorted => vv_av
          ENDIF
          flag_nr = 2

          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'vw_xy', 'vw_xz', 'vw_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => vw
          ELSE
             IF ( .NOT. ALLOCATED( vw_av ) )  THEN
                ALLOCATE( vw_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vw_av = 0.0_wp
             ENDIF
             to_be_resorted => vw_av
          ENDIF
          flag_nr = 2

          IF ( mode == 'xy' )  grid = 'zu'


       CASE ( 'wu_xy', 'wu_xz', 'wu_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wu
          ELSE
             IF ( .NOT. ALLOCATED( wu_av ) )  THEN
                ALLOCATE( wu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wu_av = 0.0_wp
             ENDIF
             to_be_resorted => wu_av
          ENDIF
          flag_nr = 3

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'wv_xy', 'wv_xz', 'wv_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wv
          ELSE
             IF ( .NOT. ALLOCATED( wv_av ) )  THEN
                ALLOCATE( wv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wv_av = 0.0_wp
             ENDIF
             to_be_resorted => wv_av
          ENDIF
          flag_nr = 3

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'ww_xy', 'ww_xz', 'ww_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => ww
          ELSE
             IF ( .NOT. ALLOCATED( ww_av ) )  THEN
                ALLOCATE( ww_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ww_av = 0.0_wp
             ENDIF
             to_be_resorted => ww_av
          ENDIF
          flag_nr = 3

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'wtheta_xy', 'wtheta_xz', 'wtheta_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wtheta
          ELSE
             IF ( .NOT. ALLOCATED( wtheta_av ) )  THEN
                ALLOCATE( wtheta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wtheta_av = 0.0_wp
             ENDIF
             to_be_resorted => wtheta_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'wq_xy', 'wq_xz', 'wq_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wq
          ELSE
             IF ( .NOT. ALLOCATED( wq_av ) )  THEN
                ALLOCATE( wq_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wq_av = 0.0_wp
             ENDIF
             to_be_resorted => wq_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'ws_xy', 'ws_xz', 'ws_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => ws
          ELSE
             IF ( .NOT. ALLOCATED( ws_av ) )  THEN
                ALLOCATE( ws_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ws_av = 0.0_wp
             ENDIF
             to_be_resorted => ws_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'theta_2m*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = pt_2m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( pt_2m_av ) )  THEN
                ALLOCATE( pt_2m_av(nysg:nyng,nxlg:nxrg) )
                pt_2m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = pt_2m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'qv_2m*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = qv_2m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( qv_2m_av ) )  THEN
                ALLOCATE( qv_2m_av(nysg:nyng,nxlg:nxrg) )
                qv_2m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = qv_2m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'ta_2m*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = ta_2m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( ta_2m_av ) )  THEN
                ALLOCATE( ta_2m_av(nysg:nyng,nxlg:nxrg) )
                ta_2m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = ta_2m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'wspeed_10m*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = uv_10m(j,i)
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( uv_10m_av ) )  THEN
                ALLOCATE( uv_10m_av(nysg:nyng,nxlg:nxrg) )
                uv_10m_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = uv_10m_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          two_d    = .TRUE.
          grid     = 'zu1'

       CASE ( 'wspeed_xy', 'wspeed_xz', 'wspeed_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wspeed
          ELSE
             IF ( .NOT. ALLOCATED( wspeed_av ) )  THEN
                ALLOCATE( wspeed_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wspeed_av = 0.0_wp
             ENDIF
             to_be_resorted => wspeed_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'wdir_xy', 'wdir_xz', 'wdir_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => wdir
          ELSE
             IF ( .NOT. ALLOCATED( wdir_av ) )  THEN
                ALLOCATE( wdir_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wdir_av = 0.0_wp
             ENDIF
             to_be_resorted => wdir_av
          ENDIF
          flag_nr = 0

          IF ( mode == 'xy' )  grid = 'zu'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT

    IF ( found  .AND.  .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                IF ( BTEST( topo_flags(k,j,i), flag_nr ) )  local_pf(i,j,k) = to_be_resorted(k,j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE doq_output_2d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av       !< index indicating averaged or instantaneous output
    INTEGER(iwp) ::  flag_nr  !< number of the topography flag (0: scalar, 1: u, 2: v, 3: w)
    INTEGER(iwp) ::  i        !< index variable along x-direction
    INTEGER(iwp) ::  j        !< index variable along y-direction
    INTEGER(iwp) ::  k        !< index variable along z-direction
    INTEGER(iwp) ::  nzb_do   !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do   !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL ::  found             !< true if variable is in list
    LOGICAL ::  resorted          !< true if array is resorted

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf        !<
    REAL(wp), DIMENSION(:,:,:), POINTER ::                 to_be_resorted  !< points to array which needs to be resorted for output


    flag_nr  = 0
    found    = .TRUE.
    resorted = .FALSE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'div_new' )
          IF ( av == 0 )  THEN
             to_be_resorted => div_new
          ELSE
             IF ( .NOT. ALLOCATED( div_new_av ) )  THEN
                ALLOCATE( div_new_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                div_new_av = 0.0_wp
             ENDIF
             to_be_resorted => div_new_av
          ENDIF
          flag_nr = 0

       CASE ( 'div_old' )
          IF ( av == 0 )  THEN
             to_be_resorted => div_old
          ELSE
             IF ( .NOT. ALLOCATED( div_old_av ) )  THEN
                ALLOCATE( div_old_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                div_old_av = 0.0_wp
             ENDIF
             to_be_resorted => div_old_av
          ENDIF
          flag_nr = 0

       CASE ( 'hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => hr
          ELSE
             IF ( .NOT. ALLOCATED( hr_av ) )  THEN
                ALLOCATE( hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                hr_av = 0.0_wp
             ENDIF
             to_be_resorted => hr_av
          ENDIF
          flag_nr = 0

       CASE ( 'rh' )
          IF ( av == 0 )  THEN
             to_be_resorted => rh
          ELSE
             IF ( .NOT. ALLOCATED( rh_av ) )  THEN
                ALLOCATE( rh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                rh_av = 0.0_wp
             ENDIF
             to_be_resorted => rh_av
          ENDIF
          flag_nr = 0

       CASE ( 'ta' )
          IF ( av == 0 )  THEN
             to_be_resorted => ta
          ELSE
             IF ( .NOT. ALLOCATED( ta_av ) )  THEN
                ALLOCATE( ta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ta_av = 0.0_wp
             ENDIF
             to_be_resorted => ta_av
          ENDIF
          flag_nr = 0

       CASE ( 'ti' )
          IF ( av == 0 )  THEN
             to_be_resorted => ti
          ELSE
             IF ( .NOT. ALLOCATED( ti_av ) )  THEN
                ALLOCATE( ti_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ti_av = 0.0_wp
             ENDIF
             to_be_resorted => ti_av
          ENDIF
          flag_nr = 0

       CASE ( 'uu_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => uu
          ELSE
             IF ( .NOT. ALLOCATED( uu_av ) )  THEN
                ALLOCATE( uu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uu_av = 0.0_wp
             ENDIF
             to_be_resorted => uu_av
          ENDIF
          flag_nr = 1

       CASE ( 'uv_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => uv
          ELSE
             IF ( .NOT. ALLOCATED( uv_av ) )  THEN
                ALLOCATE( uv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uv_av = 0.0_wp
             ENDIF
             to_be_resorted => uv_av
          ENDIF
          flag_nr = 1

       CASE ( 'uw_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => uw
          ELSE
             IF ( .NOT. ALLOCATED( uw_av ) )  THEN
                ALLOCATE( uw_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uw_av = 0.0_wp
             ENDIF
             to_be_resorted => uw_av
          ENDIF
          flag_nr = 1

       CASE ( 'vu_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => vu
          ELSE
             IF ( .NOT. ALLOCATED( vu_av ) )  THEN
                ALLOCATE( vu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vu_av = 0.0_wp
             ENDIF
             to_be_resorted => vu_av
          ENDIF
          flag_nr = 2

       CASE ( 'vv_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => vv
          ELSE
             IF ( .NOT. ALLOCATED( vv_av ) )  THEN
                ALLOCATE( vv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vv_av = 0.0_wp
             ENDIF
             to_be_resorted => vv_av
          ENDIF
          flag_nr = 2

       CASE ( 'vw_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => vw
          ELSE
             IF ( .NOT. ALLOCATED( vw_av ) )  THEN
                ALLOCATE( vw_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vw_av = 0.0_wp
             ENDIF
             to_be_resorted => vw_av
          ENDIF
          flag_nr = 2

       CASE ( 'wu_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => wu
          ELSE
             IF ( .NOT. ALLOCATED( wu_av ) )  THEN
                ALLOCATE( wu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wu_av = 0.0_wp
             ENDIF
             to_be_resorted => wu_av
          ENDIF
          flag_nr = 3

       CASE ( 'wv_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => wv
          ELSE
             IF ( .NOT. ALLOCATED( wv_av ) )  THEN
                ALLOCATE( wv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wv_av = 0.0_wp
             ENDIF
             to_be_resorted => wv_av
          ENDIF
          flag_nr = 3

       CASE ( 'ww_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => ww
          ELSE
             IF ( .NOT. ALLOCATED( ww_av ) )  THEN
                ALLOCATE( ww_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ww_av = 0.0_wp
             ENDIF
             to_be_resorted => ww_av
          ENDIF
          flag_nr = 3

       CASE ( 'wtheta_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => wtheta
          ELSE
             IF ( .NOT. ALLOCATED( wtheta_av ) )  THEN
                ALLOCATE( wtheta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wtheta_av = 0.0_wp
             ENDIF
             to_be_resorted => wtheta_av
          ENDIF
          flag_nr = 0

       CASE ( 'wq_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => wq
          ELSE
             IF ( .NOT. ALLOCATED( wq_av ) )  THEN
                ALLOCATE( wq_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wq_av = 0.0_wp
             ENDIF
             to_be_resorted => wq_av
          ENDIF
          flag_nr = 0

       CASE ( 'ws_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => ws
          ELSE
             IF ( .NOT. ALLOCATED( ws_av ) )  THEN
                ALLOCATE( ws_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ws_av = 0.0_wp
             ENDIF
             to_be_resorted => ws_av
          ENDIF
          flag_nr = 0

       CASE ( 'wspeed' )
          IF ( av == 0 )  THEN
             to_be_resorted => wspeed
          ELSE
             IF ( .NOT. ALLOCATED( wspeed_av ) )  THEN
                ALLOCATE( wspeed_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wspeed_av = 0.0_wp
             ENDIF
             to_be_resorted => wspeed_av
          ENDIF
          flag_nr = 0

       CASE ( 'wdir' )
          IF ( av == 0 )  THEN
             to_be_resorted => wdir
          ELSE
             IF ( .NOT. ALLOCATED( wdir_av ) )  THEN
                ALLOCATE( wdir_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wdir_av = 0.0_wp
             ENDIF
             to_be_resorted => wdir_av
          ENDIF
          flag_nr = 0

       CASE DEFAULT
          found = .FALSE.

    END SELECT

    IF ( found  .AND.  .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                IF ( BTEST( topo_flags(k,j,i), flag_nr ) )  local_pf(i,j,k) = to_be_resorted(k,j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE doq_output_3d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k) for masked data output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_output_mask( av, variable, found, local_pf, mid )

    CHARACTER (LEN=*) ::  variable   !<
    CHARACTER (LEN=5) ::  grid       !< flag to distinquish between staggered grids

    INTEGER(iwp) ::  av              !< index indicating averaged or instantaneous output
    INTEGER(iwp) ::  flag_nr         !< number of the topography flag (0: scalar, 1: u, 2: v, 3: w)
    INTEGER(iwp) ::  i               !< index variable along x-direction
    INTEGER(iwp) ::  j               !< index variable along y-direction
    INTEGER(iwp) ::  k               !< index variable along z-direction
    INTEGER(iwp) ::  im              !< loop index for masked variables
    INTEGER(iwp) ::  jm              !< loop index for masked variables
    INTEGER(iwp) ::  kk              !< masked output index variable along z-direction
    INTEGER(iwp) ::  mid             !< masked output running index
    INTEGER(iwp) ::  ktt             !< k index of highest horizontal surface

    LOGICAL      ::  found           !< true if variable is in list
    LOGICAL      ::  resorted        !< true if array is resorted

    REAL(wp), DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  local_pf   !<

    REAL(wp), DIMENSION(:,:,:), POINTER  ::  to_be_resorted  !< points to array which needs to be resorted for output


    flag_nr  = 0
    found    = .TRUE.
    resorted = .FALSE.
    grid     = 's'

    SELECT CASE ( TRIM( variable ) )

      CASE ( 'div_new' )
          IF ( av == 0 )  THEN
             to_be_resorted => div_new
          ELSE
             to_be_resorted => div_new_av
          ENDIF
          grid = 's'
          flag_nr = 0

      CASE ( 'div_old' )
          IF ( av == 0 )  THEN
             to_be_resorted => div_old
          ELSE
             to_be_resorted => div_old_av
          ENDIF
          grid = 's'
          flag_nr = 0

      CASE ( 'hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => hr
          ELSE
             to_be_resorted => hr_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'rh' )
          IF ( av == 0 )  THEN
             to_be_resorted => rh
          ELSE
             to_be_resorted => rh_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'ta' )
          IF ( av == 0 )  THEN
             to_be_resorted => ta
          ELSE
             to_be_resorted => ta_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'ti' )
          IF ( av == 0 )  THEN
             to_be_resorted => ti
          ELSE
             to_be_resorted => ti_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'uu_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => uu
          ELSE
             to_be_resorted => uu_av
          ENDIF
          grid = 'u'
          flag_nr = 1

       CASE ( 'uv_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => uv
          ELSE
             to_be_resorted => uv_av
          ENDIF
          grid = 'u'
          flag_nr = 1

       CASE ( 'uw_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => uw
          ELSE
             to_be_resorted => uw_av
          ENDIF
          grid = 'u'
          flag_nr = 1

       CASE ( 'vu_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => vu
          ELSE
             to_be_resorted => vu_av
          ENDIF
          grid = 'v'
          flag_nr = 2

       CASE ( 'vv_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => vv
          ELSE
             to_be_resorted => vv_av
          ENDIF
          grid = 'v'
          flag_nr = 2

       CASE ( 'vw_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => vw
          ELSE
             to_be_resorted => vw_av
          ENDIF
          grid = 'v'
          flag_nr = 2

       CASE ( 'wu_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => wu
          ELSE
             to_be_resorted => wu_av
          ENDIF
          grid = 'w'
          flag_nr = 3

       CASE ( 'wv_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => wv
          ELSE
             to_be_resorted => wv_av
          ENDIF
          grid = 'w'
          flag_nr = 3

       CASE ( 'ww_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => ww
          ELSE
             to_be_resorted => ww_av
          ENDIF
          grid = 'w'
          flag_nr = 3

       CASE ( 'wtheta_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => wtheta
          ELSE
             to_be_resorted => wtheta_av
          ENDIF
          grid = 'w'
          flag_nr = 0

       CASE ( 'wq_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => wq
          ELSE
             to_be_resorted => wq_av
          ENDIF
          grid = 'w'
          flag_nr = 0

       CASE ( 'ws_product' )
          IF ( av == 0 )  THEN
             to_be_resorted => ws
          ELSE
             to_be_resorted => ws_av
          ENDIF
          grid = 'w'
          flag_nr = 0

       CASE ( 'wspeed' )
          IF ( av == 0 )  THEN
             to_be_resorted => wspeed
          ELSE
             to_be_resorted => wspeed_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE ( 'wdir' )
          IF ( av == 0 )  THEN
             to_be_resorted => wdir
          ELSE
             to_be_resorted => wdir_av
          ENDIF
          grid = 's'
          flag_nr = 0

       CASE DEFAULT
          found = .FALSE.

    END SELECT

    IF ( found  .AND.  .NOT. resorted )  THEN
       IF ( .NOT. mask_surface(mid) )  THEN
!
!--       Default masked output.
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
                DO  k = 1, mask_size_l(mid,3)
                   IF ( BTEST( topo_flags(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i)), flag_nr ) )  &
                   THEN
                      local_pf(i,j,k) = to_be_resorted(mask_k(mid,k),mask_j(mid,j),mask_i(mid,i))
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ELSE
!
!--       Terrain-following masked output.
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
!
!--             Get k index of the lowest non-terrain grid point, else set fill value.
                im = mask_i(mid,i)
                jm = mask_j(mid,j)
                ktt = MINLOC( MERGE( 1, 0, BTEST( topo_flags(:,jm,im), 5 ) ), DIM=1 ) - 1
                DO  k = 1, mask_size_l(mid,3)
                   kk = MIN( ktt + mask_k(mid,k) - 1, nzt+1 )
!
!--                Set value if not in building.
                   IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
                      local_pf(i,j,k) = to_be_resorted(kk,jm,im)
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ENDIF
    ENDIF

 END SUBROUTINE doq_output_mask


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate required arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_init

    IMPLICIT NONE

    INTEGER(iwp) ::  ivar   !< loop index over all 2d/3d/mask output quantities


!
!-- Next line is to avoid compiler warnings about unused variables
    IF ( timestep_number_at_prev_calc == 0 )  CONTINUE
!
!-- Preparatory steps and initialization of output arrays
    IF ( .NOT.  prepared_diagnostic_output_quantities )  CALL doq_prepare

    initialized_diagnostic_output_quantities = .FALSE.

    ivar = 1

    DO  WHILE ( ivar <= SIZE( do_all ) )

       SELECT CASE ( TRIM( do_all(ivar) ) )
!
!--       Allocate arrays for 'flow field divergence'.
          CASE ( 'div_new' )
             IF ( .NOT. ALLOCATED( div_new ) )  THEN
                ALLOCATE( div_new(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                div_new = 0.0_wp
             ENDIF

          CASE ( 'div_old' )
             IF ( .NOT. ALLOCATED( div_old ) )  THEN
                ALLOCATE( div_old(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                div_old = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic flow depth'.
          CASE ( 'kfd*' )
             IF ( .NOT. ALLOCATED( kfd ) )  THEN
                ALLOCATE( kfd(nysg:nyng,nxlg:nxrg) )
                kfd = 0.0_wp
             ENDIF
             IF ( .NOT. ALLOCATED( count_kfd ) )  THEN
!
!--             Allocate cound_kfd without ghost layers since this has not been implemented
!--             restart IO yet.
                ALLOCATE( count_kfd(nys:nyn,nxl:nxr) )
                count_kfd = 0
             ENDIF
!
!--       Allocate array for 'heating rate'.
          CASE ( 'hr' )
             IF ( .NOT. ALLOCATED( hr ) )  THEN
                ALLOCATE( hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                hr = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux' (integrated up to 25m).
          CASE ( 'vf25m*' )
             IF ( .NOT. ALLOCATED( vf_25m ) )  THEN
                ALLOCATE( vf_25m(nysg:nyng,nxlg:nxrg) )
                vf_25m = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux' (integrated up to 50m).
          CASE ( 'vf50m*' )
             IF ( .NOT. ALLOCATED( vf_50m ) )  THEN
                ALLOCATE( vf_50m(nysg:nyng,nxlg:nxrg) )
                vf_50m = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux' (integrated up to 75m).
          CASE ( 'vf75m*' )
             IF ( .NOT. ALLOCATED( vf_75m ) )  THEN
                ALLOCATE( vf_75m(nysg:nyng,nxlg:nxrg) )
                vf_75m = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux' (integrated up to 100m).
          CASE ( 'vf100m*' )
             IF ( .NOT. ALLOCATED( vf_100m ) )  THEN
                ALLOCATE( vf_100m(nysg:nyng,nxlg:nxrg) )
                vf_100m = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux' (integrated up to variable height
!--       according to kfd*).
          CASE ( 'vfxxm*' )
             IF ( .NOT. ALLOCATED( vf_xxm ) )  THEN
                ALLOCATE( vf_xxm(nysg:nyng,nxlg:nxrg) )
                vf_xxm = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux density' (integrated up to 50m).
          CASE ( 'vfd25m*' )
             IF ( .NOT. ALLOCATED( vfd_25m ) )  THEN
                ALLOCATE( vfd_25m(nysg:nyng,nxlg:nxrg) )
                vfd_25m = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux density' (integrated up to 50m).
          CASE ( 'vfd50m*' )
             IF ( .NOT. ALLOCATED( vfd_50m ) )  THEN
                ALLOCATE( vfd_50m(nysg:nyng,nxlg:nxrg) )
                vfd_50m = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux density' (integrated up to 75m).
          CASE ( 'vfd75m*' )
             IF ( .NOT. ALLOCATED( vfd_75m ) )  THEN
                ALLOCATE( vfd_75m(nysg:nyng,nxlg:nxrg) )
                vfd_75m = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux density' (integrated up to 100m).
          CASE ( 'vfd100m*' )
             IF ( .NOT. ALLOCATED( vfd_100m ) )  THEN
                ALLOCATE( vfd_100m(nysg:nyng,nxlg:nxrg) )
                vfd_100m = 0.0_wp
             ENDIF
!
!--       Allocate array for 'katabatic volume flux density' (integrated up to variable height
!--       according to kfd*).
          CASE ( 'vfdxxm*' )
             IF ( .NOT. ALLOCATED( vfd_xxm ) )  THEN
                ALLOCATE( vfd_xxm(nysg:nyng,nxlg:nxrg) )
                vfd_xxm = 0.0_wp
             ENDIF
!
!--       Allocate array for 'relative humidity'
          CASE ( 'rh' )
             IF ( .NOT. ALLOCATED( rh ) )  THEN
                ALLOCATE( rh(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                rh = 0.0_wp
             ENDIF
!
!--       Allocate array for 'air temperature'
          CASE ( 'ta' )
             IF ( .NOT. ALLOCATED( ta ) )  THEN
                ALLOCATE( ta(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ta = 0.0_wp
             ENDIF
!
!--       Allocate array for 'turbulence intensity'
          CASE ( 'ti' )
             IF ( .NOT. ALLOCATED( ti ) )  THEN
                ALLOCATE( ti(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ti = 0.0_wp
             ENDIF
!
!--       Allocate array for uu
          CASE ( 'uu_product' )
             IF ( .NOT. ALLOCATED( uu ) )  THEN
                ALLOCATE( uu(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uu = 0.0_wp
             ENDIF
!
!--       Allocate array for uv
          CASE ( 'uv_product' )
             IF ( .NOT. ALLOCATED( uv ) )  THEN
                ALLOCATE( uv(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uv = 0.0_wp
             ENDIF
!
!--       Allocate array for uw
          CASE ( 'uw_product' )
             IF ( .NOT. ALLOCATED( uw ) )  THEN
                ALLOCATE( uw(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                uw = 0.0_wp
             ENDIF
!
!--       Allocate array for vu
          CASE ( 'vu_product' )
             IF ( .NOT. ALLOCATED( vu ) )  THEN
                ALLOCATE( vu(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vu = 0.0_wp
             ENDIF
!
!--       Allocate array for vv
          CASE ( 'vv_product' )
             IF ( .NOT. ALLOCATED( vv ) )  THEN
                ALLOCATE( vv(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vv = 0.0_wp
             ENDIF
!
!--       Allocate array for vw
          CASE ( 'vw_product' )
             IF ( .NOT. ALLOCATED( vw ) )  THEN
                ALLOCATE( vw(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                vw = 0.0_wp
             ENDIF
!
!--       Allocate array for wu
          CASE ( 'wu_product' )
             IF ( .NOT. ALLOCATED( wu ) )  THEN
                ALLOCATE( wu(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wu = 0.0_wp
             ENDIF
!
!--       Allocate array for wv
          CASE ( 'wv_product' )
             IF ( .NOT. ALLOCATED( wv ) )  THEN
                ALLOCATE( wv(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wv = 0.0_wp
             ENDIF
!
!--       Allocate array for ww
          CASE ( 'ww_product' )
             IF ( .NOT. ALLOCATED( ww ) )  THEN
                ALLOCATE( ww(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ww = 0.0_wp
             ENDIF
!
!--       Allocate array for wtheta
          CASE ( 'wtheta_product' )
             IF ( .NOT. ALLOCATED( wtheta ) )  THEN
                ALLOCATE( wtheta(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wtheta = 0.0_wp
             ENDIF
!
!--       Allocate array for wq
          CASE ( 'wq_product' )
             IF ( .NOT. ALLOCATED( wq ) )  THEN
                ALLOCATE( wq(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wq = 0.0_wp
             ENDIF
!
!--       Allocate array for ws
          CASE ( 'ws_product' )
             IF ( .NOT. ALLOCATED( ws ) )  THEN
                ALLOCATE( ws(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ws = 0.0_wp
             ENDIF
!
!--       Allocate arrays for 2-m potential and non-potential
!--       temperature. Both, pt_2m and ta_2m are required by
!--       calc_2m_temperature even if only one of the variables is
!--       requested.
          CASE ( 'theta_2m*', 'ta_2m*' )
             IF ( .NOT. ALLOCATED( pt_2m ) )  THEN
                ALLOCATE( pt_2m(nysg:nyng,nxlg:nxrg) )
                pt_2m = 0.0_wp
             ENDIF
             IF ( .NOT. ALLOCATED( ta_2m ) )  THEN
                ALLOCATE( ta_2m(nysg:nyng,nxlg:nxrg) )
                ta_2m = 0.0_wp
             ENDIF
!
!--       Allocate arrays for 2-m water vapor mixing ratio.
          CASE ( 'qv_2m*'  )
             IF ( .NOT. ALLOCATED( qv_2m ) )  THEN
                ALLOCATE( qv_2m(nysg:nyng,nxlg:nxrg) )
                qv_2m = 0.0_wp
             ENDIF
!
!--       Allocate array for 10-m wind speed
          CASE ( 'wspeed_10m*' )
             IF ( .NOT. ALLOCATED( uv_10m ) )  THEN
                ALLOCATE( uv_10m(nysg:nyng,nxlg:nxrg) )
                uv_10m = 0.0_wp
             ENDIF
!
!--       Allocate array for wspeed
          CASE ( 'wspeed' )
             IF ( .NOT. ALLOCATED( wspeed ) )  THEN
                ALLOCATE( wspeed(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wspeed = 0.0_wp
             ENDIF
!
!--       Allocate array for wdir
          CASE ( 'wdir' )
             IF ( .NOT. ALLOCATED( u_center ) )  THEN
                ALLOCATE( u_center(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                u_center = 0.0_wp
             ENDIF
             IF ( .NOT. ALLOCATED( v_center ) )  THEN
                ALLOCATE( v_center(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                v_center = 0.0_wp
             ENDIF
             IF ( .NOT. ALLOCATED( wdir ) )  THEN
                ALLOCATE( wdir(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wdir = 0.0_wp
             ENDIF

       END SELECT

       ivar = ivar + 1
    ENDDO

    ivar = 1
!
!-- Scan output variables again. This is required for special actions some output needs. At the
!-- moment this is needed to trigger terrain slope angle computation.
    DO  WHILE ( ivar <= SIZE( do_all ) )

       SELECT CASE ( TRIM( do_all(ivar) ) )
!
!--       In case of any of the following output variables using katabatic-flow depth
!--       detection,  the slope angle needs to be comuted.
          CASE ( 'kfd*', 'vfxxm*', 'vfdxxm*' )
             CALL doq_compute_slope_angle
       END SELECT
       ivar = ivar + 1
    ENDDO

    initialized_diagnostic_output_quantities = .TRUE.

 END SUBROUTINE doq_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of module-specific statistics, i.e. horizontally averaged profiles and time series.
!> This is called for every statistic region sr, but at least for the region "total domain" (sr=0).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_statistics( mode, sr, tn )


    CHARACTER (LEN=*) ::  mode   !<

!     INTEGER(iwp) ::  i    !<
!     INTEGER(iwp) ::  j    !<
!     INTEGER(iwp) ::  k    !<
    INTEGER(iwp) ::  sr   !<
    INTEGER(iwp) ::  tn   !<
!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( sr == 0  .OR.  tn == 0 )  CONTINUE

    IF ( mode == 'profiles' )  THEN

    ELSEIF ( mode == 'time_series' )  THEN

    ENDIF

 END SUBROUTINE doq_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate diagnostic quantities
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_calculate

    IMPLICIT NONE

    INTEGER(iwp) ::  i          !< running index x-dimension
    INTEGER(iwp) ::  ivar       !< loop index over all 2d/3d/mask output quantities
    INTEGER(iwp) ::  j          !< running index y-dimension
    INTEGER(iwp) ::  k          !< running index z-dimension
    INTEGER(iwp) ::  k_detect   !< running index z-dimension
    INTEGER(iwp) ::  k_max      !< grid index along z-dimension indicating upper integration bound or location of maximum
    INTEGER(iwp) ::  kk         !< derived grid index along z-dimension on terrain-following coordinates
    INTEGER(iwp) ::  m          !< running index surface elements

    REAL(wp) ::  ddxyz = 0.0_wp     !< dummy for ddx, ddy, or ddzw, depending on wall facing
    REAL(wp) ::  e_s                !< saturation vapor pressure
    REAL(wp) ::  dudx               !< gradient of u in x-direction
    REAL(wp) ::  dudy               !< gradient of u in y-direction
    REAL(wp) ::  dudz               !< gradient of u in z-direction
    REAL(wp) ::  dvdx               !< gradient of v in x-direction
    REAL(wp) ::  dvdy               !< gradient of v in y-direction
    REAL(wp) ::  dvdz               !< gradient of v in z-direction
    REAL(wp) ::  dwdx               !< gradient of w in x-direction
    REAL(wp) ::  dwdy               !< gradient of w in y-direction
    REAL(wp) ::  dwdz               !< gradient of w in z-direction
    REAL(wp) ::  integration_height !< level up to which volume fluxes are integrated
    REAL(wp) ::  q_s                !< saturation mixing ratio
    REAL(wp) ::  temp               !< temperature
    REAL(wp) ::  uuf                !< approximated resolved-scale flux u'u' at grid center
    REAL(wp) ::  uvf                !< approximated resolved-scale flux u'v' at grid center
    REAL(wp) ::  uwf                !< approximated resolved-scale flux u'w' at grid center
    REAL(wp) ::  vuf                !< approximated resolved-scale flux v'u' at grid center
    REAL(wp) ::  vvf                !< approximated resolved-scale flux v'v' at grid center
    REAL(wp) ::  vwf                !< approximated resolved-scale flux v'w' at grid center
    REAL(wp) ::  wuf                !< approximated resolved-scale flux w'u' at grid center
    REAL(wp) ::  wvf                !< approximated resolved-scale flux w'v' at grid center
    REAL(wp) ::  wwf                !< approximated resolved-scale flux w'w' at grid center

    REAL(wp), DIMENSION(nzb+1:nzt+1) ::  shear_production !< production of TKE by shear

    TYPE(surf_type), POINTER ::  surf     !< surf-type array, used to generalize subroutines


!     CALL cpu_log( log_point(41), 'calculate_quantities', 'start' )

!
!-- Save timestep number to check in time_integration if doq_calculate has been called already,
!-- since the CALL occurs at two locations, but the calculations need to be done only once per
!-- timestep.
    timestep_number_at_prev_calc = current_timestep_number

    ivar = 1

    DO  WHILE ( ivar <= SIZE( do_all ) )

       SELECT CASE ( TRIM( do_all(ivar) ) )
!
!--       Katabatic flow depth
          CASE ( 'kfd*' )
             kfd = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   CALL detect_katabatic_flow( k_detect )

                   IF ( k_detect /= -1 )  THEN
                      kfd(j,i) = zu(k_detect) - zw(topo_top_ind(j,i,0))
                      count_kfd(j,i) = count_kfd(j,i) + 1
                   ENDIF
                ENDDO
             ENDDO
!
!--       Heating rate caused by longwave emission of air-volume (only RRTMG), emission by
!--       plants, and by surface fluxes (in K/s).
          CASE ( 'hr' )
#if defined( __rrtmg ) || defined ( __tenstream )
!
!--          Heating of air-volume by longwave emission of air volume
             IF ( radiation_scheme == 'rrtmg' )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         hr(k,j,i) = ( rad_lw_hr(k,j,i) + rad_sw_hr(k,j,i) )                       &
                                    * d_exner(k) / seconds_per_hour                                &
                                    * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
#endif
!
!--          Heating of air-volume by emission of plant canopy
             IF ( plant_canopy )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
!
!--                   Determine topography-top index on scalar-grid
                      DO  k = topo_top_ind(j,i,0)+1, topo_top_ind(j,i,0) + pch_index_ji(j,i)
                         kk = k - topo_top_ind(j,i,0)                 ! lad arrays are defined flat
                         hr(k,j,i) = hr(k,j,i) + pcm_sensiblerate(kk,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--          Heating of air-volume by surface interaction
             IF ( .NOT. neutral )  THEN
                DO  m = 1, surf_def%ns
                   ddxyz = MERGE( ddzw(k), ddxyz, surf_def%upward(m)     .OR.  surf_def%downward(m)  )
                   ddxyz = MERGE( ddy,     ddxyz, surf_def%northward(m)  .OR.  surf_def%southward(m) )
                   ddxyz = MERGE( ddx,     ddxyz, surf_def%eastward(m)   .OR.  surf_def%westward(m)  )
                   i = surf_def%i(m)
                   j = surf_def%j(m)
                   k = surf_def%k(m)
                   hr(k,j,i) = hr(k,j,i) + surf_def%shf(m) * ddxyz
                ENDDO
                DO  m = 1, surf_lsm%ns
                   ddxyz = MERGE( ddzw(k), ddxyz, surf_lsm%upward(m)     .OR.  surf_lsm%downward(m)  )
                   ddxyz = MERGE( ddy,     ddxyz, surf_lsm%northward(m)  .OR.  surf_lsm%southward(m) )
                   ddxyz = MERGE( ddx,     ddxyz, surf_lsm%eastward(m)   .OR.  surf_lsm%westward(m)  )
                   i = surf_lsm%i(m)
                   j = surf_lsm%j(m)
                   k = surf_lsm%k(m)
                   hr(k,j,i) = hr(k,j,i) + surf_lsm%shf(m) * ddxyz
                ENDDO
                DO  m = 1, surf_usm%ns
                   ddxyz = MERGE( ddzw(k), ddxyz, surf_usm%upward(m)     .OR.  surf_usm%downward(m)  )
                   ddxyz = MERGE( ddy,     ddxyz, surf_usm%northward(m)  .OR.  surf_usm%southward(m) )
                   ddxyz = MERGE( ddx,     ddxyz, surf_usm%eastward(m)   .OR.  surf_usm%westward(m)  )
                   i = surf_usm%i(m)
                   j = surf_usm%j(m)
                   k = surf_usm%k(m)
                   hr(k,j,i) = hr(k,j,i) + surf_usm%shf(m) * ddxyz
                ENDDO
             ENDIF
!
!--       Katabatic volume flux (integrated up to 25m). Note, volume flux is not evaluated
!--       above buildings.
          CASE ( 'vf25m*' )
!
!--          Set the integration height (starting from the surface) and determine the index for
!--          the upper integration height.
             integration_height = 25.0_wp
             CALL find_k_index( integration_height )
             vf_25m = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   surf => surf_def
                   CALL integrate_volume_flux( vf_25m(j,i), integration_height )
                   surf => surf_lsm
                   CALL integrate_volume_flux( vf_25m(j,i), integration_height )
                ENDDO
             ENDDO

!
!--       Katabatic volume flux (integrated up to 50m). Note, volume flux is not evaluated
!--       above buildings.
          CASE ( 'vf50m*' )
!
!--          Set the integration height (starting from the surface) and determine the index for
!--          the upper integration height.
             integration_height = 50.0_wp
             CALL find_k_index( integration_height )
             vf_50m = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   surf => surf_def
                   CALL integrate_volume_flux( vf_50m(j,i), integration_height )
                   surf => surf_lsm
                   CALL integrate_volume_flux( vf_50m(j,i), integration_height )
                ENDDO
             ENDDO
!
!--       Katabatic volume flux (integrated up to 75m). Note, volume flux is not evaluated
!--       above buildings.
          CASE ( 'vf75m*' )
!
!--          Set the integration height (starting from the surface) and determine the index for
!--          the upper integration height.
             integration_height = 75.0_wp
             CALL find_k_index( integration_height )
             vf_75m = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   surf => surf_def
                   CALL integrate_volume_flux( vf_75m(j,i), integration_height )
                   surf => surf_lsm
                   CALL integrate_volume_flux( vf_75m(j,i), integration_height )
                ENDDO
             ENDDO
!
!--       Katabatic volume flux (integrated up to 100m). Note, volume flux is not evaluated
!--       above buildings.
          CASE ( 'vf100m*' )
!
!--          Set the integration height (starting from the surface) and determine the index for
!--          the upper integration height.
             integration_height = 100.0_wp
             CALL find_k_index( integration_height )
             vf_100m = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   surf => surf_def
                   CALL integrate_volume_flux( vf_100m(j,i), integration_height )
                   surf => surf_lsm
                   CALL integrate_volume_flux( vf_100m(j,i), integration_height )
                ENDDO
             ENDDO
!
!--       Katabatic volume flux (integrated up to variable depth). Note, volume flux
!--       is not evaluated above buildings.
          CASE ( 'vfxxm*' )
             vf_xxm = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Detect depth of the katabatic flow (if present).
                   CALL detect_katabatic_flow( k_detect )

                   IF ( k_detect /= -1 )  THEN
!
!--                   Set the integration height (starting from the surface) and determine the
!--                   index for the upper integration height.
                      integration_height = zu(k_detect) - zw(topo_top_ind(j,i,0))

                      CALL find_k_index( integration_height )
                      surf => surf_def
                      CALL integrate_volume_flux( vf_xxm(j,i), integration_height )
                      surf => surf_lsm
                      CALL integrate_volume_flux( vf_xxm(j,i), integration_height )
                   ENDIF

                ENDDO
             ENDDO
!
!--       Katabatic volume flux density (integrated up to 25m). Note, volume flux is not evaluated
!--       above buildings.
          CASE ( 'vfd25m*' )
!
!--          Set the integration height (starting from the surface) and determine the index for
!--          the upper integration height.
             integration_height = 25.0_wp
             CALL find_k_index( integration_height )
             vfd_25m = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   surf => surf_def
                   CALL integrate_volume_flux( vfd_25m(j,i), integration_height,                   &
                                               compute_volume_flux_density = .TRUE. )
                   surf => surf_lsm
                   CALL integrate_volume_flux( vfd_25m(j,i), integration_height,                   &
                                               compute_volume_flux_density = .TRUE. )
                ENDDO
             ENDDO

!
!--       Katabatic volume flux density (integrated up to 50m). Note, volume flux is not evaluated
!--       above buildings.
          CASE ( 'vfd50m*' )
!
!--          Set the integration height (starting from the surface) and determine the index for
!--          the upper integration height.
             integration_height = 50.0_wp
             CALL find_k_index( integration_height )
             vfd_50m = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   surf => surf_def
                   CALL integrate_volume_flux( vfd_50m(j,i), integration_height,                   &
                                               compute_volume_flux_density = .TRUE. )
                   surf => surf_lsm
                   CALL integrate_volume_flux( vfd_50m(j,i), integration_height,                   &
                                               compute_volume_flux_density = .TRUE. )
                ENDDO
             ENDDO
!
!--       Katabatic volume flux density (integrated up to 75m). Note, volume flux is not evaluated
!--       above buildings.
          CASE ( 'vfd75m*' )
!
!--          Set the integration height (starting from the surface) and determine the index for
!--          the upper integration height.
             integration_height = 75.0_wp
             CALL find_k_index( integration_height )
             vfd_75m = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   surf => surf_def
                   CALL integrate_volume_flux( vfd_75m(j,i), integration_height,                   &
                                               compute_volume_flux_density = .TRUE. )
                   surf => surf_lsm
                   CALL integrate_volume_flux( vfd_75m(j,i), integration_height,                   &
                                               compute_volume_flux_density = .TRUE. )
                ENDDO
             ENDDO
!
!--       Katabatic volume flux density (integrated up to 100m). Note, volume flux is not evaluated
!--       above buildings.
          CASE ( 'vfd100m*' )
!
!--          Set the integration height (starting from the surface) and determine the index for
!--          the upper integration height.
             integration_height = 100.0_wp
             CALL find_k_index( integration_height )
             vfd_100m = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   surf => surf_def
                   CALL integrate_volume_flux( vfd_100m(j,i), integration_height,                   &
                                               compute_volume_flux_density = .TRUE. )
                   surf => surf_lsm
                   CALL integrate_volume_flux( vfd_100m(j,i), integration_height,                   &
                                               compute_volume_flux_density = .TRUE. )
                ENDDO
             ENDDO
!
!--       Katabatic volume flux density (integrated up to variable depth). Note, volume flux
!--       is not evaluated above buildings.
          CASE ( 'vfdxxm*' )
             vfd_xxm = 0.0_wp
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Detect depth of the katabatic flow (if present).
                   CALL detect_katabatic_flow( k_detect )

                   IF ( k_detect /= -1 )  THEN
!
!--                   Set the integration height (starting from the surface) and determine the
!--                   index for the upper integration height.
                      integration_height = zu(k_detect) - zw(topo_top_ind(j,i,0))

                      CALL find_k_index( integration_height )
                      surf => surf_def
                      CALL integrate_volume_flux( vfd_xxm(j,i), integration_height,                &
                                                  compute_volume_flux_density = .TRUE. )
                      surf => surf_lsm
                      CALL integrate_volume_flux( vfd_xxm(j,i), integration_height,                &
                                                  compute_volume_flux_density = .TRUE. )
                   ENDIF
                ENDDO
             ENDDO
!
!--       rh (relative humidity)
          CASE ( 'rh' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt+1
                      IF ( humidity )  THEN
                         temp = exner(k) * pt(k,j,i)
                         e_s = magnus( temp )
                         q_s = rd_d_rv * e_s / ( hyp(k) - e_s )
                         rh(k,j,i) = q(k,j,i) / q_s * 100.0_wp                                     &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
!
!--       ta (air temperature)
          CASE ( 'ta' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt+1
                      ta(k,j,i) = exner(k) * pt(k,j,i) - degc_to_k
                   ENDDO
                ENDDO
             ENDDO
!
!--       Calculate 'turbulence intensity' from rot[(u,v,w)] at scalar grid point
          CASE ( 'ti' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      ti(k,j,i) = 0.25_wp * SQRT(                                                  &
                                       (   (   w(k,j+1,i) + w(k-1,j+1,i)                           &
                                             - w(k,j-1,i) - w(k-1,j-1,i) ) * ddy                   &
                                         - (   v(k+1,j,i) + v(k+1,j+1,i)                           &
                                             - v(k-1,j,i) - v(k-1,j+1,i) ) * ddzu(k) )**2          &
                                     + (   (   u(k+1,j,i) + u(k+1,j,i+1)                           &
                                             - u(k-1,j,i) - u(k-1,j,i+1) ) * ddzu(k)               &
                                         - (   w(k,j,i+1) + w(k-1,j,i+1)                           &
                                             - w(k,j,i-1) - w(k-1,j,i-1) ) * ddx     )**2          &
                                     + (   (   v(k,j,i+1) + v(k,j+1,i+1)                           &
                                             - v(k,j,i-1) - v(k,j+1,i-1) ) * ddx                   &
                                         - (   u(k,j+1,i) + u(k,j+1,i+1)                           &
                                             - u(k,j-1,i) - u(k,j-1,i+1) ) * ddy     )**2  )       &
                                  * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   ENDDO
                   ti(nzt+1,j,i) = ti(nzt,j,i)
                ENDDO
             ENDDO
!
!--       uu
          CASE ( 'uu_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         uu(k,j,i) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )                            &
                                   * 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )                            &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         uu(k,j,i) = u(k,j,i) * u(k,j,i)                                           &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--       uv
          CASE ( 'uv_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         uv(k,j,i) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )                            &
                                   * 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )                            &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         uv(k,j,i) = u(k,j,i)                                                       &
                                   * 0.25_wp * ( v(k,j,i) + v(k,j+1,i) + v(k,j,i-1) + v(k,j+1,i-1) )&
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--       uw
          CASE ( 'uw_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         uw(k,j,i) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )                            &
                                   * 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )                            &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         uw(k,j,i) = u(k,j,i)                                                       &
                                   * 0.25_wp * ( w(k,j,i) + w(k-1,j,i) + w(k,j,i-1) + w(k-1,j,i-1) )&
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--       vu
          CASE ( 'vu_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         vu(k,j,i) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )                            &
                                   * 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )                            &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         vu(k,j,i) = v(k,j,i)                                                       &
                                   * 0.25_wp * ( u(k,j,i) + u(k,j-1,i) + u(k,j,i+1) + u(k,j-1,i+1) )&
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--       vv
          CASE ( 'vv_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         vv(k,j,i) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )                            &
                                   * 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )                            &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         vv(k,j,i) = v(k,j,i) * v(k,j,i)                                           &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--       vw
          CASE ( 'vw_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         vw(k,j,i) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )                            &
                                   * 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )                            &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         vw(k,j,i) = v(k,j,i)                                                       &
                                   * 0.25_wp * ( w(k,j,i) + w(k-1,j,i) + w(k,j-1,i) + w(k-1,j-1,i) )&
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--       ww
          CASE ( 'ww_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         ww(k,j,i) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )                            &
                                   * 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )                            &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         ww(k,j,i) = w(k,j,i) * w(k,j,i)                                           &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
                      ENDDO
                      ww(nzt+1,j,i) = ww(nzt,j,i)
                   ENDDO
                ENDDO
             ENDIF
!
!--       wu
          CASE ( 'wu_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         wu(k,j,i) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )                            &
                                   * 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )                            &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         wu(k,j,i) = w(k,j,i)                                                         &
                                     * 0.25_wp * ( u(k,j,i) + u(k,j,i+1) + u(k+1,j,i) + u(k+1,j,i+1) )&
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
                      ENDDO
                      wu(nzt+1,j,i) = wu(nzt,j,i)
                   ENDDO
                ENDDO
             ENDIF
!
!--       wv
          CASE ( 'wv_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         wv(k,j,i) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) )                            &
                                   * 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )                            &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         wv(k,j,i) = w(k,j,i)                                                       &
                                   * 0.25_wp * ( v(k,j,i) + v(k,j+1,i) + v(k+1,j,i) + v(k+1,j+1,i) )&
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
                      ENDDO
                      wv(nzt+1,j,i) = wv(nzt,j,i)
                   ENDDO
                ENDDO
             ENDIF
!
!--       wtheta
          CASE ( 'wtheta_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         wtheta(k,j,i) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) ) * pt(k,j,i)            &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                        wtheta(k,j,i) = w(k,j,i) *  0.5_wp  * ( pt(k,j,i) + pt(k+1,j,i) )          &
                                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
                      ENDDO
                      wtheta(nzt+1,j,i) = wtheta(nzt,j,i)
                   ENDDO
                ENDDO
             ENDIF
!
!--       wq
          CASE ( 'wq_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         wq(k,j,i) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) ) * q(k,j,i)                 &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         wq(k,j,i) = w(k,j,i) * 0.5_wp * ( q(k,j,i) + q(k+1,j,i) )                 &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
                      ENDDO
                      wq(nzt+1,j,i) = wq(nzt,j,i)
                   ENDDO
                ENDDO
             ENDIF
!
!--       ws
          CASE ( 'ws_product' )
             IF ( interpolate_to_grid_center )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt+1
                         ws(k,j,i) = 0.5_wp * ( w(k,j,i) + w(k-1,j,i) ) * s(k,j,i)                 &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      ENDDO
                   ENDDO
                ENDDO
             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         ws(k,j,i) = w(k,j,i) * 0.5_wp * ( s(k,j,i) + s(k+1,j,i) )                 &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
                      ENDDO
                      ws(nzt+1,j,i) = ws(nzt,j,i)
                   ENDDO
                ENDDO
             ENDIF
!
!--       2-m water vapor mixing ratio
          CASE ( 'qv_2m*'  )
!
!--          2-m mixing ratio from surface arrays. In case the 2m level is below the first grid
!--          point, MOST is applied, else, linear interpolation between two vertical grid levels
!--          is applied. To access all surfaces, iterate over all horizontally-upward facing
!--          surface types.
             surf => surf_def
             CALL calc_2m_mixing_ratio
             surf => surf_lsm
             CALL calc_2m_mixing_ratio
             surf => surf_usm
             CALL calc_2m_mixing_ratio
!
!--       2-m potential and absolute temperature
          CASE ( 'theta_2m*', 'ta_2m*' )
!
!--          2-m potential temperature is caluclated from surface arrays. In case the 2m level is
!--          below the first grid point, MOST is applied, else, linear interpolation between two
!--          vertical grid levels is applied. To access all surfaces, iterate over all horizontally-
!--          upward facing surface types.
             surf => surf_def
             CALL calc_2m_temperature
             surf => surf_lsm
             CALL calc_2m_temperature
             surf => surf_usm
             CALL calc_2m_temperature
!
!--       10-m wind speed
          CASE ( 'wspeed_10m*' )
!
!--          10-m wind speed is caluclated from surface arrays. In case the 10m level is below the
!--          first grid point, MOST is applied, else, linear interpolation between two vertical grid
!--          levels is applied. To access all surfaces, iterate over all horizontally-upward facing
!--          surface types.
             surf => surf_def
             CALL calc_wind_10m
             surf => surf_lsm
             CALL calc_wind_10m
             surf => surf_usm
             CALL calc_wind_10m
!
!--       Horizontal wind speed
          CASE ( 'wspeed' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      wspeed(k,j,i) = SQRT( ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) )**2              &
                                          + ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) )**2 )            &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   ENDDO
                ENDDO
             ENDDO

!
!--       Horizontal wind direction
          CASE ( 'wdir' )
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      u_center(k,j,i) = 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )
                      v_center(k,j,i) = 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )

                      wdir(k,j,i) = ATAN2( u_center(k,j,i), v_center(k,j,i) )                      &
                                    / pi * 180.0_wp + 180.0_wp
                   ENDDO
                ENDDO
             ENDDO

       END SELECT

       ivar = ivar + 1
    ENDDO

!     CALL cpu_log( log_point(41), 'calculate_quantities', 'stop' )

!
!-- The following block contains subroutines to calculate diagnostic quantities.
 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Find k-index of the given height level above ground
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE find_k_index( level )

    REAL(wp) ::  level !< given height level
!
!-- Search for the k-index where the zw grid is still small/equal the given height level.
    k = nzb
    DO WHILE ( zw(k) <= level )
       k_max = k
       k = k + 1
    ENDDO

 END SUBROUTINE find_k_index


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Function to compare two angles (in degrees) including rollover.
!--------------------------------------------------------------------------------------------------!
 FUNCTION compare_angles( dir1, dir2 )

    REAL(wp) ::  compare_angles  !< difference between angle 1 and 2 in degrees
    REAL(wp) ::  dir1            !< angle 1 given in degrees
    REAL(wp) ::  dir2            !< angle 2 given in degrees


    compare_angles = MIN( 360.0_wp - ABS( dir1 - dir2 ), ABS( dir1 - dir2 ) )

 END FUNCTION compare_angles


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Estimate TKE-production by shear
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_shear( shear_prod )

    REAL(wp), DIMENSION(nzb+1:nzt) ::  shear_prod !< estimated shear production


    DO  k = nzb+1, nzt
       dudx =           ( u(k,j,i+1) - u(k,j,i)                                 ) * ddx
       dudy = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) - u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
       dudz = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) - u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)

       dvdx = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) - v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
       dvdy =           ( v(k,j+1,i) - v(k,j,i)                                 ) * ddy
       dvdz = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) - v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)

       dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) - w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
       dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) - w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
       dwdz =           ( w(k,j,i)   - w(k-1,j,i)                               ) * ddzw(k)
!
!--    Resolved scale fluxes u'_i u'_j. Note, this computation is not exact and
!--    is thought to be only an estimate for the shear production.
!--    u'_i u'_j is approximated by -K_(m,resolved)*(du_i/dx_j + du_j/dx_i), while
!--    K_(m,resolved), the resolved-scale diffusion coefficient is unknown. Hence,
!--    it is simply set to 1. Actually, u'_i u'_j need to be computed via the EC method
!--    (but this is not possible during a run).
       uuf = - 2.0_wp * dudx
       uvf = - ( dudy + dvdx )
       uwf = - ( dudz + dwdx )
       vuf = - ( dvdx + dudy )
       vvf = - 2.0_wp * dvdy
       vwf = - ( dvdz + dwdy )
       wuf = - ( dwdx + dudz )
       wvf = - ( dwdy + dvdz )
       wwf = - 2.0_wp * dwdz
!
!--    Resolved-scale shear production
       shear_prod(k) = - ( uuf * dudx + uvf * dudy + uwf * dudz +                                  &
                           vuf * dvdx + vvf * dvdy + vwf * dvdz +                                  &
                           wuf * dwdx + wvf * dwdy + wwf * dwdz )
!
!--    SGS-shear production. Note, near surface shear production is currently
!--    neglected. This is because the shear-production term is only used for detecting
!--    katabatic flow events where the maximum shear production indicates a jet-like
!--    near-surface flow. To detect the shear induced by the katabatic flow rather than
!--    shear induced by surface friction, shear at wall-bounded grid points is neglected.
       shear_prod(k) = shear_prod(k) + km(k,j,i) *                                                 &
                                       ( 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +                &
                                                    dudy**2 + dvdx**2 + dwdx**2 +                  &
                                                    dwdy**2 + dudz**2 + dvdz**2 +                  &
                                         2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz ) )

       shear_prod(k) = shear_prod(k) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 29 ) )

    ENDDO

 END SUBROUTINE calc_shear


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Estimate TKE-production by shear.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE detect_katabatic_flow( k_detect )

    INTEGER(iwp) ::  k                       !< running index along z-dimension
    INTEGER(iwp), INTENT(OUT) ::  k_detect   !< grid index indicating top of katabatic flow
    INTEGER(iwp) ::  k_dum                   !< grid index along z-dimension to indicate search bounds for shear maximum
    INTEGER(iwp) ::  ks                      !< grid index along z-dimension to indicate first grid point above plant canopy or surface

    LOGICAL ::  jet_found           !< flag to indicate whether a jet-like flow is present
    LOGICAL ::  wind_along_slope    !< flag to indicate whether the actual flow is along the valley slope

    REAL(wp) ::  diff_angles               !< difference between wind direction and terrain slope
    REAL(wp) ::  max_diff_angle = 120.0_wp !< maximum difference between terrain slope angle and wind direction (in degrees)
    REAL(wp) ::  uh_t                      !< horizontal wind speed at k+1
    REAL(wp) ::  uh_b                      !< horizontal wind speed at k-1
    REAL(wp) ::  wind_dir                  !< instantaneous wind-direction


!
!-- Determine starting index (in vertical direction) where to check for katabatic flows.
!-- No search is done within plant canopy as the flow within the trunk area might
!-- be de-coupled from the above-canopy flow.
    IF ( ALLOCATED( pch_index_ji ) )  THEN
       ks = topo_top_ind(j,i,0) + 1 + pch_index_ji(j,i)
    ELSE
       ks = topo_top_ind(j,i,0) + 1
    ENDIF
!
!-- Calculate TKE production by shear and determine index where its maximum
!-- occurs. Therefore, exclude in-canopy jets. Further, the maximum in shear must
!-- coincide with negative vertical shear, i.e. the wind speed must decrease with
!-- height. If this is not the case, go to the next local maximum beyond.
    CALL calc_shear( shear_production )

    jet_found = .FALSE.
    k_max = 1
    k_dum = ks
    DO WHILE ( .NOT. jet_found  .AND.  k_max < MIN(2*nzb_max,nzt) )
       k_max = MAXLOC( shear_production(k_dum:MIN(2*nzb_max,nzt)), DIM = 1 ) + k_dum

       uh_t = SQRT( ( 0.5_wp * ( u(k_max+1,j,i) + u(k_max+1,j,i+1) ) )**2 +                        &
                    ( 0.5_wp * ( v(k_max+1,j,i) + v(k_max+1,j+1,i) ) )**2 )
       uh_b = SQRT( ( 0.5_wp * ( u(k_max-1,j,i) + u(k_max-1,j,i+1) ) )**2 +                        &
                    ( 0.5_wp * ( v(k_max-1,j,i) + v(k_max-1,j+1,i) ) )**2 )

       IF ( ( uh_t - uh_b ) / ( zu(k_max+1) - zu(k_max-1) ) < 0.0_wp )  jet_found = .TRUE.
       k_dum = k_max+1
    ENDDO
!
!-- Check wind direction up to the height of maximum shear production.
!-- The flow needs to be along the terrain slope within a certain range (<120 degrees).
!-- Smaller values are too restrictive and sometimes give a too patchy picture.
!-- Note, this check is not done within any plant canopy. This is because the katabatic
!-- jet could be above the canopy and the in-canopy flow could be uncoupled from this.
    wind_along_slope = .TRUE.
    DO  k = ks, k_max
       wind_dir = ATAN2( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ),                                       &
                         0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) ) / pi * 180.0_wp                      &
                  + 180.0_wp
       diff_angles = compare_angles( wind_dir, slope_angle(j,i) )
       wind_along_slope = wind_along_slope  .AND.  diff_angles <= max_diff_angle
    ENDDO
!
!-- Finally, if a jet-like flow has been detected and the wind direction matches
!-- roughly with the terrain, determine the katabatic flow depth. This is defined
!-- where the shear-production drops to 1/e times its relevant maximum value.
!-- Further, as the katabatic flow cannot be as height as the topography itself, the
!-- considered height level is restricted to 1.5 times nzb_max (maximum topography
!-- within the model domain).
    k_detect = -1
    IF ( shear_production(k_max) /= 0.0_wp  .AND.  wind_along_slope  .AND.  jet_found )  THEN
       DO  k = k_max, INT( 1.5_wp * nzb_max )
!
!--       Check where shear-productions is dropped to 1/e times its maximum. Further,
!--       the katabatic flow cannot be as height as the topography itself.
          IF ( ( shear_production(k) / shear_production(k_max) ) > 0.36_wp  .AND.                  &
               ( zu(k) - zw(topo_top_ind(j,i,0)) ) <= zu(nzb_max) )                                &
          THEN
             k_detect = k
          ENDIF
       ENDDO
    ENDIF

 END SUBROUTINE detect_katabatic_flow


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute volume flux (m3/s, component-wise integrated volume flux for a fixed or variable depth)
!> or volume-flux density (m3/(s m), component-wise integrated volume flux for a fixed or variable
!> depth flow through a 1-m wide section). Component-wise integration means that the x- and y-
!> components are integrated individually. The volume flux then is the absolute value of the
!> vertically integrated components. An exception is made for the lower-half of the first prognostic
!> grid box, where MOST relations are used to compute the volume flux directly.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE integrate_volume_flux( vf, max_integration_height, compute_volume_flux_density )

    INTEGER(iwp) ::  k  !< grid index along z-direction
    INTEGER(iwp) ::  m  !< index for the respective surface element
    INTEGER(iwp) ::  mm !< running index for surfaces

    LOGICAL, OPTIONAL ::  compute_volume_flux_density !< flag to distinguish between volume flux and volume flux density
    LOGICAL           ::  compute_vd                  !< control flag to distinguish between volume flux and volume flux density

    REAL(wp) ::  int_dz                  !< integration width between z0 and the first prognostic grid point
    REAL(wp) ::  level                   !< actual height level
    REAL(wp) ::  max_integration_height  !< maximum integration height
    REAL(wp) ::  vf                      !< integrated volume flux
    REAL(wp) ::  vf_x                    !< volume flux along x
    REAL(wp) ::  vf_y                    !< volume flux along y
    REAL(wp) ::  width_x                 !< integration width along x-direction, either dx or 1
    REAL(wp) ::  width_y                 !< integration width along y-direction, either dy or 1


!
!-- Note, routine is called for each surface type. If the called surface type is not present
!-- at (j,i), return. This is to assure that the integration is only done once.
    IF ( .NOT.  surf%start_index(j,i) <= surf%end_index(j,i) )  RETURN
!
!-- Set component-wise volume-flux contributions to zero
    vf_x = 0.0_wp
    vf_y = 0.0_wp
!
!-- Initialize control flag to distinguish between volumen flux and volume-flux density
    compute_vd = .FALSE.
!
!-- Set control flag to distinguish between volume flux and volume-flux density
    IF ( PRESENT( compute_volume_flux_density ) )  compute_vd = compute_volume_flux_density
!
!-- Pre-calculate horizontal width (in x and y). This is required to distinguish between
!-- volume fluxes and volume-flux densities. width_x is used to integrate the y-component,
!-- width_y to integrate the x-component.
    width_x = MERGE( 1.0_wp, dx, compute_vd )
    width_y = MERGE( 1.0_wp, dy, compute_vd )
!
!-- Vertical integration of the wind velocity. The first grid box above the surface
!-- is treated separately. There, integration is performed from z0 to z_mo via MOST
!-- relations within small vertical intervals (layer between z0 and z_mo is devided into 10
!-- vertical layers), and from z_mo to z_surface + dz.
!-- Integrate over the lower half of the first prognostic grid box.
!-- First, search for the upward-facing surface.
    m = -HUGE( 1 )
    DO  mm = surf%start_index(j,i), surf%end_index(j,i)
       IF ( surf%upward(mm) )  m = mm
    ENDDO
!
!-- Exit the routine if there is no upward-facing surface of the given type. This may happen
!-- if a terrain step is right next to a building, where vertical LSM surfaces are present but
!-- no upward-facing surface. As volume fluxes are only evaluated above upward-facing
!-- LSM surfaces (and default surfaces if present), the subroutine needs to be exited.
    IF ( m == -HUGE( 1 ) )  RETURN

    int_dz = ( surf%z_mo(m) - surf%z0(m) ) * 0.1_wp
    level = surf%z0(m) + int_dz * 0.5_wp
    DO WHILE ( level <= surf%z_mo(m)  .AND.  level <= max_integration_height )
!
!--    Integrate volume flux. Therefore, compute uv(z=level) via MOST.
!--    In contrast to the following integration of the u- and v-components as a vector
!--    integration, the volume flux within the lowest layer is intregrated as a scalar.
       vf = vf + surf%us(m) / kappa                                                                &
               * ( LOG( level / surf%z0(m) ) - psi_m( level      / surf%ol(m) )                    &
                                             + psi_m( surf%z0(m) / surf%ol(m) ) )                  &
               * int_dz * width_x
       level = level + int_dz
    ENDDO
!
!-- Integrate over the upper half of the first prognostic grid box. Consider two cases:
!-- (i) the integration height is beyond the grid-cell top ( integration interval is
!-- 0.5 * dzw(k) ), and (ii) the integration height lies within the grid cell
!-- ( integration interval is 0.5 * (max_integration_height - dzw(k) ) ).
    k = topo_top_ind(j,i,0) + 1
    vf_x = vf_x + 0.5_wp * ( u(k,j,i) + u(k,j,i+1) )                                               &
                         * 0.5_wp * MIN( dzw(k), max_integration_height - dzw(k) ) * width_y
    vf_y = vf_y + 0.5_wp * ( v(k,j,i) + v(k,j+1,i) )                                               &
                         * 0.5_wp * MIN( dzw(k), max_integration_height - dzw(k) ) * width_x
!
!-- Integrate over the remaining grid boxes until k_max is reached.
    DO k = topo_top_ind(j,i,0) + 2, MIN( topo_top_ind(j,i,0) + k_max, nzt + 1 )
       vf_x = vf_x + 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) * dzw(k) * width_y
       vf_y = vf_y + 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) * dzw(k) * width_x
    ENDDO
!
!-- k_max indicates the uppermost grid box that lies fully within the integration height.
!-- Check if the final integration height has been already reached, else, also
!-- integrate over the remaining depth, i.e. within the grid cell at k_max + 1.
!-- Furthermore, note that this branch is not entered if
!-- zw(topo_top_ind(j,i,0)+k_max) - zw(topo_top_ind(j,i,0)) = 0 -- this case (integration
!-- height lies within the 1st prognostic grid box) has been treated already above.
    IF ( ( zw(topo_top_ind(j,i,0)+k_max) - zw(topo_top_ind(j,i,0)) ) <= max_integration_height     &
         .AND.  ( zw(topo_top_ind(j,i,0)+k_max) - zw(topo_top_ind(j,i,0)) ) /= 0 )                 &
    THEN

       k = MIN( topo_top_ind(j,i,0) + k_max + 1, nzt + 1 )
       vf_x = vf_x + 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) * width_y                                  &
                   * ( max_integration_height                                                      &
                     - ( zw(topo_top_ind(j,i,0)+k_max) - zw(topo_top_ind(j,i,0)) ) )
       vf_y = vf_y + 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) * width_x                                  &
                   * ( max_integration_height                                                      &
                     - ( zw(topo_top_ind(j,i,0)+k_max) - zw(topo_top_ind(j,i,0)) ) )

    ENDIF
!
!-- Finally, calculate component-wise integrated volume flux or volume flux density.
    vf = vf + SQRT( vf_x**2 + vf_y**2 )

 END SUBROUTINE integrate_volume_flux


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of 2-m mixing ratio.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_2m_mixing_ratio

    IMPLICIT NONE

    INTEGER(iwp) ::  kk         !< running index along the z-dimension
    INTEGER(iwp) ::  m          !< running index for surface elements

    REAL(wp) ::  qv   !< dummy value for water vapor mixing ratio at index kk
    REAL(wp) ::  qv_m !< dummy value for water vapor mixing ratio at index kk-1


    DO  m = 1, surf%ns
!
!--    Only compute 2-m mixing ratio at horizontally upward facing surfaces.
       IF ( surf%upward(m) )  THEN

          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)
!
!--       If 2-m level is below the first grid level, MOST is used for calculation of 2-m
!--       mixing ratio.
          IF ( surf%z_mo(m) > 2.0_wp )  THEN
             qv_2m(j,i) = surf%q_surface(m) + surf%qs(m) / kappa *                                 &
                                             ( LOG( 2.0_wp / surf%z0q(m) ) -                       &
                                               psi_h( 2.0_wp      / surf%ol(m) ) +                 &
                                               psi_h( surf%z0q(m) / surf%ol(m) ) )
!
!--       If 2-m level is above the first grid level, 2-m mixing ratio is linearly interpolated
!--       between the two nearest vertical grid levels. Note, since 2-m mixing ratio is only
!--       computed for horizontal upward-facing surfaces, only a vertical interpolation is
!--       necessary.
          ELSE
!
!--          zw(k-1) defines the height of the surface.
             kk = k
             DO WHILE ( ( zu(kk) - zw(k-1) ) < 2.0_wp  .AND.  kk <= nzt )
                kk = kk + 1
             ENDDO
!
!--          kk defines the index of the first grid level >= 2m.
!--          To obtain the water vapor mixing ratio, remove the liquid water portion in case
!--          cloud physics are considered.
             IF ( bulk_cloud_model )  THEN
                qv_m  = q(kk-1,j,i) - ql(kk-1,j,i)
                qv    = q(kk,j,i)   - ql(kk,j,i)
             ELSE
                qv_m  = q(kk-1,j,i)
                qv    = q(kk,j,i)
             ENDIF

             qv_2m(j,i) = qv_m + ( qv - qv_m ) * ( zw(k-1) + 2.0_wp - zu(kk-1) ) *                 &
                                                 ( zu(kk)           - zu(kk-1) )
          ENDIF
       ENDIF

    ENDDO

 END SUBROUTINE calc_2m_mixing_ratio


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of 2-m potential and absolute temperature.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_2m_temperature

    IMPLICIT NONE

    INTEGER(iwp) ::  kk         !< running index along the z-dimension
    INTEGER(iwp) ::  m          !< running index for surface elements

    REAL(wp)     ::  exner_2m   !< exner value for 2m about horizontal surface


    DO  m = 1, surf%ns
!
!--    Only compute 2-m temperature at horizontally upward facing surfaces.
       IF ( surf%upward(m) )  THEN

          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)
!
!--       If 2-m level is below the first grid level, MOST is used for calculation of
!--       2-m temperature.
          IF ( surf%z_mo(m) > 2.0_wp )  THEN
             pt_2m(j,i) = surf%pt_surface(m) + surf%ts(m) / kappa *                                &
                                               ( LOG( 2.0_wp / surf%z0h(m) ) -                     &
                                                 psi_h( 2.0_wp      / surf%ol(m) ) +               &
                                                 psi_h( surf%z0h(m) / surf%ol(m) ) )
!
!--       If 2-m level is above the first grid level, 2-m temperature is linearly interpolated
!--       between the two nearest vertical grid levels. Note, since 2-m temperature is only
!--       computed for horizontal upward-facing surfaces, only a vertical interpolation is
!--       necessary.
          ELSE
!
!--          zw(k-1) defines the height of the surface.
             kk = k
             DO WHILE ( ( zu(kk) - zw(k-1) ) < 2.0_wp  .AND.  kk <= nzt )
                kk = kk + 1
             ENDDO
!
!--          kk defines the index of the first grid level >= 2m.
             pt_2m(j,i) = pt(kk-1,j,i) + ( pt(kk,j,i) - pt(kk-1,j,i) ) *                           &
                                         ( zw(k-1) + 2.0_wp - zu(kk-1) ) /                         &
                                         ( zu(kk)           - zu(kk-1) )
          ENDIF

!
!--       Calculate absolute temperature with exner values 2m above surface at zw(k-1)
!--       and convert to degree_C.
          exner_2m = exner_function( barometric_formula( zw(k-1),                                  &
                                     pt_surface * exner_function( surface_pressure * 100.0_wp ),   &
                                     surface_pressure * 100.0_wp )                                 &
                                   )
          ta_2m(j,i) = exner_2m * pt_2m(j,i) - degc_to_k
       ENDIF

    ENDDO

 END SUBROUTINE calc_2m_temperature


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of 10-m wind speed.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_wind_10m

    INTEGER(iwp) ::  kk     !< running index along the z-dimension
    INTEGER(iwp) ::  m      !< running index for surface elements

    REAL(wp) ::  uv_l !< wind speed at lower grid point
    REAL(wp) ::  uv_u !< wind speed at upper grid point


    DO  m = 1, surf%ns
!
!--    Only compute 10-m wind speed at horizontally upward facing surfaces.
       IF ( surf%upward(m) )  THEN

          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)
!
!--       If 10-m level is below the first grid level, MOST is used for calculation of 10-m
!--       wind speed.
          IF ( surf%z_mo(m) > 10.0_wp )  THEN
             uv_10m(j,i) = surf%us(m) / kappa * ( LOG( 10.0_wp /  surf%z0(m) ) -                   &
                                                  psi_m( 10.0_wp    / surf%ol(m) ) +               &
                                                  psi_m( surf%z0(m) / surf%ol(m) ) )
!
!--       If 10-m level is above the first grid level, 10-m wind speed is linearly interpolated
!--       between the two nearest vertical grid levels. Note, since 10-m wind speed is only
!--       computed for horizontal upward-facing surfaces, only a vertical interpolation is
!--       necessary.
          ELSE
!
!--          zw(k-1) defines the height of the surface.
             kk = k
             DO WHILE ( ( zu(kk) - zw(k-1) ) < 10.0_wp  .AND.  kk <= nzt )
                kk = kk + 1
             ENDDO
!
!--          kk defines the index of the first grid level >= 10m.
             uv_l = SQRT( ( 0.5_wp * ( u(kk-1,j,i) + u(kk-1,j,i+1) ) )**2 +                        &
                          ( 0.5_wp * ( v(kk-1,j,i) + v(kk-1,j+1,i) ) )**2 )

             uv_u = SQRT( ( 0.5_wp * ( u(kk,j,i)   + u(kk,j,i+1)   ) )**2 +                        &
                          ( 0.5_wp * ( v(kk,j,i)   + v(kk,j+1,i)   ) )**2 )

             uv_10m(j,i) = uv_l + ( uv_u - uv_l ) * ( zw(k-1) + 10.0_wp - zu(kk-1) ) /             &
                                                    ( zu(kk)            - zu(kk-1) )

          ENDIF
       ENDIF

    ENDDO

 END SUBROUTINE calc_wind_10m

 END SUBROUTINE doq_calculate


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Preparation of the diagnostic output, counting of the module-specific output quantities and
!> gathering of the output names.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_prepare

    IMPLICIT NONE

    CHARACTER (LEN=varnamelength), DIMENSION(0:1,500) ::  do2d_var = ' '  !< label array for 2d output quantities

    INTEGER(iwp) ::  av         !< index defining type of output, av=0 instantaneous, av=1 averaged
    INTEGER(iwp) ::  ivar       !< loop index
    INTEGER(iwp) ::  ivar_all   !< loop index
    INTEGER(iwp) ::  l          !< index for cutting string
    INTEGER(iwp) ::  mid        !< masked output running index


    prepared_diagnostic_output_quantities = .FALSE.

    ivar     = 1
    ivar_all = 1

    DO  av = 0, 1
!
!--    Remove _xy, _xz, or _yz from string
       l = MAX( 3, LEN_TRIM( do2d(av,ivar) ) )
       do2d_var(av,ivar)(1:l-3) = do2d(av,ivar)(1:l-3)
!
!--    Gather 2d output quantity names.
!--    Check for double occurrence of output quantity, e.g. by _xy, _yz, _xz.
       DO  WHILE ( do2d_var(av,ivar)(1:1) /= ' ' )
          IF ( .NOT.  ANY( do_all == do2d_var(av,ivar) ) )  THEN
             do_all(ivar_all) = do2d_var(av,ivar)
          ENDIF
          ivar = ivar + 1
          ivar_all = ivar_all + 1
          l = MAX( 3, LEN_TRIM( do2d(av,ivar) ) )
          do2d_var(av,ivar)(1:l-3) = do2d(av,ivar)(1:l-3)
       ENDDO

       ivar = 1
!
!--    Gather 3d output quantity names
       DO  WHILE ( do3d(av,ivar)(1:1) /= ' ' )
          do_all(ivar_all) = do3d(av,ivar)
          ivar = ivar + 1
          ivar_all = ivar_all + 1
       ENDDO

       ivar = 1
!
!--    Gather masked output quantity names. Also check for double output e.g. by different masks.
       DO  mid = 1, masks
          DO  WHILE ( domask(mid,av,ivar)(1:1) /= ' ' )
             IF ( .NOT.  ANY( do_all == domask(mid,av,ivar) ) )  THEN
                do_all(ivar_all) = domask(mid,av,ivar)
             ENDIF

             ivar = ivar + 1
             ivar_all = ivar_all + 1
          ENDDO
          ivar = 1
       ENDDO

    ENDDO

    prepared_diagnostic_output_quantities = .TRUE.

 END SUBROUTINE doq_prepare

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine reads local (subdomain) restart data
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync,    &
                               nyn_on_file, nysf, nysc, nys_on_file, tmp_2d, tmp_3d, found )


    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  nxlc            !<
    INTEGER(iwp) ::  nxlf            !<
    INTEGER(iwp) ::  nxl_on_file     !<
    INTEGER(iwp) ::  nxrc            !<
    INTEGER(iwp) ::  nxrf            !<
    INTEGER(iwp) ::  nxr_on_file     !<
    INTEGER(iwp) ::  nync            !<
    INTEGER(iwp) ::  nynf            !<
    INTEGER(iwp) ::  nyn_on_file     !<
    INTEGER(iwp) ::  nysc            !<
    INTEGER(iwp) ::  nysf            !<
    INTEGER(iwp) ::  nys_on_file     !<

    LOGICAL, INTENT(OUT)  :: found

    REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_2d   !<

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !<


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'div_new_av' )
          IF ( .NOT. ALLOCATED( div_new_av ) )  THEN
             ALLOCATE( div_new_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          div_new_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                  &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'div_old_av' )
          IF ( .NOT. ALLOCATED( div_old_av ) )  THEN
             ALLOCATE( div_old_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          div_old_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                  &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'pt_2m_av' )
          IF ( .NOT. ALLOCATED( pt_2m_av ) )  THEN
             ALLOCATE( pt_2m_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          pt_2m_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                      &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'qv_2m_av' )
          IF ( .NOT. ALLOCATED( qv_2m_av ) )  THEN
             ALLOCATE( qv_2m_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          qv_2m_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                      &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'ta_2m_av' )
          IF ( .NOT. ALLOCATED( ta_2m_av ) )  THEN
             ALLOCATE( ta_2m_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          ta_2m_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                      &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rh_av' )
          IF ( .NOT. ALLOCATED( rh_av ) )  THEN
             ALLOCATE( rh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          rh_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'ta_av' )
          IF ( .NOT. ALLOCATED( ta_av ) )  THEN
             ALLOCATE( ta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          ta_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'ti_av' )
          IF ( .NOT. ALLOCATED( ti_av ) )  THEN
             ALLOCATE( ti_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          ti_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'u_center_av' )
          IF ( .NOT. ALLOCATED( u_center_av ) )  THEN
             ALLOCATE( u_center_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          u_center_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                 &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'uu_av' )
          IF ( .NOT. ALLOCATED( uu_av ) )  THEN
             ALLOCATE( uu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          uu_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'uv_10m_av' )
          IF ( .NOT. ALLOCATED( uv_10m_av ) )  THEN
             ALLOCATE( uv_10m_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_2d
          uv_10m_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =                                    &
                       tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'v_center_av' )
          IF ( .NOT. ALLOCATED( v_center_av ) )  THEN
             ALLOCATE( v_center_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          v_center_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                 &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'vv_av' )
          IF ( .NOT. ALLOCATED( vv_av ) )  THEN
             ALLOCATE( vv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          vv_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wtheta_av' )
          IF ( .NOT. ALLOCATED( wtheta_av ) )  THEN
             ALLOCATE( wtheta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          wtheta_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wq_av' )
          IF ( .NOT. ALLOCATED( wq_av ) )  THEN
             ALLOCATE( wq_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          wq_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'ws_av' )
          IF ( .NOT. ALLOCATED( ws_av ) )  THEN
             ALLOCATE( ws_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          ws_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wspeed_av' )
          IF ( .NOT. ALLOCATED( wspeed_av ) )  THEN
             ALLOCATE( wspeed_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          wspeed_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wu_av' )
          IF ( .NOT. ALLOCATED( wu_av ) )  THEN
             ALLOCATE( wu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          wu_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wv_av' )
          IF ( .NOT. ALLOCATED( wv_av ) )  THEN
             ALLOCATE( wv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          wv_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'ww_av' )
          IF ( .NOT. ALLOCATED( ww_av ) )  THEN
             ALLOCATE( ww_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
          ww_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)


       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE doq_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_rrd_local_mpi

    LOGICAL ::  array_found  !< control flad indicating whether the array is found in the file or not

!
!-- Restart input of time-averaged quantities is skipped in case of cyclic-fill initialization.
!-- This case, input of time-averaged data is useless and can lead to faulty averaging.
    IF ( .NOT. cyclic_fill_initialization )  THEN
       CALL rd_mpi_io_check_array( 'div_new_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( div_new_av ) )                                                     &
             ALLOCATE( div_new_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'div_new_av', div_new_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'div_old_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( div_old_av ) )                                                     &
             ALLOCATE( div_old_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'div_old_av', div_old_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'kfd_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( kfd_av ) )  ALLOCATE( kfd_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'kfd_av', kfd_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'count_kfd' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( count_kfd ) )  ALLOCATE( count_kfd(nys:nyn,nxl:nxr) )
          CALL rrd_mpi_io( 'count_kfd', count_kfd )
       ENDIF

       CALL rd_mpi_io_check_array( 'hr_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( hr_av ) )  ALLOCATE( hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'hr_av', hr_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'pt_2m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( pt_2m_av ) )  ALLOCATE( pt_2m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'pt_2m_av', pt_2m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'qv_2m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( qv_2m_av ) )  ALLOCATE( qv_2m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'qv_2m_av', qv_2m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'ta_2m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( ta_2m_av ) )  ALLOCATE( ta_2m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'ta_2m_av', ta_2m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'ta_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( ta_av ) )  ALLOCATE( ta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'ta_av', ta_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'ti_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( ti_av ) )  ALLOCATE( ti_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'ti_av', ti_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'rh_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rh_av ) )  ALLOCATE( rh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rh_av', rh_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'u_center_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( u_center_av ) )                                                    &
             ALLOCATE( u_center_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'u_center_av', u_center_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'uu_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( uu_av ) )  ALLOCATE( uu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'uu_av', uu_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'uv_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( uv_av ) )  ALLOCATE( uv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'uv_av', uv_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'uw_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( uw_av ) )  ALLOCATE( uw_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'uw_av', uw_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'uv_10m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( uv_10m_av ) )  ALLOCATE( uv_10m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'uv_10m_av', uv_10m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'v_center_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( v_center_av ) )                                                    &
             ALLOCATE( v_center_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'v_center_av', v_center_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vu_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vu_av ) )  ALLOCATE( vu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vu_av', vu_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vf_25m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vf_25m_av ) )  ALLOCATE( vf_25m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vf_25m_av', vf_25m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vf_50m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vf_50m_av ) )  ALLOCATE( vf_50m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vf_50m_av', vf_50m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vf_75m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vf_75m_av ) )  ALLOCATE( vf_75m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vf_75m_av', vf_75m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vf_100m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vf_100m_av ) )  ALLOCATE( vf_100m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vf_100m_av', vf_100m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vf_xxm_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vf_xxm_av ) )  ALLOCATE( vf_xxm_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vf_xxm_av', vf_xxm_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vfd_25m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vfd_25m_av ) )  ALLOCATE( vfd_25m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vfd_25m_av', vfd_25m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vfd_50m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vfd_50m_av ) )  ALLOCATE( vfd_50m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vfd_50m_av', vfd_50m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vfd_75m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vfd_75m_av ) )  ALLOCATE( vfd_75m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vfd_75m_av', vfd_75m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vfd_100m_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vfd_100m_av ) )  ALLOCATE( vfd_100m_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vfd_100m_av', vfd_100m_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vfd_xxm_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vfd_xxm_av ) )  ALLOCATE( vfd_xxm_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vfd_xxm_av', vfd_xxm_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vv_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vv_av ) )  ALLOCATE( vv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vv_av', vv_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'vw_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( vw_av ) )  ALLOCATE( vw_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'vw_av', vw_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'wu_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wu_av ) )  ALLOCATE( wu_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wu_av', wu_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'wv_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wv_av ) )  ALLOCATE( wv_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wv_av', wv_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'ww_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( ww_av ) )  ALLOCATE( ww_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'ww_av', ww_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'wtheta_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wtheta_av ) )  ALLOCATE( wtheta_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wtheta_av', wtheta_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'wq_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wq_av ) )  ALLOCATE( wq_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wq_av', wq_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'ws_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( ws_av ) )  ALLOCATE( ws_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'ws_av', ws_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'wspeed_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( wspeed_av ) )  ALLOCATE( wspeed_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'wspeed_av', wspeed_av )
       ENDIF

    ENDIF

 END SUBROUTINE doq_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes local (subdomain) restart data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE doq_wrd_local

    IMPLICIT NONE


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       IF ( ALLOCATED( div_new_av ) )  THEN
          CALL wrd_write_string( 'div_new_av' )
          WRITE ( 14 )  div_new_av
       ENDIF

       IF ( ALLOCATED( div_old_av ) )  THEN
          CALL wrd_write_string( 'div_old_av' )
          WRITE ( 14 )  div_old_av
       ENDIF

       IF ( ALLOCATED( pt_2m_av ) )  THEN
          CALL wrd_write_string( 'pt_2m_av' )
          WRITE ( 14 )  pt_2m_av
       ENDIF

       IF ( ALLOCATED( qv_2m_av ) )  THEN
          CALL wrd_write_string( 'qv_2m_av' )
          WRITE ( 14 )  qv_2m_av
       ENDIF

       IF ( ALLOCATED( ta_2m_av ) )  THEN
          CALL wrd_write_string( 'ta_2m_av' )
          WRITE ( 14 )  ta_2m_av
       ENDIF

       IF ( ALLOCATED( rh_av ) )  THEN
          CALL wrd_write_string( 'rh_av' )
          WRITE ( 14 )  rh_av
       ENDIF

       IF ( ALLOCATED( ta_av ) )  THEN
          CALL wrd_write_string( 'ta_av' )
          WRITE ( 14 )  ta_av
       ENDIF

       IF ( ALLOCATED( ti_av ) )  THEN
          CALL wrd_write_string( 'ti_av' )
          WRITE ( 14 )  ti_av
       ENDIF

       IF ( ALLOCATED( u_center_av ) )  THEN
          CALL wrd_write_string( 'u_center_av' )
          WRITE ( 14 )  u_center_av
       ENDIF

       IF ( ALLOCATED( uu_av ) )  THEN
          CALL wrd_write_string( 'uu_av' )
          WRITE ( 14 )  uu_av
       ENDIF

       IF ( ALLOCATED( uv_10m_av ) )  THEN
          CALL wrd_write_string( 'uv_10m_av' )
          WRITE ( 14 )  uv_10m_av
       ENDIF

       IF ( ALLOCATED( v_center_av ) )  THEN
          CALL wrd_write_string( 'v_center_av' )
          WRITE ( 14 )  v_center_av
       ENDIF

       IF ( ALLOCATED( vv_av ) )  THEN
          CALL wrd_write_string( 'vv_av' )
          WRITE ( 14 )  vv_av
       ENDIF

       IF ( ALLOCATED( ww_av ) )  THEN
          CALL wrd_write_string( 'ww_av' )
          WRITE ( 14 )  ww_av
       ENDIF

       IF ( ALLOCATED( wu_av ) )  THEN
          CALL wrd_write_string( 'wu_av' )
          WRITE ( 14 )  wu_av
       ENDIF

       IF ( ALLOCATED( wv_av ) )  THEN
          CALL wrd_write_string( 'wv_av' )
          WRITE ( 14 )  wv_av
       ENDIF

       IF ( ALLOCATED( wtheta_av ) )  THEN
          CALL wrd_write_string( 'wtheta_av' )
          WRITE ( 14 )  wtheta_av
       ENDIF

       IF ( ALLOCATED( wq_av ) )  THEN
          CALL wrd_write_string( 'wq_av' )
          WRITE ( 14 )  wq_av
       ENDIF

       IF ( ALLOCATED( ws_av ) )  THEN
          CALL wrd_write_string( 'ws_av' )
          WRITE ( 14 )  ws_av
       ENDIF

       IF ( ALLOCATED( wspeed_av ) )  THEN
          CALL wrd_write_string( 'wspeed_av' )
          WRITE ( 14 )  wspeed_av
       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       IF ( ALLOCATED( div_new_av ) )   CALL wrd_mpi_io( 'div_new_av', div_new_av )
       IF ( ALLOCATED( div_old_av ) )   CALL wrd_mpi_io( 'div_old_av', div_old_av )
       IF ( ALLOCATED( kfd_av ) )       CALL wrd_mpi_io( 'kfd_av', kfd_av )
       IF ( ALLOCATED( count_kfd ) )    CALL wrd_mpi_io( 'count_kfd', count_kfd )
       IF ( ALLOCATED( hr_av ) )        CALL wrd_mpi_io( 'hr_av', hr_av )
       IF ( ALLOCATED( pt_2m_av ) )     CALL wrd_mpi_io( 'pt_2m_av', pt_2m_av )
       IF ( ALLOCATED( qv_2m_av ) )     CALL wrd_mpi_io( 'qv_2m_av', qv_2m_av )
       IF ( ALLOCATED( ta_2m_av ) )     CALL wrd_mpi_io( 'ta_2m_av', ta_2m_av )
       IF ( ALLOCATED( rh_av ) )        CALL wrd_mpi_io( 'rh_av', rh_av )
       IF ( ALLOCATED( ta_av ) )        CALL wrd_mpi_io( 'ta_av', ta_av )
       IF ( ALLOCATED( ti_av ) )        CALL wrd_mpi_io( 'ti_av', ti_av )
       IF ( ALLOCATED( u_center_av ) )  CALL wrd_mpi_io( 'u_center_av', u_center_av )
       IF ( ALLOCATED( uu_av ) )        CALL wrd_mpi_io( 'uu_av', uu_av )
       IF ( ALLOCATED( uv_av ) )        CALL wrd_mpi_io( 'uv_av', uv_av )
       IF ( ALLOCATED( uv_10m_av ) )    CALL wrd_mpi_io( 'uv_10m_av', uv_10m_av )
       IF ( ALLOCATED( uw_av ) )        CALL wrd_mpi_io( 'uw_av', uw_av )
       IF ( ALLOCATED( vu_av ) )        CALL wrd_mpi_io( 'vu_av', vu_av )
       IF ( ALLOCATED( vv_av ) )        CALL wrd_mpi_io( 'vv_av', vv_av )
       IF ( ALLOCATED( vw_av ) )        CALL wrd_mpi_io( 'vw_av', vw_av )
       IF ( ALLOCATED( v_center_av ) )  CALL wrd_mpi_io( 'v_center_av', v_center_av )
       IF ( ALLOCATED( vf_25m_av ) )    CALL wrd_mpi_io( 'vf_25m_av', vf_25m_av )
       IF ( ALLOCATED( vf_50m_av ) )    CALL wrd_mpi_io( 'vf_50m_av', vf_50m_av )
       IF ( ALLOCATED( vf_75m_av ) )    CALL wrd_mpi_io( 'vf_75m_av', vf_75m_av )
       IF ( ALLOCATED( vf_100m_av ) )   CALL wrd_mpi_io( 'vf_100m_av', vf_100m_av )
       IF ( ALLOCATED( vf_xxm_av ) )    CALL wrd_mpi_io( 'vf_xxm_av', vf_xxm_av )
       IF ( ALLOCATED( vfd_25m_av ) )   CALL wrd_mpi_io( 'vfd_25m_av', vfd_25m_av )
       IF ( ALLOCATED( vfd_50m_av ) )   CALL wrd_mpi_io( 'vfd_50m_av', vfd_50m_av )
       IF ( ALLOCATED( vfd_75m_av ) )   CALL wrd_mpi_io( 'vfd_75m_av', vfd_75m_av )
       IF ( ALLOCATED( vfd_100m_av ) )  CALL wrd_mpi_io( 'vfd_100m_av', vfd_100m_av )
       IF ( ALLOCATED( vfd_xxm_av ) )   CALL wrd_mpi_io( 'vfd_xxm_av', vfd_xxm_av )
       IF ( ALLOCATED( wtheta_av ) )    CALL wrd_mpi_io( 'wtheta_av', wtheta_av )
       IF ( ALLOCATED( wq_av ) )        CALL wrd_mpi_io( 'wq_av', wq_av )
       IF ( ALLOCATED( ws_av ) )        CALL wrd_mpi_io( 'ws_av', ws_av )
       IF ( ALLOCATED( wspeed_av ) )    CALL wrd_mpi_io( 'wspeed_av', wspeed_av )
       IF ( ALLOCATED( wu_av ) )        CALL wrd_mpi_io( 'wu_av', wu_av )
       IF ( ALLOCATED( wv_av ) )        CALL wrd_mpi_io( 'wv_av', wv_av )
       IF ( ALLOCATED( ww_av ) )        CALL wrd_mpi_io( 'ww_av', ww_av )

    ENDIF

 END SUBROUTINE doq_wrd_local

 END MODULE diagnostic_output_quantities_mod
