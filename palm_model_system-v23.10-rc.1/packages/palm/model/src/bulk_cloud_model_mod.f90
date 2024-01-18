!> @file bulk_cloud_model_mod.f90
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
!> Calculate bulk cloud microphysics.
!--------------------------------------------------------------------------------------------------!
 MODULE bulk_cloud_model_mod


    USE advec_s_bc_mod,                                                                            &
        ONLY:  advec_s_bc

    USE advec_s_pw_mod,                                                                            &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                                            &
        ONLY:  advec_s_up

    USE advec_ws,                                                                                  &
        ONLY:  advec_s_ws

    USE arrays_3d,                                                                                 &
        ONLY:  ddzu,                                                                               &
               d_exner,                                                                            &
               diss,                                                                               &
               diss_l_nc,                                                                          &
               diss_l_ng,                                                                          &
               diss_l_ni,                                                                          &
               diss_l_nr,                                                                          &
               diss_l_ns,                                                                          &
               diss_l_qc,                                                                          &
               diss_l_qg,                                                                          &
               diss_l_qi,                                                                          &
               diss_l_qr,                                                                          &
               diss_l_qs,                                                                          &
               diss_s_qc,                                                                          &
               diss_s_qg,                                                                          &
               diss_s_qi,                                                                          &
               diss_s_qr,                                                                          &
               diss_s_qs,                                                                          &
               diss_s_nc,                                                                          &
               diss_s_ng,                                                                          &
               diss_s_ni,                                                                          &
               diss_s_nr,                                                                          &
               diss_s_ns,                                                                          &
               dzu,                                                                                &
               dzw,                                                                                &
               exner,                                                                              &
               flux_l_nc,                                                                          &
               flux_l_ng,                                                                          &
               flux_l_ni,                                                                          &
               flux_l_nr,                                                                          &
               flux_l_ns,                                                                          &
               flux_l_qc,                                                                          &
               flux_l_qg,                                                                          &
               flux_l_qi,                                                                          &
               flux_l_qr,                                                                          &
               flux_l_qs,                                                                          &
               flux_s_nc,                                                                          &
               flux_s_ng,                                                                          &
               flux_s_ni,                                                                          &
               flux_s_nr,                                                                          &
               flux_s_ns,                                                                          &
               flux_s_qc,                                                                          &
               flux_s_qg,                                                                          &
               flux_s_qi,                                                                          &
               flux_s_qr,                                                                          &
               flux_s_qs

    USE arrays_3d,                                                                                 &
        ONLY:  hyp,                                                                                &
               hyrho,                                                                              &
               nc,                                                                                 &
               nc_p,                                                                               &
               nc_1,                                                                               &
               nc_2,                                                                               &
               nc_3,                                                                               &
               ng,                                                                                 &
               ng_p,                                                                               &
               ng_1,                                                                               &
               ng_2,                                                                               &
               ng_3,                                                                               &
               ni,                                                                                 &
               ni_p,                                                                               &
               ni_1,                                                                               &
               ni_2,                                                                               &
               ni_3,                                                                               &
               nr,                                                                                 &
               nr_p,                                                                               &
               nr_1,                                                                               &
               nr_2,                                                                               &
               nr_3,                                                                               &
               ns,                                                                                 &
               ns_p,                                                                               &
               ns_1,                                                                               &
               ns_2,                                                                               &
               ns_3,                                                                               &
               precipitation_amount,                                                               &
               prr,                                                                                &
               prr_cloud,                                                                          &
               prr_graupel,                                                                        &
               prr_ice,                                                                            &
               prr_rain,                                                                           &
               prr_snow,                                                                           &
               pt,                                                                                 &
               pt_init,                                                                            &
               q,                                                                                  &
               qc,                                                                                 &
               qc_p,                                                                               &
               qc_1,                                                                               &
               qc_2,                                                                               &
               qc_3,                                                                               &
               qf,                                                                                 &
               qf_1,                                                                               &
               qg,                                                                                 &
               qg_p,                                                                               &
               qg_1,                                                                               &
               qg_2,                                                                               &
               qg_3,                                                                               &
               ql,                                                                                 &
               ql_1,                                                                               &
               qi,                                                                                 &
               qi_p,                                                                               &
               qi_1,                                                                               &
               qi_2,                                                                               &
               qi_3,                                                                               &
               qr,                                                                                 &
               qr_p,                                                                               &
               qr_1,                                                                               &
               qr_2,                                                                               &
               qr_3,                                                                               &
               qs,                                                                                 &
               qs_p,                                                                               &
               qs_1,                                                                               &
               qs_2,                                                                               &
               qs_3,                                                                               &
               rdf_sc,                                                                             &
               tend,                                                                               &
               tnc_m,                                                                              &
               tng_m,                                                                              &
               tni_m,                                                                              &
               tnr_m,                                                                              &
               tns_m,                                                                              &
               tqc_m,                                                                              &
               tqg_m,                                                                              &
               tqi_m,                                                                              &
               tqr_m,                                                                              &
               tqs_m,                                                                              &
               zu

    USE averaging,                                                                                 &
        ONLY:  nc_av,                                                                              &
               ng_av,                                                                              &
               ni_av,                                                                              &
               nr_av,                                                                              &
               ns_av,                                                                              &
               prr_av,                                                                             &
               prr_cloud_av,                                                                       &
               prr_graupel_av,                                                                     &
               prr_ice_av,                                                                         &
               prr_rain_av,                                                                        &
               prr_snow_av,                                                                        &
               qc_av,                                                                              &
               qg_av,                                                                              &
               qi_av,                                                                              &
               ql_av,                                                                              &
               qr_av,                                                                              &
               qs_av

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  barometric_formula,                                                                 &
               c_p,                                                                                &
               c_w,                                                                                &
               exner_function,                                                                     &
               exner_function_invers,                                                              &
               g,                                                                                  &
               ideal_gas_law_rho,                                                                  &
               ideal_gas_law_rho_pt,                                                               &
               lv_d_cp,                                                                            &
               lv_d_rd,                                                                            &
               l_m,                                                                                &
               l_v,                                                                                &
               magnus,                                                                             &
               magnus_ice,                                                                         &
               magnus_tl,                                                                          &
               rd_d_rv,                                                                            &
               l_s,                                                                                &
               ls_d_cp,                                                                            &
               molecular_weight_of_solute,                                                         &
               molecular_weight_of_water,                                                          &
               pi,                                                                                 &
               rho_i,                                                                              &
               rho_l,                                                                              &
               rho_s,                                                                              &
               r_d,                                                                                &
               r_v,                                                                                &
               vanthoff

    USE control_parameters,                                                                        &
        ONLY:  advanced_div_correction,                                                            &
               bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               cloud_droplets,                                                                     &
               cyclic_fill_initialization,                                                         &
               debug_output,                                                                       &
               dt_do2d_xy,                                                                         &
               dt_3d,                                                                              &
               humidity,                                                                           &
               initializing_actions,                                                               &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               large_scale_forcing,                                                                &
               loop_optimization,                                                                  &
               lsf_surf,                                                                           &
               mask_i,                                                                             &
               mask_j,                                                                             &
               mask_k,                                                                             &
               mask_size_l,                                                                        &
               mask_surface,                                                                       &
               message_string,                                                                     &
               pt_surface,                                                                         &
               restart_data_format_output,                                                         &
               rho_surface,                                                                        &
               scalar_advec,                                                                       &
               simulated_time,                                                                     &
               surface_pressure,                                                                   &
               timestep_scheme,                                                                    &
               time_do2d_xy,                                                                       &
               tsc,                                                                                &
               ws_scheme_sca

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE diffusion_s_mod,                                                                           &
        ONLY:  diffusion_s

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  advc_flags_s,                                                                       &
               nbgp,                                                                               &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_top_ind,                                                                       &
               topo_flags

    USE kinds

    USE pegrid,                                                                                    &
        ONLY:  threads_per_task

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array,                                                              &
               rrd_mpi_io,                                                                         &
               wrd_mpi_io

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               sums_wsncs_ws_l,                                                                    &
               sums_wsngs_ws_l,                                                                    &
               sums_wsnis_ws_l,                                                                    &
               sums_wsnrs_ws_l,                                                                    &
               sums_wsnss_ws_l,                                                                    &
               sums_wsqcs_ws_l,                                                                    &
               sums_wsqgs_ws_l,                                                                    &
               sums_wsqis_ws_l,                                                                    &
               sums_wsqrs_ws_l,                                                                    &
               sums_wsqss_ws_l,                                                                    &
               ts_value,                                                                           &
               weight_pres,                                                                        &
               weight_substep

    USE surface_mod,                                                                               &
        ONLY:   bc_hv,                                                                             &
                surf_bulk_cloud_model,                                                             &
                surf_def,                                                                          &
                surf_lsm,                                                                          &
                surf_microphysics_ice_phase,                                                       &
                surf_microphysics_morrison,                                                        &
                surf_microphysics_seifert,                                                         &
                surf_top,                                                                          &
                surf_usm

    IMPLICIT NONE

    CHARACTER (LEN=20)   ::  aerosol_bulk = 'nacl'                        !< namelist parameter
    CHARACTER (LEN=20)   ::  cloud_scheme = 'saturation_adjust'           !< namelist parameter

    INTEGER(iwp) ::  dots_start_index_bcm !< start index for time series of this module

    LOGICAL ::  aerosol_nacl =.TRUE.                         !< nacl aerosol for bulk scheme
    LOGICAL ::  aerosol_c3h4o4 =.FALSE.                      !< malonic acid aerosol for bulk scheme
    LOGICAL ::  aerosol_nh4no3 =.FALSE.                      !< malonic acid aerosol for bulk scheme
    LOGICAL ::  bulk_cloud_model = .FALSE.                   !< namelist parameter
    LOGICAL ::  call_microphysics_at_all_substeps = .FALSE.  !< namelist parameter
    LOGICAL ::  cloud_water_sedimentation = .FALSE.          !< cloud water sedimentation
    LOGICAL ::  collision_turbulence = .FALSE.               !< turbulence effects
    LOGICAL ::  curvature_solution_effects_bulk = .FALSE.    !< flag for considering koehler theory
    LOGICAL ::  ice_crystal_sedimentation = .FALSE.          !< flag for ice crystal sedimentation
    LOGICAL ::  limiter_sedimentation = .TRUE.               !< sedimentation limiter
    LOGICAL ::  microphysics_ice_phase = .FALSE.             !< use ice microphysics scheme
    LOGICAL ::  microphysics_kessler = .FALSE.               !< use kessler bulk scheme?
    LOGICAL ::  microphysics_morrison = .FALSE.              !< use 2-moment Morrison
                                                             !< (add. prog. eq. for nc and qc)
    LOGICAL ::  microphysics_morrison_no_rain = .FALSE.      !< use 2-moment Morrison
    LOGICAL ::  microphysics_sat_adjust = .FALSE.            !< use saturation adjust bulk scheme?
    LOGICAL ::  microphysics_seifert = .FALSE.               !< use 2-moment Seifert and Beheng
                                                             !< scheme
    LOGICAL ::  precipitation = .FALSE.                      !< namelist parameter
    LOGICAL ::  snow = .FALSE.                               !< namelist parameter
    LOGICAL ::  snow_sedimentation = .TRUE.                  !< namelist parameter
    LOGICAL ::  graupel_sedimentation = .TRUE.               !< namelist parameter
    LOGICAL ::  graupel = .FALSE.                            !< namelist parameter
    LOGICAL ::  ventilation_effect = .TRUE.                  !< ventilation effect
    LOGICAL ::  ice_multiplication = .TRUE.                  !< Turn on ice multiplication (Hallet Mossop process)
    LOGICAL ::  enhanced_melting = .TRUE.                    !< Turn on ice enhacement of melting during riming

    REAL(wp), PARAMETER ::  a_1 = 8.69E-4_wp                    !< coef. in turb. parametrization    (cm-2 s3)
    REAL(wp), PARAMETER ::  a_2 = -7.38E-5_wp                   !< coef. in turb. parametrization    (cm-2 s3)
    REAL(wp), PARAMETER ::  a_3 = -1.40E-2_wp                   !< coef. in turb. parametrization
    REAL(wp), PARAMETER ::  a_term = 9.65_wp                    !< coef. for terminal velocity       (m s-1)
    REAL(wp), PARAMETER ::  a_vent = 0.78_wp                    !< coef. for ventilation effect
    REAL(wp), PARAMETER ::  b_1 = 11.45E-6_wp                   !< coef. in turb. parametrization    (m)
    REAL(wp), PARAMETER ::  b_2 = 9.68E-6_wp                    !< coef. in turb. parametrization    (m)
    REAL(wp), PARAMETER ::  b_3 = 0.62_wp                       !< coef. in turb. parametrization
    REAL(wp), PARAMETER ::  b_term = 9.8_wp                     !< coef. for terminal velocity       (m s-1)
    REAL(wp), PARAMETER ::  b_vent = 0.308_wp                   !< coef. for ventilation effect
    REAL(wp), PARAMETER ::  beta_cc = 3.09E-4_wp                !< coef. in turb. parametrization    (cm-2 s3)
    REAL(wp), PARAMETER ::  c_1 = 4.82E-6_wp                    !< coef. in turb. parametrization    (m)
    REAL(wp), PARAMETER ::  c_2 = 4.8E-6_wp                     !< coef. in turb. parametrization    (m)
    REAL(wp), PARAMETER ::  c_3 = 0.76_wp                       !< coef. in turb. parametrization
    REAL(wp), PARAMETER ::  c_const = 0.93_wp                   !< const. in Taylor-microscale Reynolds number
    REAL(wp), PARAMETER ::  c_evap = 0.7_wp                     !< constant in evaporation
    REAL(wp), PARAMETER ::  c_term = 600.0_wp                   !< coef. for terminal velocity       (m-1)
    REAL(wp), PARAMETER ::  diff_coeff_l = 0.23E-4_wp           !< diffusivity of water vapor        (m2 s-1)
    REAL(wp), PARAMETER ::  eps_sb = 1.0E-10_wp                 !< threshold in two-moments scheme
    REAL(wp), PARAMETER ::  eps_sb_coll = 1.0E-7_wp             !< q-threshold for collection 1e-4 g/m3
    REAL(wp), PARAMETER ::  eps_sb_n = 1.0E-20_wp               !< add small psi to number concentration
    REAL(wp), PARAMETER ::  eps_mr = 0.0_wp                     !< threshold for morrison scheme
    REAL(wp), PARAMETER ::  k_cc = 9.44E09_wp                   !< const. cloud-cloud kernel         (m3 kg-2 s-1)
    REAL(wp), PARAMETER ::  k_cr0 = 4.33_wp                     !< const. cloud-rain kernel          (m3 kg-1 s-1)
    REAL(wp), PARAMETER ::  k_rr = 7.12_wp                      !< const. rain-rain kernel           (m3 kg-1 s-1)
    REAL(wp), PARAMETER ::  k_br = 1000.0_wp                    !< const. in breakup parametrization (m-1)
    REAL(wp), PARAMETER ::  k_st = 1.2E8_wp                     !< const. in drizzle parametrization (m-1 s-1)
    REAL(wp), PARAMETER ::  kin_vis_air = 1.4086E-5_wp          !< kin. viscosity of air             (m2 s-1)
    REAL(wp), PARAMETER ::  ql_crit = 0.0005_wp                 !< coef. in Kessler scheme           (kg kg-1)
    REAL(wp), PARAMETER ::  rho_0 = 1.225_wp                    !< reference air densitiy following SB2006
    REAL(wp), PARAMETER ::  schmidt_p_1d3 = 0.8921121_wp        !< Schmidt number**0.33333, 0.71**0.33333
    REAL(wp), PARAMETER ::  t3 = 273.15_wp                      !< temperature of freezing point
    REAL(wp), PARAMETER ::  thermal_conductivity_l = 2.43E-2_wp !< therm. cond. air         (J m-1 s-1 K-1)
    REAL(wp), PARAMETER ::  w_precipitation = 9.65_wp           !< maximum terminal velocity         (m s-1)
    REAL(wp), PARAMETER ::  x0 = 2.6E-10_wp                     !< separating drop mass              (kg)
    REAL(wp), PARAMETER ::  ximin = 4.42E-14_wp                 !< minimum mass of ice crystal       (kg) (D~10µm)
    REAL(wp), PARAMETER ::  ximax = 1.0E-7_wp                   !< maximum rain drop site            (kg)
    REAL(wp), PARAMETER ::  xcmax = 2.6E-10_wp                  !< maximum cloud drop size           (kg)
    REAL(wp), PARAMETER ::  xcmin = 4.18E-15_wp                 !< minimum cloud drop size           (kg) (~ 1µm)
    REAL(wp), PARAMETER ::  xgmax = 1.0E-4_wp                   !< maximum graupel drop size         (kg)
    REAL(wp), PARAMETER ::  xgmin = 2.6E-10_wp                  !< minimum graupel drop size         (kg)
    REAL(wp), PARAMETER ::  xsmax = 1.0E-7_wp                   !< minimum snow size                 (kg)
    REAL(wp), PARAMETER ::  xsmin = 1.73E-9_wp                  !< minimum snow size                 (kg)
    REAL(wp), PARAMETER ::  xrmin = 2.6E-10_wp                  !< minimum rain drop size            (kg)
    REAL(wp), PARAMETER ::  xrmax = 5.0E-6_wp                   !< maximum rain drop site            (kg)
    REAL(wp), PARAMETER ::  q_thres = 1.0E-6_wp                 !< q threshold for riming            (kg/kg)
    REAL(wp), PARAMETER ::  qr_crit = 1.0E-5_wp                 !< q-threshold for riming ice- rain
    REAL(wp), PARAMETER ::  dc_coll = 40.0E-6_wp                !< upper bound for diameter in collision efficiency
    REAL(wp), PARAMETER ::  dig_conv = 200.0E-6_wp              !< diameter threshold for conversion of ice to graupel (m)
    REAL(wp), PARAMETER ::  dsg_conv = 200.0E-6_wp              !< diameter threshold for conversion of snow to graupel (m)
    REAL(wp), PARAMETER ::  dc_crit = 10.0E-6_wp                !< critical diameter of cloud droplets for riming (m)
    REAL(wp), PARAMETER ::  df_crit = 100.0E-6_wp               !< critical diameter of graupel for riming (m)
    REAL(wp), PARAMETER ::  ecoll_min = 0.01_wp                 !< minimum efficences for riming of cloud droplets
    REAL(wp), PARAMETER ::  temp_mult_min = 265.0_wp            !< minimum temperature for splintering (K)
    REAL(wp), PARAMETER ::  temp_mult_max = 270.0_wp            !< maximum temperature for splintering (K)
    REAL(wp), PARAMETER ::  temp_mult_opt = 268.0_wp            !< optimal temperature splintering
    REAL(wp), PARAMETER ::  c_mult = 3.5E8_wp                   !< coeffiecnt for splintering
    REAL(wp), PARAMETER ::  x_conv = 0.1E-9_wp                  !< minimum mass of conversion due to riming
    REAL(wp), PARAMETER ::  alpha_spacefilling = 0.01_wp        !< space filling coefficient (max. 0.68)

    REAL(wp) ::  c_sedimentation = 2.0_wp        !< Courant number of sedimentation process
    REAL(wp) ::  a1_ven_coeff_snow               !< 1st coefficient for ventilation of snow
    REAL(wp) ::  b1_ven_coeff_snow               !< 2nd coefficient for ventilation of snow
    REAL(wp) ::  a1_ven_coeff_graupel            !< 1st coefficient for ventilation of graupel
    REAL(wp) ::  b1_ven_coeff_graupel            !< 2nd coefficient for ventilation of graupel
    REAL(wp) ::  cz_cloud                        !< coefficient for homogeneous freezing of cloud droplets
    REAL(wp) ::  coll_coeff_graupel_self         !< pre-calculated collision coefficient for graupel selfcollection
    REAL(wp) ::  coll_coeff_snow_self_delta      !< pre-calculated collision coefficient for snow selfcollection
    REAL(wp) ::  coll_coeff_snow_self_theta      !< pre-calculated collision coefficient for snow selfcollection
    REAL(wp) ::  coll_coeff_ice_self_delta_n     !< pre-calculated collision coefficient for ice selfcollection
    REAL(wp) ::  coll_coeff_ice_self_delta_q     !< pre-calculated collision coefficient for ice selfcollection
    REAL(wp) ::  coll_coeff_ice_self_theta_n     !< pre-calculated collision coefficient for ice selfcollection
    REAL(wp) ::  coll_coeff_ice_self_theta_q     !< pre-calculated collision coefficient for ice selfcollection
    REAL(wp) ::  coll_graupel_ice_delta_n_aa     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_delta_n_bb     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_delta_n_ab     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_delta_q_aa     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_delta_q_bb     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_delta_q_ab     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_theta_n_aa     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_theta_n_bb     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_theta_n_ab     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_theta_q_aa     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_theta_q_bb     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_graupel_ice_theta_q_ab     !< pre-calculated collision coefficient for graupel-ice collection
    REAL(wp) ::  coll_snow_ice_delta_n_aa        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_delta_n_bb        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_delta_n_ab        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_delta_q_aa        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_delta_q_bb        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_delta_q_ab        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_theta_n_aa        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_theta_n_bb        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_theta_n_ab        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_theta_q_aa        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_theta_q_bb        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_snow_ice_theta_q_ab        !< pre-calculated collision coefficient for snow-ice collection
    REAL(wp) ::  coll_graupel_snow_delta_n_aa    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_delta_n_bb    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_delta_n_ab    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_delta_q_aa    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_delta_q_bb    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_delta_q_ab    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_theta_n_aa    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_theta_n_bb    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_theta_n_ab    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_theta_q_aa    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_theta_q_bb    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  coll_graupel_snow_theta_q_ab    !< pre-calculated collision coefficient for graupel-snow collection
    REAL(wp) ::  rime_graupel_cloud_delta_n_aa   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_delta_n_ab   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_delta_n_bb   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_delta_q_aa   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_delta_q_ab   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_delta_q_ba   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_delta_q_bb   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_theta_n_aa   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_theta_n_ab   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_theta_n_bb   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_theta_q_aa   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_theta_q_ab   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_theta_q_ba   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_cloud_theta_q_bb   !< pre-calculated riming coefficient for graupel-cloud
    REAL(wp) ::  rime_graupel_rain_delta_n_aa    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_delta_n_ab    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_delta_n_bb    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_delta_q_aa    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_delta_q_ab    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_delta_q_ba    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_delta_q_bb    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_theta_n_aa    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_theta_n_ab    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_theta_n_bb    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_theta_q_aa    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_theta_q_ab    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_theta_q_ba    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_graupel_rain_theta_q_bb    !< pre-calculated riming coefficient for graupel-rain
    REAL(wp) ::  rime_ice_cloud_delta_n_aa       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_delta_n_ab       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_delta_n_bb       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_delta_q_aa       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_delta_q_ab       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_delta_q_ba       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_delta_q_bb       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_theta_n_aa       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_theta_n_ab       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_theta_n_bb       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_theta_q_aa       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_theta_q_ab       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_theta_q_ba       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_cloud_theta_q_bb       !< pre-calculated riming coefficient for ice-cloud
    REAL(wp) ::  rime_ice_rain_delta_n_aa        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_delta_n_ab        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_delta_n_bb        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_delta_q_aa        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_delta_q_ab        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_delta_q_ba        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_delta_q_bb        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_theta_n_aa        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_theta_n_ab        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_theta_n_bb        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_theta_q_aa        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_theta_q_ab        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_theta_q_ba        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_ice_rain_theta_q_bb        !< pre-calculated riming coefficient for ice-rain
    REAL(wp) ::  rime_snow_cloud_delta_n_aa      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_delta_n_ab      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_delta_n_bb      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_delta_q_aa      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_delta_q_ab      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_delta_q_bb      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_theta_n_aa      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_theta_n_ab      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_theta_n_bb      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_theta_q_aa      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_theta_q_ab      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_cloud_theta_q_bb      !< pre-calculated riming coefficient for snow-cloud
    REAL(wp) ::  rime_snow_rain_delta_n_aa       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_delta_n_ab       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_delta_n_bb       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_delta_q_aa       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_delta_q_ab       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_delta_q_ba       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_delta_q_bb       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_theta_n_aa       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_theta_n_ab       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_theta_n_bb       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_theta_q_aa       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_theta_q_ab       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_theta_q_ba       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  rime_snow_rain_theta_q_bb       !< pre-calculated riming coefficient for snow-rain
    REAL(wp) ::  dpirho_l                        !< 6.0 / ( pi * rho_l )
    REAL(wp) ::  dry_aerosol_radius = 0.05E-6_wp !< dry aerosol radius
    REAL(wp) ::  dt_micro                        !< microphysics time step
    REAL(wp) ::  dt_precipitation = 100.0_wp     !< timestep precipitation (s)
    REAL(wp) ::  e_s                             !< saturation water vapor pressure
    REAL(wp) ::  e_si                            !< saturation water vapor pressure over ice
    REAL(wp) ::  in_init = 1000.0_wp             !< initial number of potential ice nucleii
    REAL(wp) ::  na_init = 100.0E6_wp            !< Total particle/aerosol concentration (cm-3)
    REAL(wp) ::  nc_const = 70.0E6_wp            !< cloud droplet concentration
    REAL(wp) ::  pirho_l                         !< pi * rho_l / 6.0
    REAL(wp) ::  prec_time_const = 0.001_wp      !< coef. in Kessler scheme           (s-1)
    REAL(wp) ::  precipitation_amount_interval = 9999999.9_wp     !< namelist parameter
    REAL(wp) ::  q_s                             !< saturation mixing ratio
    REAL(wp) ::  q_si                            !< saturation mixing ratio over ice
    REAL(wp) ::  sat                             !< supersaturation
    REAL(wp) ::  sat_ice                         !< supersaturation over ice
    REAL(wp) ::  start_ice_microphysics = 0.0_wp !< time from which on ice microhysics are executed
    REAL(wp) ::  sed_qc_const                    !< const. for sedimentation of cloud water
    REAL(wp) ::  sigma_bulk = 2.0_wp             !< width of aerosol spectrum
    REAL(wp) ::  sigma_gc = 1.3_wp               !< geometric standard deviation cloud droplets
    REAL(wp) ::  t_l                             !< liquid-(ice) water temperature

    TYPE ::  cloud_coefficients                !< coefficients for each cloud species
       REAL(wp) ::  a                          !< a coefficient (pre-factor in diameter-mass relation, Table 1 in SB2006)
       REAL(wp) ::  b                          !< b coefficient (exponent in diameter-mass relation, Table 1 in SB2006)
       REAL(wp) ::  alpha                      !< alpha coefficient (pre-factor in power law fall speed, Table 1 in SB2006)
       REAL(wp) ::  beta                       !< beta coefficient (exponent in power law fall speed, Table 1 in SB2006)
       REAL(wp) ::  nu                         !< nu coefficient (first shape parameter of size distribution, Table 1 in SB2006)
       REAL(wp) ::  mu                         !< mu coefficient (second shape parameter of size distribution, Table 1in SB2006)
       REAL(wp) ::  sigma_v                    !< variance of fall speed (see SB2006, List of symbols)
       REAL(wp) ::  a_ven                      !< a ventilation coefficient
       REAL(wp) ::  b_ven                      !< b ventilation coefficient
       REAL(wp) ::  cap                        !< capacity of particle
       REAL(wp) ::  x_min                      !< minimum (mean) mass of particle in SB-scheme
       REAL(wp) ::  x_max                      !< maximum (mean) mass of particle in SB-scheme
       REAL(wp) ::  coll_eff                   !< maximum collision efficiency with cloud droplets
    END TYPE cloud_coefficients

    TYPE ::  cloud_species_def                  !< definition class for each ice species
       TYPE(cloud_coefficients) ::  cloud       !< cloud
       TYPE(cloud_coefficients) ::  rain        !< rain
       TYPE(cloud_coefficients) ::  graupel     !< graupel
       TYPE(cloud_coefficients) ::  ice         !< ice
       TYPE(cloud_coefficients) ::  snow        !< snow
    END TYPE

    TYPE(cloud_species_def) ::  cloud_species  !< ice species

    REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  dep_rate_ice       !< 3D-field of ice deposition rates
    REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  dep_rate_snow      !< 3D-field of ice deposition rates
    REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  rime_ice_cloud     !< 3D-field of ice-cloud riming rates
    REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  rime_ice_rain      !< 3D-field of ice rain riming rates
    REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  rime_snow_cloud    !< 3D-field of ice-cloud riming rates
    REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  rime_snow_rain     !< 3D-field of ice rain riming rates

    SAVE

    PRIVATE

    PUBLIC bcm_actions,                                                                            &
           bcm_boundary_conditions,                                                                &
           bcm_check_data_output,                                                                  &
           bcm_check_data_output_pr,                                                               &
           bcm_check_data_output_ts,                                                               &
           bcm_check_parameters,                                                                   &
           bcm_data_output_mask,                                                                   &
           bcm_data_output_2d,                                                                     &
           bcm_data_output_3d,                                                                     &
           bcm_exchange_horiz,                                                                     &
           bcm_header,                                                                             &
           bcm_init,                                                                               &
           bcm_init_arrays,                                                                        &
           bcm_non_advective_processes,                                                            &
           bcm_parin,                                                                              &
           bcm_prognostic_equations,                                                               &
           bcm_rrd_global,                                                                         &
           bcm_rrd_local,                                                                          &
           bcm_statistics,                                                                         &
           bcm_swap_timelevel,                                                                     &
           bcm_wrd_global,                                                                         &
           bcm_wrd_local,                                                                          &
           bcm_3d_data_averaging

    PUBLIC bulk_cloud_model,                                                                       &
           calc_liquid_water_content,                                                              &
           call_microphysics_at_all_substeps,                                                      &
           cloud_scheme,                                                                           &
           cloud_water_sedimentation,                                                              &
           collision_turbulence,                                                                   &
           dt_precipitation,                                                                       &
           graupel,                                                                                &
           ice_crystal_sedimentation,                                                              &
           in_init,                                                                                &
           microphysics_morrison,                                                                  &
           microphysics_morrison_no_rain,                                                          &
           microphysics_sat_adjust,                                                                &
           microphysics_seifert,                                                                   &
           microphysics_ice_phase,                                                                 &
           nc_const,                                                                               &
           na_init,                                                                                &
           precipitation,                                                                          &
           sigma_gc,                                                                               &
           snow,                                                                                   &
           start_ice_microphysics

    INTERFACE bcm_actions
       MODULE PROCEDURE bcm_actions
       MODULE PROCEDURE bcm_actions_ij
    END INTERFACE bcm_actions

    INTERFACE bcm_boundary_conditions
       MODULE PROCEDURE bcm_boundary_conditions
    END INTERFACE bcm_boundary_conditions

    INTERFACE bcm_check_data_output
       MODULE PROCEDURE bcm_check_data_output
    END INTERFACE bcm_check_data_output

    INTERFACE bcm_check_data_output_pr
       MODULE PROCEDURE bcm_check_data_output_pr
    END INTERFACE bcm_check_data_output_pr

    INTERFACE bcm_check_data_output_ts
       MODULE PROCEDURE bcm_check_data_output_ts
    END INTERFACE bcm_check_data_output_ts

    INTERFACE bcm_check_parameters
       MODULE PROCEDURE bcm_check_parameters
    END INTERFACE bcm_check_parameters

    INTERFACE bcm_data_output_mask
       MODULE PROCEDURE bcm_data_output_mask
    END INTERFACE bcm_data_output_mask

    INTERFACE bcm_data_output_2d
       MODULE PROCEDURE bcm_data_output_2d
    END INTERFACE bcm_data_output_2d

    INTERFACE bcm_data_output_3d
       MODULE PROCEDURE bcm_data_output_3d
    END INTERFACE bcm_data_output_3d

    INTERFACE bcm_exchange_horiz
       MODULE PROCEDURE bcm_exchange_horiz
    END INTERFACE bcm_exchange_horiz

    INTERFACE bcm_header
       MODULE PROCEDURE bcm_header
    END INTERFACE bcm_header

    INTERFACE bcm_init
       MODULE PROCEDURE bcm_init
    END INTERFACE bcm_init

    INTERFACE bcm_init_arrays
       MODULE PROCEDURE bcm_init_arrays
    END INTERFACE bcm_init_arrays

    INTERFACE bcm_non_advective_processes
       MODULE PROCEDURE bcm_non_advective_processes
       MODULE PROCEDURE bcm_non_advective_processes_ij
    END INTERFACE bcm_non_advective_processes

    INTERFACE bcm_parin
       MODULE PROCEDURE bcm_parin
    END INTERFACE bcm_parin

    INTERFACE bcm_prognostic_equations
       MODULE PROCEDURE bcm_prognostic_equations
       MODULE PROCEDURE bcm_prognostic_equations_ij
    END INTERFACE bcm_prognostic_equations

    INTERFACE bcm_rrd_global
       MODULE PROCEDURE bcm_rrd_global_ftn
       MODULE PROCEDURE bcm_rrd_global_mpi
    END INTERFACE bcm_rrd_global

    INTERFACE bcm_rrd_local
       MODULE PROCEDURE bcm_rrd_local_ftn
       MODULE PROCEDURE bcm_rrd_local_mpi
    END INTERFACE bcm_rrd_local

    INTERFACE bcm_swap_timelevel
       MODULE PROCEDURE bcm_swap_timelevel
    END INTERFACE bcm_swap_timelevel

    INTERFACE bcm_statistics
       MODULE PROCEDURE bcm_statistics
    END INTERFACE bcm_statistics

    INTERFACE bcm_wrd_global
       MODULE PROCEDURE bcm_wrd_global
    END INTERFACE bcm_wrd_global

    INTERFACE bcm_wrd_local
       MODULE PROCEDURE bcm_wrd_local
    END INTERFACE bcm_wrd_local

    INTERFACE bcm_3d_data_averaging
       MODULE PROCEDURE bcm_3d_data_averaging
    END INTERFACE bcm_3d_data_averaging

    INTERFACE calc_liquid_water_content
       MODULE PROCEDURE calc_liquid_water_content
    END INTERFACE calc_liquid_water_content

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &bulk_cloud_parameters for the bulk cloud module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bcm_parin


    IMPLICIT NONE

    CHARACTER(LEN=100)  ::  line  !< Dummy string that contains the current line of the parameter
                                  !< file

    INTEGER(iwp)  ::  io_status   !< Status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /bulk_cloud_parameters/                                                               &
       aerosol_bulk,                                                                               &
       c_sedimentation,                                                                            &
       call_microphysics_at_all_substeps,                                                          &
       cloud_scheme,                                                                               &
       cloud_water_sedimentation,                                                                  &
       collision_turbulence,                                                                       &
       curvature_solution_effects_bulk,                                                            &
       dry_aerosol_radius,                                                                         &
       graupel,                                                                                    &
       ice_crystal_sedimentation,                                                                  &
       in_init,                                                                                    &
       limiter_sedimentation,                                                                      &
       microphysics_ice_phase,                                                                     &
       na_init,                                                                                    &
       nc_const,                                                                                   &
       precipitation_amount_interval,                                                              &
       sigma_bulk,                                                                                 &
       snow,                                                                                       &
       start_ice_microphysics,                                                                     &
       switch_off_module,                                                                          &
       ventilation_effect

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, bulk_cloud_parameters, IOSTAT=io_status )
!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    bulk_cloud_parameters namelist was found and read correctly. Set flag that
!--    bulk_cloud_model_mod is switched on.
       IF ( .NOT. switch_off_module )  bulk_cloud_model = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    bulk_cloud_parameters namelist was found, but contained errors. Print an error message
!--    containing the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'bulk_cloud_parameters', line )

    ENDIF

 END SUBROUTINE bcm_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for bulk cloud module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bcm_check_parameters


    IMPLICIT NONE


!
!-- BCM and LCM can not be used simultaneously.
    IF ( bulk_cloud_model  .AND.  cloud_droplets )  THEN
       message_string = 'bulk_cloud_model = .TRUE. is not allowed with cloud_droplets = .TRUE.'
       CALL message( 'check_parameters', 'BCM0010', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( bulk_cloud_model  .AND.  .NOT. humidity )  THEN
       WRITE( message_string, * ) 'bulk_cloud_model = ', bulk_cloud_model,                         &
                                  ' is not allowed with humidity = .FALSE.'
       CALL message( 'check_parameters', 'BCM0011', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check cloud scheme
!-- This scheme considers only saturation adjustment, i.e. water vapor surplus is converted into
!-- liquid water. No other microphysical processes are considered.
    IF ( cloud_scheme == 'saturation_adjust' )  THEN
       microphysics_sat_adjust = .TRUE.
       microphysics_seifert    = .FALSE.
       microphysics_kessler    = .FALSE.
       precipitation           = .FALSE.
       microphysics_morrison_no_rain = .FALSE.
!
!-- This scheme includes all process of the seifert beheng scheme (2001,2006). Especially rain
!-- processes are considered with prognostic quantities of qr and nr.
!-- Cloud droplet concentration is assumed to be constant and qc is diagnostic.
!-- Technical remark: The switch 'microphysics_seifert' allocates fields of qr and nr and enables
!-- all rain processes.
    ELSEIF ( cloud_scheme == 'seifert_beheng' )  THEN
       microphysics_sat_adjust = .FALSE.
       microphysics_seifert    = .TRUE.
       microphysics_kessler    = .FALSE.
       microphysics_morrison  = .FALSE.
       precipitation           = .TRUE.
       microphysics_morrison_no_rain = .FALSE.
!
!-- The kessler scheme is a simplified scheme without any prognostic quantities for microphyical
!-- variables but considering autoconversion.
    ELSEIF ( cloud_scheme == 'kessler' )  THEN
       microphysics_sat_adjust = .FALSE.
       microphysics_seifert    = .FALSE.
       microphysics_kessler    = .TRUE.
       microphysics_morrison   = .FALSE.
       precipitation           = .TRUE.
       microphysics_morrison_no_rain = .FALSE.
!
!-- The morrison scheme is an extension of the seifer beheng scheme including also relevant
!-- processes for cloud droplet size particles such as activation and an diagnostic mehtod for
!-- diffusional growth.
!-- I.e. here all processes of Seifert and Beheng as well as of the morrision scheme are used.
!-- Therefore, ztis includes prognostic quantities for qc and nc.
!-- Technical remark: The switch 'microphysics_morrison' allocates fields of qc and nc and
!-- enables diagnostic diffusional growth and activation.
    ELSEIF ( cloud_scheme == 'morrison' )  THEN
       microphysics_sat_adjust = .FALSE.
       microphysics_seifert    = .TRUE.
       microphysics_kessler    = .FALSE.
       microphysics_morrison   = .TRUE.
       precipitation           = .TRUE.
       microphysics_morrison_no_rain = .FALSE.
!
!-- The 'morrision_no_rain' scheme includes only processes of morrision scheme without the rain
!-- processes of seifert beheng. Therfore, the prog. quantities of qr and nr remain unallocated.
!-- This might be appropiate for cloud in which the size distribution is narrow, e.g. fog.
    ELSEIF ( cloud_scheme == 'morrison_no_rain' )  THEN
       microphysics_sat_adjust = .FALSE.
       microphysics_seifert    = .FALSE.
       microphysics_kessler    = .FALSE.
       microphysics_morrison   = .TRUE.
       microphysics_morrison_no_rain = .TRUE.
       precipitation           = .FALSE.
    ELSE
       message_string = 'unknown cloud microphysics scheme cloud_scheme ="' //                     &
                        TRIM( cloud_scheme ) // '"'
       CALL message( 'check_parameters', 'BCM0001', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- The use of snow and graupel are only implemented together
    IF ( snow  .AND.  .NOT. graupel  .OR.  .NOT.  snow  .AND. graupel )  THEN
       message_string = 'snow and graupel must be both switched on or off'
       CALL message( 'check_parameters', 'BCM0002', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Mixed phase microphysics can only be used with the cloud scheme 'seifert_beheng' or 'morrison'.
    IF ( microphysics_ice_phase  .AND.  cloud_scheme == 'saturation_adjust'  .OR.                  &
         microphysics_ice_phase  .AND.  cloud_scheme == 'kessler'            .OR.                  &
         microphysics_ice_phase  .AND.  cloud_scheme == 'morrison_no_rain' )                       &
    THEN
       message_string = 'cloud scheme = "' // TRIM( cloud_scheme ) // '" does not work for' //     &
                        'microphysics_ice_phase = .TRUE.'
       CALL message( 'check_parameters', 'BCM0003', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default value for the integration interval of precipitation amount
    IF ( microphysics_seifert  .OR.  microphysics_kessler )  THEN
       IF ( precipitation_amount_interval == 9999999.9_wp )  THEN
          precipitation_amount_interval = dt_do2d_xy
       ELSE
          IF ( precipitation_amount_interval > dt_do2d_xy )  THEN
             WRITE( message_string, * )  'precipitation_amount_interval = ',                       &
                precipitation_amount_interval, ' must not be larger than ',                        &
                'dt_do2d_xy = ', dt_do2d_xy
             CALL message( 'check_parameters', 'BCM0004', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

    ! TODO: find better solution for circular dependency problem
    surf_bulk_cloud_model = bulk_cloud_model
    surf_microphysics_morrison = microphysics_morrison
    surf_microphysics_seifert = microphysics_seifert
    surf_microphysics_ice_phase = microphysics_ice_phase
!
!-- Check aerosol
    IF ( aerosol_bulk == 'nacl' )  THEN
       aerosol_nacl   = .TRUE.
       aerosol_c3h4o4 = .FALSE.
       aerosol_nh4no3 = .FALSE.
    ELSEIF ( aerosol_bulk == 'c3h4o4' )  THEN
       aerosol_nacl   = .FALSE.
       aerosol_c3h4o4 = .TRUE.
       aerosol_nh4no3 = .FALSE.
    ELSEIF ( aerosol_bulk == 'nh4no3' )  THEN
       aerosol_nacl   = .FALSE.
       aerosol_c3h4o4 = .FALSE.
       aerosol_nh4no3 = .TRUE.
    ELSE
       message_string = 'unknown aerosol species "' // TRIM( aerosol_bulk ) // '" for bulk ' //    &
                        'microphysics'
       CALL message( 'check_parameters', 'BCM0005', 1, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE bcm_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for bulk cloud module
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_check_data_output( var, unit )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  unit  !<
       CHARACTER (LEN=*) ::  var   !<

       SELECT CASE ( TRIM( var ) )

          CASE ( 'nc' )
             IF ( .NOT.  microphysics_morrison )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'cloud_scheme = "morrison"'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '1/m3'

          CASE ( 'ng' )
             IF ( .NOT.  microphysics_ice_phase  .OR.  .NOT.  graupel )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'microphysics_ice_phase = ".TRUE."'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '1/m3'

          CASE ( 'ni' )
             IF ( .NOT.  microphysics_ice_phase )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'microphysics_ice_phase = ".TRUE."'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '1/m3'

          CASE ( 'nr' )
             IF ( .NOT.  microphysics_seifert )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'cloud_scheme = "seifert_beheng"'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '1/m3'

          CASE ( 'ns' )
             IF ( .NOT.  microphysics_ice_phase  .OR.  .NOT. snow )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'microphysics_ice_phase = ".TRUE."'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '1/m3'

          CASE ( 'pra*' )
             IF ( .NOT. microphysics_kessler  .AND.  .NOT. microphysics_seifert )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'cloud_scheme = "kessler" or "seifert_beheng"'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
! TODO: find solution (maybe connected to flow_statistics redesign?)
!             IF ( j == 1 )  THEN
!                message_string = 'temporal averaging of precipitation ' //                        &
!                                 'amount "' // TRIM( var ) // '" is not possible'
!                CALL message( 'check_parameters', 'BCM0113', 1, 2, 0, 6, 0 )
!             ENDIF
             unit = 'mm'

          CASE ( 'prr' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'is not available for ' //&
                                 'cloud_scheme = "saturation_adjust"'
                CALL message( 'check_parameters', 'BCM0007', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg m/s'

          CASE ( 'prr_cloud' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'is not available for ' //&
                                 'cloud_scheme = "saturation_adjust"'
                CALL message( 'check_parameters', 'BCM0007', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg m/s'

          CASE ( 'prr_graupel' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'is not available for ' //&
                                 'cloud_scheme = "saturation_adjust"'
                CALL message( 'check_parameters', 'BCM0007', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg m/s'

          CASE ( 'prr_ice' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'is not available for ' //&
                                 'cloud_scheme = "saturation_adjust"'
                CALL message( 'check_parameters', 'BCM0007', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg m/s'

          CASE ( 'prr_rain' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'is not available for ' //&
                                 'cloud_scheme = "saturation_adjust"'
                CALL message( 'check_parameters', 'BCM0007', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg m/s'

          CASE ( 'prr_snow' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'is not available for ' //&
                                 'cloud_scheme = "saturation_adjust"'
                CALL message( 'check_parameters', 'BCM0007', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg m/s'

          CASE ( 'qc' )
             unit = 'kg/kg'

          CASE ( 'qg' )
             IF ( .NOT.  microphysics_ice_phase  .OR.  .NOT.  graupel  ) THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'microphysics_ice_phase = ".TRUE."'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'qi' )
             IF ( .NOT.  microphysics_ice_phase ) THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'microphysics_ice_phase = ".TRUE."'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'qr' )
             IF ( .NOT.  microphysics_seifert ) THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'cloud_scheme = "seifert_beheng"'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE ( 'qs' )
             IF ( .NOT.  microphysics_ice_phase  .OR.  .NOT. snow ) THEN
                message_string = 'output of "' // TRIM( var ) // '" ' // 'requires ' //            &
                                 'microphysics_ice_phase = ".TRUE."'
                CALL message( 'check_parameters', 'BCM0006', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/kg'

          CASE DEFAULT
             unit = 'illegal'

       END SELECT


    END SUBROUTINE bcm_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for bulk cloud module
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_check_data_output_pr( variable, var_count, unit, dopr_unit )

       USE arrays_3d,                                                                              &
           ONLY: zu

       USE control_parameters,                                                                     &
           ONLY: data_output_pr

       USE profil_parameter,                                                                       &
           ONLY: dopr_index

       USE statistics,                                                                             &
           ONLY: hom,                                                                              &
                 statistic_regions

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  unit      !<
       CHARACTER (LEN=*) ::  variable  !<
       CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit

       INTEGER(iwp) ::  var_count      !<
       INTEGER(iwp) ::  pr_index       !<

       SELECT CASE ( TRIM( variable ) )

! TODO: make index generic: pr_index = pr_palm+1

          CASE ( 'nc' )
             IF ( .NOT.  microphysics_morrison )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not implemented for' // ' cloud_scheme /= morrison'
                CALL message( 'check_parameters', 'BCM0008', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 123
             dopr_index(var_count) = pr_index
             dopr_unit     = '1/m3'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'ni' )
             IF ( .NOT.  microphysics_ice_phase )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not implemented for' // ' microphysics_ice_phase = ".F."'
                CALL message( 'check_parameters', 'BCM0008', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 124
             dopr_index(var_count) = pr_index
             dopr_unit     = '1/m3'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'ng' )
             IF ( .NOT.  microphysics_ice_phase  .OR.  .NOT.  graupel )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not implemented for' // ' microphysics_ice_phase = ".F."' // &
                                 ' or graupel = ".F." '
                CALL message( 'check_parameters', 'BCM0008', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 126
             dopr_index(var_count) = pr_index
             dopr_unit     = '1/m3'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'ns' )
             IF ( .NOT.  microphysics_ice_phase  .OR.  .NOT.  snow )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not implemented for' // ' microphysics_ice_phase = ".F."' // &
                                 ' or snow = ".F." '
                CALL message( 'check_parameters', 'BCM0008', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 128
             dopr_index(var_count) = pr_index
             dopr_unit     = '1/m3'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'nr' )
             IF ( .NOT.  microphysics_seifert )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not implemented for' // ' cloud_scheme /= seifert_beheng'
                CALL message( 'check_parameters', 'BCM0008', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 73
             dopr_index(var_count) = pr_index
             dopr_unit     = '1/m3'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'prr' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not available for' // ' cloud_scheme = saturation_adjust'
                CALL message( 'check_parameters', 'BCM0009', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 76
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg m/s'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'prr_cloud' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not available for' // ' cloud_scheme = saturation_adjust'
                CALL message( 'check_parameters', 'BCM0009', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 131
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg m/s'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'prr_graupel' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not available for' // ' cloud_scheme = saturation_adjust'
                CALL message( 'check_parameters', 'BCM0009', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 132
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg m/s'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'prr_ice' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not available for' // ' cloud_scheme = saturation_adjust'
                CALL message( 'check_parameters', 'BCM0009', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 133
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg m/s'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'prr_rain' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not available for' // ' cloud_scheme = saturation_adjust'
                CALL message( 'check_parameters', 'BCM0009', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 134
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg m/s'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'prr_snow' )
             IF ( microphysics_sat_adjust )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not available for' // ' cloud_scheme = saturation_adjust'
                CALL message( 'check_parameters', 'BCM0009', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 135
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg m/s'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'qc' )
             pr_index = 75
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'qi' )
             IF ( .NOT.  microphysics_ice_phase )  THEN
                message_string = 'data_output_pr = ' //  TRIM( data_output_pr(var_count) ) //      &
                                 ' is not implemented for' // ' microphysics_ice_phase = ".F."'
                CALL message( 'check_parameters', 'BCM0008', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 125
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'qg' )
             IF ( .NOT.  microphysics_ice_phase  .OR.  .NOT.  graupel )  THEN
                message_string = 'data_output_pr = ' //  TRIM( data_output_pr(var_count) ) //      &
                                 ' is not implemented for' // ' microphysics_ice_phase = ".F."' // &
                                 ' or graupel = ".F." '
                CALL message( 'check_parameters', 'BCM0008', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 127
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'qs' )
             IF ( .NOT.  microphysics_ice_phase  .OR.  .NOT.  snow )  THEN
                message_string = 'data_output_pr = ' //  TRIM( data_output_pr(var_count) ) //      &
                                 ' is not implemented for' // ' microphysics_ice_phase = ".F."' // &
                                 ' or snow = ".F." '
                CALL message( 'check_parameters', 'BCM0008', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 129
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE ( 'qr' )
             IF ( .NOT.  microphysics_seifert )  THEN
                message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) //       &
                                 ' is not implemented for' // ' cloud_scheme /= seifert_beheng'
                CALL message( 'check_parameters', 'BCM0008', 1, 2, 0, 6, 0 )
             ENDIF
             pr_index = 74
             dopr_index(var_count) = pr_index
             dopr_unit     = 'kg/kg'
             unit = dopr_unit
             hom(:,2,pr_index,:)  = SPREAD( zu, 2, statistic_regions+1 )

          CASE DEFAULT
             unit = 'illegal'

       END SELECT

    END SUBROUTINE bcm_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of timeseries for bulk cloud module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bcm_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )

    INTEGER(iwp),      INTENT(IN)     ::  dots_max  !<
    INTEGER(iwp),      INTENT(INOUT)  ::  dots_num  !<

    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_label  !<
    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_unit   !<


!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( dots_num == 0  .OR.  dots_label(1)(1:1) == ' '  .OR.  dots_unit(1)(1:1) == ' ' )  CONTINUE

!
!-- Sample for user-defined time series:
!-- For each time series quantity you have to give a label and a unit, which will be used for the
!-- NetCDF file. They must not contain more than seven characters. The value of dots_num has to be
!-- increased by the number of new time series quantities. The start index for BCM time series is
!-- stored in dots_start_index_bcm and later used to address the module specific time series values.

    dots_start_index_bcm = dots_num + 1

    dots_num = dots_num + 1
    dots_label(dots_num) = 'lwp'
    dots_unit(dots_num)  = 'g/m2'

    dots_num = dots_num + 1
    dots_label(dots_num) = 'cwp'
    dots_unit(dots_num)  = 'g/m2'

    IF ( microphysics_seifert )  THEN
       dots_num = dots_num + 1
       dots_label(dots_num) = 'rwp'
       dots_unit(dots_num)  = 'g/m2'
    ENDIF

    IF ( microphysics_ice_phase )  THEN
       dots_num = dots_num + 1
       dots_label(dots_num) = 'iwp'
       dots_unit(dots_num)  = 'g/m2'

       IF ( snow  .AND.  graupel )  THEN
          dots_num = dots_num + 1
          dots_label(dots_num) = 'gwp'
          dots_unit(dots_num)  = 'g/m2'

          dots_num = dots_num + 1
          dots_label(dots_num) = 'swp'
          dots_unit(dots_num)  = 'g/m2'
       ENDIF

    ENDIF

 END SUBROUTINE bcm_check_data_output_ts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate bulk cloud module arrays and define pointers
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_init_arrays

       IMPLICIT NONE

!
!--    Liquid water content
       ALLOCATE ( ql_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!--    Frozen water content
       ALLOCATE ( qf_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!--    3D-cloud water content
       IF ( .NOT. microphysics_morrison )  THEN
          ALLOCATE( qc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF
!
!--    Precipitation amount and rate (only needed if output is switched)
       ALLOCATE( precipitation_amount(nysg:nyng,nxlg:nxrg) )

!
!--    3d-precipitation rate
       ALLOCATE( prr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( prr_cloud(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( prr_graupel(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( prr_ice(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( prr_rain(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( prr_snow(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )


       IF ( microphysics_morrison )  THEN
!
!--       3D-cloud drop water content, cloud drop concentration arrays
          ALLOCATE( nc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    nc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    nc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    qc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    qc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    qc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       IF ( microphysics_seifert )  THEN
!
!--       3D-rain water content, rain drop concentration arrays
          ALLOCATE( nr_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    nr_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    nr_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    qr_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    qr_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    qr_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       IF ( microphysics_ice_phase )  THEN
!
!--       3D-ice crystal content, ice crystal concentration arrays
          ALLOCATE( ni_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    ni_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    ni_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    qi_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    qi_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                           &
                    qi_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

          IF ( snow )  THEN
!
!--          3D-snow content, snow concentration arrays
             ALLOCATE( ns_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       ns_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       ns_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       qs_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       qs_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       qs_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

          ENDIF

          IF ( graupel )  THEN
!
!--          3D-graupel content, graupel concentration arrays
             ALLOCATE( ng_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       ng_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       ng_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       qg_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       qg_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                        &
                       qg_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

          ENDIF

       ENDIF

       IF ( ws_scheme_sca )  THEN
          IF ( microphysics_morrison )  THEN
             ALLOCATE( sums_wsqcs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             ALLOCATE( sums_wsncs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             sums_wsqcs_ws_l = 0.0_wp
             sums_wsncs_ws_l = 0.0_wp
          ENDIF
          IF ( microphysics_seifert )  THEN
             ALLOCATE( sums_wsqrs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             ALLOCATE( sums_wsnrs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             sums_wsqrs_ws_l = 0.0_wp
             sums_wsnrs_ws_l = 0.0_wp
          ENDIF
          IF ( microphysics_ice_phase )  THEN
             ALLOCATE( sums_wsqis_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             ALLOCATE( sums_wsnis_ws_l(nzb:nzt+1,0:threads_per_task-1) )
             sums_wsqis_ws_l = 0.0_wp
             sums_wsnis_ws_l = 0.0_wp
             IF ( snow )  THEN
                ALLOCATE( sums_wsqss_ws_l(nzb:nzt+1,0:threads_per_task-1) )
                ALLOCATE( sums_wsnss_ws_l(nzb:nzt+1,0:threads_per_task-1) )
                sums_wsqss_ws_l = 0.0_wp
                sums_wsnss_ws_l = 0.0_wp
             ENDIF
             IF ( graupel )  THEN
                ALLOCATE( sums_wsqgs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
                ALLOCATE( sums_wsngs_ws_l(nzb:nzt+1,0:threads_per_task-1) )
                sums_wsqgs_ws_l = 0.0_wp
                sums_wsngs_ws_l = 0.0_wp
             ENDIF
          ENDIF
       ENDIF

!
!--    Arrays needed for reasons of speed optimization for cache version.
!--    For the vector version the buffer arrays are not necessary, because the the fluxes can
!--    swapped directly inside the loops of the advection routines.
       IF ( loop_optimization /= 'vector' )  THEN
          IF ( ws_scheme_sca )  THEN
             IF ( microphysics_morrison )  THEN
                ALLOCATE( flux_s_qc(nzb+1:nzt,0:threads_per_task-1),                               &
                          diss_s_qc(nzb+1:nzt,0:threads_per_task-1),                               &
                          flux_s_nc(nzb+1:nzt,0:threads_per_task-1),                               &
                          diss_s_nc(nzb+1:nzt,0:threads_per_task-1) )
                ALLOCATE( flux_l_qc(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                       &
                          diss_l_qc(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                       &
                          flux_l_nc(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                       &
                          diss_l_nc(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
             ENDIF
             IF ( microphysics_seifert )  THEN
                ALLOCATE( flux_s_qr(nzb+1:nzt,0:threads_per_task-1),                               &
                          diss_s_qr(nzb+1:nzt,0:threads_per_task-1),                               &
                          flux_s_nr(nzb+1:nzt,0:threads_per_task-1),                               &
                          diss_s_nr(nzb+1:nzt,0:threads_per_task-1) )
                ALLOCATE( flux_l_qr(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                       &
                          diss_l_qr(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                       &
                          flux_l_nr(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                       &
                          diss_l_nr(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
             ENDIF
             IF ( microphysics_ice_phase )  THEN
                ALLOCATE( flux_s_qi(nzb+1:nzt,0:threads_per_task-1),                               &
                          diss_s_qi(nzb+1:nzt,0:threads_per_task-1),                               &
                          flux_s_ni(nzb+1:nzt,0:threads_per_task-1),                               &
                          diss_s_ni(nzb+1:nzt,0:threads_per_task-1) )
                ALLOCATE( flux_l_qi(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                       &
                          diss_l_qi(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                       &
                          flux_l_ni(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                       &
                          diss_l_ni(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
                IF ( snow )  THEN
                   ALLOCATE( flux_s_qs(nzb+1:nzt,0:threads_per_task-1),                            &
                             diss_s_qs(nzb+1:nzt,0:threads_per_task-1),                            &
                             flux_s_ns(nzb+1:nzt,0:threads_per_task-1),                            &
                             diss_s_ns(nzb+1:nzt,0:threads_per_task-1) )
                   ALLOCATE( flux_l_qs(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                    &
                             diss_l_qs(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                    &
                             flux_l_ns(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                    &
                             diss_l_ns(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
                ENDIF
                IF ( graupel ) THEN
                   ALLOCATE( flux_s_qg(nzb+1:nzt,0:threads_per_task-1),                            &
                             diss_s_qg(nzb+1:nzt,0:threads_per_task-1),                            &
                             flux_s_ng(nzb+1:nzt,0:threads_per_task-1),                            &
                             diss_s_ng(nzb+1:nzt,0:threads_per_task-1) )
                   ALLOCATE( flux_l_qg(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                    &
                             diss_l_qg(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                    &
                             flux_l_ng(nzb+1:nzt,nys:nyn,0:threads_per_task-1),                    &
                             diss_l_ng(nzb+1:nzt,nys:nyn,0:threads_per_task-1) )
                ENDIF
             ENDIF
          ENDIF
       ENDIF

!
!--   Allocate arrays needed for riming process
      ALLOCATE( dep_rate_ice(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                       &
                dep_rate_snow(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                      &
                rime_ice_cloud(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                     &
                rime_ice_rain(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                      &
                rime_snow_cloud(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                    &
                rime_snow_rain(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!--    Initial assignment of the pointers
       ql => ql_1
       qf => qf_1
       IF ( .NOT. microphysics_morrison )  THEN
          qc => qc_1
       ENDIF
       IF ( microphysics_morrison )  THEN
          qc => qc_1;  qc_p  => qc_2;  tqc_m  => qc_3
          nc => nc_1;  nc_p  => nc_2;  tnc_m  => nc_3
       ENDIF
       IF ( microphysics_seifert )  THEN
          qr => qr_1;  qr_p  => qr_2;  tqr_m  => qr_3
          nr => nr_1;  nr_p  => nr_2;  tnr_m  => nr_3
       ENDIF
       IF ( microphysics_ice_phase )  THEN
          qi => qi_1;  qi_p  => qi_2;  tqi_m  => qi_3
          ni => ni_1;  ni_p  => ni_2;  tni_m  => ni_3
          IF ( snow )  THEN
             qs => qs_1;  qs_p  => qs_2;  tqs_m  => qs_3
             ns => ns_1;  ns_p  => ns_2;  tns_m  => ns_3
          ENDIF
          IF ( graupel )  THEN
             qg => qg_1;  qg_p  => qg_2;  tqg_m  => qg_3
             ng => ng_1;  ng_p  => ng_2;  tng_m  => ng_3
          ENDIF
       ENDIF

    END SUBROUTINE bcm_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the bulk cloud module
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_init

       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<

       IF ( debug_output )  CALL debug_message( 'bcm_init', 'start' )

       IF ( bulk_cloud_model )  THEN
          IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
!
!--          Initialize the remaining quantities
             IF ( microphysics_morrison )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      qc(:,j,i) = 0.0_wp
                      nc(:,j,i) = 0.0_wp
                   ENDDO
                ENDDO
             ENDIF

             IF ( microphysics_seifert )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      qr(:,j,i) = 0.0_wp
                      nr(:,j,i) = 0.0_wp
                   ENDDO
                ENDDO
             ENDIF
!
!--          Initialize the remaining quantities
             IF ( microphysics_ice_phase )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      qi(:,j,i) = 0.0_wp
                      ni(:,j,i) = 0.0_wp
                   ENDDO
                ENDDO
                IF ( snow )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         qs(:,j,i) = 0.0_wp
                         ns(:,j,i) = 0.0_wp
                      ENDDO
                   ENDDO
                ENDIF
                IF ( graupel ) THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         qg(:,j,i) = 0.0_wp
                         ng(:,j,i) = 0.0_wp
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF
!
!--          Liquid water content and precipitation amount are zero at beginning of the simulation.
             ql = 0.0_wp
             qf = 0.0_wp
             qc = 0.0_wp
             precipitation_amount = 0.0_wp
             prr = 0.0_wp
             prr_cloud = 0.0_wp
             prr_graupel = 0.0_wp
             prr_ice = 0.0_wp
             prr_rain = 0.0_wp
             prr_snow = 0.0_wp
!
!--          Initialize old and new time levels.
             IF ( microphysics_morrison )  THEN
                tqc_m = 0.0_wp
                tnc_m = 0.0_wp
                qc_p  = qc
                nc_p  = nc
             ENDIF
             IF ( microphysics_seifert )  THEN
                tqr_m = 0.0_wp
                tnr_m = 0.0_wp
                qr_p  = qr
                nr_p  = nr
             ENDIF
             IF ( microphysics_ice_phase )  THEN
                tqi_m = 0.0_wp
                tni_m = 0.0_wp
                qi_p  = qi
                ni_p  = ni
                IF ( snow )  THEN
                   tqs_m = 0.0_wp
                   tns_m = 0.0_wp
                   qs_p  = qs
                   ns_p  = ns
                ENDIF
                IF ( graupel ) THEN
                   tqg_m = 0.0_wp
                   tng_m = 0.0_wp
                   qg_p  = qg
                   ng_p  = ng
                ENDIF
             ENDIF
          ENDIF ! Only if not read_restart_data
!
!--       Constant for the sedimentation of cloud water (2-moment cloud physics)
          sed_qc_const = k_st * ( 3.0_wp / ( 4.0_wp * pi * rho_l ) )**( 2.0_wp / 3.0_wp ) *        &
                         EXP( 5.0_wp * LOG( sigma_gc )**2 )

!
!--       Calculate timestep according to precipitation
          IF ( microphysics_seifert )  THEN
             dt_precipitation = c_sedimentation * MINVAL( dzu(nzb+2:nzt) ) / w_precipitation
          ENDIF

!
!--       Set constants for certain aerosol type
          IF ( microphysics_morrison )  THEN
             IF ( aerosol_nacl ) THEN
                molecular_weight_of_solute = 0.05844_wp
                rho_s                      = 2165.0_wp
                vanthoff                   = 2.0_wp
             ELSEIF ( aerosol_c3h4o4 ) THEN
                molecular_weight_of_solute = 0.10406_wp
                rho_s                      = 1600.0_wp
                vanthoff                   = 1.37_wp
             ELSEIF ( aerosol_nh4no3 ) THEN
                molecular_weight_of_solute = 0.08004_wp
                rho_s                      = 1720.0_wp
                vanthoff                   = 2.31_wp
             ENDIF
          ENDIF
!
!--       Pre-calculate frequently calculated fractions of pi and rho_l
          pirho_l  = pi * rho_l / 6.0_wp
          dpirho_l = 1.0_wp / pirho_l
!
!--       Calculate collision integrals: First, define the coefficients of the species (Here values of
!--       SB2006, Table 1 are used). Second, calculate the collision integrals as described in SB2006
!--       (Eq. 90-93, see their Appendix C).
          CALL define_cloud_coefficients
          CALL init_collision_integrals
          CALL init_coefficients

          IF ( debug_output )  CALL debug_message( 'bcm_init', 'end' )

       ELSE

          IF ( debug_output )  CALL debug_message( 'bcm_init skipped', 'end' )

       ENDIF

    END SUBROUTINE bcm_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for bulk cloud module
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_header ( io )


       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) ::  io  !< Unit of the output file

!
!--    Write bulk cloud module header
       WRITE ( io, 1 )

       WRITE ( io, 2 )
       WRITE ( io, 3 )

       IF ( microphysics_kessler )  THEN
          WRITE ( io, 4 ) 'Kessler-Scheme'
       ENDIF

       IF ( microphysics_seifert )  THEN
          WRITE ( io, 4 ) 'Seifert-Beheng-Scheme'
          IF ( cloud_water_sedimentation )  WRITE ( io, 5 )
          IF ( collision_turbulence )  WRITE ( io, 6 )
          IF ( ventilation_effect )  WRITE ( io, 7 )
          IF ( limiter_sedimentation )  WRITE ( io, 8 )
       ENDIF

       WRITE ( io, 20 )
       WRITE ( io, 21 ) surface_pressure
       WRITE ( io, 22 ) r_d
       WRITE ( io, 23 ) rho_surface
       WRITE ( io, 24 ) c_p
       WRITE ( io, 25 ) l_v

       IF ( microphysics_seifert )  THEN
          WRITE ( io, 26 ) 1.0E-6_wp * nc_const
          WRITE ( io, 27 ) c_sedimentation
       ENDIF


1   FORMAT ( //' Bulk cloud module information:'/ ' ------------------------------------------'/ )
2   FORMAT ( '--> Bulk scheme with liquid water potential temperature and'/                        &
             '    total water content is used.' )
3   FORMAT ( '--> Condensation is parameterized via 0% - or 100% scheme.' )
4   FORMAT ( '--> Precipitation parameterization via ', A )

5   FORMAT ( '--> Cloud water sedimentation parameterization via Stokes law' )
6   FORMAT ( '--> Turbulence effects on precipitation process' )
7   FORMAT ( '--> Ventilation effects on evaporation of rain drops' )
8   FORMAT ( '--> Slope limiter used for sedimentation process' )

20  FORMAT ( '--> Essential parameters:' )
21  FORMAT ( '       Surface pressure             :   p_0   = ', F7.2, ' hPa')
22  FORMAT ( '       Gas constant                 :   R     = ', F5.1, ' J/(kg K)')
23  FORMAT ( '       Density of air               :   rho_0 = ', F6.3, ' kg/m**3')
24  FORMAT ( '       Specific heat cap.           :   c_p   = ', F6.1, ' J/(kg K)')
25  FORMAT ( '       Vapourization heat           :   L_v   = ', E9.2, ' J/kg')
26  FORMAT ( '       Droplet density              :   N_c   = ', F6.1, ' 1/cm**3' )
27  FORMAT ( '       Sedimentation Courant number :   C_s   = ', F4.1 )


    END SUBROUTINE bcm_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_actions( location )


    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( ws_scheme_sca )  THEN

             IF ( microphysics_morrison )  THEN
                sums_wsqcs_ws_l = 0.0_wp
                sums_wsncs_ws_l = 0.0_wp
             ENDIF
             IF ( microphysics_seifert )  THEN
                sums_wsqrs_ws_l = 0.0_wp
                sums_wsnrs_ws_l = 0.0_wp
             ENDIF
             IF ( microphysics_ice_phase )  THEN
                sums_wsqis_ws_l = 0.0_wp
                sums_wsnis_ws_l = 0.0_wp
                IF ( snow )  THEN
                   sums_wsqss_ws_l = 0.0_wp
                   sums_wsnss_ws_l = 0.0_wp
                ENDIF
                IF ( graupel )  THEN
                   sums_wsqgs_ws_l = 0.0_wp
                   sums_wsngs_ws_l = 0.0_wp
                ENDIF
             ENDIF

          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT

    END SUBROUTINE bcm_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for grid points i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_actions_ij( i, j, location )


    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string

    INTEGER(iwp)  ::  dummy                !< call location string

    INTEGER(iwp), INTENT(IN) ::  i         !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  j         !< grid index in y-direction


    IF ( bulk_cloud_model )   dummy = i + j

    SELECT CASE ( location )

       CASE ( 'before_timestep' )

          IF ( ws_scheme_sca )  THEN
             IF ( microphysics_morrison )  THEN
                sums_wsqcs_ws_l = 0.0_wp
                sums_wsncs_ws_l = 0.0_wp
             ENDIF
             IF ( microphysics_seifert )  THEN
                sums_wsqrs_ws_l = 0.0_wp
                sums_wsnrs_ws_l = 0.0_wp
             ENDIF
             IF ( microphysics_ice_phase )  THEN
                sums_wsqis_ws_l = 0.0_wp
                sums_wsnis_ws_l = 0.0_wp
                IF ( snow )  THEN
                   sums_wsqss_ws_l = 0.0_wp
                   sums_wsnss_ws_l = 0.0_wp
                ENDIF
                IF ( graupel )  THEN
                   sums_wsqgs_ws_l = 0.0_wp
                   sums_wsngs_ws_l = 0.0_wp
                ENDIF
             ENDIF
          ENDIF

       CASE DEFAULT
          CONTINUE

    END SELECT


    END SUBROUTINE bcm_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_non_advective_processes


       CALL cpu_log( log_point(51), 'microphysics', 'start' )

       IF ( .NOT. microphysics_sat_adjust  .AND.  ( intermediate_timestep_count == 1  .OR.         &
            call_microphysics_at_all_substeps ) )                                                  &
       THEN

          IF ( large_scale_forcing  .AND.  lsf_surf ) THEN
!
!--          Calculate vertical profile of the hydrostatic pressure (hyp)
             hyp    = barometric_formula(zu, pt_surface *                                          &
                      exner_function(surface_pressure * 100.0_wp), surface_pressure * 100.0_wp)
             d_exner = exner_function_invers(hyp)
             exner = 1.0_wp / exner_function_invers(hyp)
             hyrho  = ideal_gas_law_rho_pt(hyp, pt_init)
!
!--          Compute reference density
             rho_surface = ideal_gas_law_rho(surface_pressure * 100.0_wp,                          &
                           pt_surface * exner_function(surface_pressure * 100.0_wp))
          ENDIF

!
!--       Compute length of time step
          IF ( call_microphysics_at_all_substeps )  THEN
             dt_micro = dt_3d * weight_pres(intermediate_timestep_count)
          ELSE
             dt_micro = dt_3d
          ENDIF
!
!--       Reset precipitation rate
          IF ( intermediate_timestep_count == 1 )  THEN
             prr = 0.0_wp
             prr_cloud = 0.0_wp
             prr_graupel = 0.0_wp
             prr_ice = 0.0_wp
             prr_rain = 0.0_wp
             prr_snow = 0.0_wp
          ENDIF
!
!--       Compute cloud physics.
!--       Here the the simple kessler scheme is used.
          IF ( microphysics_kessler )  THEN
             CALL autoconversion_kessler
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud
!
!--       Here the seifert beheng scheme is used. Cloud concentration is assumed to a constant value
!--       an qc a diagnostic value.
          ELSEIF ( microphysics_seifert  .AND.  .NOT. microphysics_morrison )  THEN
             CALL adjust_cloud
             IF ( microphysics_ice_phase  .AND.  simulated_time > start_ice_microphysics )  THEN
                 CALL adjust_ice
                 CALL nucleation_ice
                 CALL homogeneous_freezing_cloud
                 CALL deposition_ice
                 IF ( snow  .AND.  graupel ) THEN
                    CALL deposition_snow
                    CALL deposition_graupel
                    CALL selfcollection_ice
                    CALL selfcollection_snow
                    CALL selfcollection_graupel
                    CALL collection_graupel_ice
                    CALL collection_snow_ice
                    CALL collection_graupel_snow
                    CALL riming_graupel_cloud
                    CALL riming_graupel_rain
                    CALL riming_ice_cloud
                    CALL riming_ice_rain
                    CALL riming_snow_cloud
                    CALL riming_snow_rain
                    CALL heterogeneous_freezing_rain
                    CALL melting_ice
                    CALL melting_snow
                    CALL melting_graupel
                    CALL evaporation_graupel
                    CALL evaporation_snow
                 ENDIF
             ENDIF
             CALL autoconversion
             CALL accretion
             CALL selfcollection_breakup_rain
             CALL evaporation_rain
             CALL sedimentation_rain
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud
             IF ( microphysics_ice_phase  .AND.  simulated_time > start_ice_microphysics )  THEN
                 IF ( ice_crystal_sedimentation )             CALL sedimentation_ice
                 IF ( snow_sedimentation  .AND.  snow )       CALL sedimentation_snow
                 IF ( graupel_sedimentation  .AND.  graupel ) CALL sedimentation_graupel
             ENDIF
!
!--       Here the morrison scheme is used. No rain processes are considered and qr and nr are not
!--       allocated.
          ELSEIF ( microphysics_morrison_no_rain  .AND.  .NOT. microphysics_seifert )  THEN
             CALL activation_cloud
             CALL condensation_cloud
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud
!
!--       Here the full morrison scheme is used and all processes of Seifert and Beheng are included
          ELSEIF ( microphysics_morrison  .AND.  microphysics_seifert )  THEN
             CALL adjust_cloud
             CALL activation_cloud
             CALL condensation_cloud
             IF ( microphysics_ice_phase  .AND.  simulated_time > start_ice_microphysics )  THEN
                 CALL adjust_ice
                 CALL nucleation_ice
                 CALL homogeneous_freezing_cloud
                 CALL deposition_ice
                 IF ( snow  .AND.  graupel ) THEN
                    CALL deposition_snow
                    CALL deposition_graupel
                    CALL selfcollection_ice
                    CALL selfcollection_snow
                    CALL selfcollection_graupel
                    CALL collection_graupel_ice
                    CALL collection_snow_ice
                    CALL collection_graupel_snow
                    CALL riming_graupel_cloud
                    CALL riming_graupel_rain
                    CALL riming_ice_cloud
                    CALL riming_ice_rain
                    CALL riming_snow_cloud
                    CALL riming_snow_rain
                    CALL heterogeneous_freezing_rain
                    CALL melting_ice
                    CALL melting_snow
                    CALL melting_graupel
                    CALL evaporation_graupel
                    CALL evaporation_snow
                 ENDIF
             ENDIF
             CALL autoconversion
             CALL accretion
             CALL selfcollection_breakup_rain
             CALL evaporation_rain
             CALL sedimentation_rain
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud
             IF ( microphysics_ice_phase  .AND.  simulated_time > start_ice_microphysics )  THEN
                 IF ( ice_crystal_sedimentation )             CALL sedimentation_ice
                 IF ( snow_sedimentation  .AND.  snow )       CALL sedimentation_snow
                 IF ( graupel_sedimentation  .AND.  graupel ) CALL sedimentation_graupel
             ENDIF

          ENDIF

          CALL calc_precipitation_amount

       ENDIF

       CALL cpu_log( log_point(51), 'microphysics', 'stop' )

    END SUBROUTINE bcm_non_advective_processes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for grid points i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_non_advective_processes_ij( i, j )


       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<

       IF ( .NOT. microphysics_sat_adjust  .AND.  ( intermediate_timestep_count == 1  .OR.         &
            call_microphysics_at_all_substeps ) )                                                  &
       THEN

          IF ( large_scale_forcing  .AND.  lsf_surf ) THEN
!
!--          Calculate vertical profile of the hydrostatic pressure (hyp)
             hyp    = barometric_formula(zu, pt_surface *                                          &
                      exner_function(surface_pressure * 100.0_wp), surface_pressure * 100.0_wp)
             d_exner = exner_function_invers(hyp)
             exner = 1.0_wp / exner_function_invers(hyp)
             hyrho  = ideal_gas_law_rho_pt(hyp, pt_init)
!
!--          Compute reference density
             rho_surface = ideal_gas_law_rho(surface_pressure * 100.0_wp,                          &
                           pt_surface * exner_function(surface_pressure * 100.0_wp))
          ENDIF

!
!--       Compute length of time step
          IF ( call_microphysics_at_all_substeps )  THEN
             dt_micro = dt_3d * weight_pres(intermediate_timestep_count)
          ELSE
             dt_micro = dt_3d
          ENDIF
!
!--       Reset precipitation rate
          IF ( intermediate_timestep_count == 1 )  THEN
             prr(:,j,i) = 0.0_wp
             prr_cloud(:,j,i) = 0.0_wp
             prr_graupel(:,j,i) = 0.0_wp
             prr_ice(:,j,i) = 0.0_wp
             prr_rain(:,j,i) = 0.0_wp
             prr_snow(:,j,i) = 0.0_wp
          ENDIF
!
!--       Compute cloud physics
!--       Here the the simple kessler scheme is used.
          IF( microphysics_kessler )  THEN
             CALL autoconversion_kessler_ij( i, j )
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud_ij( i, j )
!
!--       Here the seifert beheng scheme is used. Cloud concentration is assumed to a constant value
!--       an qc a diagnostic value.
          ELSEIF ( microphysics_seifert  .AND.  .NOT. microphysics_morrison )  THEN
             CALL adjust_cloud_ij( i, j )
             IF ( microphysics_ice_phase  .AND.  simulated_time > start_ice_microphysics )  THEN
                 CALL adjust_ice_ij( i, j )
                 CALL nucleation_ice_ij( i, j )
                 CALL homogeneous_freezing_cloud_ij( i, j )
                 CALL deposition_ice_ij( i, j )
                 IF ( snow  .AND.  graupel )  THEN
                    CALL deposition_snow_ij( i, j )
                    CALL deposition_graupel_ij( i, j )
                    CALL selfcollection_ice_ij( i, j )
                    CALL selfcollection_snow_ij( i, j )
                    CALL selfcollection_graupel_ij( i, j )
                    CALL collection_graupel_ice_ij( i, j )
                    CALL collection_snow_ice_ij( i, j )
                    CALL collection_graupel_snow_ij( i, j )
                    CALL riming_graupel_cloud_ij( i, j )
                    CALL riming_graupel_rain_ij( i, j )
                    CALL riming_ice_cloud_ij( i, j )
                    CALL riming_ice_rain_ij( i, j )
                    CALL riming_snow_cloud_ij( i, j )
                    CALL riming_snow_rain_ij( i, j )
                    CALL heterogeneous_freezing_rain_ij( i, j )
                    CALL melting_ice_ij( i, j )
                    CALL melting_snow_ij( i, j )
                    CALL melting_graupel_ij( i, j )
                    CALL evaporation_graupel_ij( i, j )
                    CALL evaporation_snow_ij( i, j )
                 ENDIF
             ENDIF
             CALL autoconversion_ij( i, j )
             CALL accretion_ij( i, j )
             CALL selfcollection_breakup_rain_ij( i, j )
             CALL evaporation_rain_ij( i, j )
             CALL sedimentation_rain_ij( i, j )
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud_ij( i ,j )
             IF ( microphysics_ice_phase  .AND.  simulated_time > start_ice_microphysics )  THEN
                 IF ( ice_crystal_sedimentation )             CALL sedimentation_ice_ij( i, j )
                 IF ( snow_sedimentation  .AND.  snow )       CALL sedimentation_snow_ij( i, j )
                 IF ( graupel_sedimentation .AND.  graupel )  CALL sedimentation_graupel_ij( i, j )
             ENDIF
!
!--       Here the morrison scheme is used. No rain processes are considered and qr and nr are not
!--       allocated.
          ELSEIF ( microphysics_morrison_no_rain  .AND.  .NOT. microphysics_seifert )  THEN
             CALL activation_cloud_ij( i, j )
             CALL condensation_cloud_ij( i, j )
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud_ij( i, j )
!
!--       Here the full morrison scheme is used and all processes of Seifert and Beheng are included
          ELSEIF ( microphysics_morrison  .AND.  microphysics_seifert )  THEN
             CALL adjust_cloud_ij( i, j )
             CALL activation_cloud_ij( i, j )
             CALL condensation_cloud_ij( i, j )
             IF ( microphysics_ice_phase  .AND.  simulated_time > start_ice_microphysics )  THEN
                 CALL adjust_ice_ij( i, j )
                 CALL nucleation_ice_ij( i, j )
                 CALL homogeneous_freezing_cloud_ij( i, j )
                 CALL deposition_ice_ij( i, j )
                 IF ( snow  .AND.  graupel )  THEN
                    CALL deposition_snow_ij( i, j )
                    CALL deposition_graupel_ij( i, j )
                    CALL selfcollection_ice_ij( i, j )
                    CALL selfcollection_snow_ij( i, j )
                    CALL selfcollection_graupel_ij( i, j )
                    CALL collection_graupel_ice_ij( i, j )
                    CALL collection_snow_ice_ij( i, j )
                    CALL collection_graupel_snow_ij( i, j )
                    CALL riming_graupel_cloud_ij( i, j )
                    CALL riming_graupel_rain_ij( i, j )
                    CALL riming_ice_cloud_ij( i, j )
                    CALL riming_ice_rain_ij( i, j )
                    CALL riming_snow_cloud_ij( i, j )
                    CALL riming_snow_rain_ij( i, j )
                    CALL heterogeneous_freezing_rain_ij( i, j )
                    CALL melting_ice_ij( i, j )
                    CALL melting_snow_ij( i, j )
                    CALL melting_graupel_ij( i, j )
                    CALL evaporation_graupel_ij( i, j )
                    CALL evaporation_snow_ij( i, j )
                 ENDIF
             ENDIF
             CALL autoconversion_ij( i, j )
             CALL accretion_ij( i, j )
             CALL selfcollection_breakup_rain_ij( i, j )
             CALL evaporation_rain_ij( i, j )
             CALL sedimentation_rain_ij( i, j )
             IF ( cloud_water_sedimentation )  CALL sedimentation_cloud_ij( i, j )
             IF ( microphysics_ice_phase  .AND.  simulated_time > start_ice_microphysics )  THEN
                 IF ( ice_crystal_sedimentation )              CALL sedimentation_ice_ij( i, j )
                 IF ( snow_sedimentation  .AND.  snow )        CALL sedimentation_snow_ij( i, j )
                 IF ( graupel_sedimentation  .AND.  graupel )  CALL sedimentation_graupel_ij( i, j )
             ENDIF

          ENDIF

          CALL calc_precipitation_amount_ij( i, j )

       ENDIF

    END SUBROUTINE bcm_non_advective_processes_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange ghostpoints
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_exchange_horiz( location )

       USE exchange_horiz_mod,                                                                     &
           ONLY:  exchange_horiz

       CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

       SELECT CASE ( location )

          CASE ( 'before_prognostic_equation' )

             IF ( .NOT. microphysics_sat_adjust  .AND.  ( intermediate_timestep_count == 1  .OR.   &
                  call_microphysics_at_all_substeps ) )                                            &
             THEN
                IF ( microphysics_morrison )  THEN
                   CALL exchange_horiz( nc, nbgp )
                   CALL exchange_horiz( qc, nbgp )
                ENDIF
                IF ( microphysics_seifert ) THEN
                   CALL exchange_horiz( qr, nbgp )
                   CALL exchange_horiz( nr, nbgp )
                ENDIF
                IF ( microphysics_ice_phase ) THEN
                   CALL exchange_horiz( qi, nbgp )
                   CALL exchange_horiz( ni, nbgp )
                   IF ( snow )  THEN
                      CALL exchange_horiz( qs, nbgp )
                      CALL exchange_horiz( ns, nbgp )
                   ENDIF
                   IF ( graupel )  THEN
                      CALL exchange_horiz( qg, nbgp )
                      CALL exchange_horiz( ng, nbgp )
                   ENDIF
                ENDIF
                CALL exchange_horiz( q, nbgp )
                CALL exchange_horiz( pt, nbgp )
             ENDIF

          CASE ( 'after_prognostic_equation' )

             IF ( microphysics_morrison )  THEN
                CALL exchange_horiz( qc_p, nbgp )
                CALL exchange_horiz( nc_p, nbgp )
             ENDIF
             IF ( microphysics_seifert )  THEN
                CALL exchange_horiz( qr_p, nbgp )
                CALL exchange_horiz( nr_p, nbgp )
             ENDIF
             IF ( microphysics_ice_phase )  THEN
                CALL exchange_horiz( qi_p, nbgp )
                CALL exchange_horiz( ni_p, nbgp )
                IF ( snow )  THEN
                   CALL exchange_horiz( qs_p, nbgp )
                   CALL exchange_horiz( ns_p, nbgp )
                ENDIF
                IF ( graupel )  THEN
                   CALL exchange_horiz( qg_p, nbgp )
                   CALL exchange_horiz( ng_p, nbgp )
                ENDIF
             ENDIF

          CASE ( 'after_anterpolation' )

             IF ( microphysics_morrison )  THEN
                CALL exchange_horiz( qc, nbgp )
                CALL exchange_horiz( nc, nbgp )
             ENDIF
             IF ( microphysics_seifert )  THEN
                CALL exchange_horiz( qr, nbgp )
                CALL exchange_horiz( nr, nbgp )
             ENDIF
             IF ( microphysics_ice_phase )  THEN
                CALL exchange_horiz( qi, nbgp )
                CALL exchange_horiz( ni, nbgp )
                IF ( snow )  THEN
                   CALL exchange_horiz( qs_p, nbgp )
                   CALL exchange_horiz( ns_p, nbgp )
                ENDIF
                IF ( graupel )  THEN
                   CALL exchange_horiz( qg_p, nbgp )
                   CALL exchange_horiz( ng_p, nbgp )
                ENDIF
             ENDIF

       END SELECT

    END SUBROUTINE bcm_exchange_horiz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_prognostic_equations


       INTEGER(iwp) ::  i         !< grid index in x-direction
       INTEGER(iwp) ::  j         !< grid index in y-direction
       INTEGER(iwp) ::  k         !< grid index in z-direction

       REAL(wp)     ::  sbt  !<

       CALL cpu_log( log_point(67), 'all_bmc-prog_equation', 'start' )
!
!--    If required, calculate prognostic equations for cloud water content and cloud drop
!--    concentration.
       IF ( microphysics_morrison )  THEN
!
!--       Calculate prognostic equation for cloud water content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN
             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qc, 'qc' )

          ENDIF
!
!--       qc-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, qc, 'qc',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, qc, 'qc',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( qc )
                ENDIF
             ELSE
                CALL advec_s_up( qc )
             ENDIF
          ENDIF

          CALL diffusion_s( qc, surf_top%qcsws, surf_def%qcsws, surf_lsm%qcsws, surf_usm%qcsws )

!
!--       Prognostic equation for cloud water content
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   qc_p(k,j,i) = qc(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tqc_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * qc(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( qc_p(k,j,i) < 0.0_wp )  qc_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqc_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqc_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqc_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
!
!--       Calculate prognostic equation for cloud drop concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( nc, 'nc' )

          ENDIF

!
!--       nc-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, nc, 'nc',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, nc, 'nc',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( nc )
                ENDIF
             ELSE
                CALL advec_s_up( nc )
             ENDIF
          ENDIF

          CALL diffusion_s( nc, surf_top%ncsws, surf_def%ncsws, surf_lsm%ncsws, surf_usm%ncsws )
!
!--       Prognostic equation for cloud drop concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   nc_p(k,j,i) = nc(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tnc_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * nc(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( nc_p(k,j,i) < 0.0_wp )  nc_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnc_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnc_m(k,j,i) =  -9.5625_wp * tend(k,j,i) + 5.3125_wp * tnc_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
       ENDIF

!
!--    If required, calculate prognostic equations for ice crystal content and ice crystal
!--    concentration
       IF ( microphysics_ice_phase )  THEN
!
!--       Calculate prognostic equation for ice crystal content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qi, 'qi' )

          ENDIF

!
!--       qi-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, qi, 'qi',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, qi, 'qi',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( qi )
                ENDIF
             ELSE
                CALL advec_s_up( qi )
             ENDIF
          ENDIF

          CALL diffusion_s( qi, surf_top%qisws, surf_def%qisws, surf_lsm%qisws, surf_usm%qisws )
!
!--       Prognostic equation for ice crystal mixing ratio
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   qi_p(k,j,i) = qi(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tqi_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * qi(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( qi_p(k,j,i) < 0.0_wp )  qi_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqi_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqi_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqi_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
!
!--       Calculate prognostic equation for ice crystal concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( ni, 'ni' )

          ENDIF

!
!--       ni-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, ni, 'ni',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, ni, 'ni',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( ni )
                ENDIF
             ELSE
                CALL advec_s_up( ni )
             ENDIF
          ENDIF

          CALL diffusion_s( ni, surf_top%nisws, surf_def%nisws, surf_lsm%nisws, surf_usm%nisws )
!
!--       Prognostic equation for ice crystal concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   ni_p(k,j,i) = ni(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tni_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * ni(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( ni_p(k,j,i) < 0.0_wp )  ni_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tni_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tni_m(k,j,i) =  -9.5625_wp * tend(k,j,i) + 5.3125_wp * tni_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF

       ENDIF

!
!--    If required, calculate prognostic equations for snow content and snow
!--    concentration
       IF ( microphysics_ice_phase  .AND.  snow )  THEN
!
!--       Calculate prognostic equation for ice crystal content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qs, 'qs' )

          ENDIF

!
!--       qs-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, qs, 'qs',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, qs, 'qs',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( qs )
                ENDIF
             ELSE
                CALL advec_s_up( qs )
             ENDIF
          ENDIF

          CALL diffusion_s( qs, surf_top%qisws, surf_def%qisws, surf_lsm%qisws, surf_usm%qisws )
!
!--       Prognostic equation for snow mixing ratio
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   qs_p(k,j,i) = qs(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tqs_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * qs(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( qs_p(k,j,i) < 0.0_wp )  qs_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqs_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqs_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqs_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
!
!--       Calculate prognostic equation for ice crystal concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( ns, 'ns' )

          ENDIF

!
!--       ni-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, ns, 'ns',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, ns, 'ns',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( ns )
                ENDIF
             ELSE
                CALL advec_s_up( ns )
             ENDIF
          ENDIF

          CALL diffusion_s( ns, surf_top%nisws, surf_def%nisws, surf_lsm%nisws, surf_usm%nisws )
!
!--       Prognostic equation for snow number concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   ns_p(k,j,i) = ns(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tns_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * ns(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( ns_p(k,j,i) < 0.0_wp )  ns_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tns_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tns_m(k,j,i) =  -9.5625_wp * tend(k,j,i) + 5.3125_wp * tns_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
       ENDIF

!
!--    If required, calculate prognostic equations for ice crystal content and ice crystal
!--    concentration
       IF ( microphysics_ice_phase   .AND.  graupel )  THEN
!
!--       Calculate prognostic equation for graupel content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qg, 'qg' )

          ENDIF

!
!--       qg-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, qg, 'qg',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, qg, 'qg',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( qg )
                ENDIF
             ELSE
                CALL advec_s_up( qg )
             ENDIF
          ENDIF

          CALL diffusion_s( qg, surf_top%qisws, surf_def%qisws, surf_lsm%qisws, surf_usm%qisws )
!
!--       Prognostic equation for graupel mixing ratio
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   qg_p(k,j,i) = qg(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tqg_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * qg(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( qg_p(k,j,i) < 0.0_wp )  qg_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqg_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqg_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqg_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
!
!--       Calculate prognostic equation for graupel number concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( ng, 'ng ' )

          ENDIF

!
!--       ng-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, ng, 'ng',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, ng, 'ng',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( ng )
                ENDIF
             ELSE
                CALL advec_s_up( ng )
             ENDIF
          ENDIF

          CALL diffusion_s( ng, surf_top%nisws, surf_def%nisws, surf_lsm%nisws, surf_usm%nisws )
!
!--       Prognostic equation for ice crystal concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   ng_p(k,j,i) = ng(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tng_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * ng(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( ng_p(k,j,i) < 0.0_wp )  ng_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tng_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tng_m(k,j,i) =  -9.5625_wp * tend(k,j,i) + 5.3125_wp * tng_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
       ENDIF


!
!--    If required, calculate prognostic equations for rain water content and rain drop
!--    concentration.
       IF ( microphysics_seifert )  THEN
!
!--       Calculate prognostic equation for rain water content
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( qr, 'qr' )

          ENDIF

!
!--       qr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, qr, 'qr',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, qr, 'qr',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( qr )
                ENDIF
             ELSE
                CALL advec_s_up( qr )
             ENDIF
          ENDIF

          CALL diffusion_s( qr, surf_top%qrsws, surf_def%qrsws, surf_lsm%qrsws, surf_usm%qrsws )
!
!--       Prognostic equation for rain water content
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   qr_p(k,j,i) = qr(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tqr_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * qr(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( qr_p(k,j,i) < 0.0_wp )  qr_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tqr_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqr_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
!
!--       Calculate prognostic equation for rain drop concentration
          sbt = tsc(2)
          IF ( scalar_advec == 'bc-scheme' )  THEN

             IF ( timestep_scheme(1:5) /= 'runge' )  THEN
!
!--             Bott-Chlond scheme always uses Euler time step. Thus:
                sbt = 1.0_wp
             ENDIF
             tend = 0.0_wp
             CALL advec_s_bc( nr, 'nr' )

          ENDIF

!
!--       nr-tendency terms with no communication
          IF ( scalar_advec /= 'bc-scheme' )  THEN
             tend = 0.0_wp
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( ws_scheme_sca )  THEN
                   IF ( .NOT. advanced_div_correction )  THEN
                      CALL advec_s_ws( advc_flags_s, nr, 'nr',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s )
                   ELSE
                      CALL advec_s_ws( advc_flags_s, nr, 'nr',                                     &
                                       bc_dirichlet_l  .OR.  bc_radiation_l,                       &
                                       bc_dirichlet_n  .OR.  bc_radiation_n,                       &
                                       bc_dirichlet_r  .OR.  bc_radiation_r,                       &
                                       bc_dirichlet_s  .OR.  bc_radiation_s,                       &
                                       advanced_div_correction )
                   ENDIF
                ELSE
                   CALL advec_s_pw( nr )
                ENDIF
             ELSE
                CALL advec_s_up( nr )
             ENDIF
          ENDIF

          CALL diffusion_s( nr, surf_top%nrsws, surf_def%nrsws, surf_lsm%nrsws, surf_usm%nrsws )
!
!--       Prognostic equation for rain drop concentration
          DO  i = nxl, nxr
             DO  j = nys, nyn
                !following directive is required to vectorize on Intel19
                !DIR$ IVDEP
                DO  k = nzb+1, nzt
                   nr_p(k,j,i) = nr(k,j,i) + ( dt_3d * ( sbt * tend(k,j,i) +                       &
                                                         tsc(3) * tnr_m(k,j,i) )                   &
                                               - tsc(5) * rdf_sc(k) * nr(k,j,i)                    &
                                             )                                                     &
                                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                   IF ( nr_p(k,j,i) < 0.0_wp )  nr_p(k,j,i) = 0.0_wp
                ENDDO
             ENDDO
          ENDDO

!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnr_m(k,j,i) = tend(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb+1, nzt
                         tnr_m(k,j,i) =  -9.5625_wp * tend(k,j,i)             &
                                         + 5.3125_wp * tnr_m(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
       ENDIF

       CALL cpu_log( log_point(67), 'all_bmc-prog_equation', 'stop' )

    END SUBROUTINE bcm_prognostic_equations


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Control of microphysics for grid points i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_prognostic_equations_ij( i, j, i_omp_start, tn )


       INTEGER(iwp), INTENT(IN) ::  i            !< grid index in x-direction
       INTEGER(iwp), INTENT(IN) ::  i_omp_start  !< first loop index of i-loop in
                                                 !< prognostic_equations
       INTEGER(iwp), INTENT(IN) ::  j            !< grid index in y-direction
       INTEGER(iwp)             ::  k            !< grid index in z-direction
       INTEGER(iwp), INTENT(IN) ::  tn           !< task number of openmp task

!
!--    If required, calculate prognostic equations for cloud water content and cloud drop
!--    concentration.
       IF ( microphysics_morrison )  THEN
!
!--       Calculate prognostic equation for cloud water content
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, qc, 'qc', flux_s_qc,                          &
                                 diss_s_qc, flux_l_qc, diss_l_qc,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, qc )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, qc )
          ENDIF
          CALL diffusion_s( i, j, qc, surf_top%qcsws, surf_def%qcsws, surf_lsm%qcsws,              &
                            surf_usm%qcsws )
!
!--       Prognostic equation for cloud water content
          DO  k = nzb+1, nzt
             qc_p(k,j,i) = qc(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tqc_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * qc(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( qc_p(k,j,i) < 0.0_wp )  qc_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tqc_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tqc_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqc_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

!
!--       Calculate prognostic equation for cloud drop concentration.
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, nc, 'nc', flux_s_nc,                          &
                                 diss_s_nc, flux_l_nc, diss_l_nc,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, nc )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, nc )
          ENDIF
          CALL diffusion_s( i, j, nc, surf_top%ncsws, surf_def%ncsws, surf_lsm%ncsws,              &
                            surf_usm%ncsws )
!
!--       Prognostic equation for cloud drop concentration
          DO  k = nzb+1, nzt
             nc_p(k,j,i) = nc(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tnc_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * nc(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( nc_p(k,j,i) < 0.0_wp )  nc_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tnc_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tnc_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tnc_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

       ENDIF

!
!--    If required, calculate prognostic equations for ice crystal mixing ratio and ice crystal
!--    concentration.
       IF ( microphysics_ice_phase )  THEN
!
!--       Calculate prognostic equation for ice crystal mixing ratio
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, qi, 'qi', flux_s_qi,                          &
                                 diss_s_qi, flux_l_qi, diss_l_qi,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, qi )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, qi )
          ENDIF
          CALL diffusion_s( i, j, qi, surf_top%qisws, surf_def%qisws, surf_lsm%qisws,              &
                            surf_usm%qisws )
!
!--       Prognostic equation for ice crystal mixing ratio
          DO  k = nzb+1, nzt
             qi_p(k,j,i) = qi(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tqi_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * qi(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( qi_p(k,j,i) < 0.0_wp )  qi_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tqi_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tqi_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqi_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

!
!--       Calculate prognostic equation for ice crystal concentration.
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, ni, 'ni', flux_s_ni,                          &
                                 diss_s_ni, flux_l_ni, diss_l_ni,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, ni )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, ni )
          ENDIF
          CALL diffusion_s( i, j, ni, surf_top%nisws, surf_def%nisws, surf_lsm%nisws,              &
                            surf_usm%nisws )
!
!--       Prognostic equation for ice crystal concentration
          DO  k = nzb+1, nzt
             ni_p(k,j,i) = ni(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tni_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * ni(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( ni_p(k,j,i) < 0.0_wp )  ni_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tni_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tni_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tni_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

       ENDIF

!
!--    If required, calculate prognostic equations for ice crystal mixing ratio and ice crystal
!--    concentration.
       IF ( microphysics_ice_phase  .AND.  snow )  THEN
!
!--       Calculate prognostic equation for ice crystal mixing ratio
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, qs, 'qs', flux_s_qs,                          &
                                 diss_s_qs, flux_l_qs, diss_l_qs,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, qs )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, qs )
          ENDIF
          CALL diffusion_s( i, j, qs, surf_top%qisws, surf_def%qisws, surf_lsm%qisws,              &
                            surf_usm%qisws )
!
!--       Prognostic equation for snow  mixing ratio
          DO  k = nzb+1, nzt
             qs_p(k,j,i) = qs(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tqs_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * qs(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( qs_p(k,j,i) < 0.0_wp )  qs_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tqs_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tqs_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqs_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

!
!--       Calculate prognostic equation for snow number concentration.
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, ns, 'ns', flux_s_ns,                          &
                                 diss_s_ns, flux_l_ns, diss_l_ns,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, ns )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, ns )
          ENDIF
          CALL diffusion_s( i, j, ns, surf_top%nisws, surf_def%nisws, surf_lsm%nisws,              &
                            surf_usm%nisws )
!
!--       Prognostic equation for snow number concentration
          DO  k = nzb+1, nzt
             ns_p(k,j,i) = ns(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tns_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * ns(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( ns_p(k,j,i) < 0.0_wp )  ns_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tns_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tns_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tns_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

       ENDIF

!
!--    If required, calculate prognostic equations for graupel mixing ratio and graupel number
!--    concentration.
       IF ( microphysics_ice_phase  .AND.  graupel )  THEN
!
!--       Calculate prognostic equation for graupel mixing ratio
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, qg, 'qg', flux_s_qg,                          &
                                 diss_s_qg, flux_l_qg, diss_l_qg,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, qg )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, qg )
          ENDIF
          CALL diffusion_s( i, j, qg, surf_top%qisws, surf_def%qisws, surf_lsm%qisws,              &
                            surf_usm%qisws )
!
!--       Prognostic equation for ice crystal mixing ratio
          DO  k = nzb+1, nzt
             qg_p(k,j,i) = qg(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tqg_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * qg(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( qg_p(k,j,i) < 0.0_wp )  qg_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tqg_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tqg_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqg_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

!
!--       Calculate prognostic equation for graupel number concentration.
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, ng, 'ng', flux_s_ng,                          &
                                 diss_s_ng, flux_l_ng, diss_l_ng,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, ng )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, ng )
          ENDIF
          CALL diffusion_s( i, j, ng, surf_top%nisws, surf_def%nisws, surf_lsm%nisws,              &
                            surf_usm%nisws )
!
!--       Prognostic equation for graupel number concentration
          DO  k = nzb+1, nzt
             ng_p(k,j,i) = ng(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tng_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * ng(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( ng_p(k,j,i) < 0.0_wp )  ng_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tng_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tng_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tng_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

       ENDIF

!
!--    If required, calculate prognostic equations for rain water content and rain drop
!--    concentration
       IF ( microphysics_seifert )  THEN
!
!--       Calculate prognostic equation for rain water content
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, qr, 'qr', flux_s_qr,                          &
                                 diss_s_qr, flux_l_qr, diss_l_qr,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, qr )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, qr )
          ENDIF
          CALL diffusion_s( i, j, qr, surf_top%qrsws, surf_def%qrsws, surf_lsm%qrsws,              &
                            surf_usm%qrsws )
!
!--       Prognostic equation for rain water content
          DO  k = nzb+1, nzt
             qr_p(k,j,i) = qr(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tqr_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * qr(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( qr_p(k,j,i) < 0.0_wp )  qr_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tqr_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tqr_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tqr_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

!
!--       Calculate prognostic equation for rain drop concentration.
          tend(:,j,i) = 0.0_wp
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( ws_scheme_sca )  THEN
                CALL advec_s_ws( advc_flags_s, i, j, nr, 'nr', flux_s_nr,                          &
                                 diss_s_nr, flux_l_nr, diss_l_nr,                                  &
                                 i_omp_start, tn,                                                  &
                                 bc_dirichlet_l  .OR.  bc_radiation_l,                             &
                                 bc_dirichlet_n  .OR.  bc_radiation_n,                             &
                                 bc_dirichlet_r  .OR.  bc_radiation_r,                             &
                                 bc_dirichlet_s  .OR.  bc_radiation_s )
             ELSE
                CALL advec_s_pw( i, j, nr )
             ENDIF
          ELSE
             CALL advec_s_up( i, j, nr )
          ENDIF
          CALL diffusion_s( i, j, nr, surf_top%nrsws, surf_def%nrsws, surf_lsm%nrsws,              &
                            surf_usm%nrsws )
!
!--       Prognostic equation for rain drop concentration
          DO  k = nzb+1, nzt
             nr_p(k,j,i) = nr(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) + tsc(3) * tnr_m(k,j,i) )  &
                                         - tsc(5) * rdf_sc(k) * nr(k,j,i)                          &
                                       ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             IF ( nr_p(k,j,i) < 0.0_wp )  nr_p(k,j,i) = 0.0_wp
          ENDDO
!
!--       Calculate tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
             IF ( intermediate_timestep_count == 1 )  THEN
                DO  k = nzb+1, nzt
                   tnr_m(k,j,i) = tend(k,j,i)
                ENDDO
             ELSEIF ( intermediate_timestep_count < &
                      intermediate_timestep_count_max )  THEN
                DO  k = nzb+1, nzt
                   tnr_m(k,j,i) =   -9.5625_wp * tend(k,j,i) + 5.3125_wp * tnr_m(k,j,i)
                ENDDO
             ENDIF
          ENDIF

       ENDIF

    END SUBROUTINE bcm_prognostic_equations_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_swap_timelevel ( mod_count )

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: mod_count

       IF ( bulk_cloud_model )  THEN

          SELECT CASE ( mod_count )

             CASE ( 0 )

                IF ( microphysics_morrison )  THEN
                   qc => qc_1;    qc_p => qc_2
                   nc => nc_1;    nc_p => nc_2
                ENDIF
                IF ( microphysics_seifert )  THEN
                   qr => qr_1;    qr_p => qr_2
                   nr => nr_1;    nr_p => nr_2
                ENDIF
                IF ( microphysics_ice_phase )  THEN
                   qi => qi_1;    qi_p => qi_2
                   ni => ni_1;    ni_p => ni_2
                   IF ( snow )  THEN
                      qs => qs_1;    qs_p => qs_2
                      ns => ns_1;    ns_p => ns_2
                   ENDIF
                   IF ( graupel )  THEN
                      qg => qg_1;    qg_p => qg_2
                      ng => ng_1;    ng_p => ng_2
                   ENDIF
                ENDIF


             CASE ( 1 )

                IF ( microphysics_morrison )  THEN
                   qc => qc_2;    qc_p => qc_1
                   nc => nc_2;    nc_p => nc_1
                ENDIF
                IF ( microphysics_seifert )  THEN
                   qr => qr_2;    qr_p => qr_1
                   nr => nr_2;    nr_p => nr_1
                ENDIF
                IF ( microphysics_ice_phase )  THEN
                   qi => qi_2;    qi_p => qi_1
                   ni => ni_2;    ni_p => ni_1
                   IF ( snow )  THEN
                      qs => qs_2;    qs_p => qs_1
                      ns => ns_2;    ns_p => ns_1
                   ENDIF
                   IF ( graupel )  THEN
                      qg => qg_2;    qg_p => qg_1
                      ng => ng_2;    ng_p => ng_1
                   ENDIF
                ENDIF


          END SELECT

       ENDIF

    END SUBROUTINE bcm_swap_timelevel


!--------------------------------------------------------------------------------------------------!
! Description: Boundary conditions of the bulk cloud module variables
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_boundary_conditions

       IMPLICIT NONE

       INTEGER(iwp) ::  i !<
       INTEGER(iwp) ::  j !<
       INTEGER(iwp) ::  k !<
       INTEGER(iwp) ::  m !<

       IF ( microphysics_morrison )  THEN
!
!--       Surface conditions cloud water (Dirichlet).
!--       Run loop over all non-natural and natural walls. Note, in wall-datatype the i,j,k coordinates
!--       belong to the atmospheric grid point, therefore, set qr_p and nr_p at
!--       (k+koff,j+joff,i+ioff) walls.
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m)
             j = bc_hv%j(m)
             k = bc_hv%k(m)
             qc_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
             nc_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
          ENDDO
!
!--       Top boundary condition for cloud water (Dirichlet)
          qc_p(nzt+1,:,:) = 0.0_wp
          nc_p(nzt+1,:,:) = 0.0_wp

       ENDIF

       IF ( microphysics_ice_phase )  THEN
!
!--       Surface conditions ice crysral, snow and graupel (Dirichlet).
!--       Run loop over all non-natural and natural walls. Note, in wall-datatype the i,j,k coordinates
!--       belong to the atmospheric grid point, therefore, set qr_p and nr_p at
!--       (k+koff,j+joff,i+ioff) walls.
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m)
             j = bc_hv%j(m)
             k = bc_hv%k(m)
             qi_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
             ni_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
          ENDDO
!
!--       Top boundary condition for ice crystal (Dirichlet)
          qi_p(nzt+1,:,:) = 0.0_wp
          ni_p(nzt+1,:,:) = 0.0_wp

          IF ( snow )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_hv%ns
                i = bc_hv%i(m)
                j = bc_hv%j(m)
                k = bc_hv%k(m)
                qs_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
                ns_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
             ENDDO
!
!--          Top boundary condition for snow (Dirichlet)
             qs_p(nzt+1,:,:) = 0.0_wp
             ns_p(nzt+1,:,:) = 0.0_wp
          ENDIF
          IF ( graupel )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_hv%ns
                i = bc_hv%i(m)
                j = bc_hv%j(m)
                k = bc_hv%k(m)
                qg_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
                ng_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
             ENDDO
!
!--          Top boundary condition for snow (Dirichlet)
             qg_p(nzt+1,:,:) = 0.0_wp
             ng_p(nzt+1,:,:) = 0.0_wp
          ENDIF

       ENDIF

       IF ( microphysics_seifert )  THEN
!
!--       Surface conditions rain water (Dirichlet).
!--       Run loop over all non-natural and natural walls. Note, in wall-datatype the i,j,k coordinates
!--       belong to the atmospheric grid point, therefore, set qr_p and nr_p at
!--       (k+koff,j+joff,i+ioff) walls.
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m)
             j = bc_hv%j(m)
             k = bc_hv%k(m)
             qr_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
             nr_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
          ENDDO
!
!--       Top boundary condition for rain water (Dirichlet)
          qr_p(nzt+1,:,:) = 0.0_wp
          nr_p(nzt+1,:,:) = 0.0_wp

       ENDIF

!
!--    Lateral boundary conditions for scalar quantities at the outflow.
!--    Lateral oundary conditions for TKE and dissipation are set in tcm_boundary_conds.
       IF ( bc_radiation_s )  THEN
          IF ( microphysics_morrison )  THEN
             qc_p(:,nys-1,:) = qc_p(:,nys,:)
             nc_p(:,nys-1,:) = nc_p(:,nys,:)
          ENDIF
          IF ( microphysics_seifert )  THEN
             qr_p(:,nys-1,:) = qr_p(:,nys,:)
             nr_p(:,nys-1,:) = nr_p(:,nys,:)
          ENDIF
          IF ( microphysics_ice_phase )  THEN
             qi_p(:,nys-1,:) = qi_p(:,nys,:)
             ni_p(:,nys-1,:) = ni_p(:,nys,:)
             IF ( snow )  THEN
                qs_p(:,nys-1,:) = qs_p(:,nys,:)
                ns_p(:,nys-1,:) = ns_p(:,nys,:)
             ENDIF
             IF ( graupel )  THEN
                qg_p(:,nys-1,:) = qg_p(:,nys,:)
                ng_p(:,nys-1,:) = ng_p(:,nys,:)
             ENDIF
          ENDIF
       ELSEIF ( bc_radiation_n )  THEN
          IF ( microphysics_morrison )  THEN
             qc_p(:,nyn+1,:) = qc_p(:,nyn,:)
             nc_p(:,nyn+1,:) = nc_p(:,nyn,:)
          ENDIF
          IF ( microphysics_seifert )  THEN
             qr_p(:,nyn+1,:) = qr_p(:,nyn,:)
             nr_p(:,nyn+1,:) = nr_p(:,nyn,:)
          ENDIF
          IF ( microphysics_ice_phase )  THEN
             qi_p(:,nyn+1,:) = qi_p(:,nyn,:)
             ni_p(:,nyn+1,:) = ni_p(:,nyn,:)
             IF ( snow )  THEN
                qs_p(:,nyn+1,:) = qs_p(:,nyn,:)
                ns_p(:,nyn+1,:) = ns_p(:,nyn,:)
             ENDIF
             IF ( graupel )  THEN
                qg_p(:,nyn+1,:) = qg_p(:,nyn,:)
                ng_p(:,nyn+1,:) = ng_p(:,nyn,:)
             ENDIF
          ENDIF
       ELSEIF ( bc_radiation_l )  THEN
          IF ( microphysics_morrison )  THEN
             qc_p(:,:,nxl-1) = qc_p(:,:,nxl)
             nc_p(:,:,nxl-1) = nc_p(:,:,nxl)
          ENDIF
          IF ( microphysics_seifert )  THEN
             qr_p(:,:,nxl-1) = qr_p(:,:,nxl)
             nr_p(:,:,nxl-1) = nr_p(:,:,nxl)
          ENDIF
          IF ( microphysics_ice_phase )  THEN
             qi_p(:,:,nxl-1) = qi_p(:,:,nxl)
             ni_p(:,:,nxl-1) = ni_p(:,:,nxl)
             IF ( snow )  THEN
                qs_p(:,:,nxl-1) = qs_p(:,:,nxl)
                ns_p(:,:,nxl-1) = ns_p(:,:,nxl)
             ENDIF
             IF ( graupel )  THEN
                qg_p(:,:,nxl-1) = qg_p(:,:,nxl)
                ng_p(:,:,nxl-1) = ng_p(:,:,nxl)
             ENDIF
          ENDIF
       ELSEIF ( bc_radiation_r )  THEN
          IF ( microphysics_morrison )  THEN
             qc_p(:,:,nxr+1) = qc_p(:,:,nxr)
             nc_p(:,:,nxr+1) = nc_p(:,:,nxr)
          ENDIF
          IF ( microphysics_seifert )  THEN
             qr_p(:,:,nxr+1) = qr_p(:,:,nxr)
             nr_p(:,:,nxr+1) = nr_p(:,:,nxr)
          ENDIF
          IF ( microphysics_ice_phase )  THEN
             qi_p(:,:,nxr+1) = qi_p(:,:,nxr)
             ni_p(:,:,nxr+1) = ni_p(:,:,nxr)
             IF ( snow )  THEN
                qs_p(:,:,nxr+1) = qs_p(:,:,nxr)
                ns_p(:,:,nxr+1) = ns_p(:,:,nxr)
             ENDIF
             IF ( graupel )  THEN
                qg_p(:,:,nxr+1) = qg_p(:,:,nxr)
                ng_p(:,:,nxr+1) = ng_p(:,:,nxr)
             ENDIF
          ENDIF
       ENDIF

    END SUBROUTINE bcm_boundary_conditions

!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_3d_data_averaging( mode, variable )

       USE control_parameters,                                                                     &
           ONLY:  average_count_3d

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  mode     !<
       CHARACTER (LEN=*) ::  variable !<

       INTEGER(iwp) ::  i       !< local index
       INTEGER(iwp) ::  j       !< local index
       INTEGER(iwp) ::  k       !< local index

       IF ( mode == 'allocate' )  THEN

          SELECT CASE ( TRIM( variable ) )

             CASE ( 'nc' )
                IF ( .NOT. ALLOCATED( nc_av ) )  THEN
                   ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                nc_av = 0.0_wp

             CASE ( 'ng' )
                IF ( .NOT. ALLOCATED( ng_av ) )  THEN
                   ALLOCATE( ng_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ng_av = 0.0_wp

             CASE ( 'ni' )
                IF ( .NOT. ALLOCATED( ni_av ) )  THEN
                   ALLOCATE( ni_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ni_av = 0.0_wp

             CASE ( 'nr' )
                IF ( .NOT. ALLOCATED( nr_av ) )  THEN
                   ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                nr_av = 0.0_wp

             CASE ( 'ns' )
                IF ( .NOT. ALLOCATED( ns_av ) )  THEN
                   ALLOCATE( ns_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ns_av = 0.0_wp

             CASE ( 'prr' )
                IF ( .NOT. ALLOCATED( prr_av ) )  THEN
                   ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                prr_av = 0.0_wp

             CASE ( 'prr_cloud' )
                IF ( .NOT. ALLOCATED( prr_cloud_av ) )  THEN
                   ALLOCATE( prr_cloud_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                prr_cloud_av = 0.0_wp

             CASE ( 'prr_graupel' )
                IF ( .NOT. ALLOCATED( prr_graupel_av ) )  THEN
                   ALLOCATE( prr_graupel_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                prr_graupel_av = 0.0_wp

             CASE ( 'prr_ice' )
                IF ( .NOT. ALLOCATED( prr_ice_av ) )  THEN
                   ALLOCATE( prr_ice_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                prr_ice_av = 0.0_wp

             CASE ( 'prr_rain' )
                IF ( .NOT. ALLOCATED( prr_rain_av ) )  THEN
                   ALLOCATE( prr_rain_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                prr_rain_av = 0.0_wp

             CASE ( 'prr_snow' )
                IF ( .NOT. ALLOCATED( prr_snow_av ) )  THEN
                   ALLOCATE( prr_snow_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                prr_snow_av = 0.0_wp

             CASE ( 'qc' )
                IF ( .NOT. ALLOCATED( qc_av ) )  THEN
                   ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qc_av = 0.0_wp

             CASE ( 'qg' )
                IF ( .NOT. ALLOCATED( qg_av ) )  THEN
                   ALLOCATE( qg_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qg_av = 0.0_wp

             CASE ( 'qi' )
                IF ( .NOT. ALLOCATED( qi_av ) )  THEN
                   ALLOCATE( qi_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qi_av = 0.0_wp

             CASE ( 'ql' )
                IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                   ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                ql_av = 0.0_wp

             CASE ( 'qr' )
                IF ( .NOT. ALLOCATED( qr_av ) )  THEN
                   ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qr_av = 0.0_wp

             CASE ( 'qs' )
                IF ( .NOT. ALLOCATED( qs_av ) )  THEN
                   ALLOCATE( qs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                qs_av = 0.0_wp

             CASE DEFAULT
                CONTINUE

          END SELECT

       ELSEIF ( mode == 'sum' )  THEN

          SELECT CASE ( TRIM( variable ) )

             CASE ( 'nc' )
                IF ( ALLOCATED( nc_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            nc_av(k,j,i) = nc_av(k,j,i) + nc(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ng' )
                IF ( ALLOCATED( ng_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ng_av(k,j,i) = ng_av(k,j,i) + ng(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ni' )
                IF ( ALLOCATED( ni_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ni_av(k,j,i) = ni_av(k,j,i) + ni(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'nr' )
                IF ( ALLOCATED( nr_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            nr_av(k,j,i) = nr_av(k,j,i) + nr(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ns' )
                IF ( ALLOCATED( ns_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ns_av(k,j,i) = ns_av(k,j,i) + ns(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr' )
                IF ( ALLOCATED( prr_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_av(k,j,i) = prr_av(k,j,i) + prr(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_cloud' )
                IF ( ALLOCATED( prr_cloud_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_cloud_av(k,j,i) = prr_cloud_av(k,j,i) + prr_cloud(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_graupel' )
                IF ( ALLOCATED( prr_graupel_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_graupel_av(k,j,i) = prr_graupel_av(k,j,i) + prr_graupel(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_ice' )
                IF ( ALLOCATED( prr_ice_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_ice_av(k,j,i) = prr_ice_av(k,j,i) + prr_ice(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_rain' )
                IF ( ALLOCATED( prr_rain_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_rain_av(k,j,i) = prr_rain_av(k,j,i) + prr_rain(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_snow' )
                IF ( ALLOCATED( prr_snow_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_snow_av(k,j,i) = prr_snow_av(k,j,i) + prr_snow(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qc' )
                IF ( ALLOCATED( qc_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qc_av(k,j,i) = qc_av(k,j,i) + qc(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qg' )
                IF ( ALLOCATED( qg_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qg_av(k,j,i) = qg_av(k,j,i) + qg(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qi' )
                IF ( ALLOCATED( qi_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qi_av(k,j,i) = qi_av(k,j,i) + qi(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ql' )
                IF ( ALLOCATED( ql_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ql_av(k,j,i) = ql_av(k,j,i) + ql(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qr' )
                IF ( ALLOCATED( qr_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qr_av(k,j,i) = qr_av(k,j,i) + qr(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qs' )
                IF ( ALLOCATED( qs_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qs_av(k,j,i) = qs_av(k,j,i) + qs(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE DEFAULT
                CONTINUE

          END SELECT

       ELSEIF ( mode == 'average' )  THEN

          SELECT CASE ( TRIM( variable ) )

             CASE ( 'nc' )
                IF ( ALLOCATED( nc_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            nc_av(k,j,i) = nc_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ng' )
                IF ( ALLOCATED( ng_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ng_av(k,j,i) = ng_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ni' )
                IF ( ALLOCATED( ni_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ni_av(k,j,i) = ni_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'nr' )
                IF ( ALLOCATED( nr_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            nr_av(k,j,i) = nr_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ns' )
                IF ( ALLOCATED( ns_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ns_av(k,j,i) = ns_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr' )
                IF ( ALLOCATED( prr_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_av(k,j,i) = prr_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_cloud' )
                IF ( ALLOCATED( prr_cloud_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_cloud_av(k,j,i) = prr_cloud_av(k,j,i)                              &
                                                  / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_graupel' )
                IF ( ALLOCATED( prr_graupel_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_graupel_av(k,j,i) = prr_graupel_av(k,j,i)                          &
                                                  / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_ice' )
                IF ( ALLOCATED( prr_ice_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_ice_av(k,j,i) = prr_ice_av(k,j,i)                                  &
                                                  / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_rain' )
                IF ( ALLOCATED( prr_rain_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_rain_av(k,j,i) = prr_rain_av(k,j,i)                                &
                                                  / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'prr_snow' )
                IF ( ALLOCATED( prr_snow_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            prr_snow_av(k,j,i) = prr_snow_av(k,j,i)                                &
                                                  / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qc' )
                IF ( ALLOCATED( qc_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qc_av(k,j,i) = qc_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qg' )
                IF ( ALLOCATED( qg_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qg_av(k,j,i) = qg_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qi' )
                IF ( ALLOCATED( qi_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qi_av(k,j,i) = qi_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'ql' )
                IF ( ALLOCATED( ql_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            ql_av(k,j,i) = ql_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qr' )
                IF ( ALLOCATED( qr_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qr_av(k,j,i) = qr_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

             CASE ( 'qs' )
                IF ( ALLOCATED( qs_av ) )  THEN
                   DO  i = nxlg, nxrg
                      DO  j = nysg, nyng
                         DO  k = nzb, nzt+1
                            qs_av(k,j,i) = qs_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF


             CASE DEFAULT
                CONTINUE

          END SELECT

       ENDIF

    END SUBROUTINE bcm_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 2D output variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bcm_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )


    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(INOUT) ::  grid       !< name of vertical grid
    CHARACTER (LEN=*), INTENT(IN)    ::  mode       !< either 'xy', 'xz' or 'yz'
    CHARACTER (LEN=*), INTENT(IN)    ::  variable   !< name of variable

    INTEGER(iwp), INTENT(IN) ::  av        !< flag for (non-)average output
    INTEGER(iwp), INTENT(IN) ::  nzb_do    !< vertical output index (bottom)
    INTEGER(iwp), INTENT(IN) ::  nzt_do    !< vertical output index (top)

    INTEGER(iwp) ::  flag_nr   !< number of masking flag
    INTEGER(iwp) ::  i         !< loop index along x-direction
    INTEGER(iwp) ::  j         !< loop index along y-direction
    INTEGER(iwp) ::  k         !< loop index along z-direction

    LOGICAL ::  resorted  !< flag if output is already resorted

    LOGICAL, INTENT(INOUT) ::  found   !< flag if output variable is found
    LOGICAL, INTENT(INOUT) ::  two_d   !< flag parameter that indicates 2D variables
                                       !<  (horizontal cross sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf !< local
                                                        !< array to which output data is resorted to

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to selected output variable


    found = .TRUE.
    resorted = .FALSE.

!
!-- Set masking flag for topography for not resorted arrays
    flag_nr = 0    ! 0 = scalar, 1 = u, 2 = v, 3 = w

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'nc_xy', 'nc_xz', 'nc_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => nc
          ELSE
             IF ( .NOT. ALLOCATED( nc_av ) )  THEN
                ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                nc_av = 0.0_wp
             ENDIF
             to_be_resorted => nc_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'ng_xy', 'ng_xz', 'ng_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => ng
          ELSE
             IF ( .NOT. ALLOCATED( ng_av ) )  THEN
                ALLOCATE( ng_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ng_av = 0.0_wp
             ENDIF
             to_be_resorted => ng_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'ni_xy', 'ni_xz', 'ni_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => ni
          ELSE
             IF ( .NOT. ALLOCATED( ni_av ) )  THEN
                ALLOCATE( ni_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ni_av = 0.0_wp
             ENDIF
             to_be_resorted => ni_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'nr_xy', 'nr_xz', 'nr_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => nr
          ELSE
             IF ( .NOT. ALLOCATED( nr_av ) )  THEN
                ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                nr_av = 0.0_wp
             ENDIF
             to_be_resorted => nr_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'ns_xy', 'ns_xz', 'ns_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => ns
          ELSE
             IF ( .NOT. ALLOCATED( ns_av ) )  THEN
                ALLOCATE( ns_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ns_av = 0.0_wp
             ENDIF
             to_be_resorted => ns_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'


       CASE ( 'pra*_xy' )        ! 2d-array / integral quantity => no av
          DO  i = nxl, nxr
             DO  j = nys, nyn
                local_pf(i,j,nzb+1) =  precipitation_amount(j,i)
             ENDDO
          ENDDO
          precipitation_amount = 0.0_wp   ! reset for next integ. interval
          resorted = .TRUE.
          two_d = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu1'

       CASE ( 'prr_xy', 'prr_xz', 'prr_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_av ) )  THEN
                ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_av(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'prr_cloud_xy', 'prr_cloud_xz', 'prr_cloud_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_cloud(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_cloud_av ) )  THEN
                ALLOCATE( prr_cloud_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_cloud_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_cloud_av(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'prr_graupel_xy', 'prr_graupel_xz', 'prr_graupel_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_graupel(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_graupel_av ) )  THEN
                ALLOCATE( prr_graupel_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_graupel_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_graupel_av(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'prr_ice_xy', 'prr_ice_xz', 'prr_ice_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_ice(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_ice_av ) )  THEN
                ALLOCATE( prr_ice_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_ice_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_ice_av(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'prr_rain_xy', 'prr_rain_xz', 'prr_rain_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_rain(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_rain_av ) )  THEN
                ALLOCATE( prr_rain_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_rain_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_rain_av(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'prr_snow_xy', 'prr_snow_xz', 'prr_snow_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_snow(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_snow_av ) )  THEN
                ALLOCATE( prr_snow_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_snow_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = prr_snow_av(k,j,i) * hyrho(nzb+1)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'qc_xy', 'qc_xz', 'qc_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => qc
          ELSE
             IF ( .NOT. ALLOCATED( qc_av ) )  THEN
                ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qc_av = 0.0_wp
             ENDIF
             to_be_resorted => qc_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'qg_xy', 'qg_xz', 'qg_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => qg
          ELSE
             IF ( .NOT. ALLOCATED( qg_av ) )  THEN
                ALLOCATE( qg_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qg_av = 0.0_wp
             ENDIF
             to_be_resorted => qg_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'qi_xy', 'qi_xz', 'qi_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => qi
          ELSE
             IF ( .NOT. ALLOCATED( qi_av ) )  THEN
                ALLOCATE( qi_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qi_av = 0.0_wp
             ENDIF
             to_be_resorted => qi_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'ql_xy', 'ql_xz', 'ql_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => ql
          ELSE
             IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ql_av = 0.0_wp
             ENDIF
             to_be_resorted => ql_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'qr_xy', 'qr_xz', 'qr_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => qr
          ELSE
             IF ( .NOT. ALLOCATED( qr_av ) )  THEN
                ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qr_av = 0.0_wp
             ENDIF
             to_be_resorted => qr_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE ( 'qs_xy', 'qs_xz', 'qs_yz' )
          IF ( av == 0 )  THEN
             to_be_resorted => qs
          ELSE
             IF ( .NOT. ALLOCATED( qs_av ) )  THEN
                ALLOCATE( qs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qs_av = 0.0_wp
             ENDIF
             to_be_resorted => qs_av
          ENDIF
          IF ( mode == 'xy' ) grid = 'zu'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT

    IF ( found .AND. .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                IF ( BTEST( topo_flags(k,j,i), flag_nr ) )  local_pf(i,j,k) = to_be_resorted(k,j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE bcm_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define 3D output variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bcm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) ::  variable   !< name of variable

    INTEGER(iwp), INTENT(IN) ::  av       !< flag for (non-)average output
    INTEGER(iwp), INTENT(IN) ::  nzb_do   !< lower limit of the data output (usually 0)
    INTEGER(iwp), INTENT(IN) ::  nzt_do   !< vertical upper limit of the data output
                                          !< (usually nz_do3d)

    INTEGER(iwp) ::  flag_nr   !< number of masking flag
    INTEGER(iwp) ::  i         !< loop index along x-direction
    INTEGER(iwp) ::  j         !< loop index along y-direction
    INTEGER(iwp) ::  k         !< loop index along z-direction

    LOGICAL      ::  resorted  !< flag if output is already resorted

    LOGICAL, INTENT(INOUT) ::  found     !< flag if output variable is found

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf   !< local
                                                        !< array to which output data is resorted to

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to selected output variable


    found = .TRUE.
    resorted = .FALSE.

!
!-- Set masking flag for topography for not resorted arrays
    flag_nr = 0    ! 0 = scalar, 1 = u, 2 = v, 3 = w

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'nc' )
          IF ( av == 0 )  THEN
             to_be_resorted => nc
          ELSE
             IF ( .NOT. ALLOCATED( nc_av ) ) THEN
                ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                nc_av = 0.0_wp
             ENDIF
             to_be_resorted => nc_av
          ENDIF

       CASE ( 'ng' )
          IF ( av == 0 )  THEN
             to_be_resorted => ng
          ELSE
             IF ( .NOT. ALLOCATED( ng_av ) )  THEN
                ALLOCATE( ng_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ng_av = 0.0_wp
             ENDIF
             to_be_resorted => ng_av
          ENDIF

       CASE ( 'ni' )
          IF ( av == 0 )  THEN
             to_be_resorted => ni
          ELSE
             IF ( .NOT. ALLOCATED( ni_av ) )  THEN
                ALLOCATE( ni_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ni_av = 0.0_wp
             ENDIF
             to_be_resorted => ni_av
          ENDIF

       CASE ( 'nr' )
          IF ( av == 0 )  THEN
             to_be_resorted => nr
          ELSE
             IF ( .NOT. ALLOCATED( nr_av ) )  THEN
                ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                nr_av = 0.0_wp
             ENDIF
             to_be_resorted => nr_av
          ENDIF

       CASE ( 'ns' )
          IF ( av == 0 )  THEN
             to_be_resorted => ns
          ELSE
             IF ( .NOT. ALLOCATED( ns_av ) )  THEN
                ALLOCATE( ns_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ns_av = 0.0_wp
             ENDIF
             to_be_resorted => ns_av
          ENDIF

       CASE ( 'prr' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_av ) )  THEN
                ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.

       CASE ( 'prr_cloud' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_cloud(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_cloud_av ) )  THEN
                ALLOCATE( prr_cloud_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_cloud_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_cloud_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.

       CASE ( 'prr_graupel' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_graupel(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_graupel_av ) )  THEN
                ALLOCATE( prr_graupel_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_graupel_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_graupel_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.

       CASE ( 'prr_ice' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_ice(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_ice_av ) )  THEN
                ALLOCATE( prr_ice_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_ice_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_ice_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.

       CASE ( 'prr_rain' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_rain(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_rain_av ) )  THEN
                ALLOCATE( prr_rain_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_rain_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_rain_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.

       CASE ( 'prr_snow' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_snow(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( prr_snow_av ) )  THEN
                ALLOCATE( prr_snow_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                prr_snow_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = prr_snow_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          resorted = .TRUE.

       CASE ( 'qc' )
          IF ( av == 0 )  THEN
             to_be_resorted => qc
          ELSE
             IF ( .NOT. ALLOCATED( qc_av ) )  THEN
                ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qc_av = 0.0_wp
             ENDIF
             to_be_resorted => qc_av
          ENDIF

       CASE ( 'qg' )
          IF ( av == 0 )  THEN
             to_be_resorted => qg
          ELSE
             IF ( .NOT. ALLOCATED( qg_av ) )  THEN
                ALLOCATE( qg_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qg_av = 0.0_wp
             ENDIF
             to_be_resorted => qg_av
          ENDIF

       CASE ( 'qi' )
          IF ( av == 0 )  THEN
             to_be_resorted => qi
          ELSE
             IF ( .NOT. ALLOCATED( qi_av ) )  THEN
                ALLOCATE( qi_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qi_av = 0.0_wp
             ENDIF
             to_be_resorted => qi_av
          ENDIF

       CASE ( 'ql' )
          IF ( av == 0 )  THEN
             to_be_resorted => ql
          ELSE
             IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ql_av = 0.0_wp
             ENDIF
             to_be_resorted => ql_av
          ENDIF

       CASE ( 'qr' )
          IF ( av == 0 )  THEN
             to_be_resorted => qr
          ELSE
             IF ( .NOT. ALLOCATED( qr_av ) )  THEN
                ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qr_av = 0.0_wp
             ENDIF
             to_be_resorted => qr_av
          ENDIF

       CASE ( 'qs' )
          IF ( av == 0 )  THEN
             to_be_resorted => qs
          ELSE
             IF ( .NOT. ALLOCATED( qs_av ) )  THEN
                ALLOCATE( qs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                qs_av = 0.0_wp
             ENDIF
             to_be_resorted => qs_av
          ENDIF


       CASE DEFAULT
          found = .FALSE.

    END SELECT


    IF ( found .AND. .NOT. resorted )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb_do, nzt_do
                IF ( BTEST( topo_flags(k,j,i), flag_nr ) )  local_pf(i,j,k) = to_be_resorted(k,j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE bcm_data_output_3d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining masked data output
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bcm_data_output_mask( av, variable, found, local_pf, mid )

    CHARACTER(LEN=5) ::  grid  !< flag to distinquish between staggered grids
    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  av              !<
    INTEGER(iwp) ::  i               !<
    INTEGER(iwp) ::  j               !<
    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  mid             !< masked output running index
    INTEGER(iwp) ::  topo_top_index  !< k index of highest horizontal surface

    LOGICAL ::  found     !< true if output array was found
    LOGICAL ::  resorted  !< true if array is resorted

    REAL(wp), DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  local_pf  !<

    REAL(wp), DIMENSION(:,:,:), POINTER ::  to_be_resorted  !< points to array which needs to be resorted for output


    found    = .TRUE.
    grid     = 's'
    resorted = .FALSE.

    SELECT CASE ( TRIM( variable ) )

          CASE ( 'nc' )
             IF ( av == 0 )  THEN
                to_be_resorted => nc
             ELSE
                to_be_resorted => nc_av
             ENDIF

          CASE ( 'ng' )
             IF ( av == 0 )  THEN
                to_be_resorted => ng
             ELSE
                to_be_resorted => ng_av
             ENDIF

          CASE ( 'ni' )
             IF ( av == 0 )  THEN
                to_be_resorted => ni
             ELSE
                to_be_resorted => ni_av
             ENDIF

          CASE ( 'nr' )
             IF ( av == 0 )  THEN
                to_be_resorted => nr
             ELSE
                to_be_resorted => nr_av
             ENDIF

          CASE ( 'ns' )
             IF ( av == 0 )  THEN
                to_be_resorted => ns
             ELSE
                to_be_resorted => ns_av
             ENDIF

          CASE ( 'qc' )
             IF ( av == 0 )  THEN
                to_be_resorted => qc
             ELSE
                to_be_resorted => qc_av
             ENDIF

          CASE ( 'qg' )
             IF ( av == 0 )  THEN
                to_be_resorted => qg
             ELSE
                to_be_resorted => qg_av
             ENDIF

          CASE ( 'qi' )
             IF ( av == 0 )  THEN
                to_be_resorted => qi
             ELSE
                to_be_resorted => qi_av
             ENDIF

          CASE ( 'qr' )
             IF ( av == 0 )  THEN
                to_be_resorted => qr
             ELSE
                to_be_resorted => qr_av
             ENDIF

          CASE ( 'qs' )
             IF ( av == 0 )  THEN
                to_be_resorted => qs
             ELSE
                to_be_resorted => qs_av
             ENDIF

       CASE DEFAULT
          found = .FALSE.

    END SELECT

!
!-- Resort the array to be output, if not done above
    IF ( found  .AND.  .NOT. resorted )  THEN
       IF ( .NOT. mask_surface(mid) )  THEN
!
!--       Default masked output
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
                DO  k = 1, mask_size_l(mid,3)
                   local_pf(i,j,k) =  to_be_resorted(mask_k(mid,k), mask_j(mid,j), mask_i(mid,i))
                ENDDO
             ENDDO
          ENDDO

       ELSE
!
!--       Terrain-following masked output
          DO  i = 1, mask_size_l(mid,1)
             DO  j = 1, mask_size_l(mid,2)
!
!--             Get k index of highest horizontal surface
                topo_top_index = topo_top_ind(mask_j(mid,j),mask_i(mid,i), 0 )
!
!--             Save output array
                DO  k = 1, mask_size_l(mid,3)
                   local_pf(i,j,k) = to_be_resorted( MIN( topo_top_index+mask_k(mid,k), nzt+1 ),   &
                                                     mask_j(mid,j), mask_i(mid,i) )
                ENDDO
             ENDDO
          ENDDO

       ENDIF
    ENDIF

 END SUBROUTINE bcm_data_output_mask


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine computes profile timeseries data for the bulk cloud module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bcm_statistics( mode, sr )

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode  !<

    INTEGER(iwp) ::  k          !<
    INTEGER(iwp) ::  sr         !<

    REAL(wp) ::  clwp            !< cloud water path
    REAL(wp) ::  grwp            !< graupel water path
    REAL(wp) ::  icwp            !< ice water path
    REAL(wp) ::  lqwp            !< liquid water path
    REAL(wp) ::  rawp            !< rain water path
    REAL(wp) ::  snwp            !< snow water path

    IF ( mode == 'time_series' )  THEN
!
!--    Store time series date.
       lqwp = 0.0_wp
       DO  k = nzb+1, nzt
          lqwp  = lqwp + hom(k,1,54,sr) * dzu(k)
       ENDDO
!
!--    set lqwp in units g/m2
       ts_value(dots_start_index_bcm,sr) = lqwp * 1000.0_wp

       clwp = 0.0_wp
          DO  k = nzb+1, nzt
             clwp  = clwp + hom(k,1,75,sr) * dzu(k)
          ENDDO
!
!--    set clwp in units g/m2
       ts_value(dots_start_index_bcm+1,sr) = clwp * 1000.0_wp

       IF ( microphysics_seifert )  THEN
          rawp = 0.0_wp
          DO  k = nzb+1, nzt
             rawp  = rawp + hom(k,1,74,sr) * dzu(k)
          ENDDO
!
!--       set rawp in units g/m2
          ts_value(dots_start_index_bcm+2,sr) = rawp * 1000.0_wp
       ENDIF

       IF ( microphysics_ice_phase )  THEN
          icwp = 0.0_wp
          DO  k = nzb+1, nzt
             icwp  = icwp + hom(k,1,125,sr) * dzu(k)
          ENDDO
!
!--       set icwp in units g/m2
          ts_value(dots_start_index_bcm+3,sr) = icwp * 1000.0_wp

          IF ( snow  .AND.  graupel  )  THEN
             grwp = 0.0_wp
             DO  k = nzb+1, nzt
                grwp  = grwp + hom(k,1,127,sr) * dzu(k)
             ENDDO
!
!--          set grwp in units g/m2
             ts_value(dots_start_index_bcm+4,sr) = grwp * 1000.0_wp
             snwp = 0.0_wp
             DO  k = nzb+1, nzt
                snwp  = snwp + hom(k,1,129,sr) * dzu(k)
             ENDDO
!
!--          set snwp in units g/m2
             ts_value(dots_start_index_bcm+5,sr) = snwp * 1000.0_wp
          ENDIF
       ENDIF
    ENDIF

 END SUBROUTINE bcm_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_rrd_global_ftn( found )


       USE control_parameters,                                                                     &
           ONLY: length,                                                                           &
                 restart_string


       IMPLICIT NONE

       LOGICAL, INTENT(OUT)  ::  found


       found = .TRUE.

       SELECT CASE ( restart_string(1:length) )

          CASE ( 'c_sedimentation' )
             READ ( 13 )  c_sedimentation

          CASE ( 'bulk_cloud_model' )
             READ ( 13 )  bulk_cloud_model

          CASE ( 'cloud_scheme' )
             READ ( 13 )  cloud_scheme

          CASE ( 'cloud_water_sedimentation' )
             READ ( 13 )  cloud_water_sedimentation

          CASE ( 'collision_turbulence' )
             READ ( 13 )  collision_turbulence

          CASE ( 'limiter_sedimentation' )
             READ ( 13 )  limiter_sedimentation

          CASE ( 'nc_const' )
             READ ( 13 )  nc_const

          CASE ( 'precipitation' )
             READ ( 13 ) precipitation

          CASE ( 'ventilation_effect' )
             READ ( 13 )  ventilation_effect

          CASE ( 'na_init' )
             READ ( 13 )  na_init

          CASE ( 'dry_aerosol_radius' )
             READ ( 13 )  dry_aerosol_radius

          CASE ( 'sigma_bulk' )
             READ ( 13 )  sigma_bulk

          CASE ( 'aerosol_bulk' )
             READ ( 13 )  aerosol_bulk

          CASE ( 'curvature_solution_effects_bulk' )
             READ ( 13 )  curvature_solution_effects_bulk

          CASE ( 'microphysics_ice_phase' )
             READ ( 13 )  microphysics_ice_phase

          CASE ( 'ice_crystal_sedimentation' )
             READ ( 13 )  ice_crystal_sedimentation

          CASE ( 'in_init' )
             READ ( 13 )  in_init

          CASE ( 'start_ice_microphysics' )
             READ ( 13 )  start_ice_microphysics

          CASE ( 'snow' )
             READ ( 13 )  snow

          CASE ( 'graupel' )
             READ ( 13 )  graupel

          CASE DEFAULT

             found = .FALSE.

       END SELECT


    END SUBROUTINE bcm_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_rrd_global_mpi

       CALL rrd_mpi_io( 'c_sedimentation', c_sedimentation )
       CALL rrd_mpi_io( 'bulk_cloud_model', bulk_cloud_model )
       CALL rrd_mpi_io( 'cloud_scheme', cloud_scheme )
       CALL rrd_mpi_io( 'cloud_water_sedimentation', cloud_water_sedimentation )
       CALL rrd_mpi_io( 'collision_turbulence', collision_turbulence )
       CALL rrd_mpi_io( 'limiter_sedimentation', limiter_sedimentation )
       CALL rrd_mpi_io( 'nc_const', nc_const )
       CALL rrd_mpi_io( 'precipitation', precipitation )
       CALL rrd_mpi_io( 'ventilation_effect', ventilation_effect )
       CALL rrd_mpi_io( 'na_init', na_init )
       CALL rrd_mpi_io( 'dry_aerosol_radius', dry_aerosol_radius )
       CALL rrd_mpi_io( 'sigma_bulk', sigma_bulk )
       CALL rrd_mpi_io( 'aerosol_bulk', aerosol_bulk )
       CALL rrd_mpi_io( 'curvature_solution_effects_bulk', curvature_solution_effects_bulk )
       CALL rrd_mpi_io( 'start_ice_microphysics', start_ice_microphysics )
       CALL rrd_mpi_io( 'microphysics_ice_phase', microphysics_ice_phase )
       CALL rrd_mpi_io( 'in_init', in_init )
       CALL rrd_mpi_io( 'ice_crystal_sedimentation', ice_crystal_sedimentation )
       CALL rrd_mpi_io( 'snow', snow )
       CALL rrd_mpi_io( 'graupel', graupel )


    END SUBROUTINE bcm_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync, &
                                  nyn_on_file, nysf, nysc, nys_on_file, tmp_2d, tmp_3d, found )


       USE control_parameters

       USE indices

       USE pegrid


       IMPLICIT NONE

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

       LOGICAL, INTENT(OUT)  ::  found

       REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_2d   !<

       REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !<


       found = .TRUE.

       SELECT CASE ( restart_string(1:length) )

          CASE ( 'nc' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             nc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'nc_av' )
             IF ( .NOT. ALLOCATED( nc_av ) )  THEN
                ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             nc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ng' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ng(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ng_av' )
             IF ( .NOT. ALLOCATED( ng_av ) )  THEN
                ALLOCATE( ng_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ng_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ni' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ni(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ni_av' )
             IF ( .NOT. ALLOCATED( ni_av ) )  THEN
                ALLOCATE( ni_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ni_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'nr' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             nr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'nr_av' )
             IF ( .NOT. ALLOCATED( nr_av ) )  THEN
                ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             nr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ns' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ns(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ns_av' )
             IF ( .NOT. ALLOCATED( ns_av ) )  THEN
                ALLOCATE( ns_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ns_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'precipitation_amount' )
             IF ( k == 1 )  READ ( 13 )  tmp_2d
             precipitation_amount(nysc-nbgp:nync+nbgp, nxlc-nbgp:nxrc+nbgp) =                      &
                tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr' )
             IF ( .NOT. ALLOCATED( prr ) )  THEN
                ALLOCATE( prr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                      &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_cloud' )
             IF ( .NOT. ALLOCATED( prr_cloud ) )  THEN
                ALLOCATE( prr_cloud(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_cloud(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_graupel' )
             IF ( .NOT. ALLOCATED( prr_graupel ) )  THEN
                ALLOCATE( prr_graupel(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_graupel(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_ice' )
             IF ( .NOT. ALLOCATED( prr_ice ) )  THEN
                ALLOCATE( prr_ice(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_ice(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                  &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_rain' )
             IF ( .NOT. ALLOCATED( prr_rain ) )  THEN
                ALLOCATE( prr_rain(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_rain(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                 &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_snow' )
             IF ( .NOT. ALLOCATED( prr_snow ) )  THEN
                ALLOCATE( prr_snow(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_snow(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                 &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_av' )
             IF ( .NOT. ALLOCATED( prr_av ) )  THEN
                ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                   &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_cloud_av' )
             IF ( .NOT. ALLOCATED( prr_cloud_av ) )  THEN
                ALLOCATE( prr_cloud_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_cloud_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                             &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_graupel_av' )
             IF ( .NOT. ALLOCATED( prr_graupel_av ) )  THEN
                ALLOCATE( prr_graupel_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_graupel_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_ice_av' )
             IF ( .NOT. ALLOCATED( prr_ice_av ) )  THEN
                ALLOCATE( prr_ice_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_ice_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_rain_av' )
             IF ( .NOT. ALLOCATED( prr_rain_av ) )  THEN
                ALLOCATE( prr_rain_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_rain_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'prr_snow_av' )
             IF ( .NOT. ALLOCATED( prr_snow_av ) )  THEN
                ALLOCATE( prr_snow_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             prr_snow_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qc' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qg' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qg(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qg_av' )
             IF ( .NOT. ALLOCATED( qg_av ) )  THEN
                ALLOCATE( qg_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qg_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qi' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qi(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qi_av' )
             IF ( .NOT. ALLOCATED( qi_av ) )  THEN
                ALLOCATE( qi_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qi_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qc_av' )
             IF ( .NOT. ALLOCATED( qc_av ) )  THEN
                ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qf' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qf(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ql' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ql(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'ql_av' )
             IF ( .NOT. ALLOCATED( ql_av ) )  THEN
                ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             ql_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qr' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qr_av' )
             IF ( .NOT. ALLOCATED( qr_av ) )  THEN
                ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qs' )
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qs(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                       &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

          CASE ( 'qs_av' )
             IF ( .NOT. ALLOCATED( qs_av ) )  THEN
                ALLOCATE( qs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             qs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

!
          CASE DEFAULT

             found = .FALSE.

          END SELECT


    END SUBROUTINE bcm_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_rrd_local_mpi

       LOGICAL ::  array_found  !<


       CALL rd_mpi_io_check_array( 'prr' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( prr ) )  ALLOCATE( prr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'prr', prr )
       ENDIF
       
       CALL rd_mpi_io_check_array( 'prr_cloud' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( prr_cloud ) )  ALLOCATE( prr_cloud(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'prr_cloud', prr_cloud )
       ENDIF

       CALL rd_mpi_io_check_array( 'prr_graupel' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( prr_graupel ) )  THEN
              ALLOCATE( prr_graupel(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          CALL rrd_mpi_io( 'prr_graupel', prr_graupel )
       ENDIF

       CALL rd_mpi_io_check_array( 'prr_ice' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( prr_ice ) )  ALLOCATE( prr_ice(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'prr_ice', prr_ice )
       ENDIF

       CALL rd_mpi_io_check_array( 'prr_rain' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( prr_rain ) )  ALLOCATE( prr_rain(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'prr_rain', prr_rain )
       ENDIF

       CALL rd_mpi_io_check_array( 'prr_snow' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( prr_snow ) )  ALLOCATE( prr_snow(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'prr_snow', prr_snow )
       ENDIF
!
!--    Note, restart input of time-averaged quantities is skipped in case of cyclic-fill
!--    initialization. This case, input of time-averaged data is useless and can lead to faulty
!--    averaging.
       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'prr_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( prr_av ) )  ALLOCATE( prr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'prr_av', prr_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'prr_cloud_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( prr_cloud_av ) )  THEN
                ALLOCATE( prr_cloud_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             CALL rrd_mpi_io( 'prr_cloud_av', prr_cloud_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'prr_graupel_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( prr_graupel_av ) )  THEN
                ALLOCATE( prr_graupel_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             CALL rrd_mpi_io( 'prr_graupel_av', prr_graupel_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'prr_ice_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( prr_ice_av ) )  THEN
                ALLOCATE( prr_ice_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             CALL rrd_mpi_io( 'prr_ice_av', prr_ice_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'prr_rain_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( prr_rain_av ) )  THEN
                ALLOCATE( prr_rain_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             CALL rrd_mpi_io( 'prr_rain_av', prr_rain_av )
          ENDIF

          CALL rd_mpi_io_check_array( 'prr_snow_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( prr_snow_av ) )  THEN
                ALLOCATE( prr_snow_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             CALL rrd_mpi_io( 'prr_snow_av', prr_snow_av )
          ENDIF
       ENDIF

       CALL rd_mpi_io_check_array( 'precipitation_amount' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( precipitation_amount ) )  THEN
             ALLOCATE( precipitation_amount(nysg:nyng,nxlg:nxrg) )
          ENDIF
          CALL rrd_mpi_io( 'precipitation_amount', precipitation_amount )
       ENDIF

       CALL rd_mpi_io_check_array( 'qf' , found = array_found )
       IF ( array_found )  CALL rrd_mpi_io( 'qf', qf )

       CALL rd_mpi_io_check_array( 'ql' , found = array_found )
       IF ( array_found )  CALL rrd_mpi_io( 'ql', ql )

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'ql_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( ql_av ) )  ALLOCATE( ql_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'ql_av', ql_av )
          ENDIF
       ENDIF

       CALL rd_mpi_io_check_array( 'qc' , found = array_found )
       IF ( array_found )  CALL rrd_mpi_io( 'qc', qc )

       IF ( .NOT. cyclic_fill_initialization )  THEN
          CALL rd_mpi_io_check_array( 'qc_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( qc_av ) )  ALLOCATE( qc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'qc_av', qc_av )
          ENDIF
       ENDIF

       IF ( microphysics_morrison )  THEN

          CALL rd_mpi_io_check_array( 'nc' , found = array_found )
          IF ( array_found )  CALL rrd_mpi_io( 'nc', nc )

          IF ( .NOT. cyclic_fill_initialization )  THEN
             CALL rd_mpi_io_check_array( 'nc_av' , found = array_found )
             IF ( array_found )  THEN
                IF ( .NOT. ALLOCATED( nc_av ) )  ALLOCATE( nc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                CALL rrd_mpi_io( 'nc_av', nc_av )
             ENDIF
          ENDIF

       ENDIF

       IF ( microphysics_seifert )  THEN

          CALL rd_mpi_io_check_array( 'nr' , found = array_found )
          IF ( array_found )  CALL rrd_mpi_io( 'nr', nr )

          IF ( .NOT. cyclic_fill_initialization )  THEN
             CALL rd_mpi_io_check_array( 'nr_av' , found = array_found )
             IF ( array_found )  THEN
                IF ( .NOT. ALLOCATED( nr_av ) )  ALLOCATE( nr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                CALL rrd_mpi_io( 'nr_av', nr_av )
             ENDIF
          ENDIF

          CALL rd_mpi_io_check_array( 'qr' , found = array_found )
          IF ( array_found )  CALL rrd_mpi_io( 'qr', qr )

          IF ( .NOT. cyclic_fill_initialization )  THEN
             CALL rd_mpi_io_check_array( 'qr_av' , found = array_found )
             IF ( array_found )  THEN
                IF ( .NOT. ALLOCATED( qr_av ) )  ALLOCATE( qr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                CALL rrd_mpi_io( 'qr_av', qr_av )
             ENDIF
          ENDIF

       ENDIF

       IF ( microphysics_ice_phase )  THEN

          CALL rd_mpi_io_check_array( 'ni' , found = array_found )
          IF ( array_found )  CALL rrd_mpi_io( 'ni', ni )

          IF ( .NOT. cyclic_fill_initialization )  THEN
             CALL rd_mpi_io_check_array( 'ni_av' , found = array_found )
             IF ( array_found )  THEN
                IF ( .NOT. ALLOCATED( ni_av ) )  ALLOCATE( ni_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                CALL rrd_mpi_io( 'ni_av', ni_av )
             ENDIF
          ENDIF

          CALL rd_mpi_io_check_array( 'qi' , found = array_found )
          IF ( array_found )  CALL rrd_mpi_io( 'qi', qi )

          IF ( .NOT. cyclic_fill_initialization )  THEN
             CALL rd_mpi_io_check_array( 'qi_av' , found = array_found )
             IF ( array_found )  THEN
                IF ( .NOT. ALLOCATED( qi_av ) )  ALLOCATE( qi_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                CALL rrd_mpi_io( 'qi_av', qi_av )
             ENDIF
          ENDIF

          IF ( snow .OR. graupel )  THEN
             CALL rd_mpi_io_check_array( 'ng' , found = array_found )
             IF ( array_found )  CALL rrd_mpi_io( 'ng', ng )

             CALL rd_mpi_io_check_array( 'ns' , found = array_found )
             IF ( array_found )  CALL rrd_mpi_io( 'ns', ns )

             CALL rd_mpi_io_check_array( 'qg' , found = array_found )
             IF ( array_found )  CALL rrd_mpi_io( 'qg', qg )

             CALL rd_mpi_io_check_array( 'qs' , found = array_found )
             IF ( array_found )  CALL rrd_mpi_io( 'qs', qs )

             IF ( .NOT. cyclic_fill_initialization )  THEN
                CALL rd_mpi_io_check_array( 'ng_av' , found = array_found )
                IF ( array_found )  THEN
                   IF ( .NOT. ALLOCATED( ng_av ) )  ALLOCATE( ng_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   CALL rrd_mpi_io( 'ng_av', ng_av )
                ENDIF

                CALL rd_mpi_io_check_array( 'ns_av' , found = array_found )
                IF ( array_found )  THEN
                   IF ( .NOT. ALLOCATED( ns_av ) )  ALLOCATE( ns_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   CALL rrd_mpi_io( 'ns_av', ns_av )
                ENDIF

                CALL rd_mpi_io_check_array( 'qg_av' , found = array_found )
                IF ( array_found )  THEN
                   IF ( .NOT. ALLOCATED( qg_av ) )  ALLOCATE( qg_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   CALL rrd_mpi_io( 'qg_av', qg_av )
                ENDIF

                CALL rd_mpi_io_check_array( 'qs_av' , found = array_found )
                IF ( array_found )  THEN
                   IF ( .NOT. ALLOCATED( qs_av ) )  ALLOCATE( qs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   CALL rrd_mpi_io( 'qs_av', qs_av )
                ENDIF
             ENDIF
          ENDIF

       ENDIF


    END SUBROUTINE bcm_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the bulk cloud module.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_wrd_global


       IMPLICIT NONE


       IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

          CALL wrd_write_string( 'c_sedimentation' )
          WRITE ( 14 )  c_sedimentation

          CALL wrd_write_string( 'bulk_cloud_model' )
          WRITE ( 14 )  bulk_cloud_model

          CALL wrd_write_string( 'cloud_scheme' )
          WRITE ( 14 )  cloud_scheme

          CALL wrd_write_string( 'cloud_water_sedimentation' )
          WRITE ( 14 )  cloud_water_sedimentation

          CALL wrd_write_string( 'collision_turbulence' )
          WRITE ( 14 )  collision_turbulence

          CALL wrd_write_string( 'limiter_sedimentation' )
          WRITE ( 14 )  limiter_sedimentation

          CALL wrd_write_string( 'nc_const' )
          WRITE ( 14 )  nc_const

          CALL wrd_write_string( 'precipitation' )
          WRITE ( 14 )  precipitation

          CALL wrd_write_string( 'ventilation_effect' )
          WRITE ( 14 )  ventilation_effect

          CALL wrd_write_string( 'na_init' )
          WRITE ( 14 )  na_init

          CALL wrd_write_string( 'dry_aerosol_radius' )
          WRITE ( 14 )  dry_aerosol_radius

          CALL wrd_write_string( 'sigma_bulk' )
          WRITE ( 14 )  sigma_bulk

          CALL wrd_write_string( 'aerosol_bulk' )
          WRITE ( 14 )  aerosol_bulk

          CALL wrd_write_string( 'curvature_solution_effects_bulk' )
          WRITE ( 14 )  curvature_solution_effects_bulk

          CALL wrd_write_string( 'start_ice_microphysics' )
          WRITE ( 14 )  start_ice_microphysics

          CALL wrd_write_string( 'microphysics_ice_phase' )
          WRITE ( 14 )  microphysics_ice_phase

          CALL wrd_write_string( 'in_init' )
          WRITE ( 14 )  in_init

          CALL wrd_write_string( 'ice_crystal_sedimentation' )
          WRITE ( 14 )  ice_crystal_sedimentation

          CALL wrd_write_string( 'snow' )
          WRITE ( 14 )  snow

          CALL wrd_write_string( 'graupel' )
          WRITE ( 14 )  graupel

       ELSEIF ( TRIM( restart_data_format_output(1:3) ) == 'mpi' )  THEN

          CALL wrd_mpi_io( 'c_sedimentation', c_sedimentation )
          CALL wrd_mpi_io( 'bulk_cloud_model', bulk_cloud_model )
          CALL wrd_mpi_io( 'cloud_scheme', cloud_scheme )
          CALL wrd_mpi_io( 'cloud_water_sedimentation', cloud_water_sedimentation )
          CALL wrd_mpi_io( 'collision_turbulence', collision_turbulence )
          CALL wrd_mpi_io( 'limiter_sedimentation', limiter_sedimentation )
          CALL wrd_mpi_io( 'nc_const', nc_const )
          CALL wrd_mpi_io( 'precipitation', precipitation )
          CALL wrd_mpi_io( 'ventilation_effect', ventilation_effect )
          CALL wrd_mpi_io( 'na_init', na_init )
          CALL wrd_mpi_io( 'dry_aerosol_radius', dry_aerosol_radius )
          CALL wrd_mpi_io( 'sigma_bulk', sigma_bulk )
          CALL wrd_mpi_io( 'aerosol_bulk', aerosol_bulk )
          CALL wrd_mpi_io( 'curvature_solution_effects_bulk', curvature_solution_effects_bulk )
          CALL wrd_mpi_io( 'start_ice_microphysics', start_ice_microphysics )
          CALL wrd_mpi_io( 'microphysics_ice_phase', microphysics_ice_phase )
          CALL wrd_mpi_io( 'in_init', in_init )
          CALL wrd_mpi_io( 'ice_crystal_sedimentation', ice_crystal_sedimentation )
          CALL wrd_mpi_io( 'snow', snow )
          CALL wrd_mpi_io( 'graupel', graupel )

       ENDIF

    END SUBROUTINE bcm_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the bulk cloud module.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE bcm_wrd_local


       IMPLICIT NONE


       IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

          IF ( ALLOCATED( prr ) )  THEN
             CALL wrd_write_string( 'prr' )
             WRITE ( 14 )  prr
          ENDIF

          IF ( ALLOCATED( prr_cloud ) )  THEN
             CALL wrd_write_string( 'prr_cloud' )
             WRITE ( 14 )  prr_cloud
          ENDIF

          IF ( ALLOCATED( prr_graupel ) )  THEN
             CALL wrd_write_string( 'prr_graupel' )
             WRITE ( 14 )  prr_graupel
          ENDIF

          IF ( ALLOCATED( prr_ice ) )  THEN
             CALL wrd_write_string( 'prr_ice' )
             WRITE ( 14 )  prr_ice
          ENDIF

          IF ( ALLOCATED( prr_rain ) )  THEN
             CALL wrd_write_string( 'prr_rain' )
             WRITE ( 14 )  prr_rain
          ENDIF

          IF ( ALLOCATED( prr_snow ) )  THEN
             CALL wrd_write_string( 'prr_snow' )
             WRITE ( 14 )  prr_snow
          ENDIF

          IF ( ALLOCATED( prr_av ) )  THEN
             CALL wrd_write_string( 'prr_av' )
             WRITE ( 14 )  prr_av
          ENDIF

          IF ( ALLOCATED( prr_cloud_av ) )  THEN
             CALL wrd_write_string( 'prr_cloud_av' )
             WRITE ( 14 )  prr_cloud_av
          ENDIF

          IF ( ALLOCATED( prr_graupel_av ) )  THEN
             CALL wrd_write_string( 'prr_graupel_av' )
             WRITE ( 14 )  prr_graupel_av
          ENDIF

          IF ( ALLOCATED( prr_ice_av ) )  THEN
             CALL wrd_write_string( 'prr_ice_av' )
             WRITE ( 14 )  prr_ice_av
          ENDIF

          IF ( ALLOCATED( prr_rain_av ) )  THEN
             CALL wrd_write_string( 'prr_rain_av' )
             WRITE ( 14 )  prr_rain_av
          ENDIF

          IF ( ALLOCATED( prr_snow_av ) )  THEN
             CALL wrd_write_string( 'prr_snow_av' )
             WRITE ( 14 )  prr_snow_av
          ENDIF

          IF ( ALLOCATED( precipitation_amount ) )  THEN
             CALL wrd_write_string( 'precipitation_amount' )
             WRITE ( 14 )  precipitation_amount
          ENDIF

          CALL wrd_write_string( 'qf' )
          WRITE ( 14 )  qf

          CALL wrd_write_string( 'ql' )
          WRITE ( 14 )  ql

          IF ( ALLOCATED( ql_av ) )  THEN
             CALL wrd_write_string( 'ql_av' )
             WRITE ( 14 )  ql_av
          ENDIF

          CALL wrd_write_string( 'qc' )
          WRITE ( 14 )  qc

          IF ( ALLOCATED( qc_av ) )  THEN
             CALL wrd_write_string( 'qc_av' )
             WRITE ( 14 )  qc_av
          ENDIF

          IF ( microphysics_morrison )  THEN

             CALL wrd_write_string( 'nc' )
             WRITE ( 14 )  nc

             IF ( ALLOCATED( nc_av ) )  THEN
                CALL wrd_write_string( 'nc_av' )
                WRITE ( 14 )  nc_av
             ENDIF

          ENDIF

          IF ( microphysics_ice_phase )  THEN

             CALL wrd_write_string( 'ni' )
             WRITE ( 14 )  ni

             IF ( ALLOCATED( ni_av ) )  THEN
                CALL wrd_write_string( 'ni_av' )
                WRITE ( 14 )  ni_av
             ENDIF

             CALL wrd_write_string( 'qi' )
             WRITE ( 14 )  qi

             IF ( ALLOCATED( qi_av ) )  THEN
                CALL wrd_write_string( 'qi_av' )
                WRITE ( 14 )  qi_av
             ENDIF

             IF ( snow  .OR.  graupel  )  THEN
                CALL wrd_write_string( 'ng' )
                WRITE ( 14 )  ng
                CALL wrd_write_string( 'ns' )
                WRITE ( 14 )  ns
                CALL wrd_write_string( 'qg' )
                WRITE ( 14 )  qg
                CALL wrd_write_string( 'qs' )
                WRITE ( 14 )  qs

                IF ( ALLOCATED( ng_av ) )  THEN
                   CALL wrd_write_string( 'ng_av' )
                   WRITE ( 14 )  ng_av
                ENDIF
                IF ( ALLOCATED( ns_av ) )  THEN
                   CALL wrd_write_string( 'ns_av' )
                   WRITE ( 14 )  ns_av
                ENDIF
                IF ( ALLOCATED( qg_av ) )  THEN
                   CALL wrd_write_string( 'qg_av' )
                   WRITE ( 14 )  qg_av
                ENDIF
                IF ( ALLOCATED( qs_av ) )  THEN
                   CALL wrd_write_string( 'qs_av' )
                   WRITE ( 14 )  qs_av
                ENDIF
             ENDIF

          ENDIF


          IF ( microphysics_seifert )  THEN

             CALL wrd_write_string( 'nr' )
             WRITE ( 14 )  nr

             IF ( ALLOCATED( nr_av ) )  THEN
                CALL wrd_write_string( 'nr_av' )
                WRITE ( 14 )  nr_av
             ENDIF

             CALL wrd_write_string( 'qr' )
             WRITE ( 14 )  qr

             IF ( ALLOCATED( qr_av ) )  THEN
                CALL wrd_write_string( 'qr_av' )
                WRITE ( 14 )  qr_av
             ENDIF

          ENDIF

       ELSEIF ( TRIM( restart_data_format_output(1:3) ) == 'mpi' )  THEN

          IF ( ALLOCATED( prr ) )  CALL wrd_mpi_io( 'prr', prr )
          IF ( ALLOCATED( prr_cloud ) )  CALL wrd_mpi_io( 'prr_cloud', prr_cloud )
          IF ( ALLOCATED( prr_graupel ) )  CALL wrd_mpi_io( 'prr_graupel', prr_graupel )
          IF ( ALLOCATED( prr_ice ) )  CALL wrd_mpi_io( 'prr_ice', prr_ice )
          IF ( ALLOCATED( prr_rain ) )  CALL wrd_mpi_io( 'prr_rain', prr_rain )
          IF ( ALLOCATED( prr_snow ) )  CALL wrd_mpi_io( 'prr_snow', prr_snow )
          IF ( ALLOCATED( prr_av ) )  CALL wrd_mpi_io( 'prr_av', prr_av )
          IF ( ALLOCATED( prr_cloud_av ) )  CALL wrd_mpi_io( 'prr_cloud_av', prr_cloud_av )
          IF ( ALLOCATED( prr_graupel_av ) )  CALL wrd_mpi_io( 'prr_graupel_av', prr_graupel_av )
          IF ( ALLOCATED( prr_ice_av ) )  CALL wrd_mpi_io( 'prr_ice_av', prr_ice_av )
          IF ( ALLOCATED( prr_rain_av ) )  CALL wrd_mpi_io( 'prr_rain_av', prr_rain_av )
          IF ( ALLOCATED( prr_snow_av ) )  CALL wrd_mpi_io( 'prr_snow_av', prr_snow_av )
          IF ( ALLOCATED( precipitation_amount ) )  THEN
             CALL wrd_mpi_io( 'precipitation_amount', precipitation_amount )
          ENDIF
          CALL wrd_mpi_io( 'qf', qf )
          CALL wrd_mpi_io( 'ql', ql )
          IF ( ALLOCATED( ql_av ) )  CALL wrd_mpi_io( 'ql_av', ql_av )
          CALL wrd_mpi_io( 'qc', qc )
          IF ( ALLOCATED( qc_av ) )  CALL wrd_mpi_io( 'qc_av', qc_av )
          IF ( microphysics_morrison )  THEN
             CALL wrd_mpi_io( 'nc', nc )
             IF ( ALLOCATED( nc_av ) )  CALL wrd_mpi_io( 'nc_av', nc_av )
          ENDIF
          IF ( microphysics_seifert )  THEN
             CALL wrd_mpi_io( 'nr', nr )
             IF ( ALLOCATED( nr_av ) )  CALL wrd_mpi_io( 'nr_av', nr_av )
             CALL wrd_mpi_io( 'qr', qr )
             IF ( ALLOCATED( qr_av ) )  CALL wrd_mpi_io( 'qr_av', qr_av )
          ENDIF
          IF ( microphysics_ice_phase )  THEN
             CALL wrd_mpi_io( 'ni', ni )
             IF ( ALLOCATED( ni_av ) )  CALL wrd_mpi_io( 'ni_av', ni_av )
             CALL wrd_mpi_io( 'qi', qi )
             IF ( ALLOCATED( qi_av ) )  CALL wrd_mpi_io( 'qi_av', qi_av )
             IF ( snow  .OR.  graupel  )  THEN
                CALL wrd_mpi_io( 'ng', ng )
                CALL wrd_mpi_io( 'ns', ns )
                CALL wrd_mpi_io( 'qg', qg )
                CALL wrd_mpi_io( 'qs', qs )
                IF ( ALLOCATED( ng_av ) )  CALL wrd_mpi_io( 'ng_av', ng_av )
                IF ( ALLOCATED( ns_av ) )  CALL wrd_mpi_io( 'ns_av', ns_av )
                IF ( ALLOCATED( qg_av ) )  CALL wrd_mpi_io( 'qg_av', qg_av )
                IF ( ALLOCATED( qs_av ) )  CALL wrd_mpi_io( 'qs_av', qs_av )
             ENDIF
          ENDIF

       ENDIF

    END SUBROUTINE bcm_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> In this routine the cloued coefficients as given in Table 1 in SB2006 are defined. Right now a
!> namelist steering is not implemented. Accordingly, values must be directly changed in the code.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE define_cloud_coefficients
!
!--    Define values corresponding Table 1 in SB2006 for cloud species
       cloud_species%cloud%a     = 0.124_wp
       cloud_species%cloud%b     = 1.0_wp / 3.0_wp
       cloud_species%cloud%alpha = 3.75E5_wp
       cloud_species%cloud%beta  = 2.0_wp / 3.0_wp
       cloud_species%cloud%nu    = 1.0_wp
       cloud_species%cloud%mu    = 1.0_wp
       cloud_species%cloud%a_ven = 1.0_wp
       cloud_species%cloud%b_ven = 1.0_wp
       cloud_species%cloud%cap   = 1.0_wp
       cloud_species%cloud%x_min = 4.2E-15_wp
       cloud_species%cloud%x_max = 2.6E-10_wp
!
!--    Define values corresponding Table 1 in SB2006 for cloud species
       cloud_species%rain%a     = 0.124_wp
       cloud_species%rain%b     = 1.0_wp / 3.0_wp
       cloud_species%rain%alpha = 159.0_wp
       cloud_species%rain%beta  = 0.266_wp
       cloud_species%rain%nu    = -2.0_wp / 3.0_wp
       cloud_species%rain%mu    = 1.0_wp / 3.0_wp
       cloud_species%rain%a_ven = 1.0_wp
       cloud_species%rain%b_ven = 1.0_wp
       cloud_species%rain%cap   = 1.0_wp
       cloud_species%rain%x_min = 2.6E-10_wp
       cloud_species%rain%x_max = 5.0E-6_wp
!
!--    Define values corresponding Table 1 in SB2006 for graupel species
       cloud_species%graupel%a        = 0.19_wp
       cloud_species%graupel%b        = 0.323_wp
       cloud_species%graupel%alpha    = 40.0_wp
       cloud_species%graupel%beta     = 0.23_wp
       cloud_species%graupel%nu       = 1.0_wp
       cloud_species%graupel%mu       = 1.0_wp / 3.0_wp
       cloud_species%graupel%a_ven    = 0.78_wp
       cloud_species%graupel%b_ven    = 0.308_wp
       cloud_species%graupel%cap      = 2.0_wp
       cloud_species%graupel%x_min    = 2.6E-10_wp
       cloud_species%graupel%x_max    = 1.0E-4_wp
       cloud_species%graupel%coll_eff = 0.8_wp
!
!--    Define values corresponding Table 1 in SB2006 for ice species
       cloud_species%ice%a        = 0.217_wp
       cloud_species%ice%b        = 0.302_wp
       cloud_species%ice%alpha    = 317.0_wp
       cloud_species%ice%beta     = 0.363_wp
       cloud_species%ice%nu       = 1.0_wp
       cloud_species%ice%mu       = 1.0_wp / 3.0_wp
       cloud_species%ice%sigma_v  = 0.2_wp
       cloud_species%ice%a_ven    = 0.86_wp
       cloud_species%ice%b_ven    = 0.28_wp
       cloud_species%ice%cap      = 2.0_wp
       cloud_species%ice%x_min    = 1.0E-12_wp
       cloud_species%ice%x_max    = 1.0E-7_wp
       cloud_species%ice%coll_eff = 0.8_wp
!
!--    Define values corresponding Table 1 in SB2006 for snow species
       cloud_species%snow%a        = 8.156_wp
       cloud_species%snow%b        = 0.526_wp
       cloud_species%snow%alpha    = 27.7_wp
       cloud_species%snow%beta     = 0.216_wp
       cloud_species%snow%nu       = 1.0_wp
       cloud_species%snow%mu       = 1.0_wp / 3.0_wp
       cloud_species%snow%sigma_v  = 0.2_wp
       cloud_species%snow%a_ven    = 0.78_wp
       cloud_species%snow%b_ven    = 0.308_wp
       cloud_species%snow%cap      = 2.0_wp
       cloud_species%snow%x_min    = 1.73E-9_wp
       cloud_species%snow%x_max    = 1.0E-7_wp
       cloud_species%snow%coll_eff = 0.8_wp

    END SUBROUTINE define_cloud_coefficients


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> In this routine the collection integrals are calculated at the beginning of the simulation to
!> avoid the calculation at each timestep. This dimensionless constants, depending on the diameter-
!> mass and velocity-mass relations of particle type "a" and "b". Quations are given in Seifert and
!> Beheng, 2006 (Appendix C Eq. 90 - 94 ).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE init_collision_integrals

       REAL(wp) ::  delta_n                     !< collision integral
       REAL(wp) ::  theta_n                     !< collision integral

       REAL(wp) ::  delta_n_11                  !< collision integral
       REAL(wp) ::  delta_n_12                  !< collision integral
       REAL(wp) ::  delta_n_22                  !< collision integral
       REAL(wp) ::  theta_n_11                  !< collision integral
       REAL(wp) ::  theta_n_12                  !< collision integral
       REAL(wp) ::  theta_n_22                  !< collision integral
       REAL(wp) ::  delta_q_11                  !< collision integral
       REAL(wp) ::  delta_q_12                  !< collision integral
       REAL(wp) ::  delta_q_22                  !< collision integral
       REAL(wp) ::  theta_q_11                  !< collision integral
       REAL(wp) ::  theta_q_12                  !< collision integral
       REAL(wp) ::  theta_q_22                  !< collision integral
!
!--    Calculate collision integrals, last argument is always the moment:
!--    0-th moment: number concentration
!--    1-th moment: mass mixing ratio
!
!---   For graupel selfcollection
       delta_n_11 = collision_integral_delta_b(cloud_species%graupel, 0)
       theta_n_11 = collision_integral_theta_b(cloud_species%graupel, 0)
       delta_n_12 = collision_integral_delta_ba(cloud_species%graupel, cloud_species%graupel, 0)
       theta_n_12 = collision_integral_theta_ba(cloud_species%graupel, cloud_species%graupel, 0)

       delta_n = ( 2.0_wp * delta_n_11 + delta_n_12 )
       theta_n = ( 2.0_wp * theta_n_11 - theta_n_12 )**0.5

       coll_coeff_graupel_self  = pi / 8.0_wp * delta_n * theta_n
!
!--    For snow selfcollection
       delta_n_11 = collision_integral_delta_b(cloud_species%snow, 0)
       theta_n_11 = collision_integral_theta_b(cloud_species%snow, 0)
       delta_n_12 = collision_integral_delta_ba(cloud_species%snow, cloud_species%snow, 0)
       theta_n_12 = collision_integral_theta_ba(cloud_species%snow, cloud_species%snow, 0)

       coll_coeff_snow_self_delta = ( 2.0_wp * delta_n_11 + delta_n_12 )
       coll_coeff_snow_self_theta = ( 2.0_wp * theta_n_11 - theta_n_12 )
!
!--    Ice selfcollection
       delta_n_11 = collision_integral_delta_b(cloud_species%ice, 0)
       delta_n_22 = collision_integral_delta_b(cloud_species%ice, 0)
       delta_q_11 = collision_integral_delta_b(cloud_species%ice, 0)
       delta_q_22 = collision_integral_delta_b(cloud_species%ice, 1)
       delta_n_12 = collision_integral_delta_ba(cloud_species%ice, cloud_species%ice, 0)
       delta_q_12 = collision_integral_delta_ba(cloud_species%ice, cloud_species%ice, 1)

       theta_n_11 = collision_integral_theta_b(cloud_species%ice, 0)
       theta_n_12 = collision_integral_theta_b(cloud_species%ice, 0)
       theta_n_22 = collision_integral_theta_b(cloud_species%ice, 0)
       theta_q_11 = collision_integral_theta_b(cloud_species%ice, 1)
       theta_q_12 = collision_integral_theta_ba(cloud_species%ice, cloud_species%ice, 0)
       theta_q_22 = collision_integral_theta_ba(cloud_species%ice, cloud_species%ice, 1)

       coll_coeff_ice_self_delta_n = delta_n_11 + delta_n_12 + delta_n_22
       coll_coeff_ice_self_delta_q = delta_q_11 + delta_q_12 + delta_q_22
       coll_coeff_ice_self_theta_n = theta_n_11 - theta_n_12 + theta_n_22
       coll_coeff_ice_self_theta_q = theta_q_11 - theta_q_12 + theta_q_22
!
!--    For particle particle collection ( a + b -> a ), Here: graupel + ice -> graupel
       coll_graupel_ice_delta_n_aa = collision_integral_delta_b(cloud_species%graupel, 0)
       coll_graupel_ice_delta_n_bb = collision_integral_delta_b(cloud_species%ice, 0)
       coll_graupel_ice_delta_n_ab = collision_integral_delta_ba(cloud_species%graupel, cloud_species%ice, 0)
       coll_graupel_ice_delta_q_aa = collision_integral_delta_b(cloud_species%graupel, 0) !@TODO: unclear why 0 here, but checked with ICON code
       coll_graupel_ice_delta_q_bb = collision_integral_delta_b(cloud_species%ice, 1)
       coll_graupel_ice_delta_q_ab = collision_integral_delta_ba(cloud_species%graupel, cloud_species%ice, 1)
       coll_graupel_ice_theta_n_aa = collision_integral_theta_b(cloud_species%graupel, 0)
       coll_graupel_ice_theta_n_bb = collision_integral_theta_b(cloud_species%ice, 0)
       coll_graupel_ice_theta_n_ab = collision_integral_theta_ba(cloud_species%graupel, cloud_species%ice, 0)
       coll_graupel_ice_theta_q_aa = collision_integral_theta_b(cloud_species%graupel, 0) !@TODO: unclear why 0 here, but checked with ICON code
       coll_graupel_ice_theta_q_bb = collision_integral_theta_b(cloud_species%ice, 1)
       coll_graupel_ice_theta_q_ab = collision_integral_theta_ba(cloud_species%graupel, cloud_species%ice, 1)
!
!--    For particle particle collection ( a + b -> a ), Here: snow + ice -> snow
       coll_snow_ice_delta_n_aa = collision_integral_delta_b(cloud_species%snow, 0)
       coll_snow_ice_delta_n_bb = collision_integral_delta_b(cloud_species%ice, 0)
       coll_snow_ice_delta_n_ab = collision_integral_delta_ba(cloud_species%snow, cloud_species%ice, 0)
       coll_snow_ice_delta_q_aa = collision_integral_delta_b(cloud_species%snow, 0) !@TODO: unclear why 0 here, but checked with ICON code
       coll_snow_ice_delta_q_bb = collision_integral_delta_b(cloud_species%ice, 1)
       coll_snow_ice_delta_q_ab = collision_integral_delta_ba(cloud_species%snow, cloud_species%ice, 1)
       coll_snow_ice_theta_n_aa = collision_integral_theta_b(cloud_species%snow, 0)
       coll_snow_ice_theta_n_bb = collision_integral_theta_b(cloud_species%ice, 0)
       coll_snow_ice_theta_n_ab = collision_integral_theta_ba(cloud_species%snow, cloud_species%ice, 0)
       coll_snow_ice_theta_q_aa = collision_integral_theta_b(cloud_species%snow, 0) !@TODO: unclear why 0 here, but checked with ICON code
       coll_snow_ice_theta_q_bb = collision_integral_theta_b(cloud_species%ice, 1)
       coll_snow_ice_theta_q_ab = collision_integral_theta_ba(cloud_species%snow, cloud_species%ice, 1)
!
!--    For particle particle collection ( a + b -> a ), Here: graupel + snow -> graupel
       coll_graupel_snow_delta_n_aa = collision_integral_delta_b(cloud_species%graupel, 0)
       coll_graupel_snow_delta_n_bb = collision_integral_delta_b(cloud_species%snow, 0)
       coll_graupel_snow_delta_n_ab = collision_integral_delta_ba(cloud_species%graupel, cloud_species%snow, 0)
       coll_graupel_snow_delta_q_aa = collision_integral_delta_b(cloud_species%graupel, 0) !@TODO: unclear why 0 here, but checked with ICON code
       coll_graupel_snow_delta_q_bb = collision_integral_delta_b(cloud_species%snow, 1)
       coll_graupel_snow_delta_q_ab = collision_integral_delta_ba(cloud_species%graupel, cloud_species%snow, 1)
       coll_graupel_snow_theta_n_aa = collision_integral_theta_b(cloud_species%graupel, 0)
       coll_graupel_snow_theta_n_bb = collision_integral_theta_b(cloud_species%snow, 0)
       coll_graupel_snow_theta_n_ab = collision_integral_theta_ba(cloud_species%graupel, cloud_species%snow, 0)
       coll_graupel_snow_theta_q_aa = collision_integral_theta_b(cloud_species%graupel, 0) !@TODO: unclear why 0 here, but checked with ICON code
       coll_graupel_snow_theta_q_bb = collision_integral_theta_b(cloud_species%snow, 1)
       coll_graupel_snow_theta_q_ab = collision_integral_theta_ba(cloud_species%graupel, cloud_species%snow, 1)
!
!--    For graupel - cloud riming
       rime_graupel_cloud_delta_n_aa = collision_integral_delta_b(cloud_species%graupel, 0)
       rime_graupel_cloud_delta_n_ab = collision_integral_delta_ba(cloud_species%graupel, cloud_species%cloud, 0)
       rime_graupel_cloud_delta_n_bb = collision_integral_delta_b(cloud_species%cloud, 0)
       rime_graupel_cloud_delta_q_aa = collision_integral_delta_b(cloud_species%graupel, 1) ! mass weighted
       rime_graupel_cloud_delta_q_ab = collision_integral_delta_ba(cloud_species%graupel, cloud_species%cloud, 1)
       rime_graupel_cloud_delta_q_ba = collision_integral_delta_ba(cloud_species%cloud, cloud_species%graupel, 1)
       rime_graupel_cloud_delta_q_bb = collision_integral_delta_b(cloud_species%cloud, 1)

       rime_graupel_cloud_theta_n_aa = collision_integral_theta_b(cloud_species%graupel, 0)
       rime_graupel_cloud_theta_n_ab = collision_integral_theta_ba(cloud_species%graupel, cloud_species%cloud, 0)
       rime_graupel_cloud_theta_n_bb = collision_integral_theta_b(cloud_species%cloud, 0)
       rime_graupel_cloud_theta_q_aa = collision_integral_theta_b(cloud_species%graupel, 1) ! mass weighted
       rime_graupel_cloud_theta_q_ab = collision_integral_theta_ba(cloud_species%graupel, cloud_species%cloud, 1)
       rime_graupel_cloud_theta_q_ba = collision_integral_theta_ba(cloud_species%cloud, cloud_species%graupel, 1)
       rime_graupel_cloud_theta_q_bb = collision_integral_theta_b(cloud_species%cloud, 1)
!
!--    For graupel - rain riming
       rime_graupel_rain_delta_n_aa = collision_integral_delta_b(cloud_species%graupel, 0)
       rime_graupel_rain_delta_n_ab = collision_integral_delta_ba(cloud_species%graupel, cloud_species%rain, 0)
       rime_graupel_rain_delta_n_bb = collision_integral_delta_b(cloud_species%rain, 0)
       rime_graupel_rain_delta_q_aa = collision_integral_delta_b(cloud_species%graupel, 1) ! mass weighted
       rime_graupel_rain_delta_q_ab = collision_integral_delta_ba(cloud_species%graupel, cloud_species%rain, 1)
       rime_graupel_rain_delta_q_ba = collision_integral_delta_ba(cloud_species%rain, cloud_species%graupel, 1)
       rime_graupel_rain_delta_q_bb = collision_integral_delta_b(cloud_species%rain, 1)

       rime_graupel_rain_theta_n_aa = collision_integral_theta_b(cloud_species%graupel, 0)
       rime_graupel_rain_theta_n_ab = collision_integral_theta_ba(cloud_species%graupel, cloud_species%rain, 0)
       rime_graupel_rain_theta_n_bb = collision_integral_theta_b(cloud_species%rain, 0)
       rime_graupel_rain_theta_q_aa = collision_integral_theta_b(cloud_species%graupel, 1) ! mass weighted
       rime_graupel_rain_theta_q_ab = collision_integral_theta_ba(cloud_species%graupel, cloud_species%rain, 1)
       rime_graupel_rain_theta_q_ba = collision_integral_theta_ba(cloud_species%rain, cloud_species%graupel, 1)
       rime_graupel_rain_theta_q_bb = collision_integral_theta_b(cloud_species%rain, 1)
!
!--    For ice - cloud riming
       rime_ice_cloud_delta_n_aa = collision_integral_delta_b(cloud_species%ice, 0)
       rime_ice_cloud_delta_n_ab = collision_integral_delta_ba(cloud_species%ice, cloud_species%cloud, 0)
       rime_ice_cloud_delta_n_bb = collision_integral_delta_b(cloud_species%cloud, 0)
       rime_ice_cloud_delta_q_aa = collision_integral_delta_b(cloud_species%ice, 1) ! mass weighted
       rime_ice_cloud_delta_q_ab = collision_integral_delta_ba(cloud_species%ice, cloud_species%cloud, 1)
       rime_ice_cloud_delta_q_ba = collision_integral_delta_ba(cloud_species%cloud, cloud_species%ice, 1)
       rime_ice_cloud_delta_q_bb = collision_integral_delta_b(cloud_species%cloud, 1)

       rime_ice_cloud_theta_n_aa = collision_integral_theta_b(cloud_species%ice, 0)
       rime_ice_cloud_theta_n_ab = collision_integral_theta_ba(cloud_species%ice, cloud_species%cloud, 0)
       rime_ice_cloud_theta_n_bb = collision_integral_theta_b(cloud_species%cloud, 0)
       rime_ice_cloud_theta_q_aa = collision_integral_theta_b(cloud_species%ice, 1) ! mass weighted
       rime_ice_cloud_theta_q_ab = collision_integral_theta_ba(cloud_species%ice, cloud_species%cloud, 1)
       rime_ice_cloud_theta_q_ba = collision_integral_theta_ba(cloud_species%cloud, cloud_species%ice, 1)
       rime_ice_cloud_theta_q_bb = collision_integral_theta_b(cloud_species%cloud, 1)
!
!--    For ice - rain riming
       rime_ice_rain_delta_n_aa = collision_integral_delta_b(cloud_species%ice, 0)
       rime_ice_rain_delta_n_ab = collision_integral_delta_ba(cloud_species%ice, cloud_species%rain, 0)
       rime_ice_rain_delta_n_bb = collision_integral_delta_b(cloud_species%rain, 0)
       rime_ice_rain_delta_q_aa = collision_integral_delta_b(cloud_species%ice, 1) ! mass weighted
       rime_ice_rain_delta_q_ab = collision_integral_delta_ba(cloud_species%ice, cloud_species%rain, 1)
       rime_ice_rain_delta_q_ba = collision_integral_delta_ba(cloud_species%rain, cloud_species%ice, 1)
       rime_ice_rain_delta_q_bb = collision_integral_delta_b(cloud_species%rain, 1)
       rime_ice_rain_theta_n_aa = collision_integral_theta_b(cloud_species%ice, 0)
       rime_ice_rain_theta_n_ab = collision_integral_theta_ba(cloud_species%ice, cloud_species%rain, 0)
       rime_ice_rain_theta_n_bb = collision_integral_theta_b(cloud_species%rain, 0)
       rime_ice_rain_theta_q_aa = collision_integral_theta_b(cloud_species%ice, 1) ! mass weighted
       rime_ice_rain_theta_q_ab = collision_integral_theta_ba(cloud_species%ice, cloud_species%rain, 1)
       rime_ice_rain_theta_q_ba = collision_integral_theta_ba(cloud_species%rain, cloud_species%ice, 1)
       rime_ice_rain_theta_q_bb = collision_integral_theta_b(cloud_species%rain, 1)
!
!--    For snow - cloud riming
       rime_snow_cloud_delta_n_aa = collision_integral_delta_b(cloud_species%snow, 0)
       rime_snow_cloud_delta_n_ab = collision_integral_delta_ba(cloud_species%snow,cloud_species%cloud,0)
       rime_snow_cloud_delta_n_bb = collision_integral_delta_b(cloud_species%cloud,0)
       rime_snow_cloud_delta_q_aa = collision_integral_delta_b(cloud_species%snow, 0)
       rime_snow_cloud_delta_q_ab = collision_integral_delta_ba(cloud_species%snow,cloud_species%cloud,1)
       rime_snow_cloud_delta_q_bb = collision_integral_delta_b(cloud_species%cloud, 1)

       rime_snow_cloud_theta_n_aa = collision_integral_theta_b(cloud_species%snow, 0)
       rime_snow_cloud_theta_n_ab = collision_integral_theta_ba(cloud_species%snow,cloud_species%cloud,0)
       rime_snow_cloud_theta_n_bb = collision_integral_theta_b(cloud_species%cloud, 0)
       rime_snow_cloud_theta_q_aa = collision_integral_theta_b(cloud_species%snow,  0)
       rime_snow_cloud_theta_q_ab = collision_integral_theta_ba(cloud_species%snow,cloud_species%cloud,1)
       rime_snow_cloud_theta_q_bb = collision_integral_theta_b(cloud_species%cloud, 1)
!
!--    For snow - rain riming
       rime_snow_rain_delta_n_aa = collision_integral_delta_b(cloud_species%snow, 0)
       rime_snow_rain_delta_n_ab = collision_integral_delta_ba(cloud_species%snow, cloud_species%rain, 0)
       rime_snow_rain_delta_n_bb = collision_integral_delta_b(cloud_species%rain, 0)
       rime_snow_rain_delta_q_aa = collision_integral_delta_b(cloud_species%snow, 1) ! mass weighted
       rime_snow_rain_delta_q_ab = collision_integral_delta_ba(cloud_species%snow, cloud_species%rain, 1)
       rime_snow_rain_delta_q_ba = collision_integral_delta_ba(cloud_species%rain, cloud_species%snow, 1)
       rime_snow_rain_delta_q_bb = collision_integral_delta_b(cloud_species%rain, 1)

       rime_snow_rain_theta_n_aa = collision_integral_theta_b(cloud_species%snow, 0)
       rime_snow_rain_theta_n_ab = collision_integral_theta_ba(cloud_species%snow, cloud_species%rain, 0)
       rime_snow_rain_theta_n_bb = collision_integral_theta_b(cloud_species%rain, 0)
       rime_snow_rain_theta_q_aa = collision_integral_theta_b(cloud_species%snow, 1) ! mass weighted
       rime_snow_rain_theta_q_ab = collision_integral_theta_ba(cloud_species%snow, cloud_species%rain, 1)
       rime_snow_rain_theta_q_ba = collision_integral_theta_ba(cloud_species%rain, cloud_species%snow, 1)
       rime_snow_rain_theta_q_bb = collision_integral_theta_b(cloud_species%rain, 1)

    END SUBROUTINE init_collision_integrals

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize coefficients such as for ventilation and other processes
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE init_coefficients

!
!--    Calculate ventilation coefficients for snow calling function
       a1_ven_coeff_snow = ventilation_coeff_a( cloud_species%snow, 1)
       b1_ven_coeff_snow = ventilation_coeff_b( cloud_species%snow, 1) *                           &
                           schmidt_p_1d3 / SQRT(kin_vis_air)
!
!--    Calculate ventilation coefficients for graupel calling function
       a1_ven_coeff_graupel = ventilation_coeff_a( cloud_species%graupel, 1)
       b1_ven_coeff_graupel = ventilation_coeff_b( cloud_species%graupel, 1) *                     &
                           schmidt_p_1d3 / SQRT(kin_vis_air)
!
!--    calculate coefficient for homogeneous freezing proportional to 2nd moment of distribution
       cz_cloud = moment_gamma( cloud_species%cloud, 2 )

    END SUBROUTINE init_coefficients


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust number of raindrops to avoid nonlinear effects in sedimentation and evaporation of rain
!> drops due to too small or too big weights of rain drops (Stevens and Seifert, 2008).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE adjust_cloud

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp) ::  flag                  !< flag to indicate first grid level above

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                IF ( .NOT. microphysics_morrison_no_rain )  THEN
                   IF ( qr(k,j,i) <= eps_sb )  THEN
                      qr(k,j,i) = 0.0_wp
                      nr(k,j,i) = 0.0_wp
                   ELSE
                      IF ( nr(k,j,i) * xrmin > qr(k,j,i) * hyrho(k) )  THEN
                         nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmin * flag
                      ELSEIF ( nr(k,j,i) * xrmax < qr(k,j,i) * hyrho(k) )  THEN
                         nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmax * flag
                      ENDIF
                   ENDIF
                ENDIF

                IF ( microphysics_morrison )  THEN
                   IF ( qc(k,j,i) <= eps_sb )  THEN
                      qc(k,j,i) = 0.0_wp
                      nc(k,j,i) = 0.0_wp
                   ELSE
                      IF ( nc(k,j,i) * xcmin > qc(k,j,i) * hyrho(k) )  THEN
                         nc(k,j,i) = qc(k,j,i) * hyrho(k) / xcmin * flag
                      ENDIF
                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE adjust_cloud

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust number of raindrops to avoid nonlinear effects in sedimentation and evaporation of rain
!> drops due to too small or too big weights of rain drops (Stevens and Seifert, 2008).
!> The same procedure is applied to cloud droplets if they are determined prognostically. Call for
!> grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE adjust_cloud_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp) ::  flag                  !< flag to indicate first grid level above surface

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          IF ( .NOT. microphysics_morrison_no_rain )  THEN
             IF ( qr(k,j,i) <= eps_sb )  THEN
                qr(k,j,i) = 0.0_wp
                nr(k,j,i) = 0.0_wp
             ELSE
!
!--             Adjust number of raindrops to avoid nonlinear effects in sedimentation and
!--             evaporation of rain drops due to too small or too big weights of rain drops
!--            (Stevens and Seifert, 2008).
                IF ( nr(k,j,i) * xrmin > qr(k,j,i) * hyrho(k) )  THEN
                   nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmin * flag
                ELSEIF ( nr(k,j,i) * xrmax < qr(k,j,i) * hyrho(k) )  THEN
                   nr(k,j,i) = qr(k,j,i) * hyrho(k) / xrmax * flag
                ENDIF
             ENDIF
          ENDIF

          IF ( microphysics_morrison )  THEN
             IF ( qc(k,j,i) <= eps_sb )  THEN
                qc(k,j,i) = 0.0_wp
                nc(k,j,i) = 0.0_wp
             ELSE
                IF ( nc(k,j,i) * xcmin > qc(k,j,i) * hyrho(k) )  THEN
                   nc(k,j,i) = qc(k,j,i) * hyrho(k) / xcmin * flag
                ENDIF
             ENDIF
          ENDIF

       ENDDO

    END SUBROUTINE adjust_cloud_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust number of ice crystal to avoid nonlinear effects in sedimentation and evaporation of ice
!> crytals due to too small or too big weights of ice crytals (Stevens and Seifert, 2008).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE adjust_ice

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp) ::  flag                  !< flag to indicate first grid level above

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                IF ( qi(k,j,i) <= ximin )  THEN
                   qi(k,j,i) = 0.0_wp
                   ni(k,j,i) = 0.0_wp
                ELSE
                   IF ( ni(k,j,i) * ximin > qi(k,j,i) * hyrho(k) )  THEN
                      ni(k,j,i) = qi(k,j,i) * hyrho(k) / ximin * flag
                   ELSEIF ( ni(k,j,i) * ximax < qi(k,j,i) * hyrho(k) )  THEN
                      ni(k,j,i) = qi(k,j,i) * hyrho(k) / ximax * flag
                   ENDIF
                ENDIF
                IF ( snow  .OR.  graupel )  THEN
!
!--                snow adjust
                   IF ( qs(k,j,i) <= xsmin )  THEN
                      qs(k,j,i) = 0.0_wp
                      ns(k,j,i) = 0.0_wp
                   ELSE
                      IF ( ns(k,j,i) * xsmin > qs(k,j,i) * hyrho(k) )  THEN
                         ns(k,j,i) = qs(k,j,i) * hyrho(k) / xsmin * flag
                      ELSEIF ( ns(k,j,i) * xsmax < qs(k,j,i) * hyrho(k) )  THEN
                         ns(k,j,i) = qs(k,j,i) * hyrho(k) / xsmax * flag
                      ENDIF
                   ENDIF
!
!--                graupel adjust
                   IF ( qg(k,j,i) <= eps_sb )  THEN
                      qg(k,j,i) = 0.0_wp
                      ng(k,j,i) = 0.0_wp
                   ELSE
                      IF ( ng(k,j,i) * xgmin > qg(k,j,i) * hyrho(k) )  THEN
                         ng(k,j,i) = qg(k,j,i) * hyrho(k) / xgmin * flag
                      ELSEIF ( ng(k,j,i) * xgmax < qg(k,j,i) * hyrho(k) )  THEN
                         ng(k,j,i) = qg(k,j,i) * hyrho(k) / xgmax * flag
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE adjust_ice

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust number of of ice crystal to avoid nonlinear effects in sedimentation and evaporation of
!> ice crystals due to too small or too big weights of ice crytals (Stevens and Seifert, 2008).
!> The same procedure is applied to cloud droplets if they are determined prognostically. Call for
!> grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE adjust_ice_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp) ::  flag                  !< flag to indicate first grid level above surface

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       ice adjust
          IF ( qi(k,j,i) <= ximin )  THEN
             qi(k,j,i) = 0.0_wp
             ni(k,j,i) = 0.0_wp
          ELSE
             IF ( ni(k,j,i) * ximin > qi(k,j,i) * hyrho(k) )  THEN
                ni(k,j,i) = qi(k,j,i) * hyrho(k) / ximin * flag
             ELSEIF ( ni(k,j,i) * ximax < qi(k,j,i) * hyrho(k) )  THEN
                ni(k,j,i) = qi(k,j,i) * hyrho(k) / ximax * flag
             ENDIF
          ENDIF
          IF ( snow  .OR.  graupel )  THEN
!
!--          snow adjust
             IF ( qs(k,j,i) <= xsmin )  THEN
                qs(k,j,i) = 0.0_wp
                ns(k,j,i) = 0.0_wp
             ELSE
                IF ( ns(k,j,i) * xsmin > qs(k,j,i) * hyrho(k) )  THEN
                   ns(k,j,i) = qs(k,j,i) * hyrho(k) / xsmin * flag
                ELSEIF ( ns(k,j,i) * xsmax < qs(k,j,i) * hyrho(k) )  THEN
                   ns(k,j,i) = qs(k,j,i) * hyrho(k) / xsmax * flag
                ENDIF
             ENDIF
!
!--          graupel adjust
             IF ( qg(k,j,i) <= eps_sb )  THEN
                qg(k,j,i) = 0.0_wp
                ng(k,j,i) = 0.0_wp
             ELSE
                IF ( ng(k,j,i) * xgmin > qg(k,j,i) * hyrho(k) )  THEN
                   ng(k,j,i) = qg(k,j,i) * hyrho(k) / xgmin * flag
                ELSEIF ( ng(k,j,i) * xgmax < qg(k,j,i) * hyrho(k) )  THEN
                   ng(k,j,i) = qg(k,j,i) * hyrho(k) / xgmax * flag
                ENDIF
             ENDIF
          ENDIF
       ENDDO

    END SUBROUTINE adjust_ice_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate number of activated condensation nucleii after simple activation
!> scheme of Twomey, 1959.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE activation_cloud

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  activ             !<
       REAL(wp)     ::  afactor           !<
       REAL(wp)     ::  beta_act          !<
       REAL(wp)     ::  bfactor           !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above
       REAL(wp)     ::  k_act             !<
       REAL(wp)     ::  n_act             !<
       REAL(wp)     ::  n_ccn             !<
       REAL(wp)     ::  s_0               !<
       REAL(wp)     ::  sat_max           !<
       REAL(wp)     ::  sigma             !<
       REAL(wp)     ::  sigma_act         !<

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Call calculation of supersaturation located in subroutine
                CALL supersaturation ( i, j, k )
!
!--             Prescribe parameters for activation (see: Bott + Trautmann, 2002, Atm. Res., 64)
                k_act  = 0.7_wp
                activ  = 0.0_wp


                IF ( sat > 0.0  .AND.  .NOT. curvature_solution_effects_bulk )  THEN
!
!--                Compute the number of activated Aerosols
!--                (see: Twomey, 1959, Pure and applied Geophysics, 43)
                   n_act     = na_init * sat**k_act
!
!--                Compute the number of cloud droplets (see: Morrison + Grabowski, 2007, JAS, 64)
!                  activ = MAX( n_act - nc(k,j,i), 0.0_wp) / dt_micro

!
!--                Compute activation rate after Khairoutdinov and Kogan
!--                (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev., 128)
                   sat_max = 1.0_wp / 100.0_wp
                   activ   = MAX( 0.0_wp, ( (na_init + nc(k,j,i) ) *                               &
                             MIN( 1.0_wp, ( sat / sat_max )**k_act) - nc(k,j,i) ) ) / dt_micro
!
!--                Apply activation
                   nc(k,j,i) = MIN( (nc(k,j,i) + activ * dt_micro * flag), na_init)
!
!--             If curvature and solution effects are considered
                ELSEIF ( sat > 0.0  .AND.  curvature_solution_effects_bulk )  THEN
!
!--                Curvature effect (afactor) with surface tension parameterization by Straka (2009)
                   sigma = 0.0761_wp - 0.000155_wp * ( t_l - 273.15_wp )
                   afactor = 2.0_wp * sigma / ( rho_l * r_v * t_l )
!
!--                Solute effect (bfactor)
                   bfactor = vanthoff * molecular_weight_of_water * rho_s /                        &
                             ( molecular_weight_of_solute * rho_l )

!
!--                Prescribe power index that describes the soluble fraction of an aerosol particle
!--                (beta) (see: Morrison + Grabowski, 2007, JAS, 64)
                   beta_act  = 0.5_wp
                   sigma_act = sigma_bulk**( 1.0_wp + beta_act )
!
!--                Calculate mean geometric supersaturation (s_0) with parameterization by
!--                Khvorostyanov and Curry (2006)
                   s_0 = dry_aerosol_radius **(- ( 1.0_wp + beta_act ) ) * ( 4.0_wp * afactor**3 / &
                         ( 27.0_wp * bfactor ) )**0.5_wp

!
!--                Calculate number of activated CCN as a function of supersaturation and taking
!--                Koehler theory into account (see: Khvorostyanov + Curry, 2006, J. Geo. Res., 111)
                   n_ccn = ( na_init / 2.0_wp ) * ( 1.0_wp - ERF(                                  &
                           LOG( s_0 / sat ) / ( SQRT(2.0_wp) * LOG(sigma_act) ) ) )
                   activ = MAX( ( n_ccn - nc(k,j,i) ) / dt_micro, 0.0_wp )
!
!--                Apply activation
                   nc(k,j,i) = MIN( (nc(k,j,i) + activ * dt_micro * flag), na_init)
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE activation_cloud

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate number of activated condensation nucleii after simple activation scheme of
!> Twomey, 1959.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE activation_cloud_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !<
       INTEGER(iwp) ::  j                 !<
       INTEGER(iwp) ::  k                 !<

       REAL(wp)     ::  activ             !<
       REAL(wp)     ::  afactor           !<
       REAL(wp)     ::  beta_act          !<
       REAL(wp)     ::  bfactor           !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  k_act             !<
       REAL(wp)     ::  n_act             !<
       REAL(wp)     ::  n_ccn             !<
       REAL(wp)     ::  s_0               !<
       REAL(wp)     ::  sat_max           !<
       REAL(wp)     ::  sigma             !<
       REAL(wp)     ::  sigma_act         !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Call calculation of supersaturation
          CALL supersaturation ( i, j, k )
!
!--       Prescribe parameters for activation (see: Bott + Trautmann, 2002, Atm. Res., 64)
          k_act  = 0.7_wp
          activ  = 0.0_wp

          IF ( sat > 0.0 .AND. .NOT. curvature_solution_effects_bulk )  THEN
!
!--          Compute the number of activated Aerosols
!--          (see: Twomey, 1959, Pure and applied Geophysics, 43)
             n_act     = na_init * sat**k_act
!
!--          Compute activation rate after Khairoutdinov and Kogan
!--          (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev., 128)
             sat_max = 0.8_wp / 100.0_wp
             activ   = MAX( 0.0_wp, ( (na_init + nc(k,j,i) ) *                                     &
                       MIN( 1.0_wp, ( sat / sat_max )**k_act) - nc(k,j,i) ) ) / dt_micro

             nc(k,j,i) = MIN( (nc(k,j,i) + activ * dt_micro), na_init)
          ELSEIF ( sat > 0.0  .AND.  curvature_solution_effects_bulk )  THEN
!
!--          Curvature effect (afactor) with surface tension parameterization by Straka (2009)
             sigma = 0.0761_wp - 0.000155_wp * ( t_l - 273.15_wp )
             afactor = 2.0_wp * sigma / ( rho_l * r_v * t_l )
!
!--          Solute effect (bfactor)
             bfactor = vanthoff * molecular_weight_of_water * rho_s /                              &
                       ( molecular_weight_of_solute * rho_l )

!
!--          Prescribe power index that describes the soluble fraction of an aerosol particle
!--          (beta). (see: Morrison + Grabowski, 2007, JAS, 64)
             beta_act  = 0.5_wp
             sigma_act = sigma_bulk**( 1.0_wp + beta_act )
!
!--          Calculate mean geometric supersaturation (s_0) with parameterization by Khvorostyanov
!--          and Curry (2006)
             s_0   = dry_aerosol_radius **(- ( 1.0_wp + beta_act ) ) *                             &
                     ( 4.0_wp * afactor**3 / ( 27.0_wp * bfactor ) )**0.5_wp

!
!--          Calculate number of activated CCN as a function of supersaturation and taking Koehler
!--          theory into account (see: Khvorostyanov + Curry, 2006, J. Geo. Res., 111)
             n_ccn = ( na_init / 2.0_wp ) * ( 1.0_wp - ERF(                                        &
                     LOG( s_0 / sat ) / ( SQRT(2.0_wp) * LOG(sigma_act) ) ) )
             activ = MAX( ( n_ccn ) / dt_micro, 0.0_wp )

             nc(k,j,i) = MIN( (nc(k,j,i) + activ * dt_micro * flag), na_init)
          ENDIF

       ENDDO

    END SUBROUTINE activation_cloud_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate ice nucleation by applying the deposition-condensation formula as given by
!> Meyers et al 1992 and as described in Seifert and Beheng 2006
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE nucleation_ice

       INTEGER(iwp) ::  i             !< loop index
       INTEGER(iwp) ::  j             !< loop index
       INTEGER(iwp) ::  k             !< loop index

       LOGICAL :: isdac = .FALSE.

       REAL(wp) ::  a_m92 = -0.639_wp !< parameter for nucleation
       REAL(wp) ::  b_m92 = 12.96_wp  !< parameter for nucleation
       REAL(wp) ::  flag              !< flag to indicate first grid level
       REAL(wp) ::  n_in              !< number of ice nucleii
       REAL(wp) ::  nucle = 0.0_wp    !< nucleation rate

       IF  ( isdac )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
!
!--                Predetermine flag to mask topography
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--                Call calculation of supersaturation located in subroutine
                   CALL supersaturation_ice ( i, j, k )
                   nucle = 0.0_wp
                   IF ( sat_ice >= 0.05_wp  .OR.  ql(k,j,i) >= 0.001E-3_wp  )  THEN
!
!--                   Calculate ice nucleation
                      nucle = MAX( ( in_init - ni(k,j,i) ) / dt_micro, 0.0_wp )
                      ni(k,j,i) = MIN( (ni(k,j,i) + nucle * dt_micro * flag), in_init)
                   ENDIF

                ENDDO
             ENDDO
          ENDDO
       ELSE
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
!
!--                Predetermine flag to mask topography
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--                Call calculation of supersaturation located in subroutine
                   CALL supersaturation_ice ( i, j, k )
                   nucle = 0.0_wp
                   IF ( sat_ice > 0.0_wp  .AND.  sat_ice < 2.0_wp ) THEN
!
!--                   Calculate ice nucleation
                      n_in = in_init * EXP( a_m92 + b_m92 * sat_ice )
                      nucle = MAX( ( n_in - ni(k,j,i) ) / dt_micro, 0.0_wp )
                   ENDIF
                   ni(k,j,i) = MIN( (ni(k,j,i) + nucle * dt_micro * flag), in_init )
                ENDDO
             ENDDO
          ENDDO
       ENDIF

    END SUBROUTINE nucleation_ice


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate ice nucleation by applying the deposition-condensation formula as given by
!> Meyers et al 1992 and as described in Seifert and Beheng 2006.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE nucleation_ice_ij( i, j )

       INTEGER(iwp) ::  i             !< loop index
       INTEGER(iwp) ::  j             !< loop index
       INTEGER(iwp) ::  k             !< loop index

       LOGICAL :: isdac = .FALSE.

       REAL(wp) ::  a_m92 = -0.639_wp !< parameter for nucleation
       REAL(wp) ::  b_m92 = 12.96_wp  !< parameter for nucleation
       REAL(wp) ::  flag              !< flag to indicate first grid level
       REAL(wp) ::  n_in              !< number of ice nucleii
       REAL(wp) ::  nucle = 0.0_wp    !< nucleation rate

       IF ( isdac )  THEN
          DO  k = nzb+1, nzt
!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--          Call calculation of supersaturation located in subroutine
             CALL supersaturation_ice ( i, j, k )
             nucle = 0.0_wp
             IF ( sat_ice >= 0.05_wp  .OR.  ql(k,j,i) >= 0.001E-3_wp )  THEN
!
!--             Calculate ice nucleation
                nucle = MAX( ( in_init - ni(k,j,i) ) / dt_micro, 0.0_wp )
                ni(k,j,i) = MIN( (ni(k,j,i) + nucle * dt_micro * flag), in_init)
             ENDIF
          ENDDO
       ELSE
          DO  k = nzb+1, nzt
!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--          Call calculation of supersaturation located in subroutine
             CALL supersaturation_ice ( i, j, k )
             nucle = 0.0_wp
             IF ( sat_ice > 0.0_wp  .AND.  sat_ice < 2.0_wp )  THEN
!
!--             Calculate ice nucleation
                n_in = in_init * EXP( a_m92 + b_m92 * sat_ice )
                nucle = MAX( ( n_in - ni(k,j,i) ) / dt_micro, 0.0_wp )
             ENDIF
             ni(k,j,i) = MIN( (ni(k,j,i) + nucle * dt_micro * flag), in_init )
          ENDDO
       ENDIF


    END SUBROUTINE nucleation_ice_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate condensation/evaporation rate for cloud water (after Khairoutdinov and Kogan, 2000),
!> assuming a gamma size distribution for droplets, with shape nu=1. This is also applied in
!> Seifert and Stevens, 2008.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE condensation_cloud

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  alpha_rc = 1.0_wp !< Tuning parameter see (Seifert and Stevens, 2010)
       REAL(wp)     ::  cond              !< condensation rate
       REAL(wp)     ::  cond_max          !< maximum condensation rate
       REAL(wp)     ::  evap              !< evaporation rate
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  g_fac             !< factor 1 / Fk + Fd
       REAL(wp)     ::  nu = 1.0_wp       !< Shape parameter of gernerlized gama distribution
       REAL(wp)     ::  rc                !< mean cloud droplets radius assuming gamma distribution
       REAL(wp)     ::  temp              !< actual temperature

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Call calculation of supersaturation
                CALL supersaturation ( i, j, k )
!
!--             Actual temperature, t_l is calculated directly before in supersaturation
                IF ( microphysics_ice_phase ) THEN
                   temp = t_l + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
                ELSE
                   temp = t_l + lv_d_cp * ql(k,j,i)
                ENDIF

                g_fac  = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *                            &
                                      l_v / ( thermal_conductivity_l  * temp )                     &
                                    + r_v * temp / ( diff_coeff_l * e_s )                          &
                                  )
!
!--             Mean weight of cloud drops
                IF ( nc(k,j,i) <= 0.0_wp ) CYCLE
!
!--             Calculating mean radius of cloud droplets assuming gamma distribution with shape
!--             parameter nu=1 (Seifert and Beheng, 2006). Tuning factor alpha_rc (introduced
!--             in Seifert and Stevens, 2010 ) is switched off. Minimum radius is set to 1Âµm
!--             following (Seifert and Beheng, 2006, Kogan and Khairoutdinov, 2000,
!--             Seifert and Stevens, 2010)
                rc = MAX( alpha_rc  * GAMMA( nu + 1.33_wp ) / GAMMA( nu + 1.0_wp ) *               &
                          ( 3.0_wp * qc(k,j,i) /                                                   &
                                   ( 4.0_wp * pi * rho_l  * ( nu + 2.0_wp ) * nc(k,j,i) )          &
                          )**0.33_wp, 1.0E-6_wp )
!
!--             Condensation needs only to be calculated in supersaturated regions
                IF ( sat > 0.0_wp )  THEN
!
!--                Condensation rate of cloud water content after KK scheme.
!--                (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev.,128, Morrison and Grabowski,
!--                2007, and Seifert and Stevens, 2010)
                   cond      = 4.0_wp * pi * nc(k,j,i) * g_fac * sat * rc / hyrho(k)
                   IF ( microphysics_seifert )  THEN
                      cond_max  = q(k,j,i) - q_s - qc(k,j,i) - qr(k,j,i)
                   ELSEIF ( microphysics_morrison_no_rain )  THEN
                      cond_max  = q(k,j,i) - q_s - qc(k,j,i)
                   ENDIF
                   cond      = MIN( cond, cond_max / dt_micro )

                   qc(k,j,i) = qc(k,j,i) + cond * dt_micro * flag
                ELSEIF ( sat < 0.0_wp )  THEN
                   evap      = 4.0_wp * pi * nc(k,j,i) * g_fac * sat * rc / hyrho(k)
                   evap      = MAX( evap, -qc(k,j,i) / dt_micro )

                   qc(k,j,i) = qc(k,j,i) + evap * dt_micro * flag
                ENDIF
                IF ( nc(k,j,i) * xcmin > qc(k,j,i) * hyrho(k) )  THEN
                   nc(k,j,i) = qc(k,j,i) * hyrho(k) / xcmin
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE condensation_cloud

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate condensation/evaporation rate for cloud water (after Khairoutdinov and Kogan, 2000),
!> assuming a gamma size distribution for droplets, with shape nu=1. This is also applied in
!> Seifert and Stevens, 2008.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE condensation_cloud_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  alpha_rc = 1.0_wp !< Tuning parameter see (Seifert and Stevens, 2010)
       REAL(wp)     ::  cond              !< condensation rate
       REAL(wp)     ::  cond_max          !< maximum condensation rate
       REAL(wp)     ::  evap              !< evaporation rate
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  g_fac             !< factor 1 / Fk + Fd
       REAL(wp)     ::  nu = 1.0_wp       !< Shape parameter of gernerlized gama distribution
       REAL(wp)     ::  rc                !< mean cloud droplets radius assuming gamma distribution
       REAL(wp)     ::  temp              !< actual temperature

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Call calculation of supersaturation
          CALL supersaturation ( i, j, k )
!
!--       Actual temperature, t_l is calculated directly before in supersaturation
          IF ( microphysics_ice_phase ) THEN
             temp = t_l + lv_d_cp * ql(k,j,i) + ls_d_cp * qi(k,j,i)
          ELSE
             temp = t_l + lv_d_cp * ql(k,j,i)
          ENDIF

          g_fac  = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *                                  &
                                l_v / ( thermal_conductivity_l  * temp )                           &
                              + r_v * temp / ( diff_coeff_l * e_s )                                &
                            )
!
!--       Mean weight of cloud drops
          IF ( nc(k,j,i) <= 0.0_wp)  CYCLE
!
!--       Calculating mean radius of cloud droplets assuming gamma distribution with shape
!--       parameter nu=1 (Seifert and Beheng, 2006). Tuning factor alpha_rc (introduced
!--       in Seifert and Stevens, 2010 ) is switched off. Minimum radius is set to 1Âµm following
!--       (Seifert and Beheng, 2006, Kogan and Khairoutdinov, 2000, Seifert and Stevens, 2010)
          rc = MAX( alpha_rc  * GAMMA( nu + 1.33_wp ) / GAMMA( nu + 1.0_wp ) *                     &
                    ( 3.0_wp * qc(k,j,i) / ( 4.0_wp * pi * rho_l  * ( nu + 2.0_wp ) * nc(k,j,i) )  &
                  )**0.33_wp, 1.0E-6_wp )
!
!--       Condensation needs only to be calculated in supersaturated regions
          IF ( sat > 0.0_wp )  THEN
!
!--          Condensation rate of cloud water content after KK scheme.
!--          (see: Khairoutdinov + Kogan, 2000, Mon. Wea. Rev.,128, Morrison and Grabowski, 2007
!--          and Seifert and Stevens, 2010)
             cond      = 4.0_wp * pi * nc(k,j,i) * g_fac * sat * rc / hyrho(k)
             IF ( microphysics_seifert )  THEN
                cond_max  = q(k,j,i) - q_s - qc(k,j,i) - qr(k,j,i)
             ELSEIF ( microphysics_morrison_no_rain )  THEN
                cond_max  = q(k,j,i) - q_s - qc(k,j,i)
             ENDIF
             cond      = MIN( cond, cond_max / dt_micro )

             qc(k,j,i) = qc(k,j,i) + cond * dt_micro * flag
          ELSEIF ( sat < 0.0_wp )  THEN
             evap      = 4.0_wp * pi * nc(k,j,i) * g_fac * sat * rc / hyrho(k)
             evap      = MAX( evap, -qc(k,j,i) / dt_micro )

             qc(k,j,i) = qc(k,j,i) + evap * dt_micro * flag
          ENDIF
       ENDDO

    END SUBROUTINE condensation_cloud_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the growth/shrink of ice particles by water vapor deposition/sublimation
!> (after Seifert and Beheng, 2006) assuming gamam size distribution. Right now ventilation is
!> neglected.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE deposition_ice

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  ac = 0.09_wp        !< parameter for ice capacitance
       REAL(wp) ::  bc = 0.33_wp        !< parameter for ice capacitance
       REAL(wp) ::  deposition_rate = 0.0_wp !< depositions rate
       REAL(wp) ::  deposition_rate_max !< maximum deposition rate
       REAL(wp) ::  fac_gamma = 0.76_wp !< parameter to describe spectral distribution, here
                                        !< following gamma size distribution with µ =1/3 and nu=0
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  gfac_dep            !< factor
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  xi                  !< mean mass of ice crystal

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Call calculation of supersaturation over a plane ice surface
                CALL supersaturation_ice ( i, j, k )
!
!--             Actual temperature:
                temp = t_l + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--             calculating gfac_dep ( 1/ (Fk + Fd) ) see e.g. Rogers and Yau, 1989
                gfac_dep  = 1.0_wp / ( ( l_s / ( r_v * temp ) - 1.0_wp ) *                         &
                                         l_s / ( thermal_conductivity_l  * temp )                  &
                                    + r_v * temp / ( diff_coeff_l * e_si )                         &
                                     )
!
!--             If there is nothing nucleated, than there is also no deposition (above -38°C)
                IF ( ni(k,j,i) <= 0.0_wp )  CYCLE
!
!--             calculate mean mass of ice crystal
                xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
!
!--             Condensation needs only to be calculated in supersaturated regions (regarding ice)
                IF ( sat_ice > 0.0_wp )  THEN
!
!--                Calculate deposition rate assuming ice crystal shape as prescribed in
!--                Ovchinnikov et al., 2014 and a gamma size distribution according to
!--                Seifert and Beheng with µ =1/3 and nu=0
                   deposition_rate  = 4.0_wp * pi * sat_ice * gfac_dep * fac_gamma * ac * xi**bc * &
                                      ni(k,j,i)
                   IF ( microphysics_seifert )  THEN
                      deposition_rate_max  = q(k,j,i) - q_si - qr(k,j,i) - qi(k,j,i)
                   ELSEIF ( microphysics_morrison_no_rain ) THEN
                      deposition_rate_max  = q(k,j,i) - q_si - qc(k,j,i) - qi(k,j,i)
                   ENDIF
                   deposition_rate = MIN( deposition_rate, deposition_rate_max / dt_micro )
                   qi(k,j,i) = qi(k,j,i) + deposition_rate * dt_micro * flag
                ELSEIF ( sat_ice < 0.0_wp )  THEN
                   deposition_rate = 4.0_wp * pi * sat_ice * gfac_dep * fac_gamma * ac * xi**bc * &
                                      ni(k,j,i)
                   deposition_rate = MAX( deposition_rate, - qi(k,j,i) / dt_micro )
                   qi(k,j,i) = qi(k,j,i) + deposition_rate * dt_micro * flag
                ENDIF
!
!--             Store deposition rate ice for riming process
                dep_rate_ice(k,j,i) = deposition_rate
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE deposition_ice


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the growth/shrink of ice particles by water vapor deposition/sublimation
!> (after Seifert and Beheng, 2006) assuming gamam size distribution. Right now ventilation is
!> neglected.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE deposition_ice_ij( i, j )

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  ac = 0.09_wp        !< parameter for ice capacitance
       REAL(wp) ::  bc = 0.33_wp        !< parameter for ice capacitance
       REAL(wp) ::  deposition_rate = 0.0_wp !< depositions rate
       REAL(wp) ::  deposition_rate_max !< maximum deposition rate
       REAL(wp) ::  fac_gamma = 0.76_wp !< parameter to describe spectral distribution, here
                                        !< following gamma size distribution with nu=1, v=0
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  gfac_dep            !< factor
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  xi                  !< mean mass of ice crystal

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Call calculation of supersaturation over a plane ice surface
          CALL supersaturation_ice( i, j, k )
!
!--       Actual temperature, t_l is calcualted in supersaturation_ice routine:
          temp = t_l + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--       Calculating gfac_dep ( 1/ (Fk + Fd) ) see e.g.
!--       Rogers and Yau, 1989
          gfac_dep  = 1.0_wp / ( ( l_s / ( r_v * temp ) - 1.0_wp ) *                               &
                                   l_s / ( thermal_conductivity_l  * temp )                        &
                                 + r_v * temp / ( diff_coeff_l * e_si )                            &
                               )
!
!--       If there is nothing nucleated, than there is also no deposition (above -38°C)
          IF ( ni(k,j,i) > 0.0_wp )  THEN
!
!--          calculate mean mass of ice crystal
             xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
!
!--          Water vapor deposition in supersaturated regions (regarding ice)
             IF ( sat_ice > 0.0_wp )  THEN
!
!--             Calculate deposition rate assuming ice crystal shape as prescribed in
!--             Ovchinnikov et al., 2014 and a gamma size distribution according to
!--             Seifert and Beheng with µ =1/3 and nu=0
                deposition_rate  = 4.0_wp * pi * sat_ice * gfac_dep * fac_gamma * ac * xi**bc *    &
                                   ni(k,j,i)
                deposition_rate_max  = q(k,j,i) - q_si - qc(k,j,i) - qr(k,j,i) - qf(k,j,i)
                deposition_rate = MIN( deposition_rate, deposition_rate_max / dt_micro )

                qi(k,j,i) = qi(k,j,i) + deposition_rate * dt_micro * flag
!
!--          Water vapor sublimation in supersaturated regions (regarding ice)
             ELSEIF ( sat_ice < 0.0_wp )  THEN
                deposition_rate = 4.0_wp * pi * sat_ice * gfac_dep * fac_gamma * ac * xi**bc *    &
                                   ni(k,j,i)
                deposition_rate = MAX( deposition_rate, - qi(k,j,i) / dt_micro )
                qi(k,j,i) = qi(k,j,i) + deposition_rate * dt_micro * flag
             ENDIF
          ENDIF
!
!--       Store deposition rate ice for riming process
          dep_rate_ice(k,j,i) = deposition_rate
       ENDDO

    END SUBROUTINE deposition_ice_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the growth/shrink of snowflakes by water vapor deposition/sublimation
!> (after Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE deposition_snow

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       LOGICAL :: ventilation_snow = .TRUE.

       REAL(wp) ::  ds                  !< mean diameter of snow
       REAL(wp) ::  deposition_rate     !< depositions rate
       REAL(wp) ::  deposition_rate_max !< maximum deposition rate
       REAL(wp) ::  fac_gamma = 1.0_wp  !< parameter to describe spectral distribution,
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  f_vent              !< ventilation coefficient
       REAL(wp) ::  gfac_dep            !< factor
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vs                  !< mean fall velocity of snow
       REAL(wp) ::  xs                  !< mean mass of snow

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Call calculation of supersaturation over a plane ice surface
                CALL supersaturation_ice( i, j, k )
!
!--             Actual temperature, t_l is calcualted in supersaturation_ice routine:
                temp = t_l + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--             Calculating gfac_dep ( 1/ (Fk + Fd) ) see e.g.
!--             Rogers and Yau, 1989
                gfac_dep  = 1.0_wp / ( ( l_s / ( r_v * temp ) - 1.0_wp ) *                         &
                                         l_s / ( thermal_conductivity_l  * temp )                  &
                                       + r_v * temp / ( diff_coeff_l * e_si )                      &
                                     )
!
!--             If there is nothing nucleated, than there is also no deposition (above -38°C)
                IF ( ns(k,j,i) > 0.0_wp )  THEN
!
!--                calculate mean mass of snow crystal
                   xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )

!
!--                Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and
!--                mean diameter of snow
                   vs = terminal_fall_velocity( cloud_species%snow, xs, hyrho(k) )
                   ds = mean_diameter( cloud_species%snow, xs )

                   IF ( ventilation_snow )  THEN
                      f_vent  = a1_ven_coeff_snow + b1_ven_coeff_snow * SQRT( ds * vs )
                      f_vent  = MAX( f_vent, a1_ven_coeff_snow / cloud_species%snow%a_ven )
                   ELSE
                      f_vent = 1.0_wp
                   ENDIF
!
!--                Water vapor deposition in supersaturated regions (regarding ice)
                   IF ( sat_ice > 0.0_wp )  THEN
!
!--                   Calculate deposition rate assuming snow crystal shape
                      deposition_rate  = 4.0_wp * pi / cloud_species%snow%cap * sat_ice *          &
                                         gfac_dep * fac_gamma * ds * f_vent * ns(k,j,i)
                      deposition_rate_max  = q(k,j,i) - q_si - qc(k,j,i) - qr(k,j,i) - qf(k,j,i)
                      deposition_rate = MIN( deposition_rate, deposition_rate_max / dt_micro )

                      qs(k,j,i) = qs(k,j,i) + deposition_rate * dt_micro * flag
!
!--                Water vapor sublimation in supersaturated regions (regarding ice)
                   ELSEIF ( sat_ice < 0.0_wp )  THEN
                      deposition_rate = 4.0_wp * pi / cloud_species%snow%cap * sat_ice *           &
                                        gfac_dep * fac_gamma * ds * f_vent * ns(k,j,i)
                      deposition_rate = MAX( deposition_rate, - qs(k,j,i) / dt_micro )
                      qs(k,j,i) = qs(k,j,i) + deposition_rate * dt_micro * flag
                   ENDIF
!
!--                Store deposition rate ice for riming process
                   dep_rate_snow(k,j,i) = deposition_rate
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE deposition_snow


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the growth/shrink of snowflakes by water vapor deposition/sublimation
!> (after Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE deposition_snow_ij( i, j )

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       LOGICAL :: ventilation_snow = .TRUE.

       REAL(wp) ::  ds                  !< mean diameter of snow
       REAL(wp) ::  deposition_rate     !< depositions rate
       REAL(wp) ::  deposition_rate_max !< maximum deposition rate
       REAL(wp) ::  fac_gamma = 1.0_wp  !< parameter to describe spectral distribution,
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  f_vent              !< ventilation coefficient
       REAL(wp) ::  gfac_dep            !< factor
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vs                  !< mean fall velocity of snow
       REAL(wp) ::  xs                  !< mean mass of snow

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Call calculation of supersaturation over a plane ice surface
          CALL supersaturation_ice( i, j, k )
!
!--       Actual temperature, t_l is calcualted in supersaturation_ice routine:
          temp = t_l + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--       Calculating gfac_dep ( 1/ (Fk + Fd) ) see e.g.
!--       Rogers and Yau, 1989
          gfac_dep  = 1.0_wp / ( ( l_s / ( r_v * temp ) - 1.0_wp ) *                               &
                                   l_s / ( thermal_conductivity_l  * temp )                        &
                                 + r_v * temp / ( diff_coeff_l * e_si )                            &
                               )
!
!--       If there is nothing nucleated, than there is also no deposition (above -38°C)
          IF ( ns(k,j,i) > 0.0_wp )  THEN
!
!--          calculate mean mass of snow crystal
             xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )

!
!--          Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and mean
!--          diameter of snow
             vs = terminal_fall_velocity( cloud_species%snow, xs, hyrho(k) )
             ds = mean_diameter( cloud_species%snow, xs )

             IF ( ventilation_snow )  THEN
                f_vent  = a1_ven_coeff_snow + b1_ven_coeff_snow * SQRT( ds * vs )
                f_vent  = MAX( f_vent, a1_ven_coeff_snow / cloud_species%snow%a_ven )
             ELSE
                f_vent = 1.0_wp
             ENDIF
!
!--          Water vapor deposition in supersaturated regions (regarding ice)
             IF ( sat_ice > 0.0_wp )  THEN
!
!--             Calculate deposition rate assuming snow crystal shape
                deposition_rate  = 4.0_wp * pi / cloud_species%snow%cap * sat_ice * gfac_dep *     &
                                   fac_gamma * ds * f_vent * ns(k,j,i)
                deposition_rate_max  = q(k,j,i) - q_si - qc(k,j,i) - qr(k,j,i) - qf(k,j,i)
                deposition_rate = MIN( deposition_rate, deposition_rate_max / dt_micro )

                qs(k,j,i) = qs(k,j,i) + deposition_rate * dt_micro * flag
!
!--          Water vapor sublimation in supersaturated regions (regarding ice)
             ELSEIF ( sat_ice < 0.0_wp )  THEN
                deposition_rate = 4.0_wp * pi / cloud_species%snow%cap * sat_ice * gfac_dep *     &
                                   fac_gamma * ds * f_vent * ns(k,j,i)
                deposition_rate = MAX( deposition_rate, - qs(k,j,i) / dt_micro )
                qs(k,j,i) = qs(k,j,i) + deposition_rate * dt_micro * flag
             ENDIF
!
!--          Store deposition rate ice for riming process
             dep_rate_snow(k,j,i) = deposition_rate
          ENDIF
       ENDDO

    END SUBROUTINE deposition_snow_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the growth/shrink of graupel by water vapor deposition/sublimation
!> (after Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE deposition_graupel

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       LOGICAL :: ventilation_graupel = .TRUE.  !< Turn on/off ventilation effects during sublimation

       REAL(wp) ::  dg                  !< mean diameter of graupel
       REAL(wp) ::  deposition_rate     !< depositions rate
       REAL(wp) ::  deposition_rate_max !< maximum deposition rate
       REAL(wp) ::  fac_gamma = 1.0_wp  !< parameter to describe spectral distribution,
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  f_vent              !< ventilation coefficient
       REAL(wp) ::  gfac_dep            !< factor
       REAL(wp) ::  sublimation_rate    !< sublimations rate
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vg                  !< mean fall velocity of graupel
       REAL(wp) ::  xg                  !< mean mass of graupel

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Call calculation of supersaturation over a plane ice surface
                CALL supersaturation_ice( i, j, k )
!
!--             Actual temperature, t_l is calcualted in supersaturation_ice routine:
                temp = t_l + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--             Calculating gfac_dep ( 1/ (Fk + Fd) ) see e.g.
!--             Rogers and Yau, 1989
                gfac_dep  = 1.0_wp / ( ( l_s / ( r_v * temp ) - 1.0_wp ) *                         &
                                         l_s / ( thermal_conductivity_l  * temp )                  &
                                       + r_v * temp / ( diff_coeff_l * e_si )                      &
                                     )
!
!--             If there is nothing nucleated, than there is also no deposition (above -38°C)
                IF ( ng(k,j,i) > 0.0_wp )  THEN
!
!--                calculate mean mass of graupel
                   xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
!
!--                Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and
!--                mean diameter of graupel
                   vg = terminal_fall_velocity( cloud_species%graupel, xg, hyrho(k) )
                   dg = mean_diameter( cloud_species%graupel, xg )

                   IF ( ventilation_graupel )  THEN
                      f_vent  = a1_ven_coeff_graupel + b1_ven_coeff_graupel * SQRT( dg * vg )
                      f_vent  = MAX( f_vent, a1_ven_coeff_graupel / cloud_species%graupel%a_ven )
                   ELSE
                      f_vent = 1.0_wp
                   ENDIF
!
!--                Water vapor deposition in supersaturated regions (regarding ice)
                   IF ( sat_ice > 0.0_wp )  THEN
!
!--                   Calculate deposition rate assuming graupel crystal shape
                      deposition_rate  = 4.0_wp * pi / cloud_species%graupel%cap * sat_ice *       &
                                         gfac_dep * fac_gamma * dg * f_vent * ng(k,j,i)
                      deposition_rate_max  = q(k,j,i) - q_si - qc(k,j,i) - qr(k,j,i) - qf(k,j,i)
                      deposition_rate = MIN( deposition_rate, deposition_rate_max / dt_micro )

                      qg(k,j,i) = qg(k,j,i) + deposition_rate * dt_micro * flag
!
!--                Water vapor sublimation in supersaturated regions (regarding ice)
                   ELSEIF ( sat_ice < 0.0_wp )  THEN
                      sublimation_rate = 4.0_wp * pi / cloud_species%graupel%cap * sat_ice *       &
                                         gfac_dep * fac_gamma * dg * f_vent * ng(k,j,i)
                      sublimation_rate = MAX( sublimation_rate, - qg(k,j,i) / dt_micro )
                      qg(k,j,i) = qg(k,j,i) + sublimation_rate * dt_micro * flag
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE deposition_graupel

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the growth/shrink of graupel by water vapor deposition/sublimation
!> (after Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE deposition_graupel_ij( i, j )

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       LOGICAL :: ventilation_graupel = .TRUE.  !< Turn on/off ventilation effects during sublimation

       REAL(wp) ::  dg                  !< mean diameter of graupel
       REAL(wp) ::  deposition_rate     !< depositions rate
       REAL(wp) ::  deposition_rate_max !< maximum deposition rate
       REAL(wp) ::  fac_gamma = 1.0_wp  !< parameter to describe spectral distribution,
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  f_vent              !< ventilation coefficient
       REAL(wp) ::  gfac_dep            !< factor
       REAL(wp) ::  sublimation_rate    !< sublimations rate
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vg                  !< mean fall velocity of graupel
       REAL(wp) ::  xg                  !< mean mass of graupel

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Call calculation of supersaturation over a plane ice surface
          CALL supersaturation_ice( i, j, k )
!
!--       Actual temperature, t_l is calcualted in supersaturation_ice routine:
          temp = t_l + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--       Calculating gfac_dep ( 1/ (Fk + Fd) ) see e.g.
!--       Rogers and Yau, 1989
          gfac_dep  = 1.0_wp / ( ( l_s / ( r_v * temp ) - 1.0_wp ) *                               &
                                   l_s / ( thermal_conductivity_l  * temp )                        &
                                 + r_v * temp / ( diff_coeff_l * e_si )                            &
                               )
!
!--       If there is nothing nucleated, than there is also no deposition (above -38°C)
          IF ( ng(k,j,i) > 0.0_wp )  THEN
!
!--          calculate mean mass of graupel
             xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
!
!--          Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and mean
!--          diameter of graupel
             vg = terminal_fall_velocity( cloud_species%graupel, xg, hyrho(k) )
             dg = mean_diameter( cloud_species%graupel, xg )

             IF ( ventilation_graupel )  THEN
                f_vent  = a1_ven_coeff_graupel + b1_ven_coeff_graupel * SQRT( dg * vg )
                f_vent  = MAX( f_vent, a1_ven_coeff_graupel / cloud_species%graupel%a_ven )
             ELSE
                f_vent = 1.0_wp
             ENDIF
!
!--          Water vapor deposition in supersaturated regions (regarding ice)
             IF ( sat_ice > 0.0_wp )  THEN
!
!--             Calculate deposition rate assuming graupel crystal shape
                deposition_rate  = 4.0_wp * pi / cloud_species%graupel%cap * sat_ice * gfac_dep *  &
                                   fac_gamma * dg * f_vent * ng(k,j,i)
                deposition_rate_max  = q(k,j,i) - q_si - qc(k,j,i) - qr(k,j,i) - qf(k,j,i)
                deposition_rate = MIN( deposition_rate, deposition_rate_max / dt_micro )

                qg(k,j,i) = qg(k,j,i) + deposition_rate * dt_micro * flag
!
!--          Water vapor sublimation in supersaturated regions (regarding ice)
             ELSEIF ( sat_ice < 0.0_wp )  THEN
                sublimation_rate = 4.0_wp * pi / cloud_species%graupel%cap * sat_ice * gfac_dep *  &
                                   fac_gamma * dg * f_vent * ng(k,j,i)
                sublimation_rate = MAX( sublimation_rate, - qg(k,j,i) / dt_micro )
                qg(k,j,i) = qg(k,j,i) + sublimation_rate * dt_micro * flag
             ENDIF
          ENDIF
       ENDDO

    END SUBROUTINE deposition_graupel_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates melting of ice particles after SB2006. Ice particles melting immeaditely if temperature
! is high enough.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE melting_ice

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  melting_rate_n     !< melting rate for k = zeroth moment (i.e. number concentration)
       REAL(wp) ::  melting_rate_q     !< melting rate fof k = first moment (i.e. mixing ratio )
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  xi                  !< mean mass of ice crystal

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Actual temperature:
                temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)

                IF ( temp > t3  .AND.  qi(k,j,i) > 0.0_wp )  THEN
!
!--                Calculate mean mass of ice crystal
                   xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
!
!--                set ice to zero, as it's melting immeaditely
                   melting_rate_q = qi(k,j,i)
                   melting_rate_n = ni(k,j,i)

                   qi(k,j,i) = qi(k,j,i) - melting_rate_q
                   ni(k,j,i) = ni(k,j,i) - melting_rate_n
!
!--                Add number of ice crystals to rain or cloud in dependency of size
                   IF ( xi > xrmin )  THEN
                      qr(k,j,i) = qr(k,j,i) + melting_rate_q
                      nr(k,j,i) = nr(k,j,i) + melting_rate_n
                   ELSE
                      qc(k,j,i) = qc(k,j,i) + melting_rate_q
                      IF ( microphysics_morrison )  THEN
                         nc(k,j,i) = nc(k,j,i) + melting_rate_n
                      ENDIF
                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE melting_ice


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates melting of ice particles after SB2006. Ice particles melting immeaditely if temperature
! is high enough.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE melting_ice_ij( i, j )

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  melting_rate_n     !< melting rate for k = zeroth moment (i.e. number concentration)
       REAL(wp) ::  melting_rate_q     !< melting rate fof k = first moment (i.e. mixing ratio )
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  xi                  !< mean mass of ice crystal

       DO  k = nzb+1, nzt
!
!--       Actual temperature:
          temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)

          IF ( temp > t3  .AND.  qi(k,j,i) > 0.0_wp )  THEN
!
!--          Calculate mean mass of ice crystal
             xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
!
!--          set ice to zero, as it's melting immeaditely
             melting_rate_q = qi(k,j,i)
             melting_rate_n = ni(k,j,i)

             qi(k,j,i) = qi(k,j,i) - melting_rate_q
             ni(k,j,i) = ni(k,j,i) - melting_rate_n
!
!--          Add number of ice crystals to rain or cloud in dependency of size
             IF ( xi > xrmin )  THEN
                qr(k,j,i) = qr(k,j,i) + melting_rate_q
                nr(k,j,i) = nr(k,j,i) + melting_rate_n
             ELSE
                qc(k,j,i) = qc(k,j,i) + melting_rate_q
                IF ( microphysics_morrison )  THEN
                   nc(k,j,i) = nc(k,j,i) + melting_rate_n
                ENDIF
             ENDIF
          ENDIF

       ENDDO

    END SUBROUTINE melting_ice_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates melting of graupel after SB2006
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE melting_graupel

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  dg                  !< power law normalised mean graupel diameter
       REAL(wp) ::  e_a                 !< water vapor pressure
       REAL(wp) ::  e_s0                !< saturation vapor pressure at freezing point
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  fv_q                !< ventilation coefficient
       REAL(wp) ::  fh_q                !< ventilation coefficient
       REAL(wp) ::  melt_fac            !< melting rate factor
       REAL(wp) ::  melt_h              !< melting rate factor
       REAL(wp) ::  melt_v              !< melting rate factor
       REAL(wp) ::  melting_rate_n      !< melting rate for  number concentration
       REAL(wp) ::  melting_rate_q      !< melting rate for mixing ratio
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vg                  !< fall velocity graupel
       REAL(wp) ::  xg                  !< mean mass of graupel

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Actual temperature:
                temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)

                IF ( temp > t3  .AND.  qg(k,j,i) > 0.0_wp  .AND.  ng(k,j,i) > 0.0_wp )  THEN
!
!--                Calculate water vapor pressure at freezing point
                   e_s0 = magnus_tl ( t3 )
!
!--                Current vapor pressure
                   e_a = (   q(k,j,i) - ql(k,j,i) - qf(k,j,i) ) * hyp(k) /                         &
                         ( ( q(k,j,i) - ql(k,j,i) - qf(k,j,i) ) + rd_d_rv )
!
!--                Calculate mean mass of graupel crystal, limit to minimum and maximum value
!--                concerning SB-scheme
                   xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
!
!--                Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and
!--                mean diameter
                   vg = terminal_fall_velocity( cloud_species%graupel, xg, hyrho(k) )
                   dg = mean_diameter( cloud_species%graupel, xg )
!
!--                Calculate averaged ventilation coefficient for graupel
                   fv_q = a1_ven_coeff_graupel + b1_ven_coeff_graupel * SQRT( vg * dg )
!
!--                Regarding on Rasmussen and Heymsfield (1987) the ratio fh/fv is 1.05
                   fh_q = 1.05_wp * fv_q
!
!--                Calculate factors for melting rate
                   melt_fac = 2.0_wp * pi / l_m * dg * ng(k,j,i)
                   melt_h = melt_fac * thermal_conductivity_l * (temp - t3)
                   melt_v = melt_fac * diff_coeff_l * l_s / r_d * ( e_a / temp - e_s0 / t3 )
!
!--                Calculate melting rates
                   melting_rate_q = ( melt_h * fh_q + melt_v * fv_q )
!
!--                Assuming that xg is constant during melting
                   melting_rate_n = MIN( MAX( ( melting_rate_q - qg(k,j,i) / dt_micro ) / xg +     &
                                                ng(k,j,i), 0.0_wp ), ng(k,j,i) / dt_micro )
!
!--                Limit the melting rate to current number of graupel and graupel mass
                   melting_rate_q = MAX( MIN( qg(k,j,i) / dt_micro, melting_rate_q ), 0.0_wp )
                   melting_rate_n = MAX( MIN( ng(k,j,i) / dt_micro, melting_rate_n ), 0.0_wp )
!
!--                Add/substract gain/loss due to melting of ice crystals
                   qg(k,j,i) = qg(k,j,i) - melting_rate_q * dt_micro * flag
                   ng(k,j,i) = ng(k,j,i) - melting_rate_n * dt_micro * flag

                   qr(k,j,i) = qr(k,j,i) + melting_rate_q * dt_micro * flag
                   nr(k,j,i) = nr(k,j,i) + melting_rate_n * dt_micro * flag
!
!--                Adapt number of graupel
                   ng(k,j,i) = MAX( ng(k,j,i), qg(k,j,i) / xgmax )

                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE melting_graupel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates melting of graupel after SB2006
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE melting_graupel_ij( i, j )

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  dg                  !< power law normalised mean graupel diameter
       REAL(wp) ::  e_a                 !< water vapor pressure
       REAL(wp) ::  e_s0                !< saturation vapor pressure at freezing point
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  fv_q                !< ventilation coefficient
       REAL(wp) ::  fh_q                !< ventilation coefficient
       REAL(wp) ::  melt_fac            !< melting rate factor
       REAL(wp) ::  melt_h              !< melting rate factor
       REAL(wp) ::  melt_v              !< melting rate factor
       REAL(wp) ::  melting_rate_n      !< melting rate for  number concentration
       REAL(wp) ::  melting_rate_q      !< melting rate for mixing ratio
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vg                  !< fall velocity graupel
       REAL(wp) ::  xg                  !< mean mass of graupel

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Actual temperature:
          temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)

          IF ( temp > t3  .AND.  qg(k,j,i) > 0.0_wp  .AND.  ng(k,j,i) > 0.0_wp )  THEN
!
!--          Calculate water vapor pressure at freezing point
             e_s0 = magnus_tl ( t3 )
!
!--          Current vapor pressure
             e_a = (   q(k,j,i) - ql(k,j,i) - qf(k,j,i) ) * hyp(k) /                               &
                   ( ( q(k,j,i) - ql(k,j,i) - qf(k,j,i) ) + rd_d_rv )
!
!--          Calculate mean mass of graupel crystal, limit to minimum and maximum value concerning
!--          SB-scheme
             xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
!
!--          Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and mean
!--          diameter
             vg = terminal_fall_velocity( cloud_species%graupel, xg, hyrho(k) )
             dg = mean_diameter( cloud_species%graupel, xg )
!
!--          Calculate averaged ventilation coefficient for graupel
             fv_q = a1_ven_coeff_graupel + b1_ven_coeff_graupel * SQRT( vg * dg )
!
!--          Regarding on Rasmussen and Heymsfield (1987) the ratio fh/fv is 1.05
             fh_q = 1.05_wp * fv_q
!
!--          Calculate factors for melting rate
             melt_fac = 2.0_wp * pi / l_m * dg * ng(k,j,i)
             melt_h = melt_fac * thermal_conductivity_l * (temp - t3)
             melt_v = melt_fac * diff_coeff_l * l_s / r_d * ( e_a / temp - e_s0 / t3 )
!
!--          Calculate melting rates
             melting_rate_q = ( melt_h * fh_q + melt_v * fv_q )
!
!--          Assuming that xg is constant during melting
             melting_rate_n = MIN( MAX( ( melting_rate_q - qg(k,j,i) / dt_micro ) / xg +           &
                                          ng(k,j,i), 0.0_wp ), ng(k,j,i) / dt_micro )
!
!--          Limit the melting rate to current number of graupel and graupel mass
             melting_rate_q = MAX( MIN( qg(k,j,i) / dt_micro, melting_rate_q ), 0.0_wp )
             melting_rate_n = MAX( MIN( ng(k,j,i) / dt_micro, melting_rate_n ), 0.0_wp )
!
!--          Add/substract gain/loss due to melting of ice crystals
             qg(k,j,i) = qg(k,j,i) - melting_rate_q * dt_micro * flag
             ng(k,j,i) = ng(k,j,i) - melting_rate_n * dt_micro * flag

             qr(k,j,i) = qr(k,j,i) + melting_rate_q * dt_micro * flag
             nr(k,j,i) = nr(k,j,i) + melting_rate_n * dt_micro * flag
!
!--          Adapt number of graupel
             ng(k,j,i) = MAX( ng(k,j,i), qg(k,j,i) / xgmax )

          ENDIF

       ENDDO

    END SUBROUTINE melting_graupel_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates melting of snow after SB2006
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE melting_snow

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  ds                  !< power law normalised mean graupel diameter
       REAL(wp) ::  e_a                 !< water vapor pressure
       REAL(wp) ::  e_s0                !< saturation vapor pressure at freezing point
       REAL(wp) ::  fv_q                !<
       REAL(wp) ::  fh_q                !<
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  melt_fac            !< melting factor
       REAL(wp) ::  melt_h              !< metling factor
       REAL(wp) ::  melt_v              !< metling factor
       REAL(wp) ::  melting_rate_n      !< melting rate for mixing ratio
       REAL(wp) ::  melting_rate_q      !< melting rate for number concentration
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vs                  !< fall velocity graupel
       REAL(wp) ::  xs                  !< mean mass of graupel

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Actual temperature:
                temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--             Melting occurs above 0°C and if snow is present
                IF ( temp > t3  .AND.  qs(k,j,i) > 0.0_wp  .AND.  ns(k,j,i) > 0.0_wp )  THEN
!
!--                calculate water vapor pressure at freezing point
                   e_s0 = magnus_tl ( t3 )
!
!--                Current vapor pressure
                   e_a = (   q(k,j,i) - ql(k,j,i) - qf(k,j,i) ) * hyp(k) /                          &
                         ( ( q(k,j,i) - ql(k,j,i) - qf(k,j,i) ) + rd_d_rv )
!
!--                Calculate mean mass of graupel crystal, limit to minimum and maximum value
!--                concerning SB-scheme
                   xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )
!
!--                Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1)
                   vs = terminal_fall_velocity( cloud_species%snow, xs, hyrho(k) )
                   ds = mean_diameter( cloud_species%snow, xs )
!
!--                Calculate ventilation coefficients
                   fv_q = a1_ven_coeff_snow + b1_ven_coeff_snow * SQRT( vs * ds )
!
!--                Regarding on Rasmussen and Heymsfield (1987) the ratio fh/fv is 1.05
                   fh_q = 1.05_wp * fv_q
!
!--                Calculate melting factor and rates
                   melt_fac = 2.0_wp * pi / l_m * ds * ns(k,j,i)
                   melt_h = melt_fac * thermal_conductivity_l * (temp - t3)
                   melt_v = melt_fac * diff_coeff_l * l_s / r_d * ( e_a / temp - e_s0 / t3)
!
!--                Assuming that xs is constant during melting
                   melting_rate_q = ( melt_h * fh_q + melt_v * fv_q )
                   melting_rate_n = MIN ( MAX( ( melting_rate_q - qs(k,j,i) ) / xs + &
                                              ns(k,j,i), 0.0_wp ),  ns(k,j,i) )
!
!--                Limit melting rates
                   melting_rate_q = MIN( qs(k,j,i) / dt_micro, MAX( melting_rate_q, 0.0_wp ) )
                   melting_rate_n = MIN( ns(k,j,i) / dt_micro, MAX( melting_rate_n, 0.0_wp ) )
!
!--                snow melts instantaneously at 10°C
                   IF ( temp - t3 > 10.0_wp) THEN
                     melting_rate_q = qs(k,j,i) / dt_micro
                     melting_rate_n = ns(k,j,i) / dt_micro
                   ENDIF
!
!--                Substract snow due to melting
                   qs(k,j,i) = qs(k,j,i) - melting_rate_q * dt_micro * flag
                   qr(k,j,i) = qr(k,j,i) + melting_rate_q * dt_micro * flag

                   ns(k,j,i) = ns(k,j,i) - melting_rate_n * dt_micro * flag
                   nr(k,j,i) = nr(k,j,i) + melting_rate_n * dt_micro * flag
!
!--                Adapt number of snow crystals
                   ns(k,j,i) = MAX( ns(k,j,i), qs(k,j,i) / xsmax )

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE melting_snow


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates melting of snow after SB2006
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE melting_snow_ij( i, j )

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  ds                  !< power law normalised mean graupel diameter
       REAL(wp) ::  e_a                 !< water vapor pressure
       REAL(wp) ::  e_s0                !< saturation vapor pressure at freezing point
       REAL(wp) ::  fv_q                !<
       REAL(wp) ::  fh_q                !<
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  melt_fac            !< melting factor
       REAL(wp) ::  melt_h              !< metling factor
       REAL(wp) ::  melt_v              !< metling factor
       REAL(wp) ::  melting_rate_n      !< melting rate for mixing ratio
       REAL(wp) ::  melting_rate_q      !< melting rate for number concentration
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vs                  !< fall velocity graupel
       REAL(wp) ::  xs                  !< mean mass of graupel

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Actual temperature:
          temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--       Melting occurs above 0°C and if snow is present
          IF ( temp > t3  .AND.  qs(k,j,i) > 0.0_wp  .AND.  ns(k,j,i) > 0.0_wp )  THEN
!
!--          calculate water vapor pressure at freezing point
             e_s0 = magnus_tl ( t3 )
!
!--          Current vapor pressure
             e_a = (   q(k,j,i) - ql(k,j,i) - qf(k,j,i) ) * hyp(k) /                               &
                   ( ( q(k,j,i) - ql(k,j,i) - qf(k,j,i) ) + rd_d_rv )
!
!--          Calculate mean mass of graupel crystal, limit to minimum and maximum value concerning
!--          SB-scheme
             xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )
!
!--          Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1)
             vs = terminal_fall_velocity( cloud_species%snow, xs, hyrho(k) )
             ds = mean_diameter( cloud_species%snow, xs )
!
!--          Calculate ventilation coefficients
             fv_q = a1_ven_coeff_snow + b1_ven_coeff_snow * SQRT( vs * ds )
!
!--          Regarding on Rasmussen and Heymsfield (1987) the ratio fh/fv is 1.05
             fh_q = 1.05_wp * fv_q
!
!--          Calculate melting factor and rates
             melt_fac = 2.0_wp * pi / l_m * ds * ns(k,j,i)
             melt_h = melt_fac * thermal_conductivity_l * (temp - t3)
             melt_v = melt_fac * diff_coeff_l * l_s / r_d * ( e_a / temp - e_s0 / t3)
!
!--          Assuming that xs is constant during melting
             melting_rate_q = ( melt_h * fh_q + melt_v * fv_q )
             melting_rate_n = MIN ( MAX( ( melting_rate_q - qs(k,j,i) ) / xs + ns(k,j,i), 0.0_wp ),&
                                    ns(k,j,i) )
!
!--          Limit melting rates
             melting_rate_q = MIN( qs(k,j,i) / dt_micro, MAX( melting_rate_q, 0.0_wp ) )
             melting_rate_n = MIN( ns(k,j,i) / dt_micro, MAX( melting_rate_n, 0.0_wp ) )
!
!--          snow melts instantaneously at 10°C
             IF ( temp - t3 > 10.0_wp) THEN
               melting_rate_q = qs(k,j,i) / dt_micro
               melting_rate_n = ns(k,j,i) / dt_micro
             ENDIF
!
!--          Substract snow due to melting
             qs(k,j,i) = qs(k,j,i) - melting_rate_q * dt_micro * flag
             qr(k,j,i) = qr(k,j,i) + melting_rate_q * dt_micro * flag

             ns(k,j,i) = ns(k,j,i) - melting_rate_n * dt_micro * flag
             nr(k,j,i) = nr(k,j,i) + melting_rate_n * dt_micro * flag
!
!--          Adapt number of snow crystals
             ns(k,j,i) = MAX( ns(k,j,i), qs(k,j,i) / xsmax )

          ENDIF
       ENDDO

    END SUBROUTINE melting_snow_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates evaporation of melting graupel particles after SB2006
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE evaporation_graupel

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  dg                  !< power law normalised mean graupel diameter
       REAL(wp) ::  evap_graupel        !< evaporation rate
       REAL(wp) ::  f_v                 !< ventilation coefficient
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  g_fac               !< factor for evaporation
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vg                  !< fall velocity graupel
       REAL(wp) ::  xg                  !< mean mass of graupel

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Actual temperature:
                temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)

                IF ( qg(k,j,i) > 0.0_wp  .AND.  ng(k,j,i) > 0.0_wp  .AND.  temp > t3 )  THEN
!
!--                Call calculation of supersaturation over liquid water
                   CALL supersaturation ( i, j, k )
                   g_fac = 4.0_wp * pi / ( l_s**2 / ( thermal_conductivity_l * r_d * t3**2 ) +     &
                           r_d * t3 / ( diff_coeff_l * e_s ) )

!--                Calculate mean mass of graupel, limit to minimum and maximum value concerning
!--                SB-scheme
                   xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
!
!--                Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and
!--                mean diameter
                   vg = terminal_fall_velocity( cloud_species%graupel, xg, hyrho(k) )
                   dg = mean_diameter( cloud_species%graupel, xg )
!
!--                ventilation coefficient
                   f_v = a1_ven_coeff_graupel + b1_ven_coeff_graupel * SQRT( vg * dg )
!
!--                evapoartion rate
                   evap_graupel = g_fac * ng(k,j,i) / cloud_species%graupel%cap * dg * f_v * sat
!
!--                Limit evaporation rate
                   evap_graupel = MIN( MAX( -evap_graupel, 0.0_wp ), qg(k,j,i) / dt_micro )
!
!--                Evaporation of graupel
                   qg(k,j,i) = qg(k,j,i) - evap_graupel * dt_micro * flag

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE evaporation_graupel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates evaporation of melting graupel particles after SB2006
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE evaporation_graupel_ij( i, j )

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  dg                  !< power law normalised mean graupel diameter
       REAL(wp) ::  evap_graupel        !< evaporation rate
       REAL(wp) ::  f_v                 !< ventilation coefficient
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  g_fac               !< factor for evaporation
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vg                  !< fall velocity graupel
       REAL(wp) ::  xg                  !< mean mass of graupel

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Actual temperature:
          temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)

          IF ( qg(k,j,i) > 0.0_wp  .AND.  ng(k,j,i) > 0.0_wp  .AND.  temp > t3 )  THEN
!
!--          Call calculation of supersaturation over liquid water
             CALL supersaturation ( i, j, k )
             g_fac = 4.0_wp * pi / ( l_s**2 / ( thermal_conductivity_l * r_d * t3**2 ) +           &
                     r_d * t3 / ( diff_coeff_l * e_s ) )

!--          Calculate mean mass of graupel, limit to minimum and maximum value concerning
!--          SB-scheme
             xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
!
!--          Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and
!--          mean diameter
             vg = terminal_fall_velocity( cloud_species%graupel, xg, hyrho(k) )
             dg = mean_diameter( cloud_species%graupel, xg )
!
!--          ventilation coefficient
             f_v = a1_ven_coeff_graupel + b1_ven_coeff_graupel * SQRT( vg * dg )
!
!--          evapoartion rate
             evap_graupel = g_fac * ng(k,j,i) / cloud_species%graupel%cap * dg * f_v * sat
!
!--          Limit evaporation rate
             evap_graupel = MIN( MAX( -evap_graupel, 0.0_wp ), qg(k,j,i) / dt_micro )
!
!--          Evaporation of graupel
             qg(k,j,i) = qg(k,j,i) - evap_graupel * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE evaporation_graupel_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates evaporation of melting snow particles after SB2006
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE evaporation_snow

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  ds                  !< power law normalised mean snow diameter
       REAL(wp) ::  evap_snow           !< evaporation rate
       REAL(wp) ::  f_v                 !< ventilation coefficient
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  g_fac               !< factor for evaporation
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vs                  !< fall velocity snow
       REAL(wp) ::  xs                  !< mean mass of snow

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Actual temperature:
                temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)

                IF ( qs(k,j,i) > 0.0_wp  .AND.  ns(k,j,i) > 0.0_wp  .AND.  temp > t3 )  THEN
!
!--                Call calculation of supersaturation over liquid water
                   CALL supersaturation ( i, j, k )
                   g_fac = 4.0_wp * pi / ( l_s**2 / ( thermal_conductivity_l * r_d * t3**2 ) + &
                           r_d * t3 / ( diff_coeff_l * e_s ) )

!--                Calculate mean mass of snow crystal, limit to minimum and maximum value concerning
!--                SB-scheme
                   xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )
!
!--                Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and
!--                mean diameter
                   vs = terminal_fall_velocity( cloud_species%snow, xs, hyrho(k) )
                   ds = mean_diameter( cloud_species%snow, xs )
!
!--                ventilation coefficient
                   f_v = a1_ven_coeff_snow + b1_ven_coeff_snow * SQRT( vs * ds )
!
!--                evapoartion rate
                   evap_snow = g_fac * ns(k,j,i) / cloud_species%snow%cap * ds * f_v * sat
!
!--                Limit evaporation rate
                   evap_snow = MIN( MAX( -evap_snow, 0.0_wp ), qs(k,j,i) / dt_micro )
!
!--                Evaporation of snow
                   qs(k,j,i) = qs(k,j,i) - evap_snow * dt_micro * flag

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE evaporation_snow


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
! Calculates evaporation of melting snow particles after SB2006
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE evaporation_snow_ij( i, j )

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  ds                  !< power law normalised mean snow diameter
       REAL(wp) ::  evap_snow           !< evaporation rate
       REAL(wp) ::  f_v                 !< ventilation coefficient
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  g_fac               !< factor for evaporation
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vs                  !< fall velocity snow
       REAL(wp) ::  xs                  !< mean mass of snow

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Actual temperature:
          temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)

          IF ( qs(k,j,i) > 0.0_wp  .AND.  ns(k,j,i) > 0.0_wp  .AND.  temp > t3 )  THEN
!
!--          Call calculation of supersaturation over liquid water
             CALL supersaturation ( i, j, k )
             g_fac = 4.0_wp * pi / ( l_s**2 / ( thermal_conductivity_l * r_d * t3**2 ) + &
                     r_d * t3 / ( diff_coeff_l * e_s ) )

!--          Calculate mean mass of snow crystal, limit to minimum and maximum value concerning
!--          SB-scheme
             xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )
!
!--          Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1) and
!--          mean diameter
             vs = terminal_fall_velocity( cloud_species%snow, xs, hyrho(k) )
             ds = mean_diameter( cloud_species%snow, xs )
!
!--          ventilation coefficient
             f_v = a1_ven_coeff_snow + b1_ven_coeff_snow * SQRT( vs * ds )
!
!--          evapoartion rate
             evap_snow = g_fac * ns(k,j,i) / cloud_species%snow%cap * ds * f_v * sat
!
!--          Limit evaporation rate
             evap_snow = MIN( MAX( -evap_snow, 0.0_wp ), qs(k,j,i) / dt_micro )
!
!--          Evaporation of snow
             qs(k,j,i) = qs(k,j,i) - evap_snow * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE evaporation_snow_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Heteorogenous freezing of rain droplets (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE heterogeneous_freezing_rain

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  a_het = 0.2_wp      !< factor freezing 0.2 kg-1 s-1
       REAL(wp) ::  b_het = 0.65_wp     !< factor freezing K^-1
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  freezing_rate_n     !< freezing rate number concentratration
       REAL(wp) ::  freezing_rate_q     !< freezing rate mixing ratio
       REAL(wp) ::  j_het               !< factor of temperature dependent heteorogenous freezing function
       REAL(wp) ::  temp                !< temperature
       REAL(wp) ::  xr                  !< average mass rain droplet

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Actual temperature:
                temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--             Freezing occurs below 273.15K
                IF ( temp < t3 )  THEN

                   IF ( qr(k,j,i) >= eps_sb )  THEN
!
!--                   Mean weight of rain drops
                      xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
!
!--                   factor of temperature dependent heteorogenous freezing function
                      j_het = MAX( a_het * EXP( b_het * ( t3 - temp ) - 1.0_wp ) / rho_l , 0.0_wp )
!
!--                   Calculate freezing rate according Seifert and Beheng, 2006 Eq.44-47
                      freezing_rate_n = qr(k,j,i) * j_het
                      freezing_rate_q = 20.0_wp * qr(k,j,i) * xr * j_het

                      freezing_rate_n = MIN( freezing_rate_n, nr(k,j,i) / dt_micro )
                      freezing_rate_q = MIN( freezing_rate_q, qr(k,j,i) / dt_micro )

                      nr(k,j,i) = nr(k,j,i) - freezing_rate_n * dt_micro * flag
                      qr(k,j,i) = qr(k,j,i) - freezing_rate_q * dt_micro * flag

                      ng(k,j,i) = ng(k,j,i) + freezing_rate_n * dt_micro * flag
                      qg(k,j,i) = qg(k,j,i) + freezing_rate_q * dt_micro * flag

                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE heterogeneous_freezing_rain


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Heteorogenous freezing of rain droplets (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE heterogeneous_freezing_rain_ij( i, j )

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  a_het = 0.2_wp      !< factor freezing 0.2 kg-1 s-1
       REAL(wp) ::  b_het = 0.65_wp     !< factor freezing K^-1
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  freezing_rate_n     !< freezing rate number concentratration
       REAL(wp) ::  freezing_rate_q     !< freezing rate mixing ratio
       REAL(wp) ::  j_het               !< factor of temperature dependent heteorogenous freezing function
       REAL(wp) ::  temp                !< temperature
       REAL(wp) ::  xr                  !< average mass rain droplet

       DO  k = nzb+1, nzt

          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Actual temperature:
          temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--       Freezing occurs below 273.15K
          IF ( temp < t3 )  THEN

             IF ( qr(k,j,i) >= eps_sb )  THEN
!
!--             Mean weight of rain drops
                xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
!
!--             factor of temperature dependent heteorogenous freezing function
                j_het = MAX( a_het * EXP( b_het * ( t3 - temp ) - 1.0_wp ) / rho_l , 0.0_wp )
!
!--             Calculate freezing rate according Seifert and Beheng, 2006 Eq.44-47
                freezing_rate_n = qr(k,j,i) * j_het
                freezing_rate_q = 20.0_wp * qr(k,j,i) * xr * j_het

                freezing_rate_n = MIN( freezing_rate_n, nr(k,j,i) / dt_micro )
                freezing_rate_q = MIN( freezing_rate_q, qr(k,j,i) / dt_micro )

                nr(k,j,i) = nr(k,j,i) - freezing_rate_n * dt_micro * flag
                qr(k,j,i) = qr(k,j,i) - freezing_rate_q * dt_micro * flag

                ng(k,j,i) = ng(k,j,i) + freezing_rate_n * dt_micro * flag
                qg(k,j,i) = qg(k,j,i) + freezing_rate_q * dt_micro * flag

             ENDIF
          ENDIF
       ENDDO

    END SUBROUTINE heterogeneous_freezing_rain_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Homogeneous freezing of cloud droplets (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE homogeneous_freezing_cloud

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  nc_c                  !< local dummy number concentration clouds
       REAL(wp) ::  tc                    !< actual temperature in °C
       REAL(wp) ::  temp                  !< actual temperature
       REAL(wp) ::  j_hom                 !< Jhom temperature dependent freezing function (kg-1 s-1)
       REAL(wp) ::  freezing_rate_n      !< freezing rate for zeroth moment (m-3 s-1)
       REAL(wp) ::  freezing_rate_q      !< freezing rate for first moment (kg kg-1 s-1)
       REAL(wp) ::  flag                  !< flag to indicate first grid level above
       REAL(wp) ::  xc                    !< mean mass cloud droplets

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt

                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Actual temperature:
                temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--             Freezing occurs below 273.15K
                IF ( temp < t3 )  THEN
                   tc = temp - t3
                   IF ( microphysics_morrison )  THEN
                      nc_c = nc(k,j,i)
                   ELSE
                      nc_c = nc_const
                   ENDIF
                   IF ( qc(k,j,i) > eps_sb  .AND.  nc_c > 0.0_wp  .AND.  tc < -30.0_wp ) THEN
!
!--                   Instantaneous freezing below -50°C
                      IF ( tc < -50.0_wp )  THEN
                         qi(k,j,i) = qi(k,j,i) + qc(k,j,i)
                         ni(k,j,i) = ni(k,j,i) + nc_c
                         qc(k,j,i) = 0.0_wp
                         IF ( microphysics_morrison )  THEN
                            nc(k,j,i) = 0.0_wp
                         ENDIF
                      ELSE
!
!--                      Mean weight of cloud drops
                         xc = mean_mass( cloud_species%cloud, qc(k,j,i), nc_c, hyrho(k) )
!
!--                      Calculate temperature dependent function for homogeneous freezing cloud
!--                      after Jeffrey und Austin (1997) and Cotton und Field (2002)
                         IF ( tc > -30.0_wp )  THEN
                            j_hom = 1.0E6_wp / rho_l * 10**( -7.63_wp -2.996_wp * ( tc + 30.0_wp ))
                         ELSE
                            j_hom = 1.0E6_wp / rho_l * 10**( -243.4_wp -14.75_wp * tc -0.307_wp *  &
                                    tc**2 -0.00287_wp *tc**3 -0.0000102_wp * tc**4 )
                         ENDIF
                         j_hom = MERGE( j_hom, 0.0_wp, j_hom > 1.0E-20_wp )
!
!--                      Calculate freezing rate according Seifert and Beheng, 2006 Eq.49-50
                         freezing_rate_n = j_hom * qc(k,j,i)
                         freezing_rate_q = j_hom * qc(k,j,i) * xc * cz_cloud
!
!--                      Limit freezing rates to available water content
                         freezing_rate_n = MIN( freezing_rate_n, nc_c / dt_micro )
                         freezing_rate_q = MIN( freezing_rate_q, qc(k,j,i) / dt_micro )
!
!--                      Substract/Add liquid and ice to species
                         qc(k,j,i) = qc(k,j,i) - freezing_rate_q * dt_micro * flag
                         IF ( microphysics_morrison )  THEN
                            nc(k,j,i) = nc(k,j,i) - freezing_rate_n * dt_micro * flag
                         ENDIF
                         ni(k,j,i) = ni(k,j,i) + freezing_rate_n * dt_micro * flag
                         qi(k,j,i) = qi(k,j,i) + freezing_rate_q * dt_micro * flag
                      ENDIF

                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE homogeneous_freezing_cloud


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Homogeneous freezing of cloud droplets (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE homogeneous_freezing_cloud_ij( i, j )

       INTEGER(iwp) ::  i               !< loop index
       INTEGER(iwp) ::  j               !< loop index
       INTEGER(iwp) ::  k               !< loop index

       REAL(wp) ::  nc_c                  !< local dummy number concentration clouds
       REAL(wp) ::  j_hom                 !< Jhom temperature dependent freezing function (kg-1 s-1)
       REAL(wp) ::  freezing_rate_n       !< freezing rate for zeroth moment (number concentration, m-3 s-1)
       REAL(wp) ::  freezing_rate_q       !< freezing rate for first moment (mixing ratio, kg kg-1 s-1)
       REAL(wp) ::  flag                  !< flag to indicate first grid level above
       REAL(wp) ::  tc                    !< actual temperature in °C
       REAL(wp) ::  temp                  !< actual temperature
       REAL(wp) ::  xc                    !< mean mass cloud droplets

       DO  k = nzb+1, nzt

          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Actual temperature:
          temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--       Freezing occurs below 273.15K
          IF ( temp < t3 )  THEN
             tc = temp - t3
             IF ( microphysics_morrison )  THEN
                nc_c = nc(k,j,i)
             ELSE
                nc_c = nc_const
             ENDIF
             IF ( qc(k,j,i) > eps_sb  .AND.  nc_c > 0.0_wp  .AND.  tc < -30.0_wp )  THEN
!
!--             Instantaneous freezing below -50°C
                IF ( tc < -50.0_wp )  THEN
                   qi(k,j,i) = qi(k,j,i) + qc(k,j,i)
                   ni(k,j,i) = ni(k,j,i) + nc_c
                   qc(k,j,i) = 0.0_wp
                   IF ( microphysics_morrison )  THEN
                      nc(k,j,i) = 0.0_wp
                   ENDIF
                ELSE
!
!--                Mean weight of cloud drops
                   xc = mean_mass( cloud_species%cloud, qc(k,j,i), nc_c, hyrho(k) )
!
!--                Calculate temperature dependent function for homogeneous freezing cloud
!--                after Jeffrey und Austin (1997) and Cotton und Field (2002)
                   IF ( tc > -30.0_wp )  THEN
                      j_hom = 1.0E6_wp / rho_l * 10**( -7.63_wp -2.996_wp * ( tc + 30.0_wp ))
                   ELSE
                      j_hom = 1.0E6_wp / rho_l * 10**( -243.4_wp - 14.75_wp * tc - 0.307_wp *      &
                              tc**2 - 0.00287_wp * tc**3 -0.0000102_wp * tc**4 )
                   ENDIF
                   j_hom = MERGE( j_hom, 0.0_wp, j_hom > 1.0E-20_wp )
!
!--                Calculate freezing rate according Seifert and Beheng, 2006 Eq.49-50
                   freezing_rate_n = j_hom * qc(k,j,i)
                   freezing_rate_q = j_hom * qc(k,j,i) * xc * cz_cloud
!
!--                Limit freezing rates to available water content
                   freezing_rate_n = MIN( freezing_rate_n, nc_c / dt_micro )
                   freezing_rate_q = MIN( freezing_rate_q, qc(k,j,i) / dt_micro )
!
!--                Substract/Add liquid and ice to species
                   qc(k,j,i) = qc(k,j,i) - freezing_rate_q * dt_micro * flag
                   IF ( microphysics_morrison )  THEN
                      nc(k,j,i) = nc(k,j,i) - freezing_rate_n * dt_micro * flag
                   ENDIF
                   ni(k,j,i) = ni(k,j,i) + freezing_rate_n * dt_micro * flag
                   qi(k,j,i) = qi(k,j,i) + freezing_rate_q * dt_micro * flag
                ENDIF

             ENDIF
          ENDIF
       ENDDO

    END SUBROUTINE homogeneous_freezing_cloud_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion rate (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE autoconversion

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  alpha_cc          !<
       REAL(wp)     ::  autocon           !<
       REAL(wp)     ::  dissipation       !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  k_au              !<
       REAL(wp)     ::  l_mix             !<
       REAL(wp)     ::  nc_c              !<
       REAL(wp)     ::  nu_c              !<
       REAL(wp)     ::  phi_au            !<
       REAL(wp)     ::  r_cc              !<
       REAL(wp)     ::  rc                !<
       REAL(wp)     ::  re_lambda         !<
       REAL(wp)     ::  sigma_cc          !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( microphysics_morrison )  THEN
                   nc_c = nc(k,j,i)
                ELSE
                   nc_c = nc_const
                ENDIF

                IF ( qc(k,j,i) > eps_sb  .AND.  nc_c > eps_mr )  THEN

                   k_au = k_cc / ( 20.0_wp * x0 )
!
!--                Intern time scale of coagulation (Seifert and Beheng, 2006):
!--                (1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) ))
                   tau_cloud = MAX( 1.0_wp - qc(k,j,i) / ( qr(k,j,i) + qc(k,j,i) ), 0.0_wp )
!
!--                Universal function for autoconversion process
!--                (Seifert and Beheng, 2006):
                   phi_au = 600.0_wp * tau_cloud**0.68_wp * ( 1.0_wp - tau_cloud**0.68_wp )**3
!
!--                Shape parameter of gamma distribution (Geoffroy et al., 2010):
!--                (Use constant nu_c = 1.0_wp instead?)
                   nu_c = 1.0_wp !MAX( 0.0_wp, 1580.0_wp * hyrho(k) * qc(k,j,i) - 0.28_wp )
!
!--                Mean weight of cloud droplets:
                   xc = hyrho(k) * qc(k,j,i) / nc_c
!
!--                Parameterized turbulence effects on autoconversion
!--                (Seifert, Nuijens and Stevens, 2010)
                   IF ( collision_turbulence )  THEN
!
!--                   Weight averaged radius of cloud droplets:
                      rc = 0.5_wp * ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )

                      alpha_cc = ( a_1 + a_2 * nu_c ) / ( 1.0_wp + a_3 * nu_c )
                      r_cc     = ( b_1 + b_2 * nu_c ) / ( 1.0_wp + b_3 * nu_c )
                      sigma_cc = ( c_1 + c_2 * nu_c ) / ( 1.0_wp + c_3 * nu_c )
!
!--                   Mixing length (neglecting distance to ground and stratification)
                      l_mix = ( dx * dy * dzu(k) )**( 1.0_wp / 3.0_wp )
!
!--                   Limit dissipation rate according to Seifert, Nuijens and Stevens (2010)
                      dissipation = MIN( 0.06_wp, diss(k,j,i) )
!
!--                   Compute Taylor-microscale Reynolds number:
                      re_lambda = 6.0_wp / 11.0_wp * ( l_mix / c_const )**( 2.0_wp / 3.0_wp ) *    &
                                  SQRT( 15.0_wp / kin_vis_air ) * dissipation**( 1.0_wp / 6.0_wp )
!
!--                   The factor of 1.0E4 is needed to convert the dissipation rate from m2 s-3 to
!--                   cm2 s-3.
                      k_au = k_au * ( 1.0_wp +                                                     &
                                      dissipation * 1.0E4_wp *                                     &
                                      ( re_lambda * 1.0E-3_wp )**0.25_wp *                         &
                                      ( alpha_cc * EXP( -1.0_wp * ( ( rc -                         &
                                                                      r_cc ) /                     &
                                                        sigma_cc )**2                              &
                                                      ) + beta_cc                                  &
                                      )                                                            &
                                    )
                   ENDIF
!
!--                Autoconversion rate (Seifert and Beheng, 2006):
                   autocon = k_au * ( nu_c + 2.0_wp ) * ( nu_c + 4.0_wp )    /                     &
                             ( nu_c + 1.0_wp )**2 * qc(k,j,i)**2 * xc**2     *                     &
                             ( 1.0_wp + phi_au / ( 1.0_wp - tau_cloud )**2 ) * rho_surface
                   autocon = MIN( autocon, qc(k,j,i) / dt_micro )

                   qr(k,j,i) = qr(k,j,i) + autocon * dt_micro * flag
                   qc(k,j,i) = qc(k,j,i) - autocon * dt_micro * flag
                   nr(k,j,i) = nr(k,j,i) + autocon / x0 * hyrho(k) * dt_micro * flag
                   IF ( microphysics_morrison )  THEN
                      nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i), 2.0_wp *                             &
                                  autocon / x0 * hyrho(k) * dt_micro * flag )
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE autoconversion


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion rate (Seifert and Beheng, 2006). Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE autoconversion_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  alpha_cc          !<
       REAL(wp)     ::  autocon           !<
       REAL(wp)     ::  dissipation       !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  k_au              !<
       REAL(wp)     ::  l_mix             !<
       REAL(wp)     ::  nc_c              !<
       REAL(wp)     ::  nu_c              !<
       REAL(wp)     ::  phi_au            !<
       REAL(wp)     ::  r_cc              !<
       REAL(wp)     ::  rc                !<
       REAL(wp)     ::  re_lambda         !<
       REAL(wp)     ::  sigma_cc          !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
          IF ( microphysics_morrison ) THEN
             nc_c = nc(k,j,i)
          ELSE
             nc_c = nc_const
          ENDIF

          IF ( qc(k,j,i) > eps_sb  .AND.  nc_c > eps_mr )  THEN

             k_au = k_cc / ( 20.0_wp * x0 )
!
!--          Intern time scale of coagulation (Seifert and Beheng, 2006):
!--          (1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) ))
             tau_cloud = MAX( 1.0_wp - qc(k,j,i) / ( qr(k,j,i) + qc(k,j,i) ), 0.0_wp )
!
!--          Universal function for autoconversion process (Seifert and Beheng, 2006):
             phi_au = 600.0_wp * tau_cloud**0.68_wp * ( 1.0_wp - tau_cloud**0.68_wp )**3
!
!--          Shape parameter of gamma distribution (Geoffroy et al., 2010):
!--          (Use constant nu_c = 1.0_wp instead?)
             nu_c = 1.0_wp !MAX( 0.0_wp, 1580.0_wp * hyrho(k) * qc(k,j,i) - 0.28_wp )
!
!--          Mean weight of cloud droplets:
             xc = hyrho(k) * qc(k,j,i) / nc_c
!
!--          Parameterized turbulence effects on autoconversion (Seifert, Nuijens and Stevens, 2010)
             IF ( collision_turbulence )  THEN
!
!--             Weight averaged radius of cloud droplets:
                rc = 0.5_wp * ( xc * dpirho_l )**( 1.0_wp / 3.0_wp )

                alpha_cc = ( a_1 + a_2 * nu_c ) / ( 1.0_wp + a_3 * nu_c )
                r_cc     = ( b_1 + b_2 * nu_c ) / ( 1.0_wp + b_3 * nu_c )
                sigma_cc = ( c_1 + c_2 * nu_c ) / ( 1.0_wp + c_3 * nu_c )
!
!--             Mixing length (neglecting distance to ground and stratification)
                l_mix = ( dx * dy * dzu(k) )**( 1.0_wp / 3.0_wp )
!
!--             Limit dissipation rate according to Seifert, Nuijens and Stevens (2010)
                dissipation = MIN( 0.06_wp, diss(k,j,i) )
!
!--             Compute Taylor-microscale Reynolds number:
                re_lambda = 6.0_wp / 11.0_wp * ( l_mix / c_const )**( 2.0_wp / 3.0_wp ) *          &
                            SQRT( 15.0_wp / kin_vis_air ) * dissipation**( 1.0_wp / 6.0_wp )
!
!--             The factor of 1.0E4 is needed to convert the dissipation rate from m2 s-3 to
!--             cm2 s-3.
                k_au = k_au * ( 1.0_wp +                                                           &
                                dissipation * 1.0E4_wp *                                           &
                                ( re_lambda * 1.0E-3_wp )**0.25_wp *                               &
                                ( alpha_cc * EXP( -1.0_wp * ( ( rc - r_cc ) /                      &
                                                  sigma_cc )**2                                    &
                                                ) + beta_cc                                        &
                                )                                                                  &
                              )
             ENDIF
!
!--          Autoconversion rate (Seifert and Beheng, 2006):
             autocon = k_au * ( nu_c + 2.0_wp ) * ( nu_c + 4.0_wp )    /                           &
                       ( nu_c + 1.0_wp )**2 * qc(k,j,i)**2 * xc**2     *                           &
                       ( 1.0_wp + phi_au / ( 1.0_wp - tau_cloud )**2 ) * rho_surface
             autocon = MIN( autocon, qc(k,j,i) / dt_micro )
             qr(k,j,i) = qr(k,j,i) + autocon * dt_micro                 * flag
             qc(k,j,i) = qc(k,j,i) - autocon * dt_micro                 * flag
             nr(k,j,i) = nr(k,j,i) + autocon / x0 * hyrho(k) * dt_micro * flag
             IF ( microphysics_morrison )  THEN
                nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i), 2.0_wp *                                   &
                            autocon / x0 * hyrho(k) * dt_micro * flag )
             ENDIF

          ENDIF

       ENDDO

    END SUBROUTINE autoconversion_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion process (Kessler, 1969).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE autoconversion_kessler


       IMPLICIT NONE

       INTEGER(iwp) ::  i      !< loop index
       INTEGER(iwp) ::  j      !< loop index
       INTEGER(iwp) ::  k      !< loop index
       INTEGER(iwp) ::  k_wall !< topgraphy top index

       REAL(wp)    ::  dqdt_precip !<
       REAL(wp)    ::  flag        !< flag to mask topography grid points

       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Determine vertical index of topography top
             k_wall = topo_top_ind(j,i,0)
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( qc(k,j,i) > ql_crit )  THEN
                   dqdt_precip = prec_time_const * ( qc(k,j,i) - ql_crit )
                ELSE
                   dqdt_precip = 0.0_wp
                ENDIF

                qc(k,j,i) = qc(k,j,i) - dqdt_precip * dt_micro * flag
                q(k,j,i)  = q(k,j,i)  - dqdt_precip * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) + dqdt_precip * dt_micro * lv_d_cp *                         &
                                        d_exner(k)             * flag
!
!--             Compute the rain rate (stored on surface grid point)
                prr(k_wall,j,i) = prr(k_wall,j,i) + dqdt_precip * dzw(k) * flag

             ENDDO
          ENDDO
       ENDDO

   END SUBROUTINE autoconversion_kessler

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Autoconversion process (Kessler, 1969).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE autoconversion_kessler_ij( i, j )


       IMPLICIT NONE

       INTEGER(iwp) ::  i      !< loop index
       INTEGER(iwp) ::  j      !< loop index
       INTEGER(iwp) ::  k      !< loop index
       INTEGER(iwp) ::  k_wall !< topography top index

       REAL(wp)    ::  dqdt_precip       !<
       REAL(wp)    ::  flag              !< flag to indicate first grid level above surface

!
!--    Determine vertical index of topography top
       k_wall = topo_top_ind(j,i,0)
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          IF ( qc(k,j,i) > ql_crit )  THEN
             dqdt_precip = prec_time_const * ( qc(k,j,i) - ql_crit )
          ELSE
             dqdt_precip = 0.0_wp
          ENDIF

          qc(k,j,i) = qc(k,j,i) - dqdt_precip * dt_micro * flag
          q(k,j,i)  = q(k,j,i)  - dqdt_precip * dt_micro * flag
          pt(k,j,i) = pt(k,j,i) + dqdt_precip * dt_micro * lv_d_cp * d_exner(k) * flag
!
!--       Compute the rain rate (stored on surface grid point)
          prr(k_wall,j,i) = prr(k_wall,j,i) + dqdt_precip * dzw(k) * flag

       ENDDO

    END SUBROUTINE autoconversion_kessler_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Accretion rate (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE accretion

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  accr              !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  k_cr              !<
       REAL(wp)     ::  nc_c              !<
       REAL(wp)     ::  phi_ac            !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( microphysics_morrison )  THEN
                   nc_c = nc(k,j,i)
                ELSE
                   nc_c = nc_const
                ENDIF

                IF ( ( qc(k,j,i) > eps_sb )  .AND.  ( qr(k,j,i) > eps_sb )                         &
                                             .AND.  ( nc_c > eps_mr   ) )  THEN
!
!--                Intern time scale of coagulation (Seifert and Beheng, 2006):
                   tau_cloud = 1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) )
!
!--                Universal function for accretion process (Seifert and
!--                Beheng, 2001):
                   phi_ac = ( tau_cloud / ( tau_cloud + 5.0E-5_wp ) )**4
!--                Mean weight of cloud droplets:
                   xc = MAX( hyrho(k) * qc(k,j,i) / nc_c, xcmin)

!
!--                Parameterized turbulence effects on autoconversion
!--                (Seifert, Nuijens and Stevens, 2010). The factor of 1.0E4 is needed to convert
!--                the dissipation rate (diss) from m2 s-3 to cm2 s-3.
                   IF ( collision_turbulence )  THEN
                      k_cr = k_cr0 * ( 1.0_wp + 0.05_wp *                                          &
                                       MIN( 600.0_wp, diss(k,j,i) * 1.0E4_wp )**0.25_wp            &
                                     )
                   ELSE
                      k_cr = k_cr0
                   ENDIF
!
!--                Accretion rate (Seifert and Beheng, 2006):
                   accr = k_cr * qc(k,j,i) * qr(k,j,i) * phi_ac * SQRT( rho_surface * hyrho(k) )
                   accr = MIN( accr, qc(k,j,i) / dt_micro )
                   qr(k,j,i) = qr(k,j,i) + accr * dt_micro * flag
                   qc(k,j,i) = qc(k,j,i) - accr * dt_micro * flag
                   IF ( microphysics_morrison )  THEN
                      nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i),                                      &
                                                   accr / xc * hyrho(k) * dt_micro * flag)
                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE accretion

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Accretion rate (Seifert and Beheng, 2006). Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE accretion_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  accr              !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  k_cr              !<
       REAL(wp)     ::  nc_c              !<
       REAL(wp)     ::  phi_ac            !<
       REAL(wp)     ::  tau_cloud         !<
       REAL(wp)     ::  xc                !<


       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
          IF ( microphysics_morrison )  THEN
             nc_c = nc(k,j,i)
          ELSE
             nc_c = nc_const
          ENDIF

          IF ( ( qc(k,j,i) > eps_sb )  .AND.  ( qr(k,j,i) > eps_sb )  .AND.                        &
               ( nc_c > eps_mr   ) )  THEN
!
!--          Intern time scale of coagulation (Seifert and Beheng, 2006):
             tau_cloud = 1.0_wp - qc(k,j,i) / ( qc(k,j,i) + qr(k,j,i) )
!
!--          Universal function for accretion process (Seifert and Beheng, 2001):
             phi_ac = ( tau_cloud / ( tau_cloud + 5.0E-5_wp ) )**4

!
!--          Mean weight of cloud drops
             !xc = mean_mass( cloud_species%cloud, qc(k,j,i), nc_c, hyrho(k) )
             xc = MAX( hyrho(k) * qc(k,j,i) / nc_c, xcmin)

!
!--          Parameterized turbulence effects on autoconversion
!--          (Seifert, Nuijens and Stevens, 2010). The factor of 1.0E4 is needed to convert the
!--          dissipation rate (diss) from m2 s-3 to cm2 s-3.
             IF ( collision_turbulence )  THEN
                k_cr = k_cr0 * ( 1.0_wp + 0.05_wp *                                                &
                                 MIN( 600.0_wp, diss(k,j,i) * 1.0E4_wp )**0.25_wp                  &
                               )
             ELSE
                k_cr = k_cr0
             ENDIF
!
!--          Accretion rate (Seifert and Beheng, 2006):
             accr = k_cr * qc(k,j,i) * qr(k,j,i) * phi_ac * SQRT( rho_surface * hyrho(k) )
             accr = MIN( accr, qc(k,j,i) / dt_micro )
             qr(k,j,i) = qr(k,j,i) + accr * dt_micro * flag
             qc(k,j,i) = qc(k,j,i) - accr * dt_micro * flag
             IF ( microphysics_morrison )  THEN
                nc(k,j,i) = nc(k,j,i) - MIN( nc(k,j,i), accr / xc * hyrho(k) * dt_micro * flag )
             ENDIF
          ENDIF

       ENDDO

    END SUBROUTINE accretion_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collisional breakup rate of rain droplets (Seifert, 2008) and rain selfcollection.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_breakup_rain

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  breakup           !<
       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  phi_br            !<
       REAL(wp)     ::  selfcoll          !<

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( qr(k,j,i) > eps_sb  .AND.  nr(k,j,i) > 0.0_wp  )  THEN
!
!--                Selfcollection rate (Seifert and Beheng, 2001):
                   selfcoll = k_rr * nr(k,j,i) * qr(k,j,i) * SQRT( hyrho(k) * rho_surface )
!
!--                Weight averaged diameter of rain drops:
                   dr = ( hyrho(k) * qr(k,j,i) / nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--                Collisional breakup rate (Seifert, 2008):
                   IF ( dr >= 0.3E-3_wp )  THEN
                      phi_br  = k_br * ( dr - 1.1E-3_wp )
                      breakup = selfcoll * ( phi_br + 1.0_wp )
                   ELSE
                      breakup = 0.0_wp
                   ENDIF

                   selfcoll = MAX( breakup - selfcoll, -nr(k,j,i) / dt_micro )
                   nr(k,j,i) = nr(k,j,i) + selfcoll * dt_micro * flag
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE selfcollection_breakup_rain


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collisional breakup rate of rain droplets (Seifert, 2008) and rain selfcollection.
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_breakup_rain_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  breakup           !<
       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  phi_br            !<
       REAL(wp)     ::  selfcoll          !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          IF ( qr(k,j,i) > eps_sb  .AND.  nr(k,j,i) > 0.0_wp )  THEN
!
!--          Selfcollection rate (Seifert and Beheng, 2001):
             selfcoll = k_rr * nr(k,j,i) * qr(k,j,i) * SQRT( hyrho(k) * rho_surface )
!
!--          Weight averaged diameter of rain drops:
             dr = ( hyrho(k) * qr(k,j,i) / nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--          Collisional breakup rate (Seifert, 2008):
             IF ( dr >= 0.3E-3_wp )  THEN
                phi_br  = k_br * ( dr - 1.1E-3_wp )
                breakup = selfcoll * ( phi_br + 1.0_wp )
             ELSE
                breakup = 0.0_wp
             ENDIF

             selfcoll = MAX( breakup - selfcoll, -nr(k,j,i) / dt_micro )
             nr(k,j,i) = nr(k,j,i) + selfcoll * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE selfcollection_breakup_rain_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Selfcollection graupel (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_graupel

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  dg                  !< power law normalised mean graupel diameter
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  selfcollection_rate !< selfcollection rate
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vg                  !< fall velocity graupel
       REAL(wp) ::  xg                  !< mean mass of graupel

       REAL(wp), PARAMETER ::  ecoll_gg     = 0.10_wp  !< collision efficiency for graupel selfcollection
       REAL(wp), PARAMETER ::  ecoll_gg_wet = 0.40_wp  !< in case of wet graupel (T>273.15K)

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( qg(k,j,i) > eps_sb  .AND.  ng(k,j,i) > 0.0_wp )  THEN

!--                Calculate mean mass of graupel crystal, limit to minimum and maximum value
!--                SB-scheme
                   xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
!
!--                Calculate temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1)
                   vg = terminal_fall_velocity(cloud_species%graupel, xg, hyrho(k) )
                   dg = mean_diameter( cloud_species%graupel, xg )
!
!                  graupel_coeffs%sc_coll_n  = pi / 8.0_wp * delta_n * theta_n
                   selfcollection_rate = coll_coeff_graupel_self * ng(k,j,i)**2 * dg**2 * vg
!
!--                sticking efficiency does only distinguish dry and wet based on T_3
                   selfcollection_rate = selfcollection_rate *                                     &
                                            MERGE( ecoll_gg_wet, ecoll_gg, temp > t3 )
!
!--                Limit selfcollection rate
                   selfcollection_rate = MIN( selfcollection_rate, ng(k,j,i) / dt_micro )
!
!--                Reduction of number concentration due to selfcollection
                   ng(k,j,i) = ng(k,j,i) - selfcollection_rate * dt_micro * flag

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE selfcollection_graupel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Selfcollection graupel (Seifert and Beheng, 2006). Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_graupel_ij( i, j )

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  dg                  !< power law normalised mean graupel diameter
       REAL(wp) ::  flag                !< flag to indicate first grid level above
       REAL(wp) ::  selfcollection_rate !< selfcollection rate
       REAL(wp) ::  temp                !< actual temperature
       REAL(wp) ::  vg                  !< fall velocity graupel
       REAL(wp) ::  xg                  !< mean mass of graupel

       REAL(wp), PARAMETER ::  ecoll_gg     = 0.10_wp  !< collision efficiency for graupel selfcollection
       REAL(wp), PARAMETER ::  ecoll_gg_wet = 0.40_wp  !< in case of wet graupel (T>273.15K)

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          IF ( qg(k,j,i) > eps_sb  .AND.  ng(k,j,i) > 0.0_wp )  THEN

!--          Calculate mean mass of graupel crystal, limit to minimum and maximum value concerning
!--          SB-scheme
             xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
!
!--          Calculate temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Calculate mean fall velocity (Seifert and Beheng, 2006: Eq.:33 and Table 1)
             vg = terminal_fall_velocity(cloud_species%graupel, xg, hyrho(k) )
             dg = mean_diameter( cloud_species%graupel, xg )
!
!            graupel_coeffs%sc_coll_n  = pi / 8.0_wp * delta_n * theta_n
             selfcollection_rate = coll_coeff_graupel_self * ng(k,j,i)**2 * dg**2 * vg
!
!--          sticking efficiency does only distinguish dry and wet based on T_3
             selfcollection_rate = selfcollection_rate * MERGE( ecoll_gg_wet, ecoll_gg, temp > t3 )
!
!--          Limit selfcollection rate
             selfcollection_rate = MIN( selfcollection_rate, ng(k,j,i) / dt_micro )
!
!--          Reduction of number concentration due to selfcollection
             ng(k,j,i) = ng(k,j,i) - selfcollection_rate * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE selfcollection_graupel_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Selfcollection ice (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_ice

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  di                      !< diameter ice
       REAL(wp) ::  di_crit = 100.0E-6_wp   !< critical diameter ice for selfcollection
       REAL(wp) ::  d_conv_ii = 75.0E-6_wp  !< diameter for conversion to snow
       REAL(wp) ::  e_coll                  !< collision efficences
       REAL(wp) ::  flag                    !< flag to indicate first grid level above topography
       REAL(wp) ::  q_crit_ii = 1.0E-6_wp   !< critical mass
       REAL(wp) ::  selfcollection_rate_n   !< selfcollection rate number concentration
       REAL(wp) ::  selfcollection_rate_q   !< selfcollection rate mixing ratio
       REAL(wp) ::  temp                    !< actual temperature
       REAL(wp) ::  vi                      !< terminal fall velocity ice crystal
       REAL(wp) ::  xi                      !< mean mass ice
       REAL(wp) ::  x_conv_ii               !< mass for conversion to snow

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( qi(k,j,i) > q_crit_ii  .AND.  ni(k,j,i) > 0.0_wp )  THEN
!
!--                Calculate mean mass
                   xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
!
!--                Calculate mean diameter of ice particle
                   di = mean_diameter( cloud_species%ice, xi )
!
!--                Selfcollection only takes place if ice particles are large enough
                   IF ( di < di_crit )  CYCLE

                   x_conv_ii = ( d_conv_ii / cloud_species%snow%a )**(1.0_wp / cloud_species%snow%b )
!
!--                Calculate temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Calculate temperaute dependent collision efficiency after Cotton et al., 1986
                   e_coll = MIN( 10**( 0.035_wp * ( temp - t3 ) - 0.7_wp ), 0.2_wp )
!
!--                Calculate terminal fall velocity of ice particles
                   vi = terminal_fall_velocity(cloud_species%ice, xi, hyrho(k) )
!
!--                calcualte self collection rates
                   selfcollection_rate_n = pi / 4.0_wp * e_coll * coll_coeff_ice_self_delta_n *    &
                           ni(k,j,i)**2 * di**2  * SQRT( coll_coeff_ice_self_theta_n * vi**2  +    &
                           2.0_wp * cloud_species%ice%sigma_v**2 )

                   selfcollection_rate_q = pi / 4.0_wp * e_coll * coll_coeff_ice_self_delta_q *    &
                           ni(k,j,i)* qi(k,j,i) * di**2  * SQRT( coll_coeff_ice_self_theta_q *     &
                           vi**2  +  2.0_wp * cloud_species%ice%sigma_v**2 )
!
!--                Limit selfcollection rates
                   selfcollection_rate_q = MIN( selfcollection_rate_q, qi(k,j,i) / dt_micro )
                   selfcollection_rate_n = MIN( MIN( selfcollection_rate_n,                        &
                                                     qi(k,j,i) / x_conv_ii ),                      &
                                                ni(k,j,i) / dt_micro )
!
!--                Mass transformation to snow and number concentration reduction due to
!--                selfcollection
                   qi(k,j,i) = qi(k,j,i) - selfcollection_rate_q * dt_micro * flag
                   qs(k,j,i) = qs(k,j,i) + selfcollection_rate_q * dt_micro * flag

                   ni(k,j,i) = ni(k,j,i) - selfcollection_rate_n * dt_micro * flag
                   ns(k,j,i) = ns(k,j,i) + selfcollection_rate_n * dt_micro * flag / 2.0_wp

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE selfcollection_ice


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Selfcollection ice (Seifert and Beheng, 2006). Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_ice_ij( i, j )

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  di                      !< diameter ice
       REAL(wp) ::  di_crit = 100.0E-6_wp   !< critical diameter ice for selfcollection
       REAL(wp) ::  d_conv_ii = 75.0E-6_wp  !< diameter for conversion to snow
       REAL(wp) ::  e_coll                  !< collision efficences
       REAL(wp) ::  flag                    !< flag to indicate first grid level above topography
       REAL(wp) ::  q_crit_ii = 1.0E-6_wp   !< critical mass
       REAL(wp) ::  selfcollection_rate_n   !< selfcollection rate number concentration
       REAL(wp) ::  selfcollection_rate_q   !< selfcollection rate mixing ratio
       REAL(wp) ::  temp                    !< actual temperature
       REAL(wp) ::  vi                      !< terminal fall velocity ice crystal
       REAL(wp) ::  xi                      !< mean mass ice
       REAL(wp) ::  x_conv_ii               !< mass for conversion to snow

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          IF ( qi(k,j,i) > q_crit_ii  .AND.  ni(k,j,i) > 0.0_wp )  THEN
!
!--          Calculate mean mass
             xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
!
!--          Calculate mean diameter of ice particle
             di = mean_diameter( cloud_species%ice, xi )
!
!--          Selfcollection only takes place if ice particles are large enough
             IF ( di < di_crit )  CYCLE

             x_conv_ii = ( d_conv_ii / cloud_species%snow%a )**(1.0_wp / cloud_species%snow%b )
!
!--          Calculate temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Calculate temperaute dependent collision efficiency after Cotton et al., 1986
             e_coll = MIN( 10**( 0.035_wp * ( temp - t3 ) - 0.7_wp ), 0.2_wp )
!
!--          Calculate terminal fall velocity of ice particles
             vi = terminal_fall_velocity(cloud_species%ice, xi, hyrho(k) )
!
!--          calcualte self collection rates
             selfcollection_rate_n = pi / 4.0_wp * e_coll * coll_coeff_ice_self_delta_n *          &
                     ni(k,j,i)**2 * di**2  * SQRT( coll_coeff_ice_self_theta_n * vi**2  +          &
                     2.0_wp * cloud_species%ice%sigma_v**2 )

             selfcollection_rate_q = pi / 4.0_wp * e_coll * coll_coeff_ice_self_delta_q *          &
                     ni(k,j,i)* qi(k,j,i) * di**2  * SQRT( coll_coeff_ice_self_theta_q * vi**2  +  &
                     2.0_wp * cloud_species%ice%sigma_v**2 )
!
!--          Limit selfcollection rates
             selfcollection_rate_q = MIN( selfcollection_rate_q, qi(k,j,i) / dt_micro )
             selfcollection_rate_n = MIN( MIN( selfcollection_rate_n, qi(k,j,i) / x_conv_ii ),     &
                                          ni(k,j,i) / dt_micro )
!
!--          Mass transformation to snow and number concentration reduction due to selfcollection
             qi(k,j,i) = qi(k,j,i) - selfcollection_rate_q * dt_micro * flag
             qs(k,j,i) = qs(k,j,i) + selfcollection_rate_q * dt_micro * flag

             ni(k,j,i) = ni(k,j,i) - selfcollection_rate_n * dt_micro * flag
             ns(k,j,i) = ns(k,j,i) + selfcollection_rate_n * dt_micro * flag / 2.0_wp

          ENDIF
       ENDDO

    END SUBROUTINE selfcollection_ice_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Selfcollection snow (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_snow

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  ds                      !< diameter snow
       REAL(wp) ::  e_coll                  !< collision efficiency
       REAL(wp) ::  flag                    !< flag to indicate first grid level above topography
       REAL(wp) ::  selfcollection_rate_n   !< selfcollection rate number concentration
       REAL(wp) ::  temp                    !< actual temperature
       REAL(wp) ::  vs                      !< terminal fall velocites snow
       REAL(wp) ::  xs                      !< mean mass snow

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( qs(k,j,i) > eps_sb  .AND.  ns(k,j,i) > 0.0_wp )  THEN
!
!--                Calculate temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Temperature depending collision efficiency
                   e_coll = MAX( 0.1_wp, MIN( EXP( 0.09_wp * (temp - t3) ), 1.0_wp ) )
!
!--                Calculate mean mass
                   xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )
!
!--                Calculate mean diameter of ice particle
                   ds = mean_diameter( cloud_species%snow, xs )
!
!--                Calculate terminal fall velocity of ice particles
                   vs = terminal_fall_velocity(cloud_species%snow, xs, hyrho(k) )
!
!--                Calculate selfcollection
                   selfcollection_rate_n = pi / 8.0_wp * e_coll * ns(k,j,i)**2 *                   &
                      coll_coeff_snow_self_delta * ds**2 *                                         &
                      SQRT( coll_coeff_snow_self_theta * vs**2 + 2.0_wp *                          &
                            cloud_species%snow%sigma_v**2 )
!
!--                Limit selfcollection rate
                   selfcollection_rate_n = MIN( selfcollection_rate_n, ns(k,j,i) / dt_micro )
!
!--                Reduction of number concentration due to selfcollection
                   ns(k,j,i) = ns(k,j,i) - selfcollection_rate_n * dt_micro * flag

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE selfcollection_snow


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Selfcollection snow (Seifert and Beheng, 2006). Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE selfcollection_snow_ij( i, j )

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp) ::  ds                      !< diameter snow
       REAL(wp) ::  e_coll                  !< collision efficiency
       REAL(wp) ::  flag                    !< flag to indicate first grid level above topography
       REAL(wp) ::  selfcollection_rate_n   !< selfcollection rate number concentration
       REAL(wp) ::  temp                    !< actual temperature
       REAL(wp) ::  vs                      !< terminal fall velocites snow
       REAL(wp) ::  xs                      !< mean mass snow

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          IF ( qs(k,j,i) > eps_sb  .AND.  ns(k,j,i) > 0.0_wp )  THEN
!
!--          Calculate temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Temperature depending collision efficiency
             e_coll = MAX( 0.1_wp, MIN( EXP( 0.09_wp * (temp - t3) ), 1.0_wp ) )
!
!--          Calculate mean mass
             xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )
!
!--          Calculate mean diameter of ice particle
             ds = mean_diameter( cloud_species%snow, xs )
!
!--          Calculate terminal fall velocity of ice particles
             vs = terminal_fall_velocity(cloud_species%snow, xs, hyrho(k) )
!
!--          Calculate selfcollection
             selfcollection_rate_n = pi / 8.0_wp * e_coll * ns(k,j,i)**2 *                         &
                coll_coeff_snow_self_delta * ds**2 *                                               &
                SQRT( coll_coeff_snow_self_theta * vs**2 + 2.0_wp * cloud_species%snow%sigma_v**2 )
!
!--          Limit selfcollection rate
             selfcollection_rate_n = MIN( selfcollection_rate_n, ns(k,j,i) / dt_micro )
!
!--          Reduction of number concentration due to selfcollection
             ns(k,j,i) = ns(k,j,i) - selfcollection_rate_n * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE selfcollection_snow_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collection of type a + b -> a. Here: graupel + ice -> graupel (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE collection_graupel_ice

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  coll_rate_n   !< collection rate for number concentration
       REAL(wp) ::  coll_rate_q   !< collection factor for mixing ratio
       REAL(wp) ::  di            !< diameter ice
       REAL(wp) ::  dg            !< diameter graupel
       REAL(wp) ::  e_coll        !< collection efficiency
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vg            !< terminal fall velocitiy graupel
       REAL(wp) ::  vi            !< terminal fall velocitiy ice
       REAL(wp) ::  xg            !< mean mass graupel
       REAL(wp) ::  xi            !< mean mass ice

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Apply collection only above threshold
                IF ( qi(k,j,i) > eps_sb_coll  .AND.  qg(k,j,i) > eps_sb_coll )  THEN
!
!--                Calculate temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Sticking efficiency of Lin (1983)
                   e_coll = MIN( EXP( 0.09_wp * (temp - t3) ), 1.0_wp )
!
!--                Calculate mean mass, mean diameter and mean terminal fall velocites for both ice
!--                and graupel
                   xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
                   xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )

                   di = mean_diameter( cloud_species%ice, xi )
                   dg = mean_diameter( cloud_species%graupel, xg )

                   vi = terminal_fall_velocity(cloud_species%ice, xi, hyrho(k) )
                   vg = terminal_fall_velocity(cloud_species%graupel, xg, hyrho(k) )
!
!--                Calculate collection rate for number concentration
                   coll_rate_n = pi / 4.0_wp * ng(k,j,i) * ni(k,j,i) * e_coll *                    &
                             ( coll_graupel_ice_delta_n_aa * dg**2 +                               &
                               coll_graupel_ice_delta_n_ab * dg * di +                             &
                               coll_graupel_ice_delta_n_bb * di**2                                 &
                             ) *                                                                   &
                         SQRT( coll_graupel_ice_theta_n_aa * vg**2 -                               &
                               coll_graupel_ice_theta_n_ab * vg * vi +                             &
                               coll_graupel_ice_theta_n_bb * vi**2 +                               &
                               cloud_species%ice%sigma_v**2                                        &
                             )
!
!--                Calculate collection rate for mixing ratio
                   coll_rate_q = pi / 4.0_wp * ng(k,j,i) * qi(k,j,i) * e_coll *                    &
                             ( coll_graupel_ice_delta_q_aa * dg**2 +                               &
                               coll_graupel_ice_delta_q_ab * dg * di +                             &
                               coll_graupel_ice_delta_q_bb * di**2                                 &
                             ) *                                                                   &
                         SQRT( coll_graupel_ice_theta_q_aa * vg**2 -                               &
                               coll_graupel_ice_theta_q_ab * vg * vi +                             &
                               coll_graupel_ice_theta_q_bb * vi**2 +                               &
                               cloud_species%ice%sigma_v**2                                        &
                             )
!
!--                Limit collection rates
                   coll_rate_n = MIN( ni(k,j,i) / dt_micro, coll_rate_n )
                   coll_rate_q = MIN( qi(k,j,i) / dt_micro, coll_rate_q )
!
!--                Add/substract rates from species
                   qg(k,j,i) = qg(k,j,i) + coll_rate_q * dt_micro * flag
                   qi(k,j,i) = qi(k,j,i) - coll_rate_q * dt_micro * flag
                   ni(k,j,i) = ni(k,j,i) - coll_rate_n * dt_micro * flag

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE collection_graupel_ice


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collection of type a + b -> a. Here: graupel + ice -> graupel (Seifert and Beheng, 2006).
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE collection_graupel_ice_ij( i, j )

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  coll_rate_n   !< collection rate for number concentration
       REAL(wp) ::  coll_rate_q   !< collection factor for mixing ratio
       REAL(wp) ::  di            !< diameter ice
       REAL(wp) ::  dg            !< diameter graupel
       REAL(wp) ::  e_coll        !< collection efficiency
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vg            !< terminal fall velocitiy graupel
       REAL(wp) ::  vi            !< terminal fall velocitiy ice
       REAL(wp) ::  xg            !< mean mass graupel
       REAL(wp) ::  xi            !< mean mass ice

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Apply collection only above threshold
          IF ( qi(k,j,i) > eps_sb_coll  .AND.  qg(k,j,i) > eps_sb_coll )  THEN
!
!--          Calculate temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Sticking efficiency of Lin (1983)
             e_coll = MIN( EXP( 0.09_wp * (temp - t3) ), 1.0_wp )
!
!--          Calculate mean mass, mean diameter and mean terminal fall velocites for both ice
!--          and graupel
             xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
             xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )

             di = mean_diameter( cloud_species%ice, xi )
             dg = mean_diameter( cloud_species%graupel, xg )

             vi = terminal_fall_velocity(cloud_species%ice, xi, hyrho(k) )
             vg = terminal_fall_velocity(cloud_species%graupel, xg, hyrho(k) )
!
!--          Calculate collection rate for number concentration
             coll_rate_n = pi / 4.0_wp * ng(k,j,i) * ni(k,j,i) * e_coll *                          &
                       ( coll_graupel_ice_delta_n_aa * dg**2 +                                     &
                         coll_graupel_ice_delta_n_ab * dg * di +                                   &
                         coll_graupel_ice_delta_n_bb * di**2                                       &
                       ) *                                                                         &
                   SQRT( coll_graupel_ice_theta_n_aa * vg**2 -                                     &
                         coll_graupel_ice_theta_n_ab * vg * vi +                                   &
                         coll_graupel_ice_theta_n_bb * vi**2 +                                     &
                         cloud_species%ice%sigma_v**2                                              &
                       )
!
!--          Calculate collection rate for mixing ratio
             coll_rate_q = pi / 4.0_wp * ng(k,j,i) * qi(k,j,i) * e_coll *                          &
                       ( coll_graupel_ice_delta_q_aa * dg**2 +                                     &
                         coll_graupel_ice_delta_q_ab * dg * di +                                   &
                         coll_graupel_ice_delta_q_bb * di**2                                       &
                       ) *                                                                         &
                   SQRT( coll_graupel_ice_theta_q_aa * vg**2 -                                     &
                         coll_graupel_ice_theta_q_ab * vg * vi +                                   &
                         coll_graupel_ice_theta_q_bb * vi**2 +                                     &
                         cloud_species%ice%sigma_v**2                                              &
                       )
!
!--          Limit collection rates
             coll_rate_n = MIN( ni(k,j,i) / dt_micro, coll_rate_n )
             coll_rate_q = MIN( qi(k,j,i) / dt_micro, coll_rate_q )
!
!--          Add/substract rates from species
             qg(k,j,i) = qg(k,j,i) + coll_rate_q * dt_micro * flag
             qi(k,j,i) = qi(k,j,i) - coll_rate_q * dt_micro * flag
             ni(k,j,i) = ni(k,j,i) - coll_rate_n * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE collection_graupel_ice_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collection of type a + b -> a. Here: snow + ice -> snow (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE collection_snow_ice

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  coll_rate_n   !< collection rate for number concentration
       REAL(wp) ::  coll_rate_q   !< collection factor for mixing ratio
       REAL(wp) ::  di            !< diameter ice
       REAL(wp) ::  ds            !< diameter snow
       REAL(wp) ::  e_coll        !< collection efficiency
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vs            !< terminal fall velocitiy snow
       REAL(wp) ::  vi            !< terminal fall velocitiy ice
       REAL(wp) ::  xs            !< mean mass snow
       REAL(wp) ::  xi            !< mean mass ice

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Apply collection only above threshold
                IF ( qi(k,j,i) > eps_sb_coll  .AND.  qs(k,j,i) > eps_sb_coll )  THEN
!
!--                Calculate temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Sticking efficiency of Lin (1983)
                   e_coll = MIN( EXP( 0.09_wp * (temp - t3) ), 1.0_wp )
!
!--                Calculate mean mass, mean diameter and mean terminal fall velocites for both ice
!--                and snow
                   xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
                   xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )

                   di = mean_diameter( cloud_species%ice,  xi )
                   ds = mean_diameter( cloud_species%snow, xs )

                   vi = terminal_fall_velocity(cloud_species%ice,  xi, hyrho(k) )
                   vs = terminal_fall_velocity(cloud_species%snow, xs, hyrho(k) )
!
!--                Calculate collection rate for number concentration
                   coll_rate_n = pi / 4.0_wp * ns(k,j,i) * ni(k,j,i) * e_coll *                    &
                             ( coll_snow_ice_delta_n_aa * ds**2 +                                  &
                               coll_snow_ice_delta_n_ab * ds * di +                                &
                               coll_snow_ice_delta_n_bb * di**2                                    &
                             ) *                                                                   &
                         SQRT( coll_snow_ice_theta_n_aa * vs**2 -                                  &
                               coll_snow_ice_theta_n_ab * vs * vi +                                &
                               coll_snow_ice_theta_n_bb * vi**2 +                                  &
                               cloud_species%ice%sigma_v**2                                        &
                             )
!
!--                Calculate collection rate for mixing ratio
                   coll_rate_q = pi / 4.0_wp * ns(k,j,i) * qi(k,j,i) * e_coll *                    &
                             ( coll_snow_ice_delta_q_aa * ds**2 +                                  &
                               coll_snow_ice_delta_q_ab * ds * di +                                &
                               coll_snow_ice_delta_q_bb * di**2                                    &
                             ) *                                                                   &
                         SQRT( coll_snow_ice_theta_q_aa * vs**2 -                                  &
                               coll_snow_ice_theta_q_ab * vs * vi +                                &
                               coll_snow_ice_theta_q_bb * vi**2 +                                  &
                               cloud_species%ice%sigma_v**2                                        &
                             )
!
!--                Limit collection rates
                   coll_rate_n = MIN( ni(k,j,i) / dt_micro, coll_rate_n )
                   coll_rate_q = MIN( qi(k,j,i) / dt_micro, coll_rate_q )
!
!--                Add/substract rates from species
                   qs(k,j,i) = qs(k,j,i) + coll_rate_q * dt_micro * flag
                   qi(k,j,i) = qi(k,j,i) - coll_rate_q * dt_micro * flag
                   ni(k,j,i) = ni(k,j,i) - coll_rate_n * dt_micro * flag

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE collection_snow_ice


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collection of type a + b -> a. Here: snow + ice -> snow (Seifert and Beheng, 2006).
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE collection_snow_ice_ij( i, j )

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  coll_rate_n   !< collection rate for number concentration
       REAL(wp) ::  coll_rate_q   !< collection factor for mixing ratio
       REAL(wp) ::  di            !< diameter ice
       REAL(wp) ::  ds            !< diameter snow
       REAL(wp) ::  e_coll        !< collection efficiency
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vs            !< terminal fall velocitiy snow
       REAL(wp) ::  vi            !< terminal fall velocitiy ice
       REAL(wp) ::  xs            !< mean mass snow
       REAL(wp) ::  xi            !< mean mass ice

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Apply collection only above threshold
          IF ( qi(k,j,i) > eps_sb_coll  .AND.  qs(k,j,i) > eps_sb_coll )  THEN
!
!--          Calculate temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Sticking efficiency of Lin (1983)
             e_coll = MIN( EXP( 0.09_wp * (temp - t3) ), 1.0_wp )
!
!--          Calculate mean mass, mean diameter and mean terminal fall velocites for both ice
!--          and snow
             xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
             xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )

             di = mean_diameter( cloud_species%ice,  xi )
             ds = mean_diameter( cloud_species%snow, xs )

             vi = terminal_fall_velocity(cloud_species%ice,  xi, hyrho(k) )
             vs = terminal_fall_velocity(cloud_species%snow, xs, hyrho(k) )
!
!--          Calculate collection rate for number concentration
             coll_rate_n = pi / 4.0_wp * ns(k,j,i) * ni(k,j,i) * e_coll *                          &
                       ( coll_snow_ice_delta_n_aa * ds**2 +                                        &
                         coll_snow_ice_delta_n_ab * ds * di +                                      &
                         coll_snow_ice_delta_n_bb * di**2                                          &
                       ) *                                                                         &
                   SQRT( coll_snow_ice_theta_n_aa * vs**2 -                                        &
                         coll_snow_ice_theta_n_ab * vs * vi +                                      &
                         coll_snow_ice_theta_n_bb * vi**2 +                                        &
                         cloud_species%ice%sigma_v**2                                              &
                       )
!
!--          Calculate collection rate for mixing ratio
             coll_rate_q = pi / 4.0_wp * ns(k,j,i) * qi(k,j,i) * e_coll *                          &
                       ( coll_snow_ice_delta_q_aa * ds**2 +                                        &
                         coll_snow_ice_delta_q_ab * ds * di +                                      &
                         coll_snow_ice_delta_q_bb * di**2                                          &
                       ) *                                                                         &
                   SQRT( coll_snow_ice_theta_q_aa * vs**2 -                                        &
                         coll_snow_ice_theta_q_ab * vs * vi +                                      &
                         coll_snow_ice_theta_q_bb * vi**2 +                                        &
                         cloud_species%ice%sigma_v**2                                              &
                       )
!
!--          Limit collection rates
             coll_rate_n = MIN( ni(k,j,i) / dt_micro, coll_rate_n )
             coll_rate_q = MIN( qi(k,j,i) / dt_micro, coll_rate_q )
!
!--          Add/substract rates from species
             qs(k,j,i) = qs(k,j,i) + coll_rate_q * dt_micro * flag
             qi(k,j,i) = qi(k,j,i) - coll_rate_q * dt_micro * flag
             ni(k,j,i) = ni(k,j,i) - coll_rate_n * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE collection_snow_ice_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collection of type a + b -> a. Here: graupel + snow -> snow (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE collection_graupel_snow

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  coll_rate_n   !< collection rate for number concentration
       REAL(wp) ::  coll_rate_q   !< collection factor for mixing ratio
       REAL(wp) ::  dg            !< diameter graupel
       REAL(wp) ::  ds            !< diameter snow
       REAL(wp) ::  e_coll        !< collection efficiency
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vg            !< terminal fall velocitiy graupel
       REAL(wp) ::  vs            !< terminal fall velocitiy snow
       REAL(wp) ::  xg            !< mean mass graupel
       REAL(wp) ::  xs            !< mean mass snow

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Apply collection only above threshold
                IF ( qg(k,j,i) > eps_sb_coll  .AND.  qs(k,j,i) > eps_sb_coll )  THEN
!
!--                Calculate temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Sticking efficiency of Lin (1983)
                   e_coll = MIN( EXP( 0.09_wp * (temp - t3) ), 1.0_wp )
!
!--                Calculate mean mass, mean diameter and mean terminal fall velocites for both
!--                graupel and snow
                   xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
                   xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )

                   dg = mean_diameter( cloud_species%graupel,  xg )
                   ds = mean_diameter( cloud_species%snow,     xs )

                   vg = terminal_fall_velocity(cloud_species%graupel,  xg, hyrho(k) )
                   vs = terminal_fall_velocity(cloud_species%snow,     xs, hyrho(k) )
!
!--                Calculate collection rate for number concentration
                   coll_rate_n = pi / 4.0_wp * ng(k,j,i) * ns(k,j,i) * e_coll *                    &
                             ( coll_graupel_snow_delta_n_aa * dg**2 +                              &
                               coll_graupel_snow_delta_n_ab * dg * ds +                            &
                               coll_graupel_snow_delta_n_bb * ds**2                                &
                             ) *                                                                   &
                         SQRT( coll_graupel_snow_theta_n_aa * vg**2 -                              &
                               coll_graupel_snow_theta_n_ab * vg * vs +                            &
                               coll_graupel_snow_theta_n_bb * vs**2 +                              &
                               cloud_species%snow%sigma_v**2                                       &
                             )
!
!--                Calculate collection rate for mixing ratio
                   coll_rate_q = pi / 4.0_wp * ng(k,j,i) * qs(k,j,i) * e_coll *                    &
                             ( coll_graupel_snow_delta_q_aa * dg**2 +                              &
                               coll_graupel_snow_delta_q_ab * dg * ds +                            &
                               coll_graupel_snow_delta_q_bb * ds**2                                &
                             ) *                                                                   &
                         SQRT( coll_graupel_snow_theta_q_aa * vg**2 -                              &
                               coll_graupel_snow_theta_q_ab * vg * vs +                            &
                               coll_graupel_snow_theta_q_bb * vs**2 +                              &
                               cloud_species%snow%sigma_v**2                                       &
                             )
!
!--                Limit collection rates
                   coll_rate_n = MIN( ns(k,j,i) / dt_micro, coll_rate_n )
                   coll_rate_q = MIN( qs(k,j,i) / dt_micro, coll_rate_q )
!
!--                Add/substract rates from species
                   qg(k,j,i) = qg(k,j,i) + coll_rate_q * dt_micro * flag
                   qs(k,j,i) = qs(k,j,i) - coll_rate_q * dt_micro * flag
                   ns(k,j,i) = ns(k,j,i) - coll_rate_n * dt_micro * flag

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE collection_graupel_snow


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collection of type a + b -> a. Here: graupel + snow -> snow (Seifert and Beheng, 2006).
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE collection_graupel_snow_ij( i, j )

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  coll_rate_n   !< collection rate for number concentration
       REAL(wp) ::  coll_rate_q   !< collection factor for mixing ratio
       REAL(wp) ::  dg            !< diameter graupel
       REAL(wp) ::  ds            !< diameter snow
       REAL(wp) ::  e_coll        !< collection efficiency
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vg            !< terminal fall velocitiy graupel
       REAL(wp) ::  vs            !< terminal fall velocitiy snow
       REAL(wp) ::  xg            !< mean mass graupel
       REAL(wp) ::  xs            !< mean mass snow

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Apply collection only above threshold
          IF ( qg(k,j,i) > eps_sb_coll  .AND.  qs(k,j,i) > eps_sb_coll )  THEN
!
!--          Calculate temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Sticking efficiency of Lin (1983)
             e_coll = MIN( EXP( 0.09_wp * (temp - t3) ), 1.0_wp )
!
!--          Calculate mean mass, mean diameter and mean terminal fall velocites for both graupel
!--          and snow
             xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
             xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )

             dg = mean_diameter( cloud_species%graupel,  xg )
             ds = mean_diameter( cloud_species%snow,     xs )

             vg = terminal_fall_velocity(cloud_species%graupel,  xg, hyrho(k) )
             vs = terminal_fall_velocity(cloud_species%snow,     xs, hyrho(k) )
!
!--          Calculate collection rate for number concentration
             coll_rate_n = pi / 4.0_wp * ng(k,j,i) * ns(k,j,i) * e_coll *                              &
                       ( coll_graupel_snow_delta_n_aa * dg**2 +                                        &
                         coll_graupel_snow_delta_n_ab * dg * ds +                                      &
                         coll_graupel_snow_delta_n_bb * ds**2                                          &
                       ) *                                                                             &
                   SQRT( coll_graupel_snow_theta_n_aa * vg**2 -                                        &
                         coll_graupel_snow_theta_n_ab * vg * vs +                                      &
                         coll_graupel_snow_theta_n_bb * vs**2 +                                        &
                         cloud_species%snow%sigma_v**2                                                 &
                       )
!
!--          Calculate collection rate for mixing ratio
             coll_rate_q = pi / 4.0_wp * ng(k,j,i) * qs(k,j,i) * e_coll *                              &
                       ( coll_graupel_snow_delta_q_aa * dg**2 +                                        &
                         coll_graupel_snow_delta_q_ab * dg * ds +                                      &
                         coll_graupel_snow_delta_q_bb * ds**2                                          &
                       ) *                                                                             &
                   SQRT( coll_graupel_snow_theta_q_aa * vg**2 -                                        &
                         coll_graupel_snow_theta_q_ab * vg * vs +                                      &
                         coll_graupel_snow_theta_q_bb * vs**2 +                                        &
                         cloud_species%snow%sigma_v**2                                                 &
                       )
!
!--          Limit collection rates
             coll_rate_n = MIN( ns(k,j,i) / dt_micro, coll_rate_n )
             coll_rate_q = MIN( qs(k,j,i) / dt_micro, coll_rate_q )
!
!--          Add/substract rates from species
             qg(k,j,i) = qg(k,j,i) + coll_rate_q * dt_micro * flag
             qs(k,j,i) = qs(k,j,i) - coll_rate_q * dt_micro * flag
             ns(k,j,i) = ns(k,j,i) - coll_rate_n * dt_micro * flag

          ENDIF
       ENDDO

    END SUBROUTINE collection_graupel_snow_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of cloud and graupel (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_graupel_cloud

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  dc            !< diameter cloud
       REAL(wp) ::  dg            !< diameter graupel
       REAL(wp) ::  const1        !< constant for riming efficiency
       REAL(wp) ::  e_coll_n      !< collision efficences for number concentration
       REAL(wp) ::  e_coll_q      !< collision efficences for mixing ratio
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  melt_n        !< melting rate number concentration
       REAL(wp) ::  melt_q        !< melting rate mixing ratio
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  nc_c          !< dummy number concentration cloud
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_q !< riming rate mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vc            !< terminal fall velocitiy cloud
       REAL(wp) ::  vg            !< terminal fall velocitiy graupel
       REAL(wp) ::  xc            !< mean mass cloud
       REAL(wp) ::  xg            !< mean mass graupel

       REAL(wp), PARAMETER ::  const0 = 1.0 / ( dc_coll - dc_crit )             !< constant for riming
       REAL(wp), PARAMETER ::  const2 = 1.0 / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplaction
       REAL(wp), PARAMETER ::  const3 = 1.0 / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = c_w / l_m                               !< constant for melting

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                IF ( microphysics_morrison )  THEN
                   nc_c = nc(k,j,i)
                ELSE
                   nc_c = nc_const
                ENDIF
!
!--             Calculate mean mass, mean diameter for both graupel and cloud
                xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
                xc = mean_mass( cloud_species%cloud,   qc(k,j,i), nc_c, hyrho(k) )

                dg = mean_diameter( cloud_species%graupel,  xg )
                dc = mean_diameter( cloud_species%cloud,    xc )
!
!--             Apply riming only above threshold
                IF ( qc(k,j,i) > q_thres  .AND.  qg(k,j,i) > q_thres  .AND.                        &
                     dc > dc_crit  .AND.  dg > df_crit )  THEN

                   const1 = const0 * cloud_species%graupel%coll_eff
!
!--                Calculate terminal velocites
                   vg = terminal_fall_velocity(cloud_species%graupel,  xg, hyrho(k) )
                   vc = terminal_fall_velocity(cloud_species%cloud,    xc, hyrho(k) )
!
!--                Calcualte collision_effciens
                   e_coll_n = MIN( cloud_species%graupel%coll_eff, MAX( const1 * ( dc - dc_crit ), &
                                                                        ecoll_min ) )
                   e_coll_q = e_coll_n
!
!--                Calculate riming rates
                   riming_rate_n = pi / 4.0_wp * e_coll_n * ng(k,j,i) * nc_c *                     &
                            ( rime_graupel_cloud_delta_n_aa * dg**2 +                              &
                              rime_graupel_cloud_delta_n_ab * dg * dc +                            &
                              rime_graupel_cloud_delta_n_bb * dc**2 ) *                            &
                        SQRT( rime_graupel_cloud_theta_n_aa * vg**2 -                              &
                              rime_graupel_cloud_theta_n_ab * vg * vc +                            &
                              rime_graupel_cloud_theta_n_bb * vc**2 )

                   riming_rate_q = pi / 4.0_wp * e_coll_q * ng(k,j,i) * qc(k,j,i) *                &
                            ( rime_graupel_cloud_delta_q_aa * dg**2 +                              &
                              rime_graupel_cloud_delta_q_ab * dg * dc +                            &
                              rime_graupel_cloud_delta_q_bb * dc**2) *                             &
                        SQRT( rime_graupel_cloud_theta_q_aa * vg**2 -                              &
                              rime_graupel_cloud_theta_q_ab * vg * vc +                            &
                              rime_graupel_cloud_theta_q_bb * vc**2 )
!
!--                Limit riming rates
                   riming_rate_q = MIN( qc(k,j,i) / dt_micro, riming_rate_q )
                   riming_rate_n = MIN( nc_c / dt_micro, riming_rate_n )
!
!--                Add/substract riming rates
                   qg(k,j,i) = qg(k,j,i) + riming_rate_q * dt_micro * flag
                   qc(k,j,i) = qc(k,j,i) - riming_rate_q * dt_micro * flag
                   IF ( microphysics_morrison )  THEN
                       nc(k,j,i) = nc(k,j,i) - riming_rate_n * dt_micro * flag
                   ENDIF
!
!--                Calcualte temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Ice multiplication based on Hallet and Mossop
                   IF ( temp < t3  .AND.  ice_multiplication )  THEN
                      mult_1 = const2 * ( temp - temp_mult_min )
                      mult_2 = const3 * ( temp - temp_mult_max )
                      mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                      mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                      mult_n = c_mult * mult_1 * mult_2 * riming_rate_q
                      mult_q = mult_n * cloud_species%ice%x_min
                      mult_q = MIN( riming_rate_q, mult_q )
!
!--                   Apply Hallet Mossop process
                      ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                      qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
                      qg(k,j,i) = qg(k,j,i) - mult_q * dt_micro * flag
                   ENDIF
!
!--                Enhancement of melting
                   IF ( temp > t3  .AND.  enhanced_melting )  THEN
                      melt_q = const4 * ( temp - t3 ) * riming_rate_q
                      melt_n = melt_q / xg
                      melt_q = MIN( qg(k,j,i) / dt_micro, melt_q )
                      melt_n = MIN( ng(k,j,i) / dt_micro, melt_n )
!
!--                   Apply enhancement due to melting
                      qg(k,j,i) = qg(k,j,i) - melt_q * dt_micro * flag
                      qr(k,j,i) = qr(k,j,i) + melt_q * dt_micro * flag
                      ng(k,j,i) = ng(k,j,i) - melt_n * dt_micro * flag
                      nr(k,j,i) = nr(k,j,i) + melt_n * dt_micro * flag
                   ENDIF

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE riming_graupel_cloud


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of cloud and graupel (Seifert and Beheng, 2006).
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_graupel_cloud_ij( i, j )

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  dc            !< diameter cloud
       REAL(wp) ::  dg            !< diameter graupel
       REAL(wp) ::  const1        !< constant for riming efficiency
       REAL(wp) ::  e_coll_n      !< collision efficences for number concentration
       REAL(wp) ::  e_coll_q      !< collision efficences for mixing ratio
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  melt_n        !< melting rate number concentration
       REAL(wp) ::  melt_q        !< melting rate mixing ratio
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  nc_c          !< dummy number concentration cloud
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_q !< riming rate mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vc            !< terminal fall velocitiy cloud
       REAL(wp) ::  vg            !< terminal fall velocitiy graupel
       REAL(wp) ::  xc            !< mean mass cloud
       REAL(wp) ::  xg            !< mean mass graupel

       REAL(wp), PARAMETER ::  const0 = 1.0 / ( dc_coll - dc_crit )             !< constant for riming
       REAL(wp), PARAMETER ::  const2 = 1.0 / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplaction
       REAL(wp), PARAMETER ::  const3 = 1.0 / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = c_w / l_m                               !< constant for melting

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
          IF ( microphysics_morrison )  THEN
             nc_c = nc(k,j,i)
          ELSE
             nc_c = nc_const
          ENDIF
!
!--       Calculate mean mass, mean diameter for both graupel and cloud
          xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
          xc = mean_mass( cloud_species%cloud,   qc(k,j,i), nc_c, hyrho(k) )

          dg = mean_diameter( cloud_species%graupel,  xg )
          dc = mean_diameter( cloud_species%cloud,    xc )
!
!--       Apply riming only above threshold
          IF ( qc(k,j,i) > q_thres  .AND.  qg(k,j,i) > q_thres  .AND.                              &
               dc > dc_crit  .AND.  dg > df_crit )  THEN

             const1 = const0 * cloud_species%graupel%coll_eff
!
!--          Calculate terminal velocites
             vg = terminal_fall_velocity(cloud_species%graupel,  xg, hyrho(k) )
             vc = terminal_fall_velocity(cloud_species%cloud,    xc, hyrho(k) )
!
!--          Calcualte collision_effciens
             e_coll_n = MIN( cloud_species%graupel%coll_eff, MAX( const1 * ( dc - dc_crit ),       &
                                                                  ecoll_min ) )
             e_coll_q = e_coll_n
!
!--          Calculate riming rates
             riming_rate_n = pi / 4.0_wp * e_coll_n * ng(k,j,i) * nc_c *                           &
                      ( rime_graupel_cloud_delta_n_aa * dg**2 +                                    &
                        rime_graupel_cloud_delta_n_ab * dg * dc +                                  &
                        rime_graupel_cloud_delta_n_bb * dc**2 ) *                                  &
                  SQRT( rime_graupel_cloud_theta_n_aa * vg**2 -                                    &
                        rime_graupel_cloud_theta_n_ab * vg * vc +                                  &
                        rime_graupel_cloud_theta_n_bb * vc**2 )

             riming_rate_q = pi / 4.0_wp * e_coll_q * ng(k,j,i) * qc(k,j,i) *                      &
                      ( rime_graupel_cloud_delta_q_aa * dg**2 +                                    &
                        rime_graupel_cloud_delta_q_ab * dg * dc +                                  &
                        rime_graupel_cloud_delta_q_bb * dc**2) *                                   &
                  SQRT( rime_graupel_cloud_theta_q_aa * vg**2 -                                    &
                        rime_graupel_cloud_theta_q_ab * vg * vc +                                  &
                        rime_graupel_cloud_theta_q_bb * vc**2 )
!
!--          Limit riming rates
             riming_rate_q = MIN( qc(k,j,i) / dt_micro, riming_rate_q )
             riming_rate_n = MIN( nc_c / dt_micro, riming_rate_n )
!
!--          Add/substract riming rates
             qg(k,j,i) = qg(k,j,i) + riming_rate_q * dt_micro * flag
             qc(k,j,i) = qc(k,j,i) - riming_rate_q * dt_micro * flag
             IF ( microphysics_morrison )  THEN
                 nc(k,j,i) = nc(k,j,i) - riming_rate_n * dt_micro * flag
             ENDIF
!
!--          Calcualte temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Ice multiplication based on Hallet and Mossop
             IF ( temp < t3  .AND.  ice_multiplication )  THEN
                mult_1 = const2 * ( temp - temp_mult_min )
                mult_2 = const3 * ( temp - temp_mult_max )
                mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                mult_n = c_mult * mult_1 * mult_2 * riming_rate_q
                mult_q = mult_n * cloud_species%ice%x_min
                mult_q = MIN( riming_rate_q, mult_q )
!
!--             Apply Hallet Mossop process
                ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
                qg(k,j,i) = qg(k,j,i) - mult_q * dt_micro * flag
             ENDIF
!
!--          Enhancement of melting
             IF ( temp > t3  .AND.  enhanced_melting )  THEN
                melt_q = const4 * ( temp - t3 ) * riming_rate_q
                melt_n = melt_q / xg
                melt_q = MIN( qg(k,j,i) / dt_micro, melt_q )
                melt_n = MIN( ng(k,j,i) / dt_micro, melt_n )
!
!--             Apply enhancement due to melting
                qg(k,j,i) = qg(k,j,i) - melt_q * dt_micro * flag
                qr(k,j,i) = qr(k,j,i) + melt_q * dt_micro * flag
                ng(k,j,i) = ng(k,j,i) - melt_n * dt_micro * flag
                nr(k,j,i) = nr(k,j,i) + melt_n * dt_micro * flag
             ENDIF

          ENDIF
       ENDDO

    END SUBROUTINE riming_graupel_cloud_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of rain and graupel (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_graupel_rain

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  dg            !< diameter graupel
       REAL(wp) ::  dr            !< diameter rain
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  melt_n        !< melting rate number concentration
       REAL(wp) ::  melt_q        !< melting rate mixing ratio
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_q !< riming rate mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vg            !< terminal fall velocitiy graupel
       REAL(wp) ::  vr            !< terminal fall velocitiy rain
       REAL(wp) ::  xg            !< mean mass graupel
       REAL(wp) ::  xr            !< mean mass rain

       REAL(wp), PARAMETER ::  const2 = 1.0 / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const3 = 1.0 / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = c_w / l_m                               !< constant for melting

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Apply riming only above threshold
                IF ( qr(k,j,i) > eps_sb_coll  .AND.  qg(k,j,i) > eps_sb_coll                       &
                                              .AND.  ng(k,j,i) > 0.0_wp )  THEN
!
!--                Calculate mean mass, mean diameter for both graupel and rain
                   xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
                   xr = mean_mass( cloud_species%rain,    qr(k,j,i), nr(k,j,i), hyrho(k) )

                   dg = mean_diameter( cloud_species%graupel,  xg )
                   dr = mean_diameter( cloud_species%rain,     xr )
!
!--                Calculate terminal velocites
                   vg = terminal_fall_velocity(cloud_species%graupel,  xg, hyrho(k) )
                   vr = terminal_fall_velocity(cloud_species%rain,     xr, hyrho(k) )

!--                Calculate riming rates
                   riming_rate_n = pi / 4.0_wp *  ng(k,j,i) * nr(k,j,i) *                          &
                            ( rime_graupel_rain_delta_n_aa * dg**2 +                               &
                              rime_graupel_rain_delta_n_ab * dg * dr +                             &
                              rime_graupel_rain_delta_n_bb * dr**2 ) *                             &
                        SQRT( rime_graupel_rain_theta_n_aa * vg**2 -                               &
                              rime_graupel_rain_theta_n_ab * vg * vr +                             &
                              rime_graupel_rain_theta_n_bb * vr**2 )

                   riming_rate_q = pi / 4.0_wp *  ng(k,j,i) * qr(k,j,i) *                          &
                            ( rime_graupel_rain_delta_q_aa * dg**2 +                               &
                              rime_graupel_rain_delta_q_ab * dg * dr +                             &
                              rime_graupel_rain_delta_q_bb * dr**2) *                              &
                        SQRT( rime_graupel_rain_theta_q_aa * vg**2 -                               &
                              rime_graupel_rain_theta_q_ab * vg * vr +                             &
                              rime_graupel_rain_theta_q_bb * vr**2 )

!
!--                Limit riming rates
                   riming_rate_q = MIN( qr(k,j,i) / dt_micro, MAX( riming_rate_q, 0.0_wp ) )
                   riming_rate_n = MIN( nr(k,j,i) / dt_micro, MAX( riming_rate_n, 0.0_wp ) )
!
!--                Add/substract riming rates
                   qg(k,j,i) = qg(k,j,i) + riming_rate_q * dt_micro * flag
                   qr(k,j,i) = qr(k,j,i) - riming_rate_q * dt_micro * flag
                   nr(k,j,i) = nr(k,j,i) - riming_rate_n * dt_micro * flag
!
!--                Calcualte temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Ice multiplication based on Hallet and Mossop
                   IF ( temp < t3  .AND.  ice_multiplication )  THEN
                      mult_1 = const2 * ( temp - temp_mult_min )
                      mult_2 = const3 * ( temp - temp_mult_max )
                      mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                      mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                      mult_n = c_mult * mult_1 * mult_2 * riming_rate_q
                      mult_q = mult_n * cloud_species%ice%x_min
                      mult_q = MIN( riming_rate_q, mult_q )
!
!--                   Apply Hallet Mossop process
                      ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                      qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
                      qg(k,j,i) = qg(k,j,i) - mult_q * dt_micro * flag
                   ENDIF
!
!--                Enhancement of melting
                   IF ( temp > t3  .AND.  enhanced_melting )  THEN
                      melt_q = const4 * ( temp - t3 ) * riming_rate_q
                      melt_n = melt_q / xg
                      melt_q = MIN( qg(k,j,i) / dt_micro, melt_q )
                      melt_n = MIN( ng(k,j,i) / dt_micro, melt_n )
!
!--                   Apply enhancement due to melting
                      qg(k,j,i) = qg(k,j,i) - melt_q * dt_micro * flag
                      qr(k,j,i) = qr(k,j,i) + melt_q * dt_micro * flag
                      ng(k,j,i) = ng(k,j,i) - melt_n * dt_micro * flag
                      nr(k,j,i) = nr(k,j,i) + melt_n * dt_micro * flag
                   ENDIF

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE riming_graupel_rain


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of rain and graupel (Seifert and Beheng, 2006).
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_graupel_rain_ij( i, j )

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  dg            !< diameter graupel
       REAL(wp) ::  dr            !< diameter rain
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  melt_n        !< melting rate number concentration
       REAL(wp) ::  melt_q        !< melting rate mixing ratio
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_q !< riming rate mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vg            !< terminal fall velocitiy graupel
       REAL(wp) ::  vr            !< terminal fall velocitiy rain
       REAL(wp) ::  xg            !< mean mass graupel
       REAL(wp) ::  xr            !< mean mass rain

       REAL(wp), PARAMETER ::  const2 = 1.0 / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const3 = 1.0 / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = c_w / l_m                               !< constant for melting

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Apply riming only above threshold
          IF ( qr(k,j,i) > eps_sb_coll  .AND.  qg(k,j,i) > eps_sb_coll                             &
                                        .AND.  ng(k,j,i) > 0.0_wp )  THEN
!
!--          Calculate mean mass, mean diameter for both graupel and rain
             xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
             xr = mean_mass( cloud_species%rain,    qr(k,j,i), nr(k,j,i), hyrho(k) )

             dg = mean_diameter( cloud_species%graupel,  xg )
             dr = mean_diameter( cloud_species%rain,     xr )
!
!--          Calculate terminal velocites
             vg = terminal_fall_velocity(cloud_species%graupel,  xg, hyrho(k) )
             vr = terminal_fall_velocity(cloud_species%rain,     xr, hyrho(k) )

!--          Calculate riming rates
             riming_rate_n = pi / 4.0_wp *  ng(k,j,i) * nr(k,j,i) *                                &
                      ( rime_graupel_rain_delta_n_aa * dg**2 +                                     &
                        rime_graupel_rain_delta_n_ab * dg * dr +                                   &
                        rime_graupel_rain_delta_n_bb * dr**2 ) *                                   &
                  SQRT( rime_graupel_rain_theta_n_aa * vg**2 -                                     &
                        rime_graupel_rain_theta_n_ab * vg * vr +                                   &
                        rime_graupel_rain_theta_n_bb * vr**2 )

             riming_rate_q = pi / 4.0_wp *  ng(k,j,i) * qr(k,j,i) *                                &
                      ( rime_graupel_rain_delta_q_aa * dg**2 +                                     &
                        rime_graupel_rain_delta_q_ab * dg * dr +                                   &
                        rime_graupel_rain_delta_q_bb * dr**2) *                                    &
                  SQRT( rime_graupel_rain_theta_q_aa * vg**2 -                                     &
                        rime_graupel_rain_theta_q_ab * vg * vr +                                   &
                        rime_graupel_rain_theta_q_bb * vr**2 )

!
!--          Limit riming rates
             riming_rate_q = MIN( qr(k,j,i) / dt_micro, MAX( riming_rate_q, 0.0_wp ) )
             riming_rate_n = MIN( nr(k,j,i) / dt_micro, MAX( riming_rate_n, 0.0_wp ) )
!
!--          Add/substract riming rates
             qg(k,j,i) = qg(k,j,i) + riming_rate_q * dt_micro * flag
             qr(k,j,i) = qr(k,j,i) - riming_rate_q * dt_micro * flag
             nr(k,j,i) = nr(k,j,i) - riming_rate_n * dt_micro * flag
!
!--          Calcualte temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Ice multiplication based on Hallet and Mossop
             IF ( temp < t3  .AND.  ice_multiplication )  THEN
                mult_1 = const2 * ( temp - temp_mult_min )
                mult_2 = const3 * ( temp - temp_mult_max )
                mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                mult_n = c_mult * mult_1 * mult_2 * riming_rate_q
                mult_q = mult_n * cloud_species%ice%x_min
                mult_q = MIN( riming_rate_q, mult_q )
!
!--             Apply Hallet Mossop process
                ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
                qg(k,j,i) = qg(k,j,i) - mult_q * dt_micro * flag
             ENDIF
!
!--          Enhancement of melting
             IF ( temp > t3  .AND.  enhanced_melting )  THEN
                melt_q = const4 * ( temp - t3 ) * riming_rate_q
                melt_n = melt_q / xg
                melt_q = MIN( qg(k,j,i) / dt_micro, melt_q )
                melt_n = MIN( ng(k,j,i) / dt_micro, melt_n )
!
!--             Apply enhancement due to melting
                qg(k,j,i) = qg(k,j,i) - melt_q * dt_micro * flag
                qr(k,j,i) = qr(k,j,i) + melt_q * dt_micro * flag
                ng(k,j,i) = ng(k,j,i) - melt_n * dt_micro * flag
                nr(k,j,i) = nr(k,j,i) + melt_n * dt_micro * flag
             ENDIF

          ENDIF
       ENDDO

    END SUBROUTINE riming_graupel_rain_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of cloud and ice (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_ice_cloud

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wP) ::  coll_eff      !< collision efficiency
       REAL(wp) ::  dc            !< diameter cloud
       REAL(wp) ::  di            !< diameter ice
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  nc_c          !< dummy number concentration cloud droplets
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_q !< riming rate mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vc            !< terminal fall velocitiy cloud
       REAL(wp) ::  vi            !< terminal fall velocitiy ice
       REAL(wp) ::  xc            !< mean mass cloud
       REAL(wp) ::  xi            !< mean mass ice
       REAL(wp) ::  conv_n        !< conversion rate ice graupel number concentration
       REAL(wp) ::  conv_q        !< conversion rate ice graupel mixing ratio
       REAL(wp) ::  const1        !< constant for riming

       REAL(wp), PARAMETER ::  const0 = 1.0_wp / ( dc_coll - dc_crit )  !< constant for riming calculation
       REAL(wp), PARAMETER ::  const3 = 1.0_wp / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = 1.0_wp / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const5 = alpha_spacefilling * rho_l / rho_i

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                IF ( microphysics_morrison )  THEN
                   nc_c = nc(k,j,i)
                ELSE
                   nc_c = nc_const
                ENDIF
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Calculate mean mass, mean diameter for both graupel and rain
                xc = mean_mass( cloud_species%cloud, qc(k,j,i), nc_c, hyrho(k) )
                xi = mean_mass( cloud_species%ice,   qi(k,j,i), ni(k,j,i), hyrho(k) )

                dc = mean_diameter( cloud_species%cloud,  xc )
                di = mean_diameter( cloud_species%ice,    xi )
!
!--             Apply riming only above threshold
                IF ( qc(k,j,i) > q_thres  .AND.  qi(k,j,i) > eps_sb_coll  .AND.                    &
                     dc > dc_crit  .AND.  di > df_crit )  THEN
!
!--                Calcualte temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Calculate terminal velocites
                   vc = terminal_fall_velocity(cloud_species%cloud,   xc, hyrho(k) )
                   vi = terminal_fall_velocity(cloud_species%ice,     xi, hyrho(k) )
!
!--                Calculate collision efficiency
                   const1 = const0 * cloud_species%ice%coll_eff
                   coll_eff = MIN( cloud_species%ice%coll_eff,                                     &
                                   MAX( const1 * ( dc - dc_crit ), ecoll_min) )
!
!--                Calcualte riming rates
                   riming_rate_n = pi / 4.0_wp * coll_eff * ni(k,j,i) * nc_c *                     &
                                ( rime_ice_cloud_delta_n_aa * di**2   +                            &
                                  rime_ice_cloud_delta_n_ab * di * dc +                            &
                                  rime_ice_cloud_delta_n_bb * dc**2)  *                            &
                            SQRT( rime_ice_cloud_theta_n_aa * vi**2   -                            &
                                  rime_ice_cloud_theta_n_ab * vi * vc +                            &
                                  rime_ice_cloud_theta_n_bb * vc**2   +                            &
                                  cloud_species%ice%sigma_v**2 )

                   riming_rate_q = pi / 4.0_wp * coll_eff * ni(k,j,i) * qc(k,j,i) *                &
                                ( rime_ice_cloud_delta_q_aa * di**2 +                              &
                                  rime_ice_cloud_delta_q_ab * di * dc +                            &
                                  rime_ice_cloud_delta_q_bb * dc**2) *                             &
                            SQRT( rime_ice_cloud_theta_q_aa * vi**2 -                              &
                                  rime_ice_cloud_theta_q_ab * vi * vc +                            &
                                  rime_ice_cloud_theta_q_bb * vc**2 +                              &
                                  cloud_species%ice%sigma_v**2 )
!
!--                Limit riming rates
                   riming_rate_q = MIN( qc(k,j,i) / dt_micro, MAX( riming_rate_q, 0.0_wp ) )
                   riming_rate_n = MIN( nc_c / dt_micro, MAX( riming_rate_n, 0.0_wp ) )
!
!--                Apply riming ice - cloud
                   qi(k,j,i) = qi(k,j,i) + riming_rate_q * dt_micro * flag
                   qc(k,j,i) = qc(k,j,i) - riming_rate_q * dt_micro * flag
                   IF ( microphysics_morrison )  THEN
                      nc(k,j,i) = nc(k,j,i) - riming_rate_n * dt_micro * flag
                   ENDIF
!
!--                Store riming rate for rain riming
                   rime_ice_cloud(k,j,i) = riming_rate_q
!
!--                Ice multiplication during riming
                   IF ( temp < t3 .AND. ice_multiplication )  THEN
                      mult_1 = ( temp - temp_mult_min ) * const3
                      mult_2 = ( temp - temp_mult_max ) * const4
                      mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                      mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                      mult_n = c_mult * mult_1 * mult_2 * riming_rate_q

                      ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                   ENDIF
!
!--                If deposition is negative or deposition is smaller than riming, conversion from
!--                ice to graupel takes place. Note, that here the riming rate of ice-rain is from
!--                the last time-step. However, the error is assumed to be small
                   IF ( dep_rate_ice(k,j,i) < 0.0_wp  .OR.                                         &
                        dep_rate_ice(k,j,i) < rime_ice_cloud(k,j,i) + rime_ice_rain(k,j,i) )  THEN
!
!--                    Conversion from ice to graupel if diameter small enough and riming rates are
!--                    larger than deposition. If deposition is larger ice stay
!--                    factor alpha_spacefilling after SB 2006 and icon model
                       IF ( di > dig_conv ) THEN
                         conv_q =  riming_rate_q /                                                 &
                                  ( const5 * ( pi / 6.0_wp * rho_i * di**3 / xi - 1.0_wp ) )
                         conv_q = MIN( qi(k,j,i) / dt_micro , conv_q )
                         conv_n = conv_q / MAX( xi, x_conv )
                         conv_n = MIN( ni(k,j,i) / dt_micro, conv_n )
                      ELSE
                         conv_q = 0.0_wp
                         conv_n = 0.0_wp
                      ENDIF
!
!--                   Apply conversion from ice to grauppel
                      qi(k,j,i) = qi(k,j,i) - conv_q * dt_micro * flag
                      qg(k,j,i) = qg(k,j,i) + conv_q * dt_micro * flag
                      ni(k,j,i) = ni(k,j,i) - conv_n * dt_micro * flag
                      ng(k,j,i) = ng(k,j,i) + conv_n * dt_micro * flag
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE riming_ice_cloud


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of cloud and ice (Seifert and Beheng, 2006).
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_ice_cloud_ij( i, j )

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wP) ::  coll_eff      !< collision efficiency
       REAL(wp) ::  dc            !< diameter cloud
       REAL(wp) ::  di            !< diameter ice
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  nc_c          !< dummy number concentration cloud droplets
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_q !< riming rate mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vc            !< terminal fall velocitiy cloud
       REAL(wp) ::  vi            !< terminal fall velocitiy ice
       REAL(wp) ::  xc            !< mean mass cloud
       REAL(wp) ::  xi            !< mean mass ice
       REAL(wp) ::  conv_n        !< conversion rate ice graupel number concentration
       REAL(wp) ::  conv_q        !< conversion rate ice graupel mixing ratio
       REAL(wp) ::  const1        !< constant for riming

       REAL(wp), PARAMETER ::  const0 = 1.0_wp / ( dc_coll - dc_crit )  !< constant for riming calculation
       REAL(wp), PARAMETER ::  const3 = 1.0_wp / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = 1.0_wp / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const5 = alpha_spacefilling * rho_l / rho_i

       DO  k = nzb+1, nzt
          IF ( microphysics_morrison )  THEN
             nc_c = nc(k,j,i)
          ELSE
             nc_c = nc_const
          ENDIF
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Calculate mean mass, mean diameter for both graupel and rain
          xc = mean_mass( cloud_species%cloud, qc(k,j,i), nc_c, hyrho(k) )
          xi = mean_mass( cloud_species%ice,   qi(k,j,i), ni(k,j,i), hyrho(k) )

          dc = mean_diameter( cloud_species%cloud,  xc )
          di = mean_diameter( cloud_species%ice,    xi )
!
!--       Apply riming only above threshold
          IF ( qc(k,j,i) > q_thres  .AND.  qi(k,j,i) > eps_sb_coll  .AND.                          &
               dc > dc_crit  .AND.  di > df_crit )  THEN
!
!--          Calcualte temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Calculate terminal velocites
             vc = terminal_fall_velocity(cloud_species%cloud,   xc, hyrho(k) )
             vi = terminal_fall_velocity(cloud_species%ice,     xi, hyrho(k) )
!
!--          Calculate collision efficiency
             const1 = const0 * cloud_species%ice%coll_eff
             coll_eff = MIN( cloud_species%ice%coll_eff,                                           &
                             MAX( const1 * ( dc - dc_crit ), ecoll_min) )
!
!--          Calcualte riming rates
             riming_rate_n = pi / 4.0_wp * coll_eff * ni(k,j,i) * nc_c *                           &
                          ( rime_ice_cloud_delta_n_aa * di**2   +                                  &
                            rime_ice_cloud_delta_n_ab * di * dc +                                  &
                            rime_ice_cloud_delta_n_bb * dc**2)  *                                  &
                      SQRT( rime_ice_cloud_theta_n_aa * vi**2   -                                  &
                            rime_ice_cloud_theta_n_ab * vi * vc +                                  &
                            rime_ice_cloud_theta_n_bb * vc**2   +                                  &
                            cloud_species%ice%sigma_v**2 )

             riming_rate_q = pi / 4.0_wp * coll_eff * ni(k,j,i) * qc(k,j,i) *                      &
                          ( rime_ice_cloud_delta_q_aa * di**2 +                                    &
                            rime_ice_cloud_delta_q_ab * di * dc +                                  &
                            rime_ice_cloud_delta_q_bb * dc**2) *                                   &
                      SQRT( rime_ice_cloud_theta_q_aa * vi**2 -                                    &
                            rime_ice_cloud_theta_q_ab * vi * vc +                                  &
                            rime_ice_cloud_theta_q_bb * vc**2 +                                    &
                            cloud_species%ice%sigma_v**2 )
!
!--          Limit riming rates
             riming_rate_q = MIN( qc(k,j,i) / dt_micro, MAX( riming_rate_q, 0.0_wp ) )
             riming_rate_n = MIN( nc_c / dt_micro, MAX( riming_rate_n, 0.0_wp ) )
!
!--          Apply riming ice - cloud
             qi(k,j,i) = qi(k,j,i) + riming_rate_q * dt_micro * flag
             qc(k,j,i) = qc(k,j,i) - riming_rate_q * dt_micro * flag
             IF ( microphysics_morrison )  THEN
                nc(k,j,i) = nc(k,j,i) - riming_rate_n * dt_micro * flag
             ENDIF
!
!--          Store riming rate for rain riming
             rime_ice_cloud(k,j,i) = riming_rate_q
!
!--          Ice multiplication during riming
             IF ( temp < t3 .AND. ice_multiplication )  THEN
                mult_1 = ( temp - temp_mult_min ) * const3
                mult_2 = ( temp - temp_mult_max ) * const4
                mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                mult_n = c_mult * mult_1 * mult_2 * riming_rate_q

                ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
             ENDIF
!
!--          If deposition is negative or deposition is smaller than riming, conversion from ice
!--          ice to graupel takes place. Note, that here the riming rate of ice-rain is from the
!--          last time-step. However, the error is assumed to be small
             IF ( dep_rate_ice(k,j,i) < 0.0_wp  .OR.                                               &
                  dep_rate_ice(k,j,i) < rime_ice_cloud(k,j,i) + rime_ice_rain(k,j,i) )  THEN
!
!--              Conversion from ice to graupel if diameter small enough and riming rates are larger
!--              than deposition. If deposition is larger ice stay
!--              factor alpha_spacefilling after SB 2006 and icon model
                 IF ( di > dig_conv ) THEN
                   conv_q =  riming_rate_q /                                                       &
                            ( const5 * ( pi / 6.0_wp * rho_i * di**3 / xi - 1.0_wp ) )
                   conv_q = MIN( qi(k,j,i) / dt_micro , conv_q )
                   conv_n = conv_q / MAX( xi, x_conv )
                   conv_n = MIN( ni(k,j,i) / dt_micro, conv_n )
                ELSE
                   conv_q = 0.0_wp
                   conv_n = 0.0_wp
                ENDIF
!
!--             Apply conversion from ice to grauppel
                qi(k,j,i) = qi(k,j,i) - conv_q * dt_micro * flag
                qg(k,j,i) = qg(k,j,i) + conv_q * dt_micro * flag
                ni(k,j,i) = ni(k,j,i) - conv_n * dt_micro * flag
                ng(k,j,i) = ng(k,j,i) + conv_n * dt_micro * flag
             ENDIF
          ENDIF

       ENDDO

    END SUBROUTINE riming_ice_cloud_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of rain and ice (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_ice_rain

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  dr            !< diameter rain
       REAL(wp) ::  di            !< diameter ice
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_qi!< riming rate ice mixing ratio
       REAL(wp) ::  riming_rate_qr!< riming rate rain mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vr            !< terminal fall velocitiy rain
       REAL(wp) ::  vi            !< terminal fall velocitiy ice
       REAL(wp) ::  xr            !< mean mass rain
       REAL(wp) ::  xi            !< mean mass ice

       REAL(wp), PARAMETER ::  const3 = 1.0_wp / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = 1.0_wp / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Calculate mean mass, mean diameter for both graupel and rain
                xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
                xi = mean_mass( cloud_species%ice,  qi(k,j,i), ni(k,j,i), hyrho(k) )

                dr = mean_diameter( cloud_species%rain,  xr )
                di = mean_diameter( cloud_species%ice,   xi )
!
!--             Apply ice-rain riming only above thresholds
                IF ( qr(k,j,i) > eps_sb_coll  .AND.  qi(k,j,i) > qr_crit  .AND.  di > df_crit      &
                     .AND.  nr(k,j,i) > 0.0_wp  .AND.  ni(k,j,i) > 0.0_wp  )  THEN
!
!--                Calcualte temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Calculate terminal velocites
                   vr = terminal_fall_velocity(cloud_species%rain, xr, hyrho(k) )
                   vi = terminal_fall_velocity(cloud_species%ice,  xi, hyrho(k) )
!
!--                Calcualte riming rates
                   riming_rate_n = pi / 4.0_wp * ni(k,j,i) * nr(k,j,i) *                           &
                                ( rime_ice_rain_delta_n_aa * di**2   +                             &
                                  rime_ice_rain_delta_n_ab * di * dr +                             &
                                  rime_ice_rain_delta_n_bb * dr**2 ) *                             &
                            SQRT( rime_ice_rain_theta_n_aa * vi**2   -                             &
                                  rime_ice_rain_theta_n_ab * vi * vr +                             &
                                  rime_ice_rain_theta_n_bb * vr**2   +                             &
                                  cloud_species%ice%sigma_v**2 )

                   riming_rate_qr = pi / 4.0_wp * ni(k,j,i) * qr(k,j,i) *                          &
                                ( rime_ice_rain_delta_n_aa * di**2 +                               &
                                  rime_ice_rain_delta_q_ab * di * dr +                             &
                                  rime_ice_rain_delta_q_bb * dr**2 ) *                             &
                            SQRT( rime_ice_rain_theta_n_aa * vi**2 -                               &
                                  rime_ice_rain_theta_q_ab * vi * vr +                             &
                                  rime_ice_rain_theta_q_bb * vr**2 +                               &
                                  cloud_species%ice%sigma_v**2 )

                   riming_rate_qi = pi / 4.0_wp * nr(k,j,i) * qi(k,j,i) *                          &
                                ( rime_ice_rain_delta_q_aa * di**2 +                               &
                                  rime_ice_rain_delta_q_ab * di * dr +                             &
                                  rime_ice_rain_delta_q_bb * dr**2 ) *                             &
                            SQRT( rime_ice_rain_theta_n_aa * vi**2 -                               &
                                  rime_ice_rain_theta_q_ab * vi * vr +                             &
                                  rime_ice_rain_theta_q_bb * vr**2 +                               &
                                  cloud_species%ice%sigma_v**2 )

!Q TODO hi      er muss noch ein fehler behoben werden, evtl wird die wurzel negativ
!                    riming_rate_qi = pi / 4.0_wp * nr(k,j,i) * qi(k,j,i) *                                &
!                                 ( rime_ice_rain_delta_q_aa * di**2 +                                     &
!                                   rime_ice_rain_delta_q_ab * di * dr +                                   &
!                                   rime_ice_rain_delta_n_bb * dr**2 ) *                                   &
!                             SQRT( rime_ice_rain_theta_q_aa * vi**2 -                                     &
!                                   rime_ice_rain_theta_q_ab * vi * vr +                                   &
!                                   rime_ice_rain_theta_n_bb * vr**2 +                                     &
!                                   cloud_species%ice%sigma_v**2 )
!
!                 rime_qi = pi4 * n_r * q_a * dt &
!                      &  *     (  coeffs%delta_q_aa * D_a * D_a &
!                      &         + coeffs%delta_q_ba * D_a * D_r &
!                      &         + coeffs%delta_n_bb * D_r * D_r) &
!                      &  * SQRT(  coeffs%theta_q_aa * v_a * v_a &
!                      &         - coeffs%theta_q_ba * v_a * v_r &
!                      &         + coeffs%theta_n_bb * v_r * v_r &
!                      &         + ptype%s_vel**2)

!--                Limit riming rates
                   riming_rate_qr = MIN( qr(k,j,i) / dt_micro, MAX( riming_rate_qr, 0.0_wp ) )
                   riming_rate_qi = MIN( qi(k,j,i) / dt_micro, MAX( riming_rate_qi, 0.0_wp ) )
                   riming_rate_n  = MIN( nr(k,j,i) / dt_micro, MAX( riming_rate_n,  0.0_wp ) )
!
!--                store ice-rain riming rate for ice-cloud riming
                   rime_ice_rain(k,j,i) = riming_rate_qr
!
!--                If deposition is larger than riming ice stays ice
                   IF ( dep_rate_ice(k,j,i) > 0.0_wp  .AND.                                        &
                        dep_rate_ice(k,j,i) < rime_ice_rain(k,j,i) + rime_ice_cloud(k,j,i) ) THEN
!
!--                   Apply riming ice - rain
                      qi(k,j,i) = qi(k,j,i) + riming_rate_qr * dt_micro * flag
                      qr(k,j,i) = qr(k,j,i) - riming_rate_qr * dt_micro * flag
                      nr(k,j,i) = nr(k,j,i) - riming_rate_n  * dt_micro * flag
!
!--                   Ice multiplication during riming
                      IF ( temp < t3 .AND. ice_multiplication )  THEN
                         mult_1 = ( temp - temp_mult_min ) * const3
                         mult_2 = ( temp - temp_mult_max ) * const4
                         mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                         mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                         mult_n = c_mult * mult_1 * mult_2 * riming_rate_qr

                         ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                      ENDIF
!
!--                If deposition is negative (i.e. sublimation or riming is stronger, a different
!--                treatment is conducted: ice might convert to graupel
                   ELSE
!
!--                   Ensure that riming rate is also not larger than ni
                      riming_rate_n  = MIN( MIN( nr(k,j,i) / dt_micro, ni(k,j,i) / dt_micro ),     &
                                             MAX( riming_rate_n,  0.0_wp ) )
!
!--                   Substract all riming rates, those riming rates are added to species (either
!--                   ice again or graupel furhter down)
                      qi(k,j,i) = qi(k,j,i) - riming_rate_qi * dt_micro * flag
                      qr(k,j,i) = qr(k,j,i) - riming_rate_qr * dt_micro * flag
                      ni(k,j,i) = ni(k,j,i) - riming_rate_n  * dt_micro * flag
                      nr(k,j,i) = nr(k,j,i) - riming_rate_n  * dt_micro * flag
!
!--                   Ice multiplication
                      mult_q = 0.0_wp
                      mult_n = 0.0_wp
                      IF ( temp < t3  .AND.  ice_multiplication )  THEN
                         mult_1 = ( temp - temp_mult_min ) * const3
                         mult_2 = ( temp - temp_mult_max ) * const4
                         mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp) )
                         mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp) )
                         mult_n = c_mult * mult_1 * mult_2 * riming_rate_qr
                         mult_q = mult_n * cloud_species%ice%x_min
                         mult_q = MIN( riming_rate_qr, mult_q )
                      ENDIF
!
!--                   shedding of rain at warm temperatures
                      IF ( temp >= t3 )  THEN
!
!--                      but with modified nr
                         xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
                         nr(k,j,i) = nr(k,j,i) + riming_rate_qr / xr * dt_micro * flag
                         ni(k,j,i) = ni(k,j,i) + riming_rate_n  * dt_micro * flag
                         qi(k,j,i) = qi(k,j,i) + riming_rate_qi * dt_micro * flag
                         qr(k,j,i) = qr(k,j,i) + riming_rate_qr * dt_micro * flag
!
!--                   at tempatures below freezing point graupel is generated
                      ELSE
!
!--                      new ice particles from multiplication
                         ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                         qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
!
!--                      and riming to graupel
                         ng(k,j,i) = ng(k,j,i) + riming_rate_n * dt_micro * flag
                         qg(k,j,i) = qg(k,j,i) + ( riming_rate_qi + riming_rate_qr - mult_q ) *    &
                                                                                     dt_micro * flag
                      END IF

                   ENDIF ! loop if deposition or riming is stronger
                ENDIF ! loop if thresholds for riming are exceeded, otherwise nothings happen

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE riming_ice_rain


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of rain and ice (Seifert and Beheng, 2006).
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_ice_rain_ij( i, j )

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  dr            !< diameter rain
       REAL(wp) ::  di            !< diameter ice
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_qi!< riming rate ice mixing ratio
       REAL(wp) ::  riming_rate_qr!< riming rate rain mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vr            !< terminal fall velocitiy rain
       REAL(wp) ::  vi            !< terminal fall velocitiy ice
       REAL(wp) ::  xr            !< mean mass rain
       REAL(wp) ::  xi            !< mean mass ice

       REAL(wp), PARAMETER ::  const3 = 1.0_wp / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = 1.0_wp / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Calculate mean mass, mean diameter for both graupel and rain
          xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
          xi = mean_mass( cloud_species%ice,  qi(k,j,i), ni(k,j,i), hyrho(k) )

          dr = mean_diameter( cloud_species%rain,  xr )
          di = mean_diameter( cloud_species%ice,   xi )
!
!--       Apply ice-rain riming only above thresholds
          IF ( qr(k,j,i) > eps_sb_coll  .AND.  qi(k,j,i) > qr_crit  .AND.  di > df_crit            &
               .AND.  nr(k,j,i) > 0.0_wp  .AND.  ni(k,j,i) > 0.0_wp  )  THEN
!
!--          Calcualte temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Calculate terminal velocites
             vr = terminal_fall_velocity(cloud_species%rain, xr, hyrho(k) )
             vi = terminal_fall_velocity(cloud_species%ice,  xi, hyrho(k) )
!
!--          Calcualte riming rates
             riming_rate_n = pi / 4.0_wp * ni(k,j,i) * nr(k,j,i) *                                 &
                          ( rime_ice_rain_delta_n_aa * di**2   +                                   &
                            rime_ice_rain_delta_n_ab * di * dr +                                   &
                            rime_ice_rain_delta_n_bb * dr**2 ) *                                   &
                      SQRT( rime_ice_rain_theta_n_aa * vi**2   -                                   &
                            rime_ice_rain_theta_n_ab * vi * vr +                                   &
                            rime_ice_rain_theta_n_bb * vr**2   +                                   &
                            cloud_species%ice%sigma_v**2 )

             riming_rate_qr = pi / 4.0_wp * ni(k,j,i) * qr(k,j,i) *                                &
                          ( rime_ice_rain_delta_n_aa * di**2 +                                     &
                            rime_ice_rain_delta_q_ab * di * dr +                                   &
                            rime_ice_rain_delta_q_bb * dr**2 ) *                                   &
                      SQRT( rime_ice_rain_theta_n_aa * vi**2 -                                     &
                            rime_ice_rain_theta_q_ab * vi * vr +                                   &
                            rime_ice_rain_theta_q_bb * vr**2 +                                     &
                            cloud_species%ice%sigma_v**2 )

             riming_rate_qi = pi / 4.0_wp * nr(k,j,i) * qi(k,j,i) *                                &
                          ( rime_ice_rain_delta_q_aa * di**2 +                                     &
                            rime_ice_rain_delta_q_ab * di * dr +                                   &
                            rime_ice_rain_delta_q_bb * dr**2 ) *                                   &
                      SQRT( rime_ice_rain_theta_n_aa * vi**2 -                                     &
                            rime_ice_rain_theta_q_ab * vi * vr +                                   &
                            rime_ice_rain_theta_q_bb * vr**2 +                                     &
                            cloud_species%ice%sigma_v**2 )

!Q TODO hier muss noch ein fehler behoben werden, evtl wird die wurzel negativ
!              riming_rate_qi = pi / 4.0_wp * nr(k,j,i) * qi(k,j,i) *                                &
!                           ( rime_ice_rain_delta_q_aa * di**2 +                                     &
!                             rime_ice_rain_delta_q_ab * di * dr +                                   &
!                             rime_ice_rain_delta_n_bb * dr**2 ) *                                   &
!                       SQRT( rime_ice_rain_theta_q_aa * vi**2 -                                     &
!                             rime_ice_rain_theta_q_ab * vi * vr +                                   &
!                             rime_ice_rain_theta_n_bb * vr**2 +                                     &
!                             cloud_species%ice%sigma_v**2 )
!
!           rime_qi = pi4 * n_r * q_a * dt &
!                &  *     (  coeffs%delta_q_aa * D_a * D_a &
!                &         + coeffs%delta_q_ba * D_a * D_r &
!                &         + coeffs%delta_n_bb * D_r * D_r) &
!                &  * SQRT(  coeffs%theta_q_aa * v_a * v_a &
!                &         - coeffs%theta_q_ba * v_a * v_r &
!                &         + coeffs%theta_n_bb * v_r * v_r &
!                &         + ptype%s_vel**2)

!--          Limit riming rates
             riming_rate_qr = MIN( qr(k,j,i) / dt_micro, MAX( riming_rate_qr, 0.0_wp ) )
             riming_rate_qi = MIN( qi(k,j,i) / dt_micro, MAX( riming_rate_qi, 0.0_wp ) )
             riming_rate_n  = MIN( nr(k,j,i) / dt_micro, MAX( riming_rate_n,  0.0_wp ) )
!
!--          store ice-rain riming rate for ice-cloud riming
             rime_ice_rain(k,j,i) = riming_rate_qr
!
!--          If deposition is larger than riming ice stays ice
             IF ( dep_rate_ice(k,j,i) > 0.0_wp  .AND.                                             &
                  dep_rate_ice(k,j,i) < rime_ice_rain(k,j,i) + rime_ice_cloud(k,j,i) ) THEN
!
!--             Apply riming ice - rain
                qi(k,j,i) = qi(k,j,i) + riming_rate_qr * dt_micro * flag
                qr(k,j,i) = qr(k,j,i) - riming_rate_qr * dt_micro * flag
                nr(k,j,i) = nr(k,j,i) - riming_rate_n  * dt_micro * flag
!
!--             Ice multiplication during riming
                IF ( temp < t3 .AND. ice_multiplication )  THEN
                   mult_1 = ( temp - temp_mult_min ) * const3
                   mult_2 = ( temp - temp_mult_max ) * const4
                   mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                   mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                   mult_n = c_mult * mult_1 * mult_2 * riming_rate_qr

                   ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                ENDIF
!
!--          If deposition is negative (i.e. sublimation or riming is stronger, a different
!--          treatment is conducted: ice might convert to graupel
             ELSE
!
!--             Ensure that riming rate is also not larger than ni
                riming_rate_n  = MIN( MIN( nr(k,j,i) / dt_micro, ni(k,j,i) / dt_micro ),          &
                                       MAX( riming_rate_n,  0.0_wp ) )
!
!--             Substract all riming rates, those riming rates are added to species (either ice
!--             again or graupel furhter down)
                qi(k,j,i) = qi(k,j,i) - riming_rate_qi * dt_micro * flag
                qr(k,j,i) = qr(k,j,i) - riming_rate_qr * dt_micro * flag
                ni(k,j,i) = ni(k,j,i) - riming_rate_n  * dt_micro * flag
                nr(k,j,i) = nr(k,j,i) - riming_rate_n  * dt_micro * flag
!
!--             Ice multiplication
                mult_q = 0.0_wp
                mult_n = 0.0_wp
                IF ( temp < t3  .AND.  ice_multiplication )  THEN
                   mult_1 = ( temp - temp_mult_min ) * const3
                   mult_2 = ( temp - temp_mult_max ) * const4
                   mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp) )
                   mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp) )
                   mult_n = c_mult * mult_1 * mult_2 * riming_rate_qr
                   mult_q = mult_n * cloud_species%ice%x_min
                   mult_q = MIN( riming_rate_qr, mult_q )
                ENDIF
!
!--             shedding of rain at warm temperatures
                IF ( temp >= t3 )  THEN
!
!--                but with modified nr
                   xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
                   nr(k,j,i) = nr(k,j,i) + riming_rate_qr / xr * dt_micro * flag
                   ni(k,j,i) = ni(k,j,i) + riming_rate_n  * dt_micro * flag
                   qi(k,j,i) = qi(k,j,i) + riming_rate_qi * dt_micro * flag
                   qr(k,j,i) = qr(k,j,i) + riming_rate_qr * dt_micro * flag
!
!--             at tempatures below freezing point graupel is generated
                ELSE
!
!--                new ice particles from multiplication
                   ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                   qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
!
!--                and riming to graupel
                   ng(k,j,i) = ng(k,j,i) + riming_rate_n * dt_micro * flag
                   qg(k,j,i) = qg(k,j,i) + ( riming_rate_qi + riming_rate_qr - mult_q ) *          &
                                                                                     dt_micro * flag
                END IF

             ENDIF ! loop if deposition or riming is stronger
          ENDIF ! loop if thresholds for riming are exceeded, otherwise nothings happen

       ENDDO

    END SUBROUTINE riming_ice_rain_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of cloud and snow (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_snow_cloud

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wP) ::  coll_eff      !< collision efficiency
       REAL(wp) ::  dc            !< diameter cloud
       REAL(wp) ::  ds            !< diameter snow
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  nc_c          !< dummy number concentration cloud dro
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_q !< riming rate mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vc            !< terminal fall velocitiy cloud
       REAL(wp) ::  vs            !< terminal fall velocitiy snow
       REAL(wp) ::  xc            !< mean mass cloud
       REAL(wp) ::  xs            !< mean mass snow
       REAL(wp) ::  conv_n        !< conversion rate snow number concentration
       REAL(wp) ::  conv_q        !< conversion rate snow mixing ratio
       REAL(wp) ::  const1        !< constant for riming

       REAL(wp), PARAMETER ::  const0 = 1.0_wp / ( dc_coll - dc_crit )  !< constant for riming calculation
       REAL(wp), PARAMETER ::  const3 = 1.0_wp / ( temp_mult_opt - temp_mult_min ) !< constant for snow multiplication
       REAL(wp), PARAMETER ::  const4 = 1.0_wp / ( temp_mult_opt - temp_mult_max ) !< constant for snow multiplication
       REAL(wp), PARAMETER ::  const5 = alpha_spacefilling * rho_l / rho_i

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                IF ( microphysics_morrison )  THEN
                   nc_c = nc(k,j,i)
                ELSE
                   nc_c = nc_const
                ENDIF
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Calculate mean mass, mean diameter for both graupel and rain
                xc = mean_mass( cloud_species%cloud, qc(k,j,i), nc_c     , hyrho(k) )
                xs = mean_mass( cloud_species%snow,  qs(k,j,i), ns(k,j,i), hyrho(k) )

                dc = mean_diameter( cloud_species%cloud,  xc )
                ds = mean_diameter( cloud_species%ice,    xs )
!
!--             Apply riming only above threshold
                IF ( qc(k,j,i) > q_thres  .AND.  qs(k,j,i) > eps_sb_coll  .AND.                    &
                     dc > dc_crit  .AND.  ds > df_crit )  THEN
!
!--                Calcualte temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Calculate terminal velocites
                   vc = terminal_fall_velocity(cloud_species%cloud, xc, hyrho(k) )
                   vs = terminal_fall_velocity(cloud_species%snow,  xs, hyrho(k) )
!
!--                Calculate collision efficiency
                   const1 = const0 * cloud_species%snow%coll_eff
                   coll_eff = MIN( cloud_species%snow%coll_eff,                                    &
                                   MAX( const1 * ( dc - dc_crit ), ecoll_min) )
!
!--                Calcualte riming rates
                   riming_rate_n = pi / 4.0_wp * coll_eff * ns(k,j,i) * nc_c *                     &
                                ( rime_snow_cloud_delta_n_aa * ds**2   +                           &
                                  rime_snow_cloud_delta_n_ab * ds * dc +                           &
                                  rime_snow_cloud_delta_n_bb * dc**2)  *                           &
                            SQRT( rime_snow_cloud_theta_n_aa * vs**2   -                           &
                                  rime_snow_cloud_theta_n_ab * vs * vc +                           &
                                  rime_snow_cloud_theta_n_bb * vc**2   +                           &
                                  cloud_species%snow%sigma_v**2 )

                   riming_rate_q = pi / 4.0_wp * coll_eff * ns(k,j,i) * qc(k,j,i) *                &
                                ( rime_snow_cloud_delta_q_aa * ds**2 +                             &
                                  rime_snow_cloud_delta_q_ab * ds * dc +                           &
                                  rime_snow_cloud_delta_q_bb * dc**2) *                            &
                            SQRT( rime_snow_cloud_theta_q_aa * vs**2 -                             &
                                  rime_snow_cloud_theta_q_ab * vs * vc +                           &
                                  rime_snow_cloud_theta_q_bb * vc**2 +                             &
                                  cloud_species%snow%sigma_v**2 )
!
!--                Limit riming rates
                   riming_rate_q = MIN( qc(k,j,i) / dt_micro, MAX( riming_rate_q, 0.0_wp ) )
                   riming_rate_n = MIN( nc_c      / dt_micro, MAX( riming_rate_n, 0.0_wp ) )
!
!--                Apply riming snow - cloud
                   qs(k,j,i) = qs(k,j,i) + riming_rate_q * dt_micro * flag
                   qc(k,j,i) = qc(k,j,i) - riming_rate_q * dt_micro * flag
                   IF ( microphysics_morrison )  THEN
                      nc(k,j,i) = nc(k,j,i) - riming_rate_n * dt_micro * flag
                   ENDIF
!
!--                Store riming rate for rain riming
                   rime_snow_cloud(k,j,i) = riming_rate_q
!
!--                Ice multiplication during riming
                   mult_q = 0.0
                   IF ( temp <  t3  .AND.  ice_multiplication )  THEN
                      mult_1 = ( temp - temp_mult_min ) * const3
                      mult_2 = ( temp - temp_mult_max ) * const4
                      mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                      mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                      mult_n = c_mult * mult_1 * mult_2 * riming_rate_q
                      mult_q = mult_n * cloud_species%snow%x_min
                      mult_q = MIN( riming_rate_q, mult_q )

                      ni(k,j,i) = ni(k,j,i) + mult_n
                      qi(k,j,i) = qi(k,j,i) + mult_q
                      qs(k,j,i) = qs(k,j,i) - mult_q
                    ENDIF
!
!--                If deposition is negative or deposition is smaller than riming, conversion from snow
!--                to graupel takes place. Note, that here the riming rate of snow-rain is from the
!--                last time-step. However, the error is assumed to be small.
                   IF ( dep_rate_snow(k,j,i) < 0.0_wp  .OR.                                        &
                       dep_rate_snow(k,j,i) < rime_snow_cloud(k,j,i) + rime_snow_rain(k,j,i) )  THEN
!
!--                    Conversion from ice to graupel if diameter small enough and riming rates are larger
!--                    than deposition. If deposition is larger ice stay
!--                    factor alpha_spacefilling after SB 2006 and icon model
                       IF ( ds > dsg_conv ) THEN
                         conv_q = ( riming_rate_q - mult_q ) /                                     &
                                  ( const5 * ( pi / 6.0_wp * rho_i * ds**3 / xs - 1.0_wp ) )
                         conv_q = MIN( qs(k,j,i) / dt_micro , conv_q )
                         conv_n = conv_q / MAX( xs, x_conv )
                         conv_n = MIN( ns(k,j,i) / dt_micro, conv_n )
                      ELSE
                         conv_q = 0.0_wp
                         conv_n = 0.0_wp
                      ENDIF
!
!--                   Apply conversion from ice to grauppel
                      qs(k,j,i) = qs(k,j,i) - conv_q * dt_micro * flag
                      qg(k,j,i) = qg(k,j,i) + conv_q * dt_micro * flag
                      ns(k,j,i) = ns(k,j,i) - conv_n * dt_micro * flag
                      ng(k,j,i) = ng(k,j,i) + conv_n * dt_micro * flag
                   ENDIF

                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE riming_snow_cloud


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of cloud and snow (Seifert and Beheng, 2006).
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_snow_cloud_ij( i, j )

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wP) ::  coll_eff      !< collision efficiency
       REAL(wp) ::  dc            !< diameter cloud
       REAL(wp) ::  ds            !< diameter snow
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  mult_n        !< ice multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  nc_c          !< dummy number concentration cloud dro
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_q !< riming rate mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vc            !< terminal fall velocitiy cloud
       REAL(wp) ::  vs            !< terminal fall velocitiy snow
       REAL(wp) ::  xc            !< mean mass cloud
       REAL(wp) ::  xs            !< mean mass snow
       REAL(wp) ::  conv_n        !< conversion rate snow number concentration
       REAL(wp) ::  conv_q        !< conversion rate snow mixing ratio
       REAL(wp) ::  const1        !< constant for riming

       REAL(wp), PARAMETER ::  const0 = 1.0_wp / ( dc_coll - dc_crit )  !< constant for riming calculation
       REAL(wp), PARAMETER ::  const3 = 1.0_wp / ( temp_mult_opt - temp_mult_min ) !< constant for snow multiplication
       REAL(wp), PARAMETER ::  const4 = 1.0_wp / ( temp_mult_opt - temp_mult_max ) !< constant for snow multiplication
       REAL(wp), PARAMETER ::  const5 = alpha_spacefilling * rho_l / rho_i

       DO  k = nzb+1, nzt
          IF ( microphysics_morrison )  THEN
             nc_c = nc(k,j,i)
          ELSE
             nc_c = nc_const
          ENDIF
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Calculate mean mass, mean diameter for both graupel and rain
          xc = mean_mass( cloud_species%cloud, qc(k,j,i), nc_c     , hyrho(k) )
          xs = mean_mass( cloud_species%snow,  qs(k,j,i), ns(k,j,i), hyrho(k) )

          dc = mean_diameter( cloud_species%cloud,  xc )
          ds = mean_diameter( cloud_species%ice,    xs )
!
!--       Apply riming only above threshold
          IF ( qc(k,j,i) > q_thres  .AND.  qs(k,j,i) > eps_sb_coll  .AND.                          &
               dc > dc_crit  .AND.  ds > df_crit )  THEN
!
!--          Calcualte temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Calculate terminal velocites
             vc = terminal_fall_velocity(cloud_species%cloud, xc, hyrho(k) )
             vs = terminal_fall_velocity(cloud_species%snow,  xs, hyrho(k) )
!
!--          Calculate collision efficiency
             const1 = const0 * cloud_species%snow%coll_eff
             coll_eff = MIN( cloud_species%snow%coll_eff,                                          &
                             MAX( const1 * ( dc - dc_crit ), ecoll_min) )
!
!--          Calcualte riming rates
             riming_rate_n = pi / 4.0_wp * coll_eff * ns(k,j,i) * nc_c *                           &
                          ( rime_snow_cloud_delta_n_aa * ds**2   +                                 &
                            rime_snow_cloud_delta_n_ab * ds * dc +                                 &
                            rime_snow_cloud_delta_n_bb * dc**2)  *                                 &
                      SQRT( rime_snow_cloud_theta_n_aa * vs**2   -                                 &
                            rime_snow_cloud_theta_n_ab * vs * vc +                                 &
                            rime_snow_cloud_theta_n_bb * vc**2   +                                 &
                            cloud_species%snow%sigma_v**2 )

             riming_rate_q = pi / 4.0_wp * coll_eff * ns(k,j,i) * qc(k,j,i) *                      &
                          ( rime_snow_cloud_delta_q_aa * ds**2 +                                   &
                            rime_snow_cloud_delta_q_ab * ds * dc +                                 &
                            rime_snow_cloud_delta_q_bb * dc**2) *                                  &
                      SQRT( rime_snow_cloud_theta_q_aa * vs**2 -                                   &
                            rime_snow_cloud_theta_q_ab * vs * vc +                                 &
                            rime_snow_cloud_theta_q_bb * vc**2 +                                   &
                            cloud_species%snow%sigma_v**2 )
!
!--          Limit riming rates
             riming_rate_q = MIN( qc(k,j,i) / dt_micro, MAX( riming_rate_q, 0.0_wp ) )
             riming_rate_n = MIN( nc_c      / dt_micro, MAX( riming_rate_n, 0.0_wp ) )
!
!--          Apply riming snow - cloud
             qs(k,j,i) = qs(k,j,i) + riming_rate_q * dt_micro * flag
             qc(k,j,i) = qc(k,j,i) - riming_rate_q * dt_micro * flag
             IF ( microphysics_morrison )  THEN
                nc(k,j,i) = nc(k,j,i) - riming_rate_n * dt_micro * flag
             ENDIF
!
!--          Store riming rate for rain riming
             rime_snow_cloud(k,j,i) = riming_rate_q
!
!--          Ice multiplication during riming
             mult_q = 0.0
             IF ( temp <  t3  .AND.  ice_multiplication )  THEN
                mult_1 = ( temp - temp_mult_min ) * const3
                mult_2 = ( temp - temp_mult_max ) * const4
                mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp ) )
                mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp ) )
                mult_n = c_mult * mult_1 * mult_2 * riming_rate_q
                mult_q = mult_n * cloud_species%snow%x_min
                mult_q = MIN( riming_rate_q, mult_q )

                ni(k,j,i) = ni(k,j,i) + mult_n
                qi(k,j,i) = qi(k,j,i) + mult_q
                qs(k,j,i) = qs(k,j,i) - mult_q
              ENDIF
!
!--          If deposition is negative or deposition is smaller than riming, conversion from snow
!--          to graupel takes place. Note, that here the riming rate of snow-rain is from the
!--          last time-step. However, the error is assumed to be small.
             IF ( dep_rate_snow(k,j,i) < 0.0_wp  .OR.                                              &
                  dep_rate_snow(k,j,i) < rime_snow_cloud(k,j,i) + rime_snow_rain(k,j,i) )  THEN
!
!--              Conversion from ice to graupel if diameter small enough and riming rates are larger
!--              than deposition. If deposition is larger ice stay
!--              factor alpha_spacefilling after SB 2006 and icon model
                 IF ( ds > dsg_conv ) THEN
                   conv_q = ( riming_rate_q - mult_q ) /                                           &
                            ( const5 * ( pi / 6.0_wp * rho_i * ds**3 / xs - 1.0_wp ) )
                   conv_q = MIN( qs(k,j,i) / dt_micro , conv_q )
                   conv_n = conv_q / MAX( xs, x_conv )
                   conv_n = MIN( ns(k,j,i) / dt_micro, conv_n )
                ELSE
                   conv_q = 0.0_wp
                   conv_n = 0.0_wp
                ENDIF
!
!--             Apply conversion from ice to grauppel
                qs(k,j,i) = qs(k,j,i) - conv_q * dt_micro * flag
                qg(k,j,i) = qg(k,j,i) + conv_q * dt_micro * flag
                ns(k,j,i) = ns(k,j,i) - conv_n * dt_micro * flag
                ng(k,j,i) = ng(k,j,i) + conv_n * dt_micro * flag
             ENDIF

          ENDIF

       ENDDO

    END SUBROUTINE riming_snow_cloud_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of rain and snow (Seifert and Beheng, 2006).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_snow_rain

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  dr            !< diameter rain
       REAL(wp) ::  ds            !< diameter snow
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  mult_n        !< snow multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_qs!< riming rate snow mixing ratio
       REAL(wp) ::  riming_rate_qr!< riming rate rain mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vr            !< terminal fall velocitiy rain
       REAL(wp) ::  vs            !< terminal fall velocitiy snow
       REAL(wp) ::  xr            !< mean mass rain
       REAL(wp) ::  xs            !< mean mass snow

       REAL(wp), PARAMETER ::  const3 = 1.0_wp / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = 1.0_wp / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Calculate mean mass, mean diameter for both graupel and rain
                xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
                xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )

                dr = mean_diameter( cloud_species%rain, xr )
                ds = mean_diameter( cloud_species%snow, xs )
!
!--             Apply ice-rain riming only above thresholds
                IF ( qr(k,j,i) > eps_sb_coll  .AND.  qs(k,j,i) > qr_crit  .AND.  ds > df_crit      &
                     .AND.  nr(k,j,i) > 0.0_wp  .AND.  ns(k,j,i) > 0.0_wp  )  THEN
!
!--                Calcualte temperature
                   temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--                Calculate terminal velocites
                   vr = terminal_fall_velocity(cloud_species%rain, xr, hyrho(k) )
                   vs = terminal_fall_velocity(cloud_species%snow, xs, hyrho(k) )
!
!--                Calcualte riming rates
                   riming_rate_n = pi / 4.0_wp * ns(k,j,i) * nr(k,j,i) *                           &
                                ( rime_snow_rain_delta_n_aa * ds**2   +                            &
                                  rime_snow_rain_delta_n_ab * ds * dr +                            &
                                  rime_snow_rain_delta_n_bb * dr**2 ) *                            &
                            SQRT( rime_snow_rain_theta_n_aa * vs**2   -                            &
                                  rime_snow_rain_theta_n_ab * vs * vr +                            &
                                  rime_snow_rain_theta_n_bb * vr**2   +                            &
                                  cloud_species%snow%sigma_v**2 )

                   riming_rate_qr = pi / 4.0_wp * ns(k,j,i) * qr(k,j,i) *                          &
                                ( rime_snow_rain_delta_n_aa * ds**2 +                              &
                                  rime_snow_rain_delta_q_ab * ds * dr +                            &
                                  rime_snow_rain_delta_q_bb * dr**2 ) *                            &
                            SQRT( rime_snow_rain_theta_n_aa * vs**2 -                              &
                                  rime_snow_rain_theta_q_ab * vs * vr +                            &
                                  rime_snow_rain_theta_q_bb * vr**2 +                              &
                                  cloud_species%snow%sigma_v**2 )

                   riming_rate_qs = pi / 4.0_wp * nr(k,j,i) * qs(k,j,i) *                          &
                                ( rime_snow_rain_delta_q_aa * ds**2 +                              &
                                  rime_snow_rain_delta_q_ab * ds * dr +                            &
                                  rime_snow_rain_delta_q_bb * dr**2 ) *                            &
                            SQRT( rime_snow_rain_theta_n_aa * vs**2 -                              &
                                  rime_snow_rain_theta_q_ab * vs * vr +                            &
                                  rime_snow_rain_theta_q_bb * vr**2 +                              &
                                  cloud_species%snow%sigma_v**2 )

!                    riming_rate_qs = pi / 4.0_wp * nr(k,j,i) * qs(k,j,i) *                          &
!                                 ( rime_snow_rain_delta_q_aa * ds**2 +                              &
!                                   rime_snow_rain_delta_q_ab * ds * dr +                            &
!                                   rime_snow_rain_delta_n_bb * dr**2 ) *                            &
!                             SQRT( rime_snow_rain_theta_q_aa * vs**2 -                              &
!                                   rime_snow_rain_theta_q_ab * vs * vr +                            &
!                                   rime_snow_rain_theta_n_bb * vr**2 +                              &
!                                   cloud_species%ice%sigma_v**2 )
!
!--                Limit riming rates
                   riming_rate_qr = MIN( qr(k,j,i) / dt_micro, MAX( riming_rate_qr, 0.0_wp ) )
                   riming_rate_qs = MIN( qs(k,j,i) / dt_micro, MAX( riming_rate_qs, 0.0_wp ) )
                   riming_rate_n  = MIN( nr(k,j,i) / dt_micro, MAX( riming_rate_n,  0.0_wp ) )
!
!--                store snow-rain riming rate for ice-cloud riming
                   rime_snow_rain(k,j,i) = riming_rate_qr
!
!--                If deposition is larger than riming snow stays snow
                   IF ( dep_rate_ice(k,j,i) > 0.0_wp  .AND.                                        &
                        dep_rate_ice(k,j,i) < rime_snow_rain(k,j,i) + rime_snow_cloud(k,j,i) ) THEN
!
!--                   Apply riming snow - rain
                      qs(k,j,i) = qs(k,j,i) + riming_rate_qr * dt_micro * flag
                      qr(k,j,i) = qr(k,j,i) - riming_rate_qr * dt_micro * flag
                      nr(k,j,i) = nr(k,j,i) - riming_rate_n  * dt_micro * flag
!
!--                   Ice multiplication during riming snow - rain
                      mult_q = 0.0
                      IF ( temp < t3  .AND.  ice_multiplication )  THEN
                         mult_1 = ( temp - temp_mult_min ) * const3
                         mult_2 = ( temp - temp_mult_max ) * const4
                         mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp) )
                         mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp) )
                         mult_n = c_mult * mult_1 * mult_2 * riming_rate_qr
                         mult_q = mult_n * cloud_species%ice%x_min
                         mult_q = MIN( riming_rate_qr, mult_q )

                         ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                         qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
                         qs(k,j,i) = qs(k,j,i) - mult_q * dt_micro * flag
                      ENDIF
!
!--                If deposition is negative (i.e. sublimation or riming is stronger, a different
!--                treatment is conducted: ice might convert to graupel
                   ELSE
!
!--                   Ensure that riming rate is also not larger than ni
                      riming_rate_n  = MIN( MIN( nr(k,j,i) / dt_micro, ns(k,j,i) / dt_micro ),     &
                                             MAX( riming_rate_n,  0.0_wp ) )
!
!--                   Substract all riming rates, those riming rates are added to species (either
!--                   snow again or graupel furhter down)
                      qs(k,j,i) = qs(k,j,i) - riming_rate_qs * dt_micro * flag
                      qr(k,j,i) = qr(k,j,i) - riming_rate_qr * dt_micro * flag
                      ns(k,j,i) = ns(k,j,i) - riming_rate_n  * dt_micro * flag
                      nr(k,j,i) = nr(k,j,i) - riming_rate_n  * dt_micro * flag
!
!--                   Ice multiplication
                      mult_q = 0.0_wp
                      mult_n = 0.0_wp
                      IF ( temp < t3  .AND.  ice_multiplication )  THEN
                         mult_1 = ( temp - temp_mult_min ) * const3
                         mult_2 = ( temp - temp_mult_max ) * const4
                         mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp) )
                         mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp) )
                         mult_n = c_mult * mult_1 * mult_2 * riming_rate_qr
                         mult_q = mult_n * cloud_species%ice%x_min
                         mult_q = MIN( riming_rate_qr, mult_q )
                      ENDIF
!
!--                   shedding of rain at warm temperatures
                      IF ( temp >= t3 )  THEN
!
!--                      but with modified nr
                         xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
                         nr(k,j,i) = nr(k,j,i) + riming_rate_qr / xr * dt_micro * flag
                         ns(k,j,i) = ns(k,j,i) + riming_rate_n  * dt_micro * flag
                         qs(k,j,i) = qs(k,j,i) + riming_rate_qs * dt_micro * flag
                         qr(k,j,i) = qr(k,j,i) + riming_rate_qr * dt_micro * flag
!
!--                   at tempatures below freezing point graupel is generated
                      ELSE
!
!--                      new ice particles from multiplication
                         ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                         qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
!
!--                      and riming to graupel
                         ng(k,j,i) = ng(k,j,i) + riming_rate_n * dt_micro * flag
                         qg(k,j,i) = qg(k,j,i) + ( riming_rate_qs + riming_rate_qr - mult_q ) *     &
                                                                                      dt_micro * flag
                      END IF

                   ENDIF ! loop if deposition or riming is stronger
                ENDIF ! loop if thresholds for riming are exceeded, otherwise nothings happen
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE riming_snow_rain


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Riming of rain and snow (Seifert and Beheng, 2006).
!> Call for grid point i,j.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE riming_snow_rain_ij( i, j )

       INTEGER(iwp) ::  i         !< loop index
       INTEGER(iwp) ::  j         !< loop index
       INTEGER(iwp) ::  k         !< loop index

       REAL(wp) ::  dr            !< diameter rain
       REAL(wp) ::  ds            !< diameter snow
       REAL(wp) ::  flag          !< flag to mask topography grid points
       REAL(wp) ::  mult_n        !< snow multiplaction for number concentration
       REAL(wp) ::  mult_q        !< ice multiplaction
       REAL(wp) ::  mult_1        !< factor for ice multiplaction
       REAL(wp) ::  mult_2        !< factor for ice multiplaction
       REAL(wp) ::  riming_rate_n !< riming rate number concentration
       REAL(wp) ::  riming_rate_qs!< riming rate snow mixing ratio
       REAL(wp) ::  riming_rate_qr!< riming rate rain mixing ratio
       REAL(wp) ::  temp          !< air temperature
       REAL(wp) ::  vr            !< terminal fall velocitiy rain
       REAL(wp) ::  vs            !< terminal fall velocitiy snow
       REAL(wp) ::  xr            !< mean mass rain
       REAL(wp) ::  xs            !< mean mass snow

       REAL(wp), PARAMETER ::  const3 = 1.0_wp / ( temp_mult_opt - temp_mult_min ) !< constant for ice multiplication
       REAL(wp), PARAMETER ::  const4 = 1.0_wp / ( temp_mult_opt - temp_mult_max ) !< constant for ice multiplication

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Calculate mean mass, mean diameter for both graupel and rain
          xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
          xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )

          dr = mean_diameter( cloud_species%rain, xr )
          ds = mean_diameter( cloud_species%snow, xs )
!
!--       Apply ice-rain riming only above thresholds
          IF ( qr(k,j,i) > eps_sb_coll  .AND.  qs(k,j,i) > qr_crit  .AND.  ds > df_crit            &
               .AND.  nr(k,j,i) > 0.0_wp  .AND.  ns(k,j,i) > 0.0_wp  )  THEN
!
!--          Calcualte temperature
             temp = exner(k) * pt(k,j,i) + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--          Calculate terminal velocites
             vr = terminal_fall_velocity(cloud_species%rain, xr, hyrho(k) )
             vs = terminal_fall_velocity(cloud_species%snow, xs, hyrho(k) )
!
!--          Calcualte riming rates
             riming_rate_n = pi / 4.0_wp * ns(k,j,i) * nr(k,j,i) *                                 &
                          ( rime_snow_rain_delta_n_aa * ds**2   +                                  &
                            rime_snow_rain_delta_n_ab * ds * dr +                                  &
                            rime_snow_rain_delta_n_bb * dr**2 ) *                                  &
                      SQRT( rime_snow_rain_theta_n_aa * vs**2   -                                  &
                            rime_snow_rain_theta_n_ab * vs * vr +                                  &
                            rime_snow_rain_theta_n_bb * vr**2   +                                  &
                            cloud_species%snow%sigma_v**2 )

             riming_rate_qr = pi / 4.0_wp * ns(k,j,i) * qr(k,j,i) *                                &
                          ( rime_snow_rain_delta_n_aa * ds**2 +                                    &
                            rime_snow_rain_delta_q_ab * ds * dr +                                  &
                            rime_snow_rain_delta_q_bb * dr**2 ) *                                  &
                      SQRT( rime_snow_rain_theta_n_aa * vs**2 -                                    &
                            rime_snow_rain_theta_q_ab * vs * vr +                                  &
                            rime_snow_rain_theta_q_bb * vr**2 +                                    &
                            cloud_species%snow%sigma_v**2 )

             riming_rate_qs = pi / 4.0_wp * nr(k,j,i) * qs(k,j,i) *                                &
                          ( rime_snow_rain_delta_q_aa * ds**2 +                                    &
                            rime_snow_rain_delta_q_ab * ds * dr +                                  &
                            rime_snow_rain_delta_q_bb * dr**2 ) *                                  &
                      SQRT( rime_snow_rain_theta_n_aa * vs**2 -                                    &
                            rime_snow_rain_theta_q_ab * vs * vr +                                  &
                            rime_snow_rain_theta_q_bb * vr**2 +                                    &
                            cloud_species%snow%sigma_v**2 )

!              riming_rate_qs = pi / 4.0_wp * nr(k,j,i) * qs(k,j,i) *                                &
!                           ( rime_snow_rain_delta_q_aa * ds**2 +                                    &
!                             rime_snow_rain_delta_q_ab * ds * dr +                                  &
!                             rime_snow_rain_delta_n_bb * dr**2 ) *                                  &
!                       SQRT( rime_snow_rain_theta_q_aa * vs**2 -                                    &
!                             rime_snow_rain_theta_q_ab * vs * vr +                                  &
!                             rime_snow_rain_theta_n_bb * vr**2 +                                    &
!                             cloud_species%ice%sigma_v**2 )
!
!--          Limit riming rates
             riming_rate_qr = MIN( qr(k,j,i) / dt_micro, MAX( riming_rate_qr, 0.0_wp ) )
             riming_rate_qs = MIN( qs(k,j,i) / dt_micro, MAX( riming_rate_qs, 0.0_wp ) )
             riming_rate_n  = MIN( nr(k,j,i) / dt_micro, MAX( riming_rate_n,  0.0_wp ) )
!
!--          store snow-rain riming rate for ice-cloud riming
             rime_snow_rain(k,j,i) = riming_rate_qr
!
!--          If deposition is larger than riming snow stays snow
             IF ( dep_rate_ice(k,j,i) > 0.0_wp  .AND.                                             &
                  dep_rate_ice(k,j,i) < rime_snow_rain(k,j,i) + rime_snow_cloud(k,j,i) ) THEN
!
!--             Apply riming snow - rain
                qs(k,j,i) = qs(k,j,i) + riming_rate_qr * dt_micro * flag
                qr(k,j,i) = qr(k,j,i) - riming_rate_qr * dt_micro * flag
                nr(k,j,i) = nr(k,j,i) - riming_rate_n  * dt_micro * flag
!
!--             Ice multiplication during riming snow - rain
                mult_q = 0.0
                IF ( temp < t3  .AND.  ice_multiplication )  THEN
                   mult_1 = ( temp - temp_mult_min ) * const3
                   mult_2 = ( temp - temp_mult_max ) * const4
                   mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp) )
                   mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp) )
                   mult_n = c_mult * mult_1 * mult_2 * riming_rate_qr
                   mult_q = mult_n * cloud_species%ice%x_min
                   mult_q = MIN( riming_rate_qr, mult_q )

                   ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                   qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
                   qs(k,j,i) = qs(k,j,i) - mult_q * dt_micro * flag
                ENDIF
!
!--          If deposition is negative (i.e. sublimation or riming is stronger, a different
!--          treatment is conducted: ice might convert to graupel
             ELSE
!
!--             Ensure that riming rate is also not larger than ni
                riming_rate_n  = MIN( MIN( nr(k,j,i) / dt_micro, ns(k,j,i) / dt_micro ),          &
                                       MAX( riming_rate_n,  0.0_wp ) )
!
!--             Substract all riming rates, those riming rates are added to species (either snow
!--             again or graupel furhter down)
                qs(k,j,i) = qs(k,j,i) - riming_rate_qs * dt_micro * flag
                qr(k,j,i) = qr(k,j,i) - riming_rate_qr * dt_micro * flag
                ns(k,j,i) = ns(k,j,i) - riming_rate_n  * dt_micro * flag
                nr(k,j,i) = nr(k,j,i) - riming_rate_n  * dt_micro * flag
!
!--             Ice multiplication
                mult_q = 0.0_wp
                mult_n = 0.0_wp
                IF ( temp < t3  .AND.  ice_multiplication )  THEN
                   mult_1 = ( temp - temp_mult_min ) * const3
                   mult_2 = ( temp - temp_mult_max ) * const4
                   mult_1 = MAX( 0.0_wp, MIN( mult_1, 1.0_wp) )
                   mult_2 = MAX( 0.0_wp, MIN( mult_2, 1.0_wp) )
                   mult_n = c_mult * mult_1 * mult_2 * riming_rate_qr
                   mult_q = mult_n * cloud_species%ice%x_min
                   mult_q = MIN( riming_rate_qr, mult_q )
                ENDIF
!
!--             shedding of rain at warm temperatures
                IF ( temp >= t3 )  THEN
!
!--                but with modified nr
                   xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
                   nr(k,j,i) = nr(k,j,i) + riming_rate_qr / xr * dt_micro * flag
                   ns(k,j,i) = ns(k,j,i) + riming_rate_n  * dt_micro * flag
                   qs(k,j,i) = qs(k,j,i) + riming_rate_qs * dt_micro * flag
                   qr(k,j,i) = qr(k,j,i) + riming_rate_qr * dt_micro * flag
!
!--             at tempatures below freezing point graupel is generated
                ELSE
!
!--                new ice particles from multiplication
                   ni(k,j,i) = ni(k,j,i) + mult_n * dt_micro * flag
                   qi(k,j,i) = qi(k,j,i) + mult_q * dt_micro * flag
!
!--                and riming to graupel
                   ng(k,j,i) = ng(k,j,i) + riming_rate_n * dt_micro * flag
                   qg(k,j,i) = qg(k,j,i) + ( riming_rate_qs + riming_rate_qr - mult_q ) *          &
                                                                                     dt_micro * flag
                END IF

             ENDIF ! loop if deposition or riming is stronger
          ENDIF ! loop if thresholds for riming are exceeded, otherwise nothings happen

       ENDDO

    END SUBROUTINE riming_snow_rain_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Evaporation of precipitable water. Condensation is neglected for precipitable water.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE evaporation_rain

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  evap_nr           !<
       REAL(wp)     ::  f_vent            !<
       REAL(wp)     ::  flag              !< flag to mask topography grid points
       REAL(wp)     ::  g_evap            !<
       REAL(wp)     ::  lambda_r          !<
       REAL(wp)     ::  mu_r              !<
       REAL(wp)     ::  mu_r_2            !<
       REAL(wp)     ::  mu_r_5d2          !<
       REAL(wp)     ::  nr_0              !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xr                !<

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( qr(k,j,i) > eps_sb  .AND.  nr(k,j,i) > 0.0_wp )  THEN
!
!--                Call calculation of supersaturation
                   CALL supersaturation ( i, j, k )
!
!--                Evaporation needs only to be calculated in subsaturated regions
                   IF ( sat < 0.0_wp )  THEN
!
!--                   Actual temperature:
                      temp = t_l + lv_d_cp * ( qc(k,j,i) + qr(k,j,i) )

                      g_evap = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *                      &
                                            l_v / ( thermal_conductivity_l  * temp )               &
                                          + r_v * temp / ( diff_coeff_l * e_s )                    &
                                        )
!
!--                   Mean weight of rain drops
                      !xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
                      xr = hyrho(k) * qr(k,j,i) / nr(k,j,i)
!
!--                   Weight averaged diameter of rain drops:
                      dr = ( xr * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--                   Compute ventilation factor and intercept parameter
!--                   (Seifert and Beheng, 2006; Seifert, 2008):
                      IF ( ventilation_effect )  THEN
!
!--                      Shape parameter of gamma distribution (Milbrandt and Yau,
!--                      2005; Stevens and Seifert, 2008):
                         mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp * ( dr - 1.4E-3_wp ) ) )
!
!--                      Slope parameter of gamma distribution (Seifert, 2008):
                         lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) * ( mu_r + 1.0_wp )    &
                                    )**( 1.0_wp / 3.0_wp ) / dr

                         mu_r_2   = mu_r + 2.0_wp
                         mu_r_5d2 = mu_r + 2.5_wp

                         f_vent = a_vent * gamm( mu_r_2 ) *                                        &
                                  lambda_r**( -mu_r_2 ) + b_vent *                                 &
                                  schmidt_p_1d3 * SQRT( a_term / kin_vis_air ) *                   &
                                  gamm( mu_r_5d2 ) * lambda_r**( -mu_r_5d2 ) *                     &
                                  ( 1.0_wp -                                                       &
                                    0.5_wp       * ( b_term / a_term ) *                           &
                                    ( lambda_r / ( c_term + lambda_r ) )**mu_r_5d2          -      &
                                    0.125_wp     * ( b_term / a_term )**2 *                        &
                                    ( lambda_r / ( 2.0_wp * c_term + lambda_r ) )**mu_r_5d2 -      &
                                    0.0625_wp    * ( b_term / a_term )**3 *                        &
                                    ( lambda_r / ( 3.0_wp * c_term + lambda_r ) )**mu_r_5d2 -      &
                                    0.0390625_wp * ( b_term / a_term )**4 *                        &
                                    ( lambda_r / ( 4.0_wp * c_term + lambda_r ) )**mu_r_5d2        &
                                  )

                         nr_0   = nr(k,j,i) * lambda_r**( mu_r + 1.0_wp ) /    &
                                  gamm( mu_r + 1.0_wp )
                      ELSE
                         f_vent = 1.0_wp
                         nr_0   = nr(k,j,i) * dr
                      ENDIF
!
!--                   Evaporation rate of rain water content (Seifert and Beheng, 2006):
                      evap    = 2.0_wp * pi * nr_0 * g_evap * f_vent * sat / hyrho(k)
                      evap    = MAX( evap, -qr(k,j,i) / dt_micro )
                      evap_nr = MAX( c_evap * evap / xr * hyrho(k), - nr(k,j,i) / dt_micro )
                      qr(k,j,i) = qr(k,j,i) + evap    * dt_micro * flag
                      nr(k,j,i) = nr(k,j,i) + evap_nr * dt_micro * flag

                   ENDIF
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE evaporation_rain


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Evaporation of precipitable water. Condensation is neglected for precipitable water. Call for
!> grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE evaporation_rain_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i                 !< loop index
       INTEGER(iwp) ::  j                 !< loop index
       INTEGER(iwp) ::  k                 !< loop index

       REAL(wp)     ::  dr                !<
       REAL(wp)     ::  evap              !<
       REAL(wp)     ::  evap_nr           !<
       REAL(wp)     ::  f_vent            !<
       REAL(wp)     ::  flag              !< flag to indicate first grid level above surface
       REAL(wp)     ::  g_evap            !<
       REAL(wp)     ::  lambda_r          !<
       REAL(wp)     ::  mu_r              !<
       REAL(wp)     ::  mu_r_2            !<
       REAL(wp)     ::  mu_r_5d2          !<
       REAL(wp)     ::  nr_0              !<
       REAL(wp)     ::  temp              !<
       REAL(wp)     ::  xr                !<

       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          IF ( qr(k,j,i) > eps_sb  .AND.  nr(k,j,i) > 0.0_wp )  THEN
!
!--          Call calculation of supersaturation
             CALL supersaturation ( i, j, k )
!
!--          Evaporation needs only to be calculated in subsaturated regions
             IF ( sat < 0.0_wp )  THEN
!
!--             Actual temperature:
                temp = t_l + lv_d_cp * ( qc(k,j,i) + qr(k,j,i) )

                g_evap = 1.0_wp / ( ( l_v / ( r_v * temp ) - 1.0_wp ) *                            &
                                      l_v / ( thermal_conductivity_l  * temp )                     &
                                    + r_v * temp / ( diff_coeff_l * e_s )                          &
                                  )
!
!--             Mean weight of rain drops
                !xr = mean_mass( cloud_species%rain, qr(k,j,i), nr(k,j,i), hyrho(k) )
                xr = hyrho(k) * qr(k,j,i) / nr(k,j,i)
!
!--             Weight averaged diameter of rain drops:
                dr = ( xr * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--             Compute ventilation factor and intercept parameter
!--             (Seifert and Beheng, 2006; Seifert, 2008):
                IF ( ventilation_effect )  THEN
!
!--                Shape parameter of gamma distribution (Milbrandt and Yau, 2005;
!--                Stevens and Seifert, 2008):
                   mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp * ( dr - 1.4E-3_wp ) ) )
!
!--                Slope parameter of gamma distribution (Seifert, 2008):
                   lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) * ( mu_r + 1.0_wp )          &
                              )**( 1.0_wp / 3.0_wp ) / dr

                   mu_r_2   = mu_r + 2.0_wp
                   mu_r_5d2 = mu_r + 2.5_wp

                   f_vent = a_vent * gamm( mu_r_2 ) * lambda_r**( -mu_r_2 ) +                      &
                            b_vent * schmidt_p_1d3 *                                               &
                            SQRT( a_term / kin_vis_air ) * gamm( mu_r_5d2 ) *                      &
                            lambda_r**( -mu_r_5d2 ) *                                              &
                            ( 1.0_wp -                                                             &
                              0.5_wp       * ( b_term / a_term ) *                                 &
                              ( lambda_r / ( c_term + lambda_r ) )**mu_r_5d2          -            &
                              0.125_wp     * ( b_term / a_term )**2 *                              &
                              ( lambda_r / ( 2.0_wp * c_term + lambda_r ) )**mu_r_5d2 -            &
                              0.0625_wp    * ( b_term / a_term )**3 *                              &
                              ( lambda_r / ( 3.0_wp * c_term + lambda_r ) )**mu_r_5d2 -            &
                              0.0390625_wp * ( b_term / a_term )**4 *                              &
                              ( lambda_r / ( 4.0_wp * c_term + lambda_r ) )**mu_r_5d2              &
                            )

                   nr_0   = nr(k,j,i) * lambda_r**( mu_r + 1.0_wp ) / gamm( mu_r + 1.0_wp )
                ELSE
                   f_vent = 1.0_wp
                   nr_0   = nr(k,j,i) * dr
                ENDIF
!
!--             Evaporation rate of rain water content (Seifert and Beheng, 2006):
                evap    = 2.0_wp * pi * nr_0 * g_evap * f_vent * sat / hyrho(k)
                evap    = MAX( evap, -qr(k,j,i) / dt_micro )
                evap_nr = MAX( c_evap * evap / xr * hyrho(k), - nr(k,j,i) / dt_micro )
                qr(k,j,i) = qr(k,j,i) + evap    * dt_micro * flag
                nr(k,j,i) = nr(k,j,i) + evap_nr * dt_micro * flag

             ENDIF
          ENDIF

       ENDDO

    END SUBROUTINE evaporation_rain_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of cloud droplets (Ackermann et al., 2009, MWR).
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_cloud


       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< loop index
       INTEGER(iwp) ::  j             !< loop index
       INTEGER(iwp) ::  k             !< loop index

       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  nc_c              !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nc !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qc !<

       sed_qc(nzt+1) = 0.0_wp
       sed_nc(nzt+1) = 0.0_wp

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzt, nzb+1, -1
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( microphysics_morrison )  THEN
                   nc_c = nc(k,j,i)
                ELSE
                   nc_c = nc_const
                ENDIF

!
!--             Sedimentation fluxes for number concentration are only calculated for cloud_scheme =
!--             'morrison'
                IF ( microphysics_morrison )  THEN
                   IF ( qc(k,j,i) > eps_sb  .AND.  nc(k,j,i) > eps_mr )  THEN
                      sed_nc(k) = sed_qc_const * ( qc(k,j,i) * hyrho(k) )**( 2.0_wp / 3.0_wp ) *   &
                                  ( nc(k,j,i) )**( 1.0_wp / 3.0_wp )
                   ELSE
                      sed_nc(k) = 0.0_wp
                   ENDIF

                   sed_nc(k) = MIN( sed_nc(k), hyrho(k) * dzu(k+1) *                               &
                                    nc(k,j,i) / dt_micro + sed_nc(k+1)                             &
                                  ) * flag

                   nc(k,j,i) = nc(k,j,i) + ( sed_nc(k+1) - sed_nc(k) ) * ddzu(k+1) /               &
                               hyrho(k) * dt_micro * flag
                ENDIF

                IF ( qc(k,j,i) > eps_sb .AND.  nc_c > eps_mr )  THEN
                   sed_qc(k) = sed_qc_const * nc_c**( -2.0_wp / 3.0_wp )  *                        &
                               ( qc(k,j,i) * hyrho(k) )**( 5.0_wp / 3.0_wp ) * flag
                ELSE
                   sed_qc(k) = 0.0_wp
                ENDIF

                sed_qc(k) = MIN( sed_qc(k), hyrho(k) * dzu(k+1) * q(k,j,i) / dt_micro + sed_qc(k+1)&
                               ) * flag

                q(k,j,i)  = q(k,j,i)  + ( sed_qc(k+1) - sed_qc(k) ) *                              &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                qc(k,j,i) = qc(k,j,i) + ( sed_qc(k+1) - sed_qc(k) ) *                              &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) - ( sed_qc(k+1) - sed_qc(k) ) *                              &
                                        ddzu(k+1) / hyrho(k) * lv_d_cp *                           &
                                        d_exner(k) * dt_micro           * flag

!
!--             Compute the precipitation rate due to cloud (fog) droplets.
                IF ( call_microphysics_at_all_substeps )  THEN
                   prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) *                                &
                                             weight_substep(intermediate_timestep_count) * flag
                   prr_cloud(k,j,i) = prr_cloud(k,j,i) + sed_qc(k) / hyrho(k) *                    &
                                             weight_substep(intermediate_timestep_count) * flag
                ELSE
                   prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) * flag
                   prr_cloud(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) * flag
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE sedimentation_cloud


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of cloud droplets (Ackermann et al., 2009, MWR).
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_cloud_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< loop index
       INTEGER(iwp) ::  j             !< loop index
       INTEGER(iwp) ::  k             !< loop index

       REAL(wp)     ::  flag    !< flag to indicate first grid level above surface
       REAL(wp)     ::  nc_c    !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nc  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qc  !<

       sed_qc(nzt+1) = 0.0_wp
       sed_nc(nzt+1) = 0.0_wp


       DO  k = nzt, nzb+1, -1
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
          IF ( microphysics_morrison )  THEN
             nc_c = nc(k,j,i)
          ELSE
             nc_c = nc_const
          ENDIF
!
!--       Sedimentation fluxes for number concentration are only calculated for cloud_scheme =
!--       'morrison'
          IF ( microphysics_morrison )  THEN
             IF ( qc(k,j,i) > eps_sb  .AND.  nc(k,j,i) > eps_mr )  THEN
                sed_nc(k) = sed_qc_const * ( qc(k,j,i) * hyrho(k) )**( 2.0_wp / 3.0_wp ) *         &
                            ( nc(k,j,i) )**( 1.0_wp / 3.0_wp )
             ELSE
                sed_nc(k) = 0.0_wp
             ENDIF

             sed_nc(k) = MIN( sed_nc(k), hyrho(k) * dzu(k+1) * nc(k,j,i) / dt_micro + sed_nc(k+1)  &
                            ) * flag

             nc(k,j,i) = nc(k,j,i) + ( sed_nc(k+1) - sed_nc(k) ) * ddzu(k+1) /                     &
                         hyrho(k) * dt_micro * flag
          ENDIF

          IF ( qc(k,j,i) > eps_sb  .AND.  nc_c > eps_mr )  THEN
             sed_qc(k) = sed_qc_const * nc_c**( -2.0_wp / 3.0_wp ) *                               &
                         ( qc(k,j,i) * hyrho(k) )**( 5.0_wp / 3.0_wp ) * flag
          ELSE
             sed_qc(k) = 0.0_wp
          ENDIF

          sed_qc(k) = MIN( sed_qc(k), hyrho(k) * dzu(k+1) * q(k,j,i) / dt_micro + sed_qc(k+1)      &
                         ) * flag

          q(k,j,i)  = q(k,j,i)  + ( sed_qc(k+1) - sed_qc(k) ) *                                    &
                                  ddzu(k+1) / hyrho(k) * dt_micro * flag
          qc(k,j,i) = qc(k,j,i) + ( sed_qc(k+1) - sed_qc(k) ) *                                    &
                                  ddzu(k+1) / hyrho(k) * dt_micro * flag
          pt(k,j,i) = pt(k,j,i) - ( sed_qc(k+1) - sed_qc(k) ) *                                    &
                                  ddzu(k+1) / hyrho(k) * lv_d_cp * d_exner(k) * dt_micro * flag

!
!--       Compute the precipitation rate of cloud (fog) droplets
          IF ( call_microphysics_at_all_substeps )  THEN
             prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) *                                      &
                                       weight_substep(intermediate_timestep_count) * flag
             prr_cloud(k,j,i) = prr_cloud(k,j,i) + sed_qc(k) / hyrho(k) *                          &
                                       weight_substep(intermediate_timestep_count) * flag
          ELSE
             prr(k,j,i) = prr(k,j,i) + sed_qc(k) / hyrho(k) * flag
             prr_cloud(k,j,i) = prr_cloud(k,j,i) + sed_qc(k) / hyrho(k) * flag
          ENDIF

       ENDDO

    END SUBROUTINE sedimentation_cloud_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of ice crystals
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_ice

       INTEGER(iwp) ::  i             !< loop index
       INTEGER(iwp) ::  j             !< loop index
       INTEGER(iwp) ::  k             !< loop index

       REAL(wp) ::  av = 6.39_wp                  !< parameter for calculating fall speed
       REAL(wp) ::  bv = 0.1666_wp                !< parameter (1/6)
       REAL(wp) ::  factor_sed_gamma_n = 0.76_wp !< factor for zeroth moment and µ =1/3 and nu=0
       REAL(wp) ::  factor_sed_gamma_q = 1.61_wp !< factor for first moment and µ =1/3 and nu=0
       REAL(wp) ::  flag                          !< flag to indicate first grid level
       REAL(wp) ::  xi = 0.0_wp                   !< mean mass of ice crystal
       REAL(wp) ::  vi = 0.0_wp                   !< mean fall speed of ice crystal

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_ni  !< sedimentation rate zeroth moment
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qi  !< sedimentation rate fist moment

       sed_qi(nzt+1) = 0.0_wp
       sed_ni(nzt+1) = 0.0_wp

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzt, nzb+1, -1

                IF ( ni(k,j,i) <= 0.0_wp )  THEN
                   xi = 0.0_wp
                ELSE
!--                Calculate mean mass of ice crystal
                   xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
                ENDIF
!
!--             Calculate fall speed of ice crystal
                vi = av * xi**bv
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Calculate sedimentation rate for each grid box
                IF ( qi(k,j,i) > eps_sb  .AND.  ni(k,j,i) >= 0.0_wp )  THEN
                   sed_qi(k) = qi(k,j,i) * vi * factor_sed_gamma_q * flag
                   sed_ni(k) = ni(k,j,i) * vi * factor_sed_gamma_n * flag
                ELSE
                   sed_qi(k) = 0.0_wp
                   sed_ni(k) = 0.0_wp
                ENDIF
!
!--             Calculate sedimentation: divergence of sedimentation flux
                sed_qi(k) = MIN( sed_qi(k), hyrho(k) * dzu(k+1) * q(k,j,i) / dt_micro + sed_qi(k+1)&
                               ) * flag

                q(k,j,i)  = q(k,j,i)  + ( sed_qi(k+1) - sed_qi(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * dt_micro * flag
                qi(k,j,i) = qi(k,j,i) + ( sed_qi(k+1) - sed_qi(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * dt_micro * flag
                ni(k,j,i) = ni(k,j,i) + ( sed_ni(k+1) - sed_ni(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * dt_micro * flag

                pt(k,j,i) = pt(k,j,i) - ( sed_qi(k+1) - sed_qi(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * l_s / c_p * d_exner(k) * dt_micro * flag
!
!--             Compute the precipitation rate of cloud (fog) droplets
                IF ( call_microphysics_at_all_substeps )  THEN
                   prr(k,j,i) = prr(k,j,i) + sed_qi(k) / hyrho(k) *                                &
                                             weight_substep(intermediate_timestep_count) * flag
                   prr_ice(k,j,i) = prr_ice(k,j,i) + sed_qi(k) / hyrho(k) *                        &
                                             weight_substep(intermediate_timestep_count) * flag
                ELSE
                   prr(k,j,i) = prr(k,j,i) + sed_qi(k) / hyrho(k) * flag
                   prr_ice(k,j,i) = prr_ice(k,j,i) + sed_qi(k) / hyrho(k) * flag
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE sedimentation_ice


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of ice crystals
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_ice_ij( i, j )

       INTEGER(iwp) ::  i             !<
       INTEGER(iwp) ::  j             !<
       INTEGER(iwp) ::  k             !<

       REAL(wp) ::  av = 6.39_wp   !< parameter for calculating fall speed
       REAL(wp) ::  bv = 0.1666_wp !<
       REAL(wp) ::  factor_sed_gamma_n = 0.76_wp
       REAL(wp) ::  factor_sed_gamma_q = 1.61_wp
       REAL(wp) ::  flag           !< flag to indicate first grid level
       REAL(wp) ::  xi = 0.0_wp    !< mean mass of ice crystal
       REAL(wp) ::  vi = 0.0_wp    !< mean fall speed of ice crystal

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_ni  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qi  !<

       sed_qi(nzt+1) = 0.0_wp
       sed_ni(nzt+1) = 0.0_wp

       DO  k = nzt, nzb+1, -1
          IF ( ni(k,j,i) <= 0.0_wp )  THEN
             xi = 0.0_wp
          ELSE
!--          Calculate mean mass of ice crystal limit to minimum and maximum mass
             xi = mean_mass( cloud_species%ice, qi(k,j,i), ni(k,j,i), hyrho(k) )
          ENDIF
!
!--       calculate fall speed of ice crystal
          vi = av * xi**bv
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Calculate sedimentation rate for each grid box
          IF ( qi(k,j,i) > eps_sb  .AND.  ni(k,j,i) >= 0.0_wp )  THEN
             sed_qi(k) = qi(k,j,i) * vi * factor_sed_gamma_q * flag
             sed_ni(k) = ni(k,j,i) * vi * factor_sed_gamma_n * flag
          ELSE
             sed_qi(k) = 0.0_wp
             sed_ni(k) = 0.0_wp
          ENDIF
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          sed_qi(k) = MIN( sed_qi(k), hyrho(k) * dzu(k+1) * qi(k,j,i) / dt_micro + sed_qi(k+1)     &
                         ) * flag
          sed_ni(k) = MIN( sed_ni(k), hyrho(k) * dzu(k+1) * ni(k,j,i) / dt_micro + sed_ni(k+1)     &
                         ) * flag

          q(k,j,i)  = q(k,j,i)  + ( sed_qi(k+1) - sed_qi(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          qi(k,j,i) = qi(k,j,i) + ( sed_qi(k+1) - sed_qi(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          ni(k,j,i) = ni(k,j,i) + ( sed_ni(k+1) - sed_ni(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          pt(k,j,i) = pt(k,j,i) - ( sed_qi(k+1) - sed_qi(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * l_s / c_p * d_exner(k) * dt_micro * flag
!
!--       Compute the precipitation rate of ice
          IF ( call_microphysics_at_all_substeps )  THEN
             prr(k,j,i) = prr(k,j,i) + sed_qi(k) / hyrho(k) *                                      &
                                       weight_substep(intermediate_timestep_count) * flag
             prr_ice(k,j,i) = prr_ice(k,j,i) + sed_qi(k) / hyrho(k) *                              &
                                       weight_substep(intermediate_timestep_count) * flag
          ELSE
             prr(k,j,i) = prr(k,j,i) + sed_qi(k) / hyrho(k) * flag
             prr_ice(k,j,i) = prr_ice(k,j,i) + sed_qi(k) / hyrho(k) * flag
          ENDIF

       ENDDO

    END SUBROUTINE sedimentation_ice_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of graupel assuming generalized gamma distibution (Seifert and Beheng, 2006)
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_graupel

       INTEGER(iwp) ::  i             !<
       INTEGER(iwp) ::  j             !<
       INTEGER(iwp) ::  k             !<

       REAL(wp) ::  av = 46.4_wp   !< parameter for calculating fall speed
       REAL(wp) ::  bv = 0.26_wp   !< constant in fall speed relation
       REAL(wp) ::  factor_sed_gamma_n = 0.88_wp !< corresponds to µ=1/3 and v=1 for gamma distribution
       REAL(wp) ::  factor_sed_gamma_q = 1.21_wp !< corresponds to µ=1/3 and v=1 for gamma distribution
       REAL(wp) ::  flag           !< flag to indicate first grid level
       REAL(wp) ::  xg = 0.0_wp    !< mean mass of graupel
       REAL(wp) ::  vg = 0.0_wp    !< mean fall speed of graupel

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_ng  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qg  !<

       sed_qg(nzt+1) = 0.0_wp
       sed_ng(nzt+1) = 0.0_wp

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzt, nzb+1, -1
                IF ( ng(k,j,i) <= 0.0_wp )  THEN
                   xg = 0.0_wp
                ELSE
!--                Calculate mean mass of ice crystal
                   xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
                ENDIF
!
!--             calculate fall speed of graupel
                vg = av * xg**bv
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Calculate sedimentation rate for each grid box, factors are
!--             calculated using
                IF ( qg(k,j,i) > eps_sb  .AND.  ng(k,j,i) >= 0.0_wp )  THEN
                   sed_qg(k) = qg(k,j,i) * vg * factor_sed_gamma_q * flag
                   sed_ng(k) = ng(k,j,i) * vg * factor_sed_gamma_n * flag
                ELSE
                   sed_qg(k) = 0.0_wp
                   sed_ng(k) = 0.0_wp
                ENDIF
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                sed_qg(k) = MIN( sed_qg(k), hyrho(k) * dzu(k+1) * qg(k,j,i) / dt_micro +           &
                                 sed_qg(k+1) ) * flag

                q(k,j,i)  = q(k,j,i)  + ( sed_qg(k+1) - sed_qg(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * dt_micro * flag
                qg(k,j,i) = qg(k,j,i) + ( sed_qg(k+1) - sed_qg(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * dt_micro * flag
                ng(k,j,i) = ng(k,j,i) + ( sed_ng(k+1) - sed_ng(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) - ( sed_qg(k+1) - sed_qg(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * l_s / c_p * d_exner(k) * dt_micro * flag
!
!--             Compute the precipitation rate of graupel
                IF ( call_microphysics_at_all_substeps )  THEN
                   prr(k,j,i) = prr(k,j,i) + sed_qg(k) / hyrho(k) *                                &
                                             weight_substep(intermediate_timestep_count) * flag
                   prr_graupel(k,j,i) = prr_graupel(k,j,i) + sed_qg(k) / hyrho(k) *                &
                                             weight_substep(intermediate_timestep_count) * flag
                ELSE
                   prr(k,j,i) = prr(k,j,i) + sed_qg(k) / hyrho(k) * flag
                   prr_graupel(k,j,i) = prr_graupel(k,j,i) + sed_qg(k) / hyrho(k) * flag
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE sedimentation_graupel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of graupel assuming generalized gamma distibution (Seifert and Beheng, 2006)
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_graupel_ij( i, j )

       INTEGER(iwp) ::  i             !<
       INTEGER(iwp) ::  j             !<
       INTEGER(iwp) ::  k             !<

       REAL(wp) ::  av = 46.4_wp   !< parameter for calculating fall speed
       REAL(wp) ::  bv = 0.26_wp   !< constant in fall speed relation
       REAL(wp) ::  factor_sed_gamma_n = 0.88_wp !< corresponds to µ=1/3 and v=1 for gamma distribution
       REAL(wp) ::  factor_sed_gamma_q = 1.21_wp !< corresponds to µ=1/3 and v=1 for gamma distribution
       REAL(wp) ::  flag           !< flag to indicate first grid level
       REAL(wp) ::  xg = 0.0_wp    !< mean mass of graupel
       REAL(wp) ::  vg = 0.0_wp    !< mean fall speed of graupel

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_ng  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qg  !<

       sed_qg(nzt+1) = 0.0_wp
       sed_ng(nzt+1) = 0.0_wp

       DO  k = nzt, nzb+1, -1
          IF ( ng(k,j,i) <= 0.0_wp )  THEN
             xg = 0.0_wp
          ELSE
!--          Calculate mean mass of ice crystal
             xg = mean_mass( cloud_species%graupel, qg(k,j,i), ng(k,j,i), hyrho(k) )
          ENDIF
!
!--       calculate fall speed of graupel
          vg = av * xg**bv
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Calculate sedimentation rate for each grid box, factors are
!--       calculated using
          IF ( qg(k,j,i) > eps_sb  .AND.  ng(k,j,i) >= 0.0_wp )  THEN
             sed_qg(k) = qg(k,j,i) * vg * factor_sed_gamma_q * flag
             sed_ng(k) = ng(k,j,i) * vg * factor_sed_gamma_n * flag
          ELSE
             sed_qg(k) = 0.0_wp
             sed_ng(k) = 0.0_wp
          ENDIF
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          sed_qg(k) = MIN( sed_qg(k), hyrho(k) * dzu(k+1) * qg(k,j,i) / dt_micro + sed_qg(k+1)     &
                         ) * flag

          q(k,j,i)  = q(k,j,i)  + ( sed_qg(k+1) - sed_qg(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          qg(k,j,i) = qg(k,j,i) + ( sed_qg(k+1) - sed_qg(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          ng(k,j,i) = ng(k,j,i) + ( sed_ng(k+1) - sed_ng(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          pt(k,j,i) = pt(k,j,i) - ( sed_qg(k+1) - sed_qg(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * l_s / c_p * d_exner(k) * dt_micro * flag
!
!--       Compute the precipitation rate of graupel
          IF ( call_microphysics_at_all_substeps )  THEN
             prr(k,j,i) = prr(k,j,i) + sed_qg(k) / hyrho(k) *                                      &
                                       weight_substep(intermediate_timestep_count) * flag
             prr_graupel(k,j,i) = prr_graupel(k,j,i) + sed_qg(k) / hyrho(k) *                      &
                                       weight_substep(intermediate_timestep_count) * flag
          ELSE
             prr(k,j,i) = prr(k,j,i) + sed_qg(k) / hyrho(k) * flag
             prr_graupel(k,j,i) = prr_graupel(k,j,i) + sed_qg(k) / hyrho(k) * flag
          ENDIF

       ENDDO

    END SUBROUTINE sedimentation_graupel_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of snow assuming generalized gamma distibution (Seifert and Beheng, 2006)
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_snow

       INTEGER(iwp) ::  i             !<
       INTEGER(iwp) ::  j             !<
       INTEGER(iwp) ::  k             !<

       REAL(wp) ::  av = 27.7_wp   !< parameter for calculating fall speed
       REAL(wp) ::  bv = 0.22_wp   !< constant in fall speed relation
       REAL(wp) ::  factor_sed_gamma_n = 0.89_wp !< corresponds to µ=1/3 and v=1 for gamma distribution
       REAL(wp) ::  factor_sed_gamma_q = 1.17_wp !< corresponds to µ=1/3 and v=1 for gamma distribution
       REAL(wp) ::  flag           !< flag to indicate first grid level
       REAL(wp) ::  xs = 0.0_wp    !< mean mass of snow
       REAL(wp) ::  vs = 0.0_wp    !< mean fall speed of snow

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_ns  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qs  !<

       sed_qs(nzt+1) = 0.0_wp
       sed_ns(nzt+1) = 0.0_wp

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzt, nzb+1, -1
                IF ( ns(k,j,i) <= 0.0_wp )  THEN
                   xs = 0.0_wp
                ELSE
!--                Calculate mean mass snow
                   xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )
                ENDIF
!
!--             calculate fall speed of snow
                vs = av * xs**bv
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Calculate sedimentation rate for each grid box, factors are
!--             calculated using
                IF ( qs(k,j,i) > eps_sb  .AND.  ns(k,j,i) >= 0.0_wp )  THEN
                   sed_qs(k) = qs(k,j,i) * vs * factor_sed_gamma_q * flag
                   sed_ns(k) = ns(k,j,i) * vs * factor_sed_gamma_n * flag
                ELSE
                   sed_qs(k) = 0.0_wp
                   sed_ns(k) = 0.0_wp
                ENDIF
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                sed_qs(k) = MIN( sed_qs(k), hyrho(k) * dzu(k+1) * qs(k,j,i) / dt_micro +           &
                                 sed_qs(k+1) ) * flag

                q(k,j,i)  = q(k,j,i)  + ( sed_qs(k+1) - sed_qs(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * dt_micro * flag
                qs(k,j,i) = qs(k,j,i) + ( sed_qs(k+1) - sed_qs(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * dt_micro * flag
                ns(k,j,i) = ns(k,j,i) + ( sed_ns(k+1) - sed_ns(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) - ( sed_qs(k+1) - sed_qs(k) ) * ddzu(k+1) /                  &
                                        hyrho(k) * l_s / c_p * d_exner(k) * dt_micro * flag
!
!--             Compute the precipitation rate of snow
                IF ( call_microphysics_at_all_substeps )  THEN
                   prr(k,j,i) = prr(k,j,i) + sed_qs(k) / hyrho(k) *                                &
                                             weight_substep(intermediate_timestep_count) * flag
                   prr_snow(k,j,i) = prr_snow(k,j,i) + sed_qs(k) / hyrho(k) *                      &
                                             weight_substep(intermediate_timestep_count) * flag
                ELSE
                   prr(k,j,i) = prr(k,j,i) + sed_qs(k) / hyrho(k) * flag
                   prr_snow(k,j,i) = prr_snow(k,j,i) + sed_qs(k) / hyrho(k) * flag
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE sedimentation_snow


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sedimentation of snow assuming generalized gamma distibution (Seifert and Beheng, 2006)
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_snow_ij( i, j )

       INTEGER(iwp) ::  i             !<
       INTEGER(iwp) ::  j             !<
       INTEGER(iwp) ::  k             !<

       REAL(wp) ::  av = 27.7_wp   !< parameter for calculating fall speed
       REAL(wp) ::  bv = 0.22_wp   !< constant in fall speed relation
       REAL(wp) ::  factor_sed_gamma_n = 0.89_wp !< corresponds to µ=1/3 and v=1 for gamma distribution
       REAL(wp) ::  factor_sed_gamma_q = 1.17_wp !< corresponds to µ=1/3 and v=1 for gamma distribution
       REAL(wp) ::  flag           !< flag to indicate first grid level
       REAL(wp) ::  xs = 0.0_wp    !< mean mass of snow
       REAL(wp) ::  vs = 0.0_wp    !< mean fall speed of snow

       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_ns  !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qs  !<

       sed_qs(nzt+1) = 0.0_wp
       sed_ns(nzt+1) = 0.0_wp

       DO  k = nzt, nzb+1, -1
          IF ( ns(k,j,i) <= 0.0_wp )  THEN
             xs = 0.0_wp
          ELSE
!--          Calculate mean mass of ice crystal
             xs = mean_mass( cloud_species%snow, qs(k,j,i), ns(k,j,i), hyrho(k) )
          ENDIF
!
!--       calculate fall speed of ice crystal
          vs = av * xs**bv
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Calculate sedimentation rate for each grid box, factors are
!--       calculated using
          IF ( qs(k,j,i) > eps_sb  .AND.  ns(k,j,i) > 0.0_wp )  THEN
             sed_qs(k) = qs(k,j,i) * vs * factor_sed_gamma_q * flag
             sed_ns(k) = ns(k,j,i) * vs * factor_sed_gamma_n * flag
          ELSE
             sed_qs(k) = 0.0_wp
             sed_ns(k) = 0.0_wp
          ENDIF
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          sed_qs(k) = MIN( sed_qs(k), hyrho(k) * dzu(k+1) * qs(k,j,i) / dt_micro + sed_qs(k+1)     &
                         ) * flag

          q(k,j,i)  = q(k,j,i)  + ( sed_qs(k+1) - sed_qs(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          qs(k,j,i) = qs(k,j,i) + ( sed_qs(k+1) - sed_qs(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          ns(k,j,i) = ns(k,j,i) + ( sed_ns(k+1) - sed_ns(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          pt(k,j,i) = pt(k,j,i) - ( sed_qs(k+1) - sed_qs(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * l_s / c_p * d_exner(k) * dt_micro * flag
!
!--       Compute the precipitation rate of cloud (fog) droplets
          IF ( call_microphysics_at_all_substeps )  THEN
             prr(k,j,i) = prr(k,j,i) + sed_qs(k) / hyrho(k) *                                      &
                                       weight_substep(intermediate_timestep_count) * flag
             prr_snow(k,j,i) = prr_snow(k,j,i) + sed_qs(k) / hyrho(k) *                            &
                                       weight_substep(intermediate_timestep_count) * flag
          ELSE
             prr(k,j,i) = prr(k,j,i) + sed_qs(k) / hyrho(k) * flag
             prr_snow(k,j,i) = prr_snow(k,j,i) + sed_qs(k) / hyrho(k) * flag
          ENDIF

       ENDDO

    END SUBROUTINE sedimentation_snow_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of sedimentation flux. Implementation according to Stevens
!> and Seifert (2008). Code is based on UCLA-LES.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_rain

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< loop index
       INTEGER(iwp) ::  j             !< loop index
       INTEGER(iwp) ::  k             !< loop index
       INTEGER(iwp) ::  k_run         !<

       REAL(wp)     ::  c_run                      !<
       REAL(wp)     ::  d_max                      !<
       REAL(wp)     ::  d_mean                     !<
       REAL(wp)     ::  d_min                      !<
       REAL(wp)     ::  dr                         !<
       REAL(wp)     ::  flux                       !<
       REAL(wp)     ::  flag                       !< flag to mask topography grid points
       REAL(wp)     ::  lambda_r                   !<
       REAL(wp)     ::  mu_r                       !<
       REAL(wp)     ::  z_run                      !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_qr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  nr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  qr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_qr     !<
!
!--    Compute velocities
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                !IF ( qr(k,j,i) > eps_sb  .AND.  nr(k,j,i) > 0.0_wp )  THEN
                IF ( qr(k,j,i) > eps_sb )  THEN

!
!--                Weight averaged diameter of rain drops:
                   dr = ( hyrho(k) * qr(k,j,i) / nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--                Shape parameter of gamma distribution (Milbrandt and Yau, 2005;
!--                Stevens and Seifert, 2008):
                   mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp * ( dr - 1.4E-3_wp ) ) )
!
!--                Slope parameter of gamma distribution (Seifert, 2008):
                   lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *                            &
                                ( mu_r + 1.0_wp ) )**( 1.0_wp / 3.0_wp ) / dr

                   w_nr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                                            &
                                               a_term - b_term * ( 1.0_wp +                        &
                                                  c_term /                                         &
                                                  lambda_r )**( -1.0_wp *                          &
                                                  ( mu_r + 1.0_wp ) )                              &
                                              )                                                    &
                                ) * flag

                   w_qr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                                            &
                                               a_term - b_term * ( 1.0_wp +                        &
                                                  c_term /                                         &
                                                  lambda_r )**( -1.0_wp *                          &
                                                  ( mu_r + 4.0_wp ) )                              &
                                             )                                                     &
                                ) * flag
                ELSE
                   w_nr(k) = 0.0_wp
                   w_qr(k) = 0.0_wp
                ENDIF
             ENDDO
!
!--          Adjust boundary values.
             k = topo_top_ind(j,i,0) + 1
             w_nr(k-1) = w_nr(k)
             w_qr(k-1) = w_qr(k)
!
!--          Model top boundary value
             w_nr(nzt+1) = 0.0_wp
             w_qr(nzt+1) = 0.0_wp
!
!--          Compute Courant number
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                c_nr(k) = 0.25_wp * ( w_nr(k-1) + 2.0_wp * w_nr(k) + w_nr(k+1) ) *                 &
                          dt_micro * ddzu(k) * flag
                c_qr(k) = 0.25_wp * ( w_qr(k-1) + 2.0_wp * w_qr(k) + w_qr(k+1) ) *                 &
                          dt_micro * ddzu(k) * flag
             ENDDO
!
!--          Limit slopes with monotonized centered (MC) limiter (van Leer, 1977):
             IF ( limiter_sedimentation )  THEN

                DO k = nzb+1, nzt
!
!--                Predetermine flag to mask topography
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                   d_mean = 0.5_wp * ( qr(k+1,j,i) - qr(k-1,j,i) )
                   d_min  = qr(k,j,i) - MIN( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) )
                   d_max  = MAX( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) ) - qr(k,j,i)

                   qr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min, 2.0_wp * d_max,      &
                                                              ABS( d_mean ) )                      &
                                                      * flag

                   d_mean = 0.5_wp * ( nr(k+1,j,i) - nr(k-1,j,i) )
                   d_min  = nr(k,j,i) - MIN( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) )
                   d_max  = MAX( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) ) - nr(k,j,i)

                   nr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min, 2.0_wp * d_max,      &
                                                              ABS( d_mean ) )
                ENDDO

             ELSE

                nr_slope = 0.0_wp
                qr_slope = 0.0_wp

             ENDIF

             sed_nr(nzt+1) = 0.0_wp
             sed_qr(nzt+1) = 0.0_wp
!
!--          Compute sedimentation flux
             DO  k = nzt, nzb+1, -1
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Sum up all rain drop number densities which contribute to the flux
!--             through k-1/2
                flux  = 0.0_wp
                z_run = 0.0_wp ! height above z(k)
                k_run = k
                c_run = MIN( 1.0_wp, c_nr(k) )
                DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )
                   flux  = flux + hyrho(k_run) * ( nr(k_run,j,i) + nr_slope(k_run) *               &
                           ( 1.0_wp - c_run ) * 0.5_wp ) * c_run * dzu(k_run) * flag
                   z_run = z_run + dzu(k_run) * flag
                   k_run = k_run + 1 * flag
                   c_run = MIN( 1.0_wp, c_nr(k_run) - z_run * ddzu(k_run) ) * flag
                ENDDO
!
!--             It is not allowed to sediment more rain drop number density than available.
                flux = MIN( flux, hyrho(k) * dzu(k+1) * nr(k,j,i) + sed_nr(k+1) * dt_micro )

                sed_nr(k) = flux / dt_micro * flag
                nr(k,j,i) = nr(k,j,i) + ( sed_nr(k+1) - sed_nr(k) ) *                              &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
!
!--             Sum up all rain water content which contributes to the flux through k-1/2
                flux  = 0.0_wp
                z_run = 0.0_wp ! height above z(k)
                k_run = k
                c_run = MIN( 1.0_wp, c_qr(k) )

                DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )

                   flux  = flux +                                                                  &
                           hyrho(k_run) *                                                          &
                           ( qr(k_run,j,i) + qr_slope(k_run) * ( 1.0_wp - c_run ) * 0.5_wp ) *     &
                           c_run * dzu(k_run) * flag
                   z_run = z_run + dzu(k_run) * flag
                   k_run = k_run + 1 * flag
                   c_run = MIN( 1.0_wp, c_qr(k_run) - z_run * ddzu(k_run) ) * flag

                ENDDO
!
!--             It is not allowed to sediment more rain water content than available.
                flux = MIN( flux, hyrho(k) * dzu(k) * qr(k,j,i) + sed_qr(k+1) * dt_micro )

                sed_qr(k) = flux / dt_micro * flag

                qr(k,j,i) = qr(k,j,i) + ( sed_qr(k+1) - sed_qr(k) ) *                              &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                q(k,j,i)  = q(k,j,i)  + ( sed_qr(k+1) - sed_qr(k) ) *                              &
                                        ddzu(k+1) / hyrho(k) * dt_micro * flag
                pt(k,j,i) = pt(k,j,i) - ( sed_qr(k+1) - sed_qr(k) ) *                              &
                                        ddzu(k+1) / hyrho(k) * lv_d_cp *                           &
                                        d_exner(k) * dt_micro * flag
!
!--             Compute the rain rate
                IF ( call_microphysics_at_all_substeps )  THEN
                   prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k) *                                &
                                             weight_substep(intermediate_timestep_count) * flag
                   prr_rain(k,j,i) = prr_rain(k,j,i) + sed_qr(k) / hyrho(k) *                      &
                                             weight_substep(intermediate_timestep_count) * flag
                ELSE
                   prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k) * flag
                   prr_rain(k,j,i) = prr_rain(k,j,i) + sed_qr(k) / hyrho(k) * flag
                ENDIF

             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE sedimentation_rain


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of sedimentation flux. Implementation according to Stevens
!> and Seifert (2008). Code is based on UCLA-LES. Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE sedimentation_rain_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< loop index
       INTEGER(iwp) ::  j             !< loop index
       INTEGER(iwp) ::  k             !< loop index
       INTEGER(iwp) ::  k_run         !<

       REAL(wp)     ::  c_run                     !<
       REAL(wp)     ::  d_max                     !<
       REAL(wp)     ::  d_mean                    !<
       REAL(wp)     ::  d_min                     !<
       REAL(wp)     ::  dr                        !<
       REAL(wp)     ::  flux                      !<
       REAL(wp)     ::  flag                      !< flag to indicate first grid level above surface
       REAL(wp)     ::  lambda_r                  !<
       REAL(wp)     ::  mu_r                      !<
       REAL(wp)     ::  z_run                     !<

       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  c_qr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  nr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  qr_slope !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_nr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  sed_qr   !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_nr     !<
       REAL(wp), DIMENSION(nzb:nzt+1) ::  w_qr     !<

!
!--    Compute velocities
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

!          IF ( qr(k,j,i) > eps_sb  .AND.  nr(k,j,i) > 0.0_wp )  THEN
           IF ( qr(k,j,i) > eps_sb )  THEN

!
!--          Weight averaged diameter of rain drops:
             dr = ( hyrho(k) * qr(k,j,i) / nr(k,j,i) * dpirho_l )**( 1.0_wp / 3.0_wp )
!
!--          Shape parameter of gamma distribution (Milbrandt and Yau, 2005;
!--          Stevens and Seifert, 2008):
             mu_r = 10.0_wp * ( 1.0_wp + TANH( 1.2E3_wp * ( dr - 1.4E-3_wp ) ) )
!
!--          Slope parameter of gamma distribution (Seifert, 2008):
             lambda_r = ( ( mu_r + 3.0_wp ) * ( mu_r + 2.0_wp ) *                                  &
                          ( mu_r + 1.0_wp ) )**( 1.0_wp / 3.0_wp ) / dr

             w_nr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                                                  &
                                         a_term - b_term * ( 1.0_wp + c_term / lambda_r            &
                                                           )**( -1.0_wp * ( mu_r + 1.0_wp ) )      &
                                       )                                                           &
                          ) * flag
             w_qr(k) = MAX( 0.1_wp, MIN( 20.0_wp,                                                  &
                                         a_term - b_term * ( 1.0_wp + c_term / lambda_r            &
                                                           )**( -1.0_wp * ( mu_r + 4.0_wp ) )      &
                                       )                                                           &
                          ) * flag
          ELSE
             w_nr(k) = 0.0_wp
             w_qr(k) = 0.0_wp
          ENDIF
       ENDDO
!
!--    Adjust boundary values.
       k = topo_top_ind(j,i,0) + 1
       w_nr(k-1) = w_nr(k)
       w_qr(k-1) = w_qr(k)
!
!--    Neumann boundary condition at model top
       w_nr(nzt+1) = 0.0_wp
       w_qr(nzt+1) = 0.0_wp
!
!--    Compute Courant number
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          c_nr(k) = 0.25_wp * ( w_nr(k-1) + 2.0_wp * w_nr(k) + w_nr(k+1) ) *                       &
                    dt_micro * ddzu(k) * flag
          c_qr(k) = 0.25_wp * ( w_qr(k-1) + 2.0_wp * w_qr(k) + w_qr(k+1) ) *                       &
                    dt_micro * ddzu(k) * flag
       ENDDO
!
!--    Limit slopes with monotonized centered (MC) limiter (van Leer, 1977):
       IF ( limiter_sedimentation )  THEN

          DO k = nzb+1, nzt
!
!--          Predetermine flag to mask topography
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

             d_mean = 0.5_wp * ( qr(k+1,j,i) - qr(k-1,j,i) )
             d_min  = qr(k,j,i) - MIN( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) )
             d_max  = MAX( qr(k+1,j,i), qr(k,j,i), qr(k-1,j,i) ) - qr(k,j,i)

             qr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min, 2.0_wp * d_max,            &
                                                        ABS( d_mean ) ) * flag

             d_mean = 0.5_wp * ( nr(k+1,j,i) - nr(k-1,j,i) )
             d_min  = nr(k,j,i) - MIN( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) )
             d_max  = MAX( nr(k+1,j,i), nr(k,j,i), nr(k-1,j,i) ) - nr(k,j,i)

             nr_slope(k) = SIGN(1.0_wp, d_mean) * MIN ( 2.0_wp * d_min,  2.0_wp * d_max,           &
                                                        ABS( d_mean ) ) * flag
          ENDDO

       ELSE

          nr_slope = 0.0_wp
          qr_slope = 0.0_wp

       ENDIF

       sed_nr(nzt+1) = 0.0_wp
       sed_qr(nzt+1) = 0.0_wp
!
!--    Compute sedimentation flux
       DO  k = nzt, nzb+1, -1
!
!--       Predetermine flag to mask topography
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Sum up all rain drop number densities which contribute to the flux through k-1/2
          flux  = 0.0_wp
          z_run = 0.0_wp ! height above z(k)
          k_run = k
          c_run = MIN( 1.0_wp, c_nr(k) )
          DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )
             flux  = flux + hyrho(k_run) *                                                         &
                     ( nr(k_run,j,i) + nr_slope(k_run) * ( 1.0_wp - c_run ) * 0.5_wp ) *           &
                     c_run * dzu(k_run) * flag
             z_run = z_run + dzu(k_run) * flag
             k_run = k_run + 1 * flag
             c_run = MIN( 1.0_wp, c_nr(k_run) - z_run * ddzu(k_run) ) * flag
          ENDDO
!
!--       It is not allowed to sediment more rain drop number density than available.
          flux = MIN( flux, hyrho(k) * dzu(k+1) * nr(k,j,i) + sed_nr(k+1) * dt_micro )

          sed_nr(k) = flux / dt_micro * flag
          nr(k,j,i)  = nr(k,j,i) + ( sed_nr(k+1) - sed_nr(k) ) * ddzu(k+1) /                       &
                                    hyrho(k) * dt_micro * flag
!
!--       Sum up all rain water content which contributes to the flux
!--       through k-1/2
          flux  = 0.0_wp
          z_run = 0.0_wp ! height above z(k)
          k_run = k
          c_run = MIN( 1.0_wp, c_qr(k) )

          DO WHILE ( c_run > 0.0_wp  .AND.  k_run <= nzt )

             flux  = flux + hyrho(k_run) *                                                         &
                     ( qr(k_run,j,i) + qr_slope(k_run) * ( 1.0_wp - c_run ) * 0.5_wp ) *           &
                     c_run * dzu(k_run) * flag
             z_run = z_run + dzu(k_run) * flag
             k_run = k_run + 1 * flag
             c_run = MIN( 1.0_wp, c_qr(k_run) - z_run * ddzu(k_run) ) * flag

          ENDDO
!
!--       It is not allowed to sediment more rain water content than available.
          flux = MIN( flux, hyrho(k) * dzu(k) * qr(k,j,i) + sed_qr(k+1) * dt_micro )

          sed_qr(k) = flux / dt_micro * flag

          qr(k,j,i) = qr(k,j,i) + ( sed_qr(k+1) - sed_qr(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          q(k,j,i)  = q(k,j,i)  + ( sed_qr(k+1) - sed_qr(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * dt_micro * flag
          pt(k,j,i) = pt(k,j,i) - ( sed_qr(k+1) - sed_qr(k) ) * ddzu(k+1) /                        &
                                  hyrho(k) * lv_d_cp * d_exner(k) * dt_micro * flag
!
!--       Compute the rain rate
          IF ( call_microphysics_at_all_substeps )  THEN
             prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k) *                                      &
                                       weight_substep(intermediate_timestep_count) * flag
             prr_rain(k,j,i) = prr_rain(k,j,i) + sed_qr(k) / hyrho(k) *                            &
                                       weight_substep(intermediate_timestep_count) * flag
          ELSE
             prr(k,j,i) = prr(k,j,i) + sed_qr(k) / hyrho(k) * flag
             prr_rain(k,j,i) = prr_rain(k,j,i) + sed_qr(k) / hyrho(k) * flag
          ENDIF

       ENDDO

    END SUBROUTINE sedimentation_rain_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the precipitation amount due to gravitational settling of rain and cloud (fog)
!> droplets
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE calc_precipitation_amount

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< loop index
       INTEGER(iwp) ::  j             !< loop index
       INTEGER(iwp) ::  k             !< loop index
       INTEGER(iwp) ::  m             !< running index surface elements

       IF ( ( dt_do2d_xy - time_do2d_xy ) < precipitation_amount_interval .AND.                    &
            ( .NOT. call_microphysics_at_all_substeps .OR.                                         &
              intermediate_timestep_count == intermediate_timestep_count_max ) )                   &
       THEN
!
!--       Run over all upward-facing surface elements koff(m) = -1, i.e. non-natural, natural and
!--       urban.
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m)
             j = bc_hv%j(m)
             k = bc_hv%k(m)
             precipitation_amount(j,i) = MERGE( precipitation_amount(j,i)                          &
                                              + prr(k,j,i) * hyrho(k) * dt_3d,                     &
                                                precipitation_amount(j,i),                         &
                                                bc_hv%koff(m) == -1 )
          ENDDO

       ENDIF

    END SUBROUTINE calc_precipitation_amount


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine computes the precipitation amount due to gravitational settling of rain and cloud
!> (fog) droplets
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE calc_precipitation_amount_ij( i, j )

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< loop index
       INTEGER(iwp) ::  j             !< loop index
       INTEGER(iwp) ::  k             !< loop index

       IF ( ( dt_do2d_xy - time_do2d_xy ) < precipitation_amount_interval .AND.                    &
            ( .NOT. call_microphysics_at_all_substeps .OR.                                         &
              intermediate_timestep_count == intermediate_timestep_count_max ) )                   &
       THEN
!
!--       Note, in contrast to the vector-optimized version, computation of precipitaton amount
!--       at the surface is done using the topo_top index instead of the surface data structure.
!--       This is because no local information from only horizontally upward facing surfaces can
!--       be inferred from bc_hv any more after putting all surfaces into one array.
          k = topo_top_ind(j,i,0) + 1
          precipitation_amount(j,i) = precipitation_amount(j,i) + prr(k,j,i) * hyrho(k) * dt_3d

       ENDIF

    END SUBROUTINE calc_precipitation_amount_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the diagnostic supersaturation sat, actual liquid water temperature t_l and
!> saturation water vapor mixing ratio q_s
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE supersaturation ( i, j, k )

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< loop index
       INTEGER(iwp) ::  j   !< loop index
       INTEGER(iwp) ::  k   !< loop index

       REAL(wp) ::  alpha   !< correction factor
!
!--    Actual liquid water temperature:
       t_l = exner(k) * pt(k,j,i)
!
!--    Calculate water vapor saturation pressure
       e_s = magnus_tl( t_l )
!
!--    Computation of saturation mixing ratio:
       q_s   = rd_d_rv * e_s / ( hyp(k) - e_s )
!
!--    Correction factor
       alpha = rd_d_rv * lv_d_rd * lv_d_cp / ( t_l * t_l )
!
!--    Correction of the approximated value (see: Cuijpers + Duynkerke, 1993, JAS, 23)
       q_s   = q_s * ( 1.0_wp + alpha * q(k,j,i) ) / ( 1.0_wp + alpha * q_s )
!
!--    Only liquid phase
       IF ( .NOT. microphysics_ice_phase )  THEN
          IF ( microphysics_seifert )  THEN
             sat = ( q(k,j,i) - qr(k,j,i) - qc(k,j,i) ) / q_s - 1.0_wp
          ELSEIF ( microphysics_morrison_no_rain )  THEN
             sat = ( q(k,j,i) - qc(k,j,i) ) / q_s - 1.0_wp
          ENDIF
!
!--    ice phase
       ELSE
          IF ( microphysics_seifert .AND.  .NOT. graupel  .OR.  .NOT. snow )  THEN
             sat = ( q(k,j,i) - qr(k,j,i) - qc(k,j,i) - qi(k,j,i) ) / q_s - 1.0_wp
!
!--       IF snow and graupel is also turned on
          ELSE
             sat = ( q(k,j,i) - qc(k,j,i) - qr(k,j,i) - qi(k,j,i) - qs(k,j,i) - qg(k,j,i) ) / &
                   q_s - 1.0_wp
          ENDIF
       ENDIF

    END SUBROUTINE supersaturation

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the diagnostic supersaturation sat, actual liquid water temperature t_l and
!> saturation water vapor mixing ratio q_s
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE supersaturation_ice ( i, j, k )

       INTEGER(iwp) ::  i   !< loop index
       INTEGER(iwp) ::  j   !< loop index
       INTEGER(iwp) ::  k   !< loop index

       REAL(wp) ::  e_a     !< water vapor pressure
       REAL(wp) ::  temp    !< actual temperature

!
!--    Actual liquid water temperature:
       t_l = exner(k) * pt(k,j,i)
!
!--    Actual temperature:
       temp = t_l + lv_d_cp * ql(k,j,i) + ls_d_cp * qf(k,j,i)
!
!--    Calculate water vapor saturation pressure with formular using actual temperature
       e_si = magnus_ice( temp )
!
!--    Computation of ice saturation mixing ratio:
       q_si = rd_d_rv * e_si / ( hyp(k) - e_si )
!
!--    Current vapor pressure
       IF ( .NOT. snow  .AND.  .NOT. graupel )  THEN
          e_a = (   q(k,j,i) - qr(k,j,i) - qc(k,j,i) - qi(k,j,i) ) * hyp(k) /                      &
                ( ( q(k,j,i) - qr(k,j,i) - qc(k,j,i) - qi(k,j,i) ) + rd_d_rv )
       ELSE
          e_a = ( q(k,j,i) - qr(k,j,i) - qc(k,j,i) - qi(k,j,i) - qg(k,j,i) - qs(k,j,i) ) * hyp(k) /&
                ( ( q(k,j,i) - qr(k,j,i) - qc(k,j,i) - qi(k,j,i) - qg(k,j,i) - qs(k,j,i) ) +       &
                   rd_d_rv )
       ENDIF
!
!--    Supersaturation:
!--    Not in case of microphysics_kessler or microphysics_sat_adjust since qr is unallocated
       sat_ice = e_a / e_si - 1.0_wp
!
!--    If temperature of ice is over 273.15°C set sat_ice <= 0.0
       IF ( temp > t3  .AND.  sat_ice > 0.0_wp  ) THEN
          sat_ice = 0.0_wp
       ENDIF

    END SUBROUTINE supersaturation_ice


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the liquid water content (0%-or-100%-scheme). This scheme is used by the one and
!> the two moment cloud physics scheme. Using the two moment scheme, this calculation results in the
!> cloud water content.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE calc_liquid_water_content

       INTEGER(iwp) ::  i !< loop index
       INTEGER(iwp) ::  j !< loop index
       INTEGER(iwp) ::  k !< loop index

       REAL(wp)     ::  flag !< flag to indicate first grid level above surface

       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Call calculation of supersaturation located
                CALL supersaturation( i, j, k )
!
!--             Compute the liquid water content
                IF ( .NOT. microphysics_ice_phase ) THEN
                   IF ( microphysics_seifert  .AND. .NOT. microphysics_morrison )  THEN
!--                   Seifert and Beheng scheme: saturation adjustment
                      IF ( ( q(k,j,i) - q_s - qr(k,j,i) ) > 0.0_wp )  THEN
                         qc(k,j,i) = ( q(k,j,i) - q_s - qr(k,j,i) ) * flag
                         ql(k,j,i) = ( qc(k,j,i) + qr(k,j,i) ) * flag
                      ELSE
                         IF ( q(k,j,i) < qr(k,j,i) )  q(k,j,i) = qr(k,j,i)
                         qc(k,j,i) = 0.0_wp
                         ql(k,j,i) = qr(k,j,i) * flag
                      ENDIF
!
!--                Morrison scheme: explicit condensation (see above)
                   ELSEIF ( microphysics_morrison  .AND.  microphysics_seifert )  THEN
                      ql(k,j,i) = qc(k,j,i) + qr(k,j,i) * flag
!
!--                Morrison without rain: explicit condensation
                   ELSEIF ( microphysics_morrison  .AND. .NOT. microphysics_seifert )  THEN
                      ql(k,j,i) = qc(k,j,i)
!
!--                Kessler and saturation adjustment scheme
                   ELSE
                      IF ( ( q(k,j,i) - q_s ) > 0.0_wp )  THEN
                         qc(k,j,i) = ( q(k,j,i) - q_s ) * flag
                         ql(k,j,i) = qc(k,j,i) * flag
                      ELSE
                         qc(k,j,i) = 0.0_wp
                         ql(k,j,i) = 0.0_wp
                      ENDIF
                   ENDIF
!--             Calculations of liquid water content in case of mixed-phase cloud microphyics
                ELSE
                   IF ( microphysics_seifert  .AND. .NOT. microphysics_morrison ) &
                   THEN
!
!--                   Seifert and Beheng scheme: saturation adjustment if mixed phase but no snow
!--                   and graupel
                      IF ( .NOT. snow  .AND.  .NOT. graupel ) THEN
                         IF ( ( q(k,j,i) - q_s - qr(k,j,i) - qi(k,j,i) ) > 0.0_wp )  THEN
                            qc(k,j,i) = ( q(k,j,i) - q_s - qr(k,j,i) - qi(k,j,i) ) * flag
                            ql(k,j,i) = ( qc(k,j,i) + qr(k,j,i) ) * flag
                            qf(k,j,i) = qi(k,j,i)
                         ELSE
                            IF ( q(k,j,i) < ( qr(k,j,i) + qi(k,j,i) ) )  THEN
                               q(k,j,i) = qr(k,j,i) + qi(k,j,i)
                            ENDIF
                            qc(k,j,i) = 0.0_wp
                            ql(k,j,i) = qr(k,j,i) * flag
                            qf(k,j,i) = qi(k,j,i) * flag
                         ENDIF
!
!--                   Seifert and Beheng scheme: saturation adjustment if mixed phase but with snow
!--                   and graupel. Note: qs = snow mixing ratio, q_s = saturation mixing ratio
                      ELSE
                         IF ( ( q(k,j,i) - q_s - qr(k,j,i) - qi(k,j,i) - qs(k,j,i) - qg(k,j,i) )   &
                                > 0.0_wp )  THEN
                            qc(k,j,i) = ( q(k,j,i) - q_s - qr(k,j,i) - qi(k,j,i) - qs(k,j,i)       &
                                                   - qg(k,j,i) ) * flag
                            ql(k,j,i) = ( qc(k,j,i) + qr(k,j,i) ) * flag
                            qf(k,j,i) = qi(k,j,i) + qs(k,j,i) + qg(k,j,i)
                         ELSE
                            IF ( q(k,j,i) < ( qr(k,j,i) + qi(k,j,i) + qs(k,j,i) + qg(k,j,i) ) )    &
                            THEN
                               q(k,j,i) = qr(k,j,i) + qi(k,j,i) + qs(k,j,i) + qg(k,j,i)
                            ENDIF
                            qc(k,j,i) = 0.0_wp
                            ql(k,j,i) = qr(k,j,i) * flag
                            qf(k,j,i) = qi(k,j,i) + qs(k,j,i) + qg(k,j,i) * flag
                         ENDIF
                      ENDIF
!
!--                Morrison scheme: explicit condensation (see above)
                   ELSEIF ( microphysics_morrison  .AND.  microphysics_seifert )  THEN
                      ql(k,j,i) = qc(k,j,i) + qr(k,j,i) * flag
                      IF ( snow  .AND.  graupel )  THEN
                         qf(k,j,i) = qi(k,j,i) + qs(k,j,i) + qg(k,j,i) * flag
                      ELSE
                         qf(k,j,i) = qi(k,j,i) * flag
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    END SUBROUTINE calc_liquid_water_content


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the gamma function (Press et al., 1992).
!> The gamma function is needed for the calculation of the evaporation of rain drops.
!--------------------------------------------------------------------------------------------------!
    FUNCTION gamm( xx )

       IMPLICIT NONE

       INTEGER(iwp) ::  j            !<

       REAL(wp), PARAMETER  ::  stp = 2.5066282746310005_wp               !<
       REAL(wp), PARAMETER  ::  cof(6) = (/ 76.18009172947146_wp,                                  &
                                           -86.50532032941677_wp,                                  &
                                            24.01409824083091_wp,                                  &
                                            -1.231739572450155_wp,                                 &
                                             0.1208650973866179E-2_wp,                             &
                                            -0.5395239384953E-5_wp /)     !<

       REAL(wp)     ::  gamm         !<
       REAL(wp)     ::  ser          !<
       REAL(wp)     ::  tmp          !<
       REAL(wp)     ::  x_gamm       !<
       REAL(wp)     ::  xx           !<
       REAL(wp)     ::  y_gamm       !<


       x_gamm = xx
       y_gamm = x_gamm
       tmp = x_gamm + 5.5_wp
       tmp = ( x_gamm + 0.5_wp ) * LOG( tmp ) - tmp
       ser = 1.000000000190015_wp

       DO  j = 1, 6
          y_gamm = y_gamm + 1.0_wp
          ser    = ser + cof( j ) / y_gamm
       ENDDO

!
!--    Until this point the algorithm computes the logarithm of the gamma function. Hence, the
!--    exponential function is used.
!      gamm = EXP( tmp + LOG( stp * ser / x_gamm ) )
       gamm = EXP( tmp ) * stp * ser / x_gamm

       RETURN

    END FUNCTION gamm


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the complete mass moment of particle size distribution, Eq (82) of SB2006
!--------------------------------------------------------------------------------------------------!
    FUNCTION moment_gamma(specie, n)

       INTEGER, INTENT(in) ::  n                         !<
       CLASS(cloud_coefficients), INTENT(IN) ::  specie  !<
       REAL(wp) ::  moment_gamma                         !<

       moment_gamma = GAMMA( ( n + specie%nu + 1.0_wp ) / specie%mu ) /                            &
                      GAMMA( ( specie%nu + 1.0_wp ) / specie%mu )     *                            &
                    ( GAMMA( ( specie%nu + 1.0_wp ) / specie%mu )     /                            &
                      GAMMA( ( specie%nu + 2.0_wp ) / specie%mu )                                  &
                    )**n

  END FUNCTION moment_gamma


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the collision integrals as given in SB2006, EQ.90
!--------------------------------------------------------------------------------------------------!
    FUNCTION collision_integral_delta_b(specie, n)

       INTEGER(iwp), INTENT(IN) :: n
       CLASS(cloud_coefficients), INTENT(IN) :: specie
       REAL(wp) ::  collision_integral_delta_b
!
!--    Eq. 90 from SB2006
       collision_integral_delta_b = GAMMA( ( 2.0_wp * specie%b + specie%nu + 1.0_wp + n ) /         &
                                             specie%mu ) /                                          &
                            GAMMA( ( specie%nu + 1.0_wp ) / specie%mu ) *                           &
                            GAMMA( ( specie%nu + 1.0_wp ) / specie%mu )**(2.0_wp * specie%b + n ) / &
                            GAMMA( ( specie%nu + 2.0_wp ) / specie%mu )**(2.0_wp * specie%b + n)
       RETURN

     END FUNCTION collision_integral_delta_b


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the collision integrals as given in SB2006, EQ.91
!--------------------------------------------------------------------------------------------------!
    FUNCTION collision_integral_delta_ba(specie1, specie2, n)

       INTEGER(iwp), INTENT(IN) :: n
       CLASS(cloud_coefficients), INTENT(IN) :: specie1
       CLASS(cloud_coefficients), INTENT(IN) :: specie2
       REAL(wp) ::  collision_integral_delta_ba
!
!--    Eq. 91 from SB2006
       collision_integral_delta_ba = 2.0_wp *                                                      &
            GAMMA( ( specie1%b + specie1%nu + 1.0_wp ) / specie1%mu ) /                            &
            GAMMA( ( specie1%nu + 1.0_wp ) / specie1%mu )              *                           &
            GAMMA( ( specie1%nu + 1.0_wp ) / specie1%mu )**(specie1%b) /                           &
            GAMMA( ( specie1%nu + 2.0_wp ) / specie1%mu )**(specie1%b) *                           &
            GAMMA( ( specie2%b + specie2%nu + 1.0_wp + n ) / specie2%mu ) /                        &
            GAMMA( ( specie2%nu + 1.0_wp ) / specie2%mu )  *                                       &
            GAMMA( ( specie2%nu + 1.0_wp ) / specie2%mu )**( specie2%b + n ) /                     &
            GAMMA( ( specie2%nu + 2.0_wp ) / specie2%mu )**( specie2%b + n )
       RETURN

     END FUNCTION collision_integral_delta_ba


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the collision integrals as given in SB2006, EQ.92
!--------------------------------------------------------------------------------------------------!
    FUNCTION collision_integral_theta_b(specie, n)

       INTEGER(iwp), INTENT(IN) :: n
       CLASS(cloud_coefficients), INTENT(IN) :: specie
       REAL(wp) ::  collision_integral_theta_b
!
!--    Eq. 92 from SB2006
       collision_integral_theta_b = &
                  GAMMA( ( 2.0_wp * specie%beta + 2.0_wp * specie%b + specie%nu + 1.0_wp + n ) /   &
                          specie%mu ) /                                                            &
                  GAMMA( ( 2.0_wp * specie%b + specie%nu + 1.0_wp + n) / specie%mu) *              &
                  GAMMA( ( specie%nu + 1.0_wp ) / specie%mu)**( 2.0_wp * specie%beta ) /           &
                  GAMMA( ( specie%nu + 2.0_wp ) / specie%mu)**( 2.0_wp * specie%beta )
       RETURN

    END FUNCTION collision_integral_theta_b


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the collision integrals as given in SB2006, EQ.93
!--------------------------------------------------------------------------------------------------!
    FUNCTION collision_integral_theta_ba(specie1, specie2, n)

       INTEGER(iwp), INTENT(IN) :: n
       CLASS(cloud_coefficients), INTENT(IN) :: specie1
       CLASS(cloud_coefficients), INTENT(IN) :: specie2
       REAL(wp) ::  collision_integral_theta_ba
!
!--    Eq. 93 from SB2006
       collision_integral_theta_ba = 2.0_wp *                                                      &
          GAMMA( ( specie1%beta + 2.0_wp * specie1%b + specie1%nu + 1.0_wp ) / specie1%mu ) /      &
          GAMMA( ( 2.0_wp * specie1%b + specie1%nu + 1.0_wp ) / specie1%mu ) *                     &
          GAMMA( ( specie1%nu + 1.0_wp ) / specie1%mu )**(specie1%beta) /                          &
          GAMMA( ( specie1%nu + 2.0_wp ) / specie1%mu )**(specie1%beta) *                          &
          GAMMA( ( specie2%beta + 2.0_wp * specie2%b + specie2%nu + 1.0_wp + n ) / specie2%mu ) /  &
          GAMMA( ( 2.0_wp * specie2%b + specie2%nu + 1.0_wp + n ) / specie2%mu ) *                 &
          GAMMA( ( specie2%nu + 1.0_wp ) / specie2%mu )**(specie2%beta) /                          &
          GAMMA( ( specie2%nu + 2.0_wp ) / specie2%mu )**(specie2%beta)
       RETURN

    END FUNCTION collision_integral_theta_ba


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the bulk ventilation coefficient as in Eq. 88 of SB2006
!--------------------------------------------------------------------------------------------------!
    FUNCTION ventilation_coeff_a( specie, n )

       REAL(wp) ::  ventilation_coeff_a  !<
       INTEGER, INTENT(IN) ::  n         !<
       CLASS(cloud_coefficients), INTENT(IN) :: specie  !<

       ventilation_coeff_a = specie%a_ven * GAMMA( ( specie%nu + n + specie%b ) / specie%mu ) /    &
                                            GAMMA((specie%nu + 1.0_wp ) / specie%mu )         *    &
                                          ( GAMMA((specie%nu + 1.0_wp ) / specie%mu )         /    &
                                            GAMMA((specie%nu + 2.0_wp ) / specie%mu )              &
                                          )**( specie%b + n - 1.0_wp)
       RETURN

    END FUNCTION ventilation_coeff_a


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the bulk ventilation coefficient as in Eq. 89 of SB2006
!--------------------------------------------------------------------------------------------------!
    FUNCTION ventilation_coeff_b( specie, n )

       REAL(wp) ::  ventilation_coeff_b !<
       INTEGER, INTENT(IN) ::  n        !<
       CLASS(cloud_coefficients), INTENT(IN) :: specie  !<

       REAL(wp), PARAMETER :: m_f = 0.500 ! see PK, S.541.

       ventilation_coeff_b = specie%b_ven *                                                        &
            GAMMA( ( specie%nu + n + ( m_f + 1.0_wp ) * specie%b + m_f * specie%beta ) /           &
                     specie%mu ) /                                                                 &
            GAMMA( ( specie%nu + 1.0_wp ) / specie%mu ) *                                          &
          ( GAMMA( ( specie%nu + 1.0_wp ) / specie%mu ) /                                          &
            GAMMA( ( specie%nu + 2.0_wp ) / specie%mu )                                            &
          )**( ( m_f + 1.0_wp ) * specie%b + m_f * specie%beta + n - 1.0_wp )

       RETURN

    END FUNCTION ventilation_coeff_b


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function gernerically computes the mean mass of specie and limit it to allowed masses from
!> scheme. A very small number of particles is added to ensure that it will not devided by zero.
!--------------------------------------------------------------------------------------------------!
    FUNCTION mean_mass( specie, mixing_ratio, number_concentration, rho_env )

       CLASS(cloud_coefficients), INTENT(IN) :: specie              !<
       REAL(wp), INTENT(IN) ::  mixing_ratio                        !<
       REAL(wp), INTENT(IN) ::  number_concentration                !<
       REAL(wp), INTENT(IN) ::  rho_env                             !<
       REAL(wp) ::  mean_mass                                       !<

       mean_mass = MIN( MAX( (rho_env * mixing_ratio / ( number_concentration + eps_sb_n ) ),      &
                            specie%x_min ), specie%x_max )
       RETURN

    END FUNCTION mean_mass


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the terminal fall velocity as in SB2006 Eq. 28
!--------------------------------------------------------------------------------------------------!
    FUNCTION terminal_fall_velocity( specie, particle_mean_mass, rho_env )

       CLASS(cloud_coefficients), INTENT(IN) :: specie              !<
       REAL(wp), INTENT(IN) ::  particle_mean_mass                  !<
       REAL(wp), INTENT(IN) ::  rho_env                             !<
       REAL(wp) ::  terminal_fall_velocity                          !<

       terminal_fall_velocity = specie%alpha * particle_mean_mass**specie%beta *                   &
                                ( rho_0 / rho_env )**0.5

       RETURN

    END FUNCTION terminal_fall_velocity


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function computes the mean diameter of the cloud specie as in SB2006 Eq. 32
!--------------------------------------------------------------------------------------------------!
    FUNCTION mean_diameter( specie, particle_mean_mass )

       CLASS(cloud_coefficients), INTENT(IN) ::  specie             !<
       REAL(wp), INTENT(IN) ::  particle_mean_mass                  !<
       REAL(wp) ::  mean_diameter                                   !<

       mean_diameter = specie%a * particle_mean_mass**specie%b
       RETURN

    END FUNCTION mean_diameter


 END MODULE bulk_cloud_model_mod
