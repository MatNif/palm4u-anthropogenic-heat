!> @file time_integration.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
! Copyright 2023 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Integration in time of the model equations, statistical analysis and graphic
!> output
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE time_integration


#if defined( __parallel )
    USE MPI
#endif

    USE advec_ws,                                                                                  &
        ONLY:  ws_statistics

    USE arrays_3d,                                                                                 &
        ONLY:  dzu,                                                                                &
               prho,                                                                               &
               pt,                                                                                 &
               pt_init,                                                                            &
               q,                                                                                  &
               q_init,                                                                             &
               ref_state,                                                                          &
               rho_ocean,                                                                          &
               tend,                                                                               &
               u,                                                                                  &
               v,                                                                                  &
               vpt

    USE biometeorology_mod,                                                                        &
        ONLY:  bio_calculate_thermal_index_maps,                                                   &
               thermal_comfort,                                                                    &
               bio_calculate_uv_exposure,                                                          &
               uv_exposure

    USE bulk_cloud_model_mod,                                                                      &
        ONLY: bulk_cloud_model,                                                                    &
              calc_liquid_water_content

    USE calc_mean_profile_mod,                                                                     &
        ONLY:  calc_mean_profile

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    USE chem_modules,                                                                              &
        ONLY:  bc_cs_t_val,                                                                        &
               chem_species

    USE control_parameters,                                                                        &
        ONLY:  advected_distance_x,                                                                &
               advected_distance_y,                                                                &
               air_chemistry,                                                                      &
               average_count_3d,                                                                   &
               averaging_interval,                                                                 &
               averaging_interval_pr,                                                              &
               bc_lr_cyc,                                                                          &
               bc_ns_cyc,                                                                          &
               bc_pt_t_val,                                                                        &
               bc_q_t_val,                                                                         &
               biometeorology,                                                                     &
               call_psolver_at_all_substeps,                                                       &
               child_domain,                                                                       &
               constant_flux_layer,                                                                &
               constant_heatflux,                                                                  &
               create_disturbances,                                                                &
               constant_diffusion,                                                                 &
               coupling_start_time,                                                                &
               current_timestep_number,                                                            &
               debug_output_timestep,                                                              &
               debug_string,                                                                       &
               dcep,                                                                               &
               disturbance_created,                                                                &
               disturbance_energy_limit,                                                           &
               dist_range,                                                                         &
               dopr_n,                                                                             &
               do_sum,                                                                             &
               dt_3d,                                                                              &
               dt_averaging_input,                                                                 &
               dt_averaging_input_pr,                                                              &
               dt_coupling,                                                                        &
               dt_data_output_av,                                                                  &
               dt_disturb,                                                                         &
               dt_do2d_xy,                                                                         &
               dt_do2d_xz,                                                                         &
               dt_do2d_yz,                                                                         &
               dt_do3d,                                                                            &
               dt_domask,                                                                          &
               dt_dopr,                                                                            &
               dt_dopr_listing,                                                                    &
               dt_dots,                                                                            &
               dt_run_control,                                                                     &
               end_time,                                                                           &
               first_call_mas,                                                                     &
               galilei_transformation,                                                             &
               humidity,                                                                           &
               indoor_model,                                                                       &
               initializing_actions,                                                               &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               land_surface,                                                                       &
               large_scale_forcing,                                                                &
               loop_optimization,                                                                  &
               lsf_surf,                                                                           &
               lsf_vert,                                                                           &
               masks,                                                                              &
               multi_agent_system_end,                                                             &
               multi_agent_system_start,                                                           &
               nesting_offline,                                                                    &
               neutral,                                                                            &
               nr_timesteps_this_run,                                                              &
               nudging,                                                                            &
               ocean_mode,                                                                         &
               pt_reference,                                                                       &
               pt_slope_offset,                                                                    &
               random_heatflux,                                                                    &
               salsa,                                                                              &
               simulated_time,                                                                     &
               simulated_time_chr,                                                                 &
               skip_time_do2d_xy,                                                                  &
               skip_time_do2d_xz,                                                                  &
               skip_time_do2d_yz,                                                                  &
               skip_time_do3d,                                                                     &
               skip_time_domask,                                                                   &
               skip_time_dopr,                                                                     &
               skip_time_data_output_av,                                                           &
               sloping_surface,                                                                    &
               stop_dt,                                                                            &
               surface_output,                                                                     &
               syn_turb_gen,                                                                       &
               terminate_coupled,                                                                  &
               terminate_run,                                                                      &
               timestep_scheme,                                                                    &
               time_coupling,                                                                      &
               time_do2d_xy,                                                                       &
               time_do2d_xz,                                                                       &
               time_do2d_yz,                                                                       &
               time_do3d,                                                                          &
               time_domask,                                                                        &
               time_dopr,                                                                          &
               time_dopr_av,                                                                       &
               time_dopr_listing,                                                                  &
               time_dosp,                                                                          &
               time_dosp_av,                                                                       &
               time_dots,                                                                          &
               time_do_av,                                                                         &
               time_do_sla,                                                                        &
               time_disturb,                                                                       &
               time_run_control,                                                                   &
               time_since_reference_point,                                                         &
               timestep_count,                                                                     &
               turbulent_inflow,                                                                   &
               turbulent_outflow,                                                                  &
               urban_surface,                                                                      &
               use_initial_profile_as_reference,                                                   &
               use_single_reference_value,                                                         &
               u_gtrans,                                                                           &
               v_gtrans,                                                                           &
               virtual_flight,                                                                     &
               virtual_measurement,                                                                &
               ws_scheme_mom,                                                                      &
               ws_scheme_sca

#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  nested_run
#endif

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE data_output_particle_mod,                                                                  &
        ONLY:  dop_data_output_ptseries

    USE dcep_mod,                                                                                  &
        ONLY:  dcep_main

    USE diagnostic_output_quantities_mod,                                                          &
        ONLY:  doq_calculate,                                                                      &
               timestep_number_at_prev_calc

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE flight_mod,                                                                                &
        ONLY:  flight_measurement

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nzb,                                                                                &
               nzt

    USE indoor_model_mod,                                                                          &
        ONLY:  dt_indoor,                                                                          &
               im_main_heatcool,                                                                   &
               time_indoor

    USE interfaces

    USE kinds

    USE land_surface_model_mod,                                                                    &
        ONLY:  lsm_boundary_condition,                                                             &
               lsm_energy_balance,                                                                 &
               skip_time_do_lsm

    USE lsf_nudging_mod,                                                                           &
        ONLY:  calc_tnudge,                                                                        &
               ls_forcing_surf,                                                                    &
               ls_forcing_vert,                                                                    &
               nudge_ref

    USE module_interface,                                                                          &
        ONLY:  module_interface_actions,                                                           &
               module_interface_swap_timelevel,                                                    &
               module_interface_boundary_conditions,                                               &
               module_interface_exchange_horiz

    USE multi_agent_system_mod,                                                                    &
        ONLY:  agents_active,                                                                      &
               multi_agent_system

    USE nesting_offl_mod,                                                                          &
        ONLY:  nesting_offl_bc,                                                                    &
               nesting_offl_input,                                                                 &
               nesting_offl_interpolation_factor,                                                  &
               nesting_offl_mass_conservation

    USE ocean_mod,                                                                                 &
        ONLY:  prho_reference

    USE particle_attributes,                                                                       &
        ONLY:  dt_dopts,                                                                           &
               first_call_lpm,                                                                     &
               particle_advection,                                                                 &
               time_dopts

    USE pegrid

#if defined( __parallel )
    USE pmc_interface,                                                                             &
        ONLY:  atmosphere_ocean_coupled_run,                                                       &
               cpl_id,                                                                             &
               homogeneous_initialization_child,                                                   &
               nesting_mode,                                                                       &
               pmci_adjust_dt_coupling,                                                            &
               pmci_atmos_ocean,                                                                   &
               pmci_boundary_conds,                                                                &
               pmci_datatrans,                                                                     &
               pmci_synchronize,                                                                   &
               pmci_ensure_nest_mass_conservation,                                                 &
               pmci_set_swaplevel
#else
    USE pmc_interface,                                                                             &
        ONLY:  atmosphere_ocean_coupled_run,                                                       &
               cpl_id
#endif

    USE progress_bar,                                                                              &
        ONLY:  finish_progress_bar,                                                                &
               output_progress_bar

    USE prognostic_equations_mod,                                                                  &
        ONLY:  prognostic_equations_cache,                                                         &
               prognostic_equations_vector

    USE radiation_model_mod,                                                                       &
        ONLY:  dt_radiation,                                                                       &
               force_radiation_call,                                                               &
               radiation,                                                                          &
               radiation_control,                                                                  &
               radiation_interaction,                                                              &
               radiation_interactions,                                                             &
               skip_time_do_radiation,                                                             &
               time_radiation

    USE salsa_mod,                                                                                 &
        ONLY:  aerosol_number,                                                                     &
               aerosol_mass,                                                                       &
               bc_am_t_val,                                                                        &
               bc_an_t_val,                                                                        &
               bc_gt_t_val,                                                                        &
               nbins_aerosol,                                                                      &
               ncomponents_mass,                                                                   &
               ngases_salsa,                                                                       &
               salsa_gas,                                                                          &
               salsa_gases_from_chem,                                                              &
               skip_time_do_salsa

    USE spectra_mod,                                                                               &
        ONLY:  average_count_sp,                                                                   &
               averaging_interval_sp,                                                              &
               calc_spectra,                                                                       &
               dt_dosp,                                                                            &
               skip_time_dosp

    USE statistics,                                                                                &
        ONLY:  flow_statistics_called,                                                             &
               hom,                                                                                &
               pr_palm,                                                                            &
               sums_ls_l

    USE surface_layer_fluxes_mod,                                                                  &
        ONLY:  surface_layer_fluxes

    USE surface_data_output_mod,                                                                   &
        ONLY:  average_count_surf,                                                                 &
               averaging_interval_surf,                                                            &
               dt_dosurf,                                                                          &
               dt_dosurf_av,                                                                       &
               surface_data_output,                                                                &
               surface_data_output_averaging,                                                      &
               skip_time_dosurf,                                                                   &
               skip_time_dosurf_av,                                                                &
               time_dosurf,                                                                        &
               time_dosurf_av

    USE surface_mod,                                                                               &
        ONLY:  surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_usm

    USE synthetic_turbulence_generator_mod,                                                        &
        ONLY:  stg_main

    USE turbulence_closure_mod,                                                                    &
        ONLY:  tcm_diffusivities

    USE turbulent_inflow_mod,                                                                      &
        ONLY:  turbulent_inflow_impose

    USE urban_surface_mod,                                                                         &
        ONLY:  usm_boundary_condition,                                                             &
               usm_energy_balance

    USE virtual_measurement_mod,                                                                   &
        ONLY:  dt_virtual_measurement_pr,                                                          &
               dt_virtual_measurement_ts,                                                          &
               dt_virtual_measurement_tr,                                                          &
               time_virtual_measurement_pr,                                                        &
               time_virtual_measurement_ts,                                                        &
               time_virtual_measurement_tr,                                                        &
               vm_data_output,                                                                     &
               vm_sampling,                                                                        &
               vm_time_start


    USE wind_turbine_model_mod,                                                                    &
        ONLY:  dt_data_output_wtm,                                                                 &
               time_wtm,                                                                           &
               wind_turbine,                                                                       &
               wtm_data_output

#if defined( _OPENACC )
    USE arrays_3d,                                                                                 &
        ONLY:  d, dd2zu, ddzu, ddzw,                                                               &
               diss,                                                                               &
               diss_p,                                                                             &
               diss_l_u,                                                                           &
               diss_l_v,                                                                           &
               diss_l_w,                                                                           &
               diss_s_u,                                                                           &
               diss_s_v,                                                                           &
               diss_s_w,                                                                           &
               drho_air, drho_air_zw, dzw, e,                                                      &
               flux_l_u,                                                                           &
               flux_l_v,                                                                           &
               flux_l_w,                                                                           &
               flux_s_u,                                                                           &
               flux_s_v,                                                                           &
               flux_s_w,                                                                           &
               e_p,                                                                                &
               heatflux_output_conversion,                                                         &
               kh,                                                                                 &
               km,                                                                                 &
               momentumflux_output_conversion,                                                     &
               nc,                                                                                 &
               ni,                                                                                 &
               nr,                                                                                 &
               odf_x,                                                                              &
               odf_y,                                                                              &
               p,                                                                                  &
               pt_p,                                                                               &
               ptdf_x,                                                                             &
               ptdf_y,                                                                             &
               qc,                                                                                 &
               qi,                                                                                 &
               qr,                                                                                 &
               rdf,                                                                                &
               rdf_sc,                                                                             &
               rho_air,                                                                            &
               rho_air_zw,                                                                         &
               s, tdiss_m,                                                                         &
               te_m,                                                                               &
               tpt_m,                                                                              &
               tu_m,                                                                               &
               tv_m,                                                                               &
               tw_m,                                                                               &
               u_p,                                                                                &
               ug,                                                                                 &
               u_init,                                                                             &
               u_stokes_zu,                                                                        &
               v_p,                                                                                &
               vg,                                                                                 &
               v_init,                                                                             &
               v_stokes_zu,                                                                        &
               w,                                                                                  &
               w_p,                                                                                &
               zu

    USE control_parameters,                                                                        &
        ONLY:  tsc

    USE indices,                                                                                   &
        ONLY:  advc_flags_m,                                                                       &
               advc_flags_s,                                                                       &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nz,                                                                                 &
               nzb_max,                                                                            &
               topo_flags

    USE statistics,                                                                                &
        ONLY:  pr_max,                                                                             &
               rmask,                                                                              &
               statistic_regions,                                                                  &
               sums_l,                                                                             &
               sums_l_l,                                                                           &
               sums_us2_ws_l,                                                                      &
               sums_wsus_ws_l,                                                                     &
               sums_vs2_ws_l,                                                                      &
               sums_wsvs_ws_l,                                                                     &
               sums_ws2_ws_l,                                                                      &
               sums_wspts_ws_l,                                                                    &
               sums_wsqs_ws_l,                                                                     &
               sums_wssas_ws_l,                                                                    &
               sums_wsqcs_ws_l,                                                                    &
               sums_wsqrs_ws_l,                                                                    &
               sums_wsncs_ws_l,                                                                    &
               sums_wsnrs_ws_l,                                                                    &
               sums_wsss_ws_l,                                                                     &
               weight_substep,                                                                     &
               sums_salsa_ws_l,                                                                    &
               sums_wsqis_ws_l,                                                                    &
               sums_wsnis_ws_l

    USE surface_mod,                                                                               &
        ONLY:  bc_hv,                                                                              &
               enter_surface_arrays,                                                               &
               exit_surface_arrays
#endif


    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string   !<

    INTEGER(iwp) ::  ib                  !< index for aerosol size bins
    INTEGER(iwp) ::  ic                  !< index for aerosol mass bins
    INTEGER(iwp) ::  icc                 !< additional index for aerosol mass bins
    INTEGER(iwp) ::  ig                  !< index for salsa gases
    INTEGER(iwp) ::  mid                 !< masked output running index
    INTEGER(iwp) ::  n                   !< loop counter for chemistry species

    REAL(wp) ::  dt_3d_old                        !< temporary storage of timestep to be used for
                                                  !< steering of run control output interval
    REAL(wp) ::  time_since_reference_point_save  !< original value of
                                                  !< time_since_reference_point


!
!-- Copy data from arrays_3d
!$ACC DATA &
!$ACC COPY(d(nzb+1:nzt,nys:nyn,nxl:nxr)) &
!$ACC COPY(diss(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(e(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(u(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(v(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(w(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(kh(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(km(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(pt(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

!$ACC DATA &
!$ACC COPYIN(diss_l_u(0:nz+1,nys:nyn,0), flux_l_u(0:nz+1,nys:nyn,0)) &
!$ACC COPYIN(diss_l_v(0:nz+1,nys:nyn,0), flux_l_v(0:nz+1,nys:nyn,0)) &
!$ACC COPYIN(diss_l_w(0:nz+1,nys:nyn,0), flux_l_w(0:nz+1,nys:nyn,0)) &
!$ACC COPYIN(diss_s_u(0:nz+1,0), flux_s_u(0:nz+1,0)) &
!$ACC COPYIN(diss_s_v(0:nz+1,0), flux_s_v(0:nz+1,0)) &
!$ACC COPYIN(diss_s_w(0:nz+1,0), flux_s_w(0:nz+1,0)) &
!$ACC COPY(diss_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(e_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(u_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(v_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(w_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(pt_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tend(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tdiss_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(te_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tu_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tv_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tw_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(tpt_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

!$ACC DATA &
!$ACC COPYIN(rho_air(nzb:nzt+1), drho_air(nzb:nzt+1)) &
!$ACC COPYIN(rho_air_zw(nzb:nzt+1), drho_air_zw(nzb:nzt+1)) &
!$ACC COPYIN(zu(nzb:nzt+1)) &
!$ACC COPYIN(dzu(1:nzt+1), dzw(1:nzt+1)) &
!$ACC COPYIN(ddzu(1:nzt+1), dd2zu(1:nzt)) &
!$ACC COPYIN(ddzw(1:nzt+1)) &
!$ACC COPYIN(heatflux_output_conversion(nzb:nzt+1)) &
!$ACC COPYIN(momentumflux_output_conversion(nzb:nzt+1)) &
!$ACC COPYIN(rdf(nzb+1:nzt), rdf_sc(nzb+1:nzt)) &
!$ACC COPYIN(ptdf_x(nxlg:nxrg), ptdf_y(nysg:nyng)) &
!$ACC COPYIN(odf_x(nxlg:nxrg), odf_y(nysg:nyng)) &
!$ACC COPYIN(ref_state(0:nz+1)) &
!$ACC COPYIN(u_init(0:nz+1), v_init(0:nz+1)) &
!$ACC COPYIN(u_stokes_zu(nzb:nzt+1), v_stokes_zu(nzb:nzt+1)) &
!$ACC COPYIN(pt_init(0:nz+1)) &
!$ACC COPYIN(ug(0:nz+1), vg(0:nz+1))

!
!-- Copy data from control_parameters
!$ACC DATA &
!$ACC COPYIN(tsc(1:5))

!
!-- Copy data from grid_variables
!$ACC DATA &
!$ACC COPYIN(ddx, ddy)

!
!-- Copy data from indices
!$ACC DATA &
!$ACC COPYIN(advc_flags_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPYIN(advc_flags_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPYIN(topo_flags(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

!
!-- Copy data from surface_mod
!$ACC DATA &
!$ACC COPYIN(bc_hv) &
!$ACC COPYIN(bc_hv%i(1:bc_hv%ns)) &
!$ACC COPYIN(bc_hv%j(1:bc_hv%ns)) &
!$ACC COPYIN(bc_hv%k(1:bc_hv%ns)) &
!$ACC COPYIN(bc_hv%ioff(1:bc_hv%ns)) &
!$ACC COPYIN(bc_hv%joff(1:bc_hv%ns)) &
!$ACC COPYIN(bc_hv%koff(1:bc_hv%ns))

!
!-- Copy data from statistics
!$ACC DATA &
!$ACC COPYIN(hom(0:nz+1,1:2,1:4,0)) &
!$ACC COPYIN(rmask(nysg:nyng,nxlg:nxrg,0:statistic_regions)) &
!$ACC COPYIN(weight_substep(1:intermediate_timestep_count_max)) &
!$ACC COPY(sums_l(nzb:nzt+1,1:pr_max,0)) &
!$ACC COPY(sums_l_l(nzb:nzt+1,0:statistic_regions,0)) &
!$ACC COPY(sums_us2_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsus_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_vs2_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsvs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_ws2_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wspts_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wssas_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsqs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsqcs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsqis_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsqrs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsncs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsnis_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsnrs_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_wsss_ws_l(nzb:nzt+1,0)) &
!$ACC COPY(sums_salsa_ws_l(nzb:nzt+1,0))

!
!-- Next statement is to avoid compiler warnings about unused variables. Please
!-- remove in case that you are using them. ddx and ddy need to be defined in
!-- time_integration because of ACC COPYIN directive.
    ddx = ddx
    ddy = ddy

#if defined( _OPENACC )
    CALL enter_surface_arrays
#endif
!
!-- If child domains are initialized with domain-averaged parent profiles
!-- run pressure solver before time integration starts.
    IF ( TRIM( initializing_actions )  /= 'read_restart_data' )  THEN
#if defined( __parallel )
       IF ( nested_run )  THEN
          IF ( child_domain  .AND.  homogeneous_initialization_child )  THEN
             CALL pmci_ensure_nest_mass_conservation
             CALL pres
          ENDIF
       ENDIF
#endif
    ENDIF
!
!-- At beginning determine the first time step
    CALL timestep

#if defined( __parallel )
!
!-- Synchronize the timestep in case of nested run.
    IF ( nested_run )  THEN
!
!--    Synchronization by unifying the time step.
!--    Global minimum of all time-steps is used for all.
       CALL pmci_synchronize
    ENDIF
#endif
!
!-- Determine and print out the run control quantities before the first time
!-- step of this run. For the initial run, some statistics (e.g. divergence)
!-- need to be determined first --> CALL flow_statistics at the beginning of
!-- run_control
    CALL run_control
!
!-- Data exchange between coupled models in case that a call has been omitted
!-- at the end of the previous run of a job chain.
    IF ( atmosphere_ocean_coupled_run )  THEN
!
!--    Adjust the time interval after which atmosphere and ocean are coupled, to be at least
!--    as large as the maximum timestep of the atmosphere and the ocean model.
#if defined( __parallel )
       CALL pmci_adjust_dt_coupling
#endif

!
!--    In case of model termination initiated by the local model the coupler
!--    must not be called because this would again cause an MPI hang.
       DO WHILE ( time_coupling >= dt_coupling  .AND.  terminate_coupled == 0 )
#if defined( __parallel )
          CALL pmci_atmos_ocean
#endif
          time_coupling = time_coupling - dt_coupling
       ENDDO
       IF ( time_coupling == 0.0_wp  .AND.  time_since_reference_point < dt_coupling )  THEN
          time_coupling = time_since_reference_point
       ENDIF
    ENDIF
!
!-- Execute necessary module actions before time-integration starts
    CALL module_interface_actions( 'before_integration' )

!
!-- Start of the time loop
    CALL location_message( 'time-stepping', 'start' )
    DO  WHILE ( simulated_time < end_time  .AND.  .NOT. stop_dt  .AND. .NOT. terminate_run )

       CALL cpu_log( log_point_s(10), 'timesteps', 'start' )

       IF ( debug_output_timestep )  THEN
           WRITE( debug_string, * ) 'time_integration', simulated_time
           CALL debug_message( debug_string, 'start' )
       ENDIF

!
!--    Determine ug, vg and w_subs in dependence on data from external file
!--    LSF_DATA
       IF ( large_scale_forcing .AND. lsf_vert )  THEN
           CALL ls_forcing_vert ( simulated_time )
           sums_ls_l = 0.0_wp
       ENDIF

!
!--    Set pt_init and q_init to the current profiles taken from
!--    NUDGING_DATA
       IF ( nudging )  THEN
           CALL nudge_ref ( simulated_time )
!
!--        Store temperature gradient at the top boundary for possible Neumann
!--        boundary condition
           bc_pt_t_val = ( pt_init(nzt+1) - pt_init(nzt) ) / dzu(nzt+1)
           bc_q_t_val  = ( q_init(nzt+1) - q_init(nzt) ) / dzu(nzt+1)
           IF ( air_chemistry )  THEN
              DO  n = 1, nvar
                 bc_cs_t_val = (  chem_species(n)%conc_pr_init(nzt+1)                              &
                                - chem_species(n)%conc_pr_init(nzt) )                              &
                               / dzu(nzt+1)
              ENDDO
           ENDIF
           IF ( salsa  .AND.  time_since_reference_point >= skip_time_do_salsa )  THEN
              DO  ib = 1, nbins_aerosol
                 bc_an_t_val = ( aerosol_number(ib)%init(nzt+1) - aerosol_number(ib)%init(nzt) ) / &
                               dzu(nzt+1)
                 DO  ic = 1, ncomponents_mass
                    icc = ( ic - 1 ) * nbins_aerosol + ib
                    bc_am_t_val = ( aerosol_mass(icc)%init(nzt+1) - aerosol_mass(icc)%init(nzt) ) /&
                                  dzu(nzt+1)
                 ENDDO
              ENDDO
              IF ( .NOT. salsa_gases_from_chem )  THEN
                 DO  ig = 1, ngases_salsa
                    bc_gt_t_val = ( salsa_gas(ig)%init(nzt+1) - salsa_gas(ig)%init(nzt) ) /        &
                                  dzu(nzt+1)
                 ENDDO
              ENDIF
           ENDIF
       ENDIF
!
!--    Input of boundary data.
       IF ( nesting_offline )  CALL nesting_offl_input
!
!--    Execute all other module actions routines
       CALL module_interface_actions( 'before_timestep' )

!
!--    Start of intermediate step loop
       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1
!
!--       Set the steering factors for the prognostic equations which depend
!--       on the timestep scheme
          CALL timestep_scheme_steering

!
!--       Calculate those variables needed in the tendency terms which need
!--       global communication
          IF ( .NOT. use_single_reference_value  .AND.  .NOT. use_initial_profile_as_reference )   &
          THEN
!
!--          Horizontally averaged profiles to be used as reference state in
!--          buoyancy terms (WARNING: only the respective last call of
!--          calc_mean_profile defines the reference state!)
             IF ( .NOT. neutral )  THEN
                CALL calc_mean_profile( pt, 4 )
                ref_state(:)  = hom(:,1,4,0) ! this is used in the buoyancy term
             ENDIF
             IF ( ocean_mode )  THEN
                CALL calc_mean_profile( rho_ocean, 64 )
                ref_state(:)  = hom(:,1,64,0)
             ENDIF
             IF ( humidity )  THEN
                CALL calc_mean_profile( vpt, 44 )
                ref_state(:)  = hom(:,1,44,0)
             ENDIF
!
!--          Assure that ref_state does not become zero at any level
!--          ( might be the case if a vertical level is completely occupied
!--            with topography ).
             ref_state = MERGE( MAXVAL(ref_state), ref_state, ref_state == 0.0_wp )
          ENDIF

          IF ( ( ws_scheme_mom .OR. ws_scheme_sca )  .AND.  intermediate_timestep_count == 1 )     &
          THEN
             CALL ws_statistics
          ENDIF
!
!--       In case of nudging calculate current nudging time scale and horizontal
!--       means of u, v, pt and q
          IF ( nudging )  THEN
             CALL calc_tnudge( simulated_time )
             CALL calc_mean_profile( u, 1 )
             CALL calc_mean_profile( v, 2 )
             CALL calc_mean_profile( pt, 4 )
             CALL calc_mean_profile( q, 41 )
          ENDIF
!
!--       Execute all other module actions routines
          CALL module_interface_actions( 'before_prognostic_equations' )
!
!--       Solve the prognostic equations. A fast cache optimized version with
!--       only one single loop is used in case of Piascek-Williams advection
!--       scheme. NEC vector machines use a different version, because
!--       in the other versions a good vectorization is prohibited due to
!--       inlining problems.
          IF ( loop_optimization == 'cache' )  THEN
             CALL prognostic_equations_cache
          ELSEIF ( loop_optimization == 'vector' )  THEN
             CALL prognostic_equations_vector
          ENDIF
!
!--       Movement of agents in multi agent system
          IF ( agents_active  .AND.  time_since_reference_point >= multi_agent_system_start  .AND. &
               time_since_reference_point <= multi_agent_system_end  .AND.                         &
               intermediate_timestep_count == 1 )                                                  &
          THEN
             CALL multi_agent_system
             first_call_mas = .FALSE.
          ENDIF

!
!--       Exchange of ghost points (lateral boundary conditions)
          CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'start' )

          CALL module_interface_exchange_horiz( 'after_prognostic_equation' )

          CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'stop' )

!
!--       Boundary conditions for the prognostic quantities (except of the
!--       velocities at the outflow in case of a non-cyclic lateral wall) and
!--       boundary conditions for module-specific variables
          CALL module_interface_boundary_conditions

!
!--       Incrementing timestep counter
          timestep_count = timestep_count + 1

          CALL cpu_log( log_point(28), 'swap_timelevel', 'start' )
!
!--       Set the swap level for all modules
          CALL module_interface_swap_timelevel( MOD( timestep_count, 2) )

#if defined( __parallel )
!
!--       Set the swap level for steering the pmc data transfer
          IF ( nested_run  .OR.  atmosphere_ocean_coupled_run )  THEN
             CALL pmci_set_swaplevel( MOD( timestep_count, 2) + 1 )  !> @todo: why the +1 ?
          ENDIF
#endif

          CALL cpu_log( log_point(28), 'swap_timelevel', 'stop' )

!
!--       Impose a turbulent inflow conditions. Depending on the chosen method, either the
!--       turbulence signal will be recycled or a turbulent inflow from a pre-run will be imposed.
!--       Inflow conditions should be set before nesting data exchange, because in case of pure
!--       vertical nesting with non-cyclic boundary conditions, parents should receive the inflow
!--       conditions from the childs. Otherwise strong horizontal flow divergence/convergence may
!--       appear at the parent inflow. Note, boundary conditions are imposed after the last
!--       time-integration substep and will apply in the next time-integration step.
          IF ( turbulent_inflow )  CALL turbulent_inflow_impose

#if defined( __parallel )
          IF ( nested_run )  THEN

             CALL cpu_log( log_point(60), 'nesting', 'start' )
!
!--          Domain nesting. The data transfer subroutines pmci_parent_datatrans
!--          and pmci_child_datatrans are called inside the wrapper
!--          subroutine pmci_datatrans according to the control parameters
!--          nesting_mode and nesting_datatransfer_mode.
!--          TO_DO: why is nesting_mode given as a parameter here?
             CALL pmci_datatrans( nesting_mode )

             IF ( TRIM( nesting_mode ) == 'two-way' )  THEN

                CALL cpu_log( log_point_s(92), 'exchange-horiz-nest', 'start' )
!
!--             Exchange_horiz is needed for all parent-domains after the
!--             anterpolation
                CALL module_interface_exchange_horiz( 'after_anterpolation' )

                CALL cpu_log( log_point_s(92), 'exchange-horiz-nest', 'stop' )

             ENDIF

!
!--          Set boundary conditions again after interpolation and anterpolation.
             CALL pmci_boundary_conds

             CALL cpu_log( log_point(60), 'nesting', 'stop' )

          ENDIF
#endif

!
!--       Temperature offset must be imposed at cyclic boundaries in x-direction
!--       when a sloping surface is used
          IF ( sloping_surface )  THEN
             IF ( nxl ==  0 )  pt(:,:,nxlg:nxl-1) = pt(:,:,nxlg:nxl-1) - pt_slope_offset
             IF ( nxr == nx )  pt(:,:,nxr+1:nxrg) = pt(:,:,nxr+1:nxrg) + pt_slope_offset
          ENDIF
!
!--       Set values at outflow boundary using the special outflow condition
          IF ( turbulent_outflow )  CALL outflow_turbulence

!
!--       Impose a random perturbation on the horizontal velocity field
          IF ( create_disturbances  .AND.  ( call_psolver_at_all_substeps  .AND.                   &
               intermediate_timestep_count == intermediate_timestep_count_max )                    &
               .OR. ( .NOT. call_psolver_at_all_substeps  .AND.  intermediate_timestep_count == 1 ) ) &
          THEN
             time_disturb = time_disturb + dt_3d
             IF ( time_disturb >= dt_disturb )  THEN
                IF ( disturbance_energy_limit /= 0.0_wp  .AND.                                     &
                     hom(nzb+5,1,pr_palm,0) < disturbance_energy_limit )  THEN
                   CALL disturb_field( 'u', tend, u )
                   CALL disturb_field( 'v', tend, v )
                ELSEIF ( ( .NOT. bc_lr_cyc  .OR.  .NOT. bc_ns_cyc )                                &
                         .AND. .NOT. child_domain  .AND.  .NOT.  nesting_offline )                 &
                THEN
!
!--                Runs with a non-cyclic lateral wall need perturbations
!--                near the inflow throughout the whole simulation
                   dist_range = 1
                   CALL disturb_field( 'u', tend, u )
                   CALL disturb_field( 'v', tend, v )
                   dist_range = 0
                ENDIF
                time_disturb = time_disturb - dt_disturb
             ENDIF
          ENDIF

!
!--       Map mesoscale model onto domain boundaries.
          IF ( nesting_offline  .AND.  intermediate_timestep_count ==                              &
                                       intermediate_timestep_count_max  )  THEN
!--          Determine interpolation factor before boundary conditions are updated.
             CALL nesting_offl_interpolation_factor
             CALL nesting_offl_bc
          ENDIF
!
!--       Impose a turbulent inflow using synthetic generated turbulence.
          IF ( syn_turb_gen  .AND.                                                                 &
               intermediate_timestep_count == intermediate_timestep_count_max )  THEN
             CALL cpu_log( log_point(57), 'synthetic_turbulence_gen', 'start' )
             CALL stg_main
             CALL cpu_log( log_point(57), 'synthetic_turbulence_gen', 'stop' )
          ENDIF
!
!--       Ensure mass conservation. This need to be done after imposing
!--       synthetic turbulence and top boundary condition for pressure is set to
!--       Neumann conditions.
!--       Is this also required in case of Dirichlet?
          IF ( nesting_offline )  CALL nesting_offl_mass_conservation
!
!--       Reduce the velocity divergence via the equation for perturbation
!--       pressure.
          IF ( intermediate_timestep_count == 1  .OR.  call_psolver_at_all_substeps )  THEN

#if defined( __parallel )
!
!--          Mass (volume) flux correction to ensure global mass conservation for child domains.
             IF ( child_domain .AND. .NOT. atmosphere_ocean_coupled_run )  THEN
                CALL pmci_ensure_nest_mass_conservation
             ENDIF
#endif
             CALL pres

          ENDIF
!
!--       Particle transport/physics with the Lagrangian particle model
!--       (only once during intermediate steps, because it uses an Euler-step)
!--       ### particle model should be moved before prognostic_equations, in order
!--       to regard droplet interactions directly

          CALL module_interface_actions( 'after_pressure_solver' )
!
!--       Interaction of droplets with temperature and mixing ratio.
!--       Droplet condensation and evaporation is calculated within
!--       advec_particles.
!
!--       If required, compute liquid water content
          IF ( bulk_cloud_model )  THEN
             CALL calc_liquid_water_content
          ENDIF
!
!--       If required, compute virtual potential temperature
          IF ( humidity )  THEN
             CALL compute_vpt
          ENDIF

!
!--       Compute the diffusion quantities
          IF ( .NOT. constant_diffusion )  THEN

!
!--          Determine surface fluxes shf and qsws and surface values
!--          pt_surface and q_surface in dependence on data from external
!--          file LSF_DATA respectively
             IF ( ( large_scale_forcing .AND. lsf_surf ) .AND.                                     &
                 intermediate_timestep_count == intermediate_timestep_count_max )                  &
             THEN
                CALL ls_forcing_surf( simulated_time )
             ENDIF

!
!--          First the vertical (and horizontal) fluxes in the surface
!--          (constant flux) layer are computed
             IF ( constant_flux_layer )  THEN
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'start' )
                CALL surface_layer_fluxes
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'stop' )
             ENDIF
!
!--          If required, solve the energy balance for the surface and run soil
!--          model. Call for horizontal as well as vertical surfaces
             IF ( land_surface .AND. time_since_reference_point >= skip_time_do_lsm)  THEN

                CALL cpu_log( log_point(54), 'land_surface', 'start' )
                CALL lsm_energy_balance( .FALSE. )

!
!--             At the end, set boundary conditons for potential temperature
!--             and humidity after running the land-surface model. This
!--             might be important for the nesting, where arrays are transfered.
                CALL lsm_boundary_condition


                CALL cpu_log( log_point(54), 'land_surface', 'stop' )
             ENDIF
!
!--          If required, solve the energy balance for urban surfaces and run
!--          the material heat model
             IF (urban_surface) THEN
                CALL cpu_log( log_point(74), 'urban_surface', 'start' )
                CALL usm_energy_balance( .FALSE. )

!
!--             At the end, set boundary conditons for potential temperature
!--             and humidity after running the urban-surface model. This
!--             might be important for the nesting, where arrays are transfered.
                CALL usm_boundary_condition

                CALL cpu_log( log_point(74), 'urban_surface', 'stop' )
             ENDIF
!
!--          Compute the diffusion coefficients
             CALL cpu_log( log_point(17), 'diffusivities', 'start' )
             IF ( .NOT. humidity ) THEN
                IF ( ocean_mode )  THEN
                   CALL tcm_diffusivities( prho, prho_reference )
                ELSE
                   CALL tcm_diffusivities( pt, pt_reference )
                ENDIF
             ELSE
                CALL tcm_diffusivities( vpt, pt_reference )
             ENDIF
             CALL cpu_log( log_point(17), 'diffusivities', 'stop' )

          ENDIF

       ENDDO   ! Intermediate step loop

!
!--    Will be used at some point by flow_statistics.
       !$ACC UPDATE &
       !$ACC HOST(sums_l_l(nzb:nzt+1,0:statistic_regions,0)) &
       !$ACC HOST(sums_us2_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsus_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_vs2_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsvs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_ws2_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wspts_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wssas_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsqs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsqcs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsqis_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsqrs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsncs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsnis_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsnrs_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_wsss_ws_l(nzb:nzt+1,0)) &
       !$ACC HOST(sums_salsa_ws_l(nzb:nzt+1,0))

!
!--    If required, calculate radiative fluxes and heating rates
       IF ( radiation  .AND.  time_since_reference_point >= skip_time_do_radiation )  THEN

          time_radiation = time_radiation + dt_3d
          IF ( time_radiation >= dt_radiation  .OR.  force_radiation_call )  THEN

             CALL cpu_log( log_point(50), 'radiation', 'start' )

             IF ( time_radiation >= dt_radiation )  THEN
                time_radiation = time_radiation - dt_radiation
             ENDIF

!
!--          Adjust the current time to the time step of the radiation model.
!--          Needed since radiation is pre-calculated and stored only on apparent
!--          solar positions
             time_since_reference_point_save = time_since_reference_point
             time_since_reference_point = time_since_reference_point -                             &
                                          MODULO(time_since_reference_point, dt_radiation)

             CALL radiation_control

             IF ( ( urban_surface  .OR.  land_surface )  .AND.  radiation_interactions )  THEN
                CALL cpu_log( log_point_s(46), 'radiation_interaction', 'start' )
                CALL radiation_interaction
                CALL cpu_log( log_point_s(46), 'radiation_interaction', 'stop' )
             ENDIF

!
!--          Return the current time to its original value
             time_since_reference_point = time_since_reference_point_save
!
!--          Reset forcing of radiation call
             force_radiation_call = .FALSE.

             CALL cpu_log( log_point(50), 'radiation', 'stop' )

          ENDIF
       ENDIF
!
!--    If DCEP set, call dcep main routine.
       IF ( dcep )  CALL dcep_main

       CALL module_interface_actions( 'update_emission_sources' )

!
!--    If required, calculate indoor temperature, waste heat, heat flux
!--    through wall, etc.
!--    dt_indoor steers the frequency of the indoor model calculations.
!--    Note, at first timestep indoor model is called, in order to provide
!--    a waste heat flux.
       IF ( indoor_model )  THEN

          time_indoor = time_indoor + dt_3d

          IF ( time_indoor >= dt_indoor  .OR.  current_timestep_number == 0 )  THEN

             IF ( time_indoor >= dt_indoor )  time_indoor = time_indoor - dt_indoor

             CALL cpu_log( log_point(76), 'indoor_model', 'start' )
             CALL im_main_heatcool
             CALL cpu_log( log_point(76), 'indoor_model', 'stop' )

          ENDIF
       ENDIF
!
!--    Increase simulation time and output times
       nr_timesteps_this_run      = nr_timesteps_this_run + 1
       current_timestep_number    = current_timestep_number + 1
       simulated_time             = simulated_time   + dt_3d
       time_since_reference_point = simulated_time - coupling_start_time
       simulated_time_chr         = time_to_string( time_since_reference_point )

       IF ( time_since_reference_point >= skip_time_data_output_av )  THEN
          time_do_av         = time_do_av       + dt_3d
       ENDIF

       IF ( time_since_reference_point >= skip_time_do2d_xy )  THEN
          time_do2d_xy       = time_do2d_xy     + dt_3d
       ENDIF

       IF ( time_since_reference_point >= skip_time_do2d_xz )  THEN
          time_do2d_xz       = time_do2d_xz     + dt_3d
       ENDIF

       IF ( time_since_reference_point >= skip_time_do2d_yz )  THEN
          time_do2d_yz       = time_do2d_yz     + dt_3d
       ENDIF

       IF ( time_since_reference_point >= skip_time_do3d    )  THEN
          time_do3d          = time_do3d        + dt_3d
       ENDIF

       DO  mid = 1, masks
          IF ( time_since_reference_point >= skip_time_domask(mid) )  THEN
             time_domask(mid) = time_domask(mid) + dt_3d
          ENDIF
       ENDDO

       IF ( time_since_reference_point >= skip_time_dosp )  THEN
          time_dosp       = time_dosp        + dt_3d
       ENDIF
       time_dots          = time_dots        + dt_3d

       IF ( time_since_reference_point >= skip_time_dopr )  THEN
          time_dopr       = time_dopr        + dt_3d
       ENDIF

       time_dopr_listing  = time_dopr_listing + dt_3d
       time_run_control   = time_run_control + dt_3d
!
!--    Increment time-counter for surface output
       IF ( surface_output )  THEN
          IF ( time_since_reference_point >= skip_time_dosurf )  THEN
             time_dosurf    = time_dosurf + dt_3d
          ENDIF

          IF ( time_since_reference_point >= skip_time_dosurf_av )  THEN
             time_dosurf_av = time_dosurf_av + dt_3d
          ENDIF
       ENDIF
!
!--    Increment time-counter for virtual measurements
       IF ( virtual_measurement  .AND.  vm_time_start <= time_since_reference_point )  THEN
          time_virtual_measurement_pr = time_virtual_measurement_pr + dt_3d
          time_virtual_measurement_ts = time_virtual_measurement_ts + dt_3d
          time_virtual_measurement_tr = time_virtual_measurement_tr + dt_3d
       ENDIF
!
!--    Increment time-counter for wind turbine data output
       IF ( wind_turbine )  time_wtm = time_wtm + dt_3d
!
!--    Increment time counter for particle time series
       IF ( particle_advection  .AND.  .NOT. first_call_lpm )  time_dopts = time_dopts + dt_3d

!
!--    Data exchange between coupled models
       IF ( atmosphere_ocean_coupled_run )  THEN
          time_coupling = time_coupling + dt_3d
!
!--       In case of model termination initiated by the atmosphere or the ocean model
!--       (terminate_coupled > 0), the coupler must be skipped because it would
!--       cause an MPI intercommunication hang.
          DO WHILE ( time_coupling >= dt_coupling  .AND.  terminate_coupled == 0 )
#if defined( __parallel )
             CALL pmci_atmos_ocean
!
!--          Adjust the time interval after which atmosphere and ocean are coupled, to be at least
!--          as large as the maximum timestep of the atmosphere and the ocean model.
             CALL pmci_adjust_dt_coupling
#endif
             time_coupling = time_coupling - dt_coupling
!
!--          Check needs to be done here because it requires global communication, so it can only
!--          be done after both atmosphere and ocean have exchanged data.
             CALL check_for_restart
          ENDDO
       ENDIF

!
!--    Biometeorology calculation of stationary thermal indices
!--    Todo (kanani): biometeorology needs own time_... treatment.
!--                   It might be that time_do2d_xy differs from time_do3d,
!--                   and then we might get trouble with the biomet output,
!--                   because we can have 2d and/or 3d biomet output!!
       IF ( biometeorology                                                                         &
            .AND. ( ( time_do3d >= dt_do3d  .AND.  time_since_reference_point >= skip_time_do3d )  &
                  .OR.                                                                             &
            ( time_do2d_xy >= dt_do2d_xy  .AND.  time_since_reference_point >= skip_time_do2d_xy ) &
                    ) )  THEN
!
!--       If required, do thermal comfort calculations
          IF ( thermal_comfort )  THEN
             CALL bio_calculate_thermal_index_maps ( .FALSE. )
          ENDIF
!
!--       If required, do UV exposure calculations
          IF ( uv_exposure )  THEN
             CALL bio_calculate_uv_exposure
          ENDIF
       ENDIF

!
!--    Execute alle other module actions routines
       CALL module_interface_actions( 'after_integration' )

!
!--    If Galilei transformation is used, determine the distance that the
!--    model has moved so far
       IF ( galilei_transformation )  THEN
          advected_distance_x = advected_distance_x + u_gtrans * dt_3d
          advected_distance_y = advected_distance_y + v_gtrans * dt_3d
       ENDIF

!
!--    Check, if restart is necessary (because cpu-time is expiring or
!--    because it is forced by user) and set stop flag
!--    This call is skipped if the remote model has already initiated a restart.
!--    In case of atmosphere - ocean coupling, check_for_restart is only called at coupling step
!--    (see further above).
       IF ( .NOT. terminate_run  .AND.  .NOT. atmosphere_ocean_coupled_run )  CALL check_for_restart

!
!--    Carry out statistical analysis and output at the requested output times.
!--    The MOD function is used for calculating the output time counters (like
!--    time_dopr) in order to regard a possible decrease of the output time
!--    interval in case of restart runs

!
!--    Set a flag indicating that so far no statistics have been created
!--    for this time step
       flow_statistics_called = .FALSE.

!
!--    If required, call flow_statistics for averaging in time
       IF ( averaging_interval_pr /= 0.0_wp  .AND.                                                 &
            ( dt_dopr - time_dopr ) <= averaging_interval_pr  .AND.                                &
            time_since_reference_point >= skip_time_dopr )  THEN
          time_dopr_av = time_dopr_av + dt_3d
          IF ( time_dopr_av >= dt_averaging_input_pr )  THEN
             do_sum = .TRUE.
             time_dopr_av = MOD( time_dopr_av, MAX( dt_averaging_input_pr, dt_3d ) )
          ENDIF
       ENDIF

       IF ( do_sum )  CALL flow_statistics

!
!--    Sum-up 3d-arrays for later output of time-averaged 2d/3d/masked data
       IF ( averaging_interval /= 0.0_wp  .AND.                                                    &
            ( dt_data_output_av - time_do_av ) <= averaging_interval  .AND.                        &
            time_since_reference_point >= skip_time_data_output_av )                               &
       THEN
          time_do_sla = time_do_sla + dt_3d
          IF ( time_do_sla >= dt_averaging_input )  THEN
             IF ( current_timestep_number > timestep_number_at_prev_calc )                         &
                CALL doq_calculate

             CALL sum_up_3d_data
             average_count_3d = average_count_3d + 1
             time_do_sla = MOD( time_do_sla, MAX( dt_averaging_input, dt_3d ) )
          ENDIF
       ENDIF
!
!--    Average surface data
       IF ( surface_output )  THEN
          IF ( averaging_interval_surf /= 0.0_wp                                                   &
                .AND.  ( dt_dosurf_av - time_dosurf_av ) <= averaging_interval_surf                &
                .AND.  time_since_reference_point >= skip_time_dosurf_av )  THEN
             IF ( time_dosurf_av >= dt_averaging_input )  THEN
                CALL surface_data_output_averaging
                average_count_surf = average_count_surf + 1
             ENDIF
          ENDIF
       ENDIF

!
!--    Calculate spectra for time averaging
       IF ( averaging_interval_sp /= 0.0_wp  .AND. ( dt_dosp - time_dosp ) <= averaging_interval_sp&
            .AND.  time_since_reference_point >= skip_time_dosp )  THEN
          time_dosp_av = time_dosp_av + dt_3d
          IF ( time_dosp_av >= dt_averaging_input_pr )  THEN
             CALL calc_spectra
             time_dosp_av = MOD( time_dosp_av, MAX( dt_averaging_input_pr, dt_3d ) )
          ENDIF
       ENDIF

!
!--    Call flight module and output data.
       IF ( virtual_flight )  THEN
          CALL flight_measurement
          CALL data_output_flight
       ENDIF
!
!--    Take virtual measurements. Depending on the type of measurement, different output intervals
!--    are considered.
       IF ( virtual_measurement  .AND.  vm_time_start <= time_since_reference_point )  THEN
!
!--       Timeseries data for point measurements.
          IF ( time_virtual_measurement_ts >= dt_virtual_measurement_ts )  THEN
             CALL vm_sampling( sample_pr = .FALSE., sample_ts = .TRUE., sample_tr = .FALSE. )
             CALL vm_data_output( output_pr = .FALSE., output_ts = .TRUE., output_tr = .FALSE. )
             time_virtual_measurement_ts = MOD( time_virtual_measurement_ts,                       &
                                                MAX( dt_virtual_measurement_ts, dt_3d ) )
          ENDIF
!
!--       Trajectory data.
          IF ( time_virtual_measurement_tr >= dt_virtual_measurement_tr )  THEN
             CALL vm_sampling( sample_pr = .FALSE., sample_ts = .FALSE., sample_tr = .TRUE. )
             CALL vm_data_output( output_pr = .FALSE., output_ts = .FALSE., output_tr = .TRUE. )
             time_virtual_measurement_tr = MOD( time_virtual_measurement_tr,                       &
                                                MAX( dt_virtual_measurement_tr, dt_3d ) )
          ENDIF
!
!--       Profile data.
          IF ( time_virtual_measurement_pr >= dt_virtual_measurement_pr )  THEN
             CALL vm_sampling( sample_pr = .TRUE., sample_ts = .FALSE., sample_tr = .FALSE. )
             CALL vm_data_output( output_pr = .TRUE., output_ts = .FALSE., output_tr = .FALSE. )
             time_virtual_measurement_pr = MOD( time_virtual_measurement_pr,                       &
                                                MAX( dt_virtual_measurement_pr, dt_3d ) )
          ENDIF
       ENDIF

!
!--    Output wind turbine data.
       IF ( wind_turbine  .AND.  time_wtm >= dt_data_output_wtm )  THEN
          CALL wtm_data_output
          time_wtm = MOD( time_wtm, MAX( dt_data_output_wtm, dt_3d ) )
       ENDIF

!
!--    Profile output (ASCII) on file
       IF ( time_dopr_listing >= dt_dopr_listing )  THEN
          CALL print_1d
          time_dopr_listing = MOD( time_dopr_listing, MAX( dt_dopr_listing, dt_3d ) )
       ENDIF

!
!--    Graphic output for PROFIL
       IF ( time_dopr >= dt_dopr  .AND.  time_since_reference_point >= skip_time_dopr )  THEN
          IF ( dopr_n /= 0 )  CALL data_output_profiles
          time_dopr = MOD( time_dopr, MAX( dt_dopr, dt_3d ) )
          time_dopr_av = 0.0_wp    ! due to averaging (see above)
       ENDIF

!
!--    Graphic output for time series
       IF ( time_dots >= dt_dots )  THEN
          CALL data_output_tseries
          time_dots = MOD( time_dots, MAX( dt_dots, dt_3d ) )
       ENDIF

!
!--    Output of spectra (formatted for use with PROFIL), in case of no
!--    time averaging, spectra has to be calculated before
       IF ( time_dosp >= dt_dosp  .AND.  time_since_reference_point >= skip_time_dosp )  THEN
          IF ( average_count_sp == 0 )  CALL calc_spectra
          CALL data_output_spectra
          time_dosp = MOD( time_dosp, MAX( dt_dosp, dt_3d ) )
       ENDIF

!
!--    2d-data output (cross-sections)
       IF ( time_do2d_xy >= dt_do2d_xy  .AND.  time_since_reference_point >= skip_time_do2d_xy )   &
       THEN
          IF ( current_timestep_number > timestep_number_at_prev_calc )  CALL doq_calculate
          CALL data_output_2d( 'xy', 0 )
          time_do2d_xy = MOD( time_do2d_xy, MAX( dt_do2d_xy, dt_3d ) )
       ENDIF
       IF ( time_do2d_xz >= dt_do2d_xz  .AND.  time_since_reference_point >= skip_time_do2d_xz )   &
       THEN
          IF ( current_timestep_number > timestep_number_at_prev_calc )  CALL doq_calculate
          CALL data_output_2d( 'xz', 0 )
          time_do2d_xz = MOD( time_do2d_xz, MAX( dt_do2d_xz, dt_3d ) )
       ENDIF
       IF ( time_do2d_yz >= dt_do2d_yz  .AND.  time_since_reference_point >= skip_time_do2d_yz )   &
       THEN
          IF ( current_timestep_number > timestep_number_at_prev_calc )  CALL doq_calculate
          CALL data_output_2d( 'yz', 0 )
          time_do2d_yz = MOD( time_do2d_yz, MAX( dt_do2d_yz, dt_3d ) )
       ENDIF

!
!--    3d-data output (volume data)
       IF ( time_do3d >= dt_do3d  .AND.  time_since_reference_point >= skip_time_do3d )  THEN
          IF ( current_timestep_number > timestep_number_at_prev_calc )                            &
             CALL doq_calculate

          CALL data_output_3d( 0 )
          time_do3d = MOD( time_do3d, MAX( dt_do3d, dt_3d ) )
       ENDIF

!
!--    Masked data output
       DO  mid = 1, masks
          IF ( time_domask(mid) >= dt_domask(mid)                                                  &
               .AND.  time_since_reference_point >= skip_time_domask(mid) )  THEN
             IF ( current_timestep_number > timestep_number_at_prev_calc )                         &
                CALL doq_calculate

             CALL data_output_mask( 0, mid )
             time_domask(mid) = MOD( time_domask(mid), MAX( dt_domask(mid), dt_3d ) )
          ENDIF
       ENDDO

!
!--    Output of time-averaged 2d/3d/masked data
       IF ( time_do_av >= dt_data_output_av                                                        &
            .AND.  time_since_reference_point >= skip_time_data_output_av )  THEN
          CALL average_3d_data
!
!--       Udate thermal comfort indices based on updated averaged input
          IF ( biometeorology  .AND.  thermal_comfort )  THEN
             CALL bio_calculate_thermal_index_maps ( .TRUE. )
          ENDIF
          CALL data_output_2d( 'xy', 1 )
          CALL data_output_2d( 'xz', 1 )
          CALL data_output_2d( 'yz', 1 )
          CALL data_output_3d( 1 )
          DO  mid = 1, masks
             CALL data_output_mask( 1, mid )
          ENDDO
          time_do_av = MOD( time_do_av, MAX( dt_data_output_av, dt_3d ) )
       ENDIF
!
!--    Output of surface data, instantaneous and averaged data
       IF ( surface_output )  THEN
          IF ( time_dosurf >= dt_dosurf  .AND.  time_since_reference_point >= skip_time_dosurf )   &
          THEN
             CALL surface_data_output( 0 )
             time_dosurf = MOD( time_dosurf, MAX( dt_dosurf, dt_3d ) )
          ENDIF
          IF ( time_dosurf_av >= dt_dosurf_av  .AND.                                               &
               time_since_reference_point >= skip_time_dosurf_av )  THEN
             CALL surface_data_output( 1 )
             time_dosurf_av = MOD( time_dosurf_av, MAX( dt_dosurf_av, dt_3d ) )
          ENDIF
       ENDIF

!
!--    Output of particle time series
       IF ( particle_advection )  THEN
          IF ( time_dopts >= dt_dopts  .AND.  .NOT. first_call_lpm )  THEN
             CALL dop_data_output_ptseries
             time_dopts = MOD( time_dopts, MAX( dt_dopts, dt_3d ) )
          ENDIF
       ENDIF

!
!--    If required, set the heat flux for the next time step to a random value
       IF ( constant_heatflux  .AND.  random_heatflux )  THEN
          IF ( surf_def%ns >= 1 )  THEN
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'start' )
             CALL disturb_heatflux( surf_def )
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'stop' )
          ENDIF
          IF ( surf_lsm%ns >= 1 )  THEN
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'start' )
             CALL disturb_heatflux( surf_lsm )
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'stop' )
          ENDIF
          IF ( surf_usm%ns >= 1 )  THEN
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'start' )
             CALL disturb_heatflux( surf_usm )
             CALL cpu_log( log_point(23), 'disturb_heatflux', 'stop' )
          ENDIF
       ENDIF

!
!--    Execute alle other module actions routines
       CALL module_interface_actions( 'after_timestep' )

!
!--    Determine size of next time step. Save timestep dt_3d because it is
!--    newly calculated in routine timestep, but required further below for
!--    steering the run control output interval
       dt_3d_old = dt_3d
       CALL timestep

#if defined( __parallel )
!
!--    Synchronize the timestep in case of nested run.
       IF ( nested_run )  THEN
!
!--       Synchronize by unifying the time step.
!--       Global minimum of all time-steps is used for all.
          CALL pmci_synchronize
       ENDIF
#endif

!
!--    Computation and output of run control parameters.
!--    This is also done whenever perturbations have been imposed
       IF ( time_run_control >= dt_run_control  .OR.                                               &
            timestep_scheme(1:5) /= 'runge'  .OR.  disturbance_created )                           &
       THEN
          CALL run_control
          IF ( time_run_control >= dt_run_control )  THEN
             time_run_control = MOD( time_run_control, MAX( dt_run_control, dt_3d_old ) )
          ENDIF
       ENDIF

!
!--    Output elapsed simulated time in form of a progress bar on stdout
       IF ( myid == 0  .AND.  cpl_id == 1 )  CALL output_progress_bar

       IF ( debug_output_timestep )  THEN
           WRITE( debug_string, * ) 'time_integration', simulated_time
           CALL debug_message( debug_string, 'end' )
       ENDIF

       CALL cpu_log( log_point_s(10), 'timesteps', 'stop' )

    ENDDO   ! time loop

    CALL module_interface_actions( 'after_time_integration' )

#if defined( _OPENACC )
    CALL exit_surface_arrays
#endif
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA
!$ACC END DATA

    IF ( myid == 0  .AND.  cpl_id == 1 )  CALL finish_progress_bar

    CALL location_message( 'time-stepping', 'finished' )

 END SUBROUTINE time_integration
