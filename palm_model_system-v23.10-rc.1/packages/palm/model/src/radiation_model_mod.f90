!> @file radiation_model_mod.f90
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
! Copyright 2015-2021 Institute of Computer Science of the Czech Academy of Sciences, Prague
! Copyright 2015-2021 Czech Technical University in Prague
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> Radiation models and interfaces:
!> Constant, simple and RRTMG models, interface to external radiation model
!> Radiative Transfer Model (RTM) version 3.0 for modelling of radiation
!> Interactions within urban canopy or other surface layer in complex terrain
!> Integrations of RTM with other PALM-4U modules:
!> Integration with RRTMG, USM, LSM, PCM, BIO modules
!>
!> @todo Move variable definitions used in radiation_init only to the subroutine as they are no
!>       longer required after initialization.
!> @todo Output of full column vertical profiles used in RRTMG
!> @todo Output of other rrtm arrays (such as volume mixing ratios)
!> @todo Optimize radiation_tendency routines
!>
!> @note Many variables have a leading dummy dimension (0:0) in order to match the assume-size shape
!>       expected by the RRTMG model.
!--------------------------------------------------------------------------------------------------!
 MODULE radiation_model_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  dzw,                                                                                &
               d_exner,                                                                            &
               exner,                                                                              &
               hyp,                                                                                &
               hyrho,                                                                              &
               nc,                                                                                 &
               pt,                                                                                 &
               p,                                                                                  &
               q,                                                                                  &
               qi,                                                                                 &
               ql,                                                                                 &
               u,                                                                                  &
               v,                                                                                  &
               w,                                                                                  &
               zu,                                                                                 &
               zw

    USE array_utilities,                                                                           &
        ONLY:  quicksort

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  barometric_formula,                                                                 &
               c_p,                                                                                &
               g,                                                                                  &
               lv_d_cp,                                                                            &
               l_v,                                                                                &
               pi,                                                                                 &
               r_d,                                                                                &
               rho_l,                                                                              &
               solar_constant,                                                                     &
               sigma_sb

    USE calc_mean_profile_mod,                                                                     &
        ONLY:  calc_mean_profile

    USE control_parameters,                                                                        &
        ONLY:  biometeorology,                                                                     &
               cloud_droplets,                                                                     &
               coupling_char,                                                                      &
               cyclic_fill_initialization,                                                         &
               dcep,                                                                               &
               debug_output,                                                                       &
               debug_output_timestep,                                                              &
               debug_string,                                                                       &
               dt_3d,                                                                              &
               dz,                                                                                 &
               dt_spinup,                                                                          &
               end_time,                                                                           &
               humidity,                                                                           &
               initializing_actions,                                                               &
               io_blocks,                                                                          &
               io_group,                                                                           &
               land_surface,                                                                       &
               large_scale_forcing,                                                                &
               latitude,                                                                           &
               length,                                                                             &
               longitude,                                                                          &
               loop_optimization,                                                                  &
               lsf_surf,                                                                           &
               message_string,                                                                     &
               plant_canopy,                                                                       &
               pt_surface,                                                                         &
               read_svf,                                                                           &
               restart_data_format_input,                                                          &
               restart_data_format_output,                                                         &
               restart_string,                                                                     &
               rho_surface,                                                                        &
               simulated_time,                                                                     &
               spinup_time,                                                                        &
               surface_pressure,                                                                   &
               time_since_reference_point,                                                         &
               urban_surface,                                                                      &
               varnamelength,                                                                      &
               write_svf

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nnx,                                                                                &
               nny,                                                                                &
               nx,                                                                                 &
               nxl,                                                                                &
               nxl_pe,                                                                             &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxr_pe,                                                                             &
               nxrg,                                                                               &
               nx_on_file,                                                                         &
               ny,                                                                                 &
               nyn,                                                                                &
               nyn_pe,                                                                             &
               nyng,                                                                               &
               nys,                                                                                &
               nys_pe,                                                                             &
               nysg,                                                                               &
               ny_on_file,                                                                         &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_top_ind,                                                                       &
               topo_flags

    USE, INTRINSIC :: iso_c_binding

    USE kinds

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model,                                                                   &
               microphysics_ice_phase,                                                             &
               microphysics_morrison,                                                              &
               na_init,                                                                            &
               nc_const,                                                                           &
               sigma_gc

#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  albedo_type_f,                                                                      &
               albedo_pars_f,                                                                      &
               building_type_f,                                                                    &
               building_surface_pars_f,                                                            &
               char_fill,                                                                          &
               char_lod,                                                                           &
               check_existence,                                                                    &
               close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               input_file_dynamic,                                                                 &
               input_pids_dynamic,                                                                 &
               num_var_pids,                                                                       &
               pavement_type_f,                                                                    &
               pids_id,                                                                            &
               open_read_file,                                                                     &
               real_1d_3d,                                                                         &
               vars_pids,                                                                          &
               vegetation_type_f,                                                                  &
               water_type_f

    USE palm_date_time_mod,                                                                        &
        ONLY:  date_time_str_len,                                                                  &
               get_date_time,                                                                      &
               hours_per_day,                                                                      &
               seconds_per_hour

    USE plant_canopy_model_mod,                                                                    &
        ONLY:  lad_s,                                                                              &
               pcm_calc_transpiration_rate,                                                        &
               pcm_latentflux,                                                                     &
               pcm_latentrate,                                                                     &
               pcm_sensibleflux,                                                                   &
               pcm_sensiblerate,                                                                   &
               pcm_transpiration_rate,                                                             &
               plant_canopy_transpiration

    USE pegrid

#if defined( __rrtmg )
    USE parrrsw,                                                                                   &
        ONLY:  naerec,                                                                             &
               nbndsw

    USE parrrtm,                                                                                   &
        ONLY:  nbndlw

    USE rrtmg_lw_init,                                                                             &
        ONLY:  rrtmg_lw_ini

    USE rrtmg_sw_init,                                                                             &
        ONLY:  rrtmg_sw_ini

    USE rrtmg_lw_rad,                                                                              &
        ONLY:  rrtmg_lw

    USE rrtmg_sw_rad,                                                                              &
        ONLY:  rrtmg_sw
#endif

#if defined( __tenstream )
    USE M_BUILDINGS,                                                                               &
        ONLY:  CHECK_BUILDINGS_CONSISTENCY,                                                        &
               CLONE_BUILDINGS,                                                                    &
               FACEIDX_BY_CELL_PLUS_OFFSET,                                                        &
               INIT_BUILDINGS,                                                                     &
               PPRTS_TOP_FACE,                                                                     &
               PPRTS_BOT_FACE,                                                                     &
               PPRTS_LEFT_FACE,                                                                    &
               PPRTS_RIGHT_FACE,                                                                   &
               PPRTS_REAR_FACE,                                                                    &
               PPRTS_FRONT_FACE,                                                                   &
               T_PPRTS_BUILDINGS

    USE M_PPRTS_BASE,                                                                              &
        ONLY:  T_SOLVER,                                                                           &
               ALLOCATE_PPRTS_SOLVER_FROM_COMMANDLINE

    USE M_DATA_PARAMETERS,                                                                         &
        ONLY:  DEFAULT_STR_LEN,                                                                    &
               IINTEGERS,                                                                          &
               INIT_MPI_DATA_PARAMETERS,                                                           &
               IREALS,                                                                             &
               MPIINT

    USE m_pprts_rrtmg,                                                                             &
        ONLY:  DESTROY_PPRTS_RRTMG,                                                                &
               PPRTS_RRTMG

    USE M_DYN_ATM_TO_RRTMG,                                                                        &
        ONLY:  DESTROY_TENSTR_ATM,                                                                 &
               SETUP_TENSTR_ATM,                                                                   &
               T_TENSTR_ATM

    USE M_HELPER_FUNCTIONS,                                                                        &
        ONLY:  CHKERR,                                                                             &
               REORDER_MPI_COMM

    USE M_TENSTR_RRTMG_SW_RAD,                                                                     &
        ONLY:  ts_earth_sun => EARTH_SUN

    USE M_TENSTR_PARRRTM,                                                                          &
        ONLY:  ts_nbndlw => NBNDLW

    USE M_TENSTR_PARRRSW,                                                                          &
        ONLY:  ts_nbndsw => NBNDSW

    USE M_TENSTR_RRLW_WVN,                                                                         &
        ONLY:  ts_ngblw => NGB,                                                                    &
               ts_ngptlw => NGPTLW

    USE M_TENSTR_RRSW_WVN,                                                                         &
        ONLY:  ts_ngbsw => NGB,                                                                    &
               ts_ngptsw => NGPTSW
#endif

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array,                                                              &
               rrd_mpi_io,                                                                         &
               rrd_mpi_io_global_array,                                                            &
               rd_mpi_io_open,                                                                     &
               rd_mpi_io_surface_filetypes,                                                        &
               rrd_mpi_io_surface,                                                                 &
               wrd_mpi_io,                                                                         &
               wrd_mpi_io_global_array,                                                            &
               wrd_mpi_io_surface,                                                                 &
               rd_mpi_io_close,                                                                    &
               tgh

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               ts_value

    USE surface_mod,                                                                               &
        ONLY:  albedop_dcep => albedop_urb,                                                        &
               emiss_dcep => emiss_urb,                                                            &
               fr_urb,                                                                             &
               ind_pav_green,                                                                      &
               ind_veg_wall,                                                                       &
               ind_wat_win,                                                                        &
               surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_out,                                                                           &
               surface_restore_elements,                                                           &
               surf_type,                                                                          &
               surf_usm,                                                                           &
               t_grad_dcep => t_grad_urb,                                                          &
               vertical_surfaces_exist

    IMPLICIT NONE

    CHARACTER(10) ::  radiation_scheme = 'clear-sky'  !< 'constant', 'clear-sky', 'rrtmg', or 'tenstream'

!
!-- Predefined Land surface classes (albedo_type) after Briegleb (1992)
    CHARACTER(37), DIMENSION(0:42), PARAMETER :: albedo_type_name = (/      &
                                   'user defined                         ', & !  0
                                   'ocean                                ', & !  1
                                   'mixed farming, tall grassland        ', & !  2
                                   'tall/medium grassland                ', & !  3
                                   'evergreen shrubland                  ', & !  4
                                   'short grassland/meadow/shrubland     ', & !  5
                                   'evergreen needleleaf forest          ', & !  6
                                   'mixed deciduous evergreen forest     ', & !  7
                                   'deciduous forest                     ', & !  8
                                   'tropical evergreen broadleaved forest', & !  9
                                   'medium/tall grassland/woodland       ', & ! 10
                                   'desert, sandy                        ', & ! 11
                                   'desert, rocky                        ', & ! 12
                                   'tundra                               ', & ! 13
                                   'land ice                             ', & ! 14
                                   'sea ice                              ', & ! 15
                                   'snow                                 ', & ! 16
                                   'bare soil                            ', & ! 17
                                   'asphalt/concrete mix                 ', & ! 18
                                   'asphalt (asphalt concrete)           ', & ! 19
                                   'concrete (Portland concrete)         ', & ! 20
                                   'sett                                 ', & ! 21
                                   'paving stones                        ', & ! 22
                                   'cobblestone                          ', & ! 23
                                   'metal                                ', & ! 24
                                   'wood                                 ', & ! 25
                                   'gravel                               ', & ! 26
                                   'fine gravel                          ', & ! 27
                                   'pebblestone                          ', & ! 28
                                   'woodchips                            ', & ! 29
                                   'tartan (sports)                      ', & ! 30
                                   'artifical turf (sports)              ', & ! 31
                                   'clay (sports)                        ', & ! 32
                                   'building (dummy)                     ', & ! 33
                                   'building wall - reflecting facade    ', & ! 34
                                   'building wall - bright facade        ', & ! 35
                                   'building wall - other materials      ', & ! 36
                                   'building window - double glazing     ', & ! 37
                                   'building window - double glazing     ', & ! 38
                                   'building window - reflecting         ', & ! 39
                                   'building roof - reflecting           ', & ! 40
                                   'building roof - bright               ', & ! 41
                                   'building roof - other materials      '  & ! 42
                                                         /)
!
!-- Indices of radiation-related input attributes in building_surface_pars
!-- (other are in urban_surface_mod)
    INTEGER(iwp), PARAMETER ::  ind_s_alb_b_wall  = 19  !< index for Broadband albedo of wall fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_l_wall  = 20  !< index for Longwave albedo of wall fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_s_wall  = 21  !< index for Shortwave albedo of wall fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_b_win   = 22  !< index for Broadband albedo of window fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_l_win   = 23  !< index for Longwave albedo of window fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_s_win   = 24  !< index for Shortwave albedo of window fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_b_green = 24  !< index for Broadband albedo of green fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_l_green = 25  !< index for Longwave albedo of green fraction
    INTEGER(iwp), PARAMETER ::  ind_s_alb_s_green = 26  !< index for Shortwave albedo of green fraction

    INTEGER(iwp) ::  albedo_type = 9999999  !< albedo surface type
    INTEGER(iwp) ::  day_of_year            !< day of the current year
    INTEGER(iwp) ::  dim_albedo_pars = 7    !< actual required dimension size of albedo_pars from static input (0..6)
    INTEGER(iwp) ::  dots_start_index_rtm   !< start index for time series of this module


    LOGICAL ::  average_radiation = .FALSE.            !< flag to set the calculation of radiation averaging for the domain
    LOGICAL ::  constant_albedo = .FALSE.              !< flag parameter indicating whether the albedo may
                                                       !< change depending on zenith
    LOGICAL ::  dcep_average_radiation = .FALSE.       !< flag to activiate average_radiation from DCEP module. It is define here
                                                       !< to avoid circular dependency.
    LOGICAL ::  force_radiation_call = .FALSE.         !< flag parameter for unscheduled radiation calls
    LOGICAL ::  lw_radiation = .TRUE.                  !< flag parameter indicating whether longwave radiation shall be calculated
    LOGICAL ::  radiation = .FALSE.                    !< flag parameter indicating whether the radiation model is used
    LOGICAL ::  radiation_interactions = .FALSE.       !< flag to activiate RTM (TRUE only if vertical
                                                       !< urban/land surface and trees exist)
    LOGICAL ::  radiation_interactions_on = .TRUE.     !< namelist flag to force RTM activiation regardless
                                                       !< to vertical urban/land surface and trees
    LOGICAL ::  radiation_only = .FALSE.               !< flag to activate radiation model (only RRTMG) without LSM or USM
                                                       !< does not work with option average_radiation
    LOGICAL ::  sun_direction = .FALSE.                !< flag parameter indicating whether solar direction shall be calculated
    LOGICAL ::  sun_up    = .TRUE.                     !< flag parameter indicating whether the sun is up or down
    LOGICAL ::  surface_reflections = .TRUE.           !< flag to switch the calculation of radiation
                                                       !< interaction between surfaces.
                                                       !< When it switched off, only the effect of buildings and trees shadow
                                                       !< will be considered. However fewer SVFs are expected.
    LOGICAL ::  sw_radiation = .TRUE.                  !< flag parameter indicating whether shortwave
                                                       !< radiation shall be calculated
    LOGICAL ::  unscheduled_radiation_calls = .FALSE.  !< flag parameter indicating whether additional calls
                                                       !< of the radiation code are allowed
#if defined( __rrtmg ) || defined( __tenstream )
    LOGICAL :: use_broadband_albedo = .FALSE.            !< namelist flag to use broadband albedo instead of diffuse/direct albedo
#endif

    REAL(wp), PARAMETER ::  emissivity_atm_clsky = 0.8_wp  !< emissivity of the clear-sky atmosphere

    REAL(wp) ::  albedo = 9999999.9_wp,           &  !< NAMELIST alpha
                 albedo_lw_dif = 9999999.9_wp,    &  !< NAMELIST aldif
                 albedo_lw_dir = 9999999.9_wp,    &  !< NAMELIST aldir
                 albedo_sw_dif = 9999999.9_wp,    &  !< NAMELIST asdif
                 albedo_sw_dir = 9999999.9_wp,    &  !< NAMELIST asdir
                 decl_1,                          &  !< declination coef. 1
                 decl_2,                          &  !< declination coef. 2
                 decl_3,                          &  !< declination coef. 3
                 dt_radiation = 0.0_wp,           &  !< radiation model timestep
                 emissivity = 9999999.9_wp,       &  !< NAMELIST surface emissivity
                 lon = 0.0_wp,                    &  !< longitude in radians
                 lat = 0.0_wp,                    &  !< latitude in radians
                 net_radiation = 0.0_wp,          &  !< net radiation at surface
                 skip_time_do_radiation = 0.0_wp, &  !< Radiation model is not called before this time
                 sky_trans,                       &  !< sky transmissivity
                 time_radiation = 0.0_wp,         &  !< time since last call of radiation code
                 trace_fluxes_above = 2000.0_wp,  &  !< NAMELIST option for debug printing of largest radiative fluxes
                                                     !< (W/m2 for surfaces, W/m3 for PC). -1=off, 0=all fluxes
                 min_stable_coszen = 0.0262_wp       !< 1.5 deg above horizon, eliminates most of circumsolar
    REAL(wp) ::  cos_zenith                          !< cosine of solar zenith angle, also z-coordinate of solar unit vector
    REAL(wp) ::  d_hours_day                         !< 1 / hours-per-day
    REAL(wp) ::  d_seconds_hour                      !< 1 / seconds-per-hour
    REAL(wp) ::  second_of_day                       !< second of the current day
    REAL(wp) ::  sun_dir_lat                         !< y-coordinate of solar unit vector
    REAL(wp) ::  sun_dir_lon                         !< x-coordinate of solar unit vector

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_net_av        !< average of net radiation (rad_net) at surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_in_xy_av   !< average of incoming longwave radiation at surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_out_xy_av  !< average of outgoing longwave radiation at surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_in_xy_av   !< average of incoming shortwave radiation at surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_out_xy_av  !< average of outgoing shortwave radiation at surface


!
!-- Land surface albedos for solar zenith angle of 60 degree after Briegleb (1992)
!-- (broadband, longwave, shortwave ):   bb,      lw,      sw,
    REAL(wp), DIMENSION(0:2,1:42), PARAMETER :: albedo_pars = RESHAPE( (/&
                                   0.06_wp, 0.06_wp, 0.06_wp,            & !  1 - ocean
                                   0.19_wp, 0.28_wp, 0.09_wp,            & !  2 - mixed farming, tall grassland
                                   0.23_wp, 0.33_wp, 0.11_wp,            & !  3 - tall/medium grassland
                                   0.23_wp, 0.33_wp, 0.11_wp,            & !  4 - evergreen shrubland
                                   0.25_wp, 0.34_wp, 0.14_wp,            & !  5 - short grassland/meadow/shrubland
                                   0.14_wp, 0.22_wp, 0.06_wp,            & !  6 - evergreen needleleaf forest
                                   0.17_wp, 0.27_wp, 0.06_wp,            & !  7 - mixed deciduous forest
                                   0.19_wp, 0.31_wp, 0.06_wp,            & !  8 - deciduous forest
                                   0.14_wp, 0.22_wp, 0.06_wp,            & !  9 - tropical evergreen broadleaved forest
                                   0.18_wp, 0.28_wp, 0.06_wp,            & ! 10 - medium/tall grassland/woodland
                                   0.43_wp, 0.51_wp, 0.35_wp,            & ! 11 - desert, sandy
                                   0.32_wp, 0.40_wp, 0.24_wp,            & ! 12 - desert, rocky
                                   0.19_wp, 0.27_wp, 0.10_wp,            & ! 13 - tundra
                                   0.77_wp, 0.65_wp, 0.90_wp,            & ! 14 - land ice
                                   0.77_wp, 0.65_wp, 0.90_wp,            & ! 15 - sea ice
                                   0.82_wp, 0.70_wp, 0.95_wp,            & ! 16 - snow
                                   0.08_wp, 0.08_wp, 0.08_wp,            & ! 17 - bare soil
                                   0.25_wp, 0.25_wp, 0.25_wp,            & ! 18 - asphalt/concrete mix / Same as concrete
                                   0.08_wp, 0.08_wp, 0.08_wp,            & ! 19 - asphalt (asphalt concrete) / Masson et al. (2002)
                                   0.35_wp, 0.35_wp, 0.35_wp,            & ! 20 - concrete (Portland concrete) / Yaghoobian et al. (2009)
                                   0.30_wp, 0.30_wp, 0.30_wp,            & ! 21 - sett / Own estimation
                                   0.25_wp, 0.25_wp, 0.25_wp,            & ! 22 - paving stone / Oke (1987)
                                   0.30_wp, 0.30_wp, 0.30_wp,            & ! 23 - cobblestone / Own estimation
                                   0.15_wp, 0.15_wp, 0.15_wp,            & ! 24 - metal / Oke (1987)
                                   0.20_wp, 0.20_wp, 0.20_wp,            & ! 25 - wood / Roberts et al. (2006)
                                   0.12_wp, 0.12_wp, 0.12_wp,            & ! 26 - gravel / Masson et al. (2002)
                                   0.12_wp, 0.12_wp, 0.12_wp,            & ! 27 - fine gravel / Same as gravel
                                   0.12_wp, 0.12_wp, 0.12_wp,            & ! 28 - pebblestone / Same as gravel
                                   0.20_wp, 0.20_wp, 0.20_wp,            & ! 29 - woodchips / Same as wood
                                   0.31_wp, 0.31_wp, 0.31_wp,            & ! 30 - tartan (sports) / Heldens (2010)
                                   0.08_wp, 0.08_wp, 0.08_wp,            & ! 31 - artificial turf (sports) / Yaghoobian et al. (2009)
                                   0.33_wp, 0.33_wp, 0.33_wp,            & ! 32 - clay (sports) / Oke (1987) + LBNL
                                   0.17_wp, 0.17_wp, 0.17_wp,            & ! 33 - building (dummy)
                                   0.60_wp, 0.60_wp, 0.60_wp,            & ! 34 - building wall type 1 - reflecting facade
                                   0.30_wp, 0.30_wp, 0.30_wp,            & ! 35 - building wall type 2 - bright facacde
                                   0.07_wp, 0.07_wp, 0.07_wp,            & ! 36 - building wall type 3) - other materials
                                   0.12_wp, 0.12_wp, 0.12_wp,            & ! 37 - building window type 1 - double glazing
                                   0.17_wp, 0.18_wp, 0.18_wp,            & ! 38 - building window type 2 - triple glazing
                                   0.48_wp, 0.48_wp, 0.48_wp,            & ! 39 - building window type 3 - reflecting
                                   0.60_wp, 0.60_wp, 0.60_wp,            & ! 40 - building roof type 1 - reflecting
                                   0.30_wp, 0.30_wp, 0.30_wp,            & ! 41 - building roof type 2 - bright
                                   0.07_wp, 0.07_wp, 0.07_wp             & ! 42 - building roof type 3 - other materials
                                 /), (/ 3, 42 /) )

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  &
                        rad_lw_cs_hr,                   &  !< longwave clear sky radiation heating rate (K/s)
                        rad_lw_cs_hr_av,                &  !< average of rad_lw_cs_hr
                        rad_lw_hr,                      &  !< longwave radiation heating rate (K/s)
                        rad_lw_hr_av,                   &  !< average of rad_sw_hr
                        rad_lw_in,                      &  !< incoming longwave radiation (W/m2)
                        rad_lw_in_av,                   &  !< average of rad_lw_in
                        rad_lw_out,                     &  !< outgoing longwave radiation (W/m2)
                        rad_lw_out_av,                  &  !< average of rad_lw_out
                        rad_sw_cs_hr,                   &  !< shortwave clear sky radiation heating rate (K/s)
                        rad_sw_cs_hr_av,                &  !< average of rad_sw_cs_hr
                        rad_sw_hr,                      &  !< shortwave radiation heating rate (K/s)
                        rad_sw_hr_av,                   &  !< average of rad_sw_hr
                        rad_sw_in,                      &  !< incoming shortwave radiation (W/m2)
                        rad_sw_in_av,                   &  !< average of rad_sw_in
                        rad_sw_out,                     &  !< outgoing shortwave radiation (W/m2)
                        rad_sw_out_av                      !< average of rad_sw_out


!
!-- Variables and parameters used in RRTMG only
#if defined( __rrtmg )
    CHARACTER(LEN=12) ::  rrtm_input_file = 'RAD_SND_DATA'  !< name of the NetCDF input file (sounding data)


!
!-- Flag parameters to be passed to RRTMG
    INTEGER(iwp), PARAMETER ::  rrtm_idrv     = 1, &  !< flag for longwave upward flux calculation option (0,1)
                                rrtm_inflglw  = 2, &  !< flag for lw cloud optical properties (0,1,2)
                                rrtm_iceflglw = 2, &  !< flag for lw ice particle specifications (0,1,2,3)
                                rrtm_liqflglw = 1, &  !< flag for lw liquid droplet specifications
                                rrtm_inflgsw  = 2, &  !< flag for sw cloud optical properties (0,1,2)
                                rrtm_iceflgsw = 2, &  !< flag for sw ice particle specifications (0,1,2,3)
                                rrtm_liqflgsw = 1     !< flag for sw liquid droplet specifications

!
!-- The following variables should only be changed with care, as this will require further setting
!-- of some variables, which is currently not implemented (aerosols, ice phase).
    INTEGER(iwp) ::  nzt_rad,           &  !< upper vertical limit for radiation calculations
                     rrtm_icld = 0,     &  !< cloud flag (0: clear sky column, 1: cloudy column)
                     rrtm_iaer = 0         !< aerosol option flag (0: no aerosol layers, for lw only: 6
                                           !< (requires setting of rrtm_sw_ecaer), 10: one or more aerosol layers (not implemented)

    INTEGER(iwp) ::  nc_stat  !< local variable for storin the result of netCDF calls for error message handling

    LOGICAL ::  snd_exists = .FALSE.  !< flag parameter to check whether a user-defined input files exists
    LOGICAL ::  sw_exists  = .FALSE.  !< flag parameter to check whether that required rrtmg sw file exists
    LOGICAL ::  lw_exists  = .FALSE.  !< flag parameter to check whether that required rrtmg lw file exists

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  hyp_snd,     &  !< hypostatic pressure from sounding data (hPa)
                                            rrtm_tsfc,   &  !< dummy array for storing surface temperature
                                            t_snd           !< actual temperature from sounding data (hPa)

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rrtm_ccl4vmr,   &  !< CCL4 volume mixing ratio (g/mol)
                                              rrtm_cfc11vmr,  &  !< CFC11 volume mixing ratio (g/mol)
                                              rrtm_cfc12vmr,  &  !< CFC12 volume mixing ratio (g/mol)
                                              rrtm_cfc22vmr,  &  !< CFC22 volume mixing ratio (g/mol)
                                              rrtm_ch4vmr,    &  !< CH4 volume mixing ratio
                                              rrtm_cicewp,    &  !< in-cloud ice water path (g/m2)
                                              rrtm_cldfr,     &  !< cloud fraction (0,1)
                                              rrtm_cliqwp,    &  !< in-cloud liquid water path (g/m2)
                                              rrtm_co2vmr,    &  !< CO2 volume mixing ratio (g/mol)
                                              rrtm_emis,      &  !< surface emissivity (0-1)
                                              rrtm_h2ovmr,    &  !< H2O volume mixing ratio
                                              rrtm_n2ovmr,    &  !< N2O volume mixing ratio
                                              rrtm_o2vmr,     &  !< O2 volume mixing ratio
                                              rrtm_o3vmr,     &  !< O3 volume mixing ratio
                                              rrtm_play,      &  !< pressure layers (hPa, zu-grid)
                                              rrtm_plev,      &  !< pressure layers (hPa, zw-grid)
                                              rrtm_reice,     &  !< cloud ice effective radius (microns)
                                              rrtm_reliq,     &  !< cloud water drop effective radius (microns)
                                              rrtm_tlay,      &  !< actual temperature (K, zu-grid)
                                              rrtm_tlev,      &  !< actual temperature (K, zw-grid)
                                              rrtm_lwdflx,    &  !< RRTM output of incoming longwave radiation flux (W/m2)
                                              rrtm_lwdflxc,   &  !< RRTM output of outgoing clear sky longwave radiation flux (W/m2)
                                              rrtm_lwuflx,    &  !< RRTM output of outgoing longwave radiation flux (W/m2)
                                              rrtm_lwuflxc,   &  !< RRTM output of incoming clear sky longwave radiation flux (W/m2)
                                              rrtm_lwuflx_dt, &  !< RRTM output of incoming clear sky longwave radiation flux (W/m2)
                                              rrtm_lwuflxc_dt,&  !< RRTM output of outgoing clear sky longwave radiation flux (W/m2)
                                              rrtm_lwhr,      &  !< RRTM output of longwave radiation heating rate (K/d)
                                              rrtm_lwhrc,     &  !< RRTM output of incoming longwave
                                                                 !< clear sky radiation heating rate (K/d)
                                              rrtm_swdflx,    &  !< RRTM output of incoming shortwave radiation flux (W/m2)
                                              rrtm_swdflxc,   &  !< RRTM output of outgoing clear sky
                                                                 !< shortwave radiation flux (W/m2)
                                              rrtm_swuflx,    &  !< RRTM output of outgoing shortwave radiation flux (W/m2)
                                              rrtm_swuflxc,   &  !< RRTM output of incoming clear sky
                                                                 !< shortwave radiation flux (W/m2)
                                              rrtm_swhr,      &  !< RRTM output of shortwave radiation heating rate (K/d)
                                              rrtm_swhrc,     &  !< RRTM output of incoming shortwave
                                                                 !< clear sky radiation heating rate (K/d)
                                              rrtm_dirdflux,  &  !< RRTM output of incoming direct shortwave (W/m2)
                                              rrtm_difdflux      !< RRTM output of incoming diffuse shortwave (W/m2)

    REAL(wp), DIMENSION(1) ::  rrtm_aldif,     &  !< surface albedo for longwave diffuse radiation
                               rrtm_aldir,     &  !< surface albedo for longwave direct radiation
                               rrtm_asdif,     &  !< surface albedo for shortwave diffuse radiation
                               rrtm_asdir         !< surface albedo for shortwave direct radiation

!
!-- Definition of arrays that are currently not used for calling RRTMG (due to setting of flag parameters)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  rad_lw_cs_in,   &  !< incoming clear sky longwave radiation (W/m2) (not used)
                                                rad_lw_cs_out,  &  !< outgoing clear sky longwave radiation (W/m2) (not used)
                                                rad_sw_cs_in,   &  !< incoming clear sky shortwave radiation (W/m2) (not used)
                                                rad_sw_cs_out,  &  !< outgoing clear sky shortwave radiation (W/m2) (not used)
                                                rrtm_lw_tauaer, &  !< lw aerosol optical depth
                                                rrtm_lw_taucld, &  !< lw in-cloud optical depth
                                                rrtm_sw_taucld, &  !< sw in-cloud optical depth
                                                rrtm_sw_ssacld, &  !< sw in-cloud single scattering albedo
                                                rrtm_sw_asmcld, &  !< sw in-cloud asymmetry parameter
                                                rrtm_sw_fsfcld, &  !< sw in-cloud forward scattering fraction
                                                rrtm_sw_tauaer, &  !< sw aerosol optical depth
                                                rrtm_sw_ssaaer, &  !< sw aerosol single scattering albedo
                                                rrtm_sw_asmaer, &  !< sw aerosol asymmetry parameter
                                                rrtm_sw_ecaer      !< sw aerosol optical detph at 0.55 microns (rrtm_iaer = 6 only)

#endif

#if defined( __rrtmg ) || defined( __tenstream )
    REAL(wp), PARAMETER ::  mol_mass_air_d_wv = 1.607793_wp  !< molecular weight dry air / water vapor
#endif

#if defined(__tenstream)
!
!-- TenStream variables
    CHARACTER(LEN=10), PARAMETER :: tenstream_solver = '2str'  !< kind of TenStream solver, default is a 1D solver can be changed via TenStream options
    CHARACTER(LEN=DEFAULT_STR_LEN), PARAMETER ::  ts_atm_filename = 'TS_BACKGROUND_ATM' !< Filename of background atmosphere file. ASCII file with columns:
                                                                                        !< z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    INTEGER(MPIINT) ::  ts_nranksx  !< number of ranks along x-axis
    INTEGER(MPIINT) ::  ts_nranksy  !< number of ranks along y-axis
    INTEGER(MPIINT) ::  ts_comm     !< number of ranks along x-axis

    INTEGER(IINTEGERS) ::  ts_icollapse = -1_IINTEGERS  !< flag to return flux results from the background atmosphere above the dynamical grid (1) or not (-1)
    INTEGER(IINTEGERS) ::  ts_xm                        !< tenstream x-direction size
    INTEGER(IINTEGERS) ::  ts_ym                        !< tenstream y-direction size
    INTEGER(IINTEGERS) ::  ts_zm                        !< tenstream z-direction size

    INTEGER(IINTEGERS), ALLOCATABLE ::  ts_nxproc(:)  !< number of ranks along x-axis
    INTEGER(IINTEGERS), ALLOCATABLE ::  ts_nyproc(:)  !< number of ranks along y-axis

    REAL(IREALS) ::  albedo_sol  !< solar albedo (global value), not used if a 2d albedo is provided
    REAL(IREALS) ::  albedo_th   !< thermal albedo (global value), not used if a 2d albedo is provided
    REAL(IREALS) ::  ts_dx       !< tenstream dx
    REAL(IREALS) ::  ts_dy       !< tenstream dy

    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ts_h2ovmr            !< H2O volume mixing ratio
    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ts_lwc               !< Liquid water cloud content [g/kg] and effective radius in micron
    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ts_plev              !< pressure on layer interfaces [hPa] in zu grid
    REAL(IREALS), DIMENSION(:)    , ALLOCATABLE, TARGET ::  ts_play              !< pressure on layer interfaces [hPa] in zw grid
    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ts_reliq             !< cloud water drop effective radius (microns)
    REAL(IREALS), DIMENSION(:,:)  , ALLOCATABLE, TARGET ::  ts_skin_temperature  !< skin temperature interfaces [K]
    REAL(IREALS), DIMENSION(:,:)  , ALLOCATABLE, TARGET ::  ts_solar_albedo_2d   !< solar albedo interfaces
    REAL(IREALS), DIMENSION(:,:)  , ALLOCATABLE, TARGET ::  ts_thermal_albedo_2d !< thermal albedo interfaces
    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ts_tlay              !< temperature on layer interfaces [K]
    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ts_tlev              !< temperature on level interfaces [K]
!
!-  Fluxes from TenStream solver
    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE ::  ts_abso  !< absorption in W/m3 [nlev_merged(-1), nxp, nyp]
    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE ::  ts_edir  !< direct SW flux in W/m2 [nlev_merged(-1), nxp, nyp]
    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE ::  ts_edn   !< incoming diffuse radiation in W/m2 [nlev_merged(-1), nxp, nyp]
    REAL(IREALS), DIMENSION(:,:,:), ALLOCATABLE ::  ts_eup   !< outgoing diffuse radiation in W/m2 [nlev_merged(-1), nxp, nyp]
!
!-  Optical properties due to vegetation
    REAL(IREALS), DIMENSION(:,:,:,:), ALLOCATABLE :: tree_tau_solar    !< optical properties due to vegetation, tau for solar
    REAL(IREALS), DIMENSION(:,:,:,:), ALLOCATABLE :: tree_w0_solar     !< w0, and
    REAL(IREALS), DIMENSION(:,:,:,:), ALLOCATABLE :: tree_tau_thermal  !< tau for thermal
!
!-  Pressure and temperature fields for TenStream
    REAL(IREALS), POINTER, DIMENSION(:,:) ::  pplev  !< pressure on layer interfaces [hPa] in zu grid
    REAL(IREALS), POINTER, DIMENSION(:,:) ::  ptlay  !< temperature on layer interfaces [K]
    REAL(IREALS), POINTER, DIMENSION(:,:) ::  ptlev  !< temperature on level interfaces [K]

!
!-  Types for buildings data structure
    TYPE(T_PPRTS_BUILDINGS), ALLOCATABLE ::  buildings_solar
    TYPE(T_PPRTS_BUILDINGS), ALLOCATABLE ::  buildings_thermal

!
!-  Number of facad
    INTEGER(iwp) ::  nfacad         !< number of facad
    INTEGER(iwp) ::  nfacad_east    !< number of facad at east boundary
    INTEGER(iwp) ::  nfacad_eastg   !< number of facad at east neighbour PE
    INTEGER(iwp) ::  nfacad_north   !< number of facad at north boundary
    INTEGER(iwp) ::  nfacad_northg  !< number of facad at neighbour north PE
    INTEGER(iwp) ::  nfacad_south   !< number of facad at south boundary
    INTEGER(iwp) ::  nfacad_southg  !< number of facad at neighbour south PE
    INTEGER(iwp) ::  nfacad_west    !< number of facad at west boundary
    INTEGER(iwp) ::  nfacad_westg   !< number of facad at west neighbour PE
!
!-  Variables for send/receive signals of the data exchange
    INTEGER(iwp) ::  requests(8)
    INTEGER(iwp) ::  request_count
!
!-  IDs for surface direction
    INTEGER(iwp), PARAMETER ::  iup_l    = 0  !< ID for land up-surface
    INTEGER(iwp), PARAMETER ::  idown_l  = 1  !< ID for land down-surface
    INTEGER(iwp), PARAMETER ::  ieast_l  = 2  !< ID for land east-surface
    INTEGER(iwp), PARAMETER ::  iwest_l  = 3  !< ID for land west-surface
    INTEGER(iwp), PARAMETER ::  inorth_l = 4  !< ID for land north-surface
    INTEGER(iwp), PARAMETER ::  isouth_l = 5  !< ID for land south-surface
    INTEGER(iwp), PARAMETER ::  iup_u    = 6  !< ID for urban up-surface
    INTEGER(iwp), PARAMETER ::  idown_u  = 1  !< ID for urban down-surface
    INTEGER(iwp), PARAMETER ::  ieast_u  = 7  !< ID for urban east-surface
    INTEGER(iwp), PARAMETER ::  iwest_u  = 8  !< ID for urban west-surface
    INTEGER(iwp), PARAMETER ::  inorth_u = 9  !< ID for urban north-surface
    INTEGER(iwp), PARAMETER ::  isouth_u = 10 !< ID for urban south-surface
!
!-  Variables for send/receive signals of the data exchange
    INTEGER(iwp), PARAMETER ::  tag_e = 1  !< tage for east  side
    INTEGER(iwp), PARAMETER ::  tag_n = 3  !< tage for north side
    INTEGER(iwp), PARAMETER ::  tag_s = 4  !< tage for east  side
    INTEGER(iwp), PARAMETER ::  tag_w = 2  !< tage for east  side
!
!-  Arrays for faces IDs at boarders
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  surf_ids     !< faces belonging to this PE (without those located at boarders)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  surf_ids_e   !< faces at east  boarder belonging to this PE
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  surf_ids_eg  !< faces at east  boarder belonging to east neighbour PE
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  surf_ids_n   !< faces at north boarder belonging to this PE
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  surf_ids_ng  !< faces at north boarder belonging to north neighbour PE
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  surf_ids_s   !< faces at south boarder belonging to this PE
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  surf_ids_sg  !< faces at south boarder belonging to south neighbour PE
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  surf_ids_w   !< faces at west  boarder belonging to this PE
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  surf_ids_wg  !< faces at west  boarder belonging to west neighbour PE
!
!-  Arrays for faces properties at boarders
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  surf_e   !< faces at east  boarder belonging to this PE
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  surf_eg  !< faces at east  boarder belonging to east neighbour PE
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  surf_n   !< faces at north boarder belonging to this PE
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  surf_ng  !< faces at north boarder belonging to north neighbour PE
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  surf_s   !< faces at south boarder belonging to this PE
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  surf_sg  !< faces at south boarder belonging to south neighbour PE
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  surf_w   !< faces at west  boarder belonging to this PE
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  surf_wg  !< faces at west  boarder belonging to west neighbour PE
!
!-  Data structure
    CLASS(T_SOLVER), ALLOCATABLE ::  ts_solver
    TYPE(T_TENSTR_ATM) ::  ts_atm
#endif
!
!-- Parameters of urban and land surface models
    INTEGER(iwp) ::  nz_urban    !< number of layers of urban surface (will be calculated)
    INTEGER(iwp) ::  nz_plant    !< number of layers of plant canopy (will be calculated)
    INTEGER(iwp) ::  nz_urban_b  !< bottom layer of urban surface (will be calculated)
    INTEGER(iwp) ::  nz_urban_t  !< top layer of urban surface (will be calculated)
    INTEGER(iwp) ::  nz_plant_t  !< top layer of plant canopy (will be calculated)
!
!-- Parameters of urban and land surface models
    INTEGER(iwp), PARAMETER ::  nzut_free = 3   !< number of free layers above top of of topography
    INTEGER(iwp), PARAMETER ::  ndsvf = 2       !< number of dimensions of real values in SVF
    INTEGER(iwp), PARAMETER ::  idsvf = 2       !< number of dimensions of integer values in SVF
    INTEGER(iwp), PARAMETER ::  ndcsf = 1       !< number of dimensions of real values in CSF
    INTEGER(iwp), PARAMETER ::  idcsf = 2       !< number of dimensions of integer values in CSF
    INTEGER(iwp), PARAMETER ::  kdcsf = 4       !< number of dimensions of integer values in CSF calculation array
    INTEGER(iwp), PARAMETER ::  id = 1          !< position of d-index in surfl and surf
    INTEGER(iwp), PARAMETER ::  iz = 2          !< position of k-index in surfl and surf
    INTEGER(iwp), PARAMETER ::  iy = 3          !< position of j-index in surfl and surf
    INTEGER(iwp), PARAMETER ::  ix = 4          !< position of i-index in surfl and surf
    INTEGER(iwp), PARAMETER ::  nidx_surf = 4   !< number of indices in surfl and surf
    INTEGER(iwp), PARAMETER ::  nsurf_type = 5  !< number of surf types = surface directions
    INTEGER(iwp), PARAMETER ::  iup    = 0      !< 0 - index of upward surface (ground or roof)
    INTEGER(iwp), PARAMETER ::  idown  = 1      !< 1 - index of downward surface (overhanging)
    INTEGER(iwp), PARAMETER ::  inorth = 2      !< 2 - index of northward facing wall
    INTEGER(iwp), PARAMETER ::  isouth = 3      !< 3 - index of southward facing wall
    INTEGER(iwp), PARAMETER ::  ieast  = 4      !< 4 - index of eastward facing wall
    INTEGER(iwp), PARAMETER ::  iwest  = 5      !< 5 - index of westward facing wall

    INTEGER(iwp), DIMENSION(0:nsurf_type), PARAMETER ::  idir = (/0, 0, 0, 0, 1,-1/)  !< surface normal direction x indices
    INTEGER(iwp), DIMENSION(0:nsurf_type), PARAMETER ::  jdir = (/0, 0, 1,-1, 0, 0/)  !< surface normal direction y indices
    INTEGER(iwp), DIMENSION(0:nsurf_type), PARAMETER ::  kdir = (/1,-1, 0, 0, 0, 0/)  !< surface normal direction z indices

    REAL(wp), DIMENSION(0:nsurf_type) ::  facearea  !< area of single face in respective, direction (will be calc'd)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::    discr_azim_cent  !< centres of discretized azimuths
    REAL(wp), DIMENSION(:), ALLOCATABLE ::    discr_azim_bdry  !< boundaries of discretized azimuths
    REAL(wp), DIMENSION(:), ALLOCATABLE ::    discr_elev_cent  !< centres of discretized elevation angles
    REAL(wp), DIMENSION(:), ALLOCATABLE ::    discr_elev_bdry  !< boundaries of discretized elevation angles
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  discr_azim_yxdir !< horizontal directions (normal vector in grid
                                                               !< units) of discretized azimuths

!
!-- Indices needed for RTM netcdf output subroutines
    INTEGER(iwp), PARAMETER ::  nd = 6  !< number of directions

    CHARACTER(LEN=6), DIMENSION(0:nd-1), PARAMETER ::  dirname = (/ '_up   ', '_down ', '_south', '_north', '_west ', '_east ' /) !<

    INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER ::  dirint = (/ iup, idown, isouth, inorth, iwest, ieast /)  !< direction integers
    INTEGER(iwp), DIMENSION(0:nd-1), PARAMETER ::  diridx =  (/  0,     1,      1,      0,     3,     2 /)  !< mapping to surf_h
                                                                                                            !< and surf_v
!
!-- Indices and sizes of urban and land surface models
    INTEGER(iwp) ::  nsurfl  !< number of all surfaces in local processor
    INTEGER(iwp) ::  nsurf   !< global number of surfaces in index array of surfaces (nsurf = proc nsurfs)

    INTEGER(iwp), DIMENSION(:,:), POINTER ::  surf   !< coordinates of i-th surface in grid - surf[:,k] = [d, z, y, x, m]
    INTEGER(iwp), DIMENSION(:,:), POINTER ::  surfl  !< coordinates of i-th local surface in local grid - surfl[:,k] =
                                                     !< [d, z, y, x, m]

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  surfl_linear     !< dtto (linearly allocated array)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  surf_linear      !< dtto (linearly allocated array)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  nsurfs           !< array of number of all surfaces in individual processors
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  surfstart        !< starts of blocks of surfaces for individual
                                                                        !< processors in array surf (indexed from 1)
                                                                        !< respective block for particular processor is
                                                                        !< surfstart[iproc+1]+1 : surfstart[iproc+1]+nsurfs[iproc+1]
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  surfl_col_start  !< start of surfaces in surfl
                                                                        !< for each x,y column (local surfaces)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  surfg_col_start  !< start of surfaces in surfl
                                                                        !< for each x,y column (all surfaces)
!
!-- Block variables needed for calculation of the plant canopy model inside the urban surface model
    INTEGER(iwp) ::  npcbl = 0  !< number of the plant canopy gridboxes in local processor

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  pct   !< top layer of the plant canopy
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  pch   !< heights of the plant canopy
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  pcbl  !< k,j,i coordinates of l-th local plant canopy box pcbl[:,l] = [k, j, i]

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcinsw      !< array of received sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcinswdir   !< array of received direct sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcinswdif   !< array of received diffuse sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcbinsw     !< array of absorbed sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcbinswdir  !< array of absorbed direct sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcbinswdif  !< array of absorbed diffusion sw radiation for local plant canopy box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcbinlw     !< array of absorbed lw radiation for local plant canopy box
!
!-- block of indices used during MPI calls
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE  ::  nnxy             !< numbers of PE subdomans xy grids
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE  ::  nnxyd            !< displacements of gathered values
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE  ::  ipx              !< index of gridcell processor along x axis
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE  ::  ipy              !< index of gridcell processor along y axis

!
!-- Configuration parameters (they can be setup in PALM config)
    INTEGER(iwp), PARAMETER ::  rad_version_len = 10  !< length of identification string of rad version

    CHARACTER(rad_version_len), PARAMETER ::  rad_version = 'RAD v. 5.1'  !< identification of version of binary svf and
                                                                          !< restart files

    INTEGER(iwp) ::  bufsize_alltoall = 0          !< max no. of items to send in MPI_ALLTOALL at once (0=infinite)
    INTEGER(iwp) ::  mrt_minlevel = 0              !< minumum vertical box above surface for which to calculate MRT
    INTEGER(iwp) ::  mrt_nlevels = 0               !< number of vertical boxes above surface for which to calculate MRT
    INTEGER(iwp) ::  nrefsteps = 3                 !< number of reflection steps to perform
    INTEGER(iwp) ::  raytrace_discrete_elevs = 40  !< number of discretization steps for elevation (nadir to zenith)
    INTEGER(iwp) ::  raytrace_discrete_azims = 80  !< number of discretization steps for azimuth (out of 360 degrees)

    INTEGER(wp) ::  mrt_geom = 1  !< method for MRT direction weights simulating a sphere or a human body

    LOGICAL ::  localized_raytracing = .FALSE.       !< new parallel raytracing algorithm with localized computation using
                                                     !< MPI-based request queue, avoiding one-sided MPI operations
    LOGICAL ::  mrt_skip_roof = .TRUE.               !< do not calculate MRT above roof surfaces
    LOGICAL ::  mrt_include_sw = .TRUE.              !< should MRT calculation include SW radiation as well?
    LOGICAL ::  plant_lw_interact = .TRUE.           !< whether plant canopy interacts with LW radiation (in addition to SW)
    LOGICAL ::  raytrace_mpi_rma = .TRUE.            !< use MPI RMA to access LAD and gridsurf from remote processes
                                                     !< during raytracing
    LOGICAL ::  radiation_volumetric_flux = .FALSE.  !< flag indicating whether volumetric radiative fluxes will be calculated

    REAL(wp), PARAMETER ::  ext_coef = 0.6_wp        !< extinction coefficient (a.k.a. alpha)
    REAL(wp), PARAMETER ::  min_opaque_lad = 0.5_wp  !< minimum value of LAD where trees are considered opaque for the
                                                     !< purpose of volumetric fluxes

    REAL(wp), DIMENSION(2) ::  mrt_geom_params = (/ .12_wp, .88_wp /)  !< parameters for the selected method

!
!-- Radiation related arrays to be used in radiation_interaction routine.
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_in_dir   !< direct sw radiation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_sw_in_diff  !< diffusion sw radiation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rad_lw_in_diff  !< diffusion lw radiation

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  rad_shade_h   !< Height of shadow
!
!-- Parameters required for average_radiation lower boundary conditions. Set reasonable
!-- default values which will be overwritten later on.
    REAL(wp) ::  albedo_eff = 0.1_wp      !< effective albedo value for average_radiation
    REAL(wp) ::  emissivity_eff = 0.9_wp  !< effective emissivity value for average_radiation
    REAL(wp) ::  t_rad_eff = 300.0_wp     !< effective radiative surface temperature for average_radiation
!
!-- Type for calculation of svf
    TYPE t_svf
        INTEGER(iwp) :: isurflt  !<
        INTEGER(iwp) :: isurfs   !<
        REAL(wp) :: rsvf     !<
        REAL(wp) :: rtransp  !<
    END TYPE
!
!-- Type for calculation of csf
    TYPE t_csf
        INTEGER(iwp) :: ip      !<
        INTEGER(iwp) :: itx     !<
        INTEGER(iwp) :: ity     !<
        INTEGER(iwp) :: itz     !<
        INTEGER(iwp) :: isurfs  !< Idx of source face / -1 for sky
        REAL(wp) :: rcvf  !< Canopy view factor for faces / canopy sink factor for sky (-1)
    END TYPE
!
!-- Localized raytracing IPC message; record for the whole raytraced angular section
    TYPE t_traced_section
       CHARACTER ::  message  !< Message code, uppercase=finalize the respective task:
                              !< f = VF raytrace forward
                              !< b = VF raytrace backward
                              !< F = finalize 'f' (also when 'b' is skipped)
                              !< p = plant canopy direct raytrace (forward only)
                              !< P = finalize 'p'
                              !< m = MRT raytrace (forward only)
                              !< M = finalize 'm'
                              !< c = process completed
                              !< t = terminate

       INTEGER(iwp) ::  iaz               !< discretized azimuth index
       INTEGER(iwp) ::  iorig             !< index of origin face (global) or PCGB/MRTB (local)
       INTEGER(iwp) ::  lowest_free_ray   !< index of the lowest ray that is free (above mixed horizon)
       INTEGER(iwp) ::  lowest_mixed_ray  !< index of the lowest ray that is mixed (above full horizon)
       INTEGER(iwp) ::  nrays             !< number of rays (z directions) to raytrace

       INTEGER(iwp), DIMENSION(2) ::  dimnext  !< next dimension increments along path

       REAL(wp) ::  aorig     !< origin face area for csf
       REAL(wp) ::  lastdist  !< beginning of current crossing

       REAL(wp), DIMENSION(2) ::  dimnextdist  !< distance for each dimension increments
       REAL(wp), DIMENSION(3) ::  origin       !< z,y,x coordinates of ray origin
    END TYPE
!
!-- Localized raytracing IPC message; record for each traced ray within the angular section. Rays
!-- are indexed zenith to nadir top to bottom, so r(1) is the highest ray
    TYPE t_traced_ray
       INTEGER(iwp) ::  itarget  !< global index of the target face or <0 for sky

       REAL(wp) ::  vffrac     !< view factor fraction of each ray for csf
       REAL(wp) ::  transp     !< transparency of the whole path
       REAL(wp) ::  zdir       !< z directions to raytrace (z/hdist in grid coords)
    END TYPE
!
!-- Arrays storing the values of USM
    INTEGER(iwp) ::  ndsidir  !< number of apparent solar directions used
    INTEGER(iwp) ::  nmrtbl   !< No. of local grid boxes for which MRT is calculated
    INTEGER(iwp) ::  nmrtf    !< number of MRT factors for local processor

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  dsidir_rev  !< dsidir_rev[ielev,iazim] = i for dsidir or -1 if not present
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  mrtbl       !< coordinates of i-th local MRT box - surfl[:,i] = [z, y, x]
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  mrtfsurf    !< mrtfsurf[:,imrtf] = index of target MRT box and
                                                              !< source surface for mrtf[imrtf]
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  svfsurf     !< svfsurf[:,isvf] = index of target and source surface for svf[isvf]

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  shadow_top  !< shadow_top(j,i,idir) = k; k and below are shaded

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrtf         !< array of MRT factors for each local MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrtft        !< array of MRT factors including transparency for each local MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrtsky       !< array of sky view factor for each local MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrtskyt      !< array of sky view factor including transparency for each local MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrtinsw      !< mean SW radiant flux for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrtinlw      !< mean LW radiant flux for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrt          !< mean radiant temperature for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrtinsw_av   !< time average mean SW radiant flux for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrtinlw_av   !< time average mean LW radiant flux for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mrt_av       !< time average mean radiant temperature for each MRT box
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsw    !< array of total sw radiation outgoing from nonvirtual surfaces
                                                         !< surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfins      !< array of sw radiation falling to local surface after i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinl      !< array of lw radiation for local surface after i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  skyvf        !< array of sky view factor for each local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  skyvft       !< array of sky view factor including transparency for each local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinsw     !< array of sw radiation falling to local surface including
                                                         !< radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlw     !< array of lw radiation falling to local surface including
                                                         !< radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswdir  !< array of direct sw radiation falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswdif  !< array of diffuse sw radiation from sky and model boundary
                                                         !< falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlwdif  !< array of diffuse lw radiation from sky and model boundary
                                                         !< falling to local surface
                                                         !< Outward radiation is only valid for nonvirtual surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsl    !< array of reflected sw radiation for local surface in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutll    !< array of reflected + emitted lw radiation for local
                                                         !< surface in i-th reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutlw    !< array of total lw radiation outgoing from nonvirtual surfaces
                                                         !< surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfemitlwl  !< array of emitted lw radiation for local surface used to calculate
                                                         !<effective surface temperature for radiation model

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dsitrans   !< dsidir[isvfl,i] = path transmittance of i-th
                                                         !< direction of direct solar irradiance per target surface
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dsitransc  !< dtto per plant canopy box
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dsidir     !< dsidir[:,i] = unit vector of i-th
                                                         !< direction of direct solar irradiance
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  mrtdsit    !< array of direct solar transparencies for each local MRT box
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  svf        !< array of shape view factors+direct irradiation factors
                                                         !< for local surfaces

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  skyvf_vol  !< volumetric sky view factor at (k,j,i)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  swflux_vol !< omnidirectional volumetric SW rad flux through (k,j,i)
                                                          !< as average W/m2 onto imaginary sphere

!
!-- MPI Optimization global variables
#if defined( __parallel )
    INTEGER(iwp) ::  niters_radx    !< number of partial MPI_ALLTOALLV loops for radiation exchange
    INTEGER(iwp) ::  nmaxsend_radx  !< maximum size of MPI_ALLTOALLV sent array per processor
    INTEGER(iwp) ::  nrecv_radx     !< total number of received surfaces in MPI_ALLTOALLV radiation exchange
    INTEGER(iwp) ::  nsend_radx     !< total number of sent surfaces in MPI_ALLTOALLV radiation exchange

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  disp_recvbuf_radx  !< displacement of receive buffers in MPI_ALLTOALLV
                                                                   !< radiation exchange
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  disp_sendbuf_radx  !< displacement of send buffers in MPI_ALLTOALLV
                                                                   !< radiation exchange
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  isurf_recv_radx    !< surface ids received in MPI_ALLTOALLV rad exchange
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  isurf_send_radx    !< surface ids sent in MPI_ALLTOALLV rad exchange
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  num_recvbuf_radx   !< size of receive buffers in MPI_ALLTOALLV radiation
                                                                   !< exchange
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  num_sendbuf_radx   !< size of send buffers in MPI_ALLTOALLV rad exchange
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  disp_recv_radx     !< displacements of surfaces received from other
                                                                   !< procs in MPI_ALLTOALLV rad exchange
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  disp_send_radx     !< displacements of surfaces sent to other procs in
                                                                   !< MPI_ALLTOALLV rad exchange

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  radx_send          !< send buffer for MPI_ALLTOALLV rad exchange
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  radx_send_surfinl  !< send buffer for incoming LW radiation from plant
                                                               !< canopy in MPI_ALLTOALLV rad exchange
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinl_recv       !< surface incoming LW radiation from plant canopy
                                                               !< received from other procs
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutl_recv      !< surface outgoing LW rad received from other procs
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfouts_recv      !< surface outgoing SW rad received from other procs
#endif

!
!-- Block variables needed for calculation of the plant canopy model inside the urban surface model
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  csfsurf  !< csfsurf[:,icsf] = index of target surface and
                                                           !< csf grid index for csf[icsf]

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  csf  !< array of plant canopy sink fators + direct irradiation factors (transparency)

    REAL(wp), DIMENSION(:,:,:), POINTER ::  sub_lad  !< subset of lad_s within urban surface, transformed to plain Z coordinate
#if defined( __parallel )
    REAL(wp), DIMENSION(:), POINTER ::  sub_lad_g  !< sub_lad globalized (used to avoid MPI RMA calls in raytracing)
#endif
    INTEGER(iwp) ::  plantt_max  !<
!
!-- Temporary global arrays for raytracing
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nzterrt     !< global 2-d terrain top (full-3d)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nzterrb     !< global 2-d terrain bottom (full-3d)
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  plantt      !< global 2-d plant top
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  opaque_top  !< global 2-d top of terrain + opaque plant canopy
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  opaque_top_l_lin  !< local linear array of opaque_top
    INTEGER(iwp), DIMENSION(:,:), POINTER ::  opaque_top_l  !< local 2-d pointer to opaque_top_l_lin

    REAL(wp) ::  prototype_lad  !< prototype leaf area density for computing effective optical depth

!
!-- Arrays and variables for calculation of svf and csf
    INTEGER(iwp), PARAMETER ::  gasize = 100000   !< initial size of growing arrays
    INTEGER(iwp), PARAMETER ::  nsurf_type_u = 6  !< number of urban surf types (used in gridsurf)

    INTEGER(iwp) ::  nsvfl              !< number of svf for local processor
    INTEGER(iwp) ::  nsvf_ins           !< number of svf for inserting currently in itarget/vffrac/ztransp
    INTEGER(iwp) ::  ncsfl              !< no. of csf in local processor needed only during calc_svf but must be here because it is
                                        !< shared between subroutines calc_svf and raytrace
    INTEGER(iwp) ::  nsvfla             !< dimmension of array allocated for storage of svf in local processor
    INTEGER(iwp) ::  ncsfla             !< dimmension of array allocated for storage of csf in local processor
    INTEGER(iwp) ::  nmrtfa             !< dimmension of array allocated for storage of mrt
    INTEGER(iwp) ::  msvf, mcsf, mmrtf  !< mod for swapping the growing array

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  itarget  !< face indices of detected obstacles

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  gridpcbl  !< reverse index of local pcbl[k,j,i]

    INTEGER(iwp), DIMENSION(:,:,:,:), POINTER ::  gridsurf  !< reverse index of local surfl[d,k,j,i]

    REAL(wp), PARAMETER ::  grow_factor = 1.4_wp  !< growth factor of growing arrays

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  vffrac   !< view factor fractions for individual rays
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ztransp  !< array of transparency in z steps

    TYPE(t_svf), DIMENSION(:), POINTER ::  asvf   !< pointer to growing svc array
    TYPE(t_csf), DIMENSION(:), POINTER ::  acsf   !< pointer to growing csf array
    TYPE(t_svf), DIMENSION(:), POINTER ::  amrtf  !< pointer to growing mrtf array

    TYPE(t_svf), DIMENSION(:), ALLOCATABLE, TARGET ::  asvf1, asvf2    !< realizations of svf array
    TYPE(t_csf), DIMENSION(:), ALLOCATABLE, TARGET ::  acsf1, acsf2    !< realizations of csf array
    TYPE(t_svf), DIMENSION(:), ALLOCATABLE, TARGET ::  amrtf1, amrtf2  !< realizations of mftf array

#if defined( __parallel )
    CHARACTER, DIMENSION(:), ALLOCATABLE, TARGET ::  lrt_msg_buf !< buffer for outgoing messages
#endif

    INTEGER(iwp)       ::  lrt_msg_type_incoming
    INTEGER(iwp)       ::  lrt_msg_type_processing
#if defined( __parallel )
    INTEGER(iwp)       ::  lrt_req_incoming         !< MPI request for incoming message
    INTEGER(iwp)       ::  lrt_wait_time            !< total time spent waiting for LRT messages (rapid updates)
    INTEGER, PARAMETER ::  lrt_msg_tag = 200        !< MPI message tag for lrt messages

    LOGICAL ::  lrt_prev_process_complete !< previous process has completed all raytracing
    LOGICAL ::  lrt_local_complete        !< local process has completed all raytracing
    LOGICAL ::  lrt_unfinished            !< not yet received termination message (still something to do)
#endif

    TYPE(t_traced_section), POINTER ::  lrt_msg_incoming_section    !<
    TYPE(t_traced_section), POINTER ::  lrt_msg_processing_section  !<

    TYPE(t_traced_ray), DIMENSION(:), POINTER ::  lrt_msg_incoming_rays    !<
    TYPE(t_traced_ray), DIMENSION(:), POINTER ::  lrt_msg_processing_rays  !<

!
!-- Temporary arrays for calculation of csf in raytracing
    INTEGER(iwp) ::  maxboxesg  !< max number of boxes ray can cross in the domain

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  boxes  !< coordinates of gridboxes being crossed by ray

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  crlens  !< array of crossing lengths of ray for particular grid boxes

#if defined( __parallel )
    INTEGER(iwp) ::  win_lad       !< MPI RMA window for leaf area density
    INTEGER(iwp) ::  win_gridsurf  !< MPI RMA window for reverse grid surface index

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  lad_ip  !< array of numbers of process where lad is stored

    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(:), ALLOCATABLE  ::  lad_disp  !< array of displaycements of lad
                                                                             !< in local array of proc lad_ip

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  lad_s_ray  !< array of received lad_s for appropriate gridboxes crossed by ray
#endif
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  target_surfl  !<

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  rt2_track  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rt2_track_dist  !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rt2_dist        !<

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  rt2_track_lad  !<
!
!-- Arrays for time averages
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcbinlw_av      !< Average of pcbinlw
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcbinsw_av      !< Average of pcbinsw
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcbinswdir_av   !< Average of pcbinswdir
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcbinswdif_av   !< Average of pcbinswdif
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcbinswref_av   !< Average of pcbinswref
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcinsw_av       !< Average of pcinsw
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcinswdir_av    !< Average of pcinswdir
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pcinswdif_av    !< Average of pcinswdif
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfradnet_av   !< average of net radiation to local surface including radiation
                                                            !< from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinsw_av     !< average of sw radiation falling to local surface including
                                                            !< radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlw_av     !< average of lw radiation falling to local surface including
                                                            !< radiation from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswdir_av  !< average of direct sw radiation falling to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswdif_av  !< average of diffuse sw radiation from sky and model boundary falling
                                                            !< to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlwdif_av  !< average of diffuse lw radiation from sky and model boundary falling
                                                            !< to local surface
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinswref_av  !< average of sw radiation falling to surface from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinlwref_av  !< average of lw radiation falling to surface from reflections
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutsw_av    !< average of total sw radiation outgoing from nonvirtual surfaces
                                                            !<  surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfoutlw_av    !< average of total lw radiation outgoing from nonvirtual surfaces
                                                            !< surfaces after all reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfins_av      !< average of array of residua of sw radiation absorbed in surface
                                                            !< after last reflection
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surfinl_av      !< average of array of residua of lw radiation absorbed in surface
                                                            !< after last reflection


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- Energy balance variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-- Parameters of the land, roof and wall surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  albedo_surf  !< albedo of the surface
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  emiss_surf   !< emissivity of the wall surface
!
!-- External radiation. Depending on the given level of detail either a 1D or a 3D array will be
!-- allocated.
    TYPE( real_1d_3d ) ::  rad_lw_in_f      !< external incoming longwave radiation, from observation or model
    TYPE( real_1d_3d ) ::  rad_sw_in_f      !< external incoming shortwave radiation, from observation or model
    TYPE( real_1d_3d ) ::  rad_sw_in_dif_f  !< external incoming shortwave radiation, diffuse part, from observation or model
    TYPE( real_1d_3d ) ::  time_rad_f       !< time dimension for external radiation, from observation or model
!
!-- Variables required for optimized(vectorized)-version of radiation_interaction.
    INTEGER(iwp)                              ::  nr_blocks = 0  !< number of ipcgb blocks
    INTEGER(iwp)                              ::  nr_surf_bl     !< number of blocsk of sorted plant-surface interaction
    INTEGER(iwp)                              ::  nr_surf_bl_2   !< number of blocks of sorted surface-surface interaction
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:)   ::  ipcgb_start    !< start indices of ipcgb blocks
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:)   ::  ipcgb_end      !< end indices of ipcgb blocks
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:)   ::  surf_end       !< end indices of surface index in sorted plant-canopy surface interaction
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:)   ::  surf_end_2     !< end indices of surface index in sorted surface-surface interaction
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:)   ::  surf_start     !< start indices of surface index in sorted plant-canopy surface interaction
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:)   ::  surf_start_2   !< start indices of surface index in sorted surface-surface interaction
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  sorted_ipcgb   !< array to sort ipcgb values
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  sorted_surf    !< array to sort receiving surface-surface radiation
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  sorted_surf_2  !< array to sort surfaces according to plant-canopy surface interaction


    INTERFACE radiation_calc_diffusion_radiation
       MODULE PROCEDURE radiation_calc_diffusion_radiation
    END INTERFACE radiation_calc_diffusion_radiation

#if defined( __tenstream )
    INTERFACE radiation_calc_sundir
       MODULE PROCEDURE radiation_calc_sundir
    END INTERFACE radiation_calc_sundir
#endif

    INTERFACE radiation_calc_svf
       MODULE PROCEDURE radiation_calc_svf
    END INTERFACE radiation_calc_svf

    INTERFACE radiation_check_data_output
       MODULE PROCEDURE radiation_check_data_output
    END INTERFACE radiation_check_data_output

    INTERFACE radiation_check_data_output_pr
       MODULE PROCEDURE radiation_check_data_output_pr
    END INTERFACE radiation_check_data_output_pr

    INTERFACE radiation_check_data_output_surf
       MODULE PROCEDURE radiation_check_data_output_surf
    END INTERFACE radiation_check_data_output_surf

    INTERFACE radiation_check_data_output_ts
       MODULE PROCEDURE radiation_check_data_output_ts
    END INTERFACE radiation_check_data_output_ts

    INTERFACE radiation_check_parameters
       MODULE PROCEDURE radiation_check_parameters
    END INTERFACE radiation_check_parameters

    INTERFACE radiation_clearsky
       MODULE PROCEDURE radiation_clearsky
    END INTERFACE radiation_clearsky

    INTERFACE radiation_constant
       MODULE PROCEDURE radiation_constant
    END INTERFACE radiation_constant

    INTERFACE radiation_control
       MODULE PROCEDURE radiation_control
    END INTERFACE radiation_control

    INTERFACE radiation_surface_data_averaging
       MODULE PROCEDURE radiation_surface_data_averaging
    END INTERFACE radiation_surface_data_averaging

    INTERFACE radiation_data_output_2d
       MODULE PROCEDURE radiation_data_output_2d
    END INTERFACE radiation_data_output_2d

    INTERFACE radiation_data_output_3d
       MODULE PROCEDURE radiation_data_output_3d
    END INTERFACE radiation_data_output_3d

    INTERFACE radiation_data_output_mask
       MODULE PROCEDURE radiation_data_output_mask
    END INTERFACE radiation_data_output_mask

    INTERFACE radiation_data_output_surf
       MODULE PROCEDURE radiation_data_output_surf
    END INTERFACE radiation_data_output_surf

    INTERFACE radiation_define_netcdf_grid
       MODULE PROCEDURE radiation_define_netcdf_grid
    END INTERFACE radiation_define_netcdf_grid

    INTERFACE radiation_header
       MODULE PROCEDURE radiation_header
    END INTERFACE radiation_header

    INTERFACE radiation_init
       MODULE PROCEDURE radiation_init
    END INTERFACE radiation_init

    INTERFACE radiation_interaction
       MODULE PROCEDURE radiation_interaction
    END INTERFACE radiation_interaction

    INTERFACE radiation_interaction_init
       MODULE PROCEDURE radiation_interaction_init
    END INTERFACE radiation_interaction_init

    INTERFACE radiation_parin
       MODULE PROCEDURE radiation_parin
    END INTERFACE radiation_parin

    INTERFACE radiation_presimulate_solar_pos
       MODULE PROCEDURE radiation_presimulate_solar_pos
    END INTERFACE radiation_presimulate_solar_pos

    INTERFACE radiation_read_svf
       MODULE PROCEDURE radiation_read_svf
    END INTERFACE radiation_read_svf

    INTERFACE radiation_rrd_global
       MODULE PROCEDURE radiation_rrd_global_ftn
       MODULE PROCEDURE radiation_rrd_global_mpi
    END INTERFACE radiation_rrd_global

    INTERFACE radiation_rrd_local
       MODULE PROCEDURE radiation_rrd_local_ftn
       MODULE PROCEDURE radiation_rrd_local_mpi
    END INTERFACE radiation_rrd_local

    INTERFACE radiation_rrtmg
       MODULE PROCEDURE radiation_rrtmg
    END INTERFACE radiation_rrtmg

#if defined( __rrtmg ) || defined( __tenstream )
    INTERFACE radiation_tendency
       MODULE PROCEDURE radiation_tendency
       MODULE PROCEDURE radiation_tendency_ij
    END INTERFACE radiation_tendency
#endif

    INTERFACE radiation_tenstream
       MODULE PROCEDURE radiation_tenstream
    END INTERFACE radiation_tenstream

    INTERFACE radiation_tenstream_init
       MODULE PROCEDURE radiation_tenstream_init
    END INTERFACE radiation_tenstream_init

    INTERFACE radiation_statistics
       MODULE PROCEDURE radiation_statistics
    END INTERFACE radiation_statistics

    INTERFACE radiation_vm_sampling
       MODULE PROCEDURE radiation_vm_sampling
    END INTERFACE radiation_vm_sampling

    INTERFACE radiation_wrd_global
       MODULE PROCEDURE radiation_wrd_global
    END INTERFACE radiation_wrd_global

    INTERFACE radiation_wrd_local
       MODULE PROCEDURE radiation_wrd_local
    END INTERFACE radiation_wrd_local

    INTERFACE radiation_write_svf
       MODULE PROCEDURE radiation_write_svf
    END INTERFACE radiation_write_svf

    INTERFACE radiation_3d_data_averaging
       MODULE PROCEDURE radiation_3d_data_averaging
    END INTERFACE radiation_3d_data_averaging

    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC radiation_calc_diffusion_radiation,                                                     &
           radiation_calc_svf,                                                                     &
           radiation_check_data_output,                                                            &
           radiation_check_data_output_pr,                                                         &
           radiation_check_data_output_surf,                                                       &
           radiation_check_data_output_ts,                                                         &
           radiation_check_parameters,                                                             &
           radiation_control,                                                                      &
           radiation_data_output_2d,                                                               &
           radiation_data_output_3d,                                                               &
           radiation_data_output_mask,                                                             &
           radiation_data_output_surf,                                                             &
           radiation_define_netcdf_grid,                                                           &
           radiation_header,                                                                       &
           radiation_init,                                                                         &
           radiation_interaction,                                                                  &
           radiation_interaction_init,                                                             &
           radiation_parin,                                                                        &
           radiation_presimulate_solar_pos,                                                        &
           radiation_read_svf,                                                                     &
           radiation_rrd_global,                                                                   &
           radiation_rrd_local,                                                                    &
           radiation_surface_data_averaging,                                                       &
           radiation_statistics,                                                                   &
           radiation_vm_sampling,                                                                  &
           radiation_wrd_global,                                                                   &
           radiation_wrd_local,                                                                    &
           radiation_write_svf,                                                                    &
           radiation_3d_data_averaging

!
!-- Public variables and constants: some are only used in the radiation - NEEDS REVISION
    PUBLIC albedo,                                                                                 &
           albedo_type,                                                                            &
           albedo_eff,                                                                             &
           average_radiation,                                                                      &
           cos_zenith,                                                                             &
           calc_zenith,                                                                            &
           dcep_average_radiation,                                                                 &
           decl_1,                                                                                 &
           decl_2,                                                                                 &
           decl_3,                                                                                 &
           dirname,                                                                                &
           diridx,                                                                                 &
           dirint,                                                                                 &
           dt_radiation,                                                                           &
           emissivity,                                                                             &
           emissivity_eff,                                                                         &
           force_radiation_call,                                                                   &
           idir,                                                                                   &
           id,                                                                                     &
           iz,                                                                                     &
           iy,                                                                                     &
           ix,                                                                                     &
           iup,                                                                                    &
           idown,                                                                                  &
           inorth,                                                                                 &
           isouth,                                                                                 &
           ieast,                                                                                  &
           iwest,                                                                                  &
           idsvf,                                                                                  &
           idcsf,                                                                                  &
           jdir,                                                                                   &
           kdir,                                                                                   &
           kdcsf,                                                                                  &
           lat,                                                                                    &
           lon,                                                                                    &
           mrt_geom,                                                                               &
           mrt_geom_params,                                                                        &
           mrt_include_sw,                                                                         &
           mrt_minlevel,                                                                           &
           mrt_nlevels,                                                                            &
           mrtbl,                                                                                  &
           mrtinsw,                                                                                &
           mrtinlw,                                                                                &
           nd,                                                                                     &
           nmrtbl,                                                                                 &
           nsurf_type,                                                                             &
           nz_urban_b,                                                                             &
           nz_urban_t,                                                                             &
           nz_urban,                                                                               &
           nsurf,                                                                                  &
           ndsvf,                                                                                  &
           ndcsf,                                                                                  &
           pch,                                                                                    &
           pct,                                                                                    &
           pcinsw,                                                                                 &
           pcinswdir,                                                                              &
           pcinswdif,                                                                              &
           pcbl,                                                                                   &
           npcbl,                                                                                  &
           rad_net_av,                                                                             &
           radiation,                                                                              &
           radiation_scheme,                                                                       &
           radiation_interactions,                                                                 &
           radiation_interactions_on,                                                              &
           radiation_volumetric_flux,                                                              &
           rad_shade_h,                                                                            &
           rad_sw_in_diff,                                                                         &
           rad_sw_in_dir,                                                                          &
           rad_lw_in,                                                                              &
           rad_lw_in_av,                                                                           &
           rad_lw_in_diff,                                                                         &
           rad_lw_out,                                                                             &
           rad_lw_out_av,                                                                          &
           rad_lw_cs_hr,                                                                           &
           rad_lw_cs_hr_av,                                                                        &
           rad_lw_hr,                                                                              &
           rad_lw_hr_av,                                                                           &
           rad_sw_in,                                                                              &
           rad_sw_in_av,                                                                           &
           rad_sw_out,                                                                             &
           rad_sw_out_av,                                                                          &
           rad_sw_cs_hr,                                                                           &
           rad_sw_cs_hr_av,                                                                        &
           rad_sw_hr,                                                                              &
           rad_sw_hr_av,                                                                           &
           solar_constant,                                                                         &
           skip_time_do_radiation,                                                                 &
           skyvf,                                                                                  &
           skyvft,                                                                                 &
           sun_direction,                                                                          &
           sun_dir_lat,                                                                            &
           sun_dir_lon,                                                                            &
           time_radiation,                                                                         &
           t_rad_eff,                                                                              &
           unscheduled_radiation_calls

#if defined( __rrtmg )
    PUBLIC radiation_tendency,                                                                     &
           rrtm_aldif,                                                                             &
           rrtm_aldir,                                                                             &
           rrtm_asdif,                                                                             &
           rrtm_asdir
#endif

#if defined( __tenstream ) && !defined( __rrtmg )
    PUBLIC radiation_tendency
#endif

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine controls the calls of the radiation schemes
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_control


    IMPLICIT NONE


    IF ( debug_output_timestep )  CALL debug_message( 'radiation_control', 'start' )


    SELECT CASE ( TRIM( radiation_scheme ) )

       CASE ( 'constant' )
          CALL radiation_constant

       CASE ( 'clear-sky' )
          CALL radiation_clearsky

       CASE ( 'rrtmg' )
          CALL radiation_rrtmg

       CASE ( 'tenstream' )
          CALL radiation_tenstream

       CASE ( 'external' )
!
!--       During spinup apply clear-sky model
          IF ( time_since_reference_point < 0.0_wp )  THEN
             CALL radiation_clearsky
          ELSE
             CALL radiation_external
          ENDIF

       CASE DEFAULT

    END SELECT

    IF ( debug_output_timestep )  CALL debug_message( 'radiation_control', 'end' )

 END SUBROUTINE radiation_control


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for radiation model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_check_data_output( variable, unit, i, ilen, k )


    USE control_parameters,                                                                        &
        ONLY:  data_output,                                                                        &
               message_string

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  unit      !<
    CHARACTER(LEN=*) ::  variable  !<

    CHARACTER(LEN=varnamelength) ::   var  !< TRIM(variable)

    INTEGER(iwp) ::  i, k        !<
    INTEGER(iwp) ::  ilast_word  !<
    INTEGER(iwp) ::  ilen        !<
    INTEGER(iwp) ::  id          !<

    LOGICAL ::  directional  !<

    var = TRIM( variable )
!
!-- Identify directional variables
    ilast_word = SCAN( var, '_', back = .TRUE. )
    directional = .FALSE.
    IF ( ilast_word > 0 )  THEN
       DO  id = 0, nd-1
          IF ( TRIM( var(ilast_word:) ) == TRIM( dirname(id) ) )  THEN
             directional = .TRUE.
             var = var(1:ilast_word-1)
             EXIT
          ENDIF
       ENDDO
    ENDIF

    IF ( directional )  THEN
       IF ( var(1:8) == 'rtm_svf_'  .OR.  var(1:8) == 'rtm_dif_' )  THEN
          IF ( .NOT.  radiation )  THEN
             message_string = 'output of "' // var // '" requires radiation = .TRUE.'
             CALL message( 'radiation_check_data_output', 'RAD0001', 1, 2, 0, 6, 0 )
          ENDIF
          unit = '1'
       ELSE
          SELECT CASE ( TRIM( var ) )
             CASE ( 'rtm_rad_net', 'rtm_rad_insw', 'rtm_rad_inlw', 'rtm_rad_inswdir',              &
                    'rtm_rad_inswdif', 'rtm_rad_inswref', 'rtm_rad_inlwdif', 'rtm_rad_inlwref',    &
                    'rtm_rad_outsw', 'rtm_rad_outlw', 'rtm_rad_ressw', 'rtm_rad_reslw' )
                IF ( .NOT.  radiation )  THEN
                   message_string = 'output of "' // var // '" requires radiation = .TRUE.'
                   CALL message( 'radiation_check_data_output', 'RAD0001', 1, 2, 0, 6, 0 )
                ENDIF
                unit = 'W/m2'

             CASE ( 'rtm_surfalb', 'rtm_surfemis' )
                IF ( .NOT.  radiation )  THEN
                   message_string = 'output of "' // var // '" requires radiation = .TRUE.'
                   CALL message( 'radiation_check_data_output', 'RAD0001', 1, 2, 0, 6, 0 )
                ENDIF
                unit = '1'

             CASE DEFAULT
                unit = 'illegal'
          END SELECT
       ENDIF

    ELSE
       SELECT CASE ( var )
          CASE ( 'rad_lw_cs_hr', 'rad_sw_cs_hr' )
             IF ( .NOT.  radiation  .OR.  radiation_scheme /= 'rrtmg' )  THEN
                message_string = '"output of "' // var // '" requires radiation = .TRUE. and ' //  &
                                 'radiation_scheme = "rrtmg"'
                CALL message( 'radiation_check_data_output', 'RAD0002', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'K/h'

          CASE ( 'rad_lw_hr', 'rad_sw_hr' )
             IF ( .NOT.  radiation  .OR.  ( radiation_scheme /= 'rrtmg'  .AND.                     &
                                            radiation_scheme /= 'tenstream' ) )  THEN
                message_string = '"output of "' // var // '" requires radiation = .TRUE. and ' //  &
                                 'radiation_scheme = "rrtmg" or radiation_scheme = "tenstream"'
                CALL message( 'radiation_check_data_output', 'RAD0003', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'K/h'

          CASE ( 'rad_lw_in', 'rad_lw_out', 'rad_sw_in', 'rad_sw_out' )
             IF ( .NOT.  radiation  .OR.  ( radiation_scheme /= 'rrtmg'  .AND.                     &
                                            radiation_scheme /= 'tenstream' ) )  THEN
                message_string = '"output of "' // var // '" requires radiation = .TRUE. and ' //  &
                                 'radiation_scheme = "rrtmg" or radiation_scheme = "tenstream"'
                CALL message( 'radiation_check_data_output', 'RAD0003', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'W/m2'

          CASE ( 'rad_net*', 'rad_lw_in*', 'rad_lw_out*', 'rad_sw_in*', 'rad_sw_out*' )
             IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
                message_string = 'illegal value for data_output: "' // var //                      &
                                 '" & only 2d-horizontal cross sections are allowed for this value'
                CALL message( 'radiation_check_data_output', 'RAD0004', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( radiation_only )  THEN
                message_string = 'output of surface radiation cross sections of: "' // var //      &
                                 '&is not allowed in radiation-only mode.'
                CALL message( 'radiation_check_data_output', 'RAD0005', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'W/m2'

          CASE ( 'rrtm_aldif*', 'rrtm_aldir*', 'rrtm_asdif*', 'rrtm_asdir*' )
             IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
                message_string = 'illegal value for data_output: "' // var //                      &
                                 '" & only 2d-horizontal cross sections are allowed for this value'
                CALL message( 'radiation_check_data_output', 'RAD0004', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( .NOT.  radiation  .OR.  radiation_scheme /= "rrtmg" )  THEN
                message_string = 'output of "' // var // '" requires radiation = .TRUE. and ' //   &
                                 'radiation_scheme = "rrtmg"'
                CALL message( 'radiation_check_data_output', 'RAD0002', 1, 2, 0, 6, 0 )
             ENDIF
             unit = ''

          CASE ( 'rtm_rad_pc_inlw', 'rtm_rad_pc_insw', 'rtm_rad_pc_inswdir', 'rtm_rad_pc_inswdif', &
                 'rtm_rad_pc_inswref' )
             IF ( .NOT.  radiation )  THEN
                message_string = 'output of "' // var // '" requires radiation = .TRUE.'
                CALL message( 'radiation_check_data_output', 'RAD0001', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'W m-3'

          CASE ( 'rtm_rad_pc_sw_in', 'rtm_rad_pc_sw_dir', 'rtm_rad_pc_sw_dif' )
             IF ( .NOT.  radiation )  THEN
                message_string = 'output of "' // var // '" requires radiation = .TRUE.'
                CALL message( 'radiation_check_data_output', 'RAD0001', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'W/m2'

          CASE ( 'rtm_rad_vol_sw' )
             IF ( .NOT.  radiation )  THEN
                message_string = 'output of "' // var // '" requires radiation = .TRUE.'
                CALL message( 'radiation_check_data_output', 'RAD0001', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'W m-2'

          CASE ( 'rtm_mrt', 'rtm_mrt_sw', 'rtm_mrt_lw'  )
             IF ( .NOT.  radiation )  THEN
                message_string = 'output of "' // var // '" requires radiation = .TRUE.'
                CALL message( 'radiation_check_data_output', 'RAD0001', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( mrt_nlevels == 0 )  THEN
                message_string = 'output of "' // var // '" requires mrt_nlevels > 0'
                CALL message( 'radiation_check_data_output', 'RAD0006', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( var == 'rtm_mrt_sw'  .AND.  .NOT. mrt_include_sw )  THEN
                message_string = 'output of "' // var // '" requires rtm_mrt_sw = .TRUE.'
                CALL message( 'radiation_check_data_output', 'RAD0007', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( var == 'rtm_mrt' )  THEN
                unit = 'K'
             ELSE
                unit = 'W m-2'
             ENDIF

          CASE DEFAULT
             unit = 'illegal'

       END SELECT
    ENDIF

 END SUBROUTINE radiation_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set module-specific timeseries units and labels
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )


    INTEGER(iwp), INTENT(IN)    ::  dots_max  !<
    INTEGER(iwp), INTENT(INOUT) ::  dots_num  !<

    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_label  !<
    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_unit   !<


!
!-- RTM time series:
!-- For each time series quantity you have to give a label and a unit, which will be used for the
!-- NetCDF file. The value of dots_num has to be increased by the number of new time series
!-- quantities. The start index for RTM time series is stored in dots_start_index_rtm and later used
!-- to address the module specific time series values.
    dots_start_index_rtm = dots_num + 1

    dots_num = dots_num + 1
    dots_label(dots_num) = 'rad_net'
    dots_unit(dots_num)  = 'W/m2'

    dots_num = dots_num + 1
    dots_label(dots_num) = 'rad_lw_in'
    dots_unit(dots_num)  = 'W/m2'

    dots_num = dots_num + 1
    dots_label(dots_num) = 'rad_lw_out'
    dots_unit(dots_num)  = 'W/m2'

    dots_num = dots_num + 1
    dots_label(dots_num) = 'rad_sw_in'
    dots_unit(dots_num)  = 'W/m2'

    dots_num = dots_num + 1
    dots_label(dots_num) = 'rad_sw_out'
    dots_unit(dots_num)  = 'W/m2'

    IF ( radiation_scheme /= 'tenstream' )  THEN

       IF ( average_radiation )  THEN

          dots_num = dots_num + 1
          dots_label(dots_num) = 't_rad_eff'
          dots_unit(dots_num)  = 'K'

          dots_num = dots_num + 1
          dots_label(dots_num) = 'emiss_eff'
          dots_unit(dots_num)  = ''

          dots_num = dots_num + 1
          dots_label(dots_num) = 'albedo_eff'
          dots_unit(dots_num)  = ''

       ENDIF

       IF ( radiation_scheme == 'rrtmg' )  THEN

          dots_num = dots_num + 1
          dots_label(dots_num) = 'rrtm_aldif'
          dots_unit(dots_num)  = ''

          dots_num = dots_num + 1
          dots_label(dots_num) = 'rrtm_aldir'
          dots_unit(dots_num)  = ''

          dots_num = dots_num + 1
          dots_label(dots_num) = 'rrtm_asdif'
          dots_unit(dots_num)  = ''

          dots_num = dots_num + 1
          dots_label(dots_num) = 'rrtm_asdir'
          dots_unit(dots_num)  = ''

       ENDIF

    ENDIF

 END SUBROUTINE radiation_check_data_output_ts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for radiation model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_check_data_output_pr( variable, var_count, unit, dopr_unit )

    USE arrays_3d,                                                                                 &
        ONLY:  zu

    USE control_parameters,                                                                        &
        ONLY:  data_output_pr,                                                                     &
               message_string

    USE indices

    USE profil_parameter

    USE statistics

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  dopr_unit  !< local value of dopr_unit
    CHARACTER(LEN=*) ::  unit       !<
    CHARACTER(LEN=*) ::  variable   !<

    INTEGER(iwp) ::  var_count  !<

    SELECT CASE ( TRIM( variable ) )

      CASE ( 'rad_net' )
          IF ( (  .NOT.  radiation )  .OR.  radiation_scheme == 'constant' )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not available for radiation = .FALSE. or ' //                       &
                              'radiation_scheme = "constant"'
             CALL message( 'radiation_check_data_output_pr', 'RAD0008', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 99
             dopr_unit  = 'W/m2'
             hom(:,2,99,:)  = SPREAD( zw, 2, statistic_regions+1 )
             unit = dopr_unit
          ENDIF

       CASE ( 'rad_lw_in' )
          IF ( ( .NOT.  radiation)  .OR.  radiation_scheme == 'constant' )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not available for radiation = .FALSE. or ' //                       &
                              'radiation_scheme = "constant"'
             CALL message( 'radiation_check_data_output_pr', 'RAD0008', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 100
             dopr_unit  = 'W/m2'
             hom(:,2,100,:)  = SPREAD( zw, 2, statistic_regions+1 )
             unit = dopr_unit
          ENDIF

       CASE ( 'rad_lw_out' )
          IF ( ( .NOT. radiation )  .OR.  radiation_scheme == 'constant' )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not available for radiation = .FALSE. or ' //                       &
                              'radiation_scheme = "constant"'
             CALL message( 'radiation_check_data_output_pr', 'RAD0008', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 101
             dopr_unit  = 'W/m2'
             hom(:,2,101,:)  = SPREAD( zw, 2, statistic_regions+1 )
             unit = dopr_unit
          ENDIF

       CASE ( 'rad_sw_in' )
          IF ( (  .NOT. radiation )  .OR.  radiation_scheme == 'constant' )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not available for radiation = .FALSE. or ' //                       &
                              'radiation_scheme = "constant"'
             CALL message( 'radiation_check_data_output_pr', 'RAD0008', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 102
             dopr_unit  = 'W/m2'
             hom(:,2,102,:)  = SPREAD( zw, 2, statistic_regions+1 )
             unit = dopr_unit
          ENDIF

       CASE ( 'rad_sw_out')
          IF ( ( .NOT.  radiation )  .OR.  radiation_scheme == 'constant' )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not available for radiation = .FALSE. or ' //                       &
                              'radiation_scheme = "constant"'
             CALL message( 'radiation_check_data_output_pr', 'RAD0008', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 103
             dopr_unit  = 'W/m2'
             hom(:,2,103,:)  = SPREAD( zw, 2, statistic_regions+1 )
             unit = dopr_unit
          ENDIF

       CASE ( 'rad_lw_cs_hr' )
          IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not available for radiation = .FALSE. or ' //                       &
                              'radiation_scheme /= "rrtmg"'
             CALL message( 'radiation_check_data_output_pr', 'RAD0009', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 104
             dopr_unit  = 'K/h'
             hom(:,2,104,:)  = SPREAD( zu, 2, statistic_regions+1 )
             unit = dopr_unit
          ENDIF

       CASE ( 'rad_lw_hr' )
          IF ( ( .NOT.  radiation )  .OR.  ( radiation_scheme /= 'rrtmg'  .AND.                    &
                                            radiation_scheme /= 'tenstream' ) )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not available for radiation = .FALSE. or ' //                       &
                              'radiation_scheme /= "rrtmg" or radiation_scheme = "tenstream"'
             CALL message( 'radiation_check_data_output_pr', 'RAD0009', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 105
             dopr_unit  = 'K/h'
             hom(:,2,105,:)  = SPREAD( zu, 2, statistic_regions+1 )
             unit = dopr_unit
          ENDIF

       CASE ( 'rad_sw_cs_hr' )
          IF ( (  .NOT.  radiation )  .OR.  radiation_scheme /= 'rrtmg' )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not available for radiation = .FALSE. or ' //                       &
                              'radiation_scheme /= "rrtmg"'
             CALL message( 'radiation_check_data_output_pr', 'RAD0009', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 106
             dopr_unit  = 'K/h'
             hom(:,2,106,:)  = SPREAD( zu, 2, statistic_regions+1 )
             unit = dopr_unit
          ENDIF

       CASE ( 'rad_sw_hr' )
          IF ( (  .NOT.  radiation )  .OR.  ( radiation_scheme /= 'rrtmg'  .AND.                   &
                                            radiation_scheme /= 'tenstream' ) )  THEN
             message_string = 'data_output_pr = ' // TRIM( data_output_pr(var_count) ) // ' is' // &
                              'not available for radiation = .FALSE. or ' //                       &
                              'radiation_scheme /= "rrtmg" or radiation_scheme = "tenstream"'
             CALL message( 'radiation_check_data_output_pr', 'RAD0009', 1, 2, 0, 6, 0 )
          ELSE
             dopr_index(var_count) = 107
             dopr_unit  = 'K/h'
             hom(:,2,107,:)  = SPREAD( zu, 2, statistic_regions+1 )
             unit = dopr_unit
          ENDIF


       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE radiation_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check surface data output variables from the radiation model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_check_data_output_surf( trimvar, unit, av )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)    ::  trimvar  !< dummy for single output variable
    CHARACTER(LEN=*), INTENT(INOUT) ::  unit     !< dummy for unit of output variable

    INTEGER(iwp), INTENT(IN) ::  av  !< id indicating average or non-average data output


    SELECT CASE ( TRIM( trimvar ) )

       CASE ( 'rtm_skyvf', 'rtm_skyvft' )
          IF ( av == 1 )  THEN
             message_string = 'time averaging of static quantity "' // TRIM( trimvar ) //          &
                              '" is not provided'
             CALL message( 'radiation_check_data_output_surf', 'RAD0010', 1, 2, 0, 6, 0 )
          ENDIF
          unit = '1'

       CASE DEFAULT
           unit = 'illegal'

    END SELECT

 END SUBROUTINE radiation_check_data_output_surf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for radiation model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_check_parameters

    USE control_parameters,                                                                        &
        ONLY:  land_surface,                                                                       &
               message_string,                                                                     &
               urban_surface

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  input_pids_static

    IMPLICIT NONE

!
!-- In case no urban-surface or land-surface model is applied, usage of a radiation model makes
!-- no sense, except radiation-only is employed.
    IF ( .NOT. land_surface  .AND.  .NOT. urban_surface  .AND.  .NOT. radiation_only )  THEN
       message_string = 'usage of radiation model is only allowed if ' //                          &
                        'land-surface and/or urban-surface model is applied, &or if ' //           &
                        '"radiation_only = .T."'
       CALL message( 'radiation_check_parameters', 'RAD0011', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Further checks for radiation-only mode.
    IF ( ( land_surface  .OR.  urban_surface )  .AND.  radiation_only )  THEN
       message_string = 'usage of radiation model in radiation-only mode in combination with ' //  &
                        'an energy-balance model (LSM/USM) is not allowed.'
       CALL message( 'radiation_check_parameters', 'RAD0012', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( radiation_only  .AND.  radiation_interactions_on )  THEN
       message_string = 'usage of radiation model in radiation-only mode in combination with ' //  &
                        '"radiation_interactions_on = .T." &is not allowed.'
       CALL message( 'radiation_check_parameters', 'RAD0013', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( radiation_only  .AND.  radiation_scheme /= 'rrtmg' )  THEN
       message_string = 'usage of radiation model in radiation-only mode is only allowed ' //      &
                        'in combination with "radiation_scheme = rrtmg".'
       CALL message( 'radiation_check_parameters', 'RAD0014', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( radiation_only  .AND.  albedo_type == 9999999 )  THEN
       message_string = 'radiation_only = .T. requires explicit setting of albedo_type'
       CALL message( 'radiation_check_parameters', 'RAD0064', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( radiation_scheme /= 'constant'   .AND.  radiation_scheme /= 'clear-sky'  .AND.            &
         radiation_scheme /= 'rrtmg'      .AND.  radiation_scheme /= 'tenstream'  .AND.            &
         radiation_scheme /= 'external' )  THEN
       message_string = 'unknown radiation_scheme = '//  TRIM( radiation_scheme )
       CALL message( 'radiation_check_parameters', 'RAD0015', 1, 2, 0, 6, 0 )
    ELSEIF ( radiation_scheme == 'rrtmg' )  THEN
#if ! defined( __rrtmg )
       message_string = 'radiation_scheme = "rrtmg" requires compilation of PALM with ' //         &
                        'pre-processor directive -D__rrtmg'
       CALL message( 'radiation_check_parameters', 'RAD0016', 1, 2, 0, 6, 0 )
#endif
#if defined( __rrtmg ) && ! defined( __netcdf )
       message_string = 'radiation_scheme = "rrtmg" requires the use of NetCDF preprocessor ' //   &
                        'directive -D__netcdf'
       CALL message( 'radiation_check_parameters', 'RAD0017', 1, 2, 0, 6, 0 )
#endif
    ELSEIF ( radiation_scheme == 'tenstream' )  THEN
#if ! defined( __tenstream )
       message_string = 'radiation_scheme = "tenstream" requires compilation of PALM with ' //     &
                        'pre-processor directive -D__tenstream'
       CALL message( 'radiation_check_parameters', 'RAD0018', 1, 2, 0, 6, 0 )
#endif

    ENDIF
!
!-- Checks performed only if data is given via namelist only.
    IF ( .NOT. input_pids_static )  THEN
       IF ( albedo_type == 0  .AND.  albedo == 9999999.9_wp  .AND.          &
            radiation_scheme == 'clear-sky')  THEN
          message_string = 'radiation_scheme = "clear-sky" in combination with albedo_type = 0 ' //&
                           'requires setting of albedo /= 9999999.9'
          CALL message( 'radiation_check_parameters', 'RAD0019', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( albedo_type == 0  .AND.  ( radiation_scheme == 'rrtmg'  .OR.                           &
                                       radiation_scheme == 'tenstream' )  .AND.                    &
          ( albedo_lw_dif == 9999999.9_wp .OR. albedo_lw_dir == 9999999.9_wp  .OR.                 &
            albedo_sw_dif == 9999999.9_wp .OR. albedo_sw_dir == 9999999.9_wp ) )  THEN
          message_string = 'radiation_scheme = "rrtmg"/"tenstream" in combination with ' //        &
                           'albedo_type = 0 requires setting of albedo_lw_dif /= 9999999.9' //     &
                           'albedo_lw_dir /= 9999999.9 albedo_sw_dif /= 9999999.9 and' //          &
                           'albedo_sw_dir /= 9999999.9'
          CALL message( 'radiation_check_parameters', 'RAD0020', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- Check if albedo_pars is correctly dimensioned.
    IF ( albedo_pars_f%from_file )  THEN
       IF ( SIZE( albedo_pars_f%pars ) /= dim_albedo_pars )  THEN
          WRITE( message_string, * ) 'Dimension size of static input variable albedo_pars is ',    &
                                     SIZE( albedo_pars_f%pars ), '.&',                             &
                                     'Dimension size of ', dim_albedo_pars, 'is required.'

          CALL message( 'radiation_check_parameters', 'RAD0021', 2, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- Parallel angular discretization without raytrace_mpi_rma is not possible
!-- Serial mode does not allow mpi_rma
#if defined( __parallel )
    IF ( surface_reflections  .AND.  .NOT. raytrace_mpi_rma )  THEN
       message_string = 'surface_reflections can only be used together with ' //            &
                        'raytrace_mpi_rma or when no parallelization is applied.'
       CALL message( 'radiation_check_parameters', 'RAD0022', 1, 2, 0, 6, 0 )
    ENDIF
#else
    IF ( raytrace_mpi_rma )  THEN
       message_string = 'raytrace_mpi_rma = .T. not allowed in serial mode'
       CALL message( 'radiation_check_parameters', 'RAD0023', 1, 2, 0, 6, 0 )
    ENDIF
#endif

    IF ( cloud_droplets  .AND.   ( radiation_scheme == 'rrtmg'  .OR.                               &
         radiation_scheme == 'tenstream' )  .AND.  average_radiation )  THEN
       message_string = 'average_radiation = .T. with radiation_scheme = "rrtmg"/"tenstream" ' //  &
                        'in combination with cloud_droplets = .T. is not implemented'
       CALL message( 'radiation_check_parameters', 'RAD0024', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check for dt_radiation
    IF ( dt_radiation <= 0.0 )  THEN
       message_string = 'dt_radiation must be > 0.0'
       CALL message( 'radiation_check_parameters', 'RAD0025', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check for the angular discretization parameters
!-- Error message when too small values are set
    IF ( raytrace_discrete_elevs < 4  .OR.  raytrace_discrete_azims < 8 )  THEN
       message_string = 'Too coarse angular discretization settings: ' //                          &
                        'raytrace_discrete_elevs < 4 and/or raytrace_discrete_elevs < 8'
       CALL message( 'radiation_check_parameters', 'RAD0026', 1, 2, 0, 6, 0 )
    ENDIF
!-- Warning message when small values are set
    IF ( raytrace_discrete_elevs < 9  .OR.  raytrace_discrete_azims < 18 )  THEN
       message_string = 'Relatively coarse angular discretization settings are set: ' //           &
                        'raytrace_discrete_elevs < 9 and/or raytrace_discrete_elevs < 18'
       CALL message( 'radiation_check_parameters', 'RAD0027', 0, 1, 0, 6, 0 )
    ENDIF
    !TODO: add a check that radiation_volumetric_flux requires radiation_interactions
#if defined( __tenstream )
!
!-- Error/warning messages for TenStream
    IF ( ts_icollapse /= -1_IINTEGERS  .AND.  ts_icollapse /= 1_IINTEGERS ) THEN
       message_string = 'Invalid ts_icollapse value. ts_icollapse should be either -1 (collapse'// &
                        ' atmosphere above the dynamic domain) or 1 (no collapse)'
       CALL message( 'radiation_check_parameters', 'RAD0028', 1, 2, 0, 6, 0 )
    ENDIF
#endif

 END SUBROUTINE radiation_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the radiation model and Radiative Transfer Model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_init

#if defined( __rrtmg )
    USE control_parameters,                                                                        &
        ONLY:  run_identifier
#endif
    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc

    IMPLICIT NONE

    INTEGER(iwp) ::  i          !< running index x-direction
    INTEGER(iwp) ::  is         !< running index for input surface elements
    INTEGER(iwp) ::  j          !< running index y-direction
    INTEGER(iwp) ::  k          !< running index z-direction
    INTEGER(iwp) ::  m          !< running index for surface elements
    INTEGER(iwp) ::  ntime = 0  !< number of available external radiation timesteps
#if defined( __rrtmg )
    INTEGER(iwp) ::  ind_type  !< running index for subgrid-surface tiles
#endif
    LOGICAL ::  radiation_input_root_domain  !< flag indicating the existence of a dynamic input file for the root domain


    IF ( debug_output )  CALL debug_message( 'radiation_init', 'start' )
    CALL cpu_log( log_point_s(93), 'radiation_init', 'start' )
!
!-- Activate radiation_interactions according to the existence of vertical surfaces and/or trees
!   or if biometeorology output is required for flat surfaces.
!-- The namelist parameter radiation_interactions_on can override this behavior (this check cannot
!-- be performed in check_parameters, because vertical_surfaces_exist is first set in
!-- init_surface_arrays).
    IF ( radiation_interactions_on )  THEN

       IF ( radiation_scheme == 'tenstream' ) THEN
          radiation_interactions_on = .FALSE.

          IF ( biometeorology ) THEN
             message_string = 'The tenstream radiation scheme does not support the ' //            &
                              'biometeorology, which requires RTM. Please use different' //        &
                              'radiation scheme, e.g. RRTMG.'
             CALL message( 'radiation_init', 'RAD0029', 1, 2, 0, 6, 0 )
          ENDIF

       ELSEIF ( ( vertical_surfaces_exist  .OR.  plant_canopy  .OR.  biometeorology )              &
                .AND.  .NOT. dcep )                                                                &
       THEN
          radiation_interactions = .TRUE.
          average_radiation      = .TRUE.
       ELSE
          radiation_interactions_on = .FALSE.   !< Reset namelist parameter: no interactions
                                                !< calculations necessary in case of flat surface
       ENDIF
    ELSEIF ( vertical_surfaces_exist  .OR.  plant_canopy  .OR.  biometeorology )  THEN
       message_string = 'radiation_interactions_on is set to .FALSE. although vertical ' //        &
                        'surfaces and/or trees or biometeorology exist is ON. The model will ' //  &
                        'run without RTM (no shadows, no radiation reflections)'
       CALL message( 'radiation_init', 'RAD0030', 0, 1, 0, 6, 0 )
    ENDIF
!
!-- Activiate average radiation in case of dcep average radiation is required
    IF ( dcep_average_radiation )  average_radiation = .TRUE.
!
!-- Warning message when cyclic boundary conditions are set
    IF ( radiation_interactions  .AND.  ( bc_lr_cyc  .OR.  bc_ns_cyc ) )  THEN
       message_string = 'The current raytracing algorithm in the Radiative Transfer Model does' // &
                        '&NOT support explicitly cyclic boundary conditions. Surface radiation' // &
                        '&fluxes near the boundaries should be evaluated in this view point.'
       CALL message( 'radiation_init', 'RAD0031', 0, 1, 0, 6, 0 )
    ENDIF
!
!-- In case of radiation_only runs, initialize default surface properties. The final albedo is
!-- calculated later like for LSM surfaces.
    IF ( radiation_only ) THEN
       IF ( .NOT. ALLOCATED( surf_def%albedo_type ) )  THEN
          ALLOCATE( surf_def%albedo_type(1:surf_def%ns,0:0) )
       ENDIF
       IF ( .NOT. ALLOCATED( surf_def%emissivity ) )  THEN
          ALLOCATE( surf_def%emissivity(1:surf_def%ns,0:0) )
       ENDIF
       surf_def%albedo_type = albedo_type
       surf_def%emissivity  = emissivity
    ENDIF
!
!-- Precalculate some time constants
    d_hours_day    = 1.0_wp / REAL( hours_per_day, KIND = wp )
    d_seconds_hour = 1.0_wp / seconds_per_hour

!
!-- If required, initialize radiation interactions between surfaces via sky-view factors. This must
!-- be done before radiation is initialized.
    IF ( radiation_interactions )  CALL radiation_interaction_init
!
!-- Allocate surface radiation arrays that are part of the restart mechanism.
!-- Arrays need to be allocated always in case of restarts, even if the number of surfaces
!-- on the subdomain is zero.
!-- Natural surfaces.
    IF ( .NOT. ALLOCATED( surf_lsm%rad_sw_in ) )  THEN
       ALLOCATE( surf_lsm%rad_sw_in(1:surf_lsm%ns) )
       IF ( surf_lsm%ns > 0 )  surf_lsm%rad_sw_in = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_lsm%rad_sw_out ) )  THEN
       ALLOCATE( surf_lsm%rad_sw_out(1:surf_lsm%ns) )
       IF ( surf_lsm%ns > 0 )  surf_lsm%rad_sw_out = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_lsm%rad_lw_in ) )  THEN
       ALLOCATE( surf_lsm%rad_lw_in(1:surf_lsm%ns) )
       IF ( surf_lsm%ns > 0 )  surf_lsm%rad_lw_in = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_lsm%rad_lw_out ) )  THEN
       ALLOCATE( surf_lsm%rad_lw_out(1:surf_lsm%ns) )
       IF ( surf_lsm%ns > 0 )  surf_lsm%rad_lw_out = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_lsm%rad_net ) )  THEN
       ALLOCATE( surf_lsm%rad_net(1:surf_lsm%ns) )
       IF ( surf_lsm%ns > 0 )  surf_lsm%rad_net = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_lsm%rad_lw_out_change_0 ) )  THEN
       ALLOCATE( surf_lsm%rad_lw_out_change_0(1:surf_lsm%ns) )
       IF ( surf_lsm%ns > 0 )  surf_lsm%rad_lw_out_change_0 = 0.0_wp
    ENDIF
!
!-- Urban surfaces.
    IF ( .NOT. ALLOCATED( surf_usm%rad_sw_in ) )  THEN
       ALLOCATE( surf_usm%rad_sw_in(1:surf_usm%ns) )
       IF ( surf_usm%ns > 0 )  surf_usm%rad_sw_in = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_usm%rad_sw_out ) )  THEN
       ALLOCATE( surf_usm%rad_sw_out(1:surf_usm%ns) )
       IF ( surf_usm%ns > 0 )  surf_usm%rad_sw_out = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_usm%rad_lw_in ) )  THEN
       ALLOCATE( surf_usm%rad_lw_in(1:surf_usm%ns) )
       IF ( surf_usm%ns > 0 )  surf_usm%rad_lw_in = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_usm%rad_lw_out ) )  THEN
       ALLOCATE( surf_usm%rad_lw_out(1:surf_usm%ns) )
       IF ( surf_usm%ns > 0 )  surf_usm%rad_lw_out = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_usm%rad_net ) )  THEN
       ALLOCATE( surf_usm%rad_net(1:surf_usm%ns) )
       IF ( surf_usm%ns > 0 )  surf_usm%rad_net = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( surf_usm%rad_lw_out_change_0 ) )  THEN
       ALLOCATE( surf_usm%rad_lw_out_change_0(1:surf_usm%ns) )
       IF ( surf_usm%ns > 0 )  surf_usm%rad_lw_out_change_0 = 0.0_wp
    ENDIF
!
!-- Allocate further arrays that are not part of the restart mechanism.
    IF ( surf_lsm%ns > 0 )  THEN
       ALLOCATE( surf_lsm%rad_sw_dir(1:surf_lsm%ns) )
       ALLOCATE( surf_lsm%rad_sw_dif(1:surf_lsm%ns) )
       ALLOCATE( surf_lsm%rad_sw_ref(1:surf_lsm%ns) )
       ALLOCATE( surf_lsm%rad_sw_res(1:surf_lsm%ns) )
       ALLOCATE( surf_lsm%rad_lw_dif(1:surf_lsm%ns) )
       ALLOCATE( surf_lsm%rad_lw_ref(1:surf_lsm%ns) )
       ALLOCATE( surf_lsm%rad_lw_res(1:surf_lsm%ns) )
       surf_lsm%rad_sw_dir = 0.0_wp
       surf_lsm%rad_sw_dif = 0.0_wp
       surf_lsm%rad_sw_ref = 0.0_wp
       surf_lsm%rad_sw_res = 0.0_wp
       surf_lsm%rad_lw_dif = 0.0_wp
       surf_lsm%rad_lw_ref = 0.0_wp
       surf_lsm%rad_lw_res = 0.0_wp
    ENDIF

    IF ( surf_usm%ns > 0 )  THEN
       ALLOCATE( surf_usm%rad_sw_dir(1:surf_usm%ns) )
       ALLOCATE( surf_usm%rad_sw_dif(1:surf_usm%ns) )
       ALLOCATE( surf_usm%rad_sw_ref(1:surf_usm%ns) )
       ALLOCATE( surf_usm%rad_sw_res(1:surf_usm%ns) )
       ALLOCATE( surf_usm%rad_lw_dif(1:surf_usm%ns) )
       ALLOCATE( surf_usm%rad_lw_ref(1:surf_usm%ns) )
       ALLOCATE( surf_usm%rad_lw_res(1:surf_usm%ns) )
       surf_usm%rad_sw_dir = 0.0_wp
       surf_usm%rad_sw_dif = 0.0_wp
       surf_usm%rad_sw_ref = 0.0_wp
       surf_usm%rad_sw_res = 0.0_wp
       surf_usm%rad_lw_dif = 0.0_wp
       surf_usm%rad_lw_ref = 0.0_wp
       surf_usm%rad_lw_res = 0.0_wp
    ENDIF
!
!-- Fix net radiation in case of radiation_scheme = 'constant'
    IF ( radiation_scheme == 'constant' )  THEN
!
!--    @Todo: weight with inclination angle
       IF ( ALLOCATED( surf_lsm%rad_net ) )  surf_lsm%rad_net = net_radiation
       IF ( ALLOCATED( surf_usm%rad_net ) )  surf_usm%rad_net = net_radiation
!
!-- Calculate orbital constants
    ELSE
       decl_1 = SIN( 23.45_wp * pi / 180.0_wp )
       decl_2 = 2.0_wp * pi / 365.0_wp
       decl_3 = decl_2 * 81.0_wp
       lat    = latitude * pi / 180.0_wp
       lon    = longitude * pi / 180.0_wp
    ENDIF
!
!-- Allocate direct and diffuse incoming radiation in case of DCEP model.
    IF ( dcep )  THEN
       IF ( .NOT. ALLOCATED( rad_sw_in_dir ) )  THEN
          ALLOCATE( rad_sw_in_dir(nysg:nyng,nxlg:nxrg) )
          rad_sw_in_dir  = 0.0_wp
       ENDIF
       IF ( .NOT. ALLOCATED( rad_sw_in_diff ) )  THEN
          ALLOCATE( rad_sw_in_diff(nysg:nyng,nxlg:nxrg) )
          rad_sw_in_diff = 0.0_wp
       ENDIF
       IF ( .NOT. ALLOCATED( rad_lw_in_diff ) )  THEN
          ALLOCATE( rad_lw_in_diff(nysg:nyng,nxlg:nxrg) )
          rad_lw_in_diff = 0.0_wp
       ENDIF
    ENDIF
!
!-- Allocate arrays based on the radiation scheme.
    IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'   .OR.              &
         radiation_scheme == 'external' )                                                          &
    THEN
!
!--    Allocate arrays for incoming/outgoing short/longwave radiation
       IF ( .NOT. ALLOCATED( rad_sw_in ) )  THEN
          ALLOCATE( rad_sw_in(0:0,nysg:nyng,nxlg:nxrg) )
          rad_sw_in = 0.0_wp
       ENDIF
       IF ( .NOT. ALLOCATED( rad_sw_out ) )  THEN
          ALLOCATE( rad_sw_out(0:0,nysg:nyng,nxlg:nxrg) )
          rad_sw_out = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_in ) )  THEN
          ALLOCATE( rad_lw_in(0:0,nysg:nyng,nxlg:nxrg) )
          rad_lw_in = 0.0_wp
       ENDIF
       IF ( .NOT. ALLOCATED( rad_lw_out ) )  THEN
          ALLOCATE( rad_lw_out(0:0,nysg:nyng,nxlg:nxrg) )
          rad_lw_out = 0.0_wp
       ENDIF

!
!--    Allocate average arrays for incoming/outgoing short/longwave radiation
       IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
          ALLOCATE( rad_sw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
       ENDIF
       IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
          ALLOCATE( rad_sw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
          ALLOCATE( rad_lw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
       ENDIF
       IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
          ALLOCATE( rad_lw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
       ENDIF
!
!--    Allocate arrays for broadband albedo, and level 1 initialization via namelist paramter,
!--    unless not already allocated.
       IF ( .NOT. ALLOCATED(surf_lsm%albedo) )  THEN
          ALLOCATE( surf_lsm%albedo(1:surf_lsm%ns,0:2) )
          surf_lsm%albedo = albedo
       ENDIF
       IF ( .NOT. ALLOCATED(surf_usm%albedo) )  THEN
          ALLOCATE( surf_usm%albedo(1:surf_usm%ns,0:2) )
          surf_usm%albedo = albedo
       ENDIF
!
!--    Level 2 initialization of broadband albedo via given albedo_type.
!--    Only if albedo_type is non-zero.
       DO  m = 1, surf_lsm%ns
          IF ( surf_lsm%albedo_type(m,ind_veg_wall) /= 0 )                                         &
             surf_lsm%albedo(m,ind_veg_wall) = albedo_pars(0,surf_lsm%albedo_type(m,ind_veg_wall))
          IF ( surf_lsm%albedo_type(m,ind_pav_green) /= 0 )                                        &
             surf_lsm%albedo(m,ind_pav_green) = albedo_pars(0,surf_lsm%albedo_type(m,ind_pav_green))
          IF ( surf_lsm%albedo_type(m,ind_wat_win) /= 0 )                                          &
             surf_lsm%albedo(m,ind_wat_win) = albedo_pars(0,surf_lsm%albedo_type(m,ind_wat_win))
       ENDDO
       DO  m = 1, surf_usm%ns
          IF ( surf_usm%albedo_type(m,ind_veg_wall) /= 0 )                                         &
             surf_usm%albedo(m,ind_veg_wall) = albedo_pars(0,surf_usm%albedo_type(m,ind_veg_wall))
          IF ( surf_usm%albedo_type(m,ind_pav_green) /= 0 )                                        &
             surf_usm%albedo(m,ind_pav_green) = albedo_pars(0,surf_usm%albedo_type(m,ind_pav_green))
          IF ( surf_usm%albedo_type(m,ind_wat_win) /= 0 )                                          &
             surf_usm%albedo(m,ind_wat_win) = albedo_pars(0,surf_usm%albedo_type(m,ind_wat_win))
       ENDDO

!
!--    Level 3 initialization at grid points where albedo type is zero.
!--    In this case, albedo is taken from file. In case of constant radiation or clear sky, only
!--    broadband albedo is given.
       IF ( albedo_pars_f%from_file )  THEN
          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m) + surf_lsm%ioff(m)
             j = surf_lsm%j(m) + surf_lsm%joff(m)
             IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                surf_lsm%albedo(m,ind_veg_wall)  = albedo_pars_f%pars_xy(0,j,i)
                surf_lsm%albedo(m,ind_pav_green) = albedo_pars_f%pars_xy(0,j,i)
                surf_lsm%albedo(m,ind_wat_win)   = albedo_pars_f%pars_xy(0,j,i)
             ENDIF
          ENDDO
          DO  m = 1, surf_usm%ns
             i = surf_usm%i(m) + surf_usm%ioff(m)
             j = surf_usm%j(m) + surf_usm%joff(m)
             IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )  THEN
                surf_usm%albedo(m,ind_veg_wall)  = albedo_pars_f%pars_xy(0,j,i)
                surf_usm%albedo(m,ind_pav_green) = albedo_pars_f%pars_xy(0,j,i)
                surf_usm%albedo(m,ind_wat_win)   = albedo_pars_f%pars_xy(0,j,i)
             ENDIF
          ENDDO
       ENDIF
!
!-- Read explicit albedo values from building surface pars. If present, they override all less
!-- specific albedo values and force an albedo_type to zero in order to take effect.
    IF ( building_surface_pars_f%from_file )  THEN
       DO  m = 1, surf_usm%ns
          i = surf_usm%i(m)
          j = surf_usm%j(m)
          k = surf_usm%k(m)
!
!-        Iterate over surfaces in column, check height and orientation
          DO  is = building_surface_pars_f%index_ji(1,j,i),                                        &
                   building_surface_pars_f%index_ji(2,j,i)
             IF ( ( ( ( building_surface_pars_f%coords(4,is) == -surf_usm%koff(m) )  .AND.         &
                      ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) ) )  .OR.                   &
                    ( building_surface_pars_f%coords(5,is) == -surf_usm%joff(m)    .AND.           &
                      building_surface_pars_f%coords(6,is) == -surf_usm%ioff(m) )  .AND.           &
                    .NOT. ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) ) )  .AND.              &
                  building_surface_pars_f%coords(1,is) == k )                                      &
             THEN

                IF ( building_surface_pars_f%pars(ind_s_alb_b_wall,is) /=                          &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%albedo(m,ind_veg_wall) =                                               &
                        building_surface_pars_f%pars(ind_s_alb_b_wall,is)
                   surf_usm%albedo_type(m,ind_veg_wall) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_b_win,is) /=                           &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%albedo(m,ind_wat_win) =                                                &
                        building_surface_pars_f%pars(ind_s_alb_b_win,is)
                   surf_usm%albedo_type(m,ind_wat_win) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_b_green,is) /=                         &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%albedo(m,ind_pav_green) =                                              &
                        building_surface_pars_f%pars(ind_s_alb_b_green,is)
                   surf_usm%albedo_type(m,ind_pav_green) = 0
                ENDIF
!
!--             Surface was found and processed
                EXIT
             ENDIF
          ENDDO
       ENDDO
    ENDIF
!
!-- Initialization actions for RRTMG
    ELSEIF ( radiation_scheme == 'rrtmg' )  THEN
#if defined( __rrtmg ) && defined( __netcdf )
!
!--    Allocate albedos for short/longwave radiation, horizontal surfaces for wall/green/window
!--    (USM), vegetation/pavement/water surfaces (LSM), or default surfaces
       IF ( radiation_only )  THEN
          ALLOCATE( surf_def%aldif(1:surf_def%ns,0:0) )
          ALLOCATE( surf_def%aldir(1:surf_def%ns,0:0) )
          ALLOCATE( surf_def%asdif(1:surf_def%ns,0:0) )
          ALLOCATE( surf_def%asdir(1:surf_def%ns,0:0) )
          ALLOCATE( surf_def%rrtm_aldif(1:surf_def%ns,0:0) )
          ALLOCATE( surf_def%rrtm_aldir(1:surf_def%ns,0:0) )
          ALLOCATE( surf_def%rrtm_asdif(1:surf_def%ns,0:0) )
          ALLOCATE( surf_def%rrtm_asdir(1:surf_def%ns,0:0) )
       ENDIF

       ALLOCATE( surf_lsm%aldif(1:surf_lsm%ns,0:2) )
       ALLOCATE( surf_lsm%aldir(1:surf_lsm%ns,0:2) )
       ALLOCATE( surf_lsm%asdif(1:surf_lsm%ns,0:2) )
       ALLOCATE( surf_lsm%asdir(1:surf_lsm%ns,0:2) )
       ALLOCATE( surf_lsm%rrtm_aldif(1:surf_lsm%ns,0:2) )
       ALLOCATE( surf_lsm%rrtm_aldir(1:surf_lsm%ns,0:2) )
       ALLOCATE( surf_lsm%rrtm_asdif(1:surf_lsm%ns,0:2) )
       ALLOCATE( surf_lsm%rrtm_asdir(1:surf_lsm%ns,0:2) )

       ALLOCATE( surf_usm%aldif(1:surf_usm%ns,0:2) )
       ALLOCATE( surf_usm%aldir(1:surf_usm%ns,0:2) )
       ALLOCATE( surf_usm%asdif(1:surf_usm%ns,0:2) )
       ALLOCATE( surf_usm%asdir(1:surf_usm%ns,0:2) )
       ALLOCATE( surf_usm%rrtm_aldif(1:surf_usm%ns,0:2) )
       ALLOCATE( surf_usm%rrtm_aldir(1:surf_usm%ns,0:2) )
       ALLOCATE( surf_usm%rrtm_asdif(1:surf_usm%ns,0:2) )
       ALLOCATE( surf_usm%rrtm_asdir(1:surf_usm%ns,0:2) )
!
!--    Allocate broadband albedo (temporary for the current radiation implementations)
       IF ( radiation_only )  THEN
          IF ( .NOT. ALLOCATED(surf_def%albedo) )  ALLOCATE( surf_def%albedo(1:surf_def%ns,0:0) )
       ENDIF
       IF ( .NOT. ALLOCATED(surf_lsm%albedo) )  ALLOCATE( surf_lsm%albedo(1:surf_lsm%ns,0:2) )
       IF ( .NOT. ALLOCATED(surf_usm%albedo) )  ALLOCATE( surf_usm%albedo(1:surf_usm%ns,0:2) )
!
!--    Level 1 initialization of spectral albedos via namelist paramters. Please note, in this case
!--    all surface tiles are initialized the same.
       IF ( radiation_only )  THEN
          IF ( surf_def%ns > 0 )  THEN
             surf_def%aldif  = albedo_lw_dif
             surf_def%aldir  = albedo_lw_dir
             surf_def%asdif  = albedo_sw_dif
             surf_def%asdir  = albedo_sw_dir
             surf_def%albedo = albedo_sw_dif
          ENDIF
       ENDIF
       IF ( surf_lsm%ns > 0 )  THEN
          surf_lsm%aldif  = albedo_lw_dif
          surf_lsm%aldir  = albedo_lw_dir
          surf_lsm%asdif  = albedo_sw_dif
          surf_lsm%asdir  = albedo_sw_dir
          surf_lsm%albedo = albedo_sw_dif
       ENDIF
       IF ( surf_usm%ns > 0 )  THEN
          surf_usm%aldif  = albedo_lw_dif
          surf_usm%aldir  = albedo_lw_dir
          surf_usm%asdif  = albedo_sw_dif
          surf_usm%asdir  = albedo_sw_dir
          surf_usm%albedo = albedo_sw_dif
       ENDIF

!
!--    Level 2 initialization of spectral albedos via albedo_type.
!--    Please note, for natural- and urban-type surfaces, a tile approach is applied so that the
!--    resulting albedo is calculated via the weighted average of respective surface fractions.
!--    For default surfaces, which become only relevant in radiation-only mode, no tile approach
!--    is employed.
       IF ( radiation_only )  THEN
          DO  m = 1, surf_def%ns
!
!--          Spectral albedos for default surfaces.
             IF ( surf_def%albedo_type(m,0) /= 0 )  THEN
                surf_def%aldif(m,0)  = albedo_pars(1,surf_def%albedo_type(m,0))
                surf_def%asdif(m,0)  = albedo_pars(2,surf_def%albedo_type(m,0))
                surf_def%aldir(m,0)  = albedo_pars(1,surf_def%albedo_type(m,0))
                surf_def%asdir(m,0)  = albedo_pars(2,surf_def%albedo_type(m,0))
                surf_def%albedo(m,0) = albedo_pars(0,surf_def%albedo_type(m,0))
             ENDIF
          ENDDO
       ENDIF

       DO  m = 1, surf_lsm%ns
!
!--       Spectral albedos for vegetation/pavement/water surfaces
          DO  ind_type = 0, 2
             IF ( surf_lsm%albedo_type(m,ind_type) /= 0 )  THEN
                surf_lsm%aldif(m,ind_type)  = albedo_pars(1,surf_lsm%albedo_type(m,ind_type))
                surf_lsm%asdif(m,ind_type)  = albedo_pars(2,surf_lsm%albedo_type(m,ind_type))
                surf_lsm%aldir(m,ind_type)  = albedo_pars(1,surf_lsm%albedo_type(m,ind_type))
                surf_lsm%asdir(m,ind_type)  = albedo_pars(2,surf_lsm%albedo_type(m,ind_type))
                surf_lsm%albedo(m,ind_type) = albedo_pars(0,surf_lsm%albedo_type(m,ind_type))
             ENDIF
          ENDDO

       ENDDO

       DO  m = 1, surf_usm%ns
!
!--       Spectral albedos for wall/green/window surfaces
          DO  ind_type = 0, 2
             IF ( surf_usm%albedo_type(m,ind_type) /= 0 )  THEN
                surf_usm%aldif(m,ind_type)  = albedo_pars(1,surf_usm%albedo_type(m,ind_type))
                surf_usm%asdif(m,ind_type)  = albedo_pars(2,surf_usm%albedo_type(m,ind_type))
                surf_usm%aldir(m,ind_type)  = albedo_pars(1,surf_usm%albedo_type(m,ind_type))
                surf_usm%asdir(m,ind_type)  = albedo_pars(2,surf_usm%albedo_type(m,ind_type))
                surf_usm%albedo(m,ind_type) = albedo_pars(0,surf_usm%albedo_type(m,ind_type))
             ENDIF
          ENDDO
       ENDDO
!
!--    Level 3 initialization at grid points where albedo type is zero.
!--    This case, spectral albedos are taken from file if available.
       IF ( albedo_pars_f%from_file )  THEN
          IF ( radiation_only )  THEN
             DO  m = 1, surf_def%ns
                i = surf_def%i(m) + surf_def%ioff(m)
                j = surf_def%j(m) + surf_def%joff(m)
!
!--             Spectral albedos for default surfaces.
                IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )                          &
                   surf_def%albedo(m,0) = albedo_pars_f%pars_xy(0,j,i)
                IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )                          &
                   surf_def%aldir(m,0) = albedo_pars_f%pars_xy(1,j,i)
                IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )                          &
                   surf_def%aldif(m,0) = albedo_pars_f%pars_xy(1,j,i)
                IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )                          &
                   surf_def%asdir(m,0) = albedo_pars_f%pars_xy(2,j,i)
                IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )                          &
                   surf_def%asdif(m,0) = albedo_pars_f%pars_xy(2,j,i)
             ENDDO
          ENDIF

          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m) + surf_lsm%ioff(m)
             j = surf_lsm%j(m) + surf_lsm%joff(m)
!
!--          Spectral albedos for vegetation/pavement/water surfaces
             DO  ind_type = 0, 2
                IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )                          &
                   surf_lsm%albedo(m,ind_type) = albedo_pars_f%pars_xy(0,j,i)
                IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )                          &
                   surf_lsm%aldir(m,ind_type) = albedo_pars_f%pars_xy(1,j,i)
                IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )                          &
                   surf_lsm%aldif(m,ind_type) = albedo_pars_f%pars_xy(1,j,i)
                IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )                          &
                   surf_lsm%asdir(m,ind_type) = albedo_pars_f%pars_xy(2,j,i)
                IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )                          &
                   surf_lsm%asdif(m,ind_type) = albedo_pars_f%pars_xy(2,j,i)
             ENDDO
          ENDDO

          DO  m = 1, surf_usm%ns
             i = surf_usm%i(m) + surf_usm%ioff(m)
             j = surf_usm%j(m) + surf_usm%joff(m)
!
!--          Broadband albedos for wall/green/window surfaces
             DO  ind_type = 0, 2
                IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )                          &
                   surf_usm%albedo(m,ind_type) = albedo_pars_f%pars_xy(0,j,i)
             ENDDO
!
!--          Spectral albedos especially for building wall surfaces
             IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )  THEN
                surf_usm%aldir(m,ind_veg_wall) = albedo_pars_f%pars_xy(1,j,i)
                surf_usm%aldif(m,ind_veg_wall) = albedo_pars_f%pars_xy(1,j,i)
             ENDIF
             IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )  THEN
                surf_usm%asdir(m,ind_veg_wall) = albedo_pars_f%pars_xy(2,j,i)
                surf_usm%asdif(m,ind_veg_wall) = albedo_pars_f%pars_xy(2,j,i)
             ENDIF
!
!--          Spectral albedos especially for building green surfaces
             IF ( albedo_pars_f%pars_xy(3,j,i) /= albedo_pars_f%fill )  THEN
                surf_usm%aldir(m,ind_pav_green) = albedo_pars_f%pars_xy(3,j,i)
                surf_usm%aldif(m,ind_pav_green) = albedo_pars_f%pars_xy(3,j,i)
             ENDIF
             IF ( albedo_pars_f%pars_xy(4,j,i) /= albedo_pars_f%fill )  THEN
                surf_usm%asdir(m,ind_pav_green) = albedo_pars_f%pars_xy(4,j,i)
                surf_usm%asdif(m,ind_pav_green) = albedo_pars_f%pars_xy(4,j,i)
             ENDIF
!
!--          Spectral albedos especially for building window surfaces
             IF ( albedo_pars_f%pars_xy(5,j,i) /= albedo_pars_f%fill )  THEN
                surf_usm%aldir(m,ind_wat_win) = albedo_pars_f%pars_xy(5,j,i)
                surf_usm%aldif(m,ind_wat_win) = albedo_pars_f%pars_xy(5,j,i)
             ENDIF
             IF ( albedo_pars_f%pars_xy(6,j,i) /= albedo_pars_f%fill )  THEN
                surf_usm%asdir(m,ind_wat_win) = albedo_pars_f%pars_xy(6,j,i)
                surf_usm%asdif(m,ind_wat_win) = albedo_pars_f%pars_xy(6,j,i)
             ENDIF
          ENDDO
       ENDIF
!
!--    Read explicit albedo values from building surface pars. If present, they override all less
!--    specific albedo values and force an albedo_type to zero in order to take effect.
       IF ( building_surface_pars_f%from_file )  THEN
          DO  m = 1, surf_usm%ns
             i = surf_usm%i(m)
             j = surf_usm%j(m)
             k = surf_usm%k(m)
!
!-           Iterate over surfaces in column, check height and orientation
             DO  is = building_surface_pars_f%index_ji(1,j,i),                                     &
                      building_surface_pars_f%index_ji(2,j,i)
                IF ( ( ( ( building_surface_pars_f%coords(4,is) == -surf_usm%koff(m) )  .AND.      &
                         ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) ) )  .OR.                &
                       ( building_surface_pars_f%coords(5,is) == -surf_usm%joff(m)    .AND.        &
                         building_surface_pars_f%coords(6,is) == -surf_usm%ioff(m) )  .AND.        &
                       .NOT. ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) ) )  .AND.           &
                     building_surface_pars_f%coords(1,is) == k )                                   &
                THEN

                   IF ( building_surface_pars_f%pars(ind_s_alb_b_wall,is) /=                       &
                        building_surface_pars_f%fill )                                             &
                   THEN
                      surf_usm%albedo(m,ind_veg_wall) =                                            &
                           building_surface_pars_f%pars(ind_s_alb_b_wall,is)
                      surf_usm%albedo_type(m,ind_veg_wall) = 0
                   ENDIF

                   IF ( building_surface_pars_f%pars(ind_s_alb_l_wall,is) /=                       &
                        building_surface_pars_f%fill )                                             &
                   THEN
                      surf_usm%aldir(m,ind_veg_wall) =                                             &
                           building_surface_pars_f%pars(ind_s_alb_l_wall,is)
                      surf_usm%aldif(m,ind_veg_wall) =                                             &
                           building_surface_pars_f%pars(ind_s_alb_l_wall,is)
                      surf_usm%albedo_type(m,ind_veg_wall) = 0
                   ENDIF

                   IF ( building_surface_pars_f%pars(ind_s_alb_s_wall,is) /=                       &
                        building_surface_pars_f%fill )                                             &
                   THEN
                      surf_usm%asdir(m,ind_veg_wall) =                                             &
                           building_surface_pars_f%pars(ind_s_alb_s_wall,is)
                      surf_usm%asdif(m,ind_veg_wall) =                                             &
                           building_surface_pars_f%pars(ind_s_alb_s_wall,is)
                      surf_usm%albedo_type(m,ind_veg_wall) = 0
                   ENDIF

                   IF ( building_surface_pars_f%pars(ind_s_alb_b_win,is) /=                        &
                        building_surface_pars_f%fill )                                             &
                   THEN
                      surf_usm%albedo(m,ind_wat_win) =                                             &
                           building_surface_pars_f%pars(ind_s_alb_b_win,is)
                      surf_usm%albedo_type(m,ind_wat_win) = 0
                   ENDIF

                   IF ( building_surface_pars_f%pars(ind_s_alb_l_win,is) /=                        &
                        building_surface_pars_f%fill )                                             &
                   THEN
                      surf_usm%aldir(m,ind_wat_win) =                                              &
                           building_surface_pars_f%pars(ind_s_alb_l_win,is)
                      surf_usm%aldif(m,ind_wat_win) =                                              &
                           building_surface_pars_f%pars(ind_s_alb_l_win,is)
                      surf_usm%albedo_type(m,ind_wat_win) = 0
                   ENDIF

                   IF ( building_surface_pars_f%pars(ind_s_alb_s_win,is) /=                        &
                        building_surface_pars_f%fill )                                             &
                   THEN
                      surf_usm%asdir(m,ind_wat_win) =                                              &
                           building_surface_pars_f%pars(ind_s_alb_s_win,is)
                      surf_usm%asdif(m,ind_wat_win) =                                              &
                           building_surface_pars_f%pars(ind_s_alb_s_win,is)
                      surf_usm%albedo_type(m,ind_wat_win) = 0
                   ENDIF

                   IF ( building_surface_pars_f%pars(ind_s_alb_b_green,is) /=                      &
                        building_surface_pars_f%fill )                                             &
                   THEN
                      surf_usm%albedo(m,ind_pav_green) =                                           &
                           building_surface_pars_f%pars(ind_s_alb_b_green,is)
                      surf_usm%albedo_type(m,ind_pav_green) = 0
                   ENDIF

                   IF ( building_surface_pars_f%pars(ind_s_alb_l_green,is) /=                      &
                        building_surface_pars_f%fill )                                             &
                   THEN
                      surf_usm%aldir(m,ind_pav_green) =                                            &
                           building_surface_pars_f%pars(ind_s_alb_l_green,is)
                      surf_usm%aldif(m,ind_pav_green) =                                            &
                           building_surface_pars_f%pars(ind_s_alb_l_green,is)
                      surf_usm%albedo_type(m,ind_pav_green) = 0
                   ENDIF

                   IF ( building_surface_pars_f%pars(ind_s_alb_s_green,is) /=                      &
                        building_surface_pars_f%fill )                                             &
                   THEN
                      surf_usm%asdir(m,ind_pav_green) =                                            &
                           building_surface_pars_f%pars(ind_s_alb_s_green,is)
                      surf_usm%asdif(m,ind_pav_green) =                                            &
                           building_surface_pars_f%pars(ind_s_alb_s_green,is)
                      surf_usm%albedo_type(m,ind_pav_green) = 0
                   ENDIF
!
!--                Surface was found and processed.
                   EXIT
                ENDIF
             ENDDO
          ENDDO
       ENDIF
!
!--    Calculate initial values of current (cosine of) the zenith angle and whether the sun is up
       CALL get_date_time( time_since_reference_point, day_of_year=day_of_year,                    &
                           second_of_day=second_of_day )
       CALL calc_zenith( day_of_year, second_of_day )
!
!--    Initialize spectral albedos.
       IF ( radiation_only  .AND.  surf_def%ns > 0 )  THEN
          surf_def%rrtm_aldir = surf_def%aldir
          surf_def%rrtm_asdir = surf_def%asdir
          surf_def%rrtm_aldif = surf_def%aldif
          surf_def%rrtm_asdif = surf_def%asdif
       ENDIF
       IF ( surf_lsm%ns > 0 )  THEN
          surf_lsm%rrtm_aldir = surf_lsm%aldir
          surf_lsm%rrtm_asdir = surf_lsm%asdir
          surf_lsm%rrtm_aldif = surf_lsm%aldif
          surf_lsm%rrtm_asdif = surf_lsm%asdif
       ENDIF
       IF ( surf_usm%ns > 0 )  THEN
          surf_usm%rrtm_aldir = surf_usm%aldir
          surf_usm%rrtm_asdir = surf_usm%asdir
          surf_usm%rrtm_aldif = surf_usm%aldif
          surf_usm%rrtm_asdif = surf_usm%asdif
       ENDIF
!
!--    Calculate sun-inclination dependent surface albedo.
       IF ( .NOT. constant_albedo )  THEN
          IF ( radiation_only )  CALL calc_albedo( surf_def )
          CALL calc_albedo( surf_lsm )
          CALL calc_albedo( surf_usm )
       ENDIF
!
!--    Allocate 3d arrays of radiative fluxes and heating rates
       IF ( .NOT. ALLOCATED( rad_sw_in ) )  THEN
          ALLOCATE( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_sw_in = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
          ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       IF ( .NOT. ALLOCATED( rad_sw_out ) )  THEN
          ALLOCATE( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_sw_out = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
          ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       IF ( .NOT. ALLOCATED( rad_sw_hr ) )  THEN
          ALLOCATE( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_sw_hr = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
          ALLOCATE( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_sw_hr_av = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_sw_cs_hr ) )  THEN
          ALLOCATE( rad_sw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_sw_cs_hr = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
          ALLOCATE( rad_sw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_sw_cs_hr_av = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_in ) )  THEN
          ALLOCATE( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_lw_in = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
          ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_out ) )  THEN
          ALLOCATE( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
         rad_lw_out = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
          ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_hr ) )  THEN
          ALLOCATE( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_lw_hr = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
          ALLOCATE( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_lw_hr_av = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_cs_hr ) )  THEN
          ALLOCATE( rad_lw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_lw_cs_hr = 0.0_wp
       ENDIF

       IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
          ALLOCATE( rad_lw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          rad_lw_cs_hr_av = 0.0_wp
       ENDIF

       ALLOCATE( rad_sw_cs_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( rad_sw_cs_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_sw_cs_in  = 0.0_wp
       rad_sw_cs_out = 0.0_wp

       ALLOCATE( rad_lw_cs_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       ALLOCATE( rad_lw_cs_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_lw_cs_in  = 0.0_wp
       rad_lw_cs_out = 0.0_wp

!
!--    Allocate 1-element array for surface temperature
!--    (RRTMG anticipates an array as passed argument).
       ALLOCATE( rrtm_tsfc(1) )
!
!--    Allocate surface emissivity.
!--    Values will be given directly before calling rrtm_lw.
       ALLOCATE( rrtm_emis(0:0,1:nbndlw+1) )

!
!--    Initialize RRTMG, before check if files are existent
       INQUIRE( FILE = 'RRTMG_LW', EXIST = lw_exists )
       IF ( .NOT. lw_exists )  THEN
          message_string = 'Input file RRTMG_LW for rrtmg model missing.& Please provide ' //      &
                           TRIM( run_identifier ) // '_rlw file in the INPUT directory.'
          CALL message( 'radiation_init', 'RAD0032', 1, 2, 0, 6, 0 )
       ENDIF
       INQUIRE( FILE = 'RRTMG_SW', EXIST = sw_exists )
       IF ( .NOT. sw_exists )  THEN
          message_string = 'Input file RRTMG_SW for rrtmg model missing.& Please provide ' //      &
                           TRIM( run_identifier ) // '_rsw file in the INPUT directory.'
          CALL message( 'radiation_init', 'RAD0033', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( lw_radiation )  CALL rrtmg_lw_ini ( c_p )
       IF ( sw_radiation )  CALL rrtmg_sw_ini ( c_p )

!
!--    Set input files for RRTMG
       INQUIRE( FILE = 'RAD_SND_DATA', EXIST = snd_exists )
       IF ( .NOT. snd_exists )  THEN
          rrtm_input_file = 'RRTMG_LW'
       ENDIF

!
!--    Read vertical layers for RRTMG from sounding data
!--    The routine provides nzt_rad, hyp_snd(1:nzt_rad), t_snd(nzt+2:nzt_rad), rrtm_play(1:nzt_rad),
!--    rrtm_plev(1_nzt_rad+1), rrtm_tlay(nzt+2:nzt_rad), rrtm_tlev(nzt+2:nzt_rad+1)
       CALL read_sounding_data

!
!--    Read trace gas profiles from file. This routine provides the rrtm_ arrays (1:nzt_rad+1)
       CALL read_trace_gas_data
#endif
!
!-- Initialization actions for TenStream
    ELSEIF ( radiation_scheme == 'tenstream' ) THEN

       CALL radiation_tenstream_init

    ENDIF
!
!-- Initializaion actions exclusively required for external radiation forcing
    IF ( radiation_scheme == 'external' )  THEN
!
!--    Open the radiation input file. Note, for child domain, a dynamic input file is often not
!--    provided. In order to not need to duplicate the dynamic input file just for the radiation
!--    input, take it from the dynamic file for the parent if not available for the child domain(s).
!--    In this case this is possible because radiation input should be the same for each model.
       INQUIRE( FILE = TRIM( input_file_dynamic ), EXIST = radiation_input_root_domain )

       IF ( .NOT. input_pids_dynamic  .AND.  .NOT.  radiation_input_root_domain )  THEN
          message_string = 'In case of external radiation forcing a dynamic input file for ' //    &
                           'the root domain is needed'
          CALL message( 'radiation_init', 'RAD0034', 1, 2, 0, 6, 0 )
       ENDIF
#if defined( __netcdf )
!
!--    Open dynamic input file for child domain if available, else, open dynamic input file for the
!--    root domain.
       IF ( input_pids_dynamic )  THEN
       CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ), pids_id )
       ELSEIF ( radiation_input_root_domain )  THEN
          CALL open_read_file( TRIM( input_file_dynamic ), pids_id )
       ENDIF

       CALL inquire_num_variables( pids_id, num_var_pids )
!
!--    Allocate memory to store variable names and read them
       ALLOCATE( vars_pids(1:num_var_pids) )
       CALL inquire_variable_names( pids_id, vars_pids )
!
!--    Input time dimension.
       IF ( check_existence( vars_pids, 'time_rad' ) )  THEN
          CALL get_dimension_length( pids_id, ntime, 'time_rad' )

          ALLOCATE( time_rad_f%var1d(0:ntime-1) )
!
!--       Read variable
          CALL get_variable( pids_id, 'time_rad', time_rad_f%var1d )

          time_rad_f%from_file = .TRUE.
       ENDIF
!
!--    Input shortwave downwelling.
       IF ( check_existence( vars_pids, 'rad_sw_in' ) )  THEN
!
!--       Get _FillValue attribute
          CALL get_attribute( pids_id, char_fill, rad_sw_in_f%fill, .FALSE., 'rad_sw_in' )
!
!--       Get level-of-detail
          CALL get_attribute( pids_id, char_lod, rad_sw_in_f%lod, .FALSE., 'rad_sw_in' )
!
!--       Level-of-detail 1 - radiation depends only on time_rad
          IF ( rad_sw_in_f%lod == 1 )  THEN
             ALLOCATE( rad_sw_in_f%var1d(0:ntime-1) )
             CALL get_variable( pids_id, 'rad_sw_in', rad_sw_in_f%var1d )
             rad_sw_in_f%from_file = .TRUE.
!
!--       Level-of-detail 2 - radiation depends on time_rad, y, x
          ELSEIF ( rad_sw_in_f%lod == 2 )  THEN
             ALLOCATE( rad_sw_in_f%var3d(0:ntime-1,nys:nyn,nxl:nxr) )

             CALL get_variable( pids_id, 'rad_sw_in', rad_sw_in_f%var3d, nxl, nxr, nys, nyn, 0,    &
                                ntime-1 )

             rad_sw_in_f%from_file = .TRUE.
          ELSE
             message_string = '"rad_sw_in" has no valid lod attribute'
             CALL message( 'radiation_init', 'RAD0035', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Input longwave downwelling.
       IF ( check_existence( vars_pids, 'rad_lw_in' ) )  THEN
!
!--       Get _FillValue attribute
          CALL get_attribute( pids_id, char_fill, rad_lw_in_f%fill, .FALSE., 'rad_lw_in' )
!
!--       Get level-of-detail
          CALL get_attribute( pids_id, char_lod, rad_lw_in_f%lod, .FALSE., 'rad_lw_in' )
!
!--       Level-of-detail 1 - radiation depends only on time_rad
          IF ( rad_lw_in_f%lod == 1 )  THEN
             ALLOCATE( rad_lw_in_f%var1d(0:ntime-1) )
             CALL get_variable( pids_id, 'rad_lw_in', rad_lw_in_f%var1d )
             rad_lw_in_f%from_file = .TRUE.
!
!--       Level-of-detail 2 - radiation depends on time_rad, y, x
          ELSEIF ( rad_lw_in_f%lod == 2 )  THEN
             ALLOCATE( rad_lw_in_f%var3d(0:ntime-1,nys:nyn,nxl:nxr) )

             CALL get_variable( pids_id, 'rad_lw_in', rad_lw_in_f%var3d, nxl, nxr, nys, nyn, 0,    &
                                ntime-1 )

             rad_lw_in_f%from_file = .TRUE.
          ELSE
             message_string = '"rad_lw_in" has no valid lod attribute'
             CALL message( 'radiation_init', 'RAD0035', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Input shortwave downwelling, diffuse part.
       IF ( check_existence( vars_pids, 'rad_sw_in_dif' ) )  THEN
!
!--       Read _FillValue attribute
          CALL get_attribute( pids_id, char_fill, rad_sw_in_dif_f%fill, .FALSE., 'rad_sw_in_dif' )
!
!--       Get level-of-detail
          CALL get_attribute( pids_id, char_lod, rad_sw_in_dif_f%lod, .FALSE., 'rad_sw_in_dif' )
!
!--       Level-of-detail 1 - radiation depends only on time_rad
          IF ( rad_sw_in_dif_f%lod == 1 )  THEN
             ALLOCATE( rad_sw_in_dif_f%var1d(0:ntime-1) )
             CALL get_variable( pids_id, 'rad_sw_in_dif', rad_sw_in_dif_f%var1d )
             rad_sw_in_dif_f%from_file = .TRUE.
!
!--       Level-of-detail 2 - radiation depends on time_rad, y, x
          ELSEIF ( rad_sw_in_dif_f%lod == 2 )  THEN
             ALLOCATE( rad_sw_in_dif_f%var3d(0:ntime-1,nys:nyn,nxl:nxr) )

             CALL get_variable( pids_id, 'rad_sw_in_dif', rad_sw_in_dif_f%var3d, nxl, nxr, nys,    &
                                nyn, 0, ntime-1 )

             rad_sw_in_dif_f%from_file = .TRUE.
          ELSE
             message_string = '"rad_sw_in_dif" has no valid lod attribute'
             CALL message( 'radiation_init', 'RAD0035', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Finally, close the input file and deallocate temporary arrays
       DEALLOCATE( vars_pids )

       CALL close_input_file( pids_id )
#endif
!
!--    Make some consistency checks.
       IF ( .NOT. rad_sw_in_f%from_file  .OR.  .NOT. rad_lw_in_f%from_file )  THEN
          message_string = 'In case of external radiation forcing both, rad_sw_in and ' //         &
                            'rad_lw_in are required.'
          CALL message( 'radiation_init', 'RAD0036', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( .NOT. time_rad_f%from_file )  THEN
          message_string = 'In case of external radiation forcing dimension time_rad is required.'
          CALL message( 'radiation_init', 'RAD0037', 1, 2, 0, 6, 0 )
       ENDIF

       CALL get_date_time( 0.0_wp, second_of_day=second_of_day )

       IF ( end_time - spinup_time > time_rad_f%var1d(ntime-1) )  THEN
          message_string = 'External radiation forcing does not cover the entire simulation time.'
          CALL message( 'radiation_init', 'RAD0039', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check for fill values in radiation
       IF ( ALLOCATED( rad_sw_in_f%var1d ) )  THEN
          IF ( ANY( rad_sw_in_f%var1d == rad_sw_in_f%fill ) )  THEN
             message_string = 'External radiation array "rad_sw_in" must not contain any ' //      &
                              'fill values.'
             CALL message( 'radiation_init', 'RAD0040', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( ALLOCATED( rad_lw_in_f%var1d ) )  THEN
          IF ( ANY( rad_lw_in_f%var1d == rad_lw_in_f%fill ) )  THEN
             message_string = 'External radiation array "rad_lw_in" must not contain any ' //      &
                              'fill values.'
             CALL message( 'radiation_init', 'RAD0041', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( ALLOCATED( rad_sw_in_dif_f%var1d ) )  THEN
          IF ( ANY( rad_sw_in_dif_f%var1d == rad_sw_in_dif_f%fill ) )  THEN
             message_string = 'External radiation array "rad_sw_in_dif" must not contain any ' //  &
                              'fill values.'
             CALL message( 'radiation_init', 'RAD0042', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( ALLOCATED( rad_sw_in_f%var3d ) )  THEN
          IF ( ANY( rad_sw_in_f%var3d == rad_sw_in_f%fill ) )  THEN
             message_string = 'External radiation array "rad_sw_in" must not contain any ' //      &
                              'fill values.'
             CALL message( 'radiation_init', 'RAD0040', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( ALLOCATED( rad_lw_in_f%var3d ) )  THEN
          IF ( ANY( rad_lw_in_f%var3d == rad_lw_in_f%fill ) )  THEN
             message_string = 'External radiation array "rad_lw_in" must not contain any ' //      &
                              'fill values.'
             CALL message( 'radiation_init', 'RAD0041', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( ALLOCATED( rad_sw_in_dif_f%var3d ) )  THEN
          IF ( ANY( rad_sw_in_dif_f%var3d == rad_sw_in_dif_f%fill ) )  THEN
             message_string = 'External radiation array "rad_sw_in_dif" must not contain any ' //  &
                              'fill values.'
             CALL message( 'radiation_init', 'RAD0042', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Currently, 2D external radiation input is not possible in combination with topography where
!--    average radiation is used.
       IF ( ( rad_lw_in_f%lod == 2  .OR.  rad_sw_in_f%lod == 2  .OR.                               &
              rad_sw_in_dif_f%lod == 2  )  .AND. average_radiation )  THEN
          message_string = 'External radiation with lod = 2 is currently not possible with ' //    &
                           'average_radiation = .T..'
             CALL message( 'radiation_init', 'RAD0043', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    All radiation input should have the same level of detail. The sum of lods divided by the
!--    number of available radiation arrays must be 1 (if all are lod = 1) or 2 (if all are lod = 2).
       IF ( REAL( MERGE( rad_lw_in_f%lod, 0, rad_lw_in_f%from_file ) +                             &
                  MERGE( rad_sw_in_f%lod, 0, rad_sw_in_f%from_file ) +                             &
                  MERGE( rad_sw_in_dif_f%lod, 0, rad_sw_in_dif_f%from_file ), KIND = wp ) /        &
                ( MERGE( 1.0_wp, 0.0_wp, rad_lw_in_f%from_file ) +                                 &
                  MERGE( 1.0_wp, 0.0_wp, rad_sw_in_f%from_file ) +                                 &
                  MERGE( 1.0_wp, 0.0_wp, rad_sw_in_dif_f%from_file ) ) /= 1.0_wp  .AND.            &
            REAL( MERGE( rad_lw_in_f%lod, 0, rad_lw_in_f%from_file ) +                             &
                  MERGE( rad_sw_in_f%lod, 0, rad_sw_in_f%from_file ) +                             &
                  MERGE( rad_sw_in_dif_f%lod, 0, rad_sw_in_dif_f%from_file ), KIND = wp ) /        &
                ( MERGE( 1.0_wp, 0.0_wp, rad_lw_in_f%from_file ) +                                 &
                  MERGE( 1.0_wp, 0.0_wp, rad_sw_in_f%from_file ) +                                 &
                  MERGE( 1.0_wp, 0.0_wp, rad_sw_in_dif_f%from_file ) )                             &
                  /= 2.0_wp )  THEN
          message_string = 'External radiation input should have the same lod.'
          CALL message( 'radiation_init', 'RAD0044', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF
!
!-- Perform user actions if required.
    CALL user_init_radiation

!
!-- Calculate radiative fluxes at model start. Not for restart runs.
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
       SELECT CASE ( TRIM( radiation_scheme ) )

          CASE ( 'rrtmg' )
             CALL radiation_rrtmg

          CASE ( 'clear-sky' )
             CALL radiation_clearsky

          CASE ( 'constant' )
             CALL radiation_constant

          CASE ( 'tenstream' )
             CALL radiation_tenstream

          CASE ( 'external' )
!
!--          During spinup apply clear-sky model
             IF ( time_since_reference_point < 0.0_wp )  THEN
                CALL radiation_clearsky
             ELSE
                CALL radiation_external
             ENDIF

          CASE DEFAULT

       END SELECT
    ENDIF

!
!-- If required, read or calculate and write out the SVF.
    IF ( radiation_interactions )  THEN

!
!--    Find all discretized apparent solar positions for radiation interaction.
       CALL radiation_presimulate_solar_pos
       ! TODO: this should be moved to the beginning of radiation_read_svf and the positions should
       ! be saved and loaded together with svf

       IF ( read_svf )  THEN
!
!--       Read sky-view factors and further required data from file
          CALL radiation_read_svf()
       ENDIF
!
!--    read_svf can be set .FALSE. in radiation_read_svf, therefore another IF and no ELSE
       IF ( .NOT. read_svf )  THEN
!
!--       Calculate svf and csf.
          CALL radiation_calc_svf()
       ENDIF

       IF ( write_svf )  THEN
!
!--       Write svf, csf svfsurf and csfsurf data to file.
          CALL radiation_write_svf()
       ENDIF
!
!--    Further specific initialization actions are required in case the optimized version of
!--    radiation_interaction is employed.
       IF ( loop_optimization == 'vector' )  THEN
          CALL radiation_interaction_vector_prepare
       ENDIF
!
!--    Adjust radiative fluxes.
!--    In case of urban and land surfaces, also call an initial interaction when it is not a
!--    restart run.
       IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          CALL radiation_interaction
       ELSE
!
!--       Just to calculate mrt quantities.
          CALL radiation_interaction( 'at_restart' )
       ENDIF

    ENDIF

    CALL cpu_log( log_point_s(93), 'radiation_init', 'stop' )
    IF ( debug_output )  CALL debug_message( 'radiation_init', 'end' )

 END SUBROUTINE radiation_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the arrays and the variables for the radiation scheme TenStream
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_tenstream_init
#if defined( __tenstream )

    USE control_parameters,                                                                        &
        ONLY:  topo_no_distinct

    IMPLICIT NONE

    INTEGER(iwp) ::  m         !< loop index
    INTEGER(iwp) ::  m_e       !< index for surface at east boundary
    INTEGER(iwp) ::  m_w       !< index for surface at west boundary
    INTEGER(iwp) ::  m_n       !< index for surface at north boundary
    INTEGER(iwp) ::  m_s       !< index for surface at south boundary
    INTEGER(iwp) ::  ncells    !< number of cells of buildings and/or orography
    INTEGER(iwp) ::  i         !< loop index
    INTEGER(iwp) ::  icell     !< loop index for building cell
    INTEGER(iwp) ::  ifacad    !< loop index for building facad
    INTEGER(iwp) ::  ind_type  !< index type
    INTEGER(iwp) ::  is        !< running index for input surface elements
    INTEGER(iwp) ::  j         !< loop index
    INTEGER(iwp) ::  k         !< loop index in z-direction
    INTEGER(iwp) ::  k_topo    !< topography top index
    INTEGER(iwp) ::  pc_k_top  !< k index of top plant canopy box

    INTEGER(IINTEGERS) ::  iface   !< loop index for face
    INTEGER(IINTEGERS) ::  nfaces  !< number of faces
    INTEGER(IINTEGERS) ::  ts_i    !< loop index for tenstream in x-direction
    INTEGER(IINTEGERS) ::  ts_j    !< loop index for tenstream in j-direction
    INTEGER(IINTEGERS) ::  ts_k    !< loop index for tenstream in z-direction
    INTEGER(IINTEGERS) ::  ts_m    !< loop index for tenstream for surface

    INTEGER(IINTEGERS) ::  ts_da_sizes(4)  !< size of tenstream surfaces

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: tmp_e   !< temporarily array to mark the surfaces at the east  border
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: tmp_eg  !< temporarily array to mark the surfaces at the east  border (received)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: tmp_n   !< temporarily array to mark the surfaces at the north border
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: tmp_ng  !< temporarily array to mark the surfaces at the north border (received)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: tmp_s   !< temporarily array to mark the surfaces at the south border
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: tmp_sg  !< temporarily array to mark the surfaces at the south border (received)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: tmp_w   !< temporarily array to mark the surfaces at the west  border
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE :: tmp_wg  !< temporarily array to mark the surfaces at the east  border (received)

    LOGICAL ::  building             !< flag indicating building grid point
    LOGICAL ::  terrain              !< flag indicating natural terrain grid point
    LOGICAL ::  ts_albedo_check      !< check if ts_aldif, ts_aldir, ts_asdif, and ts_asdir are the same
    LOGICAL ::  unresolved_building  !< flag indicating a grid point where actually a building is
                                     !< defined but not resolved by the vertical grid

    REAL(wp), PARAMETER ::  eps = 1.0E-10_wp  !< epsilon for value comparison

    REAL(IREALS) ::  sundir(3)  !< sun direction vector

    REAL(IREALS), PARAMETER ::  tree_albedo = 0.15  !< albedo of tree (assumed constant here for all bands)


    IF ( debug_output )  THEN
       WRITE( debug_string, * ) 'radiation_tenstream_init', time_since_reference_point
       CALL debug_message( debug_string, 'start' )
    ENDIF
!
!-- Allocate short/longwave and broadband albedos.
!-- Land
    ALLOCATE( surf_lsm%ts_albedo(1:surf_lsm%ns,0:2) )
    ALLOCATE( surf_lsm%ts_aldif(1:surf_lsm%ns,0:2)  )
    ALLOCATE( surf_lsm%ts_aldir(1:surf_lsm%ns,0:2)  )
    ALLOCATE( surf_lsm%ts_asdif(1:surf_lsm%ns,0:2)  )
    ALLOCATE( surf_lsm%ts_asdir(1:surf_lsm%ns,0:2)  )
!
!-- Urban
    ALLOCATE( surf_usm%ts_albedo(1:surf_usm%ns,0:2) )
    ALLOCATE( surf_usm%ts_aldif(1:surf_usm%ns,0:2)  )
    ALLOCATE( surf_usm%ts_aldir(1:surf_usm%ns,0:2)  )
    ALLOCATE( surf_usm%ts_asdif(1:surf_usm%ns,0:2)  )
    ALLOCATE( surf_usm%ts_asdir(1:surf_usm%ns,0:2)  )
!
!-- Allocate broadband albedo.
    IF ( .NOT. ALLOCATED( surf_lsm%albedo) )  ALLOCATE( surf_lsm%albedo(1:surf_lsm%ns,0:2) )
    IF ( .NOT. ALLOCATED( surf_usm%albedo) )  ALLOCATE( surf_usm%albedo(1:surf_usm%ns,0:2) )
!
!-- Level 1 initialization of spectral albedos via namelist paramters.
!-- Please note, this case all surface tiles are initialized the same.
    IF ( surf_lsm%ns > 0 )  THEN
       surf_lsm%ts_aldif = albedo_lw_dif
       surf_lsm%ts_aldir = albedo_lw_dir
       surf_lsm%ts_asdif = albedo_sw_dif
       surf_lsm%ts_asdir = albedo_sw_dir
       surf_lsm%albedo   = albedo_sw_dif
    ENDIF
    IF ( surf_usm%ns > 0 )  THEN
       surf_usm%ts_aldif = albedo_lw_dif
       surf_usm%ts_aldir = albedo_lw_dir
       surf_usm%ts_asdif = albedo_sw_dif
       surf_usm%ts_asdir = albedo_sw_dir
       surf_usm%albedo   = albedo_sw_dif
    ENDIF
!
!-- Level 2 initialization of spectral albedos via albedo_type.
!-- Please note, for natural- and urban-type surfaces, a tile approach
!-- is applied so that the resulting albedo is calculated via the weighted
!-- average of respective surface fractions.
!-- Spectral albedos for vegetation/pavement/water surfaces
    DO  m = 1, surf_lsm%ns
       DO  ind_type = 0, 2
          IF ( surf_lsm%albedo_type(m,ind_type) /= 0 )  THEN
             surf_lsm%ts_aldif(m,ind_type) = albedo_pars(1,surf_lsm%albedo_type(m,ind_type))
             surf_lsm%ts_asdif(m,ind_type) = albedo_pars(2,surf_lsm%albedo_type(m,ind_type))
             surf_lsm%ts_aldir(m,ind_type) = albedo_pars(1,surf_lsm%albedo_type(m,ind_type))
             surf_lsm%ts_asdir(m,ind_type) = albedo_pars(2,surf_lsm%albedo_type(m,ind_type))
             surf_lsm%albedo(m,ind_type)   = albedo_pars(0,surf_lsm%albedo_type(m,ind_type))
          ENDIF
       ENDDO
    ENDDO
!
!-- Spectral albedos for wall/green/window surfaces.
    DO  m = 1, surf_usm%ns
       DO  ind_type = 0, 2
          IF ( surf_usm%albedo_type(m,ind_type) /= 0 )  THEN
             surf_usm%ts_aldif(m,ind_type) = albedo_pars(1,surf_usm%albedo_type(m,ind_type))
             surf_usm%ts_asdif(m,ind_type) = albedo_pars(2,surf_usm%albedo_type(m,ind_type))
             surf_usm%ts_aldir(m,ind_type) = albedo_pars(1,surf_usm%albedo_type(m,ind_type))
             surf_usm%ts_asdir(m,ind_type) = albedo_pars(2,surf_usm%albedo_type(m,ind_type))
             surf_usm%albedo(m,ind_type)   = albedo_pars(0,surf_usm%albedo_type(m,ind_type))
          ENDIF
       ENDDO
    ENDDO
!
!-- Level 3 initialization at grid points where albedo type is zero.
!-- This case, spectral albedos are taken from file if available
    IF ( albedo_pars_f%from_file )  THEN
!
!--    Spectral albedos for vegetation/pavement/water surfaces
       DO  m = 1, surf_lsm%ns
          i = surf_lsm%i(m) + surf_lsm%ioff(m)
          j = surf_lsm%j(m) + surf_lsm%joff(m)
          DO  ind_type = 0, 2
             IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )                             &
                 surf_lsm%albedo(m,ind_type) = albedo_pars_f%pars_xy(0,j,i)
             IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )                             &
                 surf_lsm%aldir(m,ind_type) = albedo_pars_f%pars_xy(1,j,i)
             IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )                             &
                 surf_lsm%aldif(m,ind_type) = albedo_pars_f%pars_xy(1,j,i)
             IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )                             &
                 surf_lsm%asdir(m,ind_type) = albedo_pars_f%pars_xy(2,j,i)
             IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )                             &
                 surf_lsm%asdif(m,ind_type) = albedo_pars_f%pars_xy(2,j,i)
          ENDDO
       ENDDO
!
!--    Spectral albedos for wall/green/window surfaces
       DO  m = 1, surf_usm%ns
          i = surf_usm%i(m) + surf_usm%ioff(m)
          j = surf_usm%j(m) + surf_usm%joff(m)
!
!--       Broadband albedos for wall/green/window surfaces
          DO  ind_type = 0, 2
             IF ( albedo_pars_f%pars_xy(0,j,i) /= albedo_pars_f%fill )                             &
                surf_usm%albedo(m,ind_type) = albedo_pars_f%pars_xy(0,j,i)
          ENDDO
!
!--       Spectral albedos especially for building wall surfaces
          IF ( albedo_pars_f%pars_xy(1,j,i) /= albedo_pars_f%fill )  THEN
             surf_usm%aldir(m,ind_veg_wall) = albedo_pars_f%pars_xy(1,j,i)
             surf_usm%aldif(m,ind_veg_wall) = albedo_pars_f%pars_xy(1,j,i)
          ENDIF
          IF ( albedo_pars_f%pars_xy(2,j,i) /= albedo_pars_f%fill )  THEN
             surf_usm%asdir(m,ind_veg_wall) = albedo_pars_f%pars_xy(2,j,i)
             surf_usm%asdif(m,ind_veg_wall) = albedo_pars_f%pars_xy(2,j,i)
          ENDIF
!
!--       Spectral albedos especially for building green surfaces
          IF ( albedo_pars_f%pars_xy(3,j,i) /= albedo_pars_f%fill )  THEN
             surf_usm%aldir(m,ind_pav_green) = albedo_pars_f%pars_xy(3,j,i)
             surf_usm%aldif(m,ind_pav_green) = albedo_pars_f%pars_xy(3,j,i)
          ENDIF
          IF ( albedo_pars_f%pars_xy(4,j,i) /= albedo_pars_f%fill )  THEN
             surf_usm%asdir(m,ind_pav_green) = albedo_pars_f%pars_xy(4,j,i)
             surf_usm%asdif(m,ind_pav_green) = albedo_pars_f%pars_xy(4,j,i)
          ENDIF
!
!--       Spectral albedos especially for building window surfaces
          IF ( albedo_pars_f%pars_xy(5,j,i) /= albedo_pars_f%fill )  THEN
             surf_usm%aldir(m,ind_wat_win) = albedo_pars_f%pars_xy(5,j,i)
             surf_usm%aldif(m,ind_wat_win) = albedo_pars_f%pars_xy(5,j,i)
          ENDIF
          IF ( albedo_pars_f%pars_xy(6,j,i) /= albedo_pars_f%fill )  THEN
             surf_usm%asdir(m,ind_wat_win) = albedo_pars_f%pars_xy(6,j,i)
             surf_usm%asdif(m,ind_wat_win) = albedo_pars_f%pars_xy(6,j,i)
          ENDIF

       ENDDO
    ENDIF
!
!-- Read explicit albedo values from building surface pars. If present, they override all less
!-- specific albedo values and force an albedo_type to zero in order to take effect.
    IF ( building_surface_pars_f%from_file )  THEN
       DO  m = 1, surf_usm%ns
          i = surf_usm%i(m)
          j = surf_usm%j(m)
          k = surf_usm%k(m)
!
!-        Iterate over surfaces in column, check height and orientation
          DO  is = building_surface_pars_f%index_ji(1,j,i),                                        &
                   building_surface_pars_f%index_ji(2,j,i)
             IF ( ( ( ( building_surface_pars_f%coords(4,is) == -surf_usm%koff(m) )  .AND.         &
                      ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) ) )  .OR.                   &
                    ( building_surface_pars_f%coords(5,is) == -surf_usm%joff(m)    .AND.           &
                      building_surface_pars_f%coords(6,is) == -surf_usm%ioff(m) )  .AND.           &
                    .NOT. ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) ) )  .AND.              &
                  building_surface_pars_f%coords(1,is) == k )                                      &
             THEN

                IF ( building_surface_pars_f%pars(ind_s_alb_b_wall,is) /=                          &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%albedo(m,ind_veg_wall) =                                               &
                        building_surface_pars_f%pars(ind_s_alb_b_wall,is)
                   surf_usm%albedo_type(m,ind_veg_wall) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_l_wall,is) /=                          &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%aldir(m,ind_veg_wall) =                                                &
                        building_surface_pars_f%pars(ind_s_alb_l_wall,is)
                   surf_usm%aldif(m,ind_veg_wall) =                                                &
                        building_surface_pars_f%pars(ind_s_alb_l_wall,is)
                   surf_usm%albedo_type(m,ind_veg_wall) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_s_wall,is) /=                          &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%asdir(m,ind_veg_wall) =                                                &
                        building_surface_pars_f%pars(ind_s_alb_s_wall,is)
                   surf_usm%asdif(m,ind_veg_wall) =                                                &
                        building_surface_pars_f%pars(ind_s_alb_s_wall,is)
                   surf_usm%albedo_type(m,ind_veg_wall) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_b_win,is) /=                           &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%albedo(m,ind_wat_win) =                                                &
                        building_surface_pars_f%pars(ind_s_alb_b_win,is)
                   surf_usm%albedo_type(m,ind_wat_win) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_l_win,is) /=                           &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%aldir(m,ind_wat_win) =                                                 &
                        building_surface_pars_f%pars(ind_s_alb_l_win,is)
                   surf_usm%aldif(m,ind_wat_win) =                                                 &
                        building_surface_pars_f%pars(ind_s_alb_l_win,is)
                   surf_usm%albedo_type(m,ind_wat_win) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_s_win,is) /=                           &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%asdir(m,ind_wat_win) =                                                 &
                        building_surface_pars_f%pars(ind_s_alb_s_win,is)
                   surf_usm%asdif(m,ind_wat_win) =                                                 &
                        building_surface_pars_f%pars(ind_s_alb_s_win,is)
                   surf_usm%albedo_type(m,ind_wat_win) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_b_green,is) /=                         &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%albedo(m,ind_pav_green) =                                              &
                        building_surface_pars_f%pars(ind_s_alb_b_green,is)
                   surf_usm%albedo_type(m,ind_pav_green) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_l_green,is) /=                         &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%aldir(m,ind_pav_green) =                                               &
                        building_surface_pars_f%pars(ind_s_alb_l_green,is)
                   surf_usm%aldif(m,ind_pav_green) =                                               &
                        building_surface_pars_f%pars(ind_s_alb_l_green,is)
                   surf_usm%albedo_type(m,ind_pav_green) = 0
                ENDIF

                IF ( building_surface_pars_f%pars(ind_s_alb_s_green,is) /=                         &
                     building_surface_pars_f%fill )                                                &
                THEN
                   surf_usm%asdir(m,ind_pav_green) =                                               &
                        building_surface_pars_f%pars(ind_s_alb_s_green,is)
                   surf_usm%asdif(m,ind_pav_green) =                                               &
                        building_surface_pars_f%pars(ind_s_alb_s_green,is)
                   surf_usm%albedo_type(m,ind_pav_green) = 0
                ENDIF
!
!--             Surface was found and processed
                EXIT
             ENDIF
          ENDDO
       ENDDO
    ENDIF
!
!-- Check if ts_aldif, ts_aldir, ts_asdif, and ts_asdir are the same
    ts_albedo_check = .FALSE.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             IF ( ABS( SUM( surf_lsm%frac(m,:) * surf_lsm%ts_asdir(m,:) ) -                        &
                       SUM( surf_lsm%frac(m,:) * surf_lsm%ts_asdif(m,:) )                          &
                     )  >  eps  .OR.                                                               &
                  ABS( SUM( surf_lsm%frac(m,:) * surf_lsm%ts_asdir(m,:) ) -                        &
                       SUM( surf_lsm%frac(m,:) * surf_lsm%ts_aldir(m,:) )                          &
                     )  >  eps  .OR.                                                               &
                  ABS( SUM( surf_lsm%frac(m,:) * surf_lsm%ts_asdir(m,:) ) -                        &
                       SUM( surf_lsm%frac(m,:) * surf_lsm%ts_aldif(m,:) )                          &
                     )  >  eps )                                                                   &
              THEN
                 ts_albedo_check = .TRUE.
              ENDIF
          ENDDO

          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             IF ( ABS( SUM( surf_usm%frac(m,:) * surf_usm%ts_asdir(m,:) ) -                        &
                       SUM( surf_usm%frac(m,:) * surf_usm%ts_asdif(m,:) )                          &
                     )  >  eps  .OR.                                                               &
                  ABS( SUM( surf_usm%frac(m,:) * surf_usm%ts_asdir(m,:) ) -                        &
                       SUM( surf_usm%frac(m,:) * surf_usm%ts_aldir(m,:) )                          &
                     )  >  eps  .OR.                                                               &
                  ABS( SUM( surf_usm%frac(m,:) * surf_usm%ts_asdir(m,:) ) -                        &
                       SUM( surf_usm%frac(m,:) * surf_usm%ts_aldif(m,:) )                          &
                     )  >  eps )                                                                   &
             THEN
                ts_albedo_check = .TRUE.
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF ( ts_albedo_check  .AND.  .NOT. use_broadband_albedo ) THEN
       WRITE( message_string, * ) 'Surface albedo for diffuse/direct long-/shortwave are ' //      &
                                   'different. To continue using the broadband albedo please ' //  &
                                   'set use_broadband_albedo to TRUE.'
       CALL message( 'radiation_tenstream_init', 'RAD0045', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( use_broadband_albedo ) THEN
!
!--    Use the broadband albedo for the TS albedo
       IF ( surf_lsm%ns > 0 )  surf_lsm%ts_albedo = surf_lsm%albedo
       IF ( surf_usm%ns > 0 )  surf_usm%ts_albedo = surf_usm%albedo
    ELSE
!
!--    ts_aldif, ts_aldir, ts_asdif, and ts_asdir are the same. We copy ts_asdir to the TS albedo
       IF ( surf_lsm%ns > 0 )  surf_lsm%ts_albedo = surf_lsm%ts_asdir
       IF ( surf_usm%ns > 0 )  surf_usm%ts_albedo = surf_usm%ts_asdir

    ENDIF
!
!-- Calculate initial values of current (cosine of) the zenith angle and whether the sun is up
    CALL get_date_time( time_since_reference_point, day_of_year=day_of_year,                       &
                        second_of_day=second_of_day )
    CALL calc_zenith( day_of_year, second_of_day )
!
!-- Calculate initial surface albedo for different surfaces
    IF ( .NOT. constant_albedo )  THEN
       message_string = 'variable albedo is not implemented for tenstream'
       CALL message( 'radiation_tenstream_init', 'RAD0046',1, 2, 0, 6, 0 )
    ENDIF
!
!-- Allocate 3d arrays of radiative fluxes and heating rates
    IF ( .NOT. ALLOCATED ( rad_sw_in ) )  THEN
       ALLOCATE ( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_sw_in = 0.0_wp
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_sw_in_av ) )  THEN
       ALLOCATE ( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_sw_out ) )  THEN
       ALLOCATE ( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_sw_out = 0.0_wp
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_sw_out_av ) )  THEN
       ALLOCATE ( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_sw_hr ) )  THEN
       ALLOCATE ( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_sw_hr = 0.0_wp
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_sw_hr_av ) )  THEN
       ALLOCATE ( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_sw_hr_av = 0.0_wp
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_lw_in ) )  THEN
       ALLOCATE ( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_lw_in = 0.0_wp
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_lw_in_av ) )  THEN
       ALLOCATE ( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_lw_out ) )  THEN
       ALLOCATE ( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_lw_out = 0.0_wp
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_lw_out_av ) )  THEN
       ALLOCATE ( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_lw_hr ) )  THEN
       ALLOCATE ( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_lw_hr = 0.0_wp
    ENDIF

    IF ( .NOT. ALLOCATED ( rad_lw_hr_av ) )  THEN
       ALLOCATE ( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       rad_lw_hr_av = 0.0_wp
    ENDIF

!
!-- TenStream related arrays
    ts_nranksx = npex
    ts_nranksy = npey

    ALLOCATE( ts_nxproc(ts_nranksx) )
    ts_nxproc = nnx
    ALLOCATE( ts_nyproc(ts_nranksy) )
    ts_nyproc = nny

    ALLOCATE( ts_skin_temperature(nxl:nxr,nys:nyn) )
    ALLOCATE( ts_solar_albedo_2d(nxl:nxr,nys:nyn) )
    ALLOCATE( ts_thermal_albedo_2d(nxl:nxr,nys:nyn) )

    ALLOCATE( ts_plev(nzb:nzt+1,nxl:nxr,nys:nyn) )
    ALLOCATE( ts_tlev(nzb:nzt+1,nxl:nxr,nys:nyn) )
    ALLOCATE( ts_play(nzb+1:nzt+1) )
    ALLOCATE( ts_tlay(nzb+1:nzt+1,nxl:nxr,nys:nyn) )
    ALLOCATE( ts_lwc(nzb+1:nzt+1,nxl:nxr,nys:nyn) )
    ALLOCATE( ts_reliq(nzb+1:nzt+1,nxl:nxr,nys:nyn) )
    ALLOCATE( ts_h2ovmr(nzb+1:nzt+1,nxl:nxr,nys:nyn) )

    ts_dx = REAL( dx, IREALS )
    ts_dy = REAL( dy, IREALS )

    IF ( albedo_lw_dif == 9999999.9_wp  .OR.  albedo_sw_dif == 9999999.9_wp )  THEN
       nfaces = 0
       albedo_th  = 0.0_IREALS
       albedo_sol = 0.0_IREALS

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                IF ( surf_lsm%upward(m) )  THEN
                   nfaces = nfaces + 1
                   albedo_th  = albedo_th + REAL( 1.0_wp -                                         &
                                   SUM( surf_lsm%frac(m,:) * surf_lsm%emissivity(m,:) ), IREALS )
                   albedo_sol = albedo_sol + REAL(                                                 &
                                    SUM( surf_lsm%frac(m,:) * surf_lsm%ts_albedo(m,:) ), IREALS )
                ENDIF
             ENDDO
             DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                IF ( surf_usm%upward(m) )  THEN
                   nfaces = nfaces + 1
                   albedo_th  = albedo_th  + REAL( 1.0_wp -                                        &
                                   SUM( surf_usm%frac(m,:) * surf_usm%emissivity(m,:) ), IREALS )
                   albedo_sol = albedo_sol + REAL(                                                 &
                                    SUM( surf_usm%frac(m,:) * surf_usm%ts_albedo(m,:) ), IREALS )
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       albedo_th  = albedo_th  / REAL( nfaces, KIND=IREALS )
       albedo_sol = albedo_sol / REAL( nfaces, KIND=IREALS )
    ELSE
       albedo_th  = REAL( albedo_lw_dif, IREALS )
       albedo_sol = REAL( albedo_sw_dif, IREALS )
    ENDIF

!
!-- Change from PETSC (C) domain splitting to MPI(Fortran) domain, i.e. row major vs col-major
!-- ordering of subdomains in x/y.
    CALL REORDER_MPI_COMM( comm2d, ts_nranksx, ts_nranksy, ts_comm )
    CALL INIT_MPI_DATA_PARAMETERS( ts_comm )

!
!-- Allocate the tenstream solver, by default it uses the provided solver (tenstream_solver) but
!-- can be changed via runtime option with e.g. `-solver 8_16`.
    CALL ALLOCATE_PPRTS_SOLVER_FROM_COMMANDLINE( ts_solver, tenstream_solver, ierr )
    CALL CHKERR( ierr )

!
!-- Setup building structure:
!
!-- Pressure level
    ts_plev(nzb,:,:) = REAL( surface_pressure, IREALS )
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt+1
             ts_plev(k,i,j) = REAL( barometric_formula( zw(k), pt_surface * exner(k),              &
                                                        surface_pressure ),                        &
                                    IREALS )
             ts_tlay(k,i,j) = REAL( pt(k,j,i) * exner(k), IREALS )
          ENDDO
       ENDDO
    ENDDO
!
!-- Here we use upward surfaces only.
    DO  i = nxl, nxr
       DO  j = nys, nyn

          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             ts_tlev(nzb,i,j) = MERGE( REAL( surf_lsm%pt_surface(m) * exner(nzb), IREALS ),        &
                                       ts_tlev(nzb,i,j),                                           &
                                       surf_lsm%upward(m) )
          ENDDO

          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             ts_tlev(nzb,i,j) = MERGE( REAL( surf_usm%pt_surface(m) * exner(nzb), IREALS ),        &
                                       ts_tlev(nzb,i,j),                                           &
                                       surf_usm%upward(m) )
          ENDDO

          DO k = nzb+1, nzt
             ts_tlev(k,i,j) = 0.5_IREALS * ( ts_tlay(k,i,j) + ts_tlay(k+1,i,j) )
          ENDDO
          ts_tlev(nzt+1,i,j) = ts_tlev(nzt,i,j)

       ENDDO
    ENDDO

    pplev(1:SIZE(ts_plev,1)  ,1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_plev
    ptlay(1:SIZE(ts_plev,1)-1,1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_tlay
    ptlev(1:SIZE(ts_plev,1)  ,1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_tlev
!
!-- Setup tenstream
    CALL SETUP_TENSTR_ATM( ts_comm, .False., ts_atm_filename, pplev, ptlev, ts_atm )

    CALL radiation_calc_sundir( sundir )

    CALL PPRTS_RRTMG( ts_comm, ts_solver, ts_atm, INT( nnx, IINTEGERS ), INT( nny, IINTEGERS ),    &
                      ts_dx, ts_dy, sundir, albedo_th, albedo_sol, lw_radiation, sw_radiation,     &
                      ts_edir, ts_edn, ts_eup, ts_abso, icollapse = ts_icollapse,                  &
                      nxproc = ts_nxproc, nyproc = ts_nyproc, lonly_initialize = .TRUE. )

    ts_xm = ts_solver%c_one%xm
    ts_ym = ts_solver%c_one%ym
    ts_zm = ts_solver%c_one%zm
!
!-- Setup building structure
!
!-- Count the number of obstacle cells
    ncells = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
!
!--          Check if current gridpoint belongs to obstacle (buildings or terrain)
             IF ( .NOT. BTEST( topo_flags(k,j,i), 0 ) )  ncells = ncells + 1
          ENDDO
       ENDDO
    ENDDO

    nfaces = ncells * 6_IINTEGERS

    IF ( debug_output )  THEN
       WRITE( debug_string, * )  'TenStream: Number building/orography boxes = ', ncells
       CALL debug_message( debug_string, 'info' )
    ENDIF

    ts_da_sizes = [6_IINTEGERS, ts_zm, ts_xm, ts_ym]

    CALL INIT_BUILDINGS( buildings_solar, ts_da_sizes, nfaces, ierr )
    CALL CHKERR( ierr )
!
!-- Count the number of facades and orography surfaces:
!-- facades and orography surfaces = all surfaces - surfaces at ground level.
    nfacad = 0
!
!-- Land surfaces
    DO  m = 1, surf_lsm%ns
       nfacad = nfacad + MERGE( 1, 0, surf_lsm%k(m) + surf_lsm%koff(m) > 0 )
    ENDDO
!
!-- Urban surfaces
    DO  m = 1, surf_usm%ns
       nfacad = nfacad + MERGE( 1, 0, surf_usm%k(m) + surf_usm%koff(m) > 0 )
    ENDDO
!
!-- Prepare to exchange boundaries:
!-- Due to the surface module strategy to define vertical surfaces, surfaces located at the
!-- boundaries need to be exchanged because TenStream define the surfaces of each building/orography
!-- box locally. Here we can exchange (a) the whole boundary faces (not implemented) or (b) exchange
!-- only the existing facads (implemented below).
!
!-- Count the number of faces located at boundaries. Generally there are two types of these faces:
!-- (1) Faces belong to current PE, which should be sent to the neighbour PE and removed from nfacad
!-- (2) Faces belong to neighbour PEs, which should be received by this PE
!
!-- East
    nfacad_east  = 0
    nfacad_eastg = 0
!
!-- West
    nfacad_west  = 0
    nfacad_westg = 0
!
!-- North
    nfacad_north  = 0
    nfacad_northg = 0
!
!-- South
    nfacad_south  = 0
    nfacad_southg = 0

    DO  j = nys, nyn
!
!--    East
!--    Surfaces belongs to this PE but defined for TS at east neighbour PE
       DO  k = nzb+1,nzt
          IF ( BTEST( topo_flags(k,j,nxr), 0 )  .AND.  .NOT. BTEST( topo_flags(k,j,nxr+1), 0 ) )   &
          THEN
             nfacad_east = nfacad_east + 1
          ENDIF
       ENDDO
!
!--    Surfaces belongs to east neighbour PE but defined for TS at this PE
       DO  k = nzb+1,nzt
          IF ( .NOT. BTEST( topo_flags(k,j,nxr), 0 )  .AND.  BTEST( topo_flags(k,j,nxr+1), 0 ) )   &
          THEN
             nfacad_eastg = nfacad_eastg + 1
          ENDIF
       ENDDO
!
!--    West
!--    Surfaces belongs to this PE but defined for TS at west neighbour PE
       DO  k = nzb+1,nzt
          IF ( BTEST( topo_flags(k,j,nxl), 0 )  .AND.  .NOT. BTEST( topo_flags(k,j,nxl-1), 0 ) )   &
          THEN
             nfacad_west = nfacad_west + 1
          ENDIF
       ENDDO
!
!--    Surfaces belongs to west neighbour PE but defined for TS at this PE
       DO k = nzb+1,nzt
          IF ( .NOT. BTEST( topo_flags(k,j,nxl), 0 )  .AND.  BTEST( topo_flags(k,j,nxl-1), 0 ) )   &
          THEN
             nfacad_westg = nfacad_westg + 1
          ENDIF
       ENDDO
    ENDDO

    DO  i = nxl, nxr
!
!--    North
!--    Surfaces belongs to this PE but defined for TS at north neighbour PE
       DO  k = nzb+1,nzt
          IF ( BTEST( topo_flags(k,nyn,i), 0 )  .AND.  .NOT. BTEST( topo_flags(k,nyn+1,i), 0 ) )   &
          THEN
             nfacad_north = nfacad_north + 1
          ENDIF
       ENDDO
!
!--    Surfaces belongs to north neighbour PE but defined for TS at this PE
       DO  k = nzb+1,nzt
          IF ( .NOT. BTEST( topo_flags(k,nyn,i), 0 )  .AND.  BTEST( topo_flags(k,nyn+1,i), 0 ) )   &
          THEN
             nfacad_northg = nfacad_northg + 1
          ENDIF
       ENDDO
!
!--    South
!--    Surfaces belongs to this PE but defined for TS at south neighbour PE
       DO  k = nzb+1,nzt
          IF ( BTEST( topo_flags(k,nys,i), 0 )  .AND.  .NOT. BTEST( topo_flags(k,nys-1,i), 0 ) )   &
          THEN
             nfacad_south = nfacad_south + 1
          ENDIF
       ENDDO
!
!--    Surfaces belongs to south neighbour PE but defined for TS at this PE
       DO  k = nzb+1,nzt
          IF ( .NOT. BTEST( topo_flags(k,nys,i), 0 )  .AND.  BTEST( topo_flags(k,nys-1,i), 0 ) )   &
          THEN
             nfacad_southg = nfacad_southg + 1
          ENDIF
       ENDDO
    ENDDO

!
!-- Report the number of exchanged surfaces
    IF ( debug_output  ) THEN
       WRITE( 9, * ) 'TenStream: Number of exchange surfaces:'
       WRITE( 9, * ) '           East=> send: ',nfacad_east,' receive: ',nfacad_eastg
       WRITE( 9, * ) '           West=> send: ',nfacad_west,' receive: ',nfacad_westg
       WRITE( 9, * ) '           North=> send: ',nfacad_north,' receive: ',nfacad_northg
       WRITE( 9, * ) '           South=> send: ',nfacad_south,' receive: ',nfacad_southg
       FLUSH( 9 )
    ENDIF

!
!-- Send receive signals
    request_count = 0
    IF ( nfacad_eastg > 0 )  THEN
       ALLOCATE( surf_ids_eg(2,nfacad_eastg) )
       surf_ids_eg = -99999
!
!--    Add m,i,j,k
       ALLOCATE( tmp_eg(4,nfacad_eastg) )
       request_count = request_count + 1
       CALL MPI_IRECV( tmp_eg, 4*nfacad_eastg, MPI_INTEGER, pright, tag_w, MPI_COMM_WORLD,         &
                       requests(request_count), ierr )
    ENDIF

    IF ( nfacad_westg > 0 )  THEN
       ALLOCATE( surf_ids_wg(2,nfacad_westg) )
       surf_ids_wg = -99999
!
!--    Add m,i,j,k
       ALLOCATE( tmp_wg(4,nfacad_westg) )
       request_count = request_count + 1
       CALL MPI_IRECV( tmp_wg, 4*nfacad_westg, MPI_INTEGER, pleft, tag_e, MPI_COMM_WORLD,          &
                       requests(request_count), ierr )
    ENDIF

    IF ( nfacad_northg > 0 )  THEN
       ALLOCATE( surf_ids_ng(2,nfacad_northg) )
       surf_ids_ng = -99999
!
!--    Add m,i,j,k
       ALLOCATE( tmp_ng(4,nfacad_northg) )
       request_count = request_count + 1
       CALL MPI_IRECV( tmp_ng, 4*nfacad_northg, MPI_INTEGER, pnorth, tag_s, MPI_COMM_WORLD,        &
                       requests(request_count), ierr )
    ENDIF

    IF ( nfacad_southg > 0 )  THEN
       ALLOCATE( surf_ids_sg(2,nfacad_southg) )
       surf_ids_sg = -99999
!
!--    Add m,i,j,k
       ALLOCATE( tmp_sg(4,nfacad_southg) )
       request_count = request_count + 1
       CALL MPI_IRECV( tmp_sg, 4*nfacad_southg, MPI_INTEGER, psouth, tag_n, MPI_COMM_WORLD,        &
                       requests(request_count), ierr )
    ENDIF

!
!-- Construct the boarder surface arrays
!-- East
    IF ( nfacad_east > 0 )  THEN
       CALL fill_surf_ids_tmp_array( nxr, nxr, nys, nyn, "e", nfacad_east, surf_ids_e, tmp_e )
       request_count = request_count + 1
       CALL MPI_ISEND( tmp_e, 4*nfacad_east, MPI_INTEGER, pright, tag_e, MPI_COMM_WORLD,           &
                       requests(request_count), ierr )
    ENDIF
!
!-- West
    IF ( nfacad_west > 0 )  THEN
       CALL fill_surf_ids_tmp_array( nxl, nxl, nys, nyn, "w", nfacad_west, surf_ids_w, tmp_w )
       request_count = request_count + 1
       CALL MPI_ISEND( tmp_w, 4*nfacad_west, MPI_INTEGER, pleft, tag_w, MPI_COMM_WORLD,            &
                       requests(request_count), ierr )
    ENDIF
!
!-- North
    IF ( nfacad_north > 0 )  THEN
       CALL fill_surf_ids_tmp_array( nxl, nxr, nyn, nyn, "n", nfacad_north, surf_ids_n, tmp_n )
       request_count = request_count + 1
       CALL MPI_ISEND( tmp_n, 4*nfacad_north, MPI_INTEGER, pnorth, tag_n, MPI_COMM_WORLD,          &
                       requests(request_count), ierr )
    ENDIF
!
!-- South
    IF ( nfacad_south > 0 )  THEN
       CALL fill_surf_ids_tmp_array( nxl, nxr, nys, nys, "s", nfacad_south, surf_ids_s, tmp_s )
       request_count = request_count + 1
       CALL MPI_ISEND( tmp_s, 4*nfacad_south, MPI_INTEGER, psouth, tag_s, MPI_COMM_WORLD,         &
                       requests(request_count), ierr )
    ENDIF

    CALL MPI_WAITALL( request_count, requests(1:request_count), wait_stat, ierr )

!
!-- Deduct the boarder surfaces from nfacad.
    nfacad = nfacad - nfacad_east - nfacad_west - nfacad_north - nfacad_south

    IF ( debug_output )  THEN
       WRITE( debug_string, * ) 'Number of facads/surf. for TenStream = ', nfacad
       CALL debug_message( debug_string, 'info' )
    ENDIF

    ALLOCATE( surf_ids(3,nfacad) )
    surf_ids = -99999

    icell  = 0
    ifacad = 0
    m_e    = 0
    m_w    = 0
    m_n    = 0
    m_s    = 0

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO k = nzb+1, topo_top_ind(j,i,0)
!
!--          Check if a cell is obstacle or terrain
             IF ( .NOT. BTEST( topo_flags(k,j,i), 0 ) )  THEN

                terrain  = BTEST( topo_flags(k,j,i), 5 )  .OR.  topo_no_distinct
                building = BTEST( topo_flags(k,j,i), 6 )  .OR.  topo_no_distinct

                unresolved_building = BTEST( topo_flags(k,j,i), 5 )  .AND.                         &
                                      BTEST( topo_flags(k,j,i), 6 )

                icell = icell + 1

                ts_k = ts_solver%c_one%zm - ( k - ( nzb + 1 ) )
                ts_i = 1 + ( i - nxl )
                ts_j = 1 + ( j - nys )

                DO  iface = 1, 6

                   ts_m = INT( ( icell - 1 ) * 6 + iface, IINTEGERS )

                   SELECT CASE ( iface )
!
!--                   Upward surface
                      CASE ( PPRTS_TOP_FACE )
                         IF ( BTEST( topo_flags(k+1,j,i), 0 ) )  THEN
!
!--                         Check if land surface
                            IF ( terrain  .AND.  .NOT. unresolved_building )  THEN
                               CALL find_surface( surf_lsm, i, j, k, iup_l, .FALSE. )
!
!--                         Check if urban surface
                            ELSEIF ( building )  THEN
                               CALL find_surface( surf_usm, i, j, k, iup_u, .FALSE. )
                            ELSE
                               WRITE( message_string, * ) 'undefined upward surface, not urban/land'
                               CALL message( 'radiation_tenstream_init', 'RAD0047', 1, 2, 0, 6, 0 )
                            ENDIF
                         ENDIF
!
!--                   Downward surface - attention, not fully implemented yet.
                      CASE ( PPRTS_BOT_FACE )
                         IF ( BTEST( topo_flags(k-1,j,i), 0 ) )  THEN
                            CALL find_surface( surf_lsm, i, j, k, idown_l, .FALSE. )
                            CALL find_surface( surf_usm, i, j, k, idown_u, .FALSE. )
                         ENDIF
!
!--                   Eastward-facing surface
                      CASE ( PPRTS_RIGHT_FACE )
                         IF ( BTEST( topo_flags(k,j,i+1), 0 )  )  THEN
                            IF ( i+1 <= nxr )  THEN
!
!--                            Check if land surface
                               IF ( terrain  .AND.  .NOT. unresolved_building )  THEN
                                  CALL find_surface( surf_lsm, i+1, j, k, ieast_l, .TRUE. )
!
!--                            Check if urban surface
                               ELSEIF ( building )  THEN
                                  CALL find_surface( surf_usm, i+1, j, k, ieast_u, .TRUE. )
                               ELSE
                                  WRITE( message_string, * ) 'undefined eastward surface,',        &
                                                             'not urban/land'
                                  CALL message( 'radiation_tenstream_init', 'RAD0047', 1, 2, 0, 6, &
                                                0 )
                               ENDIF
                            ELSE
!
!--                            Westward facing surface at the eastern boarder but belongs to the
!--                            east PE
                               CALL find_surface_g( surf_ids_eg, tmp_eg, i+1, j, k, m_e, ts_m )
                            ENDIF
                         ENDIF
!
!--                   Westward-facing surface
                      CASE ( PPRTS_LEFT_FACE )
                         IF ( BTEST( topo_flags(k,j,i-1), 0 )  )  THEN
                            IF ( i-1 >= nxl )  THEN
!
!--                            Check if land surface
                               IF ( terrain  .AND.  .NOT. unresolved_building )  THEN
                                  CALL find_surface( surf_lsm, i-1, j, k, iwest_l, .TRUE. )
!
!--                            Check if urban surface
                               ELSEIF ( building )  THEN
                                  CALL find_surface( surf_usm, i-1, j, k, iwest_u, .TRUE. )
                               ELSE
                                  WRITE( message_string, * ) 'undefined westward surface,',        &
                                                             'not urban/land'
                                  CALL message( 'radiation_tenstream_init', 'RAD0047', 1, 2, 0, 6, &
                                                0 )
                               ENDIF
                            ELSE
!
!--                            Eastward facing surface at the western boarder but belongs to the
!--                            west PE
                               CALL find_surface_g( surf_ids_wg, tmp_wg, i-1, j, k, m_w, ts_m )
                            ENDIF
                         ENDIF
!
!--                   Northward-facing surface
                      CASE ( PPRTS_FRONT_FACE )
                         IF ( BTEST( topo_flags(k,j+1,i), 0 )  )  THEN
                            IF ( j+1 <= nyn )  THEN
!
!--                            Check if land surface
                               IF ( terrain  .AND.  .NOT. unresolved_building )  THEN
                                  CALL find_surface( surf_lsm, i, j+1, k, inorth_l, .TRUE. )
!
!--                            Check if urban surface
                               ELSEIF ( building )  THEN
                                  CALL find_surface( surf_usm, i, j+1, k, inorth_u, .TRUE. )
                               ELSE
                                  WRITE( message_string, * ) 'undefined northward surface, ',      &
                                                             'not urban/land'
                                  CALL message( 'radiation_tenstream_init', 'RAD0047', 1, 2, 0, 6, &
                                                0 )
                               ENDIF
                            ELSE
!
!--                            Southward facing surface at the northern boarder but belongs to the
!--                            south PE
                               CALL find_surface_g( surf_ids_ng, tmp_ng, i, j+1, k, m_n, ts_m )
                            ENDIF
                         ENDIF
!
!--                   Southward-facing surface
                      CASE ( PPRTS_REAR_FACE )
                         IF ( BTEST( topo_flags(k,j-1,i), 0 )  )  THEN
                            IF ( j-1 >= nys )  THEN
!
!--                            Check if land surface
                               IF ( terrain  .AND.  .NOT. unresolved_building )  THEN
                                  CALL find_surface( surf_lsm, i, j-1, k, isouth_l, .TRUE. )
!
!--                            Check if urban surface
                               ELSEIF ( building )  THEN
                                  CALL find_surface( surf_usm, i, j-1, k, isouth_u, .TRUE. )
                               ELSE
                                  WRITE( message_string, * ) 'undefined southward surface, ',      &
                                                             'not urban/land'
                                  CALL message( 'radiation_tenstream_init', 'RAD0047', 1, 2, 0, 6, &
                                                0 )
                               ENDIF
                            ELSE
!
!--                            Northward facing surface at the southern boarder but belongs to the
!--                            north PE
                               CALL find_surface_g( surf_ids_sg, tmp_sg, i, j-1, k, m_s, ts_m )
                            ENDIF
                         ENDIF

                   END SELECT

                   buildings_solar%iface(ts_m) =                                                   &
                   FACEIDX_BY_CELL_PLUS_OFFSET(buildings_solar%da_offsets, ts_k, ts_i, ts_j, iface )
                ENDDO

             ENDIF ! if obstacle or terrain

          ENDDO
       ENDDO
    ENDDO

    CALL CHECK_BUILDINGS_CONSISTENCY( buildings_solar, ts_solver%c_one%zm, ts_solver%c_one%xm,     &
                                      ts_solver%c_one%ym, ierr )
    CALL CHKERR( ierr )

    CALL CLONE_BUILDINGS( buildings_solar, buildings_thermal, l_copy_data=.True., ierr=ierr )
    CALL CHKERR( ierr )

    ALLOCATE( buildings_thermal%temp(nfaces) )

    buildings_thermal%temp(:)   = 300.0_IREALS
    buildings_solar%albedo(:)   =   0.1_IREALS
    buildings_thermal%albedo(:) =   0.1_IREALS

!
!-- Allocate surface arrays
    IF ( nfacad_eastg  > 0 )  ALLOCATE( surf_eg(3,nfacad_eastg ) )
    IF ( nfacad_east   > 0 )  ALLOCATE( surf_e (3,nfacad_east  ) )
    IF ( nfacad_westg  > 0 )  ALLOCATE( surf_wg(3,nfacad_westg ) )
    IF ( nfacad_west   > 0 )  ALLOCATE( surf_w (3,nfacad_west  ) )
    IF ( nfacad_northg > 0 )  ALLOCATE( surf_ng(3,nfacad_northg) )
    IF ( nfacad_north  > 0 )  ALLOCATE( surf_n (3,nfacad_north ) )
    IF ( nfacad_southg > 0 )  ALLOCATE( surf_sg(3,nfacad_southg) )
    IF ( nfacad_south  > 0 )  ALLOCATE( surf_s (3,nfacad_south ) )

!
!-- Plant canopy
!-- Consider the resolved vegetation using LAD field and the tree albedo. The optical properties tau and w0
!-- are calculated by LAD and albedo, respectively.
    IF ( plant_canopy )  THEN
!
!--    Find the highest PC box in the sub-domain
       pc_k_top = 0
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzt+1, 1, -1
                IF ( lad_s(k,j,i) > 0.0_wp )  THEN
                   pc_k_top = MAX( k, pc_k_top )
                   EXIT
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       IF ( pc_k_top > 0 )  THEN
          ALLOCATE( tree_tau_solar(pc_k_top, ts_xm, ts_ym, ts_ngptsw) )
          ALLOCATE( tree_w0_solar(pc_k_top, ts_xm, ts_ym, ts_ngptsw) )
          ALLOCATE( tree_tau_thermal(pc_k_top, ts_xm, ts_ym, ts_ngptlw) )
          tree_tau_thermal = 0.0_IREALS
       ENDIF
!
!--    Fill the optical properties arrays with lad and albedo values.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzt+1, 1, -1

                IF ( lad_s(k,j,i) > 0.0_wp )  THEN

                   ts_i = 1 + ( i - nxl )
                   ts_j = 1 + ( j - nys )
                   ts_k = k
                   k_topo = topo_top_ind(j,i,0)

                   tree_w0_solar(ts_k,ts_i,ts_j,:)    = tree_albedo
                   tree_tau_solar(ts_k,ts_i,ts_j,:)   = REAL( lad_s(k-k_topo,j,i) / dz(1), IREALS )
                   tree_tau_thermal(ts_k,ts_i,ts_j,:) = REAL( lad_s(k-k_topo,j,i) / dz(1), IREALS )

                ENDIF

             ENDDO
          ENDDO
       ENDDO

    ENDIF

    IF ( debug_output )  CALL debug_message( 'radiation_tenstream_init', 'end' )

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> find_surface fills the surface id array for all surfaces (surface type, surface PALM id,
!> surface TenStream id
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE find_surface( surface, ii, jj, kk, isurf, vertical_surface  )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  ii     !< index for x-direction
    INTEGER(iwp), INTENT(IN) ::  isurf  !< surface type
    INTEGER(iwp), INTENT(IN) ::  jj     !< index for y-direction
    INTEGER(iwp)             ::  kk     !< index for z-direction
    INTEGER(iwp)             ::  ksurf  !< local k index of the surface

    LOGICAL             ::  surface_found     !< surface is successfully found
    LOGICAL, INTENT(IN) ::  vertical_surface  !< flag indicating vertical surfaces

    TYPE( surf_type )   ::  surface           !< respective surface type


    surface_found = .FALSE.

    DO  m = surface%start_index(jj,ii), surface%end_index(jj,ii)

       ksurf = surface%k(m) + surface%koff(m)

       IF ( kk == ksurf  .AND.                                                                     &
            ( surface%upward(m)     .AND.  ( isurf == iup_l     .OR.  isurf == iup_u    ) )  .OR.  &
            ( surface%downward(m)   .AND.  ( isurf == idown_l   .OR.  isurf == idown_u  ) )  .OR.  &
            ( surface%eastward(m)   .AND.  ( isurf == ieast_l   .OR.  isurf == ieast_u  ) )  .OR.  &
            ( surface%westward(m)   .AND.  ( isurf == iwest_l   .OR.  isurf == iwest_u  ) )  .OR.  &
            ( surface%southward(m)  .AND.  ( isurf == isouth_l  .OR.  isurf == isouth_u ) )  .OR.  &
            ( surface%northward(m)  .AND.  ( isurf == inorth_l  .OR.  isurf == inorth_u ) ) )      &
       THEN
          ifacad = ifacad + 1

          surf_ids(1,ifacad) = isurf
          surf_ids(2,ifacad) = m
          surf_ids(3,ifacad) = INT( ts_m )

          surface_found = .TRUE.
          IF ( vertical_surface )  EXIT
       ENDIF

    ENDDO

    IF ( .NOT. surface_found )  THEN
       WRITE( message_string, '(A,3I7,A,I7,A,I7,A)' ) 'Finding surface ids at i,j,k= ', ii, jj, kk,&
                                                      'for surface type', isurf, ' failed.'
       CALL message( 'radiation_tenstream_init', 'RAD0048', 1, 2, 0, 6, 0 )
    ENDIF

  END SUBROUTINE find_surface


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>
!--------------------------------------------------------------------------------------------------!
  SUBROUTINE find_surface_g( surface_ids, surface, i, j, k, ifacad, ts_m )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::       i       !<
    INTEGER(iwp), INTENT(INOUT) ::    ifacad  !<
    INTEGER(iwp) ::                   isurf   !< surface index
    INTEGER(iwp) ::                   i_surf  !< surface i-index
    INTEGER(iwp), INTENT(IN) ::       j       !<
    INTEGER(iwp) ::                   j_surf  !< surface j-index
    INTEGER(iwp), INTENT(IN) ::       k       !<
    INTEGER(iwp) ::                   k_surf  !< surface k-idex
    INTEGER(iwp) ::                   m_surf  !< surface ID

    INTEGER(IINTEGERS), INTENT(IN) :: ts_m  !<

    INTEGER(iwp), DIMENSION(:,:), INTENT(IN) ::  surface  !< respective surface array

    INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) ::  surface_ids

    LOGICAL ::  surface_found !< surface is successfully found


    surface_found = .FALSE.

    DO  isurf = 1, SIZE( surface, 2 )

       m_surf = surface(1,isurf)
       i_surf = surface(2,isurf)
       j_surf = surface(3,isurf)
       k_surf = surface(4,isurf)

       IF ( i == i_surf  .AND.  j == j_surf  .AND.  k == k_surf )  THEN
          ifacad = ifacad + 1
          surface_ids(1,ifacad) = m_surf
          surface_ids(2,ifacad) = INT( ts_m )
          surface_found = .TRUE.
          EXIT
       ENDIF

    ENDDO

    IF ( .NOT. surface_found ) THEN
       WRITE( message_string, '(A,3I7,A,I7,A,I7,A)' ) 'Finding boarder surface id at i,j,k= ',     &
                                                      i, j, k, ' failed.'
       CALL message( 'radiation_tenstream_init', 'RAD0049', 1, 2, 0, 6, 0 )
    ENDIF

  END SUBROUTINE find_surface_g


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> fill_surf_id_tmp_array fills two arrayes: surf_ids_X and tmp_X which are used to mark the
!> surfaces at the borders
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fill_surf_ids_tmp_array( nx1, nx2, ny1, ny2, dir, nfacad_b, surf_ids_b, tmp )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  dir   !< char to indicate the facing of the surface
    INTEGER(iwp) ::              m         !< surface ID
    INTEGER(iwp) ::              m_b       !< tmp array index
    INTEGER(iwp), INTENT(IN) ::  nfacad_b
    INTEGER(iwp), INTENT(IN) ::  nx1
    INTEGER(iwp), INTENT(IN) ::  nx2
    INTEGER(iwp), INTENT(IN) ::  ny1
    INTEGER(iwp), INTENT(IN) ::  ny2

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) ::  surf_ids_b
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) ::  tmp


    ALLOCATE( surf_ids_b(2,nfacad_b) )
    ALLOCATE( tmp(4,nfacad_b) )

    m_b = 0
    DO  i = nx1, nx2
       DO  j = ny1, ny2
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             IF ( ( surf_lsm%eastward(m)   .AND.  TRIM( dir ) == "w" )  .OR.                       &
                  ( surf_lsm%westward(m)   .AND.  TRIM( dir ) == "e" )  .OR.                       &
                  ( surf_lsm%southward(m)  .AND.  TRIM( dir ) == "n" )  .OR.                       &
                  ( surf_lsm%northward(m)  .AND.  TRIM( dir ) == "s" )  )                          &
             THEN
                m_b               = m_b + 1
                tmp(1,m_b)        = m
                tmp(2,m_b)        = surf_lsm%i(m)
                tmp(3,m_b)        = surf_lsm%j(m)
                tmp(4,m_b)        = surf_lsm%k(m)
                surf_ids_b(1,m_b) = m
                surf_ids_b(2,m_b) = 1
             ENDIF
          ENDDO
          DO m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             IF ( ( surf_usm%eastward(m)   .AND.  TRIM( dir ) == "w" )  .OR.                       &
                  ( surf_usm%westward(m)   .AND.  TRIM( dir ) == "e" )  .OR.                       &
                  ( surf_usm%southward(m)  .AND.  TRIM( dir ) == "n" )  .OR.                       &
                  ( surf_usm%northward(m)  .AND.  TRIM( dir ) == "s" )  )                          &
             THEN
                m_b               = m_b + 1
                tmp(1,m_b)        = m
                tmp(2,m_b)        = surf_usm%i(m)
                tmp(3,m_b)        = surf_usm%j(m)
                tmp(4,m_b)        = surf_usm%k(m)
                surf_ids_b(1,m_b) = m
                surf_ids_b(2,m_b) = 2
             ENDIF
          ENDDO
       ENDDO
    ENDDO
!
!-- The indices i and j need to be re-adjusted due to either the cyclic BC or in case of single
!-- domain decomposition. At right border i is reset to -1, at left border i is reset to nx+1,
!-- at north border j is reset to -1, and at south border j is reset to ny+1.
    IF ( ( npex == 1  .OR.  right_border_pe )  .AND.  TRIM( dir ) == "e" )  tmp(2,:) = -1
    IF ( ( npex == 1  .OR.  left_border_pe  )  .AND.  TRIM( dir ) == "w" )  tmp(2,:) = nx+1
    IF ( ( npey == 1  .OR.  north_border_pe )  .AND.  TRIM( dir ) == "n" )  tmp(3,:) = -1
    IF ( ( npey == 1  .OR.  south_border_pe )  .AND.  TRIM( dir ) == "s" )  tmp(3,:) = ny+1

 END SUBROUTINE fill_surf_ids_tmp_array
#endif

 END SUBROUTINE radiation_tenstream_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Implementation of coupling PALM/TenStream to enable 3D radiation calculations in the domain. The
!> TenStream model is a 3D parallel radiative transfer model for atmospheric heating rates for use
!> in cloud resolving modelsdescribed in details in:
!>     http://dx.doi.org/10.1016/j.jqsrt.2015.05.003
!>     http://dx.doi.org/10.5194/gmd-9-1413-2016
!>
!> The coupling of TenStream to PALM is realized by performing the following steps each time
!> TenStream is called:
!>    1) Provide vertical profiles of actual temperature and water vapor mixing ratio for the PALM
!>       column at full and half levels of the TenStream grid
!>    2) Provide vertical profiles of the in-cloud liquid water path and the effective droplet
!>       radius for each grid volume for each PALM column in case of clouds
!>    3) Provide surface data and surface temperature for the land surfaces as well as the building
!>       surfaces to TenStream
!>    4) Call TenStream shortwave and longwave radiation routines
!>    5) map back the radiative fluxes for each PALM grid and all surfaces
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_tenstream
#if defined( __tenstream )

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  rd_d_cp

    USE indices,                                                                                   &
        ONLY:  nnz

    USE particle_attributes,                                                                       &
        ONLY:  grid_particles,                                                                     &
               number_of_particles,                                                                &
               particles,                                                                          &
               prt_count

    IMPLICIT NONE

    INTEGER(iwp) ::  i           !< loop index in x-direction
    INTEGER(iwp) ::  ifacad      !< loop index for facad
    INTEGER(iwp) ::  isurf_type  !< surface type
    INTEGER(iwp) ::  j           !< loop index in y-direction
    INTEGER(iwp) ::  k           !< loop index in z-direction
    INTEGER(iwp) ::  k_topo      !< topography top index
    INTEGER(iwp) ::  m           !< surface ID
    INTEGER(iwp) ::  n           !< loop index for particles
    INTEGER(iwp) ::  nlay        !< number of layers
    INTEGER(iwp) ::  nlev        !< number of levels
    INTEGER(iwp) ::  ts_i        !< loop index in x-direction for tenstream
    INTEGER(iwp) ::  ts_j        !< loop index in y-direction for tenstream
    INTEGER(iwp) ::  ts_k        !< loop index in z-direction for tenstream
    INTEGER(iwp) ::  ts_m        !< surface index for tenstream

    INTEGER(iwp), DIMENSION(0:10), PARAMETER ::  facad_type = (/ 1, &   !<  0: upward land surf
                                                                 3, &   !<  1: downward default surf (not implemented)
                                                                 1, &   !<  2: eastward land  surface
                                                                 1, &   !<  3: westward land  surface
                                                                 1, &   !<  4: northward land surface
                                                                 1, &   !<  5: southward land  surface
                                                                 2, &   !<  6: upward urban surface
                                                                 2, &   !<  7: eastward urban surface
                                                                 2, &   !<  8: westward urban surface
                                                                 2, &   !<  9: northward urban surface
                                                                 2  /)  !< 10: urban surf, south )

    REAL(wp)     ::  nc_rad     !< aerosol number concentration
    REAL(IREALS) ::  stime      !< current simulation time
    REAL(IREALS) ::  sundir(3)  !< sun direction vector
    REAL(wp)     ::  s_r2       !< effective particles area to calculate  ts_reliq
    REAL(wp)     ::  s_r3       !< effective particles volume to calculate  ts_reliq
    REAL(IREALS) ::  ts_cliqwp  !< in-cloud liquid water path (g/m2)

    REAL(IREALS), DIMENSION(:), POINTER ::  psktmp  !< pointer for skin temperature interfaces [K]

    REAL(IREALS), DIMENSION(:,:), POINTER ::  ph2ovmr !< pointer for H2O volume mixing ratio
    REAL(IREALS), DIMENSION(:,:), POINTER ::  plwc    !< pointer for Liquid water cloud content [g/kg] and effective radius in micron
    REAL(IREALS), DIMENSION(:,:), POINTER ::  preliq  !< pointer for cloud water drop effective radius (microns)
    REAL(IREALS), DIMENSION(:,:), POINTER ::  psalb2d !< pointer for solar albedo interfaces
    REAL(IREALS), DIMENSION(:,:), POINTER ::  ptalb2d !< pointer for thermal albedo interfaces


    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'radiation_tenstream', time_since_reference_point
       CALL debug_message( debug_string, 'start' )
    ENDIF
!
!-- Calculate the sun direction and the current time
    CALL radiation_calc_sundir( sundir )

!
!-- Calculate the atmosphere inputs for TenStream:
!-- Pressure level
    ts_plev(nzb,:,:) = REAL( surface_pressure, IREALS )
    DO  k = nzb+1, nzt+1
       ts_plev(k,:,:) = REAL( barometric_formula( zw(k), pt_surface * exner(k), surface_pressure ),&
                              IREALS )
       ts_play(k)     = REAL( barometric_formula( zu(k), pt_surface * exner(k), surface_pressure ),&
                              IREALS )
    ENDDO
!
!-- Other quantities
    ts_lwc   = 0.0_IREALS
    ts_reliq = 0.0_IREALS

    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Temperature and H2O volume mixing ratio fields
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             IF ( surf_lsm%upward(m) )                                                             &
                ts_tlev(nzb,i,j) = REAL( surf_lsm%pt_surface(m) * exner(nzb), IREALS )
          ENDDO
          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             IF ( surf_usm%upward(m) )                                                             &
                ts_tlev(nzb,i,j) = REAL( surf_usm%pt_surface(m) * exner(nzb), IREALS )
          ENDDO

          IF ( bulk_cloud_model )  THEN
             DO  k = nzb+1, nzt+1
                ts_tlay(k,i,j)   = REAL( pt(k,j,i) * exner(k) + lv_d_cp * ql(k,j,i), IREALS )
                ts_h2ovmr(k,i,j) = REAL( mol_mass_air_d_wv * ( q(k,j,i) - ql(k,j,i) ), IREALS )
             ENDDO
          ELSEIF ( cloud_droplets )  THEN
             DO  k = nzb+1, nzt+1
                ts_tlay(k,i,j)   = REAL( pt(k,j,i) * exner(k) + lv_d_cp * ql(k,j,i), IREALS )
                ts_h2ovmr(k,i,j) = REAL( mol_mass_air_d_wv * q(k,j,i), IREALS )
             ENDDO
          ELSE
             DO  k = nzb+1, nzt+1
                ts_tlay(k,i,j) = REAL( pt(k,j,i) * exner(k), IREALS )
             ENDDO

             IF ( humidity )  THEN
                DO  k = nzb+1, nzt+1
                   ts_h2ovmr(k,i,j) = REAL( mol_mass_air_d_wv * q(k,j,i), IREALS )
                ENDDO
             ELSE
!
!--             @todo: Actually such an error message should be located in radiation_check_parameters.
                WRITE( message_string, * ) 'tenstream is not supported yet when humidity is false'
                CALL message( 'radiation_tenstream_init', 'RAD0050', 0, 1, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       Linear interpolate to zw grid
          DO  k = nzb+1, nzt
             ts_tlev(k,i,j) = ts_tlay(k,i,j) + ( ts_tlay(k+1,i,j) - ts_tlay(k,i,j) )               &
                                             / ( ts_play(k+1)     - ts_play(k)     )               &
                                             * ( ts_plev(k+1,i,j) - ts_play(k)     )
          ENDDO
          ts_tlev(nzt+1,i,j) = 2.0_IREALS * ts_tlay(nzt+1,i,j) - ts_tlev(nzt,i,j)
!
!--       Calculate liquid water path and cloud fraction for each column.
!--       Note that LWP is required in g/m2 instead of kg/kg m.
          IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN

             CALL message( 'radiation_tenstream_init', 'RAD0050', 0, 1, 0, 6, 0 )

             DO  k = nzb+1, nzt

                IF ( ql(k,j,i) > 1.0E-6 )  ts_lwc(k,i,j) = REAL( ql(k,j,i) * 1000.0_wp, IREALS )
                ts_cliqwp = ts_lwc(k,i,j) * (ts_plev(k,i,j) - ts_plev(k+1,i,j))                    &
                            * REAL( 100.0_wp / g, IREALS )

                IF ( ts_cliqwp > 0.0_wp )  THEN
!
!--                Calculate cloud droplet effective radius
                   IF ( bulk_cloud_model )  THEN
!
!--                   Calculate effective droplet radius. In case of using cloud_scheme = 'morrison'
!--                   and a non reasonable number of cloud droplets the inital aerosol number
!--                   concentration is considered.
                      IF ( microphysics_morrison )  THEN
                         IF ( nc(k,j,i) > 1.0E-20_wp )  THEN
                            nc_rad = nc(k,j,i)
                         ELSE
                            nc_rad = na_init
                         ENDIF
                      ELSE
                         nc_rad = nc_const
                      ENDIF

                      ts_reliq(k,i,j) = REAL( 1.0E6_wp * ( 3.0_wp * ql(k,j,i) * rho_surface /      &
                                                           ( 4.0_wp * pi * nc_rad * rho_l )        &
                                                         )**0.3333333333333_wp                     &
                                              * EXP( LOG( sigma_gc )**2 ), IREALS )

                   ELSEIF ( cloud_droplets )  THEN

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

                      IF ( s_r2 > 0.0_wp )  ts_reliq(k,i,j) = REAL( s_r3 / s_r2, IREALS )

                   ENDIF
!
!--                Limit effective radius
                   IF ( ts_reliq(k,i,j) > 0.0_wp )  THEN
                      ts_reliq(k,i,j) = MAX( ts_reliq(k,i,j),  2.5_IREALS )
                      ts_reliq(k,i,j) = MIN( ts_reliq(k,i,j), 60.0_IREALS )
                   ENDIF

                ENDIF

             ENDDO

          ENDIF
!
!--       Thermal and solar albedo, and skin temperature
!--       Notes:
!--           1) only horizontally aligned surfaces are used
!--           2) weighted average for the different classes is used
!--           3) skin temperature is assumed to be pt_surface
!--           4) solar albedo is set based on rrtm_asdir since no difference is
!--              considered between direct and diffuse albedo
!--           5) thermal albedo is calculated as (1 - emissivity)
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             k = surf_lsm%k(m)
!
!--          Skin temperature
             ts_skin_temperature(i,j)  = MERGE( REAL( surf_lsm%pt_surface(m) * exner(k), IREALS ), &
                                                ts_skin_temperature(i,j),                          &
                                                surf_lsm%upward(m) )
!
!--          Thermal albedo
             ts_thermal_albedo_2d(i,j) = MERGE( REAL( 1.0_wp - SUM( surf_lsm%frac(m,:) *           &
                                                      surf_lsm%emissivity(m,:) ), IREALS ),        &
                                                ts_thermal_albedo_2d(i,j),                         &
                                                surf_lsm%upward(m) )
!
!--          Solar albedo
             ts_solar_albedo_2d(i,j)   = MERGE( REAL( SUM( surf_lsm%frac(m,:) *                    &
                                                      surf_lsm%ts_albedo(m,:) ), IREALS ),         &
                                                ts_solar_albedo_2d(i,j),                           &
                                                surf_lsm%upward(m) )
          ENDDO

          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             k = surf_usm%k(m)
!
!--          Skin temperature
             ts_skin_temperature(i,j)  = MERGE( REAL( surf_usm%pt_surface(m) * exner(k), IREALS ), &
                                                ts_skin_temperature(i,j),                          &
                                                surf_usm%upward(m) )
!
!--          Thermal albedo
             ts_thermal_albedo_2d(i,j) = MERGE( REAL( 1.0_wp - SUM( surf_usm%frac(m,:) *           &
                                                      surf_usm%emissivity(m,:) ), IREALS ),        &
                                                ts_thermal_albedo_2d(i,j),                         &
                                                surf_usm%upward(m) )
!
!--          Solar albedo
             ts_solar_albedo_2d(i,j)   = MERGE( REAL( SUM( surf_usm%frac(m,:) *                    &
                                                      surf_usm%ts_albedo(m,:) ), IREALS ),         &
                                                ts_solar_albedo_2d(i,j),                           &
                                                surf_usm%upward(m) )
          ENDDO

       ENDDO ! nys
    ENDDO ! nxl

!
!-- Debug information to show the PALM settings which transferred to tenstream.
    IF ( debug_output )  THEN
       WRITE( 9, * ) '*** PALM general information transferred to tenstream:'
       WRITE( 9, * ) 'numnodes = numprocs', numprocs
       WRITE( 9, * ) 'ts_nranksx, ts_nranksy', ts_nranksx, ts_nranksy
       WRITE( 9, * ) 'ts_nxproc, ts_nyproc', ts_nxproc, ts_nyproc
       WRITE( 9, * ) 'X: dx = ', dx,',   nnx = ', nnx, ', nxl = ',nxl, ', nxr = ', nxr
       WRITE( 9, * ) 'Y: dy = ', dy,',   nny = ', nny, ', nys = ',nys, ', nyn = ', nyn
       WRITE( 9, * ) 'Z:dz(1) = ', dz(1), ', nnz =', nnz, ', nzb = ', nzb, ', nzt = ', nzt
       WRITE( 9, * ) 'size(zu)', SIZE(zu)
       WRITE( 9, * ) 'size(zw)', SIZE(zw)
       WRITE( 9, * ) 'sundir', sundir
       WRITE( 9, * ) 'albedo_th = ', albedo_th, 'albedo_sol = ', albedo_sol
       WRITE( 9, * ) '*** palm: atmosphere input to tenstream at | i:', nxl, 'j:', nys
       WRITE( 9, * ) '***       k','    ts_plev  ','      ts_tlev  ', '           ts_lwc   ',      &
                     '    ts_reliq  ','    ts_h2ovmr '
       DO k = nzb+1, nzt+1
          WRITE( 9, * ) k, ts_plev(k,nxl,nys), ts_tlev(k,nxl,nys), ts_lwc(k,nxl,nys),              &
                        ts_reliq(k,nxl,nys), ts_h2ovmr(k,nxl,nys)
       ENDDO
       WRITE( 9, * ) '*** min/max    ','    ts_plev  ','      ts_tlev  ','           ts_lwc   ',   &
                     '    ts_reliq  ','    ts_h2ovmr '
       WRITE( 9, * ) 'min',  MINVAL(ts_plev), MINVAL(ts_tlev), MINVAL(ts_lwc), MINVAL(ts_reliq),   &
                     MINVAL(ts_h2ovmr)
       WRITE( 9, * ) 'max',  MAXVAL(ts_plev), MAXVAL(ts_tlev), MAXVAL(ts_lwc), MAXVAL(ts_reliq),   &
                     MAXVAL(ts_h2ovmr)
       FLUSH( 9 )
    ENDIF
!
!-- Set pointers
    pplev(1:SIZE(ts_plev,1),1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_plev
    ptlev(1:SIZE(ts_plev,1),1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_tlev
    ptlay(1:SIZE(ts_plev,1)-1,1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_tlay
    ph2ovmr(1:SIZE(ts_plev,1)-1,1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_h2ovmr
    plwc(1:SIZE(ts_plev,1)-1,1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_lwc
    preliq(1:SIZE(ts_plev,1)-1,1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_reliq
    psktmp(1:SIZE(ts_plev,2)*SIZE(ts_plev,3)) => ts_skin_temperature
    ptalb2d(1:SIZE(ts_plev,2),1:SIZE(ts_plev,3)) => ts_thermal_albedo_2d
    psalb2d(1:SIZE(ts_plev,2),1:SIZE(ts_plev,3)) => ts_solar_albedo_2d

    CALL SETUP_TENSTR_ATM( ts_comm, .FALSE., ts_atm_filename, pplev, ptlev, ts_atm, d_tlay = ptlay,&
                           d_h2ovmr = ph2ovmr, d_lwc = plwc, d_reliq = preliq,                     &
                           d_skin_temperature = psktmp )
    stime = REAL( simulated_time, IREALS )

!
!-- Exchange the boarder data (albedo, surface temperature)
    request_count = 0
!
!-- 1) Send receive signals from boarders
    IF ( nfacad_eastg > 0 )  THEN
       request_count = request_count + 1
       CALL MPI_IRECV( surf_eg, 3*nfacad_eastg, MPI_REAL, pright, tag_w, MPI_COMM_WORLD,           &
                       requests(request_count), ierr )
    ENDIF

    IF ( nfacad_westg > 0 )  THEN
       request_count = request_count + 1
       CALL MPI_IRECV( surf_wg, 3*nfacad_westg, MPI_REAL, pleft, tag_e, MPI_COMM_WORLD,            &
                       requests(request_count), ierr )
    ENDIF

    IF ( nfacad_northg > 0 )  THEN
       request_count = request_count + 1
       CALL MPI_IRECV( surf_ng, 3*nfacad_northg, MPI_REAL, pnorth, tag_s, MPI_COMM_WORLD,          &
                       requests(request_count), ierr )
    ENDIF

    IF ( nfacad_southg > 0 )  THEN
       request_count = request_count + 1
       CALL MPI_IRECV( surf_sg, 3*nfacad_southg, MPI_REAL, psouth, tag_n, MPI_COMM_WORLD,          &
                       requests(request_count), ierr )
    ENDIF
!
!-- Send to boarders
!-- East
    IF ( nfacad_east > 0 )  THEN
       CALL fill_surf_b( nfacad_east, surf_ids_e, surf_e )
       request_count = request_count + 1
       CALL MPI_ISEND( surf_e, 3*nfacad_east, MPI_REAL, pright, tag_e, MPI_COMM_WORLD,             &
                       requests(request_count), ierr )
    ENDIF
!
!-- West
    IF ( nfacad_west > 0 )  THEN
       CALL fill_surf_b( nfacad_west, surf_ids_w, surf_w )
       request_count = request_count + 1
       CALL MPI_ISEND( surf_w, 3*nfacad_west, MPI_REAL, pleft, tag_w, MPI_COMM_WORLD,              &
                       requests(request_count), ierr )
    ENDIF
!
!-  North
    IF ( nfacad_north > 0 )  THEN
       CALL fill_surf_b( nfacad_north, surf_ids_n, surf_n )
       request_count = request_count + 1
       CALL MPI_ISEND( surf_n, 3*nfacad_north, MPI_REAL, pnorth, tag_n, MPI_COMM_WORLD,            &
                       requests(request_count), ierr )
    ENDIF
!
!-- South
    IF ( nfacad_south > 0 )  THEN
       CALL fill_surf_b( nfacad_south, surf_ids_s, surf_s )
       request_count = request_count + 1
       CALL MPI_ISEND( surf_s, 3*nfacad_south, MPI_REAL, psouth, tag_s, MPI_COMM_WORLD,            &
                       requests(request_count), ierr )
    ENDIF

    CALL MPI_WAITALL( request_count, requests(1:request_count), wait_stat, ierr )

!
!-- Calculate albedo and surface temperature for facades:
!-- 1) at boarders
    IF ( nfacad_eastg  > 0 )  CALL map_surf_bg( nfacad_eastg,  surf_ids_eg, surf_eg )
    IF ( nfacad_westg  > 0 )  CALL map_surf_bg( nfacad_westg,  surf_ids_wg, surf_wg )
    IF ( nfacad_northg > 0 )  CALL map_surf_bg( nfacad_northg, surf_ids_ng, surf_ng )
    IF ( nfacad_southg > 0 )  CALL map_surf_bg( nfacad_southg, surf_ids_sg, surf_sg )
!
!-- 2) inside domain
    DO  ifacad = 1, SIZE( surf_ids, 2 )

       i    = surf_ids(1,ifacad)
       m    = surf_ids(2,ifacad)
       ts_m = surf_ids(3,ifacad)

       isurf_type = facad_type(i)

       SELECT CASE ( isurf_type )
!
!--       Land surface
          CASE ( 1 )
             k = surf_lsm%k(m)
             buildings_solar%albedo(ts_m) = REAL( SUM( surf_lsm%frac(m,:) *                        &
                                                       surf_lsm%ts_albedo(m,:) ), IREALS )
             buildings_thermal%albedo(ts_m) = REAL( 1.0_wp - SUM( surf_lsm%frac(m,:) *             &
                                                    surf_lsm%emissivity(m,:) ), IREALS )
             buildings_thermal%temp(ts_m) = REAL( surf_lsm%pt_surface(m) * exner(k), IREALS )
!
!--       Urban surface
          CASE ( 2 )
             k = surf_usm%k(m)
             buildings_solar%albedo(ts_m) = REAL( SUM( surf_usm%frac(m,:) *                        &
                                                       surf_usm%ts_albedo(m,:) ), IREALS )
             buildings_thermal%albedo(ts_m) = REAL( 1.0_wp - SUM( surf_usm%frac(m,:) *             &
                                                    surf_usm%emissivity(m,:) ), IREALS )
             buildings_thermal%temp(ts_m) = REAL( surf_usm%pt_surface(m) * exner(k), IREALS )

       END SELECT

    ENDDO

!
!-- Calculate LW flux and heatrate
    IF ( lw_radiation )  THEN
       CALL PPRTS_RRTMG( ts_comm, ts_solver, ts_atm,                                               &
                         INT( nnx, IINTEGERS ), INT( nny, IINTEGERS ),                             &
                         ts_dx, ts_dy, sundir,                                                     &
                         albedo_th, albedo_sol,                                                    &
                         lw_radiation, .FALSE.,                                                    &
                         ts_edir, ts_edn, ts_eup, ts_abso,                                         &
                         nxproc = ts_nxproc, nyproc = ts_nyproc,                                   &
                         opt_time = stime , solar_albedo_2d = psalb2d,                             &
                         thermal_albedo_2d = ptalb2d,                                              &
                         opt_solar_constant = norm2(sundir),                                       &
                         icollapse = ts_icollapse,                                                 &
                         opt_buildings_solar   = buildings_solar,                                  &
                         opt_buildings_thermal = buildings_thermal,                                &
                         opt_tau_solar = tree_tau_solar,                                           &
                         opt_w0_solar = tree_w0_solar,                                             &
                         opt_tau_thermal = tree_tau_thermal )
!
!--    Save LW flux and heatrate
       nlev = UBOUND( ts_edn, 1 )
       nlay = UBOUND( ts_abso, 1 )

       DO  i = nxl, nxr
          ts_i = i - nxl + 1
          DO  j = nys, nyn
             k_topo = topo_top_ind(j,i,0)
             ts_j = j - nys + 1
             DO  k = k_topo+1, nzt+1
!
!--             1) atmosphere
!--             LW fluxes
                ts_k = nlev - ( k - nzb )
                rad_lw_in(k,j,i)  = ts_edn(ts_k,ts_i,ts_j)
                rad_lw_out(k,j,i) = ts_eup(ts_k,ts_i,ts_j)
!
!--             LW heatrate (convert from W/m-3 to K/h)
                ts_k = nlay - ( k - ( nzb + 1 ) )
!
!--             To convert from W/m-3 to K/h, ts_abso*3600/(cp*rho). Here rho is diagnostically
!--             calculated from local pressure and temperature values.
                rad_lw_hr(k,j,i) = ts_abso(ts_k,ts_i,ts_j) * seconds_per_hour * rd_d_cp * exner(k) &
                                   * pt(k,j,i) / hyp(k)
             ENDDO
!
!--          Flux at the first level
             ts_k = nlev - ( k_topo - nzb )
             rad_lw_in(k_topo,j,i)  = ts_edn(ts_k,ts_i,ts_j)
             rad_lw_out(k_topo,j,i) = ts_eup(ts_k,ts_i,ts_j)
!
!--          2) surface arrayes (non-buildings non orography)
!
!--          Land
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                IF ( surf_lsm%upward(m) )  THEN
                   k = surf_lsm%k(m)
                   ts_k = nlev - ( k - nzb )
                   surf_lsm%rad_lw_in(m)  = ts_edn(ts_k,ts_i,ts_j)
                   surf_lsm%rad_lw_out(m) = ts_eup(ts_k,ts_i,ts_j)
                ENDIF
             ENDDO
!
!--          Urban
             DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                IF ( surf_usm%upward(m) )  THEN
                   k = surf_usm%k(m)
                   ts_k = nlev - ( k - nzb )
                   surf_usm%rad_lw_in(m)  = ts_edn(ts_k,ts_i,ts_j)
                   surf_usm%rad_lw_out(m) = ts_eup(ts_k,ts_i,ts_j)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    3) surfaces (buildings and orography)
       DO  ifacad = 1, SIZE( surf_ids, 2 )

          i    = surf_ids(1,ifacad)
          m    = surf_ids(2,ifacad)
          ts_m = surf_ids(3,ifacad)

          isurf_type = facad_type(i)

          SELECT CASE ( isurf_type )
!
!--          Land surface
             CASE ( 1 )
                surf_lsm%rad_lw_in(m)  = buildings_thermal%incoming(ts_m)
                surf_lsm%rad_lw_out(m) = buildings_thermal%outgoing(ts_m)
!
!--          Urban surface
             CASE ( 2 )
                surf_usm%rad_lw_in(m)  = buildings_thermal%incoming(ts_m)
                surf_usm%rad_lw_out(m) = buildings_thermal%outgoing(ts_m)

          END SELECT

       ENDDO

!
!--    4) surfaces (buildings and orography) at boundaries
       request_count = 0
!
!--    4.1) Send receive signals from boarders
!--    East
       IF ( nfacad_east > 0 )  THEN
          request_count = request_count + 1
          CALL MPI_IRECV( surf_e, 3*nfacad_east, MPI_REAL, pright, tag_w, MPI_COMM_WORLD,          &
                          requests(request_count), ierr )
       ENDIF
!
!--    West
       IF ( nfacad_west > 0 )  THEN
          request_count = request_count + 1
          CALL MPI_IRECV( surf_w, 3*nfacad_west, MPI_REAL, pleft, tag_e, MPI_COMM_WORLD,           &
                          requests(request_count), ierr )
       ENDIF
!
!--    North
       IF ( nfacad_north > 0 )  THEN
          request_count = request_count + 1
          CALL MPI_IRECV( surf_n, 3*nfacad_north, MPI_REAL, pnorth, tag_s, MPI_COMM_WORLD,         &
                          requests(request_count), ierr )
       ENDIF
!
!--    South
       IF ( nfacad_south > 0 ) THEN
          request_count = request_count + 1
          CALL MPI_IRECV( surf_s, 3*nfacad_south, MPI_REAL, psouth, tag_n, MPI_COMM_WORLD,         &
                          requests(request_count), ierr )
       ENDIF
!
!--    4.2) Send to boarders
!-     East
       IF ( nfacad_eastg > 0 )  THEN
          DO  ifacad = 1, nfacad_eastg
             ts_m = surf_ids_eg(2,ifacad)
             surf_eg(1,ifacad) = buildings_thermal%incoming(ts_m)
             surf_eg(2,ifacad) = buildings_thermal%outgoing(ts_m)
          ENDDO
          request_count = request_count + 1
          CALL MPI_ISEND( surf_eg, 3*nfacad_eastg, MPI_REAL, pright, tag_e, MPI_COMM_WORLD,        &
                          requests(request_count), ierr )
       ENDIF
!
!--    West
       IF ( nfacad_westg > 0 )  THEN
          DO  ifacad = 1, nfacad_westg
             ts_m = surf_ids_wg(2,ifacad)
             surf_wg(1,ifacad) = buildings_thermal%incoming(ts_m)
             surf_wg(2,ifacad) = buildings_thermal%outgoing(ts_m)
          ENDDO
          request_count = request_count + 1
          CALL MPI_ISEND( surf_wg, 3*nfacad_westg, MPI_REAL, pleft, tag_w, MPI_COMM_WORLD,         &
                          requests(request_count), ierr )
       ENDIF
!
!--    North
       IF ( nfacad_northg > 0 )  THEN
          DO  ifacad = 1, nfacad_northg
             ts_m = surf_ids_ng(2,ifacad)
             surf_ng(1,ifacad) = buildings_thermal%incoming(ts_m)
             surf_ng(2,ifacad) = buildings_thermal%outgoing(ts_m)
          ENDDO
          request_count = request_count + 1
          CALL MPI_ISEND( surf_ng, 3*nfacad_northg, MPI_REAL, pnorth, tag_n, MPI_COMM_WORLD,       &
                          requests(request_count), ierr )
       ENDIF
!
!--    South
       IF ( nfacad_southg > 0 )  THEN
          DO  ifacad = 1, nfacad_southg
             ts_m = surf_ids_sg(2,ifacad)
             surf_sg(1,ifacad) = buildings_thermal%incoming(ts_m)
             surf_sg(2,ifacad) = buildings_thermal%outgoing(ts_m)
          ENDDO
          request_count = request_count + 1
          CALL MPI_ISEND( surf_sg, 3*nfacad_southg, MPI_REAL, psouth, tag_s, MPI_COMM_WORLD,       &
                          requests(request_count), ierr )
       ENDIF
!
!--    Wait
       CALL MPI_WAITALL( request_count, requests(1:request_count), wait_stat, ierr )
!
!--    4.3) map received surface fluxes back to respective surfaces
!--    East
       IF ( nfacad_east > 0 )  THEN
          CALL map_surf_b( nfacad_east, surf_ids_e, surf_e, .FALSE. )
       ENDIF
!
!--    West
       IF ( nfacad_west > 0 )  THEN
          CALL map_surf_b( nfacad_west, surf_ids_w, surf_w, .FALSE. )
       ENDIF
!
!--    North
       IF ( nfacad_north > 0 )  THEN
          CALL map_surf_b( nfacad_north, surf_ids_n, surf_n, .FALSE. )
       ENDIF
!
!--    South
       IF ( nfacad_south > 0 )  THEN
          CALL map_surf_b( nfacad_south, surf_ids_s, surf_s, .FALSE. )
       ENDIF

       IF ( debug_output )  THEN
          WRITE( 9, * ) '*** ts: radiation input profiles PPRTS_RRTMG at i:', nxl, 'j:', nys
          WRITE( 9, * ) '***       k','    ts_plev  ','      ts_tlev  ','           ts_lwc   ',    &
                        '    ts_reliq  ','    ts_h2ovmr ','         rad_lw_hr ',                   &
                        '           rad_lw_hr/cp'
          DO  k = nzb+1, nzt+1
             WRITE( 9, * ) k, ts_plev(k, nxl, nys), ts_tlev(k, nxl, nys), ts_lwc(k, nxl, nys),     &
                           ts_reliq(k, nxl, nys), ts_h2ovmr(k, nxl, nys), rad_lw_hr(k, nys, nxl),  &
                           rad_lw_hr(k, nys, nxl) / c_p
          ENDDO
          FLUSH( 9 )
       ENDIF

    ELSE

       rad_lw_in  = 0.0_wp
       rad_lw_out = 0.0_wp
       rad_lw_hr  = 0.0_wp
       IF ( surf_lsm%ns > 0 )  THEN
          surf_lsm%rad_lw_in  = 0.0_wp
          surf_lsm%rad_lw_out = 0.0_wp
       ENDIF
       IF ( surf_usm%ns > 0 )  THEN
          surf_usm%rad_lw_in  = 0.0_wp
          surf_usm%rad_lw_out = 0.0_wp
       ENDIF

    ENDIF ! if lw_radiation

    IF ( sw_radiation  .AND. sun_up )  THEN

       CALL PPRTS_RRTMG( ts_comm, ts_solver, ts_atm,                                               &
                         INT( nnx, IINTEGERS ), INT( nny, IINTEGERS ),                             &
                         ts_dx, ts_dy, sundir,                                                     &
                         albedo_th, albedo_sol,                                                    &
                         .FALSE., sw_radiation,                                                    &
                         ts_edir, ts_edn, ts_eup, ts_abso,                                         &
                         nxproc = ts_nxproc, nyproc = ts_nyproc,                                   &
                         opt_time = stime, solar_albedo_2d = psalb2d,                              &
                         thermal_albedo_2d = ptalb2d,                                              &
                         opt_solar_constant = norm2(sundir),                                       &
                         icollapse = ts_icollapse,                                                 &
                         opt_buildings_solar = buildings_solar,                                    &
                         opt_buildings_thermal = buildings_thermal,                                &
                         opt_tau_solar = tree_tau_solar,                                           &
                         opt_w0_solar = tree_w0_solar,                                             &
                         opt_tau_thermal = tree_tau_thermal )
!
!--    Save SW flux and heatrate
       nlev = UBOUND( ts_edn, 1 )
       nlay = UBOUND( ts_abso, 1 )

       DO  i = nxl, nxr
          ts_i = i - nxl + 1
          DO  j = nys, nyn
             k_topo = topo_top_ind(j,i,0)
             ts_j = j - nys + 1
             DO  k = k_topo+1, nzt+1
!
!--             1) atmosphere
!--             SW fluxes
                ts_k = nlev - ( k - nzb )
                rad_sw_in(k,j,i)  = ts_edn(ts_k,ts_i,ts_j) + ts_edir(ts_k,ts_i,ts_j)
                rad_sw_out(k,j,i) = ts_eup(ts_k,ts_i,ts_j)
!
!--             SW heatrate (convert from W/m-3 to K/h)
                ts_k = nlay - ( k - ( nzb + 1 ) )
!
!--             To convert from W/m-3 to K/h, ts_abso*3600/(cp*rho). Here rho is diagnostically
!--             calculated from local pressure and temperature values.
                rad_sw_hr(k,j,i) = ts_abso(ts_k,ts_i,ts_j) * seconds_per_hour * rd_d_cp * exner(k) &
                                   * pt(k,j,i) / hyp(k)
             ENDDO
!
!--          Flux at the first level
             ts_k = nlev - ( k_topo - nzb )
             rad_sw_in(k_topo,j,i)  = ts_edn(ts_k,ts_i,ts_j) + ts_edir(ts_k,ts_i,ts_j)
             rad_sw_out(k_topo,j,i) = ts_eup(ts_k,ts_i,ts_j)
!
!--          Save fluxes to surface arrays.
!--          2) surface arrayes (non-buildings non orography)
!
!--          Land
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                IF ( surf_lsm%upward(m) )  THEN
                   k = surf_lsm%k(m)
                   ts_k = nlev - ( k - nzb )
                   surf_lsm%rad_sw_in(m)  = ts_edn(ts_k,ts_i,ts_j) + ts_edir(ts_k,ts_i,ts_j)
                   surf_lsm%rad_sw_out(m) = ts_eup(ts_k,ts_i,ts_j)
                   surf_lsm%rad_sw_dir(m) = ts_edir(ts_k,ts_i,ts_j)
                   surf_lsm%rad_sw_dif(m) = ts_edn(ts_k,ts_i,ts_j)

                   surf_lsm%rad_net(m) = surf_lsm%rad_sw_in (m) -                                  &
                                         surf_lsm%rad_sw_out(m) +                                  &
                                         surf_lsm%rad_lw_in(m)  -                                  &
                                         surf_lsm%rad_lw_out(m)
                ENDIF
             ENDDO
!
!--          Urban
             DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                IF ( surf_usm%upward(m) )  THEN
                   k = surf_usm%k(m)
                   ts_k = nlev - ( k - nzb )
                   surf_usm%rad_sw_in(m)  = ts_edn(ts_k,ts_i,ts_j) + ts_edir(ts_k,ts_i,ts_j)
                   surf_usm%rad_sw_out(m) = ts_eup(ts_k,ts_i,ts_j)
                   surf_usm%rad_sw_dir(m) = ts_edir(ts_k,ts_i,ts_j)
                   surf_usm%rad_sw_dif(m) = ts_edn(ts_k,ts_i,ts_j)

                   surf_usm%rad_net(m) = surf_usm%rad_sw_in (m) -                                  &
                                         surf_usm%rad_sw_out(m) +                                  &
                                         surf_usm%rad_lw_in(m)  -                                  &
                                         surf_usm%rad_lw_out(m)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       IF ( debug_output )  THEN
          WRITE( 9, * ) '*** ts: radiation input profiles PPRTS_RRTMG at i:', nxl, 'j:', nys
          WRITE( 9, * ) '***       k','    ts_plev  ','      ts_tlev  ', '           ts_lwc   ',   &
                        '    ts_reliq  ','    ts_h2ovmr ', '         rad_sw_hr ',                  &
                        '           rad_sw_hr/cp'
          DO k = nzb+1, nzt+1
             WRITE( 9, * ) k, ts_plev(k, nxl, nys), ts_tlev(k, nxl, nys), ts_lwc(k, nxl, nys),     &
                           ts_reliq(k, nxl, nys), ts_h2ovmr(k, nxl, nys), rad_sw_hr(k, nys, nxl),  &
                           rad_sw_hr(k, nys, nxl) / c_p
          ENDDO
          FLUSH( 9 )
       ENDIF

!
!--    3) surfaces (buildings and orography)
       DO  ifacad = 1, SIZE( surf_ids, 2 )

          i    = surf_ids(1,ifacad)
          m    = surf_ids(2,ifacad)
          ts_m = surf_ids(3,ifacad)

          isurf_type = facad_type(i)

          SELECT CASE ( isurf_type )
!
!--          Land surface
             CASE ( 1 )
                surf_lsm%rad_sw_dir(m) = buildings_solar%edir(ts_m)
                surf_lsm%rad_sw_dif(m) = buildings_solar%incoming(ts_m)
                surf_lsm%rad_sw_in(m)  = buildings_solar%edir(ts_m)  +                             &
                                         buildings_solar%incoming(ts_m)
                surf_lsm%rad_sw_out(m) = buildings_solar%outgoing(ts_m)
                surf_lsm%rad_net(m)    = surf_lsm%rad_sw_in(m)  -                                  &
                                         surf_lsm%rad_sw_out(m) +                                  &
                                         surf_lsm%rad_lw_in(m)  -                                  &
                                         surf_lsm%rad_lw_out(m)
!
!--          Urban surface
             CASE ( 2 )
                surf_usm%rad_sw_dir(m) = buildings_solar%edir(ts_m)
                surf_usm%rad_sw_dif(m) = buildings_solar%incoming(ts_m)
                surf_usm%rad_sw_in(m)  = buildings_solar%edir(ts_m)  +                             &
                                         buildings_solar%incoming(ts_m)
                surf_usm%rad_sw_out(m) = buildings_solar%outgoing(ts_m)
                surf_usm%rad_net(m)    = surf_usm%rad_sw_in(m)  -                                  &
                                         surf_usm%rad_sw_out(m) +                                  &
                                         surf_usm%rad_lw_in(m)  -                                  &
                                         surf_usm%rad_lw_out(m)

          END SELECT

       ENDDO
!
!--    4) surfaces (buildings and orography) at boundaries
       request_count = 0
!
!--    4.1) Send receive signals from boarders
!--    East
       IF ( nfacad_east > 0 )  THEN
          request_count = request_count + 1
          CALL MPI_IRECV( surf_e, 3*nfacad_east, MPI_REAL, pright, tag_w, MPI_COMM_WORLD,          &
                          requests(request_count), ierr )
       ENDIF
!
!--    West
       IF ( nfacad_west > 0 )  THEN
          request_count = request_count + 1
          CALL MPI_IRECV( surf_w, 3*nfacad_west, MPI_REAL, pleft, tag_e, MPI_COMM_WORLD,           &
                          requests(request_count), ierr )
       ENDIF
!
!--    North
       IF ( nfacad_north > 0 )  THEN
          request_count = request_count + 1
          CALL MPI_IRECV( surf_n, 3*nfacad_north, MPI_REAL, pnorth, tag_s, MPI_COMM_WORLD,         &
                          requests(request_count), ierr )
       ENDIF
!
!--    South
       IF ( nfacad_south > 0 )  THEN
          request_count = request_count + 1
          CALL MPI_IRECV( surf_s, 3*nfacad_south, MPI_REAL, psouth, tag_n, MPI_COMM_WORLD,         &
                          requests(request_count), ierr )
       ENDIF
!
!--    4.2) Send to boarders
!--    East
       IF ( nfacad_eastg > 0 )  THEN
          DO  ifacad = 1, nfacad_eastg
             ts_m = surf_ids_eg(2,ifacad)
             surf_eg(1,ifacad) = buildings_solar%edir(ts_m)
             surf_eg(2,ifacad) = buildings_solar%incoming(ts_m)
             surf_eg(3,ifacad) = buildings_solar%outgoing(ts_m)
          ENDDO
          request_count = request_count + 1
          CALL MPI_ISEND( surf_eg, 3*nfacad_eastg, MPI_REAL, pright, tag_e, MPI_COMM_WORLD,        &
                          requests(request_count), ierr )
       ENDIF
!
!--    West
       IF ( nfacad_westg > 0 )  THEN
          DO  ifacad = 1, nfacad_westg
             ts_m = surf_ids_wg(2,ifacad)
             surf_wg(1,ifacad) = buildings_solar%edir(ts_m)
             surf_wg(2,ifacad) = buildings_solar%incoming(ts_m)
             surf_wg(3,ifacad) = buildings_solar%outgoing(ts_m)
          ENDDO
          request_count = request_count + 1
          CALL MPI_ISEND( surf_wg, 3*nfacad_westg, MPI_REAL, pleft, tag_w, MPI_COMM_WORLD,         &
                          requests(request_count), ierr )
       ENDIF
!
!--    North
       IF ( nfacad_northg > 0 )  THEN
          DO  ifacad = 1, nfacad_northg
             ts_m = surf_ids_ng(2,ifacad)
             surf_ng(1,ifacad) = buildings_solar%edir(ts_m)
             surf_ng(2,ifacad) = buildings_solar%incoming(ts_m)
             surf_ng(3,ifacad) = buildings_solar%outgoing(ts_m)
          ENDDO
          request_count = request_count + 1
          CALL MPI_ISEND( surf_ng, 3*nfacad_northg, MPI_REAL, pnorth, tag_n, MPI_COMM_WORLD,       &
                          requests(request_count), ierr )
       ENDIF
!
!--    South
       IF ( nfacad_southg > 0 )  THEN
          DO  ifacad = 1, nfacad_southg
             ts_m = surf_ids_sg(2,ifacad)
             surf_sg(1,ifacad) = buildings_solar%edir(ts_m)
             surf_sg(2,ifacad) = buildings_solar%incoming(ts_m)
             surf_sg(3,ifacad) = buildings_solar%outgoing(ts_m)
          ENDDO
          request_count = request_count + 1
          CALL MPI_ISEND( surf_sg, 3*nfacad_southg, MPI_REAL, psouth, tag_s, MPI_COMM_WORLD,       &
                          requests(request_count), ierr )
       ENDIF

       CALL MPI_WAITALL( request_count, requests(1:request_count), wait_stat, ierr )
!
!--    4.3) map received surface fluxes back to respective surfaces
!--    East
       IF ( nfacad_east > 0 )  THEN
          CALL map_surf_b( nfacad_east, surf_ids_e, surf_e, .TRUE. )
       ENDIF
!
!--    West
       IF ( nfacad_west > 0 )  THEN
          CALL map_surf_b( nfacad_west, surf_ids_w, surf_w, .TRUE. )
       ENDIF
!
!--    North
       IF ( nfacad_north > 0 )  THEN
          CALL map_surf_b( nfacad_north, surf_ids_n, surf_n, .TRUE. )
       ENDIF
!
!--    South
       IF ( nfacad_south > 0 )  THEN
          CALL map_surf_b( nfacad_south, surf_ids_s, surf_s, .TRUE. )
       ENDIF
!
!-- Sun is down (sun_up is false)
    ELSE

       rad_sw_in  = 0.0_wp
       rad_sw_out = 0.0_wp
       rad_sw_hr  = 0.0_wp
!
!--    Horizontal surfaces
       IF ( surf_lsm%ns > 0 )  THEN
          surf_lsm%rad_sw_dir = 0.0_wp
          surf_lsm%rad_sw_dif = 0.0_wp
          surf_lsm%rad_sw_in  = 0.0_wp
          surf_lsm%rad_sw_out = 0.0_wp
          surf_lsm%rad_net = surf_lsm%rad_lw_in - surf_lsm%rad_lw_out
       ENDIF
       IF ( surf_usm%ns > 0 )  THEN
          surf_usm%rad_sw_dir = 0.0_wp
          surf_usm%rad_sw_dif = 0.0_wp
          surf_usm%rad_sw_in  = 0.0_wp
          surf_usm%rad_sw_out = 0.0_wp
          surf_usm%rad_net = surf_usm%rad_lw_in - surf_usm%rad_lw_out
       ENDIF

    ENDIF ! if SW and sun_up

    CALL exchange_horiz( rad_lw_in, nbgp )
    CALL exchange_horiz( rad_lw_out, nbgp )
    CALL exchange_horiz( rad_lw_hr, nbgp )

    CALL exchange_horiz( rad_sw_in, nbgp )
    CALL exchange_horiz( rad_sw_out, nbgp )
    CALL exchange_horiz( rad_sw_hr, nbgp )

    IF ( ( simulated_time + dt_radiation ) > end_time )  THEN
       IF ( ALLOCATED( ts_solver ) )  THEN
          CALL DESTROY_PPRTS_RRTMG( ts_solver, lfinalizepetsc=.TRUE. )
          CALL DESTROY_TENSTR_ATM( ts_atm )
       ENDIF
    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'radiation_tenstream', 'end' )

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> fill_surf_b fills the array surf_b, where b is the side, e.g. east, west, etc
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fill_surf_b( nfacad_b, surf_ids_b, surf_b )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  nfacad_b  !<

    INTEGER(iwp), DIMENSION(:,:), INTENT(IN) ::  surf_ids_b  !<

    REAL(wp), DIMENSION(:,:), INTENT(INOUT) ::  surf_b  !<


    DO  ifacad = 1, nfacad_b

       m          = surf_ids_b(1,ifacad)
       isurf_type = surf_ids_b(2,ifacad)

       SELECT CASE ( isurf_type )
!
!--       Land surface
          CASE ( 1 )
             k = surf_lsm%k(m)
             surf_b(1,ifacad) = SUM( surf_lsm%frac(m,:) * surf_lsm%ts_albedo(m,:) )
             surf_b(2,ifacad) = 1.0_wp - SUM( surf_lsm%frac(m,:) * surf_lsm%emissivity(m,:) )
             surf_b(3,ifacad) = surf_lsm%pt_surface(m) * exner(k)
!
!--       Urban surface
          CASE ( 2 )
             k = surf_usm%k(m)
             surf_b(1,ifacad) = SUM( surf_usm%frac(m,:) * surf_usm%ts_albedo(m,:) )
             surf_b(2,ifacad) = 1.0_wp - SUM( surf_usm%frac(m,:) * surf_usm%emissivity(m,:) )
             surf_b(3,ifacad) = surf_usm%pt_surface(m) * exner(k)

       END SELECT

    ENDDO

 END SUBROUTINE fill_surf_b


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> map_surf_b map the array surf_b, where b is the side, back to the surface structure surf_l/usm_v
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE map_surf_b( nfacad_b, surf_ids_b, surf_b, l_sw )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  nfacad_b  !<

    INTEGER(iwp), DIMENSION(:,:), INTENT(IN) ::  surf_ids_b  !<

    LOGICAL, INTENT(IN) ::  l_sw  !<

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::  surf_b  !<


    IF ( l_sw )  THEN
!
!--    SW
       DO  ifacad = 1, nfacad_b
          m          = surf_ids_b(1,ifacad)
          isurf_type = surf_ids_b(2,ifacad)

          SELECT CASE ( isurf_type )
!
!--          Land surface
             CASE ( 1 )
                surf_lsm%rad_sw_dir(m) = surf_b(1,ifacad)
                surf_lsm%rad_sw_dif(m) = surf_b(2,ifacad)
                surf_lsm%rad_sw_out(m) = surf_b(3,ifacad)
                surf_lsm%rad_sw_in(m)  = surf_lsm%rad_sw_dir(m) +                                  &
                                         surf_lsm%rad_sw_dif(m)
                surf_lsm%rad_net(m)    = surf_lsm%rad_sw_in(m)  -                                  &
                                         surf_lsm%rad_sw_out(m) +                                  &
                                         surf_lsm%rad_lw_in(m)  -                                  &
                                         surf_lsm%rad_lw_out(m)
!
!--          Urban surface
             CASE ( 2 )
                surf_usm%rad_sw_dir(m) = surf_b(1,ifacad)
                surf_usm%rad_sw_dif(m) = surf_b(2,ifacad)
                surf_usm%rad_sw_out(m) = surf_b(3,ifacad)
                surf_usm%rad_sw_in(m)  = surf_usm%rad_sw_dir(m) +                                  &
                                         surf_usm%rad_sw_dif(m)
                surf_usm%rad_net(m)    = surf_usm%rad_sw_in(m)  -                                  &
                                         surf_usm%rad_sw_out(m) +                                  &
                                         surf_usm%rad_lw_in(m)  -                                  &
                                         surf_usm%rad_lw_out(m)

          END SELECT

       ENDDO

    ELSE
!
!--    LW
       DO  ifacad = 1, nfacad_b
          m          = surf_ids_b(1,ifacad)
          isurf_type = surf_ids_b(2,ifacad)

          SELECT CASE ( isurf_type )
!
!--          Land surface
             CASE ( 1 )
                surf_lsm%rad_lw_in(m)  = surf_b(1,ifacad)
                surf_lsm%rad_lw_out(m) = surf_b(2,ifacad)
!
!--          Urban surface
             CASE ( 2 )
                surf_usm%rad_lw_in(m)  = surf_b(1,ifacad)
                surf_usm%rad_lw_out(m) = surf_b(2,ifacad)

         END SELECT

      ENDDO

   ENDIF

 END SUBROUTINE map_surf_b


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> map_surf_b map the array at boarder surf_b, where b is the side, back to the TS structure
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE map_surf_bg( nfacad_b, surf_ids_b, surf_b )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  nfacad_b  !<
    INTEGER(iwp), DIMENSION(:,:), INTENT(IN) ::  surf_ids_b  !<

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::  surf_b  !<

    INTEGER(IINTEGERS) :: ts_m


    DO  ifacad = 1, nfacad_b
       ts_m = INT( surf_ids_b(2,ifacad), IINTEGERS )
       buildings_thermal%albedo(ts_m) = REAL( surf_b(2,ifacad), IREALS )
       buildings_thermal%temp(ts_m)   = REAL( surf_b(3,ifacad), IREALS )
    ENDDO

 END SUBROUTINE map_surf_bg

#endif
 END SUBROUTINE radiation_tenstream


!------------------------------------------------------------------------------!
! Description:
! ------------
!> calculate the sundirection in favor of TenStream model
!------------------------------------------------------------------------------!
#if defined( __tenstream )
 SUBROUTINE radiation_calc_sundir( sundir )

    USE control_parameters,                                                                        &
        ONLY:  rotation_angle

    IMPLICIT NONE

    REAL(IREALS), INTENT(OUT) ::  sundir(3)

    REAL(wp) ::  second_of_day  !< second of the day
    REAL(wp) ::  solar_azim     !< solar azimuth in rotated model coordinates
    REAL(wp) ::  zenith         !< solar zenith  in rotated model coordinates


!
!-- Calculate current zenith and azimuth angle, sun direction, and whether the sun is up
    sun_direction = .TRUE.
    CALL get_date_time( time_since_reference_point, second_of_day = second_of_day,                 &
                        day_of_year = day_of_year )

    CALL calc_zenith( day_of_year, second_of_day )
!
!-- To avoid numerical instability near horizon, we use a minimum value for cos_zenith.
!-- zenith and azimuth in rad.
    zenith = ACOS( MAX( min_stable_coszen, cos_zenith ) )
    solar_azim = ATAN2( sun_dir_lon, sun_dir_lat ) - rotation_angle * ( pi / 180.0_wp )

    sundir(1) = REAL( -SIN( zenith ) * SIN( solar_azim ), IREALS )
    sundir(2) = REAL( -SIN( zenith ) * COS( solar_azim ), IREALS )
    sundir(3) = REAL( -COS( zenith ), IREALS )

    sundir = sundir / NORM2( sundir ) * REAL( solar_constant, IREALS )
!
!-- Calculate Earth/Sun distance from DYOFYR, the cumulative day of the year.
    IF ( day_of_year > 0 )  THEN
      sundir = sundir * REAL( ts_earth_sun( day_of_year ), IREALS )
    ENDIf

 END SUBROUTINE radiation_calc_sundir
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Use external radiative forcing (short- and longwave downwelling radiation) from a driver
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_external

    IMPLICIT NONE

    INTEGER(iwp) ::  t   !< index of current timestep
    INTEGER(iwp) ::  tm  !< index of previous timestep

    REAL(wp) ::  fac_dt              !< interpolation factor
    REAL(wp) ::  second_of_day_init  !< second of the day at model start

    TYPE(surf_type), POINTER ::  surf  !< pointer on respective surface type, used to generalize routine

!
!-- Calculate current zenith angle
    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year,                     &
                        second_of_day = second_of_day )
    CALL calc_zenith( day_of_year, second_of_day )
!
!-- Interpolate external radiation on current timestep
    IF ( time_since_reference_point  <= 0.0_wp )  THEN
       t      = 0
       tm     = 0
       fac_dt = 0
    ELSE
       CALL get_date_time( 0.0_wp, second_of_day=second_of_day_init )
       t = 0
       DO WHILE ( time_rad_f%var1d(t) <= time_since_reference_point )
          t = t + 1
       ENDDO

       tm = MAX( t-1, 0 )

       fac_dt = ( time_since_reference_point - time_rad_f%var1d(tm) + dt_3d ) /                    &
                MAX( TINY( 1.0_wp ), ( time_rad_f%var1d(t)  - time_rad_f%var1d(tm) ) )
       fac_dt = MIN( 1.0_wp, fac_dt )
    ENDIF
!
!-- Call clear-sky calculation for each surface orientation.
    surf => surf_lsm
    CALL radiation_external_surf
    surf => surf_usm
    CALL radiation_external_surf

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Todo: Subroutine description missing!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_external_surf

    USE control_parameters

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index along x-dimension
    INTEGER(iwp) ::  j  !< grid index along y-dimension
    INTEGER(iwp) ::  k  !< grid index along z-dimension
    INTEGER(iwp) ::  m  !< running index for surface elements

    REAL(wp) ::  lw_in      !< downwelling longwave radiation, interpolated value
    REAL(wp) ::  sw_in      !< downwelling shortwave radiation, interpolated value
    REAL(wp) ::  sw_in_dif  !< downwelling diffuse shortwave radiation, interpolated value

    REAL(wp), DIMENSION(1:7) ::  combine_allreduce    !< dummy array used to combine several MPI_ALLREDUCE calls
    REAL(wp), DIMENSION(1:7) ::  combine_allreduce_l  !< dummy array used to combine several MPI_ALLREDUCE calls

    IF ( surf%ns < 1 )  RETURN
!
!-- Level-of-detail = 1. Note, here it must be distinguished between averaged radiation and
!-- non-averaged radiation for the upwelling fluxes.
    IF ( rad_sw_in_f%lod == 1 )  THEN

       sw_in = ( 1.0_wp - fac_dt ) * rad_sw_in_f%var1d(tm) + fac_dt * rad_sw_in_f%var1d(t)

       lw_in = ( 1.0_wp - fac_dt ) * rad_lw_in_f%var1d(tm) + fac_dt * rad_lw_in_f%var1d(t)
!
!--    Limit shortwave incoming radiation to positive values, in order to overcome possible
!--    observation errors.
       sw_in = MAX( 0.0_wp, sw_in )
       sw_in = MERGE( sw_in, 0.0_wp, sun_up )

       surf%rad_sw_in = sw_in
       surf%rad_lw_in = lw_in

       IF ( average_radiation )  THEN
          IF ( dcep )  THEN

             combine_allreduce_l = 0.0_wp

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  m = surf%start_index(j,i), surf%end_index(j,i)
!
!--                   Albedo
                      combine_allreduce_l(1) =                                                     &
                      combine_allreduce_l(1) + SUM( surf%frac(m,:) * surf%albedo(m,:) )            &
                                     * ( 1.0_wp - fr_urb(j,i) ) + albedop_dcep(j,i) * fr_urb(j,i)
!
!--                   Emissivity
                      combine_allreduce_l(2) =                                                     &
                      combine_allreduce_l(2) + SUM( surf%frac(m,:) * surf%emissivity(m,:) )        &
                                     * ( 1.0_wp - fr_urb(j,i) ) +  emiss_dcep(j,i) * fr_urb(j,i)
!
!--                   Flux
                      combine_allreduce_l(3) = &
                      combine_allreduce_l(3) + SUM( surf%frac(m,:) * surf%emissivity(m,:) ) *      &
                      ( surf%pt_surface(m) * exner(nzb) )**4 * ( 1.0_wp - fr_urb(j,i) ) +          &
                      fr_urb(j,i) * emiss_dcep(j,i) * t_grad_dcep(j,i)**4
                   ENDDO
                ENDDO
             ENDDO

             combine_allreduce_l(4) = REAL( surf%ns, KIND = wp )

#if defined( __parallel )
             CALL MPI_ALLREDUCE( combine_allreduce_l, combine_allreduce, SIZE( combine_allreduce ),&
                                 MPI_REAL, MPI_SUM, comm2d, ierr )
#else
             combine_allreduce = combine_allreduce_l
#endif

             albedo_eff     = combine_allreduce(1) / combine_allreduce(4)
             emissivity_eff = combine_allreduce(2) / combine_allreduce(4)
             t_rad_eff = ( combine_allreduce(3) / combine_allreduce(4) / emissivity_eff )**0.25_wp

          ENDIF ! dcep

          surf%rad_sw_out = albedo_eff * surf%rad_sw_in

          surf%rad_lw_out = emissivity_eff * sigma_sb * t_rad_eff**4 +                             &
                            ( 1.0_wp - emissivity_eff ) * surf%rad_lw_in

          surf%rad_net = surf%rad_sw_in - surf%rad_sw_out + surf%rad_lw_in - surf%rad_lw_out

          surf%rad_lw_out_change_0 = 4.0_wp * emissivity_eff * sigma_sb * t_rad_eff**3

       ELSE

          DO  m = 1, surf%ns
             k = surf%k(m)
             surf%rad_sw_out(m) = ( surf%frac(m,ind_veg_wall)  * surf%albedo(m,ind_veg_wall)       &
                                  + surf%frac(m,ind_pav_green) * surf%albedo(m,ind_pav_green)      &
                                  + surf%frac(m,ind_wat_win)   * surf%albedo(m,ind_wat_win) )      &
                                  * surf%rad_sw_in(m)

             surf%rad_lw_out(m) = ( surf%frac(m,ind_veg_wall)  * surf%emissivity(m,ind_veg_wall)   &
                                  + surf%frac(m,ind_pav_green) * surf%emissivity(m,ind_pav_green)  &
                                  + surf%frac(m,ind_wat_win)   * surf%emissivity(m,ind_wat_win) )  &
                                  * sigma_sb * ( surf%pt_surface(m) * exner(k) )**4

             surf%rad_lw_out_change_0(m) = ( surf%frac(m,ind_veg_wall) *                           &
                                             surf%emissivity(m,ind_veg_wall)                       &
                                           + surf%frac(m,ind_pav_green) *                          &
                                             surf%emissivity(m,ind_pav_green)                      &
                                           + surf%frac(m,ind_wat_win)   *                          &
                                             surf%emissivity(m,ind_wat_win) ) * 4.0_wp * sigma_sb  &
                                           * ( surf%pt_surface(m) * exner(k) )**3
          ENDDO

       ENDIF
!
!--    If diffuse shortwave radiation is available, store it on the respective files.
       IF ( rad_sw_in_dif_f%from_file )  THEN
          sw_in_dif= ( 1.0_wp - fac_dt ) * rad_sw_in_dif_f%var1d(tm)                               &
                     + fac_dt * rad_sw_in_dif_f%var1d(t)

          IF ( ALLOCATED( rad_sw_in_diff ) )  rad_sw_in_diff = sw_in_dif
          IF ( ALLOCATED( rad_sw_in_dir  ) )  rad_sw_in_dir  = sw_in - sw_in_dif
!
!--       Diffuse longwave radiation equals the total downwelling longwave radiation
          IF ( ALLOCATED( rad_lw_in_diff ) )  rad_lw_in_diff = lw_in
       ENDIF
!
!-- level-of-detail = 2
    ELSE

       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)

          surf%rad_sw_in(m) = ( 1.0_wp - fac_dt ) * rad_sw_in_f%var3d(tm,j,i)                      &
                              + fac_dt * rad_sw_in_f%var3d(t,j,i)
!
!--       Limit shortwave incoming radiation to positive values, in order to overcome possible
!--       observation errors.
          surf%rad_sw_in(m) = MAX( 0.0_wp, surf%rad_sw_in(m) )
          surf%rad_sw_in(m) = MERGE( surf%rad_sw_in(m), 0.0_wp, sun_up )

          surf%rad_lw_in(m) = ( 1.0_wp - fac_dt ) * rad_lw_in_f%var3d(tm,j,i)                      &
                              + fac_dt * rad_lw_in_f%var3d(t,j,i)
!
!--       Weighted average according to surface fraction.
          surf%rad_sw_out(m) = ( surf%frac(m,ind_veg_wall)  * surf%albedo(m,ind_veg_wall)          &
                               + surf%frac(m,ind_pav_green) * surf%albedo(m,ind_pav_green)         &
                               + surf%frac(m,ind_wat_win)   * surf%albedo(m,ind_wat_win) )         &
                               * surf%rad_sw_in(m)

          surf%rad_lw_out(m) = ( surf%frac(m,ind_veg_wall)  * surf%emissivity(m,ind_veg_wall)      &
                               + surf%frac(m,ind_pav_green) * surf%emissivity(m,ind_pav_green)     &
                               + surf%frac(m,ind_wat_win)   * surf%emissivity(m,ind_wat_win) )     &
                               * sigma_sb * ( surf%pt_surface(m) * exner(k) )**4

          surf%rad_lw_out_change_0(m) = ( surf%frac(m,ind_veg_wall) *                              &
                                          surf%emissivity(m,ind_veg_wall)                          &
                                        + surf%frac(m,ind_pav_green) *                             &
                                          surf%emissivity(m,ind_pav_green)                         &
                                        + surf%frac(m,ind_wat_win) *                               &
                                          surf%emissivity(m,ind_wat_win) ) * 4.0_wp * sigma_sb     &
                                        * ( surf%pt_surface(m) * exner(k) )**3

          surf%rad_net(m) = surf%rad_sw_in(m) - surf%rad_sw_out(m) + surf%rad_lw_in(m) -           &
                            surf%rad_lw_out(m)
!
!--       If diffuse shortwave radiation is available, store it on the respective files.
          IF ( rad_sw_in_dif_f%from_file )  THEN
             IF ( ALLOCATED( rad_sw_in_diff ) )                                                    &
                rad_sw_in_diff(j,i) = ( 1.0_wp - fac_dt ) * rad_sw_in_dif_f%var3d(tm,j,i)          &
                                      + fac_dt * rad_sw_in_dif_f%var3d(t,j,i)
!
!--          dir = sw_in - sw_in_dif.
             IF ( ALLOCATED( rad_sw_in_dir  ) )                                                    &
                rad_sw_in_dir(j,i)  = surf%rad_sw_in(m) - rad_sw_in_diff(j,i)
!
!--          Diffuse longwave radiation equals the total downwelling longwave radiation
             IF ( ALLOCATED( rad_lw_in_diff ) )  rad_lw_in_diff(j,i) = surf%rad_lw_in(m)
          ENDIF

       ENDDO

    ENDIF
!
!-- Store radiation also on 2D arrays, which are still used for direct-diffuse splitting. Note,
!-- this is only required for horizontal surfaces, which cover all (x,y)-position.
    DO  m = 1, surf%ns
       IF ( surf%upward(m) )  THEN
          i = surf%i(m)
          j = surf%j(m)

          rad_sw_in(0,j,i)  = surf%rad_sw_in(m)
          rad_lw_in(0,j,i)  = surf%rad_lw_in(m)
          rad_sw_out(0,j,i) = surf%rad_sw_out(m)
          rad_lw_out(0,j,i) = surf%rad_lw_out(m)
       ENDIF
    ENDDO

 END SUBROUTINE radiation_external_surf

 END SUBROUTINE radiation_external

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> A simple clear sky radiation model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_clearsky

    IMPLICIT NONE

    REAL(wp) ::  pt1    !< potential temperature at first grid level or mean value at urban layer top
    REAL(wp) ::  pt1_l  !< potential temperature at first grid level or mean value at urban layer top at local subdomain
    REAL(wp) ::  ql1    !< liquid water mixing ratio at first grid level or mean value at urban layer top
    REAL(wp) ::  ql1_l  !< liquid water mixing ratio at first grid level or mean value at urban layer top at local subdomain

    TYPE(surf_type), POINTER ::  surf  !< pointer on respective surface type, used to generalize routine

!
!-- Calculate current zenith angle
    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year,                     &
                        second_of_day = second_of_day )
    CALL calc_zenith( day_of_year, second_of_day )

!
!-- Calculate sky transmissivity
    sky_trans = 0.6_wp + 0.2_wp * cos_zenith

!
!-- Calculate value of the Exner function at model surface
!
!-- In case averaged radiation is used, calculate mean temperature and liquid water mixing ratio at
!-- the urban-layer top.
    IF ( average_radiation )  THEN
       pt1 = 0.0_wp
       IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1 = 0.0_wp

       pt1_l = SUM( pt(nz_urban_t,nys:nyn,nxl:nxr) )
       IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1_l = SUM( ql(nz_urban_t,nys:nyn,nxl:nxr) )

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( pt1_l, pt1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )

       IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
           CALL MPI_ALLREDUCE( ql1_l, ql1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
       ENDIF
#else
       pt1 = pt1_l
       IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1 = ql1_l
#endif

       IF ( bulk_cloud_model  .OR.  cloud_droplets )  pt1 = pt1 + lv_d_cp / exner(nz_urban_t) * ql1
!
!--    Finally, divide by number of grid points
       pt1 = pt1 / REAL( ( nx + 1 ) * ( ny + 1 ), KIND = wp )
    ENDIF
!
!-- Call clear-sky calculation for each surface orientation.
    surf => surf_lsm
    CALL radiation_clearsky_surf
    surf => surf_usm
    CALL radiation_clearsky_surf

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Todo: Subroutine description missing.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_clearsky_surf

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< index x-direction
    INTEGER(iwp) ::  j  !< index y-direction
    INTEGER(iwp) ::  k  !< index z-direction
    INTEGER(iwp) ::  m  !< running index for surface elements

    REAL(wp), DIMENSION(1:7) ::  combine_allreduce    !< dummy array used to combine several MPI_ALLREDUCE calls
    REAL(wp), DIMENSION(1:7) ::  combine_allreduce_l  !< dummy array used to combine several MPI_ALLREDUCE calls


    IF ( surf%ns < 1 )  RETURN

!
!-- Calculate radiation fluxes and net radiation (rad_net) assuming homogeneous urban radiation
!-- conditions.
    IF ( average_radiation )  THEN

       IF ( dcep )  THEN

          combine_allreduce_l = 0.0_wp

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  m = surf%start_index(j,i), surf%end_index(j,i)
!
!--                Albedo
                   combine_allreduce_l(1) =                                                        &
                   combine_allreduce_l(1) + SUM( surf%frac(m,:) * surf%rrtm_asdir(m,:) )           &
                                  * ( 1.0_wp - fr_urb(j,i) ) + albedop_dcep(j,i) * fr_urb(j,i)
                   combine_allreduce_l(2) =                                                        &
                   combine_allreduce_l(2) + SUM( surf%frac(m,:) * surf%rrtm_asdif(m,:) )           &
                                  * ( 1.0_wp - fr_urb(j,i) ) + albedop_dcep(j,i) * fr_urb(j,i)
                   combine_allreduce_l(3) =                                                        &
                   combine_allreduce_l(3) + SUM( surf%frac(m,:) * surf%rrtm_aldir(m,:) )           &
                                  * ( 1.0_wp - fr_urb(j,i) ) + albedop_dcep(j,i) * fr_urb(j,i)
                   combine_allreduce_l(4) =                                                        &
                   combine_allreduce_l(4) + SUM( surf%frac(m,:) * surf%rrtm_aldif(m,:) )           &
                                  * ( 1.0_wp - fr_urb(j,i) ) + albedop_dcep(j,i) * fr_urb(j,i)
!
!--                Emissivity
                   combine_allreduce_l(2) =                                                        &
                   combine_allreduce_l(2) + SUM( surf%frac(m,:) * surf%emissivity(m,:) )           &
                                  * ( 1.0_wp - fr_urb(j,i) ) +  emiss_dcep(j,i) * fr_urb(j,i)
!
!--                Flux
                   combine_allreduce_l(6) =                                                        &
                   combine_allreduce_l(6) + SUM( surf_lsm%frac(m,:) * surf_lsm%emissivity(m,:) )   &
                         * ( surf_lsm%pt_surface(m) * exner(nzb) )**4 * ( 1.0_wp - fr_urb(j,i) )   &
                         + fr_urb(j,i) * emiss_dcep(j,i) * t_grad_dcep(j,i)**4
                ENDDO
             ENDDO
          ENDDO

          combine_allreduce_l(7) = REAL( surf_lsm%ns, KIND = wp )

#if defined( __parallel )
          CALL MPI_ALLREDUCE( combine_allreduce_l, combine_allreduce, SIZE( combine_allreduce ),   &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
#else
          combine_allreduce = combine_allreduce_l
#endif
!
!--       For rrtmg, we do not differentiate between asdir, asdif, aldir, and aldif.
!--       Here we use only aldif.
          albedo_eff     = combine_allreduce(4) / combine_allreduce(7)
          emissivity_eff = combine_allreduce(5) / combine_allreduce(7)
          t_rad_eff = ( combine_allreduce(6) / combine_allreduce(7) / emissivity_eff )**0.25_wp

       ENDIF ! dcep

       k = nz_urban_t

       surf%rad_sw_in  = solar_constant * sky_trans * cos_zenith
       surf%rad_sw_out = albedo_eff * surf%rad_sw_in

       surf%rad_lw_in  = emissivity_atm_clsky * sigma_sb * ( pt1 * exner(k+1) )**4

       surf%rad_lw_out = emissivity_eff * sigma_sb * ( t_rad_eff )**4 *                            &
                         ( 1.0_wp - emissivity_eff ) * surf%rad_lw_in

       surf%rad_net = surf%rad_sw_in - surf%rad_sw_out + surf%rad_lw_in - surf%rad_lw_out

       surf%rad_lw_out_change_0 = 4.0_wp * emissivity_eff * sigma_sb * ( t_rad_eff )**3

!
!-- Calculate radiation fluxes and net radiation (rad_net) for each surface element.
    ELSE

       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)

          surf%rad_sw_in(m) = solar_constant * sky_trans * cos_zenith

!
!--       Weighted average according to surface fraction.
!--       ATTENTION: when radiation interactions are switched on the calculated fluxes below are not
!--       actually used as they are overwritten in radiation_interaction.
          surf%rad_sw_out(m) = ( surf%frac(m,ind_veg_wall)  * surf%albedo(m,ind_veg_wall)          &
                               + surf%frac(m,ind_pav_green) * surf%albedo(m,ind_pav_green)         &
                               + surf%frac(m,ind_wat_win)   * surf%albedo(m,ind_wat_win) )         &
                               * surf%rad_sw_in(m)

          surf%rad_lw_out(m) = ( surf%frac(m,ind_veg_wall)  * surf%emissivity(m,ind_veg_wall)      &
                               + surf%frac(m,ind_pav_green) * surf%emissivity(m,ind_pav_green)     &
                               + surf%frac(m,ind_wat_win)   * surf%emissivity(m,ind_wat_win) )     &
                               * sigma_sb * ( surf%pt_surface(m) * exner(nzb) )**4

          surf%rad_lw_out_change_0(m) = ( surf%frac(m,ind_veg_wall)  *                             &
                                          surf%emissivity(m,ind_veg_wall)                          &
                                        + surf%frac(m,ind_pav_green) *                             &
                                          surf%emissivity(m,ind_pav_green)                         &
                                        + surf%frac(m,ind_wat_win)   *                             &
                                          surf%emissivity(m,ind_wat_win) ) * 4.0_wp * sigma_sb     &
                                        * ( surf%pt_surface(m) * exner(nzb) )** 3


          IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
             pt1 = pt(k,j,i) + lv_d_cp / exner(k) * ql(k,j,i)
             surf%rad_lw_in(m) = emissivity_atm_clsky * sigma_sb * ( pt1 * exner(k) )**4
          ELSE
             surf%rad_lw_in(m) = emissivity_atm_clsky * sigma_sb * ( pt(k,j,i) * exner(k) )**4
          ENDIF

          surf%rad_net(m) = surf%rad_sw_in(m) - surf%rad_sw_out(m) + surf%rad_lw_in(m) -           &
                            surf%rad_lw_out(m)

       ENDDO

    ENDIF

!
!-- Fill out values in radiation arrays. Note, this is only required for horizontal surfaces, which
!-- covers all x,y position.
    DO  m = 1, surf%ns
       IF ( surf%upward(m) )  THEN
          i = surf%i(m)
          j = surf%j(m)
          rad_sw_in(0,j,i)  = surf%rad_sw_in(m)
          rad_sw_out(0,j,i) = surf%rad_sw_out(m)
          rad_lw_in(0,j,i)  = surf%rad_lw_in(m)
          rad_lw_out(0,j,i) = surf%rad_lw_out(m)
       ENDIF
    ENDDO

 END SUBROUTINE radiation_clearsky_surf

 END SUBROUTINE radiation_clearsky


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This scheme keeps the prescribed net radiation constant during the run
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_constant


    IMPLICIT NONE

    REAL(wp) ::  pt1    !< potential temperature at first grid level or mean value at urban layer top
    REAL(wp) ::  pt1_l  !< potential temperature at first grid level or mean value at urban layer top at local subdomain
    REAL(wp) ::  ql1    !< liquid water mixing ratio at first grid level or mean value at urban layer top
    REAL(wp) ::  ql1_l  !< liquid water mixing ratio at first grid level or mean value at urban layer top at local subdomain

    TYPE(surf_type), POINTER ::  surf  !< pointer on respective surface type, used to generalize routine

!
!-- In case averaged radiation is used, calculate mean temperature and liquid water mixing ratio at
!-- the urban-layer top.
    IF ( average_radiation )  THEN
       pt1 = 0.0_wp
       IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1 = 0.0_wp

       pt1_l = SUM( pt(nz_urban_t,nys:nyn,nxl:nxr) )
       IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1_l = SUM( ql(nz_urban_t,nys:nyn,nxl:nxr) )

#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( pt1_l, pt1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
       IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
          CALL MPI_ALLREDUCE( ql1_l, ql1, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
       ENDIF
#else
       pt1 = pt1_l
       IF ( bulk_cloud_model  .OR.  cloud_droplets )  ql1 = ql1_l
#endif
       IF ( bulk_cloud_model  .OR.  cloud_droplets )  pt1 = pt1 + lv_d_cp / exner(nz_urban_t+1) *  &
                                                            ql1
!
!--    Finally, divide by number of grid points
       pt1 = pt1 / REAL( ( nx + 1 ) * ( ny + 1 ), KIND = wp )
    ENDIF

!
!-- First, natural surfaces, then building surfaces
    surf => surf_lsm
    CALL radiation_constant_surf
    surf => surf_usm
    CALL radiation_constant_surf

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Todo: Subroutine description missing!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_constant_surf

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !< index x-direction
    INTEGER(iwp) ::  j     !< index y-direction
    INTEGER(iwp) ::  k     !< index z-direction
    INTEGER(iwp) ::  m     !< running index for surface elements

    IF ( surf%ns < 1 )  RETURN

!-- Calculate homogenoeus urban radiation fluxes
    IF ( average_radiation )  THEN

       surf%rad_net = net_radiation

       surf%rad_lw_in = emissivity_atm_clsky * sigma_sb * ( pt1 * exner(nz_urban_t+1) )**4

       surf%rad_lw_out = emissivity_eff * sigma_sb * t_rad_eff**4                                  &
                         + ( 1.0_wp - emissivity_eff ) * surf%rad_lw_in

       surf%rad_lw_out_change_0 = 4.0_wp * emissivity_eff * sigma_sb * t_rad_eff**3

       surf%rad_sw_in = ( surf%rad_net - surf%rad_lw_in + surf%rad_lw_out ) /                      &
                        ( 1.0_wp - albedo_eff )

       surf%rad_sw_out = albedo_eff * surf%rad_sw_in

!
!-- Calculate radiation fluxes for each surface element
    ELSE
!
!--    Prescribe net radiation and estimate the remaining radiative fluxes
       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)

          surf%rad_net(m) = net_radiation

          IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
             pt1 = pt(k,j,i) + lv_d_cp / exner(k) * ql(k,j,i)
             surf%rad_lw_in(m) = emissivity_atm_clsky * sigma_sb * ( pt1 * exner(k) )**4
          ELSE
             surf%rad_lw_in(m) = emissivity_atm_clsky * sigma_sb * ( pt(k,j,i) * exner(k) )**4
          ENDIF

!
!--       Weighted average according to surface fraction.
          surf%rad_lw_out(m) = ( surf%frac(m,ind_veg_wall)  * surf%emissivity(m,ind_veg_wall)      &
                               + surf%frac(m,ind_pav_green) * surf%emissivity(m,ind_pav_green)     &
                               + surf%frac(m,ind_wat_win)   * surf%emissivity(m,ind_wat_win) )     &
                               * sigma_sb * ( surf%pt_surface(m) * exner(nzb) )**4

          surf%rad_sw_in(m) = ( surf%rad_net(m) - surf%rad_lw_in(m) + surf%rad_lw_out(m) )         &
                              / ( 1.0_wp - ( surf%frac(m,ind_veg_wall) *                           &
                                             surf%albedo(m,ind_veg_wall)                           &
                                           + surf%frac(m,ind_pav_green) *                          &
                                             surf%albedo(m,ind_pav_green)                          &
                                           + surf%frac(m,ind_wat_win) *                            &
                                             surf%albedo(m,ind_wat_win)                            &
                                           )                                                       &
                                 )

          surf%rad_sw_out(m) = ( surf%frac(m,ind_veg_wall)  * surf%albedo(m,ind_veg_wall)          &
                               + surf%frac(m,ind_pav_green) * surf%albedo(m,ind_pav_green)         &
                               + surf%frac(m,ind_wat_win)   * surf%albedo(m,ind_wat_win) )         &
                               * surf%rad_sw_in(m)

       ENDDO

    ENDIF

!
!-- Fill out values in radiation arrays. Note, this is only required for horizontal surfaces, which
!-- covers all x,y position.
    DO  m = 1, surf%ns
       IF ( surf%upward(m) )  THEN
          i = surf%i(m)
          j = surf%j(m)
          rad_sw_in(0,j,i)  = surf%rad_sw_in(m)
          rad_sw_out(0,j,i) = surf%rad_sw_out(m)
          rad_lw_in(0,j,i)  = surf%rad_lw_in(m)
          rad_lw_out(0,j,i) = surf%rad_lw_out(m)
       ENDIF
    ENDDO

 END SUBROUTINE radiation_constant_surf


 END SUBROUTINE radiation_constant

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for radiation model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_header( io )


    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  io  !< Unit of the output file


!
!-- Write radiation model header.
    WRITE( io, 3 )

    IF ( radiation_scheme == 'constant' )  THEN
       WRITE( io, 4 ) net_radiation
    ELSEIF ( radiation_scheme == 'clear-sky' )  THEN
       WRITE( io, 5 )
    ELSEIF ( radiation_scheme == 'rrtmg' )  THEN
       WRITE( io, 6 )
       IF ( .NOT. lw_radiation )  WRITE( io, 10 )
       IF ( .NOT. sw_radiation )  WRITE( io, 11 )
    ELSEIF ( radiation_scheme == 'external' )  THEN
       WRITE( io, 14 )
    ENDIF

    IF ( albedo_type_f%from_file    .OR.  vegetation_type_f%from_file  .OR.                        &
         pavement_type_f%from_file  .OR.  water_type_f%from_file       .OR.                        &
         building_type_f%from_file )  THEN
          WRITE( io, 13 )
    ELSE
       IF ( albedo_type == 0 )  THEN
          WRITE( io, 7 ) albedo
       ELSE
          WRITE( io, 8 ) TRIM( albedo_type_name(albedo_type) )
       ENDIF
    ENDIF
    IF ( constant_albedo )  THEN
       WRITE( io, 9 )
    ENDIF

    WRITE( io, 12 ) dt_radiation, ndsidir, CEILING( ( end_time - spinup_time ) / dt_radiation ) + 1



 3 FORMAT ( //' Radiation model information:'/ ' ----------------------------'/ )
 4 FORMAT ( '    --> Using constant net radiation: net_radiation = ', F6.2, // 'W/m**2' )
 5 FORMAT ( '    --> Simple radiation scheme for clear sky is used (no clouds,', ' default)' )
 6 FORMAT ( '    --> RRTMG scheme is used' )
 7 FORMAT ( /'    User-specific surface albedo: albedo =', F6.3 )
 8 FORMAT ( /'    Albedo is set for land surface type: ', A )
 9 FORMAT ( /'    --> Albedo is fixed during the run' )
10 FORMAT ( /'    --> Longwave radiation is disabled' )
11 FORMAT ( /'    --> Shortwave radiation is disabled.' )
12 FORMAT ( '    Timestep: dt_radiation = ', F6.2, '  s' /                                         &
            '    Precalculated ', I4, ' solar positions from ', I5, ' timesteps')
13 FORMAT ( /'    Albedo is set individually for each xy-location, according ',                    &
                 'to given surface type.')
14 FORMAT ( '    --> External radiation forcing is used' )


 END SUBROUTINE radiation_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &radiation_parameters for radiation model and RTM
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_parin

    IMPLICIT NONE

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /radiation_parameters/ albedo,                                                        &
                                    albedo_lw_dif,                                                 &
                                    albedo_lw_dir,                                                 &
                                    albedo_sw_dif,                                                 &
                                    albedo_sw_dir,                                                 &
                                    albedo_type,                                                   &
                                    bufsize_alltoall,                                              &
                                    constant_albedo,                                               &
                                    dt_radiation,                                                  &
                                    emissivity,                                                    &
                                    lw_radiation,                                                  &
                                    localized_raytracing,                                          &
                                    mrt_geom,                                                      &
                                    mrt_geom_params,                                               &
                                    mrt_include_sw,                                                &
                                    mrt_minlevel,                                                  &
                                    mrt_nlevels,                                                   &
                                    mrt_skip_roof,                                                 &
                                    net_radiation,                                                 &
                                    nrefsteps,                                                     &
                                    plant_lw_interact,                                             &
                                    radiation_interactions_on,                                     &
                                    radiation_only,                                                &
                                    radiation_scheme,                                              &
                                    radiation_volumetric_flux,                                     &
                                    raytrace_discrete_azims,                                       &
                                    raytrace_discrete_elevs,                                       &
                                    raytrace_mpi_rma,                                              &
                                    trace_fluxes_above,                                            &
                                    skip_time_do_radiation,                                        &
                                    surface_reflections,                                           &
                                    switch_off_module,                                             &
                                    sw_radiation,                                                  &
#if defined( __tenstream )
                                    ts_icollapse,                                                  &
#endif
#if defined( __rrtmg ) || defined( __tenstream )
                                    use_broadband_albedo,                                          &
#endif
                                    unscheduled_radiation_calls

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, radiation_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    radiation_parameters namelist was found and read correctly. Set flag that indicates that the
!--    radiation model is switched on.
       IF ( .NOT. switch_off_module )  radiation = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    radiation_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'radiation_parameters', line )

    ENDIF

 END SUBROUTINE radiation_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Implementation of the RRTMG radiation_scheme
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_rrtmg

#if defined( __rrtmg ) && defined( __netcdf )
    USE palm_date_time_mod,                                                                        &
        ONLY:  hours_per_day

    USE particle_attributes,                                                                       &
        ONLY:  grid_particles,                                                                     &
               number_of_particles,                                                                &
               particles,                                                                          &
               prt_count


    IMPLICIT NONE


    INTEGER(iwp) ::  i        !< grid index in x-direction
    INTEGER(iwp) ::  j        !< grid index in y-direction
    INTEGER(iwp) ::  k        !< grid index in z-direction
    INTEGER(iwp) ::  m        !< running index for surface elements
    INTEGER(iwp) ::  n        !< loop index
    INTEGER(iwp) ::  k_topo_l !< topography top index on subdomain
    INTEGER(iwp) ::  k_topo   !< topography top index global

    REAL(wp) ::  d_hours_day     !< 1 / hours-per-day
    REAL(wp) ::  mass_xi         !< mass of cloud ice
    REAL(wp) ::  nc_rad          !< number concentration of cloud droplets
    REAL(wp) ::  rrtm_emis_save  !< saved value of rrtm_emis value
    REAL(wp) ::  s_r2            !< weighted sum over all droplets with r^2
    REAL(wp) ::  s_r3            !< weighted sum over all droplets with r^3

    REAL(wp), PARAMETER ::  a5 = 83.8_wp   !< parameter for ice effective radius (Roeckner et al., 2003)
    REAL(wp), PARAMETER ::  b5 = 0.216_wp  !< parameter for ice effective radius (Roeckner et al., 2003)

    REAL(wp), DIMENSION(0:0) ::  zenith  !< to provide indexed array

    REAL(wp), DIMENSION(0:nzt+1) ::  pt_av  !<
    REAL(wp), DIMENSION(0:nzt+1) ::  q_av   !<
    REAL(wp), DIMENSION(0:nzt+1) ::  ql_av  !<
    REAL(wp), DIMENSION(0:nzt+1) ::  qi_av  !<

    REAL(wp), DIMENSION(1:7) ::  combine_allreduce    !< dummy array used to combine several MPI_ALLREDUCE calls
    REAL(wp), DIMENSION(1:7) ::  combine_allreduce_l  !< dummy array used to combine several MPI_ALLREDUCE calls

!
!-- Just dummy arguments
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: rrtm_lw_taucld_dum,  &  !<
                                               rrtm_lw_tauaer_dum,  &  !<
                                               rrtm_sw_taucld_dum,  &  !<
                                               rrtm_sw_ssacld_dum,  &  !<
                                               rrtm_sw_asmcld_dum,  &  !<
                                               rrtm_sw_fsfcld_dum,  &  !<
                                               rrtm_sw_tauaer_dum,  &  !<
                                               rrtm_sw_ssaaer_dum,  &  !<
                                               rrtm_sw_asmaer_dum,  &  !<
                                               rrtm_sw_ecaer_dum       !<

!
!-- Pre-calculate parameters
    d_hours_day = 1.0_wp / REAL( hours_per_day, KIND = wp )

!
!-- Calculate current (cosine of) zenith angle and whether the sun is up
    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year,                     &
                        second_of_day=second_of_day )
    CALL calc_zenith( day_of_year, second_of_day )
    zenith(0) = cos_zenith
!
!-- Calculate surface albedo. In case average radiation is applied, this is not required.
    IF ( .NOT. constant_albedo )  THEN
       IF ( radiation_only )  CALL calc_albedo( surf_def )
       CALL calc_albedo( surf_lsm )
       CALL calc_albedo( surf_usm )
    ENDIF

!
!-- Prepare input data for RRTMG.
!-- In case of large scale forcing with surface data, calculate new pressure profile. nzt_rad might
!-- be modified by these calls and all required arrays will then be re-allocated.
    IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
       CALL read_sounding_data
       CALL read_trace_gas_data
    ENDIF


    IF ( average_radiation )  THEN

       IF ( dcep )  THEN

          combine_allreduce_l = 0.0_wp

          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
!
!--                Albedo
                   combine_allreduce_l(1) =                                                        &
                   combine_allreduce_l(1) + SUM( surf_lsm%frac(m,:) * surf_lsm%rrtm_asdir(m,:) )   &
                                  * ( 1.0_wp - fr_urb(j,i) ) + albedop_dcep(j,i) * fr_urb(j,i)
                   combine_allreduce_l(2) =                                                        &
                   combine_allreduce_l(2) + SUM( surf_lsm%frac(m,:) * surf_lsm%rrtm_asdif(m,:) )   &
                                  * ( 1.0_wp - fr_urb(j,i) ) + albedop_dcep(j,i) * fr_urb(j,i)
                   combine_allreduce_l(3) =                                                        &
                   combine_allreduce_l(3) + SUM( surf_lsm%frac(m,:) * surf_lsm%rrtm_aldir(m,:) )   &
                                  * ( 1.0_wp - fr_urb(j,i) ) + albedop_dcep(j,i) * fr_urb(j,i)
                   combine_allreduce_l(4) =                                                        &
                   combine_allreduce_l(4) + SUM( surf_lsm%frac(m,:) * surf_lsm%rrtm_aldif(m,:) )   &
                                  * ( 1.0_wp - fr_urb(j,i) ) + albedop_dcep(j,i) * fr_urb(j,i)
!
!--                Emissivity
                   combine_allreduce_l(5) =                                                        &
                   combine_allreduce_l(5) + SUM( surf_lsm%frac(m,:) * surf_lsm%emissivity(m,:) )   &
                                  * ( 1.0_wp - fr_urb(j,i) ) + emiss_dcep(j,i) * fr_urb(j,i)
!
!--                Flux
                   combine_allreduce_l(6) =                                                        &
                   combine_allreduce_l(6) + SUM( surf_lsm%frac(m,:) * surf_lsm%emissivity(m,:) )   &
                         * ( surf_lsm%pt_surface(m) * exner(nzb) )**4 * ( 1.0_wp - fr_urb(j,i) )   &
                         + fr_urb(j,i) * emiss_dcep(j,i) * t_grad_dcep(j,i)**4
                ENDDO
             ENDDO
          ENDDO

          combine_allreduce_l(7) = REAL( surf_lsm%ns, KIND = wp )

#if defined( __parallel )
          CALL MPI_ALLREDUCE( combine_allreduce_l, combine_allreduce, SIZE( combine_allreduce ),   &
                              MPI_REAL, MPI_SUM, comm2d, ierr )
#else
          combine_allreduce = combine_allreduce_l
#endif
!
!--       For rrtmg, we do not differentiate between asdir, asdif, aldir, and aldif.
!--       Here we use only aldif.
          albedo_eff     = combine_allreduce(4) / combine_allreduce(7)
          emissivity_eff = combine_allreduce(5) / combine_allreduce(7)
          t_rad_eff = ( combine_allreduce(6) / combine_allreduce(7) / emissivity_eff )**0.25_wp

       ENDIF ! dcep
!
!--    Determine minimum topography top index.
       k_topo_l = MINVAL( topo_top_ind(nys:nyn,nxl:nxr,0) )
#if defined( __parallel )
       CALL MPI_ALLREDUCE( k_topo_l, k_topo, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
#else
       k_topo = k_topo_l
#endif

       rrtm_asdir(1) = albedo_eff
       rrtm_asdif(1) = albedo_eff
       rrtm_aldir(1) = albedo_eff
       rrtm_aldif(1) = albedo_eff

       rrtm_emis = emissivity_eff
!
!--    Calculate mean pt profile.
       CALL calc_mean_profile( pt, 4, .TRUE. )
       pt_av = hom(:,1,4,0)

       IF ( humidity )  THEN
          CALL calc_mean_profile( q, 41, .TRUE. )
          q_av  = hom(:,1,41,0)
       ENDIF
!
!--    Prepare profiles of temperature and H2O volume mixing ratio.
       rrtm_tlev(0,k_topo+1) = t_rad_eff

       IF ( bulk_cloud_model )  THEN

          CALL calc_mean_profile( ql, 54, .TRUE. )
          ql_av = hom(:,1,54,0)

          IF ( microphysics_ice_phase )  THEN
             CALL calc_mean_profile( qi, 125, .TRUE. )
             qi_av = hom(:,1,125,0)
          ENDIF

          DO  k = nzb+1, nzt+1
             rrtm_tlay(0,k) = pt_av(k) * ( (hyp(k) ) / 100000.0_wp )**0.286_wp + lv_d_cp * ql_av(k)
             rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * ( q_av(k) - ql_av(k) )
          ENDDO
       ELSE
          DO  k = nzb+1, nzt+1
             rrtm_tlay(0,k) = pt_av(k) * ( (hyp(k) ) / 100000.0_wp )**0.286_wp
          ENDDO

          IF ( humidity )  THEN
             DO  k = nzb+1, nzt+1
                rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * q_av(k)
             ENDDO
          ELSE
             rrtm_h2ovmr(0,nzb+1:nzt+1) = 0.0_wp
          ENDIF
       ENDIF

!
!--    Avoid temperature/humidity jumps at the top of the PALM domain by linear interpolation from
!--    nzt+2 to nzt+7. Jumps are induced by discrepancies between the values in the domain and
!--    those above that are prescribed in RRTMG.
       DO  k = nzt+2, nzt+7
          rrtm_tlay(0,k) = rrtm_tlay(0,nzt+1) + ( rrtm_tlay(0,nzt+8) - rrtm_tlay(0,nzt+1) ) /      &
                           ( rrtm_play(0,nzt+8) - rrtm_play(0,nzt+1) ) *                           &
                           ( rrtm_play(0,k) - rrtm_play(0,nzt+1) )

          rrtm_h2ovmr(0,k) = rrtm_h2ovmr(0,nzt+1) + ( rrtm_h2ovmr(0,nzt+8) -                       &
                                                      rrtm_h2ovmr(0,nzt+1) ) /                     &
                             ( rrtm_play(0,nzt+8) - rrtm_play(0,nzt+1) ) *                         &
                             ( rrtm_play(0,k) - rrtm_play(0,nzt+1) )

       ENDDO

!--    Linear interpolation to zw grid. Loop reaches one level further up due to the staggered grid
!--    in RRTMG.
       DO  k = k_topo+2, nzt+8
          rrtm_tlev(0,k) = rrtm_tlay(0,k-1) + ( rrtm_tlay(0,k) - rrtm_tlay(0,k-1) ) /              &
                           ( rrtm_play(0,k) - rrtm_play(0,k-1) ) *                                 &
                           ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
       ENDDO
!
!--    Calculate liquid water path and cloud fraction for each column.
!--    Note that LWP is required in g/m2 instead of kg/kg m.
       rrtm_cldfr  = 0.0_wp
       rrtm_reliq  = 0.0_wp
       rrtm_cliqwp = 0.0_wp
       rrtm_icld   = 0
       rrtm_reice  = 0.0_wp
       rrtm_cicewp = 0.0_wp

       IF ( bulk_cloud_model )  THEN
          DO  k = k_topo+1, nzt+1
             rrtm_cliqwp(0,k) =  ql_av(k) * 1000.0_wp * ( rrtm_plev(0,k) - rrtm_plev(0,k+1) ) *    &
                                 100.0_wp / g
!
!--          Avoid rrtmg cloud calculation for very small values.
             IF ( rrtm_cliqwp(0,k) < 1.0E-20_wp )  rrtm_cliqwp(0,k) = 0.0_wp

             IF ( rrtm_cliqwp(0,k) > 0.0_wp )  THEN
                rrtm_cldfr(0,k) = 1.0_wp
                IF ( rrtm_icld == 0 )  rrtm_icld = 1

!
!--             Calculate cloud droplet effective radius
                rrtm_reliq(0,k) = 1.0E6_wp * ( 3.0_wp * ql_av(k) * rho_surface /                   &
                                  ( 4.0_wp * pi * nc_const * rho_l ) )**0.33333333333333_wp        &
                                  * EXP( LOG( sigma_gc )**2 )
!
!--             Limit effective radius
                IF ( rrtm_reliq(0,k) > 0.0_wp )  THEN
                   rrtm_reliq(0,k) = MAX( rrtm_reliq(0,k), 2.5_wp )
                   rrtm_reliq(0,k) = MIN( rrtm_reliq(0,k), 60.0_wp )
                ENDIF
             ENDIF
!
!--          Calculation of effective radius of ice (made for RRTM) based on ECHAM5
!--          documentation (Roeckner et al, MPI report 349). In this scheme graupel and
!--          snow is neglected.
             IF ( microphysics_ice_phase )  THEN
!
!--             Calculate ice water path in g/m2
                rrtm_cicewp(0,k) =  qi_av(k) * 1000.0_wp *                                         &
                         ( rrtm_plev(0,k) - rrtm_plev(0,k+1) ) * 100.0_wp / g
!
!--             Avoid rrtmg cloud calculation for very small values.
                IF ( rrtm_cicewp(0,k) < 1.0E-20_wp )  rrtm_cicewp(0,k) = 0.0_wp

                IF ( rrtm_cicewp(0,k) > 0.0_wp )  THEN

                   rrtm_cldfr(0,k) = 1.0_wp
                   IF ( rrtm_icld == 0 )  rrtm_icld = 1
!
!--                Calculate mean mass of particle.
                   mass_xi = hyrho(k) * qi_av(k) * 1000.0_wp
                   rrtm_reice(0,k) = a5 * EXP( b5 * LOG(mass_xi) )
!
!--                Limit ice effective radius to allowed range for this parameterization
                   rrtm_reice(0,k) =  MAX( MIN( rrtm_reice(0,k) , 131.0_wp ), 5.0_wp )

                ENDIF

             ENDIF

          ENDDO
       ENDIF

!
!--    Set surface temperature
       rrtm_tsfc = t_rad_eff

       IF ( lw_radiation )  THEN
!
!--       Due to technical reasons, copy optical depth to dummy arguments which are allocated on the
!--       exact size as the rrtmg_lw is called. As one dimension is allocated with zero size,
!--       compiler complains that rank of the array does not match that of the assumed-shaped
!--       arguments in the RRTMG library. In order to avoid this, write to dummy arguments and
!--       pass the entire dummy array. Seems to be the only existing work-around.
          ALLOCATE( rrtm_lw_taucld_dum(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1) )
          ALLOCATE( rrtm_lw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1) )

          rrtm_lw_taucld_dum = rrtm_lw_taucld(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1)
          rrtm_lw_tauaer_dum = rrtm_lw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1)

          CALL rrtmg_lw( 1, nzt_rad-k_topo, rrtm_icld, rrtm_idrv, rrtm_play(:,k_topo+1:),          &
                         rrtm_plev(:,k_topo+1:), rrtm_tlay(:,k_topo+1:), rrtm_tlev(:,k_topo+1:),   &
                         rrtm_tsfc, rrtm_h2ovmr(:,k_topo+1:), rrtm_o3vmr(:,k_topo+1:),             &
                         rrtm_co2vmr(:,k_topo+1:), rrtm_ch4vmr(:,k_topo+1:),                       &
                         rrtm_n2ovmr(:,k_topo+1:), rrtm_o2vmr(:,k_topo+1:),                        &
                         rrtm_cfc11vmr(:,k_topo+1:), rrtm_cfc12vmr(:,k_topo+1:),                   &
                         rrtm_cfc22vmr(:,k_topo+1:), rrtm_ccl4vmr(:,k_topo+1:), rrtm_emis,         &
                         rrtm_inflglw, rrtm_iceflglw, rrtm_liqflglw, rrtm_cldfr(:,k_topo+1:),      &
                         rrtm_lw_taucld_dum, rrtm_cicewp(:,k_topo+1:), rrtm_cliqwp(:,k_topo+1:),   &
                         rrtm_reice(:,k_topo+1:), rrtm_reliq(:,k_topo+1:), rrtm_lw_tauaer_dum,     &
                         rrtm_lwuflx(:,k_topo:), rrtm_lwdflx(:,k_topo:), rrtm_lwhr(:,k_topo+1:),   &
                         rrtm_lwuflxc(:,k_topo:), rrtm_lwdflxc(:,k_topo:), rrtm_lwhrc(:,k_topo+1:),&
                         rrtm_lwuflx_dt(:,k_topo:), rrtm_lwuflxc_dt(:,k_topo:) )

          DEALLOCATE( rrtm_lw_taucld_dum )
          DEALLOCATE( rrtm_lw_tauaer_dum )
!
!--       Save fluxes
          DO  i = nxl, nxr
             DO  j = nys, nyn
                k_topo_l = topo_top_ind(j,i,0)
                DO  k = k_topo_l, nzt+1
                   rad_lw_in(k,j,i)  = rrtm_lwdflx(0,k-k_topo_l+k_topo)
                   rad_lw_out(k,j,i) = rrtm_lwuflx(0,k-k_topo_l+k_topo)
                ENDDO
             ENDDO
          ENDDO
          rad_lw_in_diff(:,:) = rrtm_lwdflx(0,k_topo)
!
!--       Save heating rates (convert from K/d to K/h).
!--       Further, even though an aggregated radiation is computed, map signle-column profiles on
!--       top of any topography, in order to obtain correct near surface radiation heating/cooling
!--       rates.
          DO  i = nxl, nxr
             DO  j = nys, nyn
                k_topo_l = topo_top_ind(j,i,0)
                DO  k = k_topo_l+1, nzt+1
                   rad_lw_hr(k,j,i)    = rrtm_lwhr(0,k-k_topo_l+k_topo)  * d_hours_day
                   rad_lw_cs_hr(k,j,i) = rrtm_lwhrc(0,k-k_topo_l+k_topo) * d_hours_day
                ENDDO
             ENDDO
          ENDDO
!
!--       Save radiation flux in case of dcep average radiation (only upward LSM surfaces).
          IF ( dcep_average_radiation )  THEN
             surf_lsm%rad_lw_in           = rrtm_lwdflx(0,k_topo)
             surf_lsm%rad_lw_out          = rrtm_lwuflx(0,k_topo)
             surf_lsm%rad_lw_out_change_0 = rrtm_lwuflx_dt(0,k_topo)
          ENDIF

       ENDIF

       IF ( sw_radiation .AND. sun_up )  THEN
!
!--       Due to technical reasons, copy optical depths and other to dummy arguments which are
!--       allocated on the exact size as the rrtmg_sw is called. As one dimesion is allocated with
!--       zero size, compiler complains that rank of the array does not match that of the
!--       assumed-shaped arguments in the RRTMG library. In order to avoid this, write to dummy
!--       arguments and pass the entire dummy array. Seems to be the only existing work-around.
          ALLOCATE( rrtm_sw_taucld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
          ALLOCATE( rrtm_sw_ssacld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
          ALLOCATE( rrtm_sw_asmcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
          ALLOCATE( rrtm_sw_fsfcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
          ALLOCATE( rrtm_sw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
          ALLOCATE( rrtm_sw_ssaaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
          ALLOCATE( rrtm_sw_asmaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
          ALLOCATE( rrtm_sw_ecaer_dum(0:0,k_topo+1:nzt_rad+1,1:naerec+1)  )

          rrtm_sw_taucld_dum = rrtm_sw_taucld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
          rrtm_sw_ssacld_dum = rrtm_sw_ssacld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
          rrtm_sw_asmcld_dum = rrtm_sw_asmcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
          rrtm_sw_fsfcld_dum = rrtm_sw_fsfcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
          rrtm_sw_tauaer_dum = rrtm_sw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
          rrtm_sw_ssaaer_dum = rrtm_sw_ssaaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
          rrtm_sw_asmaer_dum = rrtm_sw_asmaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
          rrtm_sw_ecaer_dum  = rrtm_sw_ecaer(0:0,k_topo+1:nzt_rad+1,1:naerec+1)

          CALL rrtmg_sw( 1, nzt_rad-k_topo, rrtm_icld, rrtm_iaer, rrtm_play(:,k_topo+1:nzt_rad+1), &
                         rrtm_plev(:,k_topo+1:nzt_rad+2), rrtm_tlay(:,k_topo+1:nzt_rad+1),         &
                         rrtm_tlev(:,k_topo+1:nzt_rad+2), rrtm_tsfc,                               &
                         rrtm_h2ovmr(:,k_topo+1:nzt_rad+1), rrtm_o3vmr(:,k_topo+1:nzt_rad+1),      &
                         rrtm_co2vmr(:,k_topo+1:nzt_rad+1), rrtm_ch4vmr(:,k_topo+1:nzt_rad+1),     &
                         rrtm_n2ovmr(:,k_topo+1:nzt_rad+1), rrtm_o2vmr(:,k_topo+1:nzt_rad+1),      &
                         rrtm_asdir, rrtm_asdif, rrtm_aldir, rrtm_aldif, zenith, 0.0_wp,           &
                         day_of_year, solar_constant, rrtm_inflgsw, rrtm_iceflgsw, rrtm_liqflgsw,  &
                         rrtm_cldfr(:,k_topo+1:nzt_rad+1), rrtm_sw_taucld_dum, rrtm_sw_ssacld_dum, &
                         rrtm_sw_asmcld_dum, rrtm_sw_fsfcld_dum, rrtm_cicewp(:,k_topo+1:nzt_rad+1),&
                         rrtm_cliqwp(:,k_topo+1:nzt_rad+1), rrtm_reice(:,k_topo+1:nzt_rad+1),      &
                         rrtm_reliq(:,k_topo+1:nzt_rad+1), rrtm_sw_tauaer_dum,                     &
                         rrtm_sw_ssaaer_dum, rrtm_sw_asmaer_dum, rrtm_sw_ecaer_dum,                &
                         rrtm_swuflx(:,k_topo:nzt_rad+1), rrtm_swdflx(:,k_topo:nzt_rad+1),         &
                         rrtm_swhr(:,k_topo+1:nzt_rad+1), rrtm_swuflxc(:,k_topo:nzt_rad+1),        &
                         rrtm_swdflxc(:,k_topo:nzt_rad+1), rrtm_swhrc(:,k_topo+1:nzt_rad+1),       &
                         rrtm_dirdflux(:,k_topo:nzt_rad+1), rrtm_difdflux(:,k_topo:nzt_rad+1) )

          DEALLOCATE( rrtm_sw_taucld_dum )
          DEALLOCATE( rrtm_sw_ssacld_dum )
          DEALLOCATE( rrtm_sw_asmcld_dum )
          DEALLOCATE( rrtm_sw_fsfcld_dum )
          DEALLOCATE( rrtm_sw_tauaer_dum )
          DEALLOCATE( rrtm_sw_ssaaer_dum )
          DEALLOCATE( rrtm_sw_asmaer_dum )
          DEALLOCATE( rrtm_sw_ecaer_dum )

!
!--       Save radiation fluxes for the entire depth of the model domain
          DO  i = nxl, nxr
             DO  j = nys, nyn
                k_topo_l = topo_top_ind(j,i,0)
                DO  k = k_topo_l, nzt+1
                   rad_sw_in(k,j,i)  = rrtm_swdflx(0,k-k_topo_l+k_topo)
                   rad_sw_out(k,j,i) = rrtm_swuflx(0,k-k_topo_l+k_topo)
                ENDDO
             ENDDO
          ENDDO
!--       Save direct and diffuse SW radiation at the surface (required by RTM)
          rad_sw_in_dir(:,:) = rrtm_dirdflux(0,k_topo)
          rad_sw_in_diff(:,:) = rrtm_difdflux(0,k_topo)

!
!--       Save heating rates (convert from K/d to K/s)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                k_topo_l = topo_top_ind(j,i,0)
                DO  k = k_topo_l+1, nzt+1
                   rad_sw_hr(k,j,i)    = rrtm_swhr(0,k-k_topo_l+k_topo)  * d_hours_day
                   rad_sw_cs_hr(k,j,i) = rrtm_swhrc(0,k-k_topo_l+k_topo) * d_hours_day
                ENDDO
             ENDDO
          ENDDO
!
!--       Save radiation flux in case of dcep average radiation.
          IF ( dcep_average_radiation )  THEN
!
!--          Only upward LSM surfaces are defined in case of dcep.
             surf_lsm%rad_sw_in  = rrtm_swdflx(0,k_topo)
             surf_lsm%rad_sw_out = rrtm_swuflx(0,k_topo)
          ENDIF
!
!--    Solar radiation is zero during night.
       ELSE
          rad_sw_in  = 0.0_wp
          rad_sw_out = 0.0_wp
          rad_sw_in_dir(:,:) = 0.0_wp
          rad_sw_in_diff(:,:) = 0.0_wp
!
!--       Save radiation flux in case of dcep average radiation.
          IF ( dcep_average_radiation )  THEN
             surf_lsm%rad_sw_in  = 0.0_wp
             surf_lsm%rad_sw_out = 0.0_wp
          ENDIF

       ENDIF
!
!-- RRTMG is called for each (j,i) grid point separately, starting at the highest topography level.
!-- Here no RTM is used since average_radiation is false.
    ELSE
!
!--    Loop over all grid points
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Obtain topography top index (lower bound of RRTMG).
             k_topo = topo_top_ind(j,i,0)
!
!--          Prepare profiles of temperature and H2O volume mixing ratio.
             IF ( radiation_only )  THEN
                DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                   rrtm_tlev(0,k_topo+1) = surf_def%pt_surface(m) * exner(k_topo)
                ENDDO
             ENDIF
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                rrtm_tlev(0,k_topo+1) = surf_lsm%pt_surface(m) * exner(k_topo)
             ENDDO
             DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                rrtm_tlev(0,k_topo+1) = surf_usm%pt_surface(m) * exner(k_topo)
             ENDDO


             IF ( bulk_cloud_model )  THEN
                DO  k = k_topo+1, nzt+1
                   rrtm_tlay(0,k) = pt(k,j,i) * exner(k) + lv_d_cp * ql(k,j,i)
                   rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * (q(k,j,i) - ql(k,j,i))
                ENDDO
             ELSEIF ( cloud_droplets )  THEN
                DO  k = k_topo+1, nzt+1
                   rrtm_tlay(0,k) = pt(k,j,i) * exner(k) + lv_d_cp * ql(k,j,i)
                   rrtm_h2ovmr(0,k) = mol_mass_air_d_wv * q(k,j,i)
                ENDDO
             ELSE
                DO  k = k_topo+1, nzt+1
                   rrtm_tlay(0,k) = pt(k,j,i) * exner(k)
                ENDDO

                IF ( humidity )  THEN
                   DO  k = k_topo+1, nzt+1
                      rrtm_h2ovmr(0,k) =  mol_mass_air_d_wv * q(k,j,i)
                   ENDDO
                ELSE
                   rrtm_h2ovmr(0,k_topo+1:nzt+1) = 0.0_wp
                ENDIF
             ENDIF

!
!--          Avoid temperature/humidity jumps at the top of the LES domain by linear interpolation
!--          from nzt+2 to nzt+7
             DO  k = nzt+2, nzt+7
                rrtm_tlay(0,k) = rrtm_tlay(0,nzt+1) + ( rrtm_tlay(0,nzt+8) - rrtm_tlay(0,nzt+1) ) /&
                               ( rrtm_play(0,nzt+8) - rrtm_play(0,nzt+1) ) *                       &
                               ( rrtm_play(0,k)     - rrtm_play(0,nzt+1) )

                rrtm_h2ovmr(0,k) = rrtm_h2ovmr(0,nzt+1) +                                          &
                                 ( rrtm_h2ovmr(0,nzt+8) - rrtm_h2ovmr(0,nzt+1) ) /                 &
                                 ( rrtm_play(0,nzt+8)   - rrtm_play(0,nzt+1)   ) *                 &
                                 ( rrtm_play(0,k)       - rrtm_play(0,nzt+1) )

             ENDDO

!--          Linear interpolation to zw grid
             DO  k = k_topo+2, nzt+8
                rrtm_tlev(0,k) = rrtm_tlay(0,k-1) + (rrtm_tlay(0,k) - rrtm_tlay(0,k-1) ) /         &
                               ( rrtm_play(0,k) - rrtm_play(0,k-1) ) *                             &
                               ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
             ENDDO


!
!--          Calculate liquid water path and cloud fraction for each column.
!--          Note that LWP is required in g/m2 instead of kg/kg m.
             rrtm_cldfr  = 0.0_wp
             rrtm_reliq  = 0.0_wp
             rrtm_reice  = 0.0_wp
             rrtm_cliqwp = 0.0_wp
             rrtm_cicewp = 0.0_wp
             rrtm_icld   = 0

             IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
                DO  k = k_topo+1, nzt+1
                   rrtm_cliqwp(0,k) =  ql(k,j,i) * 1000.0_wp *                                     &
                                     ( rrtm_plev(0,k) - rrtm_plev(0,k+1) ) * 100.0_wp / g
!
!--                Avoid rrtmg cloud calculation for very small values.
                   IF ( rrtm_cliqwp(0,k) < 1.0E-20_wp )  rrtm_cliqwp(0,k) = 0.0_wp

                   IF ( rrtm_cliqwp(0,k) > 0.0_wp )  THEN
                      rrtm_cldfr(0,k) = 1.0_wp
                      IF ( rrtm_icld == 0 )  rrtm_icld = 1

!
!--                   Calculate cloud droplet effective radius
                      IF ( bulk_cloud_model )  THEN
!
!--                      Calculate effective droplet radius. In case of using cloud_scheme =
!--                      'morrison' and a non reasonable number of cloud droplets the inital aerosol
!--                      number concentration is considered.
                         IF ( microphysics_morrison )  THEN
                            IF ( nc(k,j,i) > 1.0E-20_wp )  THEN
                               nc_rad = nc(k,j,i)
                            ELSE
                               nc_rad = na_init
                            ENDIF
                         ELSE
                            nc_rad = nc_const
                         ENDIF

                         rrtm_reliq(0,k) = 1.0E6_wp * ( 3.0_wp * ql(k,j,i) * rho_surface /         &
                                         ( 4.0_wp * pi * nc_rad * rho_l ) )**0.33333333333333_wp * &
                                           EXP( LOG( sigma_gc )**2 )
                      ELSEIF ( cloud_droplets )  THEN
                         number_of_particles = prt_count(k,j,i)

                         IF (number_of_particles <= 0)  CYCLE
                         particles => grid_particles(k,j,i)%particles(1:number_of_particles)
                         s_r2 = 0.0_wp
                         s_r3 = 0.0_wp

                         DO  n = 1, number_of_particles
                            IF ( particles(n)%particle_mask )  THEN
                               s_r2 = s_r2 + particles(n)%radius**2 * particles(n)%weight_factor
                               s_r3 = s_r3 + particles(n)%radius**3 * particles(n)%weight_factor
                            ENDIF
                         ENDDO

                         IF ( s_r2 > 0.0_wp )  rrtm_reliq(0,k) = s_r3 / s_r2

                      ENDIF
!
!--                   Limit effective radius.
                      IF ( rrtm_reliq(0,k) > 0.0_wp )  THEN
                         rrtm_reliq(0,k) = MAX( rrtm_reliq(0,k), 2.5_wp )
                         rrtm_reliq(0,k) = MIN( rrtm_reliq(0,k), 60.0_wp )
                     ENDIF
                   ENDIF
!
!--                Calculation of effective radius of ice (made for RRTM) based on ECHAM5
!--                documentation (Roeckner et al, MPI report 349). In this scheme graupel and
!--                snow is neglected.
                   IF ( microphysics_ice_phase )  THEN
!
!--                   Calculate ice water path in g/m2
                      rrtm_cicewp(0,k) =  qi(k,j,i) * 1000.0_wp *                            &
                               ( rrtm_plev(0,k) - rrtm_plev(0,k+1) ) * 100.0_wp / g
!
!--                   Avoid rrtmg cloud calculation for very small values.
                      IF ( rrtm_cicewp(0,k) < 1.0E-20_wp )  rrtm_cicewp(0,k) = 0.0_wp

                      IF ( rrtm_cicewp(0,k) > 0.0_wp )  THEN

                         rrtm_cldfr(0,k) = 1.0_wp
                         IF ( rrtm_icld == 0 )  rrtm_icld = 1
!
!--                      Calculate mean mass of particle.
                         mass_xi = hyrho(k) * qi(k,j,i) * 1000.0_wp
                         rrtm_reice(0,k) = a5 * EXP( b5 * LOG(mass_xi) )
!
!--                      Limit ice effective radius to allowed range for this parameterization
                         rrtm_reice(0,k) =  MAX( MIN( rrtm_reice(0,k) , 131.0_wp ), 5.0_wp )

                      ENDIF

                   ENDIF

                ENDDO
             ENDIF

!
!--          Write surface emissivity and surface temperature at current surface element on
!--          RRTMG-shaped array. Please note, as RRTMG is a single column model, surface attributes
!--          are only obtained from upward facing horizontally aligned surfaces (for simplicity).
!--          Taking surface attributes from horizontal and vertical walls would lead to multiple
!--          solutions. Moreover, for default, natural- and urban-type surfaces, several surface
!--          classes can exist at a surface element next to each other. To obtain bulk parameters,
!--          apply a weighted average for these surfaces.
             IF ( radiation_only )  THEN
                DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                   IF ( surf_def%upward(m) )  THEN
                      rrtm_emis = surf_def%emissivity(m,0)
                      rrtm_tsfc = surf_def%pt_surface(m) * exner(k_topo)
                   ENDIF
                ENDDO
             ENDIF
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                IF ( surf_lsm%upward(m) )  THEN
                   rrtm_emis = surf_lsm%frac(m,ind_veg_wall) *                                     &
                               surf_lsm%emissivity(m,ind_veg_wall) +                               &
                               surf_lsm%frac(m,ind_pav_green) *                                    &
                               surf_lsm%emissivity(m,ind_pav_green) +                              &
                               surf_lsm%frac(m,ind_wat_win) *                                      &
                               surf_lsm%emissivity(m,ind_wat_win)
                   rrtm_tsfc = surf_lsm%pt_surface(m) * exner(k_topo)
                ENDIF
             ENDDO
             DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                IF ( surf_usm%upward(m) )  THEN
                   rrtm_emis = surf_usm%frac(m,ind_veg_wall) *                                     &
                               surf_usm%emissivity(m,ind_veg_wall) +                               &
                               surf_usm%frac(m,ind_pav_green) *                                    &
                               surf_usm%emissivity(m,ind_pav_green) +                              &
                               surf_usm%frac(m,ind_wat_win) *                                      &
                               surf_usm%emissivity(m,ind_wat_win)
                   rrtm_tsfc = surf_usm%pt_surface(m) * exner(k_topo)
                ENDIF
             ENDDO

!
!--          DCEP effect on emis and pt.
             IF ( dcep )  THEN
!
!--             Adjust emissivity and surface temperature (see Schubert 2013 Eq. 3.68 and 3.69).
                rrtm_emis_save = rrtm_emis(0,1)
                rrtm_emis = rrtm_emis_save * ( 1.0_wp - fr_urb(j,i) ) +                            &
                            emiss_dcep(j,i) * fr_urb(j,i)
                rrtm_tsfc = ( ( emiss_dcep(j,i) * fr_urb(j,i) * t_grad_dcep(j,i)**4 +              &
                                rrtm_emis_save * ( 1.0_wp - fr_urb(j,i) ) * rrtm_tsfc**4           &
                              )                                                                    &
                              / rrtm_emis(0,1) )**0.25_wp
             ENDIF

             IF ( lw_radiation )  THEN
!
!--             Due to technical reasons, copy optical depth to dummy arguments which are allocated
!--             on the exact size as the rrtmg_lw is called. As one dimension is allocated with zero
!--             size, compiler complains that rank of the array does not match that of the
!--             assumed-shaped arguments in the RRTMG library. In order to avoid this, write to
!--             dummy arguments and pass the entire dummy array. Seems to be the only existing
!--             work-around.
                ALLOCATE( rrtm_lw_taucld_dum(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1) )
                ALLOCATE( rrtm_lw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1) )

                rrtm_lw_taucld_dum = rrtm_lw_taucld(1:nbndlw+1,0:0,k_topo+1:nzt_rad+1)
                rrtm_lw_tauaer_dum = rrtm_lw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndlw+1)

                CALL rrtmg_lw( 1, nzt_rad-k_topo, rrtm_icld, rrtm_idrv,                            &
                               rrtm_play(:,k_topo+1:nzt_rad+1), rrtm_plev(:,k_topo+1:nzt_rad+2),   &
                               rrtm_tlay(:,k_topo+1:nzt_rad+1), rrtm_tlev(:,k_topo+1:nzt_rad+2),   &
                               rrtm_tsfc, rrtm_h2ovmr(:,k_topo+1:nzt_rad+1),                       &
                               rrtm_o3vmr(:,k_topo+1:nzt_rad+1), rrtm_co2vmr(:,k_topo+1:nzt_rad+1),&
                               rrtm_ch4vmr(:,k_topo+1:nzt_rad+1),                                  &
                               rrtm_n2ovmr(:,k_topo+1:nzt_rad+1),                                  &
                               rrtm_o2vmr(:,k_topo+1:nzt_rad+1),                                   &
                               rrtm_cfc11vmr(:,k_topo+1:nzt_rad+1),                                &
                               rrtm_cfc12vmr(:,k_topo+1:nzt_rad+1),                                &
                               rrtm_cfc22vmr(:,k_topo+1:nzt_rad+1),                                &
                               rrtm_ccl4vmr(:,k_topo+1:nzt_rad+1), rrtm_emis, rrtm_inflglw,        &
                               rrtm_iceflglw, rrtm_liqflglw, rrtm_cldfr(:,k_topo+1:nzt_rad+1),     &
                               rrtm_lw_taucld_dum, rrtm_cicewp(:,k_topo+1:nzt_rad+1),              &
                               rrtm_cliqwp(:,k_topo+1:nzt_rad+1), rrtm_reice(:,k_topo+1:nzt_rad+1),&
                               rrtm_reliq(:,k_topo+1:nzt_rad+1), rrtm_lw_tauaer_dum,               &
                               rrtm_lwuflx(:,k_topo:nzt_rad+1), rrtm_lwdflx(:,k_topo:nzt_rad+1),   &
                               rrtm_lwhr(:,k_topo+1:nzt_rad+1), rrtm_lwuflxc(:,k_topo:nzt_rad+1),  &
                               rrtm_lwdflxc(:,k_topo:nzt_rad+1), rrtm_lwhrc(:,k_topo+1:nzt_rad+1), &
                               rrtm_lwuflx_dt(:,k_topo:nzt_rad+1),                                 &
                               rrtm_lwuflxc_dt(:,k_topo:nzt_rad+1) )

                DEALLOCATE( rrtm_lw_taucld_dum )
                DEALLOCATE( rrtm_lw_tauaer_dum )
!
!--             Save fluxes
                DO  k = k_topo, nzt+1
                   rad_lw_in(k,j,i)  = rrtm_lwdflx(0,k)
                   rad_lw_out(k,j,i) = rrtm_lwuflx(0,k)
                ENDDO

                IF ( dcep )  rad_lw_in_diff(j,i) = rad_lw_in(k_topo,j,i)
!
!--             Save heating rates (convert from K/d to K/h)
                DO  k = k_topo+1, nzt+1
                   rad_lw_hr(k,j,i)    = rrtm_lwhr(0,k-k_topo)  * d_hours_day
                   rad_lw_cs_hr(k,j,i) = rrtm_lwhrc(0,k-k_topo) * d_hours_day
                ENDDO
!
!--             Save surface radiative fluxes and change in LW heating rate onto respective surface
!--             elements
                DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                   k                               = surf_lsm%k(m)
                   surf_lsm%rad_lw_in(m)           = rrtm_lwdflx(0,k_topo)
                   surf_lsm%rad_lw_out(m)          = rrtm_lwuflx(0,k_topo)
                   surf_lsm%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k_topo)
                ENDDO
                DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                   k                               = surf_usm%k(m)
                   surf_usm%rad_lw_in(m)           = rrtm_lwdflx(0,k_topo)
                   surf_usm%rad_lw_out(m)          = rrtm_lwuflx(0,k_topo)
                   surf_usm%rad_lw_out_change_0(m) = rrtm_lwuflx_dt(0,k_topo)
                ENDDO

             ENDIF

             IF ( sw_radiation  .AND.  sun_up )  THEN
!
!--             Get albedo for direct/diffusive long/shortwave radiation at current (y,x)-location
!--             from surface variables. Only obtain it from upward facing horizontal surfaces,
!--             as RRTMG is a single column model. (Please note, only one loop will be entered,
!--             controlled by start-end index.)
                IF ( use_broadband_albedo )  THEN
                   IF ( radiation_only )  THEN
                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward(m) )  THEN
                            rrtm_asdir(1) = surf_def%albedo(m,0)
                            rrtm_asdif(1) = surf_def%albedo(m,0)
                            rrtm_aldir(1) = surf_def%albedo(m,0)
                            rrtm_aldif(1) = surf_def%albedo(m,0)
                         ENDIF
                      ENDDO
                   ENDIF
                   DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                      IF ( surf_lsm%upward(m) )  THEN
                         rrtm_asdir(1) = SUM( surf_lsm%frac(m,:) * surf_lsm%albedo(m,:) )
                         rrtm_asdif(1) = SUM( surf_lsm%frac(m,:) * surf_lsm%albedo(m,:) )
                         rrtm_aldir(1) = SUM( surf_lsm%frac(m,:) * surf_lsm%albedo(m,:) )
                         rrtm_aldif(1) = SUM( surf_lsm%frac(m,:) * surf_lsm%albedo(m,:) )
                      ENDIF
                   ENDDO
                   DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                      IF ( surf_usm%upward(m) )  THEN
                         rrtm_asdir(1) = SUM( surf_usm%frac(m,:) * surf_usm%albedo(m,:) )
                         rrtm_asdif(1) = SUM( surf_usm%frac(m,:) * surf_usm%albedo(m,:) )
                         rrtm_aldir(1) = SUM( surf_usm%frac(m,:) * surf_usm%albedo(m,:) )
                         rrtm_aldif(1) = SUM( surf_usm%frac(m,:) * surf_usm%albedo(m,:) )
                      ENDIF
                   ENDDO
                ELSE
                   IF ( radiation_only )  THEN
                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         IF ( surf_def%upward(m) )  THEN
                            rrtm_asdir(1) = surf_def%albedo(m,0)
                            rrtm_asdif(1) = surf_def%albedo(m,0)
                            rrtm_aldir(1) = surf_def%albedo(m,0)
                            rrtm_aldif(1) = surf_def%albedo(m,0)
                         ENDIF
                      ENDDO
                   ENDIF
                   DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                      IF ( surf_lsm%upward(m) )  THEN
                         rrtm_asdir(1) = SUM( surf_lsm%frac(m,:) * surf_lsm%rrtm_asdir(m,:) )
                         rrtm_asdif(1) = SUM( surf_lsm%frac(m,:) * surf_lsm%rrtm_asdif(m,:) )
                         rrtm_aldir(1) = SUM( surf_lsm%frac(m,:) * surf_lsm%rrtm_aldir(m,:) )
                         rrtm_aldif(1) = SUM( surf_lsm%frac(m,:) * surf_lsm%rrtm_aldif(m,:) )
                      ENDIF
                   ENDDO
                   DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                      IF ( surf_usm%upward(m) )  THEN
                         rrtm_asdir(1) = SUM( surf_usm%frac(m,:) * surf_usm%rrtm_asdir(m,:) )
                         rrtm_asdif(1) = SUM( surf_usm%frac(m,:) * surf_usm%rrtm_asdif(m,:) )
                         rrtm_aldir(1) = SUM( surf_usm%frac(m,:) * surf_usm%rrtm_aldir(m,:) )
                         rrtm_aldif(1) = SUM( surf_usm%frac(m,:) * surf_usm%rrtm_aldif(m,:) )
                      ENDIF
                   ENDDO
                ENDIF
!
!--             Edit value in case of DCEP.
                IF ( dcep ) THEN
                   rrtm_asdir(1) = rrtm_asdir(1) * ( 1.0_wp - fr_urb(j,i) )                        &
                                   + albedop_dcep(j,i) * fr_urb(j,i)
                   rrtm_asdif(1) = rrtm_asdif(1) * ( 1.0_wp - fr_urb(j,i) )                        &
                                   + albedop_dcep(j,i) * fr_urb(j,i)
                   rrtm_aldir(1) = rrtm_aldir(1) * ( 1.0_wp - fr_urb(j,i) )                        &
                                   + albedop_dcep(j,i) * fr_urb(j,i)
                   rrtm_aldif(1) = rrtm_aldif(1) * ( 1.0_wp - fr_urb(j,i) )                        &
                                   + albedop_dcep(j,i) * fr_urb(j,i)
                ENDIF

!
!--             Due to technical reasons, copy optical depths and other to dummy arguments which are
!--             allocated on the exact size as the rrtmg_sw is called. As one dimension is allocated
!--             with zero size, compiler complains that rank of the array does not match that of the
!--             assumed-shaped arguments in the RRTMG library. In order to avoid this, write to
!--             dummy arguments and pass the entire dummy array. Seems to be the only existing
!--             work-around.
                ALLOCATE( rrtm_sw_taucld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                ALLOCATE( rrtm_sw_ssacld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                ALLOCATE( rrtm_sw_asmcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                ALLOCATE( rrtm_sw_fsfcld_dum(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1) )
                ALLOCATE( rrtm_sw_tauaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
                ALLOCATE( rrtm_sw_ssaaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
                ALLOCATE( rrtm_sw_asmaer_dum(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1) )
                ALLOCATE( rrtm_sw_ecaer_dum(0:0,k_topo+1:nzt_rad+1,1:naerec+1)  )

                rrtm_sw_taucld_dum = rrtm_sw_taucld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                rrtm_sw_ssacld_dum = rrtm_sw_ssacld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                rrtm_sw_asmcld_dum = rrtm_sw_asmcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                rrtm_sw_fsfcld_dum = rrtm_sw_fsfcld(1:nbndsw+1,0:0,k_topo+1:nzt_rad+1)
                rrtm_sw_tauaer_dum = rrtm_sw_tauaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
                rrtm_sw_ssaaer_dum = rrtm_sw_ssaaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
                rrtm_sw_asmaer_dum = rrtm_sw_asmaer(0:0,k_topo+1:nzt_rad+1,1:nbndsw+1)
                rrtm_sw_ecaer_dum  = rrtm_sw_ecaer(0:0,k_topo+1:nzt_rad+1,1:naerec+1)

                CALL rrtmg_sw( 1, nzt_rad-k_topo, rrtm_icld, rrtm_iaer,                            &
                               rrtm_play(:,k_topo+1:nzt_rad+1), rrtm_plev(:,k_topo+1:nzt_rad+2),   &
                               rrtm_tlay(:,k_topo+1:nzt_rad+1), rrtm_tlev(:,k_topo+1:nzt_rad+2),   &
                               rrtm_tsfc, rrtm_h2ovmr(:,k_topo+1:nzt_rad+1),                       &
                               rrtm_o3vmr(:,k_topo+1:nzt_rad+1), rrtm_co2vmr(:,k_topo+1:nzt_rad+1),&
                               rrtm_ch4vmr(:,k_topo+1:nzt_rad+1),                                  &
                               rrtm_n2ovmr(:,k_topo+1:nzt_rad+1),                                  &
                               rrtm_o2vmr(:,k_topo+1:nzt_rad+1), rrtm_asdir, rrtm_asdif,           &
                               rrtm_aldir, rrtm_aldif, zenith, 0.0_wp, day_of_year,                &
                               solar_constant, rrtm_inflgsw, rrtm_iceflgsw, rrtm_liqflgsw,         &
                               rrtm_cldfr(:,k_topo+1:nzt_rad+1), rrtm_sw_taucld_dum,               &
                               rrtm_sw_ssacld_dum, rrtm_sw_asmcld_dum, rrtm_sw_fsfcld_dum,         &
                               rrtm_cicewp(:,k_topo+1:nzt_rad+1),                                  &
                               rrtm_cliqwp(:,k_topo+1:nzt_rad+1),                                  &
                               rrtm_reice(:,k_topo+1:nzt_rad+1), rrtm_reliq(:,k_topo+1:nzt_rad+1), &
                               rrtm_sw_tauaer_dum, rrtm_sw_ssaaer_dum, rrtm_sw_asmaer_dum,         &
                               rrtm_sw_ecaer_dum, rrtm_swuflx(:,k_topo:nzt_rad+1),                 &
                               rrtm_swdflx(:,k_topo:nzt_rad+1), rrtm_swhr(:,k_topo+1:nzt_rad+1),   &
                               rrtm_swuflxc(:,k_topo:nzt_rad+1), rrtm_swdflxc(:,k_topo:nzt_rad+1), &
                               rrtm_swhrc(:,k_topo+1:nzt_rad+1), rrtm_dirdflux(:,k_topo:nzt_rad+1),&
                               rrtm_difdflux(:,k_topo:nzt_rad+1) )

                DEALLOCATE( rrtm_sw_taucld_dum )
                DEALLOCATE( rrtm_sw_ssacld_dum )
                DEALLOCATE( rrtm_sw_asmcld_dum )
                DEALLOCATE( rrtm_sw_fsfcld_dum )
                DEALLOCATE( rrtm_sw_tauaer_dum )
                DEALLOCATE( rrtm_sw_ssaaer_dum )
                DEALLOCATE( rrtm_sw_asmaer_dum )
                DEALLOCATE( rrtm_sw_ecaer_dum )
!
!--             Save fluxes
                DO  k = k_topo, nzt+1
                   rad_sw_in(k,j,i)  = rrtm_swdflx(0,k)
                   rad_sw_out(k,j,i) = rrtm_swuflx(0,k)
                ENDDO

!
!--             Save direct and diffuse SW radiation at the surface (required by DCEP).
                IF ( dcep )  THEN
                   rad_sw_in_dir(:,:)  = rrtm_dirdflux(0,k_topo)
                   rad_sw_in_diff(:,:) = rrtm_difdflux(0,k_topo)
                ENDIF

!
!--             Save heating rates (convert from K/d to K/s)
                DO  k = k_topo+1, nzt+1
                   rad_sw_hr(k,j,i)    = rrtm_swhr(0,k)  * d_hours_day
                   rad_sw_cs_hr(k,j,i) = rrtm_swhrc(0,k) * d_hours_day
                ENDDO

!
!--             Save surface radiative fluxes onto respective surface elements
                DO   m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                   k                      = surf_lsm%k(m)
                   surf_lsm%rad_sw_in(m)  = rrtm_swdflx(0,k_topo)
                   surf_lsm%rad_sw_out(m) = rrtm_swuflx(0,k_topo)
                ENDDO
                DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                   k                      = surf_usm%k(m)
                   surf_usm%rad_sw_in(m)  = rrtm_swdflx(0,k_topo)
                   surf_usm%rad_sw_out(m) = rrtm_swuflx(0,k_topo)
                ENDDO
!
!--          Solar radiation is zero during night
             ELSE
                rad_sw_in  = 0.0_wp
                rad_sw_out = 0.0_wp

                IF ( dcep )  THEN
                   rad_sw_in_dir(j,i)  = 0.0_wp
                   rad_sw_in_diff(j,i) = 0.0_wp
                ENDIF
!
!--             Surface radiative fluxes should be also set to zero here to account for zero
!--             incoming radiation.
!--             Save surface radiative fluxes onto respective surface elements.
                DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                   surf_lsm%rad_sw_in(m)  = 0.0_wp
                   surf_lsm%rad_sw_out(m) = 0.0_wp
                ENDDO
                DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                   surf_usm%rad_sw_in(m)  = 0.0_wp
                   surf_usm%rad_sw_out(m) = 0.0_wp
                ENDDO
             ENDIF

          ENDDO
       ENDDO

    ENDIF
!
!-- Finally, calculate surface net radiation for surface elements.
    IF ( .NOT. radiation_interactions )  THEN
!
!--    Todo: weight with azimuth and zenith angle according to their orientation!
       DO  m = 1, surf_lsm%ns
          surf_lsm%rad_net(m) = surf_lsm%rad_sw_in(m) - surf_lsm%rad_sw_out(m) +                   &
                                surf_lsm%rad_lw_in(m) - surf_lsm%rad_lw_out(m)
       ENDDO
       DO  m = 1, surf_usm%ns
          surf_usm%rad_net(m) = surf_usm%rad_sw_in(m) - surf_usm%rad_sw_out(m) +                   &
                                surf_usm%rad_lw_in(m) - surf_usm%rad_lw_out(m)
       ENDDO
    ENDIF


    CALL exchange_horiz( rad_lw_in, nbgp )
    CALL exchange_horiz( rad_lw_out, nbgp )
    CALL exchange_horiz( rad_lw_hr, nbgp )
    CALL exchange_horiz( rad_lw_cs_hr, nbgp )

    CALL exchange_horiz( rad_sw_in, nbgp )
    CALL exchange_horiz( rad_sw_out, nbgp )
    CALL exchange_horiz( rad_sw_hr, nbgp )
    CALL exchange_horiz( rad_sw_cs_hr, nbgp )
#endif

 END SUBROUTINE radiation_rrtmg


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the cosine of the zenith angle (variable is called zenith)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_zenith( day_of_year, second_of_day )

    USE palm_date_time_mod,                                                                        &
        ONLY:  seconds_per_day

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  day_of_year  !< day of the year

    REAL(wp) ::  declination  !< solar declination angle
    REAL(wp) ::  hour_angle   !< solar hour angle

    REAL(wp),INTENT(IN) ::  second_of_day  !< current time of the day in UTC

!
!-- Calculate solar declination and hour angle
    declination = ASIN( decl_1 * SIN( decl_2 * REAL( day_of_year, KIND = wp ) - decl_3 ) )
    hour_angle  = 2.0_wp * pi * ( second_of_day / seconds_per_day ) + lon - pi

!
!-- Calculate cosine of solar zenith angle
    cos_zenith = SIN( lat ) * SIN( declination ) + COS( lat ) * COS( declination ) *               &
                 COS( hour_angle )
    cos_zenith = MAX( 0.0_wp, cos_zenith )

!
!-- Calculate solar directional vector
    IF ( sun_direction )  THEN

!
!--    Direction in longitudes equals to sin(solar_azimuth) * sin(zenith)
       sun_dir_lon = - SIN( hour_angle ) * COS( declination )

!
!--    Direction in latitues equals to cos(solar_azimuth) * sin(zenith)
       sun_dir_lat = SIN( declination ) * COS( lat ) - COS( hour_angle ) * COS( declination ) *    &
                     SIN( lat )
    ENDIF

!
!-- Check if the sun is up (otheriwse shortwave calculations can be skipped)
    IF ( cos_zenith > 0.0_wp )  THEN
       sun_up = .TRUE.
    ELSE
       sun_up = .FALSE.
    ENDIF

 END SUBROUTINE calc_zenith

#if defined( __rrtmg ) && defined( __netcdf )
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates surface albedo components based on Briegleb (1992) and Briegleb et al. (1986)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_albedo( surf )

    IMPLICIT NONE

    INTEGER(iwp) ::  ind_type  !< running index surface tiles
    INTEGER(iwp) ::  m         !< running index surface elements

    TYPE(surf_type) ::  surf  !< treated surfaces

    IF ( sun_up  .AND.  .NOT. average_radiation  .OR.  dcep )  THEN

       DO  m = 1, surf%ns
!
!--       Loop over surface elements
          DO  ind_type = 0, SIZE( surf%albedo_type, 2 ) - 1

!
!--          Ocean
             IF ( surf%albedo_type(m,ind_type) == 1 )  THEN
                surf%rrtm_aldir(m,ind_type) = 0.026_wp / ( cos_zenith**1.7_wp + 0.065_wp ) +       &
                                              0.15_wp * ( cos_zenith - 0.1_wp ) *                  &
                                              ( cos_zenith - 0.5_wp ) * ( cos_zenith - 1.0_wp )
                surf%rrtm_asdir(m,ind_type) = surf%rrtm_aldir(m,ind_type)
!
!--          Snow
             ELSEIF ( surf%albedo_type(m,ind_type) == 16 )  THEN
                IF ( cos_zenith < 0.5_wp )  THEN
                   surf%rrtm_aldir(m,ind_type) = 0.5_wp * ( 1.0_wp - surf%aldif(m,ind_type) ) *    &
                                                 ( ( 3.0_wp / ( 1.0_wp + 4.0_wp * cos_zenith ) ) - &
                                                 1.0_wp )
                   surf%rrtm_asdir(m,ind_type) = 0.5_wp * ( 1.0_wp - surf%asdif(m,ind_type) ) *    &
                                                 ( ( 3.0_wp / ( 1.0_wp + 4.0_wp * cos_zenith ) ) - &
                                                 1.0_wp )

                   surf%rrtm_aldir(m,ind_type) = MIN( 0.98_wp, surf%rrtm_aldir(m,ind_type) )
                   surf%rrtm_asdir(m,ind_type) = MIN( 0.98_wp, surf%rrtm_asdir(m,ind_type) )
                ELSE
                   surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type)
                   surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type)
                ENDIF
!
!--          Sea ice
             ELSEIF ( surf%albedo_type(m,ind_type) == 15 )  THEN
                surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type)
                surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type)

!
!--          Asphalt
             ELSEIF ( surf%albedo_type(m,ind_type) == 17 )  THEN
                surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type)
                surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type)


!
!--          Bare soil
             ELSEIF ( surf%albedo_type(m,ind_type) == 18 )  THEN
                surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type)
                surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type)

!
!--          Land surfaces
             ELSE
                SELECT CASE ( surf%albedo_type(m,ind_type) )

!
!--                Surface types with strong zenith dependence
                   CASE ( 1, 2, 3, 4, 11, 12, 13 )
                      surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type) * 1.4_wp /              &
                                                    ( 1.0_wp + 0.8_wp * cos_zenith )
                      surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type) * 1.4_wp /              &
                                                    ( 1.0_wp + 0.8_wp * cos_zenith )
!
!--                Surface types with weak zenith dependence
                   CASE ( 5, 6, 7, 8, 9, 10, 14 )
                      surf%rrtm_aldir(m,ind_type) = surf%aldif(m,ind_type) * 1.1_wp /              &
                                                    ( 1.0_wp + 0.2_wp * cos_zenith )
                      surf%rrtm_asdir(m,ind_type) = surf%asdif(m,ind_type) * 1.1_wp /              &
                                                    ( 1.0_wp + 0.2_wp * cos_zenith )

                   CASE DEFAULT

                END SELECT
             ENDIF
!
!--          Diffusive albedo is taken from Table 2
             surf%rrtm_aldif(m,ind_type) = surf%aldif(m,ind_type)
             surf%rrtm_asdif(m,ind_type) = surf%asdif(m,ind_type)
          ENDDO
       ENDDO
!
!-- Set albedo in case of average radiation
    ELSEIF ( sun_up  .AND.  average_radiation )  THEN
       surf%rrtm_asdir = albedo_eff
       surf%rrtm_asdif = albedo_eff
       surf%rrtm_aldir = albedo_eff
       surf%rrtm_aldif = albedo_eff
!
!-- Darkness
    ELSE
       surf%rrtm_aldir = 0.0_wp
       surf%rrtm_asdir = 0.0_wp
       surf%rrtm_aldif = 0.0_wp
       surf%rrtm_asdif = 0.0_wp
    ENDIF

 END SUBROUTINE calc_albedo

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read sounding data (pressure and temperature) from RADIATION_DATA.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE read_sounding_data

    IMPLICIT NONE

    INTEGER(iwp) ::  id,           &  !< NetCDF id of input file
                     id_dim_zrad,  &  !< pressure level id in the NetCDF file
                     id_var,       &  !< NetCDF variable id
                     k,            &  !< loop index
                     nz_snd,       &  !< number of vertical levels in the sounding data
                     nz_snd_start, &  !< start vertical index for sounding data to be used
                     nz_snd_end       !< end vertical index for souding data to be used

    REAL(wp) ::  t_surface  !< actual surface temperature

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  hyp_snd_tmp, &  !< temporary hydrostatic pressure profile (sounding)
                                            t_snd_tmp       !< temporary temperature profile (sounding)

!
!-- In case of updates, deallocate arrays first (sufficient to check one array as the others are
!-- automatically allocated). This is required because nzt_rad might change during the update
    IF ( ALLOCATED( hyp_snd ) )  THEN
       DEALLOCATE( hyp_snd )
       DEALLOCATE( t_snd )
       DEALLOCATE( rrtm_play )
       DEALLOCATE( rrtm_plev )
       DEALLOCATE( rrtm_tlay )
       DEALLOCATE( rrtm_tlev )

       DEALLOCATE( rrtm_cicewp )
       DEALLOCATE( rrtm_cldfr )
       DEALLOCATE( rrtm_cliqwp )
       DEALLOCATE( rrtm_reice )
       DEALLOCATE( rrtm_reliq )
       DEALLOCATE( rrtm_lw_taucld )
       DEALLOCATE( rrtm_lw_tauaer )

       DEALLOCATE( rrtm_lwdflx )
       DEALLOCATE( rrtm_lwdflxc )
       DEALLOCATE( rrtm_lwuflx )
       DEALLOCATE( rrtm_lwuflxc )
       DEALLOCATE( rrtm_lwuflx_dt )
       DEALLOCATE( rrtm_lwuflxc_dt )
       DEALLOCATE( rrtm_lwhr )
       DEALLOCATE( rrtm_lwhrc )

       DEALLOCATE( rrtm_sw_taucld )
       DEALLOCATE( rrtm_sw_ssacld )
       DEALLOCATE( rrtm_sw_asmcld )
       DEALLOCATE( rrtm_sw_fsfcld )
       DEALLOCATE( rrtm_sw_tauaer )
       DEALLOCATE( rrtm_sw_ssaaer )
       DEALLOCATE( rrtm_sw_asmaer )
       DEALLOCATE( rrtm_sw_ecaer )

       DEALLOCATE( rrtm_swdflx )
       DEALLOCATE( rrtm_swdflxc )
       DEALLOCATE( rrtm_swuflx )
       DEALLOCATE( rrtm_swuflxc )
       DEALLOCATE( rrtm_swhr )
       DEALLOCATE( rrtm_swhrc )
       DEALLOCATE( rrtm_dirdflux )
       DEALLOCATE( rrtm_difdflux )

    ENDIF

!
!-- Open file for reading
    nc_stat = NF90_OPEN( rrtm_input_file, NF90_NOWRITE, id )
    CALL netcdf_handle_error_rad( 'read_sounding_data', 549 )

!
!-- Inquire dimension of z axis and save in nz_snd
    nc_stat = NF90_INQ_DIMID( id, 'Pressure', id_dim_zrad )
    nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim_zrad, len = nz_snd )
    CALL netcdf_handle_error_rad( 'read_sounding_data', 551 )

!
!-- Allocate temporary array for storing pressure data
    ALLOCATE( hyp_snd_tmp(1:nz_snd) )
    hyp_snd_tmp = 0.0_wp


!-- Read pressure from file
    nc_stat = NF90_INQ_VARID( id, 'Pressure', id_var )
    nc_stat = NF90_GET_VAR( id, id_var, hyp_snd_tmp(:), start = (/1/), count = (/nz_snd/) )
    CALL netcdf_handle_error_rad( 'read_sounding_data', 552 )

!
!-- Allocate temporary array for storing temperature data
    ALLOCATE( t_snd_tmp(1:nz_snd) )
    t_snd_tmp = 0.0_wp

!
!-- Read temperature from file
    nc_stat = NF90_INQ_VARID( id, 'ReferenceTemperature', id_var )
    nc_stat = NF90_GET_VAR( id, id_var, t_snd_tmp(:), start = (/1/), count = (/nz_snd/) )
    CALL netcdf_handle_error_rad( 'read_sounding_data', 553 )

!
!-- Calculate start of sounding data
    nz_snd_start = nz_snd + 1
    nz_snd_end   = nz_snd + 1

!
!-- Start filling vertical dimension at 10hPa above the model domain (hyp is in Pa, hyp_snd in hPa).
    DO  k = 1, nz_snd
       IF ( hyp_snd_tmp(k) < ( hyp(nzt+1) - 1000.0_wp) * 0.01_wp )  THEN
          nz_snd_start = k
          EXIT
       ENDIF
    ENDDO

    IF ( nz_snd_start <= nz_snd )  THEN
       nz_snd_end = nz_snd
    ENDIF


!
!-- Calculate of total grid points for RRTMG calculations
    nzt_rad = nzt + nz_snd_end - nz_snd_start + 1

!
!-- Save data above LES domain in hyp_snd, t_snd
    ALLOCATE( hyp_snd(nzb+1:nzt_rad) )
    ALLOCATE( t_snd(nzb+1:nzt_rad) )
    hyp_snd = 0.0_wp
    t_snd = 0.0_wp

    hyp_snd(nzt+2:nzt_rad) = hyp_snd_tmp(nz_snd_start+1:nz_snd_end)
    t_snd(nzt+2:nzt_rad)   = t_snd_tmp(nz_snd_start+1:nz_snd_end)

    nc_stat = NF90_CLOSE( id )

!
!-- Calculate pressure levels on zu and zw grid. Sounding data is added at top of the LES domain.
!-- This routine does not consider horizontal or vertical variability of pressure and temperature
    ALLOCATE( rrtm_play(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_plev(0:0,nzb+1:nzt_rad+2) )

    t_surface = pt_surface * exner(nzb)
    DO  k = nzb+1, nzt+1
       rrtm_play(0,k) = hyp(k) * 0.01_wp
       rrtm_plev(0,k) = barometric_formula(zw(k-1), pt_surface * exner(nzb), surface_pressure )
    ENDDO

    DO  k = nzt+2, nzt_rad
       rrtm_play(0,k) = hyp_snd(k)
       rrtm_plev(0,k) = 0.5_wp * ( rrtm_play(0,k) + rrtm_play(0,k-1) )
    ENDDO
    rrtm_plev(0,nzt_rad+1) = MAX( 0.5 * hyp_snd(nzt_rad), 1.5 * hyp_snd(nzt_rad) - 0.5 *           &
                                  hyp_snd(nzt_rad-1) )
    rrtm_plev(0,nzt_rad+2) = MIN( 1.0E-4_wp, 0.25_wp * rrtm_plev(0,nzt_rad+1) )

    rrtm_play(0,nzt_rad+1) = 0.5 * rrtm_plev(0,nzt_rad+1)

!
!-- Calculate temperature/humidity levels at top of the LES domain.
!-- Currently, the temperature is taken from sounding data (might lead to a temperature jump at
!-- interface. To do: Humidity is currently not calculated above the LES domain.
    ALLOCATE( rrtm_tlay(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_tlev(0:0,nzb+1:nzt_rad+2) )

    DO  k = nzt+8, nzt_rad
       rrtm_tlay(0,k) = t_snd(k)
    ENDDO
    rrtm_tlay(0,nzt_rad+1) = 2.0_wp * rrtm_tlay(0,nzt_rad) - rrtm_tlay(0,nzt_rad-1)
    DO  k = nzt+9, nzt_rad+1
       rrtm_tlev(0,k) = rrtm_tlay(0,k-1) + ( rrtm_tlay(0,k) - rrtm_tlay(0,k-1) ) /                 &
                        ( rrtm_play(0,k) - rrtm_play(0,k-1) ) *                                    &
                        ( rrtm_plev(0,k) - rrtm_play(0,k-1) )
    ENDDO

    rrtm_tlev(0,nzt_rad+2) = 2.0_wp * rrtm_tlay(0,nzt_rad+1) - rrtm_tlev(0,nzt_rad)
!
!-- Allocate remaining RRTMG arrays
    ALLOCATE( rrtm_cicewp(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_cldfr(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_cliqwp(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_reice(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_reliq(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_lw_taucld(1:nbndlw+1,0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_lw_tauaer(0:0,nzb+1:nzt_rad+1,1:nbndlw+1) )
    ALLOCATE( rrtm_sw_taucld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_sw_ssacld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_sw_asmcld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_sw_fsfcld(1:nbndsw+1,0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_sw_tauaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
    ALLOCATE( rrtm_sw_ssaaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
    ALLOCATE( rrtm_sw_asmaer(0:0,nzb+1:nzt_rad+1,1:nbndsw+1) )
    ALLOCATE( rrtm_sw_ecaer(0:0,nzb+1:nzt_rad+1,1:naerec+1) )

!
!-- The ice phase is currently not considered in PALM
    rrtm_cicewp = 0.0_wp
    rrtm_reice  = 0.0_wp

!
!-- Set other parameters (move to NAMELIST parameters in the future)
    rrtm_lw_tauaer = 0.0_wp
    rrtm_lw_taucld = 0.0_wp
    rrtm_sw_taucld = 0.0_wp
    rrtm_sw_ssacld = 0.0_wp
    rrtm_sw_asmcld = 0.0_wp
    rrtm_sw_fsfcld = 0.0_wp
    rrtm_sw_tauaer = 0.0_wp
    rrtm_sw_ssaaer = 0.0_wp
    rrtm_sw_asmaer = 0.0_wp
    rrtm_sw_ecaer  = 0.0_wp


    ALLOCATE( rrtm_swdflx(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_swuflx(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_swhr(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_swuflxc(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_swdflxc(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_swhrc(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_dirdflux(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_difdflux(0:0,nzb:nzt_rad+1) )

    rrtm_swdflx  = 0.0_wp
    rrtm_swuflx  = 0.0_wp
    rrtm_swhr    = 0.0_wp
    rrtm_swuflxc = 0.0_wp
    rrtm_swdflxc = 0.0_wp
    rrtm_swhrc   = 0.0_wp
    rrtm_dirdflux = 0.0_wp
    rrtm_difdflux = 0.0_wp

    ALLOCATE( rrtm_lwdflx(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_lwuflx(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_lwhr(0:0,nzb+1:nzt_rad+1) )
    ALLOCATE( rrtm_lwuflxc(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_lwdflxc(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_lwhrc(0:0,nzb+1:nzt_rad+1) )

    rrtm_lwdflx  = 0.0_wp
    rrtm_lwuflx  = 0.0_wp
    rrtm_lwhr    = 0.0_wp
    rrtm_lwuflxc = 0.0_wp
    rrtm_lwdflxc = 0.0_wp
    rrtm_lwhrc   = 0.0_wp

    ALLOCATE( rrtm_lwuflx_dt(0:0,nzb:nzt_rad+1) )
    ALLOCATE( rrtm_lwuflxc_dt(0:0,nzb:nzt_rad+1) )

    rrtm_lwuflx_dt = 0.0_wp
    rrtm_lwuflxc_dt = 0.0_wp

 END SUBROUTINE read_sounding_data


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read trace gas data from file and convert into trace gas paths / volume mixing ratios. If a
!> user-defined input file is provided it needs to follow the convections used in RRTMG (see
!> respective netCDF files shipped with RRTMG)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE read_trace_gas_data

    USE rrsw_ncpar

    IMPLICIT NONE

    INTEGER(iwp), PARAMETER ::  num_trace_gases = 10  !< number of trace gases (absorbers)

    CHARACTER(LEN=5), DIMENSION(num_trace_gases), PARAMETER ::      &  !< trace gas names
    trace_names = (/'O3   ', 'CO2  ', 'CH4  ', 'N2O  ', 'O2   ',    &
                    'CFC11', 'CFC12', 'CFC22', 'CCL4 ', 'H2O  '/)

    INTEGER(iwp) ::  id,     &  !< NetCDF id
                     k,      &  !< loop index
                     m,      &  !< loop index
                     n,      &  !< loop index
                     nabs,   &  !< number of absorbers
                     np,     &  !< number of pressure levels
                     id_abs, &  !< NetCDF id of the respective absorber
                     id_dim, &  !< NetCDF id of asborber's dimension
                     id_var     !< NetCDf id ot the absorber

    REAL(wp) ::  p_mls_l, &  !< pressure lower limit for interpolation
                 p_mls_u, &  !< pressure upper limit for interpolation
                 p_wgt_l, &  !< pressure weight lower limit for interpolation
                 p_wgt_u, &  !< pressure weight upper limit for interpolation
                 p_mls_m     !< mean pressure between upper and lower limits


    REAL(wp), DIMENSION(:), ALLOCATABLE ::  p_mls,          &  !< pressure levels for the absorbers
                                            rrtm_play_tmp,  &  !< temporary array for pressure zu-levels
                                            rrtm_plev_tmp,  &  !< temporary array for pressure zw-levels
                                            trace_path_tmp     !< temporary array for storing trace gas path data

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  trace_mls,      &  !< array for storing the absorber amounts
                                              trace_mls_path, &  !< array for storing trace gas path data
                                              trace_mls_tmp      !< temporary array for storing trace gas data


!
!-- In case of updates, deallocate arrays first (sufficient to check one array as the others are
!-- automatically allocated)
    IF ( ALLOCATED( rrtm_o3vmr ) )  THEN
       DEALLOCATE( rrtm_o3vmr )
       DEALLOCATE( rrtm_co2vmr )
       DEALLOCATE( rrtm_ch4vmr )
       DEALLOCATE( rrtm_n2ovmr )
       DEALLOCATE( rrtm_o2vmr )
       DEALLOCATE( rrtm_cfc11vmr )
       DEALLOCATE( rrtm_cfc12vmr )
       DEALLOCATE( rrtm_cfc22vmr )
       DEALLOCATE( rrtm_ccl4vmr )
       DEALLOCATE( rrtm_h2ovmr )
    ENDIF

!
!-- Allocate trace gas profiles
    ALLOCATE( rrtm_o3vmr(0:0,1:nzt_rad+1) )
    ALLOCATE( rrtm_co2vmr(0:0,1:nzt_rad+1) )
    ALLOCATE( rrtm_ch4vmr(0:0,1:nzt_rad+1) )
    ALLOCATE( rrtm_n2ovmr(0:0,1:nzt_rad+1) )
    ALLOCATE( rrtm_o2vmr(0:0,1:nzt_rad+1) )
    ALLOCATE( rrtm_cfc11vmr(0:0,1:nzt_rad+1) )
    ALLOCATE( rrtm_cfc12vmr(0:0,1:nzt_rad+1) )
    ALLOCATE( rrtm_cfc22vmr(0:0,1:nzt_rad+1) )
    ALLOCATE( rrtm_ccl4vmr(0:0,1:nzt_rad+1) )
    ALLOCATE( rrtm_h2ovmr(0:0,1:nzt_rad+1) )

!
!-- Open file for reading
    nc_stat = NF90_OPEN( rrtm_input_file, NF90_NOWRITE, id )
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 549 )
!
!-- Inquire dimension ids and dimensions
    nc_stat = NF90_INQ_DIMID( id, 'Pressure', id_dim )
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
    nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, len = np)
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )

    nc_stat = NF90_INQ_DIMID( id, 'Absorber', id_dim )
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
    nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, len = nabs )
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )


!
!-- Allocate pressure, and trace gas arrays
    ALLOCATE( p_mls(1:np) )
    ALLOCATE( trace_mls(1:num_trace_gases,1:np) )
    ALLOCATE( trace_mls_tmp(1:nabs,1:np) )


    nc_stat = NF90_INQ_VARID( id, 'Pressure', id_var )
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
    nc_stat = NF90_GET_VAR( id, id_var, p_mls )
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )

    nc_stat = NF90_INQ_VARID( id, 'AbsorberAmountMLS', id_var )
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )
    nc_stat = NF90_GET_VAR( id, id_var, trace_mls_tmp )
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 550 )


!
!-- Write absorber amounts (mls) to trace_mls
    DO  n = 1, num_trace_gases
       CALL getAbsorberIndex( TRIM( trace_names(n) ), id_abs )

       trace_mls(n,1:np) = trace_mls_tmp(id_abs,1:np)

!
!--    Replace missing values by zero
       WHERE ( trace_mls(n,:) > 2.0_wp )
          trace_mls(n,:) = 0.0_wp
       END WHERE
    ENDDO

    DEALLOCATE( trace_mls_tmp )

    nc_stat = NF90_CLOSE( id )
    CALL netcdf_handle_error_rad( 'read_trace_gas_data', 551 )

!
!-- Add extra pressure level for calculations of the trace gas paths
    ALLOCATE( rrtm_play_tmp(1:nzt_rad+1) )
    ALLOCATE( rrtm_plev_tmp(1:nzt_rad+2) )

    rrtm_play_tmp(1:nzt_rad)   = rrtm_play(0,1:nzt_rad)
    rrtm_plev_tmp(1:nzt_rad+1) = rrtm_plev(0,1:nzt_rad+1)
    rrtm_play_tmp(nzt_rad+1)   = rrtm_plev(0,nzt_rad+1) * 0.5_wp
    rrtm_plev_tmp(nzt_rad+2)   = MIN( 1.0E-4_wp, 0.25_wp * rrtm_plev(0,nzt_rad+1) )

!
!-- Calculate trace gas path (zero at surface) with interpolation to the sounding levels
    ALLOCATE( trace_mls_path(1:nzt_rad+2,1:num_trace_gases) )

    trace_mls_path(nzb+1,:) = 0.0_wp

    DO  k = nzb+2, nzt_rad+2
       DO  m = 1, num_trace_gases
          trace_mls_path(k,m) = trace_mls_path(k-1,m)

!
!--       When the pressure level is higher than the trace gas pressure level, assume that
          IF ( rrtm_plev_tmp(k-1) > p_mls(1) )  THEN

             trace_mls_path(k,m) = trace_mls_path(k,m) + trace_mls(m,1) * ( rrtm_plev_tmp(k-1) -   &
                                   MAX( p_mls(1), rrtm_plev_tmp(k) ) ) / g
          ENDIF

!
!--       Integrate for each sounding level from the contributing p_mls levels
          DO  n = 2, np
!
!--          Limit p_mls so that it is within the model level
             p_mls_u = MIN( rrtm_plev_tmp(k-1), MAX( rrtm_plev_tmp(k), p_mls(n) ) )
             p_mls_l = MIN( rrtm_plev_tmp(k-1), MAX( rrtm_plev_tmp(k), p_mls(n-1) ) )

             IF ( p_mls_l > p_mls_u )  THEN

!
!--             Calculate weights for interpolation
                p_mls_m = 0.5_wp * ( p_mls_l + p_mls_u )
                p_wgt_u = ( p_mls(n-1) - p_mls_m ) / ( p_mls(n-1) - p_mls(n) )
                p_wgt_l = ( p_mls_m - p_mls(n) )   / ( p_mls(n-1) - p_mls(n) )

!
!--             Add level to trace gas path
                trace_mls_path(k,m) = trace_mls_path(k,m) + ( p_wgt_u * trace_mls(m,n) +           &
                                      p_wgt_l * trace_mls(m,n-1) ) * (p_mls_l - p_mls_u) / g
             ENDIF
          ENDDO

          IF ( rrtm_plev_tmp(k) < p_mls(np) )  THEN
             trace_mls_path(k,m) = trace_mls_path(k,m) + trace_mls(m,np) *                         &
                                   ( MIN( rrtm_plev_tmp(k-1), p_mls(np) ) - rrtm_plev_tmp(k) ) / g
          ENDIF
       ENDDO
    ENDDO


!
!-- Prepare trace gas path profiles
    ALLOCATE( trace_path_tmp(1:nzt_rad+1) )

    DO  m = 1, num_trace_gases

       trace_path_tmp(1:nzt_rad+1) = ( trace_mls_path(2:nzt_rad+2,m) -                             &
                                       trace_mls_path(1:nzt_rad+1,m) ) * g /                       &
                                     ( rrtm_plev_tmp(1:nzt_rad+1) - rrtm_plev_tmp(2:nzt_rad+2) )

!
!--    Save trace gas paths to the respective arrays
       SELECT CASE ( TRIM( trace_names(m) ) )

          CASE ( 'O3' )

             rrtm_o3vmr(0,:) = trace_path_tmp(:)

          CASE ( 'CO2' )

             rrtm_co2vmr(0,:) = trace_path_tmp(:)

          CASE ( 'CH4' )

             rrtm_ch4vmr(0,:) = trace_path_tmp(:)

          CASE ( 'N2O' )

             rrtm_n2ovmr(0,:) = trace_path_tmp(:)

          CASE ( 'O2' )

             rrtm_o2vmr(0,:) = trace_path_tmp(:)

          CASE ( 'CFC11' )

             rrtm_cfc11vmr(0,:) = trace_path_tmp(:)

          CASE ( 'CFC12' )

             rrtm_cfc12vmr(0,:) = trace_path_tmp(:)

          CASE ( 'CFC22' )

             rrtm_cfc22vmr(0,:) = trace_path_tmp(:)

          CASE ( 'CCL4' )

             rrtm_ccl4vmr(0,:) = trace_path_tmp(:)

          CASE ( 'H2O' )

             rrtm_h2ovmr(0,:) = trace_path_tmp(:)

          CASE DEFAULT

       END SELECT

    ENDDO

    DEALLOCATE( trace_path_tmp )
    DEALLOCATE( trace_mls_path )
    DEALLOCATE( rrtm_play_tmp )
    DEALLOCATE( rrtm_plev_tmp )
    DEALLOCATE( trace_mls )
    DEALLOCATE( p_mls )

 END SUBROUTINE read_trace_gas_data


 SUBROUTINE netcdf_handle_error_rad( routine_name, errno )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    USE NETCDF

    USE pegrid

    IMPLICIT NONE

    CHARACTER(LEN=7) ::  message_identifier  !<
    CHARACTER(LEN=*) ::  routine_name        !<

    INTEGER(iwp) ::  errno  !<

    IF ( nc_stat /= NF90_NOERR )  THEN

       WRITE( message_identifier, '(''NCF'',I4.4)' )  errno
       message_string = TRIM( NF90_STRERROR( nc_stat ) )

       CALL message( routine_name, message_identifier, 2, 2, 0, 6, 1 )

    ENDIF

 END SUBROUTINE netcdf_handle_error_rad
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate temperature tendency due to radiative cooling/heating.
!> Cache-optimized version.
!--------------------------------------------------------------------------------------------------!
#if defined( __rrtmg ) || defined( __tenstream )
 SUBROUTINE radiation_tendency_ij( i, j, tend )

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index along x-direction
    INTEGER(iwp) ::  j  !< grid index along y-direction
    INTEGER(iwp) ::  k  !< grid index along z-direction

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  tend  !< pt tendency term

    IF ( radiation_scheme == 'rrtmg'  .OR.  radiation_scheme == 'tenstream' )  THEN
!
!--    Calculate tendency based on heating rate
       DO  k = nzb+1, nzt+1
          tend(k,j,i) = tend(k,j,i) + ( rad_lw_hr(k,j,i) + rad_sw_hr(k,j,i) ) *                    &
                        d_exner(k) * d_seconds_hour *                                              &
                        MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
       ENDDO

    ENDIF

 END SUBROUTINE radiation_tendency_ij
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate temperature tendency due to radiative cooling/heating.
!> Vector-optimized version
!--------------------------------------------------------------------------------------------------!
#if defined( __rrtmg ) || defined( __tenstream )
 SUBROUTINE radiation_tendency( tend )

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index along x-direction
    INTEGER(iwp) ::  j  !< grid index along y-direction
    INTEGER(iwp) ::  k  !< grid index along z-direction

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  tend  !< pt tendency term

    IF ( radiation_scheme == 'rrtmg'  .OR.  radiation_scheme == 'tenstream' )  THEN
!
!--    Calculate tendency based on heating rate
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt+1
                tend(k,j,i) = tend(k,j,i) + ( rad_lw_hr(k,j,i) + rad_sw_hr(k,j,i) ) *              &
                              d_exner(k) * d_seconds_hour *                                        &
                              MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE radiation_tendency
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Radiative Transfer Model (RTM) version 3.0 for modelling of radiation interactions within urban
!> canopy or inside of surface layer in complex terrain. This subroutine calculates interaction of
!> the solar SW and LW radiation with urban and land surfaces and updates all surface heatfluxes.
!> It also calculates interactions of SW and LW radiation with resolved plant canopy and calculates
!> the corresponding plant canopy heat fluxes. The subroutine also models spatial and temporal
!> distribution of Mean Radiant Temperature (MRT). The resulting values are provided to other
!> PALM-4U modules (RRTMG, USM, LSM, PCM and BIO).
!>
!> The new version 3.0 was radically rewritten from version 1.0. The most significant changes
!> include new angular discretization scheme, redesigned and significantly optimized raytracing
!> scheme, new processes included in modelling (e.g. intetrations of LW radiation with PC),
!> integrated calculation of Mean Radiant Temperature (MRT), and improved and enhanced output and
!> debug capabilities. This new version significantly improves effectivity of the paralelization and
!> the scalability of the model and allows simulation of extensive domain with appropriate HPC
!> resources.
!>
!> More info about RTM v.1.0. see:
!> Resler et al., GMD. 2017, https://doi.org/10.5194/gmd-10-3635-2017
!> Info about RTM v. 3.0 see: Krc et al. 2021,  https://doi.org/10.5194/gmd-14-3095-2021
!> Maronga et al. 2020, https://doi.org/10.5194/gmd-13-1335-2020
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_interaction( at_restart )

    CHARACTER(*), OPTIONAL ::  at_restart !< string indicating that this is an initial call in a restart simulation

    IF ( radiation_interactions )  THEN

       IF ( loop_optimization == 'vector' )   THEN
          CALL radiation_interaction_vector( at_restart )
       ELSE
          CALL radiation_interaction_cache( at_restart )
       ENDIF

    ENDIF

 END SUBROUTINE radiation_interaction


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Cache-optimized version of the RTM.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_interaction_cache( at_restart )

    USE control_parameters,                                                                        &
        ONLY:  rotation_angle

    IMPLICIT NONE

    CHARACTER(*), OPTIONAL ::  at_restart  !< string indicating that this is an initial call in a restart simulation

    INTEGER(iwp) ::  i, j, k, kk, d, refstep, m, mm      !<
    INTEGER(iwp) ::  isurf, isurfsrc, isvf, icsf, ipcgb  !<
    INTEGER(iwp) ::  imrt, imrtf                         !<
    INTEGER(iwp) ::  isd                                 !< solar direction number
    INTEGER(iwp) ::  pc_box_dimshift                     !< transform for best accuracy

    REAL(wp) ::  asrc                                  !< area of source face
    REAL(wp) ::  grid_volume_inverse                   !< 1./(dx * dy * dz(1))
    REAL(wp) ::  pcrad                                 !< irradiance from plant canopy
    REAL(wp) ::  pc_box_area, pc_abs_frac, pc_abs_eff  !<
    REAL(wp) ::  temp                                  !< temporary variable for calculation

    REAL(wp), DIMENSION(3) ::  sunorig       !< grid rotated solar direction unit vector (zyx)
    REAL(wp), DIMENSION(3) ::  sunorig_grid  !< grid squashed solar direction unit vector (zyx)

    REAL(wp), DIMENSION(3,3) ::  mrot  !< grid rotation matrix (zyx)

    REAL(wp), DIMENSION(0:nsurf_type) ::  costheta  !< direct irradiance factor of solar angle

    REAL(wp), DIMENSION(3,0:nsurf_type) ::  vnorm  !< face direction normal vectors (zyx)

    REAL(wp), DIMENSION(nz_urban_b:nz_urban_t) ::  pchf_prep  !< precalculated factor for canopy temperature tendency


!
!-- Variables for coupling the radiation modle (e.g. RRTMG) and RTM
    REAL(wp) ::  area_norm         !< reference horizontal area of domain in all processor
    REAL(wp) ::  pabsswl           !< total absorbed SW radiation energy in local processor (W)
    REAL(wp) ::  pabssw            !< total absorbed SW radiation energy in all processors (W)
    REAL(wp) ::  pabslwl           !< total absorbed LW radiation energy in local processor (W)
    REAL(wp) ::  pabslw            !< total absorbed LW radiation energy in all processors (W)
    REAL(wp) ::  pemitlwl          !< total emitted LW radiation energy in all processors (W)
    REAL(wp) ::  pemitlw           !< total emitted LW radiation energy in all processors (W)
    REAL(wp) ::  pinswl            !< total received SW radiation energy in local processor (W)
    REAL(wp) ::  pinsw             !< total received SW radiation energy in all processor (W)
    REAL(wp) ::  pinlwl            !< total received LW radiation energy in local processor (W)
    REAL(wp) ::  pinlw             !< total received LW radiation energy in all processor (W)
    REAL(wp) ::  pabs_surf_lwdifl  !< total absorbed LW radiation in surfaces from sky in local processor (W)
    REAL(wp) ::  pabs_surf_lwdif   !< total absorbed LW radiation in surfaces from sky in all processors (W)
    REAL(wp) ::  pabs_pc_lwdifl    !< total absorbed LW radiation in plant canopy from sky in local processor (W)
    REAL(wp) ::  pabs_pc_lwdif     !< total absorbed LW radiation in plant canopy from sky in all processors (W)
!
!-- Rotation related variables
    REAL(wp) ::  cos_rot            !< cosine of rotation_angle
    REAL(wp) ::  sun_direct_factor  !< factor for direct normal radiation from direct horizontal
    REAL(wp) ::  sin_rot            !< sine of rotation_angle
    REAL(wp) ::  solar_azim         !< solar azimuth in rotated model coordinates
#if defined( __parallel )
    INTEGER(iwp) ::  surf_start_id                    !< id of first surface in current processor

    REAL(wp), DIMENSION(1:7) ::  combine_allreduce    !< dummy array used to combine several MPI_ALLREDUCE calls
    REAL(wp), DIMENSION(1:7) ::  combine_allreduce_l  !< dummy array used to combine several MPI_ALLREDUCE calls
#endif

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'radiation_interaction', time_since_reference_point
       CALL debug_message( debug_string, 'start' )
    ENDIF

    IF ( plant_canopy )  THEN
       grid_volume_inverse = 1.0_wp / ( dx * dy * dz(1) )
!
!--    pchf_prep is equal to 1 / (rho * c_p * T)
       pchf_prep(:) = r_d * exner(nz_urban_b:nz_urban_t) / ( c_p * hyp(nz_urban_b:nz_urban_t) )
    ENDIF

    sun_direction = .TRUE.
    CALL get_date_time( time_since_reference_point, day_of_year=day_of_year,                       &
                        second_of_day = second_of_day )
    CALL calc_zenith( day_of_year, second_of_day )      ! Required also for diffusion radiation

!
!-- Prepare rotated normal vectors and irradiance factor
    sin_rot = SIN( rotation_angle * pi / 180.0_wp )
    cos_rot = COS( rotation_angle * pi / 180.0_wp )
    vnorm(1,:) = kdir(:)
    vnorm(2,:) = jdir(:)
    vnorm(3,:) = idir(:)

    mrot(1,:) = (/ 1.0_wp,  0.0_wp,   0.0_wp /)
    mrot(2,:) = (/ 0.0_wp,  cos_rot, sin_rot /)
    mrot(3,:) = (/ 0.0_wp, -sin_rot, cos_rot /)
    sunorig = (/ cos_zenith, sun_dir_lat, sun_dir_lon /)
    sunorig = MATMUL( mrot, sunorig )
!
!-- Direct irradiance factor of solar angle, avoid negative value to prevent negative direct SW
!-- values
    DO  d = 0, nsurf_type
       costheta(d) = MAX( DOT_PRODUCT( sunorig, vnorm(:,d) ), 0.0_wp )
    ENDDO

    IF ( cos_zenith > 0 )  THEN
!
!--    Now we will "squash" the sunorig vector by grid box size in each dimension, so that this
!--    new direction vector will allow us to traverse the ray path within grid coordinates directly
       sunorig_grid = (/ sunorig(1) / dz(1), sunorig(2) / dy, sunorig(3) / dx /)
!       sunorig_grid = sunorig_grid / norm2(sunorig_grid)
       sunorig_grid = sunorig_grid / SQRT( SUM( sunorig_grid**2 ) )

       IF ( npcbl > 0 )  THEN
!
!--       Precompute effective box depth with prototype Leaf Area Density
          pc_box_dimshift = MAXLOC( ABS( sunorig ), 1) - 1
          CALL box_absorb( CSHIFT( (/ dz(1), dy, dx/), pc_box_dimshift ), 60, prototype_lad,       &
                           CSHIFT( ABS( sunorig ), pc_box_dimshift ), pc_box_area, pc_abs_frac )
          pc_box_area = pc_box_area * ABS( sunorig( pc_box_dimshift + 1 ) / sunorig(1) )
          pc_abs_eff = LOG( 1.0_wp - pc_abs_frac ) / prototype_lad
       ENDIF
    ENDIF
!
!-- Split downwelling shortwave radiation into a diffuse and a direct part. Note, if radiation
!-- scheme is RRTMG or diffuse radiation is externally prescribed, this is not required. Please
!-- note, in case of external radiation, the clear-sky model is applied during spinup, so that
!-- radiation needs to be split also in this case.
    IF ( radiation_scheme == 'constant'  .OR.  radiation_scheme == 'clear-sky'  .OR.               &
         ( radiation_scheme == 'external'  .AND.  .NOT.  rad_sw_in_dif_f%from_file )  .OR.         &
         ( radiation_scheme == 'external'  .AND.  time_since_reference_point < 0.0_wp ) )  THEN
       CALL radiation_calc_diffusion_radiation
    ENDIF

!
!-- First pass of radiation interaction:
!--  1) direct and diffuse irradiance
!--  2) thermal emissions
!
!-- Initialize relavant surface flux arrays and radiation energy sum
!-- Surface flux
    surfinswdir  = 0.0_wp
    surfins      = 0.0_wp
    surfinl      = 0.0_wp
    surfoutsl(:) = 0.0_wp
    surfoutll(:) = 0.0_wp
    IF ( nmrtbl > 0 )  THEN
       mrtinsw(:) = 0.0_wp
       mrtinlw(:) = 0.0_wp
    ENDIF
!
!-- Radiation energy sum
    pinlwl           = 0.0_wp
    pinswl           = 0.0_wp
    pemitlwl         = 0.0_wp
    pabsswl          = 0.0_wp
    pabslwl          = 0.0_wp
    pabs_surf_lwdifl = 0.0_wp
    pabs_pc_lwdifl   = 0.0_wp
!
!-- Set up thermal radiation from surfaces
    mm = 1
!-- Following code depends on the order of the execution. Do not parallelize by OpenMP!
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Urban-type surfaces
          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             surfoutll(mm) = SUM( surf_usm%frac(m,:) * surf_usm%emissivity(m,:) ) *                &
                             sigma_sb * surf_usm%pt_surface(m)**4
             albedo_surf(mm) = SUM( surf_usm%frac(m,:) * surf_usm%albedo(m,:) )
             emiss_surf(mm)  = SUM( surf_usm%frac(m,:) * surf_usm%emissivity(m,:) )
             mm = mm + 1
          ENDDO
!
!--       Land surfaces
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             surfoutll(mm) = SUM( surf_lsm%frac(m,:) * surf_lsm%emissivity(m,:) ) *                &
                             sigma_sb * surf_lsm%pt_surface(m)**4
             albedo_surf(mm) = SUM( surf_lsm%frac(m,:) * surf_lsm%albedo(m,:) )
             emiss_surf(mm)  = SUM( surf_lsm%frac(m,:) * surf_lsm%emissivity(m,:) )
             mm = mm + 1
          ENDDO
       ENDDO
    ENDDO

    IF ( trace_fluxes_above >= 0.0_wp )  THEN
       CALL radiation_print_debug_surf( 'surfoutll before initial pass', surfoutll )
       CALL radiation_print_debug_horz( 'rad_lw_in_diff before initial pass', rad_lw_in_diff )
       CALL radiation_print_debug_horz( 'rad_sw_in_diff before initial pass', rad_sw_in_diff )
       CALL radiation_print_debug_horz( 'rad_sw_in_dir before initial pass', rad_sw_in_dir )
    ENDIF

#if defined( __parallel )
!
!-- Sending out flux, surfoutll
    radx_send_surfinl = 0.0_wp

    surf_start_id = surfstart(myid)
    !$OMP PARALLEL DO PRIVATE (i)
!
!-- local + surf_start = global index, which is required in svf_in_send_buf
    DO  i = 1, nsend_radx
       radx_send(i) = surfoutll(isurf_send_radx(i) - surf_start_id)
    ENDDO

    CALL rtm_alltoallv( radx_send, disp_send_radx, surfoutl_recv, disp_recv_radx )
#endif

    IF ( surface_reflections )  THEN
       !$OMP PARALLEL DO PRIVATE (isvf, isurf, isurfsrc, temp) SCHEDULE (STATIC)
       DO  isvf = 1, nsvfl
          isurf    = svfsurf(1,isvf)
          isurfsrc = svfsurf(2,isvf)
!
!--       For surface-to-surface factors we calculate thermal radiation in 1st pass
          IF ( plant_lw_interact )  THEN
#if defined( __parallel )
             temp = svf(1,isvf) * svf(2,isvf) * surfoutl_recv(isurfsrc)
#else
             temp = svf(1,isvf) * svf(2,isvf) * surfoutll(isurfsrc)
#endif
          ELSE
#if defined( __parallel )
             temp = svf(1,isvf) * surfoutl_recv(isurfsrc)
#else
             temp = svf(1,isvf) * surfoutll(isurfsrc)
#endif
          ENDIF
          !$OMP ATOMIC
          surfinl(isurf) = surfinl(isurf) + temp
       ENDDO
    ENDIF
!
!-- Diffuse radiation using sky view factor
    !$OMP PARALLEL DO PRIVATE (i, j, d, isurf) REDUCTION(+:pinswl, pinlwl) SCHEDULE (STATIC)
    DO  isurf = 1, nsurfl
       j = surfl(iy,isurf)
       i = surfl(ix,isurf)
       d = surfl(id,isurf)
       surfinswdif(isurf) = rad_sw_in_diff(j,i) * skyvft(isurf)
!
!--    Update received SW energy for RTM coupling
       pinswl = pinswl + surfinswdif(isurf) * facearea(d)
       IF ( plant_lw_interact )  THEN
          surfinlwdif(isurf) = rad_lw_in_diff(j,i) * skyvft(isurf)
       ELSE
          surfinlwdif(isurf) = rad_lw_in_diff(j,i) * skyvf(isurf)
       ENDIF
!
!--    Update received LW energy for RTM coupling
       pinlwl = pinlwl + surfinlwdif(isurf) * facearea(d)
    ENDDO
!
!-- MRT diffuse irradiance
    !$OMP PARALLEL DO PRIVATE (i, j, imrt) SCHEDULE (STATIC)
    DO  imrt = 1, nmrtbl
       j = mrtbl(iy, imrt)
       i = mrtbl(ix, imrt)
       mrtinsw(imrt) = mrtskyt(imrt) * rad_sw_in_diff(j,i)
       mrtinlw(imrt) = mrtsky(imrt) * rad_lw_in_diff(j,i)
    ENDDO
!
!-- Direct radiation
    IF ( cos_zenith > 0 )  THEN
!
!--    To avoid numerical instability near horizon depending on what direct radiation is used
!--    (slightly different zenith angle, considering circumsolar etc.), we use a minimum value for
!--    cos_zenith
       sun_direct_factor = 1.0_wp / MAX( min_stable_coszen, cos_zenith )
!
!--    Identify solar direction vector (discretized number) (1)
       solar_azim = ATAN2( sun_dir_lon, sun_dir_lat ) * ( 180.0_wp / pi ) - rotation_angle
       j = FLOOR( ACOS( cos_zenith ) / pi * REAL( raytrace_discrete_elevs, KIND = wp ) )
       i = MODULO( NINT( solar_azim / 360.0_wp * REAL( raytrace_discrete_azims, KIND = wp )        &
                         - 0.5_wp, iwp ), raytrace_discrete_azims )
       isd = dsidir_rev(j, i)
!
!-- TODO: check if isd = -1 to report that this solar position is not precalculated
       !$OMP PARALLEL DO PRIVATE (i, j, d, isurf)  REDUCTION(+:pinswl) SCHEDULE (STATIC)
       DO  isurf = 1, nsurfl
          j = surfl(iy,isurf)
          i = surfl(ix,isurf)
          d = surfl(id,isurf)
          surfinswdir(isurf) = rad_sw_in_dir(j,i) * costheta(surfl(id, isurf)) *                   &
                               dsitrans(isurf, isd) * sun_direct_factor
!
!--       Update received SW energy for RTM coupling
          pinswl = pinswl + surfinswdir(isurf) * facearea(d)
       ENDDO
!
!--    MRT direct irradiance
       !$OMP PARALLEL DO PRIVATE (i, j, imrt) SCHEDULE (STATIC)
       DO  imrt = 1, nmrtbl
          j = mrtbl(iy,imrt)
          i = mrtbl(ix,imrt)
          mrtinsw(imrt) = mrtinsw(imrt) + mrtdsit(imrt, isd) * rad_sw_in_dir(j,i) *                &
                          sun_direct_factor * 0.25_wp  ! Normal to sphere
       ENDDO
    ENDIF
!
!-- MRT first pass thermal
    !$OMP PARALLEL DO PRIVATE (imrtf, imrt, isurfsrc, temp) SCHEDULE (STATIC)
    DO  imrtf = 1, nmrtf
       imrt     = mrtfsurf(1,imrtf)
       isurfsrc = mrtfsurf(2,imrtf)
#if defined( __parallel )
       temp = mrtf(imrtf) * surfoutl_recv(isurfsrc)
#else
       temp = mrtf(imrtf) * surfoutll(isurfsrc)
#endif
       !$OMP ATOMIC
       mrtinlw(imrt) = mrtinlw(imrt) + temp
    ENDDO
!
!-- Absorption in each local plant canopy grid box from the first atmospheric pass of radiation
    IF ( npcbl > 0 )  THEN
       pcbinswdir(:) = 0.0_wp
       pcbinswdif(:) = 0.0_wp
       pcbinlw(:)    = 0.0_wp
       pcinswdir(:)  = 0.0_wp
       pcinswdif(:)  = 0.0_wp

       !$OMP PARALLEL DO PRIVATE (icsf, ipcgb, i, j, k, kk, isurfsrc, pc_abs_frac, pcrad, asrc) &
#if defined( __parallel )
       !$OMP&   REDUCTION(+:pinswl, pinlwl, pabslwl, pemitlwl, pabs_pc_lwdifl, pcbinlw, radx_send_surfinl) SCHEDULE (STATIC)
#else
       !$OMP&   REDUCTION(+:pinswl, pinlwl, pabslwl, pemitlwl, pabs_pc_lwdifl, pcbinlw) SCHEDULE (STATIC)
#endif
       DO  icsf = 1, ncsfl
          ipcgb    = csfsurf(1,icsf)
          i        = pcbl(ix,ipcgb)
          j        = pcbl(iy,ipcgb)
          k        = pcbl(iz,ipcgb)
          kk       = k - topo_top_ind(j,i,0)  ! lad arrays are defined flat
          isurfsrc = csfsurf(2, icsf)

          IF ( isurfsrc == -1 )  THEN
!
!--          Diffuse radiation from sky (in W)
             pcbinswdif(ipcgb) = csf(1,icsf) * rad_sw_in_diff(j,i)
!
!--          Calculate the received diffuse radiation from sky for biogenic NMVOC emission model
!--          pcbinswdif should be converted to W m-2
             pcinswdif(ipcgb) = pcbinswdif(ipcgb) / (lad_s(kk,j,i) * dx * dy * dz(1))
!
!--          Add to the sum of SW radiation energy
             pinswl = pinswl + pcbinswdif(ipcgb)
!
!--          Convert diffuse radiation from sky to Wm-3
             pcbinswdif(ipcgb) = pcbinswdif(ipcgb) * grid_volume_inverse
!
!--          Absorbed diffuse LW radiation from sky minus emitted to sky
             IF ( plant_lw_interact )  THEN
                pcbinlw(ipcgb) = csf(1,icsf) * ( rad_lw_in_diff(j,i) - sigma_sb *                  &
                                               ( pt(k,j,i) * exner(k) )**4 ) *                     &
                                               grid_volume_inverse

                pinlwl = pinlwl + csf(1,icsf) * rad_lw_in_diff(j,i)
                pabslwl = pabslwl + csf(1,icsf) * rad_lw_in_diff(j,i)
                pemitlwl = pemitlwl + csf(1,icsf) * sigma_sb * ( pt(k,j,i) * exner(k) )**4
                pabs_pc_lwdifl = pabs_pc_lwdifl + csf(1,icsf) * rad_lw_in_diff(j,i)
             ENDIF
!
!--          Direct solar radiation
             IF ( cos_zenith > 0 )  THEN
!
!--             Estimate directed box absorption
                pc_abs_frac = 1.0_wp - exp( pc_abs_eff * lad_s(kk,j,i) )
!
!--             isd has already been established, see (1)
                pcbinswdir(ipcgb) = rad_sw_in_dir(j,i) * pc_box_area * pc_abs_frac *               &
                                    dsitransc(ipcgb,isd)
!
!--             Received direct sw radiation
                pcinswdir(ipcgb) = rad_sw_in_dir(j,i) * dsitransc(ipcgb,isd)
!
!--             Add to the sum of SW radiation energy
                pinswl = pinswl + pcbinswdir(ipcgb)
!
!--             Convert direct radiation from sky to Wm-3
                pcbinswdir(ipcgb) = pcbinswdir(ipcgb) * grid_volume_inverse
             ENDIF
          ELSE
             IF ( plant_lw_interact )  THEN
!
!--             Thermal emission from plan canopy towards respective face
                pcrad = sigma_sb * ( pt(k,j,i) * exner(k) )**4 * csf(1,icsf)
#if defined( __parallel )
                radx_send_surfinl(isurfsrc) = radx_send_surfinl(isurfsrc) + pcrad
#else
                surfinl(isurfsrc) = surfinl(isurfsrc) + pcrad
#endif
!
!--             Remove the flux above + absorb LW from first pass from surfaces
#if defined( __parallel )
                asrc = facearea(surf(id, isurf_recv_radx(isurfsrc)))
!
!--             Second line: Absorb from first pass surf emit - remove emitted heatflux
                pcbinlw(ipcgb) = pcbinlw(ipcgb)                                                    &
                                 + ( csf(1,icsf) * surfoutl_recv(isurfsrc) - pcrad )               &
                                 * asrc * grid_volume_inverse
                pabslwl = pabslwl + csf(1,icsf) * surfoutl_recv(isurfsrc) * asrc
#else
                asrc = facearea(surf(id, isurfsrc))
!
!--             Second line: Absorb from first pass surf emit - remove emitted heatflux
                pcbinlw(ipcgb) = pcbinlw(ipcgb)                                                    &
                                 + ( csf(1,icsf) * surfoutll(isurfsrc) - pcrad )                   &
                                 * asrc * grid_volume_inverse
                pabslwl = pabslwl + csf(1,icsf) * surfoutll(isurfsrc) * asrc
#endif
                pemitlwl = pemitlwl + pcrad * asrc
             ENDIF
          ENDIF
       ENDDO

       pcbinsw(:) = pcbinswdir(:) + pcbinswdif(:)
       pcinsw(:)  = pcinswdir(:)  + pcinswdif(:)
    ENDIF

    IF ( trace_fluxes_above >= 0.0_wp )  THEN
       CALL radiation_print_debug_surf( 'surfinl after initial pass', surfinl )
       CALL radiation_print_debug_surf( 'surfinlwdif after initial pass', surfinlwdif )
       CALL radiation_print_debug_surf( 'surfinswdif after initial pass', surfinswdif )
       CALL radiation_print_debug_surf( 'surfinswdir after initial pass', surfinswdir )
       IF ( npcbl > 0 )  THEN
          CALL radiation_print_debug_pcb( 'pcbinlw after initial pass', pcbinlw )
          CALL radiation_print_debug_pcb( 'pcbinswdif after initial pass', pcbinswdif )
          CALL radiation_print_debug_pcb( 'pcbinswdir after initial pass', pcbinswdir )
          CALL radiation_print_debug_pcb( 'pcinswdif after initial pass', pcinswdif )
          CALL radiation_print_debug_pcb( 'pcinswdir after initial pass', pcinswdir )
       ENDIF
    ENDIF

    IF ( plant_lw_interact )  THEN
!
!--    Exchange incoming lw radiation from plant canopy
#if defined( __parallel )
       CALL rtm_alltoallv( radx_send_surfinl, disp_recv_radx, surfinl_recv, disp_send_radx )
       !$OMP PARALLEL DO PRIVATE (i, isurf)
       DO  i = 1, nsend_radx
          isurf = isurf_send_radx(i) - surf_start_id
          surfinl(isurf) = surfinl(isurf) + surfinl_recv(i)
       ENDDO
#endif
    ENDIF

    IF ( trace_fluxes_above >= 0.0_wp )  THEN
       CALL radiation_print_debug_surf( 'surfinl after PC emiss', surfinl )
    ENDIF

    surfins     = surfinswdir + surfinswdif
    surfinl     = surfinl + surfinlwdif
    surfinsw    = surfins
    surfinlw    = surfinl
    surfoutsw   = 0.0_wp
    surfoutlw   = surfoutll
    surfemitlwl = surfoutll

    IF ( .NOT.  surface_reflections )  THEN
!
!--    Set nrefsteps to 0 to disable reflections
       nrefsteps = 0
       surfoutsl = albedo_surf * surfins
       surfoutll = ( 1.0_wp - emiss_surf ) * surfinl
       surfoutsw = surfoutsw + surfoutsl
       surfoutlw = surfoutlw + surfoutll
    ENDIF
!
!-- Next passes of radiation interactions: Radiation reflections
    DO  refstep = 1, nrefsteps

       surfoutsl = albedo_surf * surfins
!
!--    For non-transparent surfaces, longwave albedo is 1 - emissivity
       surfoutll = ( 1.0_wp - emiss_surf ) * surfinl

       IF ( trace_fluxes_above >= 0.0_wp )  THEN
          CALL radiation_print_debug_surf( 'surfoutll before reflective pass', surfoutll, refstep )
          CALL radiation_print_debug_surf( 'surfoutsl before reflective pass', surfoutsl, refstep )
       ENDIF

#if defined( __parallel )
!
!--    Sending out flux, surfoutll (local + surf_start = global index)
       !$OMP PARALLEL DO PRIVATE (i)
       DO  i = 1, nsend_radx
          radx_send(i) = surfoutll(isurf_send_radx(i) - surf_start_id)
       ENDDO

       CALL rtm_alltoallv( radx_send, disp_send_radx, surfoutl_recv, disp_recv_radx )
!
!--    Sending out flux, surfoutsl
       !$OMP PARALLEL DO PRIVATE (i)
       DO  i = 1, nsend_radx
          radx_send(i) = surfoutsl(isurf_send_radx(i) - surf_start_id)
       ENDDO

       CALL rtm_alltoallv( radx_send, disp_send_radx, surfouts_recv, disp_recv_radx )

#endif
!
!--    Reset for the input from next reflective pass
       surfins = 0.0_wp
       surfinl = 0.0_wp
!
!--    Reflected radiation
       !$OMP PARALLEL DO PRIVATE (isvf, isurf, isurfsrc) REDUCTION(+:surfins) REDUCTION(+:surfinl) SCHEDULE (STATIC)
       DO  isvf = 1, nsvfl
          isurf = svfsurf(1,isvf)
          isurfsrc = svfsurf(2,isvf)
#if defined( __parallel )
          surfins(isurf) = surfins(isurf) + svf(1,isvf) * svf(2,isvf) * surfouts_recv(isurfsrc)
          IF ( plant_lw_interact )  THEN
             surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * svf(2,isvf) * surfoutl_recv(isurfsrc)
          ELSE
             surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutl_recv(isurfsrc)
          ENDIF
#else
          surfins(isurf) = surfins(isurf) + svf(1,isvf) * svf(2,isvf) * surfoutsl(isurfsrc)
          IF ( plant_lw_interact )  THEN
             surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * svf(2,isvf) * surfoutll(isurfsrc)
          ELSE
             surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutll(isurfsrc)
          ENDIF
#endif
       ENDDO
!
!--    NOTE: PC absorbtion and MRT from reflected can both be done at once after all reflections
!--    if we do one more MPI_ALLGATHERV on surfout.
!--    Advantage: less local computation. Disadvantage: one more collective MPI call.
!
!--    Radiation absorbed by plant canopy
       !$OMP PARALLEL DO PRIVATE (icsf, ipcgb, isurfsrc, asrc, temp) SCHEDULE (STATIC)
       DO  icsf = 1, ncsfl
          ipcgb = csfsurf(1,icsf)
          isurfsrc = csfsurf(2,icsf)
          IF ( isurfsrc == -1 )  CYCLE                    ! sky->face only in 1st pass, not here
!
!--       Calculate source surface area. If the `surf' array is removed before timestepping starts
!--       (future version), then asrc must be stored within `csf'
#if defined( __parallel )
          asrc = facearea(surf(id, isurf_recv_radx(isurfsrc)))
          temp = csf(1,icsf) * surfouts_recv(isurfsrc) * asrc * grid_volume_inverse
#else
          asrc = facearea(surf(id, isurfsrc))
          temp = csf(1,icsf) * surfoutsl(isurfsrc) * asrc * grid_volume_inverse
#endif
          !$OMP ATOMIC
          pcbinsw(ipcgb) = pcbinsw(ipcgb) + temp
          IF ( plant_lw_interact )  THEN
#if defined( __parallel )
             temp = csf(1,icsf) * surfoutl_recv(isurfsrc) * asrc * grid_volume_inverse
#else
             temp = csf(1,icsf) * surfoutll(isurfsrc) * asrc * grid_volume_inverse
#endif
             !$OMP ATOMIC
             pcbinlw(ipcgb) = pcbinlw(ipcgb) + temp
          ENDIF
       ENDDO
!
!--    MRT reflected
       !$OMP PARALLEL DO PRIVATE (imrtf, imrt, isurfsrc, temp) SCHEDULE (STATIC)
       DO  imrtf = 1, nmrtf
          imrt     = mrtfsurf(1,imrtf)
          isurfsrc = mrtfsurf(2,imrtf)
#if defined( __parallel )
          temp = mrtft(imrtf) * surfouts_recv(isurfsrc)
#else
          temp = mrtft(imrtf) * surfoutsl(isurfsrc)
#endif
          !$OMP ATOMIC
          mrtinsw(imrt) = mrtinsw(imrt) + temp
#if defined( __parallel )
          temp = mrtf(imrtf) * surfoutl_recv(isurfsrc)
#else
          temp = mrtf(imrtf) * surfoutll(isurfsrc)
#endif
          !$OMP ATOMIC
          mrtinlw(imrt) = mrtinlw(imrt) + temp
       ENDDO

       IF ( trace_fluxes_above >= 0.0_wp )  THEN
          CALL radiation_print_debug_surf( 'surfinl after reflected pass', surfinl, refstep )
          CALL radiation_print_debug_surf( 'surfins after reflected pass', surfins, refstep )
          IF ( npcbl > 0 )  THEN
             CALL radiation_print_debug_pcb( 'pcbinlw after reflected pass', pcbinlw, refstep )
             CALL radiation_print_debug_pcb( 'pcbinsw after reflected pass', pcbinsw, refstep )
          ENDIF
       ENDIF

       surfinsw  = surfinsw  + surfins
       surfinlw  = surfinlw  + surfinl
       surfoutsw = surfoutsw + surfoutsl
       surfoutlw = surfoutlw + surfoutll

    ENDDO ! refstep

!
!-- Calculate black body MRT (after all reflections)
    IF ( nmrtbl > 0 )  THEN
       IF ( mrt_include_sw )  THEN
          mrt(:) = SQRT( SQRT( (mrtinsw(:) + mrtinlw(:) ) / sigma_sb ) )
       ELSE
          mrt(:) = SQRT( SQRT( mrtinlw(:) / sigma_sb ) )
       ENDIF
    ENDIF
!
!-- Calculate volumetric radiative fluxes
    IF ( radiation_volumetric_flux )  THEN
       swflux_vol(:,:,:) = 0.0_wp
!
!--    Add direct radiative flux above shadow
       IF ( cos_zenith > 0.0_wp )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                k = shadow_top(j,i,isd)
                rad_shade_h(j,i) = k
!
!--             The ratio of circle_area / sphere_surface is 1/4.
                swflux_vol(k+1:nz_urban_t,j,i) = rad_sw_in_dir(j,i) * sun_direct_factor * 0.25_wp
              ENDDO
           ENDDO
        ELSE
           rad_shade_h(:,:) = nz_urban_t
        ENDIF
        !TODO: add diffuse from skyvf + reflected from average surf out
    ENDIF

!
!-- Return if called at the beginning of a restart run (which is done just to calculate the mrt
!-- quantities).
    IF ( TRIM( initializing_actions ) == 'read_restart_data'  .AND.  PRESENT( at_restart ) )  THEN
       RETURN
    ENDIF

!
!-- Push heat flux absorbed by plant canopy to respective 3D arrays and add absorbed SW radiation
!-- energy for RTM coupling variables
    IF ( npcbl > 0 )  THEN
        pcm_sensiblerate(:,:,:) = 0.0_wp
        pcm_sensibleflux(:,:,:) = 0.0_wp
        !$OMP PARALLEL DO PRIVATE (ipcgb, i, j, k, kk) REDUCTION(+:pabsswl) SCHEDULE (STATIC)
        DO  ipcgb = 1, npcbl
           j = pcbl(iy,ipcgb)
           i = pcbl(ix,ipcgb)
           k = pcbl(iz,ipcgb)
!
!--        Following expression equals former kk = k - nzb_s_inner(j,i)
           kk = k - topo_top_ind(j,i,0)  ! lad arrays are defined flat
!
!--        All available energy heats air, latent flux considered below
           pcm_sensiblerate(kk,j,i) = ( pcbinsw(ipcgb) + pcbinlw(ipcgb) ) * pchf_prep(k) *         &
                                       pt(k,j,i) !-- = dT/dt
           pcm_sensibleflux(kk,j,i) = pcm_sensiblerate(kk,j,i) * c_p *                             &
                                        hyp(k) / ( r_d * pt(k,j,i) * exner(k) )
!
!--        Add the absorbed SW radiation energy by plant canopy
           pabsswl = pabsswl + pcbinsw(ipcgb) / grid_volume_inverse
        ENDDO

        IF ( humidity .AND. plant_canopy_transpiration )  THEN
!
!--        Calculation of plant canopy transpiration rate and correspondidng latent heat rate
           pcm_transpiration_rate(:,:,:) = 0.0_wp
           pcm_latentrate(:,:,:) = 0.0_wp
           pcm_latentflux(:,:,:) = 0.0_wp
           !$OMP PARALLEL DO PRIVATE (ipcgb, i, j, k, kk) SCHEDULE (STATIC)
           DO  ipcgb = 1, npcbl
              i = pcbl(ix,ipcgb)
              j = pcbl(iy,ipcgb)
              k = pcbl(iz,ipcgb)
              kk = k - topo_top_ind(j,i,0)  ! lad arrays are defined flat
              CALL pcm_calc_transpiration_rate( i, j, k, kk, pcbinsw(ipcgb), pcbinlw(ipcgb),       &
                                                pcm_transpiration_rate(kk,j,i),                    &
                                                pcm_latentrate(kk,j,i),                            &
                                                pcm_latentflux(kk,j,i))
!
!--           Remove latent flux from the available energy that heats air
              pcm_sensiblerate(kk,j,i) = pcm_sensiblerate(kk,j,i) - pcm_latentrate(kk,j,i)
              pcm_sensibleflux(kk,j,i) = pcm_sensibleflux(kk,j,i) - pcm_latentflux(kk,j,i)
           ENDDO
       ENDIF
    ENDIF

!
!-- Transfer radiation arrays required for energy balance to the respective data types and
!-- claculate relevant radiation model-RTM coupling terms
    mm = 1
!
!-- Following code depends on the order of the execution. Do not parallelize by OpenMP!
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Urban
          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             surf_usm%rad_sw_in(m)  = surfinsw(mm)
             surf_usm%rad_sw_out(m) = surfoutsw(mm)
             surf_usm%rad_sw_dir(m) = surfinswdir(mm)
             surf_usm%rad_sw_dif(m) = surfinswdif(mm)
             surf_usm%rad_sw_ref(m) = surfinsw(mm) - surfinswdir(mm) - surfinswdif(mm)
             surf_usm%rad_sw_res(m) = surfins(mm)
             surf_usm%rad_lw_in(m)  = surfinlw(mm)
             surf_usm%rad_lw_out(m) = surfoutlw(mm)
             surf_usm%rad_net(m)    = surfinsw(mm) - surfoutsw(mm) + surfinlw(mm) - surfoutlw(mm)
             surf_usm%rad_net_l(m)  = surf_usm%rad_net(m)
             surf_usm%rad_lw_dif(m) = surfinlwdif(mm)
             surf_usm%rad_lw_ref(m) = surfinlw(mm) - surfinlwdif(mm)
             surf_usm%rad_lw_res(m) = surfinl(mm)
             mm = mm + 1
          ENDDO
!
!--       Land
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             surf_lsm%rad_sw_in(m)  = surfinsw(mm)
             surf_lsm%rad_sw_out(m) = surfoutsw(mm)
             surf_lsm%rad_sw_dir(m) = surfinswdir(mm)
             surf_lsm%rad_sw_dif(m) = surfinswdif(mm)
             surf_lsm%rad_sw_ref(m) = surfinsw(mm) - surfinswdir(mm) - surfinswdif(mm)
             surf_lsm%rad_sw_res(m) = surfins(mm)
             surf_lsm%rad_lw_in(m)  = surfinlw(mm)
             surf_lsm%rad_lw_out(m) = surfoutlw(mm)
             surf_lsm%rad_net(m)    = surfinsw(mm) - surfoutsw(mm) + surfinlw(mm) - surfoutlw(mm)
             surf_lsm%rad_lw_dif(m) = surfinlwdif(mm)
             surf_lsm%rad_lw_ref(m) = surfinlw(mm) - surfinlwdif(mm)
             surf_lsm%rad_lw_res(m) = surfinl(mm)
             mm = mm + 1
          ENDDO
       ENDDO
    ENDDO

    !$OMP PARALLEL DO PRIVATE (i, d) REDUCTION(+:pabsswl, pabslwl, pemitlwl, pabs_surf_lwdifl) &
    !$OMP&            SCHEDULE (STATIC)
    DO  i = 1, nsurfl
       d = surfl(id, i)
!
!--    RTM coupling terms
!--    Sum of absorbed SW & LW radiation energy
       pabsswl = pabsswl + ( 1.0_wp - albedo_surf(i) ) * surfinsw(i) * facearea(d)
       pabslwl = pabslwl + emiss_surf(i) * surfinlw(i) * facearea(d)
!
!--    Sum of emitted LW radiation energy
       pemitlwl = pemitlwl + surfemitlwl(i) * facearea(d)
!
!--    emiss1
       pabs_surf_lwdifl = pabs_surf_lwdifl + emiss_surf(i) * facearea(d) * surfinlwdif(i)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (m) SCHEDULE (STATIC)
    DO  m = 1, surf_usm%ns
       surf_usm%surfhf(m) = surf_usm%rad_sw_in(m) + surf_usm%rad_lw_in(m) -                       &
                            surf_usm%rad_sw_out(m) - surf_usm%rad_lw_out(m)
    ENDDO
    !$OMP PARALLEL DO PRIVATE (m) SCHEDULE (STATIC)
    DO  m = 1, surf_lsm%ns
       surf_lsm%surfhf(m) = surf_lsm%rad_sw_in(m) + surf_lsm%rad_lw_in(m) -                       &
                            surf_lsm%rad_sw_out(m) - surf_lsm%rad_lw_out(m)
    ENDDO
!
!-- Gather all rad flux energy in all processors. In order to reduce the number of MPI calls
!-- (to reduce latencies), combine the required quantities in one array, sum it up, and
!-- subsequently re-distribute back to the respective quantities.
#if defined( __parallel )
    combine_allreduce_l(1) = pinswl
    combine_allreduce_l(2) = pinlwl
    combine_allreduce_l(3) = pabsswl
    combine_allreduce_l(4) = pabslwl
    combine_allreduce_l(5) = pemitlwl
    combine_allreduce_l(6) = pabs_surf_lwdifl
    combine_allreduce_l(7) = pabs_pc_lwdifl

    CALL MPI_ALLREDUCE( combine_allreduce_l, combine_allreduce, SIZE( combine_allreduce ),         &
                        MPI_REAL, MPI_SUM, comm2d, ierr )

    pinsw           = combine_allreduce(1)
    pinlw           = combine_allreduce(2)
    pabssw          = combine_allreduce(3)
    pabslw          = combine_allreduce(4)
    pemitlw         = combine_allreduce(5)
    pabs_surf_lwdif = combine_allreduce(6)
    pabs_pc_lwdif   = combine_allreduce(7)
#else
    pinsw           = pinswl
    pinlw           = pinlwl
    pabssw          = pabsswl
    pabslw          = pabslwl
    pemitlw         = pemitlwl
    pabs_surf_lwdif = pabs_surf_lwdifl
    pabs_pc_lwdif   = pabs_pc_lwdifl
#endif
!
!-- Calculate the effective radiation surface parameters based on the parameterizations in Krc et
!-- al. 2021.
!-- (1) Albedo Eq. * in Krc et al. 2021
    IF ( pinsw /= 0.0_wp )  albedo_eff = ( pinsw - pabssw ) / pinsw

!
!-- (2) Emmsivity Eq. * in Krc et al. 2021.
!-- emissivity_eff weighted average of surface and PC emissivity = absorbed LW
!-- in [surfaces + plant canopy] / pinlw.
    emissivity_eff = (pabs_surf_lwdif + pabs_pc_lwdif) / pinlw
!
!-- (3) Temperature
!-- effective horizontal area to account for the effect of vertical surfaces,
!-- Eq. * in Krc et al. 2021.
    area_norm = pinlw / rad_lw_in_diff(nyn,nxl)
!
!-- Temperature, Eq. * in Krc et al. 2021.
    t_rad_eff = SQRT( SQRT( ( pemitlw - pabslw + emissivity_eff * pinlw ) /                        &
                            ( emissivity_eff * sigma_sb * area_norm ) ) )

    IF ( radiation_volumetric_flux )  THEN
!
!--    Until here, swflux_vol contains only direct radiation. We add diffuse radiation for the
!--    sky-view factor part and a simple guess of average non-sky background radiation for the
!--    remaining part, which is calculated as rad_sw_in * albedo_surf
       DO  i = nxl, nxr
          DO  j = nys, nyn
             swflux_vol(:,j,i) = swflux_vol(:,j,i) +                                               &
                                 skyvf_vol(:,j,i) * rad_sw_in_diff(j,i) +                          &
                                 (1.0_wp - skyvf_vol(:,j,i)) * rad_sw_in(0,j,i) * albedo_eff
          ENDDO
       ENDDO
    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'radiation_interaction', 'end' )

 END SUBROUTINE radiation_interaction_cache


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Vector-optimized version of the RTM.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_interaction_vector( at_restart )

    USE control_parameters,                                                                        &
        ONLY:  rotation_angle

    IMPLICIT NONE

    CHARACTER(*), OPTIONAL ::  at_restart  !< string indicating that this is an initial call in a restart simulation

    INTEGER(iwp) ::  i, j, k, kk, d, refstep, m, mm      !<
    INTEGER(iwp) ::  ii                                  !< running index over block
    INTEGER(iwp) ::  jj                                  !< running index over elements with same index

    INTEGER(iwp) ::  isurf, isurfsrc, isvf, icsf, ipcgb  !<
    INTEGER(iwp) ::  imrt, imrtf                         !<
    INTEGER(iwp) ::  isd                                 !< solar direction number
    INTEGER(iwp) ::  pc_box_dimshift                     !< transform for best accuracy

    REAL(wp) ::  asrc                                  !< area of source face
    REAL(wp) ::  grid_volume_inverse                   !< 1./(dx * dy * dz(1))
    REAL(wp) ::  pcrad                                 !< irradiance from plant canopy
    REAL(wp) ::  pc_box_area, pc_abs_frac, pc_abs_eff  !<
    REAL(wp) ::  temp                                  !< temporary variable for calculation

    REAL(wp), DIMENSION(3) ::  sunorig       !< grid rotated solar direction unit vector (zyx)
    REAL(wp), DIMENSION(3) ::  sunorig_grid  !< grid squashed solar direction unit vector (zyx)

    REAL(wp), DIMENSION(3,3) ::  mrot  !< grid rotation matrix (zyx)

    REAL(wp), DIMENSION(0:nsurf_type) ::  costheta  !< direct irradiance factor of solar angle

    REAL(wp), DIMENSION(3,0:nsurf_type) ::  vnorm  !< face direction normal vectors (zyx)

    REAL(wp), DIMENSION(nz_urban_b:nz_urban_t) ::  pchf_prep  !< precalculated factor for canopy temperature tendency


!
!-- Variables for coupling the radiation modle (e.g. RRTMG) and RTM
    REAL(wp) ::  area_norm         !< reference horizontal area of domain in all processor
    REAL(wp) ::  pabsswl           !< total absorbed SW radiation energy in local processor (W)
    REAL(wp) ::  pabssw            !< total absorbed SW radiation energy in all processors (W)
    REAL(wp) ::  pabslwl           !< total absorbed LW radiation energy in local processor (W)
    REAL(wp) ::  pabslw            !< total absorbed LW radiation energy in all processors (W)
    REAL(wp) ::  pemitlwl          !< total emitted LW radiation energy in all processors (W)
    REAL(wp) ::  pemitlw           !< total emitted LW radiation energy in all processors (W)
    REAL(wp) ::  pinswl            !< total received SW radiation energy in local processor (W)
    REAL(wp) ::  pinsw             !< total received SW radiation energy in all processor (W)
    REAL(wp) ::  pinlwl            !< total received LW radiation energy in local processor (W)
    REAL(wp) ::  pinlw             !< total received LW radiation energy in all processor (W)
    REAL(wp) ::  pabs_surf_lwdifl  !< total absorbed LW radiation in surfaces from sky in local processor (W)
    REAL(wp) ::  pabs_surf_lwdif   !< total absorbed LW radiation in surfaces from sky in all processors (W)
    REAL(wp) ::  pabs_pc_lwdifl    !< total absorbed LW radiation in plant canopy from sky in local processor (W)
    REAL(wp) ::  pabs_pc_lwdif     !< total absorbed LW radiation in plant canopy from sky in all processors (W)
!
!-- Rotation related variables
    REAL(wp) ::  cos_rot            !< cosine of rotation_angle
    REAL(wp) ::  sun_direct_factor  !< factor for direct normal radiation from direct horizontal
    REAL(wp) ::  sin_rot            !< sine of rotation_angle
    REAL(wp) ::  solar_azim         !< solar azimuth in rotated model coordinates
#if defined( __parallel )
    INTEGER(iwp) ::  surf_start_id                    !< id of first surface in current processor

    REAL(wp), DIMENSION(1:7) ::  combine_allreduce    !< dummy array used to combine several MPI_ALLREDUCE calls
    REAL(wp), DIMENSION(1:7) ::  combine_allreduce_l  !< dummy array used to combine several MPI_ALLREDUCE calls
#endif

    REAL(wp), DIMENSION(ncsfl)   :: pinswl_v          !< temporary array to allow vectorization
    REAL(wp), DIMENSION(ncsfl)   :: pinlwl_v          !< temporary array to allow vectorization
    REAL(wp), DIMENSION(ncsfl)   :: pabslwl_v         !< temporary array to allow vectorization
    REAL(wp), DIMENSION(ncsfl)   :: pemitlwl_v        !< temporary array to allow vectorization
    REAL(wp), DIMENSION(ncsfl)   :: pabs_pc_lwdifl_v  !< temporary array to allow vectorization


    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'radiation_interaction', time_since_reference_point
       CALL debug_message( debug_string, 'start' )
    ENDIF

    IF ( plant_canopy )  THEN
       grid_volume_inverse = 1.0_wp / ( dx * dy * dz(1) )
!
!--    pchf_prep is equal to 1 / (rho * c_p * T).
       pchf_prep(:) = r_d * exner(nz_urban_b:nz_urban_t) / ( c_p * hyp(nz_urban_b:nz_urban_t) )
    ENDIF

    sun_direction = .TRUE.
    CALL get_date_time( time_since_reference_point, day_of_year=day_of_year,                       &
                        second_of_day = second_of_day )
    CALL calc_zenith( day_of_year, second_of_day )      ! Required also for diffusion radiation

!
!-- Prepare rotated normal vectors and irradiance factor.
    sin_rot = SIN( rotation_angle * pi / 180.0_wp )
    cos_rot = COS( rotation_angle * pi / 180.0_wp )
    vnorm(1,:) = kdir(:)
    vnorm(2,:) = jdir(:)
    vnorm(3,:) = idir(:)

    mrot(1,:) = (/ 1.0_wp,  0.0_wp,   0.0_wp /)
    mrot(2,:) = (/ 0.0_wp,  cos_rot, sin_rot /)
    mrot(3,:) = (/ 0.0_wp, -sin_rot, cos_rot /)
    sunorig = (/ cos_zenith, sun_dir_lat, sun_dir_lon /)
    sunorig = MATMUL( mrot, sunorig )
!
!-- Direct irradiance factor of solar angle, avoid negative value to prevent negative direct SW
!-- values.
    DO  d = 0, nsurf_type
        costheta(d) = MAX( DOT_PRODUCT( sunorig, vnorm(:,d) ), 0.0_wp )
    ENDDO

    IF ( cos_zenith > 0 )  THEN
!
!--    Now we will "squash" the sunorig vector by grid box size in each dimension, so that this
!--    new direction vector will allow us to traverse the ray path within grid coordinates.
!--    directly
       sunorig_grid = (/ sunorig(1) / dz(1), sunorig(2) / dy, sunorig(3) / dx /)
!       sunorig_grid = sunorig_grid / norm2(sunorig_grid)
       sunorig_grid = sunorig_grid / SQRT( SUM( sunorig_grid**2 ) )

       IF ( npcbl > 0 )  THEN
!
!--       Precompute effective box depth with prototype Leaf Area Density.
          pc_box_dimshift = MAXLOC( ABS( sunorig ), 1) - 1
          CALL box_absorb( CSHIFT( (/ dz(1), dy, dx/), pc_box_dimshift ), 60, prototype_lad,       &
                           CSHIFT( ABS( sunorig ), pc_box_dimshift ), pc_box_area, pc_abs_frac )
          pc_box_area = pc_box_area * ABS( sunorig( pc_box_dimshift + 1 ) / sunorig(1) )
          pc_abs_eff = LOG( 1.0_wp - pc_abs_frac ) / prototype_lad
       ENDIF
    ENDIF
!
!-- Split downwelling shortwave radiation into a diffuse and a direct part. Note, if radiation
!-- scheme is RRTMG or diffuse radiation is externally prescribed, this is not required. Please
!-- note, in case of external radiation, the clear-sky model is applied during spinup, so that
!-- radiation needs to be split also in this case.
    IF ( radiation_scheme == 'constant'  .OR.  radiation_scheme == 'clear-sky'  .OR.               &
         ( radiation_scheme == 'external'  .AND.  .NOT.  rad_sw_in_dif_f%from_file )  .OR.         &
         ( radiation_scheme == 'external'  .AND.  time_since_reference_point < 0.0_wp ) )          &
    THEN
       CALL radiation_calc_diffusion_radiation
    ENDIF

!
!-- First pass of radiation interaction:
!--  1) direct and diffuse irradiance
!--  2) thermal emissions
!
!-- Initialize relavant surface flux arrays and radiation energy sum.
!-- Surface flux
    surfinswdir  = 0.0_wp
    surfins      = 0.0_wp
    surfinl      = 0.0_wp
    surfoutsl(:) = 0.0_wp
    surfoutll(:) = 0.0_wp
    IF ( nmrtbl > 0 )  THEN
       mrtinsw(:) = 0.0_wp
       mrtinlw(:) = 0.0_wp
    ENDIF
!
!-- Radiation energy sum
    pinlwl           = 0.0_wp
    pinswl           = 0.0_wp
    pemitlwl         = 0.0_wp
    pabsswl          = 0.0_wp
    pabslwl          = 0.0_wp
    pabs_surf_lwdifl = 0.0_wp
    pabs_pc_lwdifl   = 0.0_wp
!
!-- Set up thermal radiation from surfaces.
    mm = 1
!-- Following code depends on the order of the execution. Do not parallelize by OpenMP!
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Urban-type surfaces
          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             surfoutll(mm) = SUM( surf_usm%frac(m,:) * surf_usm%emissivity(m,:) ) *                &
                             sigma_sb * surf_usm%pt_surface(m)**4
             albedo_surf(mm) = SUM( surf_usm%frac(m,:) * surf_usm%albedo(m,:) )
             emiss_surf(mm)  = SUM( surf_usm%frac(m,:) * surf_usm%emissivity(m,:) )
             mm = mm + 1
          ENDDO
!
!--       Land surfaces
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             surfoutll(mm) = SUM( surf_lsm%frac(m,:) * surf_lsm%emissivity(m,:) ) *                &
                             sigma_sb * surf_lsm%pt_surface(m)**4
             albedo_surf(mm) = SUM( surf_lsm%frac(m,:) * surf_lsm%albedo(m,:) )
             emiss_surf(mm)  = SUM( surf_lsm%frac(m,:) * surf_lsm%emissivity(m,:) )
             mm = mm + 1
          ENDDO
       ENDDO
    ENDDO

    IF ( trace_fluxes_above >= 0.0_wp )  THEN
       CALL radiation_print_debug_surf( 'surfoutll before initial pass', surfoutll )
       CALL radiation_print_debug_horz( 'rad_lw_in_diff before initial pass', rad_lw_in_diff )
       CALL radiation_print_debug_horz( 'rad_sw_in_diff before initial pass', rad_sw_in_diff )
       CALL radiation_print_debug_horz( 'rad_sw_in_dir before initial pass', rad_sw_in_dir )
    ENDIF

#if defined( __parallel )
!
!-- Sending out flux, surfoutll.
    radx_send_surfinl = 0.0_wp

    surf_start_id = surfstart(myid)
    !$OMP PARALLEL DO PRIVATE (i)
!
!-- local + surf_start = global index, which is required in svf_in_send_buf
    DO  i = 1, nsend_radx
       radx_send(i) = surfoutll(isurf_send_radx(i) - surf_start_id)
    ENDDO

    CALL rtm_alltoallv( radx_send,  disp_send_radx, surfoutl_recv, disp_recv_radx )
#endif

    IF ( surface_reflections )  THEN
       !$OMP PARALLEL DO PRIVATE (isvf, isurf, isurfsrc, temp) SCHEDULE (STATIC)
       DO  isvf = 1, nsvfl
          isurf    = svfsurf(1,isvf)
          isurfsrc = svfsurf(2,isvf)
!
!--       For surface-to-surface factors we calculate thermal radiation in 1st pass.
          IF ( plant_lw_interact )  THEN
#if defined( __parallel )
             temp = svf(1,isvf) * svf(2,isvf) * surfoutl_recv(isurfsrc)
#else
             temp = svf(1,isvf) * svf(2,isvf) * surfoutll(isurfsrc)
#endif
          ELSE
#if defined( __parallel )
             temp = svf(1,isvf) * surfoutl_recv(isurfsrc)
#else
             temp = svf(1,isvf) * surfoutll(isurfsrc)
#endif
          ENDIF
          !$OMP ATOMIC
          surfinl(isurf) = surfinl(isurf) + temp
       ENDDO
    ENDIF
!
!-- Diffuse radiation using sky view factor
    !$OMP PARALLEL DO PRIVATE (i, j, d, isurf) REDUCTION(+:pinswl, pinlwl) SCHEDULE (STATIC)
    DO  isurf = 1, nsurfl
       j = surfl(iy,isurf)
       i = surfl(ix,isurf)
       d = surfl(id,isurf)
       surfinswdif(isurf) = rad_sw_in_diff(j,i) * skyvft(isurf)
!
!--    Update received SW energy for RTM coupling
       pinswl = pinswl + surfinswdif(isurf) * facearea(d)
       IF ( plant_lw_interact )  THEN
          surfinlwdif(isurf) = rad_lw_in_diff(j,i) * skyvft(isurf)
       ELSE
          surfinlwdif(isurf) = rad_lw_in_diff(j,i) * skyvf(isurf)
       ENDIF
!
!--    Update received LW energy for RTM coupling.
       pinlwl = pinlwl + surfinlwdif(isurf) * facearea(d)
    ENDDO
!
!-- MRT diffuse irradiance.
    !$OMP PARALLEL DO PRIVATE (i, j, imrt) SCHEDULE (STATIC)
    DO  imrt = 1, nmrtbl
       j = mrtbl(iy, imrt)
       i = mrtbl(ix, imrt)
       mrtinsw(imrt) = mrtskyt(imrt) * rad_sw_in_diff(j,i)
       mrtinlw(imrt) = mrtsky(imrt) * rad_lw_in_diff(j,i)
    ENDDO
!
!-- Direct radiation
    IF ( cos_zenith > 0 )  THEN
!
!--    To avoid numerical instability near horizon depending on what direct radiation is used
!--    (slightly different zenith angle, considering circumsolar etc.), we use a minimum value for
!--    cos_zenith.
       sun_direct_factor = 1.0_wp / MAX( min_stable_coszen, cos_zenith )
!
!--    Identify solar direction vector (discretized number) (1).
       solar_azim = ATAN2( sun_dir_lon, sun_dir_lat ) * ( 180.0_wp / pi ) - rotation_angle
       j = FLOOR( ACOS( cos_zenith ) / pi * REAL( raytrace_discrete_elevs, KIND = wp ) )
       i = MODULO( NINT( solar_azim / 360.0_wp * REAL( raytrace_discrete_azims, KIND = wp )        &
                         - 0.5_wp, iwp ), raytrace_discrete_azims )
       isd = dsidir_rev(j, i)
!
!--    TODO: check if isd = -1 to report that this solar position is not precalculated
       !$OMP PARALLEL DO PRIVATE (i, j, d, isurf)  REDUCTION(+:pinswl) SCHEDULE (STATIC)
       DO  isurf = 1, nsurfl
          j = surfl(iy,isurf)
          i = surfl(ix,isurf)
          d = surfl(id,isurf)
          surfinswdir(isurf) = rad_sw_in_dir(j,i) * costheta(surfl(id, isurf)) *                   &
                               dsitrans(isurf, isd) * sun_direct_factor
!
!--       Update received SW energy for RTM coupling.
          pinswl = pinswl + surfinswdir(isurf) * facearea(d)
       ENDDO
!
!--    MRT direct irradiance.
       !$OMP PARALLEL DO PRIVATE (i, j, imrt) SCHEDULE (STATIC)
       DO  imrt = 1, nmrtbl
          j = mrtbl(iy,imrt)
          i = mrtbl(ix,imrt)
          mrtinsw(imrt) = mrtinsw(imrt) + mrtdsit(imrt, isd) * rad_sw_in_dir(j,i) *                &
                          sun_direct_factor * 0.25_wp  ! Normal to sphere
       ENDDO
    ENDIF
!
!-- MRT first pass thermal.
    !$OMP PARALLEL DO PRIVATE (imrtf, imrt, isurfsrc, temp) SCHEDULE (STATIC)
    DO  imrtf = 1, nmrtf
       imrt     = mrtfsurf(1,imrtf)
       isurfsrc = mrtfsurf(2,imrtf)
#if defined( __parallel )
       temp = mrtf(imrtf) * surfoutl_recv(isurfsrc)
#else
       temp = mrtf(imrtf) * surfoutll(isurfsrc)
#endif
       !$OMP ATOMIC
       mrtinlw(imrt) = mrtinlw(imrt) + temp
    ENDDO
!
!-- Absorption in each local plant canopy grid box from the first atmospheric pass of radiation.
    IF ( npcbl > 0 )  THEN
       pcbinswdir(:) = 0.0_wp
       pcbinswdif(:) = 0.0_wp
       pcbinlw(:)    = 0.0_wp
       pcinswdir(:)  = 0.0_wp
       pcinswdif(:)  = 0.0_wp

       pinswl_v         = 0.0
       pinlwl_v         = 0.0
       pabslwl_v        = 0.0
       pemitlwl_v       = 0.0
       pabs_pc_lwdifl_v = 0.0

       DO  ii = 1, nr_surf_bl
          isurfsrc = sorted_surf(1,surf_start(ii))
!
!--       Loop over values with the same isurfsrc value.
          DO  jj = surf_start(ii), surf_end(ii)
             ipcgb = sorted_surf(2,jj)
             i     = pcbl(ix,ipcgb)
             j     = pcbl(iy,ipcgb)
             k     = pcbl(iz,ipcgb)
             icsf  = sorted_surf(3,jj)
             IF ( plant_lw_interact )  THEN
!
!--             Thermal emission from plan canopy towards respective face.
                pcrad = sigma_sb * ( pt(k,j,i) * exner(k) )**4 * csf(1,icsf)
#if defined( __parallel )
                radx_send_surfinl(isurfsrc) = radx_send_surfinl(isurfsrc) + pcrad
#else
                surfinl(isurfsrc) = surfinl(isurfsrc) + pcrad
#endif
             ENDIF
          ENDDO
       ENDDO
!
!--    Run loop with ipcgp values with isurfrc == -1. This happens only once per ipcgb value,
!--    so that the loop can run run with ipcgp as loop index.
       !$OMP PARALLEL DO PRIVATE (icsf, ipcgb, i, j, k, kk, isurfsrc, pc_abs_frac) &
       !$OMP&   REDUCTION(+:pinswl_v, pinlwl_v, pabslwl_v, pemitlwl_v, pabs_pc_lwdifl_v, pcbinlw) SCHEDULE (STATIC)
       DO  ipcgb = 1, npcbl
          i        = pcbl(ix,ipcgb)
          j        = pcbl(iy,ipcgb)
          k        = pcbl(iz,ipcgb)
          kk       = k - topo_top_ind(j,i,0)  ! lad arrays are defined flat

          isurfsrc = csfsurf(2,ipcgb)
          icsf     = ipcgb
!
!--       Diffuse radiation from sky (in W).
          pcbinswdif(ipcgb) = csf(1,icsf) * rad_sw_in_diff(j,i)
!
!--       Calculate the received diffuse radiation from sky for biogenic NMVOC emission model
!--       pcbinswdif should be converted to W m-2.
          pcinswdif(ipcgb) = pcbinswdif(ipcgb) / ( lad_s(kk,j,i) * dx * dy * dz(1) )
!
!--       Add to the sum of SW radiation energy.
          pinswl_v(ipcgb) = pinswl_v(ipcgb) + pcbinswdif(ipcgb)
!
!--       Convert diffuse radiation from sky to Wm-3.
          pcbinswdif(ipcgb) = pcbinswdif(ipcgb) * grid_volume_inverse
!
!--       Absorbed diffuse LW radiation from sky minus emitted to sky.
          IF ( plant_lw_interact )  THEN
             pcbinlw(ipcgb) = csf(1,icsf) * ( rad_lw_in_diff(j,i) - sigma_sb *                     &
                                              ( pt(k,j,i) * exner(k) )**4 ) *                      &
                                            grid_volume_inverse

             pinlwl_v(ipcgb)   = pinlwl_v(ipcgb)   + csf(1,icsf) * rad_lw_in_diff(j,i)
             pabslwl_v(ipcgb)  = pabslwl_v(ipcgb)  + csf(1,icsf) * rad_lw_in_diff(j,i)
             pemitlwl_v(ipcgb) = pemitlwl_v(ipcgb) + csf(1,icsf) * sigma_sb *                      &
                                 ( pt(k,j,i) * exner(k) )**4
             pabs_pc_lwdifl_v(ipcgb) = pabs_pc_lwdifl_v(ipcgb) + csf(1,icsf) * rad_lw_in_diff(j,i)
          ENDIF
!
!--       Direct solar radiation.
          IF ( cos_zenith > 0 )  THEN
!
!--          Estimate directed box absorption.
             pc_abs_frac = 1.0_wp - EXP( pc_abs_eff * lad_s(kk,j,i) )
!
!--          isd has already been established, see (1).
             pcbinswdir(ipcgb) = rad_sw_in_dir(j,i) * pc_box_area * pc_abs_frac *                  &
                                 dsitransc(ipcgb,isd)
!
!--          Received direct sw radiation.
             pcinswdir(ipcgb) = rad_sw_in_dir(j,i) * dsitransc(ipcgb,isd)
!
!--          Add to the sum of SW radiation energy.
             pinswl_v(ipcgb) = pinswl_v(ipcgb) + pcbinswdir(ipcgb)
!
!--          Convert direct radiation from sky to Wm-3.
             pcbinswdir(ipcgb) = pcbinswdir(ipcgb) * grid_volume_inverse
          ENDIF
       ENDDO

       pinswl         = pinswl         + SUM( pinswl_v(1:npcbl) )
       pinlwl         = pinlwl         + SUM( pinlwl_v(1:npcbl) )
       pabslwl        = pabslwl        + SUM( pabslwl_v(1:npcbl) )
       pemitlwl       = pemitlwl       + SUM( pemitlwl_v(1:npcbl) )
       pabs_pc_lwdifl = pabs_pc_lwdifl + SUM( pabs_pc_lwdifl_v(1:npcbl) )

       !$OMP PARALLEL DO PRIVATE (icsf, ipcgb, i, j, k, kk, isurfsrc, pcrad, asrc) &
       !$OMP&   REDUCTION(+:pcbinlw) SCHEDULE (STATIC)
       DO  ii = 1, nr_blocks
          ipcgb = sorted_ipcgb(1,ipcgb_start(ii))
          i     = pcbl(ix,ipcgb)
          j     = pcbl(iy,ipcgb)
          k     = pcbl(iz,ipcgb)
          kk    = k - topo_top_ind(j,i,0)  ! lad arrays are defined flat

!
!--       Loop over values with the same ipcgb value.
          DO  jj = ipcgb_start(ii), ipcgb_end(ii)
             isurfsrc = sorted_ipcgb(2,jj)
             icsf     = sorted_ipcgb(3,jj)
             IF ( plant_lw_interact )  THEN
!
!--             Thermal emission from plan canopy towards respective face.
                pcrad = sigma_sb * ( pt(k,j,i) * exner(k) )**4 * csf(1,icsf)
!
!--             Remove the flux above + absorb LW from first pass from surfaces.
#if defined( __parallel )
                asrc = facearea(surf(id, isurf_recv_radx(isurfsrc)))

                pcbinlw(ipcgb) = pcbinlw(ipcgb) +                                                  &
                                 ( csf(1,icsf) * surfoutl_recv(isurfsrc) - pcrad ) *               &
                                 asrc * grid_volume_inverse
                pabslwl_v(jj) = csf(1,icsf) * surfoutl_recv(isurfsrc) * asrc
#else
                asrc = facearea(surf(id, isurfsrc))
                pcbinlw(ipcgb) = pcbinlw(ipcgb) +                                                  &
                                 ( csf(1,icsf) * surfoutll(isurfsrc) - pcrad ) *                   &
                                 asrc * grid_volume_inverse
                pabslwl_v(jj) = csf(1,icsf) * surfoutll(isurfsrc) * asrc
#endif
                pemitlwl_v(jj) = pcrad * asrc
             ENDIF
          ENDDO

          pemitlwl = pemitlwl + SUM( pemitlwl_v(ipcgb_start(ii):ipcgb_end(ii)) )
          pabslwl  = pabslwl  + SUM( pabslwl_v(ipcgb_start(ii):ipcgb_end(ii)) )

       ENDDO

       pcbinsw(:) = pcbinswdir(:) + pcbinswdif(:)
       pcinsw (:) = pcinswdir (:) + pcinswdif (:)
    ENDIF

    IF ( trace_fluxes_above >= 0.0_wp )  THEN
       CALL radiation_print_debug_surf( 'surfinl after initial pass', surfinl )
       CALL radiation_print_debug_surf( 'surfinlwdif after initial pass', surfinlwdif )
       CALL radiation_print_debug_surf( 'surfinswdif after initial pass', surfinswdif )
       CALL radiation_print_debug_surf( 'surfinswdir after initial pass', surfinswdir )
       IF ( npcbl > 0 )  THEN
          CALL radiation_print_debug_pcb( 'pcbinlw after initial pass', pcbinlw )
          CALL radiation_print_debug_pcb( 'pcbinswdif after initial pass', pcbinswdif )
          CALL radiation_print_debug_pcb( 'pcbinswdir after initial pass', pcbinswdir )
          CALL radiation_print_debug_pcb( 'pcinswdif after initial pass', pcinswdif )
          CALL radiation_print_debug_pcb( 'pcinswdir after initial pass', pcinswdir )
       ENDIF
    ENDIF

    IF ( plant_lw_interact )  THEN
!
!--    Exchange incoming lw radiation from plant canopy.
#if defined( __parallel )
       CALL rtm_alltoallv( radx_send_surfinl, disp_recv_radx, surfinl_recv, disp_send_radx )
       !$OMP PARALLEL DO PRIVATE (i, isurf)
       DO  i = 1, nsend_radx
          isurf = isurf_send_radx(i) - surf_start_id
          surfinl(isurf) = surfinl(isurf) + surfinl_recv(i)
       ENDDO
#endif
    ENDIF

    IF ( trace_fluxes_above >= 0.0_wp )  THEN
       CALL radiation_print_debug_surf( 'surfinl after PC emiss', surfinl )
    ENDIF

    surfins     = surfinswdir + surfinswdif
    surfinl     = surfinl + surfinlwdif
    surfinsw    = surfins
    surfinlw    = surfinl
    surfoutsw   = 0.0_wp
    surfoutlw   = surfoutll
    surfemitlwl = surfoutll

    IF ( .NOT. surface_reflections )  THEN
!
!--    Set nrefsteps to 0 to disable reflections.
       nrefsteps = 0
       surfoutsl = albedo_surf * surfins
       surfoutll = ( 1.0_wp - emiss_surf ) * surfinl
       surfoutsw = surfoutsw + surfoutsl
       surfoutlw = surfoutlw + surfoutll
    ENDIF
!
!-- Next passes of radiation interactions: Radiation reflections.
    DO  refstep = 1, nrefsteps

       surfoutsl = albedo_surf * surfins
!
!--    For non-transparent surfaces, longwave albedo is 1 - emissivity.
       surfoutll = ( 1.0_wp - emiss_surf ) * surfinl

       IF ( trace_fluxes_above >= 0.0_wp )  THEN
          CALL radiation_print_debug_surf( 'surfoutll before reflective pass', surfoutll, refstep )
          CALL radiation_print_debug_surf( 'surfoutsl before reflective pass', surfoutsl, refstep )
       ENDIF

#if defined( __parallel )
!
!--    Sending out flux, surfoutll (local + surf_start = global index)
       !$OMP PARALLEL DO PRIVATE (i)
       DO  i = 1, nsend_radx
          radx_send(i) = surfoutll(isurf_send_radx(i)-surf_start_id)
       ENDDO

       CALL rtm_alltoallv( radx_send, disp_send_radx, surfoutl_recv, disp_recv_radx )
!
!--    Sending out flux, surfoutsl.
       !$OMP PARALLEL DO PRIVATE (i)
       DO  i = 1, nsend_radx
          radx_send(i) = surfoutsl(isurf_send_radx(i)-surf_start_id)
       ENDDO

       CALL rtm_alltoallv( radx_send, disp_send_radx, surfouts_recv, disp_recv_radx )

#endif
!
!--    Reset for the input from next reflective pass.
       surfins = 0.0_wp
       surfinl = 0.0_wp
!
!--    Reflected radiation.
!      !$OMP PARALLEL DO PRIVATE (isvf, isurf, isurfsrc) REDUCTION(+:surfins) REDUCTION(+:surfinl) SCHEDULE (STATIC)
       DO  ii = 1, nr_surf_bl_2
          isurf = sorted_surf_2(1,surf_start_2(ii))

          DO  jj = surf_start_2(ii), surf_end_2(ii)
             isurfsrc = sorted_surf_2(2,jj)
             isvf     = sorted_surf_2(3,jj)

#if defined( __parallel )
             surfins(isurf) = surfins(isurf) + svf(1,isvf) * svf(2,isvf) * surfouts_recv(isurfsrc)

             IF ( plant_lw_interact )  THEN
                surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * svf(2,isvf) * surfoutl_recv(isurfsrc)
             ELSE
                surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutl_recv(isurfsrc)
             ENDIF
#else
             surfins(isurf) = surfins(isurf) + svf(1,isvf) * svf(2,isvf) * surfoutsl(isurfsrc)

             IF ( plant_lw_interact )  THEN
                surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * svf(2,isvf) * surfoutll(isurfsrc)
             ELSE
                surfinl(isurf) = surfinl(isurf) + svf(1,isvf) * surfoutll(isurfsrc)
             ENDIF
#endif
          ENDDO
       ENDDO

!--    NOTE: PC absorbtion and MRT from reflected can both be done at once after all reflections
!--    if we do one more MPI_ALLGATHERV on surfout.
!--    Advantage: less local computation. Disadvantage: one more collective MPI call.
!
!--    Radiation absorbed by plant canopy.
!      !$OMP PARALLEL DO PRIVATE (icsf, ipcgb, isurfsrc, asrc, temp) SCHEDULE (STATIC)
       DO  ii = 1, nr_blocks
          ipcgb = sorted_ipcgb(1,ipcgb_start(ii))

          DO jj = ipcgb_start(ii), ipcgb_end(ii)
             isurfsrc = sorted_ipcgb(2,jj)
             icsf     = sorted_ipcgb(3,jj)
!
!--          Calculate source surface area. If the `surf' array is removed before timestepping
!--          starts (future version), then asrc must be stored within `csf'.
#if defined( __parallel )
             asrc = facearea(surf(id, isurf_recv_radx(isurfsrc)))
             temp = csf(1,icsf) * surfouts_recv(isurfsrc) * asrc * grid_volume_inverse
#else
             asrc = facearea(surf(id, isurfsrc))
             temp = csf(1,icsf) * surfoutsl(isurfsrc) * asrc * grid_volume_inverse
#endif
             !$OMP ATOMIC
             pcbinsw(ipcgb) = pcbinsw(ipcgb) + temp

             IF ( plant_lw_interact )  THEN
#if defined( __parallel )
                temp = csf(1,icsf) * surfoutl_recv(isurfsrc) * asrc * grid_volume_inverse
#else
                temp = csf(1,icsf) * surfoutll(isurfsrc) * asrc * grid_volume_inverse
#endif
                !$OMP ATOMIC
                pcbinlw(ipcgb) = pcbinlw(ipcgb) + temp
             ENDIF
          ENDDO
       ENDDO

!
!--    MRT reflected.
       !$OMP PARALLEL DO PRIVATE (imrtf, imrt, isurfsrc, temp) SCHEDULE (STATIC)
       DO  imrtf = 1, nmrtf
          imrt     = mrtfsurf(1,imrtf)
          isurfsrc = mrtfsurf(2,imrtf)
#if defined( __parallel )
          temp = mrtft(imrtf) * surfouts_recv(isurfsrc)
#else
          temp = mrtft(imrtf) * surfoutsl(isurfsrc)
#endif
          !$OMP ATOMIC
          mrtinsw(imrt) = mrtinsw(imrt) + temp
#if defined( __parallel )
          temp = mrtf(imrtf) * surfoutl_recv(isurfsrc)
#else
          temp = mrtf(imrtf) * surfoutll(isurfsrc)
#endif
          !$OMP ATOMIC
          mrtinlw(imrt) = mrtinlw(imrt) + temp
       ENDDO

       IF ( trace_fluxes_above >= 0.0_wp )  THEN
          CALL radiation_print_debug_surf( 'surfinl after reflected pass', surfinl, refstep )
          CALL radiation_print_debug_surf( 'surfins after reflected pass', surfins, refstep )
          IF ( npcbl > 0 )  THEN
             CALL radiation_print_debug_pcb( 'pcbinlw after reflected pass', pcbinlw, refstep )
             CALL radiation_print_debug_pcb( 'pcbinsw after reflected pass', pcbinsw, refstep )
          ENDIF
       ENDIF

       surfinsw  = surfinsw  + surfins
       surfinlw  = surfinlw  + surfinl
       surfoutsw = surfoutsw + surfoutsl
       surfoutlw = surfoutlw + surfoutll

    ENDDO ! refstep

!
!-- Calculate black body MRT (after all reflections).
    IF ( nmrtbl > 0 )  THEN
       IF ( mrt_include_sw )  THEN
          mrt(:) = SQRT( SQRT( (mrtinsw(:) + mrtinlw(:) ) / sigma_sb ) )
       ELSE
          mrt(:) = SQRT( SQRT( mrtinlw(:) / sigma_sb ) )
       ENDIF
    ENDIF
!
!-- Calculate volumetric radiative fluxes.
    IF ( radiation_volumetric_flux )  THEN
       swflux_vol(:,:,:) = 0.0_wp
!
!--    Add direct radiative flux above shadow.
       IF ( cos_zenith > 0.0_wp )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                k = shadow_top(j,i,isd)
                rad_shade_h(j,i) = k
!
!--             The ratio of circle_area / sphere_surface is 1/4.
                swflux_vol(k+1:nz_urban_t,j,i) = rad_sw_in_dir(j,i) * sun_direct_factor * 0.25_wp
              ENDDO
           ENDDO
        ELSE
           rad_shade_h(:,:) = nz_urban_t
        ENDIF
        !TODO: add diffuse from skyvf + reflected from average surf out
    ENDIF

!
!-- Return if called at the beginning of a restart run (which is done just to calculate the mrt
!-- quantities).
    IF ( TRIM( initializing_actions ) == 'read_restart_data'  .AND.  PRESENT( at_restart ) )  THEN
       RETURN
    ENDIF
!
!-- Push heat flux absorbed by plant canopy to respective 3D arrays and add absorbed SW radiation
!-- energy for RTM coupling variables.
    IF ( npcbl > 0 )  THEN
       pcm_sensiblerate(:,:,:) = 0.0_wp
       pcm_sensibleflux(:,:,:) = 0.0_wp
       !$OMP PARALLEL DO PRIVATE (ipcgb, i, j, k, kk) REDUCTION(+:pabsswl) SCHEDULE (STATIC)
       DO  ipcgb = 1, npcbl
          j = pcbl(iy,ipcgb)
          i = pcbl(ix,ipcgb)
          k = pcbl(iz,ipcgb)
!
!--       Following expression equals former kk = k - nzb_s_inner(j,i)
          kk = k - topo_top_ind(j,i,0)  ! lad arrays are defined flat
!
!--       All available energy heats air, latent flux considered below
          pcm_sensiblerate(kk,j,i) = ( pcbinsw(ipcgb) + pcbinlw(ipcgb) ) * pchf_prep(k) *          &
                                       pt(k,j,i) !-- = dT/dt
          pcm_sensibleflux(kk,j,i) = pcm_sensiblerate(kk,j,i) * c_p *                              &
                                       hyp(k) / ( r_d * pt(k,j,i) * exner(k) )
!
!--       Add the absorbed SW radiation energy by plant canopy
          pabsswl = pabsswl + pcbinsw(ipcgb) / grid_volume_inverse
       ENDDO

       IF ( humidity .AND. plant_canopy_transpiration )  THEN
!
!--       Calculation of plant canopy transpiration rate and correspondidng latent heat rate.
          pcm_transpiration_rate(:,:,:) = 0.0_wp
          pcm_latentrate(:,:,:) = 0.0_wp
          pcm_latentflux(:,:,:) = 0.0_wp
          !$OMP PARALLEL DO PRIVATE (ipcgb, i, j, k, kk) SCHEDULE (STATIC)
          DO  ipcgb = 1, npcbl
             i = pcbl(ix,ipcgb)
             j = pcbl(iy,ipcgb)
             k = pcbl(iz,ipcgb)
             kk = k - topo_top_ind(j,i,0)  ! lad arrays are defined flat
             CALL pcm_calc_transpiration_rate( i, j, k, kk, pcbinsw(ipcgb), pcbinlw(ipcgb),        &
                                               pcm_transpiration_rate(kk,j,i),                     &
                                               pcm_latentrate(kk,j,i),                             &
                                               pcm_latentflux(kk,j,i))
!
!--          Remove latent flux from the available energy that heats air.
             pcm_sensiblerate(kk, j, i) = pcm_sensiblerate(kk, j, i) - pcm_latentrate(kk,j,i)
             pcm_sensibleflux(kk, j, i) = pcm_sensibleflux(kk, j, i) - pcm_latentflux(kk,j,i)
          ENDDO
       ENDIF
    ENDIF
!
!-- Transfer radiation arrays required for energy balance to the respective data types and
!-- claculate relevant radiation model-RTM coupling terms.
    mm = 1
!
!-- Following code depends on the order of the execution. Do not parallelize by OpenMP!
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Urban
          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             surf_usm%rad_sw_in(m)  = surfinsw(mm)
             surf_usm%rad_sw_out(m) = surfoutsw(mm)
             surf_usm%rad_sw_dir(m) = surfinswdir(mm)
             surf_usm%rad_sw_dif(m) = surfinswdif(mm)
             surf_usm%rad_sw_ref(m) = surfinsw(mm) - surfinswdir(mm) - surfinswdif(mm)
             surf_usm%rad_sw_res(m) = surfins(mm)
             surf_usm%rad_lw_in(m)  = surfinlw(mm)
             surf_usm%rad_lw_out(m) = surfoutlw(mm)
             surf_usm%rad_net(m)    = surfinsw(mm) - surfoutsw(mm) + surfinlw(mm) - surfoutlw(mm)
             surf_usm%rad_net_l(m)  = surf_usm%rad_net(m)
             surf_usm%rad_lw_dif(m) = surfinlwdif(mm)
             surf_usm%rad_lw_ref(m) = surfinlw(mm) - surfinlwdif(mm)
             surf_usm%rad_lw_res(m) = surfinl(mm)
             mm = mm + 1
          ENDDO
!
!--       Land.
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             surf_lsm%rad_sw_in(m)  = surfinsw(mm)
             surf_lsm%rad_sw_out(m) = surfoutsw(mm)
             surf_lsm%rad_sw_dir(m) = surfinswdir(mm)
             surf_lsm%rad_sw_dif(m) = surfinswdif(mm)
             surf_lsm%rad_sw_ref(m) = surfinsw(mm) - surfinswdir(mm) - surfinswdif(mm)
             surf_lsm%rad_sw_res(m) = surfins(mm)
             surf_lsm%rad_lw_in(m)  = surfinlw(mm)
             surf_lsm%rad_lw_out(m) = surfoutlw(mm)
             surf_lsm%rad_net(m)    = surfinsw(mm) - surfoutsw(mm) + surfinlw(mm) - surfoutlw(mm)
             surf_lsm%rad_lw_dif(m) = surfinlwdif(mm)
             surf_lsm%rad_lw_ref(m) = surfinlw(mm) - surfinlwdif(mm)
             surf_lsm%rad_lw_res(m) = surfinl(mm)
             mm = mm + 1
          ENDDO
       ENDDO
    ENDDO

    !$OMP PARALLEL DO PRIVATE (i, d) REDUCTION(+:pabsswl, pabslwl, pemitlwl, pabs_surf_lwdifl) &
    !$OMP&            SCHEDULE (STATIC)
    DO  i = 1, nsurfl
       d = surfl(id,i)
!
!--    RTM coupling terms.
!--    Sum of absorbed SW & LW radiation energy.
       pabsswl = pabsswl + ( 1.0_wp - albedo_surf(i) ) * surfinsw(i) * facearea(d)
       pabslwl = pabslwl + emiss_surf(i) * surfinlw(i) * facearea(d)
!
!--    Sum of emitted LW radiation energy.
       pemitlwl = pemitlwl + surfemitlwl(i) * facearea(d)
!
!--    emiss1.
       pabs_surf_lwdifl = pabs_surf_lwdifl + emiss_surf(i) * facearea(d) * surfinlwdif(i)
    ENDDO

    !$OMP PARALLEL DO PRIVATE (m) SCHEDULE (STATIC)
    DO  m = 1, surf_usm%ns
       surf_usm%surfhf(m) = surf_usm%rad_sw_in(m) + surf_usm%rad_lw_in(m) -                        &
                            surf_usm%rad_sw_out(m) - surf_usm%rad_lw_out(m)
    ENDDO
    !$OMP PARALLEL DO PRIVATE (m) SCHEDULE (STATIC)
    DO  m = 1, surf_lsm%ns
       surf_lsm%surfhf(m) = surf_lsm%rad_sw_in(m) + surf_lsm%rad_lw_in(m) -                        &
                            surf_lsm%rad_sw_out(m) - surf_lsm%rad_lw_out(m)
    ENDDO
!
!-- Gather all rad flux energy in all processors. In order to reduce the number of MPI calls
!-- (to reduce latencies), combine the required quantities in one array, sum it up, and
!-- subsequently re-distribute back to the respective quantities.
#if defined( __parallel )
    combine_allreduce_l(1) = pinswl
    combine_allreduce_l(2) = pinlwl
    combine_allreduce_l(3) = pabsswl
    combine_allreduce_l(4) = pabslwl
    combine_allreduce_l(5) = pemitlwl
    combine_allreduce_l(6) = pabs_surf_lwdifl
    combine_allreduce_l(7) = pabs_pc_lwdifl

    CALL MPI_ALLREDUCE( combine_allreduce_l, combine_allreduce, SIZE( combine_allreduce ),         &
                        MPI_REAL, MPI_SUM, comm2d, ierr )

    pinsw           = combine_allreduce(1)
    pinlw           = combine_allreduce(2)
    pabssw          = combine_allreduce(3)
    pabslw          = combine_allreduce(4)
    pemitlw         = combine_allreduce(5)
    pabs_surf_lwdif = combine_allreduce(6)
    pabs_pc_lwdif   = combine_allreduce(7)
#else
    pinsw           = pinswl
    pinlw           = pinlwl
    pabssw          = pabsswl
    pabslw          = pabslwl
    pemitlw         = pemitlwl
    pabs_surf_lwdif = pabs_surf_lwdifl
    pabs_pc_lwdif   = pabs_pc_lwdifl
#endif
!
!-- Calculate the effective radiation surface parameters based on the parameterizations in Krc et
!-- al. 2021. See equations 24 - 35 in Krc et al. 2021.
    IF ( pinsw /= 0.0_wp )  albedo_eff = ( pinsw - pabssw ) / pinsw
!
!-- Emissivity_eff weighted average of surface and PC emissivity = absorbed LW
!-- in [surfaces + plant canopy] / pinlw.
    emissivity_eff = ( pabs_surf_lwdif + pabs_pc_lwdif ) / pinlw
!
!-- Eeffective horizontal area to account for the effect of vertical surfaces,
!-- Eq. 30 in Krc et al. 2021.
    area_norm = pinlw / rad_lw_in_diff(nyn,nxl)
!
!-- Effective radiative temperature.
    t_rad_eff = SQRT( SQRT( ( pemitlw - pabslw + emissivity_eff * pinlw ) /                        &
                            ( emissivity_eff * sigma_sb * area_norm ) ) )

    IF ( radiation_volumetric_flux )  THEN
!
!--    Until here, swflux_vol contains only direct radiation. We add diffuse radiation for the
!--    sky-view factor part and a simple guess of average non-sky background radiation for the
!--    remaining part, which is calculated as rad_sw_in * albedo_surf.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             swflux_vol(:,j,i) = swflux_vol(:,j,i) +                                               &
                               + skyvf_vol(:,j,i) * rad_sw_in_diff(j,i)                            &
                               + ( 1.0_wp - skyvf_vol(:,j,i) ) * rad_sw_in(0,j,i) * albedo_eff
          ENDDO
       ENDDO
    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'radiation_interaction', 'end' )

 END SUBROUTINE radiation_interaction_vector


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Specific initialization routine for the vector-optimized branch for the RTM.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_interaction_vector_prepare

    INTEGER(iwp) ::  ii             !< running index over number of sorted blocks
    INTEGER(iwp) ::  icsf           !< running index for canopy-sink factors
    INTEGER(iwp) ::  isvf           !< running index for surface-view factors
    INTEGER(iwp) ::  isurfsrc       !< dummy variable to check whether a surface is affected by plant-canopy
    INTEGER(iwp) ::  nr_block_val   !< number of values in ipcbg loop
    INTEGER(iwp) ::  nr_surf_val    !< running index for plant-canopy surface interaction
    INTEGER(iwp) ::  nr_surf_val_2  !< running index for surface-surface interaction


    IF ( ncsfl > 0 )  THEN

       ALLOCATE( ipcgb_start(ncsfl) )
       ALLOCATE( ipcgb_end(ncsfl) )
       ALLOCATE( sorted_ipcgb(3,ncsfl) )
!
!--    Compute additional index arrays used to vectorize radiation_interaction.
!--    Therefore, count and collect elements with isurfsrc /= -1.
       nr_block_val = 0
       DO  icsf = 1, ncsfl
          isurfsrc = csfsurf(2,icsf)
          IF ( isurfsrc /= -1 )  THEN
             nr_block_val = nr_block_val + 1
             sorted_ipcgb(1,nr_block_val) = csfsurf(1,icsf)
             sorted_ipcgb(2,nr_block_val) = csfsurf(2,icsf)
             sorted_ipcgb(3,nr_block_val) = icsf
          ENDIF
       ENDDO
!
!--    Sort array with ascending ipcgp values.
       CALL quicksort( sorted_ipcgb(:,1:nr_block_val) )
!
!--    Compute index spaces ipcgb_start and ipcgb_end assessing the same ipcgp value.
       nr_blocks = 1
       ipcgb_start(1) = 1
       DO  ii = 2, nr_block_val
          IF ( sorted_ipcgb(1,ii) /= sorted_ipcgb(1,ii-1) )  THEN
             nr_blocks = nr_blocks + 1
             ipcgb_start(nr_blocks) = ii
             ipcgb_end(nr_blocks-1) = ii - 1
          ENDIF
       ENDDO
       ipcgb_end(nr_blocks) = nr_block_val

       ALLOCATE( surf_start(ncsfl) )
       ALLOCATE( surf_end(ncsfl) )
       ALLOCATE( sorted_surf(3,ncsfl) )

       nr_surf_val = 0
       DO  icsf = 1, ncsfl
          isurfsrc = csfsurf(2,icsf)
          IF ( isurfsrc /= -1 )  THEN
             nr_surf_val = nr_surf_val + 1
             sorted_surf(1,nr_surf_val) = csfsurf(2,icsf)
             sorted_surf(2,nr_surf_val) = csfsurf(1,icsf)
             sorted_surf(3,nr_surf_val) = icsf
          ENDIF
       ENDDO
       CALL quicksort( sorted_surf(:,1:nr_surf_val) )

       nr_surf_bl = 1
       surf_start(1) = 1
       DO  ii = 2, nr_surf_val
          IF ( sorted_surf(1,ii) /= sorted_surf(1,ii-1) )  THEN
             nr_surf_bl = nr_surf_bl + 1
             surf_start(nr_surf_bl) = ii
             surf_end(nr_surf_bl-1) = ii - 1
          ENDIF
       ENDDO
       surf_end(nr_surf_bl) = nr_surf_val
    ENDIF

    IF ( nsvfl > 0 )   THEN
       ALLOCATE( surf_start_2(nsvfl) )
       ALLOCATE( surf_end_2(nsvfl) )
       ALLOCATE( sorted_surf_2(3,nsvfl) )

       nr_surf_val_2 = 0
       DO  isvf = 1, nsvfl
          nr_surf_val_2 = nr_surf_val_2 + 1
          sorted_surf_2(1,nr_surf_val_2) = svfsurf(1,isvf)
          sorted_surf_2(2,nr_surf_val_2) = svfsurf(2,isvf)
          sorted_surf_2(3,nr_surf_val_2) = isvf
       ENDDO
       CALL quicksort( sorted_surf_2 )

       nr_surf_bl_2 = 1
       surf_start_2(1) = 1
       DO  ii = 2, nr_surf_val_2
          IF ( sorted_surf_2(1,ii) /= sorted_surf_2(1,ii-1) )  THEN
             nr_surf_bl_2 = nr_surf_bl_2 + 1
             surf_start_2(nr_surf_bl_2) = ii
             surf_end_2(nr_surf_bl_2-1) = ii - 1
          ENDIF
       ENDDO
       surf_end_2(nr_surf_bl_2) = nr_surf_val_2
    ENDIF

!
!-- Check if number of values  (csfsurf(2,...) == -1) equals npcbl and if the the elements within
!-- (csfsurf(2,...) == -1) are contiguous at the lower indices of csfsurf.
    IF ( npcbl > 0 )  THEN
       IF ( COUNT( csfsurf(2,1:npcbl) == -1 ) /= npcbl  .AND.                                      &
            COUNT( csfsurf(2,:) == -1 ) /= npcbl )                                                 &
       THEN
          WRITE( message_string, * ) 'issue in sorting canopy sink factors: ',                     &
                                     COUNT( csfsurf(2,1:npcbl) == -1 ),                            &
                                     'does not match its pre-calculated number ', npcbl
          CALL message( 'radiation_interaction_vector_prepare', 'RAD0063', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

 END SUBROUTINE radiation_interaction_vector_prepare


!--------------------------------------------------------------------------------------------------!
!> Calculates radiation absorbed by box with given size and LAD.
!>
!> Simulates resol**2 rays (by equally spacing a bounding horizontal square conatining all possible
!> rays that would cross the box) and calculates average transparency per ray. Returns fraction of
!> absorbed radiation flux and area for which this fraction is effective.
!--------------------------------------------------------------------------------------------------!
 PURE SUBROUTINE box_absorb( boxsize, resol, dens, uvec, area, absorb )
    IMPLICIT NONE

    INTEGER(iwp) ::  i, j  !<

    INTEGER(iwp), INTENT(IN) ::  resol  !< No. of rays in x and y dimensions

    REAL(wp) ::  xshift, yshift,               &  !<
                 xmin, xmax, ymin, ymax,       &  !<
                 xorig, yorig,                 &  !<
                 dx1, dy1, dz1, dx2, dy2, dz2, &  !<
                 crdist,                       &  !<
                 transp                           !<

    REAL(wp), INTENT(IN) ::  dens  !< box density (e.g. Leaf Area Density)

    REAL(wp), INTENT(OUT) ::  area, &  !< horizontal area for flux absorbtion
                              absorb   !< fraction of absorbed flux

    REAL(wp), DIMENSION(3), INTENT(IN) ::  boxsize, &  !< z, y, x size of box in m
                                           uvec        !< z, y, x unit vector of incoming flux


    xshift = uvec(3) / uvec(1) * boxsize(1)
    xmin = MIN( 0.0_wp, - xshift )
    xmax = boxsize(3) + MAX( 0.0_wp, - xshift )
    yshift = uvec(2) / uvec(1) * boxsize(1)
    ymin = MIN( 0.0_wp, - yshift )
    ymax = boxsize(2) + MAX( 0.0_wp, - yshift )

    transp = 0.0_wp
    DO  i = 1, resol
       xorig = xmin + ( xmax - xmin ) * ( i - 0.5_wp ) / resol
       DO  j = 1, resol
          yorig = ymin + ( ymax - ymin ) * ( j - 0.5_wp ) / resol

          dz1 = 0.0_wp
          dz2 = boxsize(1) / uvec(1)

          IF ( uvec(2) > 0.0_wp )  THEN
             dy1 = - yorig                / uvec(2)  !< Crossing with y=0
             dy2 = ( boxsize(2) - yorig ) / uvec(2)  !< Crossing with y=boxsize(2)
          ELSE  ! uvec(2) == 0
             dy1 = - HUGE( 1.0_wp )
             dy2 = HUGE( 1.0_wp )
          ENDIF

          IF ( uvec(3) > 0.0_wp )  THEN
             dx1 = - xorig                / uvec(3)  !< Crossing with x=0
             dx2 = ( boxsize(3) - xorig ) / uvec(3)  !< Crossing with x=boxsize(3)
          ELSE  ! uvec(3) == 0
             dx1 = - HUGE( 1.0_wp )
             dx2 = HUGE( 1.0_wp )
          ENDIF

          crdist = MAX( 0.0_wp, ( MIN( dz2, dy2, dx2 ) - MAX( dz1, dy1, dx1 ) ) )
          transp = transp + EXP( - ext_coef * dens * crdist )
       ENDDO
    ENDDO
    transp = transp / resol**2
    area = ( boxsize(3) + xshift ) * ( boxsize(2) + yshift )
    absorb = 1.0_wp - transp

 END SUBROUTINE box_absorb


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print consecutive radiative extremes if requested to trace early radiation interaction
!> instabilities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_print_debug_surf( description, values, step )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  description   !<

    INTEGER(iwp), INTENT(IN), OPTIONAL ::  step    !<

    REAL(wp), DIMENSION(:), INTENT(IN) ::  values  !<

    CHARACTER(LEN=50)   ::  location      !<
    CHARACTER(LEN=1024) ::  debug_string  !<

    INTEGER ::  isurf  !<

    REAL(wp) ::  x  !<


    isurf = MAXLOC( values, DIM = 1 )
    x = values(isurf)
    IF ( x < trace_fluxes_above )  RETURN

    IF ( PRESENT( step ) )  THEN
       WRITE( location, '(A," #",I0)' ) description, step
    ELSE
       location = description
    ENDIF

    WRITE( debug_string, '("Maximum of ",A50," = ",F12.1," at coords i=",I4,", j=",I4,", ' //      &
                         'k=",I4,", d=",I1,". Alb=",F7.3,", emis=",F7.3)' )                        &
           location, x, surfl(ix,isurf), surfl(iy,isurf), surfl(iz,isurf), surfl(id,isurf),        &
           albedo_surf(isurf), emiss_surf(isurf)
    CALL debug_message( debug_string, 'info' )

 END SUBROUTINE radiation_print_debug_surf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @Todo: Missing Subroutine Description!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_print_debug_pcb( description, values, step )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  description   !<

    INTEGER(iwp), INTENT(IN), OPTIONAL ::  step    !<

    REAL(wp), DIMENSION(:), INTENT(IN) ::  values  !<

    CHARACTER(LEN=50)   ::  location      !<
    CHARACTER(LEN=1024) ::  debug_string  !<

    INTEGER ::  ipcb  !<

    REAL(wp) ::  x  !<

    IF ( npcbl <= 0 )  RETURN
    ipcb = MAXLOC( values, DIM = 1 )
    x = values(ipcb) / ( dx * dy * dz(1) )
    IF ( x < trace_fluxes_above )  RETURN

    IF ( PRESENT( step ) )  THEN
       WRITE( location, '(A," #",I0)' ) description, step
    ELSE
       location = description
    ENDIF

    WRITE( debug_string, '("Maximum of ",A50," = ",F12.1," at coords i=",I4,", j=",I4,", k=",I4)' )&
           location, x, pcbl(ix,ipcb), pcbl(iy,ipcb), pcbl(iz,ipcb)
    CALL debug_message( debug_string, 'info' )

 END SUBROUTINE radiation_print_debug_pcb


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @Todo: Missing Subroutine Description!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_print_debug_horz( description, values, step )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  description    !<

    INTEGER(iwp), INTENT(IN), OPTIONAL ::  step     !<

    REAL(wp), DIMENSION(:,:), INTENT(IN) ::  values  !<

    CHARACTER(LEN=50)   ::  location      !<
    CHARACTER(LEN=1024) ::  debug_string  !<

    INTEGER, DIMENSION(2) ::  ji  !<

    REAL(wp) ::  x  !<


    ji = MAXLOC( values )
    x = values( ji(1), ji(2) )
    IF ( x < trace_fluxes_above )  RETURN

    IF ( PRESENT( step ) )  THEN
       WRITE( location, '(A," #",I0)' ) description, step
    ELSE
       location = description
    ENDIF

    WRITE( debug_string, '("Maximum of ",A50," = ",F12.1," at coords i=",I4,", j=",I4)' )          &
           location, x, ji(2), ji(1)
    CALL debug_message( debug_string, 'info' )

 END SUBROUTINE radiation_print_debug_horz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine splits direct and diffusion dw radiation for RTM processing.
!> It sould not be called in case the radiation model already does it.
!> It follows Boland, Ridley & Brown (2008)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_calc_diffusion_radiation

    USE palm_date_time_mod,                                                                        &
        ONLY:  seconds_per_day

    IMPLICIT NONE

    INTEGER(iwp) ::  i              !< grid index x-direction
    INTEGER(iwp) ::  j              !< grid index y-direction
    INTEGER(iwp) ::  days_per_year  !< days in the current year

    REAL(wp), PARAMETER ::  lowest_solarUp = 0.1_wp  !< limit the sun elevation to protect stability of the calculation

    REAL(wp) ::  clearnessIndex     !< clearness index
    REAL(wp) ::  corrected_solarUp  !< corrected solar up radiation
    REAL(wp) ::  diff_frac          !< diffusion fraction of the radiation
    REAL(wp) ::  etr                !< extraterestrial radiation
    REAL(wp) ::  horizontalETR      !< horizontal extraterestrial radiation
    REAL(wp) ::  second_of_year     !< current second of the year
    REAL(wp) ::  year_angle         !< angle

!
!--  Calculate current day and time based on the initial values and simulation time
     CALL get_date_time( time_since_reference_point, second_of_year = second_of_year,              &
                         days_per_year = days_per_year    )
     year_angle = second_of_year / ( REAL( days_per_year, KIND = wp ) * seconds_per_day ) *        &
                  2.0_wp * pi

     etr = solar_constant * ( 1.00011_wp +  0.034221_wp * COS(year_angle) +                        &
                                            0.001280_wp * SIN(year_angle) +                        &
                                            0.000719_wp * COS(2.0_wp * year_angle) +               &
                                            0.000077_wp * SIN(2.0_wp * year_angle) )

!
!--  Under a very low angle, we keep extraterestrial radiation at the last small value, therefore
!--  the clearness index will be pushed towards 0 while keeping full continuity.
     IF ( cos_zenith <= lowest_solarUp )  THEN
         corrected_solarUp = lowest_solarUp
     ELSE
         corrected_solarUp = cos_zenith
     ENDIF

     horizontalETR = etr * corrected_solarUp

     DO  i = nxl, nxr
         DO  j = nys, nyn
            clearnessIndex = rad_sw_in(0,j,i) / horizontalETR
            diff_frac = 1.0_wp / ( 1.0_wp + EXP( -5.0033_wp + 8.6025_wp * clearnessIndex ) )
            rad_sw_in_diff(j,i) = rad_sw_in(0,j,i) * diff_frac
            rad_sw_in_dir(j,i)  = rad_sw_in(0,j,i) * ( 1.0_wp - diff_frac )
            rad_lw_in_diff(j,i) = rad_lw_in(0,j,i)
         ENDDO
     ENDDO

 END SUBROUTINE radiation_calc_diffusion_radiation


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine initializes structures needed for Radiative Transfer Model (RTM). This model
!> calculates transformation processes of the radiation inside urban and land canopy layer. The
!> module includes also the interaction of the radiation with the resolved plant canopy.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_interaction_init

    USE control_parameters,                                                                        &
        ONLY:  dz_stretch_level_start

    USE plant_canopy_model_mod,                                                                    &
        ONLY:  lad_s

    IMPLICIT NONE

    INTEGER(iwp) ::  facing               !< dummy argument for surface orientation
    INTEGER(iwp) ::  i, j, k, l, m, d     !<
    INTEGER(iwp) ::  icol                 !< flat column number (in (y,x) fortran order)
    INTEGER(iwp) ::  isurf, ipcgb, imrt   !<
    INTEGER(iwp) ::  k_topo               !< vertical index indicating topography top for given (j,i)
    INTEGER(iwp) ::  nzptl, nzubl, nzutl  !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nsurfs_nidx_surf      !< temporary array to hold nsurfs * nidx_surf
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  surfstart_nidx_surf   !< temporary array to hold surfstart * nidx_surf

!
!--  Precalculate face areas for different face directions using normal vector
     DO  d = 0, nsurf_type
        facearea(d) = 1.0_wp
        IF ( idir(d) == 0 ) facearea(d) = facearea(d) * dx
        IF ( jdir(d) == 0 ) facearea(d) = facearea(d) * dy
        IF ( kdir(d) == 0 ) facearea(d) = facearea(d) * dz(1)
     ENDDO
!
!-- Find nz_urban_b, nz_urban_t, nz_urban via topography top index.
!-- The following contruct finds the lowest / largest index for any upward-facing wall (see bit 12).
    nzubl = MINVAL( topo_top_ind(nys:nyn,nxl:nxr,0) )
    nzutl = MAXVAL( topo_top_ind(nys:nyn,nxl:nxr,5) )

    nzubl = MAX( nzubl, nzb )

    IF ( plant_canopy )  THEN
!
!--     Allocate needed arrays
        ALLOCATE( pct(nys:nyn,nxl:nxr) )
        ALLOCATE( pch(nys:nyn,nxl:nxr) )
!
!--     Calculate plant canopy height
        npcbl = 0
        pct   = 0
        pch   = 0
        DO  i = nxl, nxr
            DO  j = nys, nyn
!
!--            Find topography top index
               k_topo = topo_top_ind(j,i,0)

               DO  k = nzt+1, 1, -1
                  IF ( lad_s(k,j,i) > 0.0_wp )  THEN
!
!--                   We are at the top of the pcs
                      pct(j,i) = k + k_topo
                      pch(j,i) = k
                      npcbl = npcbl + COUNT( lad_s(1:k,j,i) > 0.0_wp )
                      EXIT
                  ENDIF
               ENDDO
            ENDDO
        ENDDO

        nzutl = MAX( nzutl, MAXVAL( pct ) )
        nzptl = MAXVAL( pct )

        prototype_lad = MAXVAL( lad_s ) * .9_wp  !< Better be *1.0 if lad is either 0 or maxval(lad) everywhere
        IF ( prototype_lad <= 0.0_wp )  prototype_lad = .3_wp
        !WRITE(message_string, '(a,f6.3)') 'Precomputing effective box optical ' &
        !    // 'depth using prototype leaf area density = ', prototype_lad
        !CALL message('radiation_interaction_init', 'RADxxxx', 0, 0, -1, 6, 0)
    ENDIF

    nzutl = MIN( nzutl + nzut_free, nzt )

#if defined( __parallel )
    CALL MPI_ALLREDUCE( nzubl, nz_urban_b, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
        WRITE( debug_string, * ) 'Error MPI_ALLREDUCE11:', ierr, nzubl, nz_urban_b
        CALL debug_message( debug_string, 'info' )
    ENDIF
    CALL MPI_ALLREDUCE( nzutl, nz_urban_t, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
        WRITE( debug_string, * ) 'Error MPI_ALLREDUCE12:', ierr, nzutl, nz_urban_t
        CALL debug_message( debug_string, 'info' )
    ENDIF
    CALL MPI_ALLREDUCE( nzptl, nz_plant_t, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
        WRITE( debug_string, * ) 'Error MPI_ALLREDUCE13:', ierr, nzptl, nz_plant_t
        CALL debug_message( debug_string, 'info' )
    ENDIF
#else
    nz_urban_b = nzubl
    nz_urban_t = nzutl
    nz_plant_t = nzptl
#endif
!
!-- Stretching (non-uniform grid spacing) is not considered in the radiation model. Therefore,
!-- vertical stretching has to be applied above the area where the parts of the radiation model
!-- which assume constant grid spacing are active. ABS (...) is required because the default value
!-- of dz_stretch_level_start is -9999999.9_wp (negative).
    IF ( ABS( dz_stretch_level_start(1) ) <= zw(nz_urban_t) )  THEN
       WRITE( message_string, * ) 'The lowest level where vertical stretching is applied have ' // &
                                  'to be greater than ', zw(nz_urban_t)
       CALL message( 'radiation_interaction_init', 'RAD0051', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Global number of urban and plant layers
    nz_urban = nz_urban_t - nz_urban_b + 1
    nz_plant = nz_plant_t - nz_urban_b + 1
!
!-- Allocate urban surfaces grid
!-- Calc number of surfaces in local proc
    IF ( debug_output )  CALL debug_message( 'calculation of indices for surfaces', 'info' )
!
!-- Number of horizontal surfaces including land- and roof surfaces in both USM and LSM. Note that
!-- all surface elements have been already counted in surface_mod.
    nsurfl = 0
    nsurfl = nsurfl + surf_usm%ns + surf_lsm%ns
!
!-- Fill gridpcbl and pcbl
    IF ( npcbl > 0 )  THEN
        ALLOCATE( pcbl(iz:ix, 1:npcbl) )
        ALLOCATE( gridpcbl(nz_urban_b:nz_plant_t,nys:nyn,nxl:nxr) )
        pcbl = -1
        gridpcbl(:,:,:) = 0
        ipcgb = 0
        DO  i = nxl, nxr
           DO  j = nys, nyn
!
!--           Find topography top index
              k_topo = topo_top_ind(j,i,0)

              DO  k = k_topo + 1, pct(j,i)
                 IF ( lad_s(k-k_topo,j,i) > 0.0_wp )  THEN
                    ipcgb = ipcgb + 1
                    gridpcbl(k,j,i) = ipcgb
                    pcbl(:,ipcgb) = (/ k, j, i /)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        ALLOCATE( pcbinsw(1:npcbl) )
        ALLOCATE( pcbinswdir(1:npcbl) )
        ALLOCATE( pcbinswdif(1:npcbl) )
        ALLOCATE( pcbinlw(1:npcbl) )
        ALLOCATE( pcinsw(1:npcbl) )
        ALLOCATE( pcinswdir(1:npcbl) )
        ALLOCATE( pcinswdif(1:npcbl) )
    ENDIF

!
!-- Allocate and calculate auxiliary indices needed for MPI exchanges
!-- Numbers of xy grid elements for individual PE and corresponding displacements
    ALLOCATE( nnxy(0:numprocs-1), nnxyd(0:numprocs-1) )
    k = 0
    DO  i = 0, npex-1
        DO  j = 0, npey-1
            nnxy(k) = (nxr_pe(i) - nxl_pe(i) + 1) * (nyn_pe(j) - nys_pe(j) + 1)
            k = k + 1
        ENDDO
    ENDDO
    nnxyd(0) = 0
    DO  i = 1, numprocs-1
        nnxyd(i) = nnxyd(i-1) + nnxy(i-1)
    ENDDO
!
!-- Indices of PE numbers along x a y axis
    ALLOCATE( ipx(0:nx), ipy(0:ny) )
    DO  i = 0, npex-1
        ipx(nxl_pe(i):nxr_pe(i)) = i
    ENDDO
    DO  j = 0, npey-1
        ipy(nys_pe(j):nyn_pe(j)) = j
    ENDDO
!
!-- Allocate and fill surfl and surfl_col_start. The ordering of local surfaces
!-- given by the following cycles must not be altered, certain file input
!-- routines may depend on it.
!
!-- We allocate the array as linear and then use a two-dimensional pointer
!-- into it, because some MPI implementations crash with 2D-allocated arrays.
    ALLOCATE( surfl_linear(nidx_surf*nsurfl) )
    surfl(1:nidx_surf,1:nsurfl) => surfl_linear(1:nidx_surf*nsurfl)
    ALLOCATE( surfl_col_start(0:nnx*nny-1) )
!
!-- Add horizontal and vertical surface elements (land and urban surfaces) ordered by x,y column
!-- (y most varying)
!-- TODO: remove the hard coding of l = 0 to l = idirection
    isurf = 0
    icol  = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Save column start
          surfl_col_start(icol) = isurf + 1
          icol = icol + 1

          DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             k = surf_usm%k(m)
             isurf = isurf + 1
!
!--          Set index indicating the facing. Note, this is only a preliminary work-around.
             IF ( surf_usm%upward(m)    )  facing = iup
             IF ( surf_usm%downward(m)  )  facing = idown
             IF ( surf_usm%northward(m) )  facing = inorth
             IF ( surf_usm%southward(m) )  facing = isouth
             IF ( surf_usm%eastward(m)  )  facing = ieast
             IF ( surf_usm%westward(m)  )  facing = iwest
             surfl(:,isurf) = (/facing,k,j,i/)
          ENDDO
          DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             k = surf_lsm%k(m)
             isurf = isurf + 1
!
!--          Set index indicating the facing. Note, this is only a preliminary work-around.
             IF ( surf_lsm%upward(m)    )  facing = iup
             IF ( surf_lsm%downward(m)  )  facing = idown
             IF ( surf_lsm%northward(m) )  facing = inorth
             IF ( surf_lsm%southward(m) )  facing = isouth
             IF ( surf_lsm%eastward(m)  )  facing = ieast
             IF ( surf_lsm%westward(m)  )  facing = iwest
             surfl(:,isurf) = (/facing,k,j,i/)
          ENDDO

       ENDDO
    ENDDO
!
!-- Add local MRT boxes for the specified number of levels
!-- !!!! NEEDS TO RETHINK AGAIN - With full 3D structure, only the one of the upward faced
!-- !!!! horizontal surfaces should be taken (the lowest one = ground?). mrt_nlevels number of air
!-- !!!! grid boxes might not be available in case of overhanging structures!
    nmrtbl = 0
    IF ( mrt_nlevels > 0 )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
!
!--             Skip roof if requested
                IF ( mrt_skip_roof  .AND.  surf_usm%isroof_surf(m) )  CYCLE
!
!--             Skip vertical and downward-facing surfaces - offset index in z /= -1.
!--             Note, later on this should be re-thought.
                IF ( .NOT. surf_usm%upward(m) )  CYCLE
!
!--             Cycle over specified no of levels
                nmrtbl = nmrtbl + mrt_nlevels
             ENDDO
!
!--          Ditto for LSM
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
!
!--             Skip vertical and downward-facing surfaces - offset index in z /= -1.
!--             Note, later on this should be re-thought.
                IF ( .NOT. surf_lsm%upward(m) )  CYCLE
                nmrtbl = nmrtbl + mrt_nlevels
             ENDDO
          ENDDO
       ENDDO

       ALLOCATE( mrtbl(iz:ix,nmrtbl), mrtsky(nmrtbl), mrtskyt(nmrtbl), mrtinsw(nmrtbl),            &
                 mrtinlw(nmrtbl), mrt(nmrtbl) )

       imrt = 0
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
!
!--             Skip roof if requested
                IF ( mrt_skip_roof  .AND.  surf_usm%isroof_surf(m) )  CYCLE
!
!--             Skip vertical and downward-facing surfaces - offset index in z /= -1.
!--             Note, later on this should be re-thought.
                IF ( .NOT. surf_usm%upward(m) )  CYCLE
!
!--             Cycle over specified no of levels
                l = surf_usm%k(m)
                DO  k = l + mrt_minlevel, l + mrt_minlevel + mrt_nlevels - 1
                   imrt = imrt + 1
                   mrtbl(:,imrt) = (/k,j,i/)
                ENDDO
             ENDDO
!
!--          Dtto for LSM
             DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
!
!--             Skip vertical and downward-facing surfaces - offset index in z /= -1.
!--             Note, later on this should be re-thought.
                IF ( .NOT. surf_lsm%upward(m) )  CYCLE

                l = surf_lsm%k(m)
                DO  k = l + mrt_minlevel, l + mrt_minlevel + mrt_nlevels - 1
                   imrt = imrt + 1
                   mrtbl(:,imrt) = (/k,j,i/)
                ENDDO
             ENDDO
          ENDDO
! MS mods merged end
       ENDDO
    ENDIF

!
!-- Broadband albedo of the land, roof and wall surface for domain border and sky set artifically
!-- to 1.0 what allows us to calculate heat flux leaving over side and top borders of the domain
    ALLOCATE( albedo_surf(nsurfl) )
    albedo_surf = 1.0_wp
!
!-- Also allocate further array for emissivity with identical order of surface elements as radiation
!-- arrays.
    ALLOCATE( emiss_surf(nsurfl)  )


!
!-- Global array surf of indices of surfaces and displacement index array surfstart
    ALLOCATE( nsurfs(0:numprocs-1) )
!
!-- Required to avoid compiler warnings about array temporaries created by Intel compiler when
!-- calling MPI_ALLGARHERV.
    ALLOCATE( nsurfs_nidx_surf(0:numprocs-1) )
    ALLOCATE( surfstart_nidx_surf(0:numprocs-1) )

#if defined( __parallel )
    CALL MPI_ALLGATHER( nsurfl, 1, MPI_INTEGER, nsurfs, 1, MPI_INTEGER, comm2d, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
      WRITE( debug_string, * ) 'Error MPI_ALLGATHER1:', ierr, nsurfl, nsurfs
      CALL debug_message( debug_string, 'info' )
  ENDIF

#else
    nsurfs(0) = nsurfl
#endif
    ALLOCATE( surfstart(0:numprocs) )
    k = 0
    DO  i = 0, numprocs-1
       surfstart(i) = k
       k = k + nsurfs(i)
    ENDDO
    surfstart(numprocs) = k
    nsurf = k
!
!-- We allocate the array as linear and then use a two-dimensional pointer into it, because some MPI
!-- implementations crash with 2D-allocated arrays.
    ALLOCATE( surf_linear(nidx_surf*nsurf) )
    surf(1:nidx_surf,1:nsurf) => surf_linear(1:nidx_surf*nsurf)

#if defined( __parallel )
    nsurfs_nidx_surf    = nsurfs * nidx_surf
    surfstart_nidx_surf = surfstart(0:numprocs-1) * nidx_surf
    CALL MPI_ALLGATHERV( surfl_linear, nsurfl * nidx_surf, MPI_INTEGER, surf_linear,               &
                         nsurfs_nidx_surf, surfstart_nidx_surf, MPI_INTEGER, comm2d, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
        WRITE( debug_string, * ) 'Error MPI_ALLGATHERV4:', ierr, SIZE( surfl_linear ),             &
                                 nsurfl * nidx_surf, SIZE( surf_linear ), nsurfs * nidx_surf,      &
                                 surfstart(0:numprocs-1) * nidx_surf
        CALL debug_message( debug_string, 'info' )
    ENDIF
#else
    surf = surfl
#endif
!
!-- Allocate and gather global column start indices surfg_col_start
    ALLOCATE( surfg_col_start(0:(nx+1)*(ny+1)) )
#if defined( __parallel )
    CALL MPI_ALLGATHERV( surfl_col_start, nnx*nny, MPI_INTEGER,              &
                         surfg_col_start, nnxy, nnxyd, MPI_INTEGER, comm2d, ierr)
    IF ( ierr /= 0  .AND.  debug_output )  THEN
       WRITE( debug_string, * ) 'Error MPI_ALLGATHER1b:', ierr, surfl_col_start, surfg_col_start
       CALL debug_message( debug_string, 'info' )
    ENDIF
!
!-- Convert local indices (->surfl) to global (->surf)
    DO  i = 0, numprocs-1
       surfg_col_start(nnxyd(i):nnxyd(i)+nnxy(i)-1) =                      &
                surfg_col_start(nnxyd(i):nnxyd(i)+nnxy(i)-1) + surfstart(i)
    ENDDO
#else
    surfg_col_start(0:(nx+1)*(ny+1)-1) = surfl_col_start(0:(nx+1)*(ny+1)-1)
#endif
    surfg_col_start((nx+1)*(ny+1)) = nsurf+1
!
!-- Allocation of the arrays for direct and diffusion radiation.
!-- rad_sw_in, rad_lw_in are computed in the radiation models. Splitting of direct
!-- and diffusion part is done in calc_diffusion_radiation, except for the RRTMG,
!-- where the direct and diffuse portions are computed directly.
    IF ( debug_output )  CALL debug_message( 'allocation of radiation arrays', 'info' )

    IF ( .NOT. ALLOCATED( rad_sw_in_dir ) )  THEN
       ALLOCATE( rad_sw_in_dir(nysg:nyng,nxlg:nxrg) )
       rad_sw_in_dir  = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( rad_sw_in_diff ) )  THEN
       ALLOCATE( rad_sw_in_diff(nysg:nyng,nxlg:nxrg) )
       rad_sw_in_diff = 0.0_wp
    ENDIF
    IF ( .NOT. ALLOCATED( rad_lw_in_diff ) )  THEN
       ALLOCATE( rad_lw_in_diff(nysg:nyng,nxlg:nxrg) )
       rad_lw_in_diff = 0.0_wp
    ENDIF
!
!-- Allocate radiation arrays
    ALLOCATE( surfins(nsurfl) )
    ALLOCATE( surfinl(nsurfl) )
    ALLOCATE( surfinsw(nsurfl) )
    ALLOCATE( surfinlw(nsurfl) )
    ALLOCATE( surfinswdir(nsurfl) )
    ALLOCATE( surfinswdif(nsurfl) )
    ALLOCATE( surfinlwdif(nsurfl) )
    ALLOCATE( surfoutsl(nsurfl) )
    ALLOCATE( surfoutll(nsurfl) )
    ALLOCATE( surfoutsw(nsurfl) )
    ALLOCATE( surfoutlw(nsurfl) )
    ALLOCATE( skyvf(nsurfl) )
    ALLOCATE( skyvft(nsurfl) )
    ALLOCATE( surfemitlwl(nsurfl) )

    IF ( radiation_volumetric_flux )  THEN
       ALLOCATE( skyvf_vol(nz_urban_b:nz_urban_t,nys:nyn,nxl:nxr) )
       ALLOCATE( swflux_vol(nz_urban_b:nz_urban_t,nys:nyn,nxl:nxr) )
       ALLOCATE( rad_shade_h(nys:nyn,nxl:nxr) )
       rad_shade_h = 0
    ENDIF

!
!-- In case of average_radiation, aggregated surface albedo and emissivity, also set initial value
!-- for t_rad_eff.
!-- For now set an arbitrary initial value. Do not overwrite values in restart runs.
    IF ( average_radiation  .AND.  initializing_actions /= 'read_restart_data' )  THEN
       albedo_eff = 0.1_wp
       emissivity_eff = 0.9_wp
       t_rad_eff = pt_surface
    ENDIF

 END SUBROUTINE radiation_interaction_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates shape view factors (SVF), plant sink canopy factors (PCSF), sky-view factors,
!> discretized path for direct solar radiation, MRT factors and other preprocessed data needed for
!> radiation_interaction inside RTM. This subroutine is called only once at the beginning of the
!> simulation. The resulting factors can be stored to files and reused with other simulations
!> utilizing the same surface and plant canopy structure.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_calc_svf

    IMPLICIT NONE

    INTEGER(iwp) ::  i                 !<
    INTEGER(iwp) ::  iaz               !< azimuth counter
    INTEGER(iwp) ::  icsf              !<
    INTEGER(iwp) ::  idzn              !< zenith counter
    INTEGER(iwp) ::  imrt              !<
    INTEGER(iwp) ::  imrtf             !<
    INTEGER(iwp) ::  ip                !<
    INTEGER(iwp) ::  ipcgb             !<
    INTEGER(iwp) ::  isd               !< solar direction
    INTEGER(iwp) ::  isurf             !<
    INTEGER(iwp) ::  isurflt           !<
    INTEGER(iwp) ::  isvf              !<
    INTEGER(iwp) ::  itarg0            !<
    INTEGER(iwp) ::  itarg1            !<
    INTEGER(iwp) ::  j                 !<
    INTEGER(iwp) ::  jp                !<
    INTEGER(iwp) ::  k                 !<
    INTEGER(iwp) ::  kcsf              !<
    INTEGER(iwp) ::  k_topo            !< terrain height
    INTEGER(iwp) ::  max_track_len     !< maximum 2d track length
    INTEGER(iwp) ::  msg_buf_size      !< localized raytracing buffer
    INTEGER(iwp) ::  naz               !< azimuth num of steps
    INTEGER(iwp) ::  npcsfl            !<
    INTEGER(iwp) ::  nzn               !< zenith num of steps
    INTEGER(idp) ::  ray_skip_maxdist  !< skipped raytracing counts
    INTEGER(idp) ::  ray_skip_minval   !< skipped raytracing counts
    INTEGER(iwp) ::  td                !<
    INTEGER(iwp) ::  udim              !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  azlist   !< list of discretized azimuth indices for raytracing
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  icsflt   !<
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  znlist   !< list of discretized zenith indices for raytracing

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  kcsflt_l   !<
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE,TARGET ::  kpcsflt_l  !<

    INTEGER(iwp), DIMENSION(:,:), POINTER ::  kcsflt,kpcsflt  !<

    REAL(wp) ::  azs     !< azimuth cycle step
    REAL(wp) ::  az1     !< relative azimuth of section borders
    REAL(wp) ::  az2     !< relative azimuth of section borders
    REAL(wp) ::  azmid   !< ray (center) azimuth
    REAL(wp) ::  coszen  !< cos(zenith) = sin(elevation)
    REAL(wp) ::  yxlen   !< |yxdir|
    REAL(wp) ::  zns     !< zenith cycle step


    REAL(wp), DIMENSION(2) ::  yxdir   !< y,x *unit* vector of ray direction (in grid units)
    REAL(wp), DIMENSION(3) ::  ta      !< real coordinates z,y,x of source and target

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  horizon        !< horizon heights per vertical level
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  solar_horizon  !< horizon heights per solar direction
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  vffrac0        !< view factor fractions for individual rays (original values)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  zbdry          !< zenith angle boundaries
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  zcent          !< zenith angle centers
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  zdirs          !< directions in z (tangent of elevation)

    REAL(wp),DIMENSION(:), ALLOCATABLE,TARGET ::  csflt_l, pcsflt_l  !<

    REAL(wp),DIMENSION(:,:), POINTER ::  csflt, pcsflt  !<

#if defined( __parallel )
    CHARACTER, DIMENSION(:), POINTER ::  msg_buf_freed !< does nothing but required by the MPI call

    INTEGER(iwp) ::  act_svf      !< auxiliary variables in aggregating
    INTEGER(iwp) ::  act_mrt      !< auxiliary variables in aggregating
    INTEGER(iwp) ::  act_csf      !< auxiliary variables in aggregating
    INTEGER(iwp) ::  index_id     !< auxiliary variables in aggregating
    INTEGER(iwp) ::  iter         !< auxiliary variables in aggregating
    INTEGER(iwp) ::  minfo        !< MPI RMA window info handle
    INTEGER(iwp) ::  niters_surf  !< auxiliary variables in aggregating
    INTEGER(iwp) ::  poz          !< auxiliary variables in aggregating
    INTEGER(iwp) ::  prev_glob    !< auxiliary variables in aggregating
    INTEGER(iwp) ::  proc_id      !< auxiliary variables in aggregating
    INTEGER(iwp) ::  q            !< auxiliary variables in aggregating
    INTEGER(iwp) ::  val          !< auxiliary variables in aggregating

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  size_lad_rma  !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE, TARGET ::  nzterrtl_l  !<
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE, TARGET ::  nzterrbl_l  !<

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  np_revc_radx   !< number of surfaces received from other procs in
                                                               !< MPI_ALLTOALLV rad exchange
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  np_send_radx   !< number of surfaces sent to other procs in
                                                               !< MPI_ALLTOALLV rad exchange
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  recv_buf_temp  !< temporary array used in allocation of recv_buf

    INTEGER(iwp), DIMENSION(:), POINTER, SAVE ::  gridsurf_rma  !< fortran pointer, but lower bounds are 1

    INTEGER(iwp), DIMENSION(:,:), POINTER ::  nzterrbl  !<
    INTEGER(iwp), DIMENSION(:,:), POINTER ::  nzterrtl  !<

    REAL(wp), DIMENSION(:), POINTER, SAVE ::  lad_s_rma  !< fortran 1D pointer

    TYPE(c_ptr) ::  lad_s_rma_p     !< allocated c pointer
    TYPE(c_ptr) ::  gridsurf_rma_p  !< allocated c pointer
#endif

!
!-- Calculation of the SVF
    CALL location_message( 'calculating view factors for radiation interaction', 'start' )
!
!-- Initialize variables and temporary arrays for calculation of svf and csf
    nsvfl  = 0
    ncsfl  = 0
    nsvfla = gasize
    msvf   = 1
    ALLOCATE( asvf1(nsvfla) )
    asvf => asvf1
    IF ( plant_canopy )  THEN
        ncsfla = gasize
        mcsf   = 1
        ALLOCATE( acsf1(ncsfla) )
        acsf => acsf1
    ENDIF
    nmrtf = 0
    IF ( mrt_nlevels > 0 )  THEN
       nmrtfa = gasize
       mmrtf = 1
       ALLOCATE( amrtf1(nmrtfa) )
       amrtf => amrtf1
    ENDIF
    ray_skip_maxdist = 0
    ray_skip_minval = 0

    IF ( localized_raytracing )  THEN
       ALLOCATE( rt2_dist(nz_plant_t-nz_urban_b+2) )
    ELSE
!
!--    Initialize temporary terrain and plant canopy height arrays (global 2D array!)
       ALLOCATE( nzterrt(0:(nx+1)*(ny+1)-1) )
       ALLOCATE( nzterrb(0:(nx+1)*(ny+1)-1) )
#if defined( __parallel )
       ALLOCATE( nzterrtl_l((nyn-nys+1)*(nxr-nxl+1)) )
       ALLOCATE( nzterrbl_l((nyn-nys+1)*(nxr-nxl+1)) )
       nzterrtl(nys:nyn,nxl:nxr) => nzterrtl_l(1:(nyn-nys+1)*(nxr-nxl+1))
       nzterrbl(nys:nyn,nxl:nxr) => nzterrbl_l(1:(nyn-nys+1)*(nxr-nxl+1))
       nzterrtl = topo_top_ind(nys:nyn,nxl:nxr,5)
       nzterrbl = topo_top_ind(nys:nyn,nxl:nxr,0)
       CALL MPI_ALLGATHERV( nzterrtl_l, nnx*nny, MPI_INTEGER, nzterrt,  nnxy, nnxyd, MPI_INTEGER,  &
                            comm2d, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_ALLGATHER1t:', ierr, SIZE( nzterrtl_l ), nnx*nny,   &
                                    SIZE( nzterrt ), nnx*nny
           CALL debug_message( debug_string, 'info' )
       ENDIF
       CALL MPI_ALLGATHERV( nzterrbl_l, nnx*nny, MPI_INTEGER, nzterrb, nnxy, nnxyd, MPI_INTEGER,   &
                            comm2d, ierr)
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_ALLGATHER1b:', ierr, SIZE( nzterrbl_l ), nnx*nny,   &
                                    SIZE( nzterrb ), nnx*nny
           CALL debug_message( debug_string, 'info' )
       ENDIF
       DEALLOCATE( nzterrtl_l )
       DEALLOCATE( nzterrbl_l )
#else
       nzterrt = RESHAPE( topo_top_ind(nys:nyn,nxl:nxr,5), (/(nx+1)*(ny+1)/) )
       nzterrb = RESHAPE( topo_top_ind(nys:nyn,nxl:nxr,0), (/(nx+1)*(ny+1)/) )
#endif
       IF ( plant_canopy )  THEN
           ALLOCATE( plantt(0:(nx+1)*(ny+1)-1) )
           maxboxesg = nx + ny + nz_plant + 1
           max_track_len = nx + ny + 1
!
!--        Temporary arrays storing values for csf calculation during raytracing
           ALLOCATE( boxes(3,maxboxesg) )
           ALLOCATE( crlens(maxboxesg) )

#if defined( __parallel )
           CALL MPI_ALLGATHERV( pct, nnx*nny, MPI_INTEGER, plantt,  nnxy, nnxyd, MPI_INTEGER,      &
                                comm2d, ierr )
           IF ( ierr /= 0  .AND.  debug_output )  THEN
               WRITE( debug_string, * ) 'Error MPI_ALLGATHER2:', ierr, SIZE( pct ), nnx*nny,       &
                                        SIZE( plantt ), nnx*nny
               CALL debug_message( debug_string, 'info' )
           ENDIF
!
!--        Temporary arrays storing values for csf calculation during raytracing
           ALLOCATE( lad_ip(maxboxesg) )
           ALLOCATE( lad_disp(maxboxesg) )

           IF ( raytrace_mpi_rma )  THEN
               ALLOCATE( lad_s_ray(maxboxesg) )
!
!--            Set conditions for RMA communication
               CALL MPI_INFO_CREATE( minfo, ierr )
               IF ( ierr /= 0  .AND.  debug_output )  THEN
                   WRITE( debug_string, * ) 'Error MPI_INFO_CREATE2:', ierr
                   CALL debug_message( debug_string, 'info' )
               ENDIF
               CALL MPI_INFO_SET( minfo, 'accumulate_ordering', 'none', ierr )
               IF ( ierr /= 0  .AND.  debug_output )  THEN
                   WRITE( debug_string, * ) 'Error MPI_INFO_SET5:', ierr
                   CALL debug_message( debug_string, 'info' )
               ENDIF
               CALL MPI_INFO_SET( minfo, 'accumulate_ops', 'same_op', ierr )
               IF ( ierr /= 0  .AND.  debug_output )  THEN
                  WRITE( debug_string, * ) 'Error MPI_INFO_SET6:', ierr
                  CALL debug_message( debug_string, 'info' )
               ENDIF
               IF ( .NOT. non_uniform_subdomain )  THEN
                  CALL MPI_INFO_SET(minfo, 'same_size', 'true', ierr)
                  IF ( ierr /= 0  .AND.  debug_output )  THEN
                     WRITE( debug_string , * ) 'Error MPI_INFO_SET7:', ierr
                     CALL debug_message( debug_string, 'info' )
                  ENDIF
               ENDIF
               CALL MPI_INFO_SET( minfo, 'same_disp_unit', 'true', ierr )
               IF ( ierr /= 0  .AND.  debug_output )  THEN
                   WRITE( debug_string, * ) 'Error MPI_INFO_SET8:', ierr
                   CALL debug_message( debug_string, 'info' )
               ENDIF

!--            Allocate and initialize the MPI RMA window, must be in accordance with allocation of
!--            lad_s in plant_canopy_model, optimization of memory should be done.
!--            Argument X of function STORAGE_SIZE(X) needs arbitrary REAL(wp) value, set to 1.0_wp
!--            for now.
               size_lad_rma = STORAGE_SIZE( 1.0_wp ) / 8 * nnx * nny * nz_plant
               CALL MPI_WIN_ALLOCATE( size_lad_rma, STORAGE_SIZE( 1.0_wp ) / 8, minfo, comm2d,     &
                                      lad_s_rma_p, win_lad, ierr )
               IF ( ierr /= 0  .AND.  debug_output )  THEN
                   WRITE( debug_string, * ) 'Error MPI_WIN_ALLOCATE2:', ierr, size_lad_rma,        &
                                            STORAGE_SIZE( 1.0_wp ) / 8, win_lad
                   CALL debug_message( debug_string, 'info' )
               ENDIF
               CALL C_F_POINTER( lad_s_rma_p, lad_s_rma, (/ nz_plant*nny*nnx /) )
               sub_lad(nz_urban_b:nz_plant_t, nys:nyn, nxl:nxr) => lad_s_rma(1:nz_plant*nny*nnx)
           ELSE
               ALLOCATE( sub_lad(nz_urban_b:nz_plant_t, nys:nyn, nxl:nxr) )
           ENDIF
#else
           plantt = RESHAPE( pct(nys:nyn,nxl:nxr), (/(nx+1)*(ny+1)/) )
           ALLOCATE( sub_lad(nz_urban_b:nz_plant_t, nys:nyn, nxl:nxr) )
#endif
           sub_lad(:,:,:) = 0.0_wp
           DO  i = nxl, nxr
              DO  j = nys, nyn
                 k = topo_top_ind(j,i,0)
                 sub_lad(k:nz_plant_t, j, i) = lad_s(0:nz_plant_t-k, j, i)
              ENDDO
           ENDDO

           plantt_max = MAXVAL( plantt )
           ALLOCATE( rt2_track(2, max_track_len),                                                  &
                     rt2_track_lad(nz_urban_b:plantt_max, max_track_len),                          &
                     rt2_track_dist(0:max_track_len), rt2_dist(plantt_max-nz_urban_b+2) )

#if defined( __parallel )
           IF ( raytrace_mpi_rma )  THEN
              CALL MPI_INFO_FREE( minfo, ierr )
              IF ( ierr /= 0  .AND.  debug_output )  THEN
                 WRITE( debug_string, * ) 'Error MPI_INFO_FREE2:', ierr
                 CALL debug_message( debug_string, 'info' )
              ENDIF
              CALL MPI_WIN_LOCK_ALL( 0, win_lad, ierr )
              IF ( ierr /= 0  .AND.  debug_output )  THEN
                 WRITE( debug_string, * ) 'Error MPI_WIN_LOCK_ALL1:', ierr, win_lad
                 CALL debug_message( debug_string, 'info' )
              ENDIF
           ELSE
              ALLOCATE( sub_lad_g(0:(nx+1)*(ny+1)*nz_plant-1) )
              CALL MPI_ALLGATHERV( sub_lad, nnx*nny*nz_plant, MPI_REAL, sub_lad_g, nnxy*nz_plant,  &
                                   nnxyd*nz_plant, MPI_REAL, comm2d, ierr )
              IF ( ierr /= 0  .AND.  debug_output )  THEN
                 WRITE( debug_string, * ) 'Error MPI_ALLGATHER3:', ierr, SIZE( sub_lad ),          &
                                          nnx*nny*nz_plant, SIZE( sub_lad_g ), nnx*nny*nz_plant
                 CALL debug_message( debug_string, 'info' )
              ENDIF
           ENDIF
#endif
       ENDIF ! plant_canopy
    ENDIF ! .NOT. localized_raytracing

    IF ( localized_raytracing )  THEN
       CALL lrt_allocate_message( lrt_msg_incoming_section, lrt_msg_incoming_rays,                 &
                                  lrt_msg_type_incoming )
       CALL lrt_allocate_message( lrt_msg_processing_section, lrt_msg_processing_rays,             &
                                  lrt_msg_type_processing, msg_buf_size )

#if defined( __parallel )
       lrt_wait_time = 0
!
!--    Allocate buffer space for outgoing messages. We send at most one message per azimuth, and at
!--    the end we send one message to each process.
!--    NOTE: If MPI_ALLOC_MEM would be more appropriate, but it does not work with Intel compilers
!--    (see description of lrt_allocate_message).
       msg_buf_size = msg_buf_size * 2 * MAX( numprocs, raytrace_discrete_azims )
       ALLOCATE( lrt_msg_buf(msg_buf_size) )
       CALL MPI_BUFFER_ATTACH( lrt_msg_buf, msg_buf_size, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
          WRITE( debug_string, * ) 'MPI_BUFFER_ATTACH error:', ierr
          CALL debug_message( debug_string, 'info' )
       ENDIF

       lrt_prev_process_complete = ( myid == 0 )
       lrt_local_complete = .FALSE.
       lrt_unfinished = .TRUE.

       CALL MPI_IRECV( MPI_BOTTOM, 1, lrt_msg_type_incoming, MPI_ANY_SOURCE, lrt_msg_tag, comm2d,  &
                       lrt_req_incoming, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
          WRITE( debug_string, * ) 'MPI Irecv Error:', ierr
          CALL debug_message( debug_string, 'info' )
       ENDIF
#endif

       ALLOCATE( gridsurf(0:nsurf_type_u-1,nz_urban_b:nz_urban_t,nys:nyn,nxl:nxr) )
    ELSE
!
!--    Prepare the MPI_WIN for collecting the surface indices from the reverse index arrays grids
!--    from processors of target surfaces
!
!--    Allocate and fill the reverse indexing array gridsurf
#if defined( __parallel )
!
!--    raytrace_mpi_rma is asserted
       CALL MPI_INFO_CREATE( minfo, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_INFO_CREATE1:', ierr
           CALL debug_message( debug_string, 'info' )
       ENDIF
       CALL MPI_INFO_SET( minfo, 'accumulate_ordering', 'none', ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_INFO_SET1:', ierr
           CALL debug_message( debug_string, 'info' )
       ENDIF
       CALL MPI_INFO_SET( minfo, 'accumulate_ops', 'same_op', ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_INFO_SET2:', ierr
           CALL debug_message( debug_string, 'info' )
       ENDIF
       IF ( .NOT. non_uniform_subdomain )  THEN
           CALL MPI_INFO_SET(minfo, 'same_size', 'true', ierr)
           IF ( ierr /= 0  .AND.  debug_output )  THEN
              WRITE( debug_string, * ) 'Error MPI_INFO_SET3:', ierr
              CALL debug_message( debug_string, 'info' )
           ENDIF
       ENDIF

       CALL MPI_INFO_SET( minfo, 'same_disp_unit', 'true', ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_INFO_SET4:', ierr
           CALL debug_message( debug_string, 'info' )
       ENDIF

       CALL MPI_WIN_ALLOCATE( INT( STORAGE_SIZE( 1_iwp ) / 8 * nsurf_type_u * nz_urban * nny * nnx,&
                              KIND = MPI_ADDRESS_KIND ), STORAGE_SIZE( 1_iwp ) / 8, minfo, comm2d, &
                              gridsurf_rma_p, win_gridsurf, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_WIN_ALLOCATE1:', ierr,                              &
                         INT( STORAGE_SIZE( 1_iwp ) / 8 * nsurf_type_u * nz_urban * nny * nnx,     &
                         KIND = MPI_ADDRESS_KIND ), STORAGE_SIZE( 1_iwp ) / 8, win_gridsurf
           CALL debug_message( debug_string, 'info' )
       ENDIF

       CALL MPI_INFO_FREE( minfo, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_INFO_FREE1:', ierr
           CALL debug_message( debug_string, 'info' )
       ENDIF

!
!--    On Intel compilers, calling C_F_POINTER to transform a C pointer directly to a
!--    multi-dimensional Fotran pointer leads to strange errors on dimension boundaries. However,
!--    transforming to a 1D pointer and then redirecting a multidimensional pointer to it works fine.
       CALL C_F_POINTER( gridsurf_rma_p, gridsurf_rma, (/ nsurf_type_u * nz_urban * nny * nnx /) )
       gridsurf(0:nsurf_type_u-1, nz_urban_b:nz_urban_t, nys:nyn, nxl:nxr) =>                      &
                                                gridsurf_rma(1:nsurf_type_u * nz_urban * nny * nnx)
#else
       ALLOCATE( gridsurf(0:nsurf_type_u-1,nz_urban_b:nz_urban_t,nys:nyn,nxl:nxr) )
#endif
    ENDIF ! localized_raytracing
!
!-- Populate gridsurf with reverse global indices (->surf)
    gridsurf(:,:,:,:) = -999
    DO  isurf = 1, nsurfl
       gridsurf(surfl(id,isurf),surfl(iz,isurf), surfl(iy,isurf),surfl(ix,isurf)) = isurf +        &
                                                                                    surfstart(myid)
    ENDDO

#if defined( __parallel )
    IF ( .NOT. localized_raytracing )  THEN
!
!--     Prepare the MPI_WIN for collecting the surface indices from the reverse index arrays gridsurf
!--     from processors of target surfaces
!--     raytrace_mpi_rma is asserted
        CALL MPI_WIN_LOCK_ALL( 0, win_gridsurf, ierr )
        IF ( ierr /= 0  .AND.  debug_output )  THEN
            WRITE( debug_string, * ) 'Error MPI_WIN_LOCK_ALL2:', ierr, win_gridsurf
            CALL debug_message( debug_string, 'info' )
        ENDIF
    ENDIF
#endif
!
!-- Directions opposite to face normals are not even calculated, they must be preset to 0
    dsitrans(:,:) = 0.0_wp
!
!-- Prepare indices for discretized angles
    ALLOCATE( discr_azim_cent(1:raytrace_discrete_azims),                                          &
              discr_azim_bdry(0:raytrace_discrete_azims),                                          &
              discr_elev_cent(1:raytrace_discrete_elevs),                                          &
              discr_elev_bdry(0:raytrace_discrete_elevs),                                          &
              azlist(1:raytrace_discrete_azims),                                                   &
              znlist(1:raytrace_discrete_elevs),                                                   &
              discr_azim_yxdir(2,1:raytrace_discrete_azims) )

    azs = 2.0_wp * pi / REAL( raytrace_discrete_azims, wp )
    discr_azim_cent(:) = (/( ( REAL( i, wp ) - 0.5_wp ) * azs, i = 1, raytrace_discrete_azims )/)
    discr_azim_bdry(:) = (/(   REAL( i, wp )            * azs, i = 0, raytrace_discrete_azims )/)

    zns = pi / REAL( raytrace_discrete_elevs, wp )
    discr_elev_cent(:) = (/( ( REAL( i, wp ) - 0.5_wp ) * zns, i = 1, raytrace_discrete_elevs )/)
    discr_elev_bdry(:) = (/(   REAL( i, wp )            * zns, i = 0, raytrace_discrete_elevs )/)

    DO  i = 1, raytrace_discrete_azims
       discr_azim_yxdir(:,i) = (/ COS( discr_azim_cent(i) ) / dy, SIN( discr_azim_cent(i) ) / dx /)
       yxlen = SQRT( SUM( discr_azim_yxdir(:,i)**2 ) )
       discr_azim_yxdir(:,i) = discr_azim_yxdir(:,i) / yxlen
    ENDDO
!
!-- Main SVF loop (local surface targets)
    CALL cpu_log( log_point_s(97), 'rtm_raytrace_surfaces', 'start' )
    DO  isurflt = 1, nsurfl
!
!--    Determine face centers
       td = surfl(id, isurflt)
       ta = (/ REAL( surfl(iz,isurflt), wp ) - 0.5_wp * kdir(td),                                  &
               REAL( surfl(iy,isurflt), wp ) - 0.5_wp * jdir(td),                                  &
               REAL( surfl(ix,isurflt), wp ) - 0.5_wp * idir(td)  /)
!
!--    Calculate sky view factor and raytrace DSI paths
       skyvf(isurflt)  = 0.0_wp
       skyvft(isurflt) = 0.0_wp
!
!--    Select a proper half-sphere for 2D raytracing
       SELECT CASE ( td )
          CASE ( iup )
             naz = raytrace_discrete_azims
             azlist(1:naz) = (/( i, i = 1, raytrace_discrete_azims )/)
             nzn = raytrace_discrete_elevs / 2
             znlist(1:nzn) = (/( i, i = 1, raytrace_discrete_elevs / 2 )/)
          CASE ( idown )
             naz = raytrace_discrete_azims
             azlist(1:naz) = (/( i, i = 1, raytrace_discrete_azims )/)
             nzn = raytrace_discrete_elevs / 2
             znlist(1:nzn) = (/( i, i = raytrace_discrete_elevs / 2 + 1, raytrace_discrete_elevs )/)
          CASE ( isouth )
             naz = raytrace_discrete_azims / 2
             azlist(1:naz) = (/( i, i = raytrace_discrete_azims / 4 + 1,                           &
                                        raytrace_discrete_azims * 3 / 4 )/)
             nzn = raytrace_discrete_elevs
             znlist(1:nzn) = (/( i, i = 1, raytrace_discrete_elevs )/)
          CASE ( inorth )
!
!--          For north-facing faces the open-side azimuths go from 270 through 360/0 to 90,
!--          therefore it is the only direction that needs the MODULO function.
             naz = raytrace_discrete_azims / 2
             azlist(1:naz) = (/( MODULO( i, raytrace_discrete_azims ) + 1,                         &
                                 i = -raytrace_discrete_azims / 4,                                 &
                                     raytrace_discrete_azims / 4 - 1 )/)
             nzn = raytrace_discrete_elevs
             znlist(1:nzn) = (/( i, i = 1, raytrace_discrete_elevs )/)
          CASE ( iwest )
             naz = raytrace_discrete_azims / 2
             azlist(1:naz) = (/( i, i = raytrace_discrete_azims / 2 + 1, raytrace_discrete_azims )/)
             nzn = raytrace_discrete_elevs
             znlist(1:nzn) = (/( i, i = 1, raytrace_discrete_elevs )/)
          CASE ( ieast )
             naz = raytrace_discrete_azims / 2
             azlist(1:naz) = (/( i, i = 1, raytrace_discrete_azims / 2 )/)
             nzn = raytrace_discrete_elevs
             znlist(1:nzn) = (/( i, i = 1, raytrace_discrete_elevs )/)
          CASE DEFAULT
             WRITE( message_string, * ) 'the surface type ', td, 'is not supported for ' // &
                                        'calculating SVF'
             CALL message( 'radiation_calc_svf', 'RAD0052', 1, 2, 0, 6, 0 )
       END SELECT

       ALLOCATE( zdirs(1:nzn), zcent(1:nzn), zbdry(0:nzn), vffrac(1:nzn*naz), ztransp(1:nzn*naz),  &
                 itarget(1:nzn*naz), vffrac0(1:nzn) ) !TODO improve allocation

       itarg0 = 1
       itarg1 = nzn
       nsvf_ins = 0
       zcent(:) = (/( discr_elev_cent(znlist(i)), i = 1, nzn )/)
       zbdry(0) = discr_elev_bdry(znlist(1) - 1)
       zbdry(1:nzn) = (/( discr_elev_bdry(znlist(i)), i = 1, nzn )/)
!
!--    For horizontal target, vf fractions are constant per azimuth
       IF ( td == iup )  THEN
          vffrac0(1:nzn) = ( COS( 2 * zbdry(0:nzn-1) ) - COS( 2 * zbdry(1:nzn) ) ) / 2.0_wp /      &
                            REAL( naz, wp )
       ELSEIF ( td == idown )  THEN
          vffrac0(1:nzn) = - ( COS( 2 * zbdry(0:nzn-1) ) - COS( 2 * zbdry(1:nzn) ) ) / 2.0_wp /    &
                              REAL( naz, wp )
       ENDIF
!
!--    Calculate sky-view factor and direct solar visibility using 2D raytracing
       DO  iaz = 1, naz
          azmid = discr_azim_cent(azlist(iaz))
          IF ( td /= iup  .AND.  td /= idown )  THEN
             az2 = REAL( iaz, wp ) * azs - pi / 2.0_wp
             az1 = az2 - azs
             vffrac0(1:nzn) = ( SIN( az2 ) - SIN( az1 ) ) *                                        &
                                ( zbdry(1:nzn) - zbdry(0:nzn-1) +                                  &
                                  SIN( zbdry(0:nzn-1) ) * COS( zbdry(0:nzn-1) ) -                  &
                                  SIN( zbdry(1:nzn) ) * COS( zbdry(1:nzn) ) ) / ( 2.0_wp * pi )

          ENDIF
          yxdir(:) = (/ COS( azmid ) / dy, SIN( azmid ) / dx /)
          yxlen = SQRT( SUM( yxdir(:)**2 ) )
          zdirs(:) = COS( zcent(:) ) / ( dz(1) * yxlen * SIN( zcent(:) ) )
          yxdir(:) = yxdir(:) / yxlen

          IF ( localized_raytracing )  THEN
             CALL raytrace_init( 'f', ta, azlist(iaz), nzn, zdirs, surfstart(myid) + isurflt,      &
                                 facearea(td), vffrac0, lrt_msg_processing_section,                &
                                 lrt_msg_processing_rays, lrt_msg_type_processing )
#if defined( __parallel )
             CALL lrt_process_pending( .FALSE. )
#endif
          ELSE
             vffrac(itarg0:itarg1) = vffrac0(1:nzn)

             CALL raytrace_2d( ta, yxdir, nzn, zdirs, surfstart(myid) + isurflt, facearea(td),     &
                               vffrac(itarg0:itarg1), .TRUE., .TRUE., .FALSE.,                     &
                               ztransp(itarg0:itarg1), itarget(itarg0:itarg1) )

             skyvf(isurflt) = skyvf(isurflt) + SUM( vffrac(itarg0:itarg1),                         &
                                                    MASK = ( itarget(itarg0:itarg1) < 0 ) )
             skyvft(isurflt) = skyvft(isurflt) + SUM( ztransp(itarg0:itarg1) *                     &
                                                      vffrac(itarg0:itarg1),                       &
                                                      MASK = ( itarget(itarg0:itarg1) < 0 ) )
!
!--          Save direct solar transparency.  For down direction there is no direct irradiance.
             IF ( td /= idown )  THEN
                DO  k = 1, raytrace_discrete_elevs / 2
                   i = dsidir_rev(k-1, azlist(iaz)-1)
                   IF ( i /= -1  .AND.  itarget(itarg0+k-1) < 0 )                                  &
                      dsitrans(isurflt,i) = ztransp(itarg0+k-1)
                ENDDO
             ENDIF
!
!--          Advance itarget indices
             itarg0 = itarg1 + 1
             itarg1 = itarg1 + nzn
          ENDIF
       ENDDO

#if defined( __parallel )
       IF ( localized_raytracing )  THEN
!
!--       Process remaining messages until all rays have been traced
          DO WHILE ( nsvf_ins < nzn*naz )
             CALL lrt_process_pending( .TRUE. )
          ENDDO
       ENDIF
#endif
!
!--    Sort itarget by face id
       CALL quicksort_itarget( itarget, vffrac, ztransp, 1, nzn*naz )
!
!--    For aggregation, we need fractions multiplied by transmissivities
       ztransp(:) = vffrac(:) * ztransp(:)
!
!--    Find the first valid position
       itarg0 = 1
       DO WHILE ( itarg0 <= nzn*naz )
          IF ( itarget(itarg0) >= 0 )  EXIT
          itarg0 = itarg0 + 1
       ENDDO

       DO  i = itarg0, nzn*naz
!
!--       For duplicate values, only sum up vf fraction value
          IF ( i < nzn*naz )  THEN
             IF ( itarget(i+1) == itarget(i) )  THEN
                vffrac(i+1) = vffrac(i+1) + vffrac(i)
                ztransp(i+1) = ztransp(i+1) + ztransp(i)
                CYCLE
             ENDIF
          ENDIF
!
!--       Write to the svf array
          nsvfl = nsvfl + 1
!
!--       Check dimmension of asvf array and enlarge it if needed
          IF ( nsvfla < nsvfl )  THEN
             k = CEILING( REAL( nsvfla, KIND = wp ) * grow_factor )
             IF ( msvf == 0 )  THEN
                msvf = 1
                ALLOCATE( asvf1(k) )
                asvf => asvf1
                asvf1(1:nsvfla) = asvf2
                DEALLOCATE( asvf2 )
             ELSE
                msvf = 0
                ALLOCATE( asvf2(k) )
                asvf => asvf2
                asvf2(1:nsvfla) = asvf1
                DEALLOCATE( asvf1 )
             ENDIF

             IF ( debug_output )  THEN
                WRITE( debug_string, '(A,3I12)' ) 'Grow asvf:', nsvfl, nsvfla, k
                CALL debug_message( debug_string, 'info' )
             ENDIF

             nsvfla = k
          ENDIF
!
!--       Write svf values into the array
          asvf(nsvfl)%isurflt = isurflt
          asvf(nsvfl)%isurfs = itarget(i)
          asvf(nsvfl)%rsvf = vffrac(i)
          asvf(nsvfl)%rtransp = ztransp(i) / vffrac(i)
       ENDDO

       DEALLOCATE( zdirs, zcent, zbdry, vffrac, ztransp, itarget, vffrac0 )

#if defined( __parallel )
       IF ( localized_raytracing )  CALL lrt_process_pending( .FALSE. )
#endif
    ENDDO
    CALL cpu_log( log_point_s(97), 'rtm_raytrace_surfaces', 'stop' )

!
!-- Raytrace to canopy boxes to fill dsitransc
!-- TODO: consider replacing by DSI rays toward surfaces
    dsitransc(:,:) = 0.0_wp
    naz = raytrace_discrete_azims
    nzn = raytrace_discrete_elevs / 2
    ALLOCATE( zdirs(1:nzn), zcent(1:nzn), vffrac(1:nzn), ztransp(1:nzn), itarget(1:nzn) )
    zcent(:) = discr_elev_cent(1:nzn)
    vffrac(:) = 0.0_wp

    DO  ipcgb = 1, npcbl
       ta = (/ REAL( pcbl(iz,ipcgb), wp ),                                                         &
               REAL( pcbl(iy,ipcgb), wp ),                                                         &
               REAL( pcbl(ix,ipcgb), wp ) /)
       nsvf_ins = 0
!
!--    Calculate direct solar visibility using 2D raytracing
       DO  iaz = 1, naz
          azmid = discr_azim_cent(iaz)
          yxdir(:) = (/ COS( azmid ) / dy, SIN( azmid ) / dx /)
          yxlen = SQRT( SUM( yxdir(:)**2 ) )
          zdirs(:) = COS( zcent(:) ) / ( dz(1) * yxlen * SIN( zcent(:) ) )
          yxdir(:) = yxdir(:) / yxlen
          IF ( localized_raytracing )  THEN
             CALL raytrace_init( 'p', ta, iaz, nzn, zdirs, ipcgb, -999.0_wp, vffrac,               &
                                 lrt_msg_processing_section, lrt_msg_processing_rays,              &
                                 lrt_msg_type_processing )
#if defined( __parallel )
             CALL lrt_process_pending( .FALSE. )
#endif
          ELSE
             CALL raytrace_2d( ta, yxdir, nzn, zdirs, -999, -999.0_wp, vffrac, .FALSE., .FALSE.,   &
                               .TRUE., ztransp, itarget )
!
!--          Save direct solar transparency
             DO  k = 1, raytrace_discrete_elevs / 2
                i = dsidir_rev(k-1, iaz-1)
                IF ( i /= -1  .AND.  itarget(k) < 0 )  dsitransc(ipcgb, i) = ztransp(k)
             ENDDO
          ENDIF
       ENDDO

#if defined( __parallel )
       IF ( localized_raytracing )  THEN
!
!--       Process remaining messages until all rays have been traced
          DO WHILE ( nsvf_ins < naz )
             CALL lrt_process_pending( .TRUE. )
          ENDDO
       ENDIF
#endif
    ENDDO
    DEALLOCATE( zdirs, zcent, vffrac, ztransp, itarget )
!
!-- Raytrace to MRT boxes
    CALL cpu_log( log_point_s(98), 'rtm_raytrace_mrt', 'start' )
    IF ( nmrtbl > 0 )  THEN
       mrtdsit(:,:) = 0.0_wp
       mrtsky(:) = 0.0_wp
       mrtskyt(:) = 0.0_wp
       naz = raytrace_discrete_azims
       nzn = raytrace_discrete_elevs
       ALLOCATE( zdirs(1:nzn), zcent(1:nzn), zbdry(0:nzn), vffrac(1:nzn*naz), vffrac0(1:nzn),      &
                 ztransp(1:nzn*naz), itarget(1:nzn*naz) )

       zcent(:) = discr_elev_cent(1:nzn)
       zbdry(:) = discr_elev_bdry(0:nzn)
       vffrac0(:) = ( COS( zbdry(0:nzn-1) ) - COS( zbdry(1:nzn) ) ) / 2.0_wp / REAL( naz, wp )
!
!--    Modify direction weights to simulate human body (lower weight for irradiance from zenith,
!--    higher from sides) depending on selection.
!--    For mrt_geom=0, no weighting is done (simulates spherical globe thermometer).
       SELECT CASE ( mrt_geom )

       CASE ( 1 )
          vffrac0(:) = vffrac0(:) * MAX( 0.0_wp, SIN( zcent(:) ) * mrt_geom_params(2)              &
                                               + COS( zcent(:) ) * mrt_geom_params(1) )
          vffrac0(:) = vffrac0(:) / ( SUM( vffrac0 ) * REAL( naz, wp ) )

       CASE ( 2 )
          vffrac0(:) = vffrac0(:) * SQRT( ( mrt_geom_params(1) * COS( zcent(:) ) )** 2 +           &
                                          ( mrt_geom_params(2) * SIN( zcent(:) ) )** 2 )
          vffrac0(:) = vffrac0(:) / ( SUM( vffrac0 ) * REAL( naz, wp ) )

       END SELECT

       DO  imrt = 1, nmrtbl
          ta = (/ REAL( mrtbl(iz,imrt), wp ),                                                      &
                  REAL( mrtbl(iy,imrt), wp ),                                                      &
                  REAL( mrtbl(ix,imrt), wp ) /)
          nsvf_ins = 0
!
!--       vf fractions are constant per azimuth
          DO  iaz = 0, naz-1
             vffrac(iaz*nzn+1:(iaz+1)*nzn) = vffrac0(:)
          ENDDO
!
!--       Sum of whole vffrac equals 1, verified
          itarg0 = 1
          itarg1 = nzn
!
!--       Calculate sky-view factor and direct solar visibility using 2D raytracing
          DO  iaz = 1, naz
             azmid = discr_azim_cent(iaz)
             yxdir(:) = (/ COS( azmid ) / dy, SIN( azmid ) / dx /)
             yxlen = SQRT( SUM( yxdir(:)**2 ) )
             zdirs(:) = COS( zcent(:) ) / ( dz(1) * yxlen * SIN( zcent(:) ) )
             yxdir(:) = yxdir(:) / yxlen

             IF ( localized_raytracing )  THEN
                CALL raytrace_init( 'm', ta, iaz, nzn, zdirs, imrt, -999.0_wp, vffrac0,            &
                                  lrt_msg_processing_section, lrt_msg_processing_rays,             &
                                  lrt_msg_type_processing )
#if defined( __parallel )
                CALL lrt_process_pending( .FALSE. )
#endif
             ELSE
                CALL raytrace_2d( ta, yxdir, nzn, zdirs, -999, -999.0_wp, vffrac(itarg0:itarg1),   &
                                  .TRUE., .FALSE., .FALSE., ztransp(itarg0:itarg1),                &
                                  itarget(itarg0:itarg1) )
!
!--             Sky view factors for MRT
                mrtsky(imrt) = mrtsky(imrt) + SUM( vffrac(itarg0:itarg1),                          &
                                                   MASK = ( itarget(itarg0:itarg1) < 0 ) )
                mrtskyt(imrt) = mrtskyt(imrt) + SUM( ztransp(itarg0:itarg1) *                      &
                                                     vffrac(itarg0:itarg1),                        &
                                                     MASK = ( itarget(itarg0:itarg1) < 0 ) )
!
!--             Direct solar transparency for MRT
                DO  k = 1, raytrace_discrete_elevs/2
                   i = dsidir_rev(k-1, iaz-1)
                   IF ( i /= -1  .AND.  itarget(itarg0+k-1) < 0 )  THEN
                      mrtdsit(imrt,i) = ztransp(itarg0+k-1)
                   ENDIF
                ENDDO
             ENDIF
!
!--          Advance itarget indices
             itarg0 = itarg1 + 1
             itarg1 = itarg1 + nzn
          ENDDO

#if defined( __parallel )
          IF ( localized_raytracing )  THEN
!
!--          Process remaining messages until all rays have been traced
             DO WHILE ( nsvf_ins < nzn*naz )
                CALL lrt_process_pending( .TRUE. )
             ENDDO
          ENDIF
#endif
!
!--       Sort itarget by face id
          CALL quicksort_itarget( itarget, vffrac, ztransp, 1, nzn * naz )
!
!--       For aggregation, we need fractions multiplied by transmissivities
          ztransp(:) = vffrac(:) * ztransp(:)
!
!--       Find the first valid position
          itarg0 = 1
          DO WHILE ( itarg0 <= nzn * naz )
             IF ( itarget(itarg0) >= 0 )  EXIT
             itarg0 = itarg0 + 1
          ENDDO

          DO  i = itarg0, nzn*naz
!
!--          For duplicate values, only sum up vf fraction value
             IF ( i < nzn * naz )  THEN
                IF ( itarget(i+1) == itarget(i) )  THEN
                   vffrac(i+1) = vffrac(i+1) + vffrac(i)
                   ztransp(i+1) = ztransp(i+1) + ztransp(i)
                   CYCLE
                ENDIF
             ENDIF
!
!--          Some MRT geometries might contain directions with zero weight, they need to be
!--          skipped here. This check may be removed if such geometries are also removed.
             IF ( vffrac(i) <= 0.0_wp )  CYCLE
!
!--          Write to the mrtf array
             nmrtf = nmrtf + 1
!
!--          Check dimmension of mrtf array and enlarge it if needed
             IF ( nmrtfa < nmrtf )  THEN
                k = CEILING( REAL( nmrtfa, KIND = wp ) * grow_factor )
                IF ( mmrtf == 0 )  THEN
                   mmrtf = 1
                   ALLOCATE( amrtf1(k) )
                   amrtf => amrtf1
                   amrtf1(1:nmrtfa) = amrtf2
                   DEALLOCATE( amrtf2 )
                ELSE
                   mmrtf = 0
                   ALLOCATE( amrtf2(k) )
                   amrtf => amrtf2
                   amrtf2(1:nmrtfa) = amrtf1
                   DEALLOCATE( amrtf1 )
                ENDIF

                IF ( debug_output )  THEN
                   WRITE( debug_string, '(A,3I12)' ) 'Grow amrtf:', nmrtf, nmrtfa, k
                   CALL debug_message( debug_string, 'info' )
                ENDIF

                nmrtfa = k
             ENDIF
!
!--          Write mrtf values into the array
             amrtf(nmrtf)%isurflt = imrt
             amrtf(nmrtf)%isurfs = itarget(i)
             amrtf(nmrtf)%rsvf = vffrac(i)
             amrtf(nmrtf)%rtransp = ztransp(i) / vffrac(i)
          ENDDO  ! itarg

#if defined( __parallel )
          IF ( localized_raytracing )  CALL lrt_process_pending( .FALSE. )
#endif
       ENDDO  ! imrt

       DEALLOCATE( zdirs, zcent, zbdry, vffrac, vffrac0, ztransp, itarget )
!
!--    Move MRT factors to final arrays
       ALLOCATE( mrtf(nmrtf), mrtft(nmrtf), mrtfsurf(2,nmrtf) )
       DO  imrtf = 1, nmrtf
          mrtf(imrtf) = amrtf(imrtf)%rsvf
          mrtft(imrtf) = amrtf(imrtf)%rsvf * amrtf(imrtf)%rtransp
          mrtfsurf(:,imrtf) = (/amrtf(imrtf)%isurflt, amrtf(imrtf)%isurfs /)
       ENDDO
       IF ( ALLOCATED( amrtf1 ) )  DEALLOCATE( amrtf1 )
       IF ( ALLOCATED( amrtf2 ) )  DEALLOCATE( amrtf2 )
    ENDIF ! nmrtbl > 0
    CALL cpu_log( log_point_s(98), 'rtm_raytrace_mrt', 'stop' )

    IF ( debug_output )  CALL debug_message( 'waiting for completion of SVF and CSF ' //           &
                                             'calculation in all processes', 'info' )

    IF ( localized_raytracing )  THEN
#if defined( __parallel )
!
!--    Mark local work as complete
       lrt_local_complete = .TRUE.
       CALL lrt_check_completion( lrt_msg_processing_section, lrt_msg_type_processing )
!
!--    Wait for termination message
       DO WHILE ( lrt_unfinished )
          CALL lrt_process_pending( .TRUE. )
       ENDDO
#endif
!
!--    Free message buffers for localized raytracing
       CALL lrt_free_message( lrt_msg_incoming_section, lrt_msg_incoming_rays,                     &
                              lrt_msg_type_incoming )
       CALL lrt_free_message( lrt_msg_processing_section, lrt_msg_processing_rays,                 &
                              lrt_msg_type_processing )

#if defined( __parallel )
       msg_buf_freed => lrt_msg_buf ! does nothing but required by interface
       CALL MPI_BUFFER_DETACH( msg_buf_freed, msg_buf_size, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_BUFFER_DETACH:', ierr
           CALL debug_message( debug_string, 'info' )
       ENDIF
       DEALLOCATE( lrt_msg_buf )

       IF ( debug_output )  THEN
          CALL SYSTEM_CLOCK( COUNT_RATE=i )
          WRITE( debug_string, '("Localized raytracing spent ",I0," / ",I0," seconds waiting.")' ) &
             lrt_wait_time, i
          CALL debug_message( debug_string, 'info' )
       ENDIF
#endif
    ELSE
#if defined( __parallel )
!
!--    Finalize MPI_RMA communication established to get global index of the surface from grid
!--    indices.
!--    Flush all MPI window pending requests.
       CALL MPI_WIN_FLUSH_ALL( win_gridsurf, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_WIN_FLUSH_ALL1:', ierr, win_gridsurf
           CALL debug_message( debug_string, 'info' )
       ENDIF
!
!--    Unlock MPI window
       CALL MPI_WIN_UNLOCK_ALL( win_gridsurf, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_WIN_UNLOCK_ALL1:', ierr, win_gridsurf
           CALL debug_message( debug_string, 'info' )
       ENDIF

!
!--    Free MPI window
       CALL MPI_WIN_FREE( win_gridsurf, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_WIN_FREE1:', ierr, win_gridsurf
           CALL debug_message( debug_string, 'info' )
       ENDIF
#else
       DEALLOCATE( gridsurf )
#endif
!
!--    Deallocate temporary global arrays
       DEALLOCATE( nzterrt, nzterrb )

       IF ( plant_canopy )  THEN
!
!--       Finalize mpi_rma communication and deallocate temporary arrays
#if defined( __parallel )
          IF ( raytrace_mpi_rma )  THEN
             CALL MPI_WIN_FLUSH_ALL( win_lad, ierr )
             IF ( ierr /= 0  .AND.  debug_output )  THEN
                WRITE( debug_string, * ) 'Error MPI_WIN_FLUSH_ALL2:', ierr, win_lad
                CALL debug_message( debug_string, 'info' )
             ENDIF
!
!--          Unlock MPI window
             CALL MPI_WIN_UNLOCK_ALL( win_lad, ierr )
             IF ( ierr /= 0  .AND.  debug_output )  THEN
                WRITE( debug_string, * ) 'Error MPI_WIN_UNLOCK_ALL2:', ierr, win_lad
                CALL debug_message( debug_string, 'info' )
             ENDIF
!
!--          Free MPI window
             CALL MPI_WIN_FREE( win_lad, ierr )
             IF ( ierr /= 0  .AND.  debug_output )  THEN
                WRITE( debug_string, * ) 'Error MPI_WIN_FREE2:', ierr, win_lad
                CALL debug_message( debug_string, 'info' )
             ENDIF
!
!--          Deallocate temporary arrays storing values for csf calculation during raytracing.
!--          sub_lad is the pointer to lad_s_rma in case of raytrace_mpi_rma and must not be
!--          deallocated here.
             DEALLOCATE( lad_s_ray )
          ELSE
             DEALLOCATE( sub_lad )
             DEALLOCATE( sub_lad_g )
          ENDIF
#else
          DEALLOCATE( sub_lad )
#endif
          DEALLOCATE( boxes )
          DEALLOCATE( crlens )
          DEALLOCATE( plantt )
          DEALLOCATE( rt2_track, rt2_track_lad, rt2_track_dist, rt2_dist )
       ENDIF ! plant_canopy
    ENDIF ! localized_raytracing
!
!-- Perform horizon tracing for volumetric fluxes
    IF ( radiation_volumetric_flux )  THEN
       IF ( debug_output )  THEN
          CALL debug_message( 'Calculating factors for volumetric fluxes', 'info' )
       ENDIF
!
!--    Initialize temporary opaque top array (global 2D array!)
       ALLOCATE( opaque_top(0:(nx+1)*(ny+1)-1) )
       ALLOCATE( opaque_top_l_lin((nyn-nys+1)*(nxr-nxl+1)) )
       opaque_top_l(nys:nyn,nxl:nxr) => opaque_top_l_lin(1:(nyn-nys+1)*(nxr-nxl+1))
!
!--    Determine opaque top from terrain and plant canopy
       opaque_top_l = topo_top_ind(nys:nyn,nxl:nxr,5)
       IF ( plant_canopy )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                k_topo = topo_top_ind(j,i,0)
                DO  k = pch(j,i), opaque_top_l(j,i)+1, -1
                   IF ( lad_s(k,j,i) >= min_opaque_lad )  THEN
!
!--                   We are at the top of opaque plant canopy
                      opaque_top_l(j,i) = k + k_topo
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF

#if defined( __parallel )
       CALL MPI_ALLGATHERV( opaque_top_l_lin, nnx*nny, MPI_INTEGER, opaque_top,  nnxy, nnxyd,      &
                            MPI_INTEGER, comm2d, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
           WRITE( debug_string, * ) 'Error MPI_ALLGATHER opaq:', ierr, SIZE( opaque_top_l_lin ),   &
                                    nnx*nny, SIZE( opaque_top ), nnx*nny
           CALL debug_message( debug_string, 'info' )
       ENDIF
#else
       opaque_top(:) = opaque_top_l_lin(:)
#endif

       ALLOCATE( shadow_top(nys:nyn,nxl:nxr,ndsidir) )
       ALLOCATE( horizon(nz_urban_b:nz_urban_t) )
       ALLOCATE( solar_horizon(ndsidir) )
!
!--    Calculate tangent of elevation angle in physical coords for each discretized solar direction
       DO  isd = 1, ndsidir
          solar_horizon(isd) = dsidir(1,isd) / sqrt( dsidir(2,isd)**2 + dsidir(3,isd)**2 )
       ENDDO

       naz = raytrace_discrete_azims ! TODO: separate setting, separate list of presim solar pos

       skyvf_vol(:,:,:) = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Cycle all discretized azimuths
             DO  iaz = 1, naz
                azmid = discr_azim_cent(iaz)
                yxdir(:) = (/ COS( azmid ), SIN( azmid ) /) ! for trace_horizons, we use physical
                                                            ! coords as opposed to grid coords
                CALL trace_horizons( i, j, yxdir, horizon(:))

                DO  idzn = 1, raytrace_discrete_elevs / 2
                   isd = dsidir_rev(idzn-1, iaz-1)
                   IF ( isd /= -1 )  THEN
                      shadow_top(j,i,isd) = opaque_top_l(j,i)
                      DO  k = opaque_top_l(j,i)+1, nz_urban_t
                         IF ( horizon(k) < solar_horizon(isd) )  EXIT
                         shadow_top(j,i,isd) = k
                      ENDDO
                   ENDIF
                ENDDO
!
!--             Calculate sky view factor
                DO  k = opaque_top_l(j,i)+1, nz_urban_t
                   coszen = SIN( ATAN( horizon(k) ) ) ! cos(zen) = sin(elev)
                       ! this could be possibly rewritten as sqrt( horz**2/(1+horz**2) ) * sgn(horz)
                   skyvf_vol(k,j,i) = skyvf_vol(k,j,i) + ( 1 - coszen )
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       skyvf_vol(:,:,:) = skyvf_vol(:,:,:) / 2.0_wp / REAL( naz, wp )

       DEALLOCATE( opaque_top, horizon, solar_horizon )
    ENDIF

    IF ( debug_output )  CALL debug_message( 'calculation of the complete SVF array', 'info' )

    IF ( debug_output )  THEN
       WRITE( debug_string, '("Load ",I0," SVFs from the structure array to plain arrays")' )   &
              nsvfl
       CALL debug_message( debug_string, 'info' )
    ENDIF
    ALLOCATE( svf(ndsvf,nsvfl) )
    ALLOCATE( svfsurf(idsvf,nsvfl) )

    DO  isvf = 1, nsvfl
       svf(:,isvf) = (/ asvf(isvf)%rsvf, asvf(isvf)%rtransp /)
       svfsurf(:,isvf) = (/ asvf(isvf)%isurflt, asvf(isvf)%isurfs /)
    ENDDO
!
!-- Deallocate temporary asvf array
!-- DEALLOCATE(asvf) - ifort has a problem with deallocation of allocatable target via pointing
!-- pointer - we need to test original targets
    IF ( ALLOCATED( asvf1 ) )  DEALLOCATE( asvf1 )
    IF ( ALLOCATED( asvf2 ) )  DEALLOCATE( asvf2 )

    npcsfl = 0
    IF ( plant_canopy )  THEN
        IF ( debug_output )  CALL debug_message( 'Calculation of the complete CSF array', 'info' )
!
!--     Sort and merge csf for the last time, keeping the array size to minimum
        CALL merge_and_grow_csf( - 1 )
!
!--     Aggregate csb among processors.
!--     Allocate necessary arrays.
        udim = MAX( ncsfl, 1 )
        ALLOCATE( csflt_l(ndcsf*udim) )
        csflt(1:ndcsf,1:udim) => csflt_l(1:ndcsf*udim)
        ALLOCATE( kcsflt_l(kdcsf*udim) )
        kcsflt(1:kdcsf,1:udim) => kcsflt_l(1:kdcsf*udim)
        ALLOCATE( icsflt(0:numprocs-1) )

!--     Fill out arrays of csf values and arrays of number of elements and displacements for
!--     particular precessors.
        icsflt = 0
        ip = -1
        j = -1
        DO  kcsf = 1, ncsfl
           j = j + 1
           IF ( acsf(kcsf)%ip /= ip )  THEN
!
!--           New block of the processor number of elements of previous block
              IF ( ip >= 0 )  icsflt(ip) = j
!
!--           Blank blocks
              DO  jp = ip+1, acsf(kcsf)%ip-1
!
!--              Number of elements is zero, displacement is equal to previous
                 icsflt(jp) = 0
              ENDDO
!
!--           The actual block
              ip = acsf(kcsf)%ip
              j = 0
           ENDIF
           csflt(1,kcsf) = acsf(kcsf)%rcvf
!
!--        Fill out integer values of itz,ity,itx,isurfs
           kcsflt(1,kcsf) = acsf(kcsf)%itz
           kcsflt(2,kcsf) = acsf(kcsf)%ity
           kcsflt(3,kcsf) = acsf(kcsf)%itx
           kcsflt(4,kcsf) = acsf(kcsf)%isurfs
        ENDDO
!
!--     Last blank blocks at the end of array
        j = j+1
        IF ( ip >= 0 )  icsflt(ip) = j
        DO  jp = ip+1, numprocs-1
!
!--        Number of elements is zero, displacement is equal to previous
           icsflt(jp) = 0
        ENDDO
!
!--     Deallocate temporary acsf array
!--     DEALLOCATE(acsf) - ifort has a problem with deallocation of allocatable target via pointing
!--     pointer - we need to test original targets
        IF ( ALLOCATED( acsf1 ) )  DEALLOCATE( acsf1 )
        IF ( ALLOCATED( acsf2 ) )  DEALLOCATE( acsf2 )

!
!--     We exchange csf fields between processors only if __parallel AND NOT localized_raytracing,
!--     otherwise we just copy.
#if defined( __parallel )
        IF ( localized_raytracing )  THEN
#else
        IF ( .TRUE. )  THEN
#endif
           npcsfl = ncsfl
           ALLOCATE( pcsflt(  ndcsf, MAX( npcsfl, 1 ) ),                                           &
                     kpcsflt( kdcsf, MAX( npcsfl, 1 ) ) )
           pcsflt = csflt
           kpcsflt = kcsflt
!
!--        Just silence compiler warning about unused variable
           IF ( bufsize_alltoall <= 0 )  bufsize_alltoall = 0
        ELSE
           IF ( debug_output )  THEN
              CALL debug_message( 'Exchange CSF fields between processors', 'start' )
           ENDIF

           CALL radiation_exchange_alltoall( icsflt, kdcsf, ndcsf, kcsflt_l, csflt_l, npcsfl,      &
                                             kpcsflt_l, pcsflt_l )
           pcsflt(1:ndcsf,1:npcsfl) => pcsflt_l(0:ndcsf*npcsfl-1)
           kpcsflt(1:kdcsf,1:npcsfl) => kpcsflt_l(0:kdcsf*npcsfl-1)

           IF ( debug_output )  THEN
              CALL debug_message( 'Exchange CSF fields between processors', 'end' )
           ENDIF
        ENDIF
!
!--     Deallocate temporary arrays
        DEALLOCATE( csflt_l )
        DEALLOCATE( kcsflt_l )
        DEALLOCATE( icsflt )
!
!--     Sort csf ( a version of quicksort )
        IF ( debug_output )  CALL debug_message( 'Sort csf', 'info' )
        CALL quicksort_csf2( kpcsflt, pcsflt, 1, npcsfl )
!
!--     Aggregate canopy sink factor records with identical box & source againg across all values
!--     from all processors
        IF ( debug_output )  CALL debug_message( 'Aggregate canopy sink factor records with ' //   &
                                                 'identical box', 'info' )

        IF ( npcsfl > 0 )  THEN
            icsf = 1 !< reading index
            kcsf = 1 !< writing index
            DO WHILE ( icsf < npcsfl )
!
!--             Here kpcsf(kcsf) already has values from kpcsf(icsf)
                IF ( kpcsflt(3,icsf) == kpcsflt(3,icsf+1)  .AND.                                   &
                     kpcsflt(2,icsf) == kpcsflt(2,icsf+1)  .AND.                                   &
                     kpcsflt(1,icsf) == kpcsflt(1,icsf+1)  .AND.                                   &
                     kpcsflt(4,icsf) == kpcsflt(4,icsf+1) )  THEN

                    pcsflt(1,kcsf) = pcsflt(1,kcsf) + pcsflt(1,icsf+1)
!
!--                 Advance reading index, keep writing index
                    icsf = icsf + 1
                ELSE
!
!--                 Not identical, just advance and copy
                    icsf = icsf + 1
                    kcsf = kcsf + 1
                    kpcsflt(:,kcsf) = kpcsflt(:,icsf)
                    pcsflt(:,kcsf) = pcsflt(:,icsf)
                ENDIF
            ENDDO
!
!--         Last written item is now also the last item in valid part of array
            npcsfl = kcsf
        ENDIF

        ncsfl = npcsfl
        IF ( ncsfl > 0 )  THEN
            ALLOCATE( csf(ndcsf,ncsfl) )
            ALLOCATE( csfsurf(idcsf,ncsfl) )
            DO  icsf = 1, ncsfl
                csf(:,icsf) = pcsflt(:,icsf)
                csfsurf(1,icsf) =  gridpcbl(kpcsflt(1,icsf),kpcsflt(2,icsf),kpcsflt(3,icsf))
                csfsurf(2,icsf) =  kpcsflt(4,icsf)
            ENDDO
        ENDIF
!
!--     Deallocation of temporary arrays
        IF ( npcbl > 0 )  DEALLOCATE( gridpcbl )

#if defined( __parallel )
        IF ( .NOT. localized_raytracing )  THEN
           DEALLOCATE( pcsflt_l )
           DEALLOCATE( kpcsflt_l )
        ENDIF
#endif

        IF ( debug_output )  THEN
           WRITE( debug_string, '("Finished aggregating ",I0," CSFs.")' ) ncsfl
           CALL debug_message( debug_string, 'info' )
        ENDIF

    ENDIF

#if defined( __parallel )
!
!-- MPI surface exchange optimization
!-- sort svf, mrt, pcb - sort by ___surf(:,2)
    CALL quicksort_target_svf( svfsurf, svf, 1, nsvfl )
    IF ( nmrtf > 0 )  THEN
       CALL quicksort_target_mrt( mrtfsurf, mrtf, mrtft, 1, nmrtf )
    ENDIF
    IF ( plant_canopy  .AND.  ncsfl > 0 )  THEN
       CALL quicksort_target_csf( csfsurf, csf, 1, ncsfl )
    ENDIF
    ALLOCATE( recv_buf_temp(nsvfl+nmrtf+ncsfl) )
!
!-- Aggregation of target surface
    prev_glob = -1;
    i = 1
    j = 1
    k = 1
    index_id = 0
!
!-- Can not be parallelized by OMP
    DO  q = 1, ( nsvfl + nmrtf + ncsfl )
!
!--    Check if i,j or k overflow its size
       IF ( i > nsvfl )  THEN
          act_svf = nsurf + 10 !< this should be the highiest value
       ELSE
          act_svf = svfsurf(2,i)
       ENDIF
       IF ( nmrtf > 0 )  THEN
          IF ( j > nmrtf )  THEN
              act_mrt = nsurf + 10
          ELSE
              act_mrt = mrtfsurf(2,j)
          ENDIF
       ELSE
          act_mrt = nsurf + 10
       ENDIF
       IF ( plant_canopy  .AND.  ncsfl > 0 )  THEN
          IF ( k > ncsfl )  THEN
             act_csf = nsurf + 10
          ELSE
             act_csf = csfsurf(2,k)
          ENDIF
       ELSE
          act_csf = nsurf + 10
       ENDIF
       poz = MINLOC( (/ act_svf, act_mrt, act_csf /), DIM = 1 )
       val = MINVAL( (/ act_svf, act_mrt, act_csf /), DIM = 1 )
       IF ( val /= prev_glob )  THEN
          IF (val /= -1)  THEN
!
!--          New value
             index_id = index_id + 1
             recv_buf_temp(index_id) = val
          ENDIF
          prev_glob = val
       ENDIF
       IF ( poz == 1 )  THEN
!
!--       Lowest value has svf, compare it with prev value
          svfsurf(2,i) = index_id
          i = i + 1
       ELSEIF ( poz == 2 )  THEN
!
!--       Mrt case
          IF ( nmrtf > 0 )  THEN
              mrtfsurf(2,j) = index_id
              j = j + 1
          ENDIF
       ELSEIF ( poz == 3 )  THEN
!
!--       Csf case, there can be -1 value, pointing to sky
          IF ( plant_canopy  .AND.  ncsfl > 0 )  THEN
             IF (val == -1)  THEN
                csfsurf(2,k) = -1
             ELSE
                csfsurf(2,k) = index_id
             ENDIF
             k = k + 1
          ENDIF
       ENDIF
    ENDDO
!
!-- Check if loop was done correctly TODO: is it always true, can it be removed now?
    IF ( ( q /= ( i+j+k-2 ) )  .OR.  ( i-1 /= nsvfl )  .OR.  ( j-1 /= nmrtf )  .OR.                &
         ( k-1 /= ncsfl ) )                                                                        &
    THEN
       IF ( debug_output )  THEN
          WRITE( 9, * ) 'Error in RTM MPI_ALLTOALL indexing', q, i, j, k, nsvfl, nmrtf, ncsfl
       ENDIF
    ENDIF
    nrecv_radx = index_id
!
!-- Move values from temp to original one
    ALLOCATE( isurf_recv_radx(nrecv_radx) )
    isurf_recv_radx = recv_buf_temp(1:nrecv_radx)
    DEALLOCATE( recv_buf_temp )
!
!-- Allocate arrays with counts and displacements for radiation exchange
    ALLOCATE( np_revc_radx(0:numprocs-1), disp_recv_radx(0:numprocs) )
    ALLOCATE( np_send_radx(0:numprocs-1), disp_send_radx(0:numprocs) )
!
!-- Assessing proc id to each surf in aggregated list
    proc_id = 0
    np_revc_radx = 0
!
!-- Can not be parallelized by OMP
    DO  i = 1, nrecv_radx
        IF ( ( isurf_recv_radx(i) > surfstart(proc_id) )  .AND.                                    &
           ( isurf_recv_radx(i) <= surfstart(proc_id+1) ) )  THEN
!
!--        Surface is between proc_id and proc_id + 1
           np_revc_radx(proc_id) = np_revc_radx(proc_id) + 1
        ELSE
!
!--        Surface is not in interval, find next processor that fits interval
           proc_id = proc_id + 1
           DO  j = proc_id, numprocs
              IF ( ( isurf_recv_radx(i) > surfstart(proc_id) )  .AND.                              &
                   ( isurf_recv_radx(i) <= surfstart(proc_id+1) ) )  THEN
!
!--               Next interval was found
                  np_revc_radx(proc_id) = np_revc_radx(proc_id) + 1
                  EXIT
               ENDIF
               proc_id = proc_id + 1
           ENDDO
        ENDIF
    ENDDO
!
!-- Check if all surfaces were assigned to target processor
    IF ( ( SUM( np_revc_radx ) /= nrecv_radx )  .AND.  debug_output )  THEN
        WRITE( debug_string, * ) 'ERROR IN SUM( np_revc_radx ) /= nrecv_radx    ',                 &
            'proc_id, numprocs, SUM( np_revc_radx ), nrecv_radx',                                  &
             proc_id, numprocs, SUM( np_revc_radx ), nrecv_radx
        CALL debug_message( debug_string, 'info' )
    ENDIF
!
!-- Send how many of target svf #proc requires from others
    CALL MPI_ALLTOALL( np_revc_radx, 1, MPI_INTEGER, np_send_radx, 1, MPI_INTEGER, comm2d, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
       WRITE( debug_string, * ) 'Error at MPI_ALLTOALL 01:', ierr, size( np_revc_radx ),           &
                                size( np_send_radx )
       CALL debug_message( debug_string, 'info' )
    ENDIF
    nsend_radx = SUM( np_send_radx )   ! send to others
!
!-- Calculate displacement to AllTOAllV routine
    disp_recv_radx(0) = 0
    disp_send_radx(0) = 0
!
!-- Can not be parallelized by OMP
    DO  i = 1, numprocs
        disp_recv_radx(i) = disp_recv_radx(i-1) + np_revc_radx(i-1)
        disp_send_radx(i) = disp_send_radx(i-1) + np_send_radx(i-1)
    ENDDO
!
!-- Send the information to other and receive info from other
!-- svf_send_buf store info about local/global index of required surface from local proc
!-- svf_recv_buf store info about local/global index of required surface index from other proc
    ALLOCATE( isurf_send_radx(nsend_radx) )
    ALLOCATE( radx_send(nsend_radx) )
    ALLOCATE( surfoutl_recv(nrecv_radx) )
    ALLOCATE( surfouts_recv(nrecv_radx) )
    ALLOCATE( radx_send_surfinl(nrecv_radx) )
    ALLOCATE( surfinl_recv(nsend_radx) )
!
!-- Determine number of iterations among all processes
!-- (e.g. this process may have nothing to send and receive, yet some other still might)
    IF ( bufsize_alltoall <= 0 )  THEN
       niters_radx = 1
       nmaxsend_radx = HUGE( niters_radx )
    ELSE
       nmaxsend_radx = bufsize_alltoall
       niters_surf = ( MAXVAL( np_send_radx(:) ) + nmaxsend_radx - 1 ) / nmaxsend_radx
       CALL MPI_ALLREDUCE( niters_surf, niters_radx, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
       IF ( niters_radx > 1 )  THEN
          WRITE( debug_string, '("The MPI_ALLTOALL call has been split to ",I8," iterations ' //   &
                               'of max. ",I12," records each.")' ) niters_radx, bufsize_alltoall
          CALL debug_message( debug_string, 'info' )
       ENDIF
    ENDIF
!
!-- Allocation of immediate send and receive buffers in MPI ALLTOALLV radiation exchange
    ALLOCATE( disp_sendbuf_radx(numprocs) )
    ALLOCATE( disp_recvbuf_radx(numprocs) )

    disp_recvbuf_radx(:) = disp_send_radx(0:numprocs-1)
    disp_sendbuf_radx(:) = disp_recv_radx(0:numprocs-1)

    ALLOCATE( num_sendbuf_radx(numprocs) )
    ALLOCATE( num_recvbuf_radx(numprocs) )
!
!-- Iterate ALLTOALLV using max-sized buffers
    DO  iter = 1, niters_radx
       num_sendbuf_radx(:) = MIN( disp_recv_radx(1:) - disp_sendbuf_radx(:), nmaxsend_radx )
       num_recvbuf_radx(:) = MIN( disp_send_radx(1:) - disp_recvbuf_radx(:), nmaxsend_radx )
!
!--    Send integer data
       CALL MPI_ALLTOALLV( isurf_recv_radx, num_sendbuf_radx(:), disp_sendbuf_radx(:), MPI_INTEGER,&
                           isurf_send_radx, num_recvbuf_radx(:), disp_recvbuf_radx(:), MPI_INTEGER,&
                           comm2d, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
          WRITE( debug_string, * ) 'Error at MPI_ALLTOALLV isurf_recv_radx sending:',              &
                                   ierr, iter, nmaxsend_radx, disp_send_radx, disp_sendbuf_radx,   &
                                   num_sendbuf_radx, disp_recv_radx, disp_recvbuf_radx,            &
                                   num_recvbuf_radx
          CALL debug_message( debug_string, 'info' )
       ENDIF
!
!--    Shift displacements for next iteration
       disp_sendbuf_radx(:) = disp_sendbuf_radx(:) + num_sendbuf_radx(:)
       disp_recvbuf_radx(:) = disp_recvbuf_radx(:) + num_recvbuf_radx(:)
    ENDDO
#endif

#if defined( __parallel )
    CALL MPI_BARRIER( comm2d, ierr )
#endif
    CALL location_message( 'calculating view factors for radiation interaction', 'finished' )

 END SUBROUTINE radiation_calc_svf


# if defined( __parallel )
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> rtm_alltoallv is used to send/recv target surface radiations between processor using
!> MPI_ALLTOALLV subroutine. And based on bufsize_alltoall can split single ALLTOALL call
!> into multiple ones. Used only for send/recv float data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rtm_alltoallv( send_float, dsend, recv_float, drecv )

    IMPLICIT NONE

    INTEGER(iwp) ::  iter      !< current iteration

    INTEGER(iwp), DIMENSION(0:), INTENT(IN) ::  drecv   !< received data displacements per proc
    INTEGER(iwp), DIMENSION(0:), INTENT(IN) ::  dsend   !< sent data displacements per process

    REAL(wp), DIMENSION(:), INTENT(OUT) ::  recv_float  !< float receive buffer
    REAL(wp), DIMENSION(:), INTENT(IN)  ::  send_float  !< send buffer with floats

    disp_recvbuf_radx(:) = drecv(0:numprocs-1)
    disp_sendbuf_radx(:) = dsend(0:numprocs-1)

!
!-- Iterate ALLTOALLv using max-sized buffers
    DO  iter = 1, niters_radx
       num_sendbuf_radx(:) = MIN( dsend(1:) - disp_sendbuf_radx(:), nmaxsend_radx )
       num_recvbuf_radx(:) = MIN( drecv(1:) - disp_recvbuf_radx(:), nmaxsend_radx )

!
!--    Send floating point data
       CALL MPI_ALLTOALLV( send_float, num_sendbuf_radx(:), disp_sendbuf_radx(:), MPI_REAL,        &
                           recv_float, num_recvbuf_radx(:), disp_recvbuf_radx(:), MPI_REAL,        &
                           comm2d, ierr )

       IF ( ierr /= 0  .AND.  debug_output )  THEN
          WRITE( debug_string, * ) 'Error at MPI_ALLTOALLV rtm_alltoallv:',                        &
                                   ierr, iter, nmaxsend_radx,                                      &
                                   dsend, disp_sendbuf_radx, num_sendbuf_radx,                     &
                                   drecv, disp_recvbuf_radx, num_recvbuf_radx
          CALL debug_message( debug_string, 'info' )
       ENDIF
!
!--    Shift displacements for next iteration
       disp_sendbuf_radx(:) = disp_sendbuf_radx(:) + num_sendbuf_radx(:)
       disp_recvbuf_radx(:) = disp_recvbuf_radx(:) + num_recvbuf_radx(:)
    ENDDO

 END SUBROUTINE rtm_alltoallv


#endif
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Quicksort algorithm for sorfting svfsurf array according to it's second row.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_target_svf( svfsurf, svf, first, last )

    IMPLICIT NONE

    INTEGER(iwp) ::  i, j  !<

    INTEGER(iwp), INTENT(IN) ::  first, last  !<

    INTEGER(iwp),DIMENSION(:,:), INTENT(INOUT) :: svfsurf !<

    INTEGER(iwp) :: x, t !<

    REAL(wp), DIMENSION(:,:), INTENT(INOUT) ::  svf !<
    REAL(wp) :: tt !<

    IF ( first >= last )  RETURN
    x = svfsurf(2, ( first + last ) / 2)
    i = first
    j = last
    DO
       DO WHILE ( target_lt(svfsurf(2, i),x) )
          i = i + 1
       ENDDO
       DO WHILE ( target_lt(x,svfsurf(2, j)) )
          j = j - 1
       ENDDO
       IF ( i >= j ) EXIT
       t = svfsurf(2, i);  svfsurf(2, i) = svfsurf(2, j);  svfsurf(2, j) = t
       t = svfsurf(1, i);  svfsurf(1, i) = svfsurf(1, j);  svfsurf(1, j) = t

       tt = svf(1, i);  svf(1, i) = svf(1, j);  svf(1, j) = tt
       tt = svf(2, i);  svf(2, i) = svf(2, j);  svf(2, j) = tt
       i = i+1
       j = j-1
    ENDDO
    IF ( first < i-1  )  CALL quicksort_target_svf( svfsurf, svf, first, i - 1 )
    IF ( j+1   < last )  CALL quicksort_target_svf( svfsurf, svf, j + 1, last  )

 END SUBROUTINE quicksort_target_svf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Quicksort algorithm for sorfting mrtfsurf array according to it's second row.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_target_mrt( mrtfsurf, mrtf, mrtft, first, last )

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  t  !<
    INTEGER(iwp) ::  x  !<

    INTEGER(iwp), INTENT(IN) ::  first  !<
    INTEGER(iwp), INTENT(IN) ::  last  !<

    INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) :: mrtfsurf  !<

    REAL(wp) ::  tt  !<

    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  mrtf  !<
    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  mrtft !<


    IF ( first >= last )  RETURN
    x = mrtfsurf(2, ( first + last ) / 2)
    i = first
    j = last
    DO
       DO WHILE ( target_lt(mrtfsurf(2, i),x) )
         i = i+1
       ENDDO
       DO WHILE ( target_lt(x,mrtfsurf(2, j)) )
          j=j-1
       ENDDO
       IF ( i >= j ) EXIT
       t = mrtfsurf(2, i);  mrtfsurf(2, i) = mrtfsurf(2, j);  mrtfsurf(2, j) = t
       t = mrtfsurf(1, i);  mrtfsurf(1, i) = mrtfsurf(1, j);  mrtfsurf(1, j) = t

       tt = mrtf(i);   mrtf(i)  = mrtf(j);   mrtf(j)  = tt
       tt = mrtft(i);  mrtft(i) = mrtft(j);  mrtft(j) = tt
       i = i + 1
       j = j - 1
    ENDDO
    IF ( first < i-1  )  CALL quicksort_target_mrt( mrtfsurf, mrtf, mrtft, first, i - 1 )
    IF ( j+1   < last )  CALL quicksort_target_mrt( mrtfsurf, mrtf, mrtft, j + 1, last  )

 END SUBROUTINE quicksort_target_mrt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Quicksort algorithm for sorfting csfsurf array according to it's second row.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_target_csf( csfsurf, csf, first, last )

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  t  !<
    INTEGER(iwp) ::  x  !<

    INTEGER(iwp), INTENT(IN) ::  first  !<
    INTEGER(iwp), INTENT(IN) ::  last   !<

    INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) ::  csfsurf  !<

    REAL(wp) :: tt  !<

    REAL(wp), DIMENSION(:,:), INTENT(INOUT) ::  csf  !<


    IF ( first >= last )  RETURN
    x = csfsurf(2, ( first + last ) / 2)
    i = first
    j = last
    DO
       DO WHILE ( target_lt(csfsurf(2, i),x) )
         i = i+1
       ENDDO
       DO WHILE ( target_lt(x,csfsurf(2, j)) )
          j=j-1
       ENDDO
       IF ( i >= j ) EXIT
       t = csfsurf(2, i);  csfsurf(2, i) = csfsurf(2, j);  csfsurf(2, j) = t
       t = csfsurf(1, i);  csfsurf(1, i) = csfsurf(1, j);  csfsurf(1, j) = t

       tt = csf(1, i);  csf(1, i) = csf(1, j);  csf(1, j) = tt
       i= i+1
       j= j-1
    ENDDO
    IF ( first < i-1 )  CALL quicksort_target_csf( csfsurf, csf, first, i - 1 )
    IF ( j+1 < last  )  CALL quicksort_target_csf( csfsurf, csf, j + 1, last )

 END SUBROUTINE quicksort_target_csf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Function for quicksort algorithm, which returns whether target1 is bigger than target2
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION target_lt( target1, target2 ) RESULT( res )

    INTEGER, INTENT(IN) :: target1  !<
    INTEGER, INTENT(IN) :: target2  !<

    LOGICAL :: res  !<

    IF ( target1 < target2  .OR. ( target1 == target2  .AND.  target1 < target2) )  THEN
       res = .TRUE.
    ELSE
       res = .FALSE.
    ENDIF

 END FUNCTION target_lt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> A new, more efficient version of ray tracing algorithm that processes a whole arc instead of a
!> single ray (new in RTM version 2.5).
!>
!> In all comments, horizon means tangent of horizon angle, i.e. vertical_delta /
!> horizontal_distance
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE raytrace_2d( origin, yxdir, nrays, zdirs, iorig, aorig, vffrac, calc_svf, create_csf,  &
                         skip_1st_pcb, transparency, itarget )
    IMPLICIT NONE

    INTEGER(iwp)             ::  d                 !< running index for directions
    INTEGER(iwp)             ::  i                 !< running index
    INTEGER(iwp)             ::  iray              !< index into zdirs
    INTEGER(iwp)             ::  isurf             !< index into surf(l)
    INTEGER(iwp)             ::  ip                !< number of processor where gridbox reside
    INTEGER(iwp)             ::  ip_last           !< previous number of processor where gridbox reside
    INTEGER(iwp)             ::  ig                !< 1D index of grid column in global 2D array
    INTEGER(iwp)             ::  ig_last           !< 1D index of previous column in global 2D array
    INTEGER(iwp), INTENT(IN) ::  iorig             !< index of origin face for csf
    INTEGER(iwp)             ::  k                 !< running index
    INTEGER(iwp)             ::  kz                !< running index for z-coordinate
    INTEGER(iwp)             ::  l                 !< running index
    INTEGER(iwp)             ::  lastdir           !< wall direction before hitting this column
    INTEGER(iwp)             ::  lowest_free_ray   !< index into zdirs
    INTEGER(iwp)             ::  lowest_mixed_ray  !< index into zdirs
    INTEGER(iwp)             ::  maxboxes          !< max no of CSF created
    INTEGER(iwp)             ::  nrays             !< number of rays (z directions) to raytrace
    INTEGER(iwp)             ::  nly               !< maximum  plant canopy height
    INTEGER(iwp)             ::  ntrack            !< number of horizontal points in the ray path
    INTEGER(iwp)             ::  nz                !< number of z coordinates for ray sub-path
    INTEGER(iwp)             ::  seldim            !< dimension to be incremented
    INTEGER(iwp)             ::  zb0               !< z boundary
    INTEGER(iwp)             ::  zb1               !< z boundary
    INTEGER(iwp)             ::  zsgn              !< sign of z increments in a column sub-path

    INTEGER(iwp), DIMENSION(2) ::  column      !< grid column being crossed
    INTEGER(iwp), DIMENSION(2) ::  dimnext     !< next dimension increments along path
    INTEGER(iwp), DIMENSION(2) ::  dimdelta    !< dimension direction = +- 1
    INTEGER(iwp), DIMENSION(2) ::  lastcolumn  !< previous xy coords

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  down_col  !< downward oriented surfaces in current column
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  up_col    !< upward oriented surfaces in current column
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  vert_col  !< vertical faces from previous to current column

    INTEGER(iwp), DIMENSION(nrays), INTENT(OUT) ::  itarget     !< global indices of target faces for zdirs or <0 for sky
    INTEGER(iwp), DIMENSION(nrays)              ::  ray_ntrack  !< number of track columns of each ray until its end

    LOGICAL, INTENT(IN) ::  calc_svf      !< whether to calculate SFV (identify obstacle surfaces)
    LOGICAL, INTENT(IN) ::  create_csf    !< whether to create canopy sink factors
    LOGICAL             ::  is_mixed_col  !< whether current column contains full-3D geometry
    LOGICAL, INTENT(IN) ::  skip_1st_pcb  !< whether to skip first plant canopy box during raytracing

    REAL(wp), INTENT(IN) ::  aorig             !< origin face area for csf
    REAL(wp)             ::  curtrans          !< transparency of current PC box crossing
    REAL(wp)             ::  dxxyy             !< square of real horizontal distance
    REAL(wp)             ::  full_horizon      !< highest full horizon found after raytracing (z/hdist)
    REAL(wp)             ::  horz_full_entry   !< full horizon at entry to column
    REAL(wp)             ::  horz_full_exit    !< full horizon at exit from column
    REAL(wp)             ::  horz_mixed_entry  !< mixed horizon at entry to column
    REAL(wp)             ::  horz_mixed_exit   !< mixed horizon at exit from column
    REAL(wp)             ::  lastdist          !< beginning of current crossing
    REAL(wp)             ::  nextdist          !< end of current crossing
    REAL(wp)             ::  qdist             !< ratio of real distance to z coord difference
    REAL(wp)             ::  zbottom, ztop     !< urban surface boundary in real numbers
    REAL(wp)             ::  zexit             !< z coordinate of ray column exit
    REAL(wp)             ::  zorig             !< z coordinate of ray column entry

    REAL(wp), DIMENSION(2)             ::  dimnextdist  !< distance for each dimension increments
    REAL(wp), DIMENSION(2)             ::  glob_max     !< global domain boundary (north y, east x)
    REAL(wp), DIMENSION(3), INTENT(IN) ::  origin       !< z,y,x coordinates of ray origin
    REAL(wp), DIMENSION(2), INTENT(IN) ::  yxdir        !< y,x *unit* vector of ray direction (in grid units)
    REAL(wp), DIMENSION(2)             ::  yxorigin     !< horizontal copy of origin (y,x)

    REAL(wp), DIMENSION(nrays), INTENT(OUT) ::  transparency  !< transparencies of zdirs paths
    REAL(wp), DIMENSION(nrays), INTENT(IN)  ::  vffrac        !< view factor fractions of each ray for csf
    REAL(wp), DIMENSION(nrays), INTENT(IN)  ::  zdirs         !< list of z directions to raytrace (z/hdist, grid, zenith->nadir)
    REAL(wp), DIMENSION(nrays)              ::  zstop         !< fractional Z grid-coordinate of ray end or -999.0 if at col bdry

#if defined( __parallel )
    INTEGER(iwp)              ::  lowest_lad  !< lowest column cell for which we need LAD
    INTEGER(iwp)              ::  wcount      !< RMA window item count
    INTEGER(MPI_ADDRESS_KIND) ::  wdisp       !< RMA window displacement
#endif


    yxorigin(:) = origin(2:3)
    transparency(:) = 1.0_wp !-- Pre-set the all rays to transparent before reducing
    full_horizon = -HUGE( 1.0_wp )
    lowest_mixed_ray = nrays
    lowest_free_ray = nrays

    ALLOCATE( target_surfl(nrays) )
    target_surfl(:) = -1
    ray_ntrack(:) = HUGE( 0_iwp )
    lastdir = -999
    lastcolumn(:) = -999
    ALLOCATE( vert_col(nz_urban_b:nz_urban_t) )
    ALLOCATE( up_col(nz_urban_b:nz_urban_t) )
    ALLOCATE( down_col(nz_urban_b:nz_urban_t) )
    glob_max(1) = ny
    glob_max(2) = nx

    IF ( plant_canopy )  THEN
       rt2_track_dist(0) = 0.0_wp
       rt2_track_lad(:,:) = 0.0_wp
       nly = plantt_max - nz_urban_b + 1
    ENDIF

    ip_last = -1
    ig_last = -1
    lastdist = 0.0_wp

!
!-- Since all face coordinates have values *.5 and we'd like to use integers, all these have
!-- 0.5 added
    DO  d = 1, 2
       IF ( yxdir(d) == 0.0_wp )  THEN
          dimnext(d) = HUGE( 1 )
          dimdelta(d) = HUGE( 1 )
          dimnextdist(d) = HUGE( 1.0_wp )
          column(d) = NINT( yxorigin(d) )
       ELSE IF ( yxdir(d) > 0.0_wp )  THEN
          dimnext(d) = FLOOR( yxorigin(d) + 0.5_wp ) + 1
          dimdelta(d) = 1
          dimnextdist(d) = ( dimnext(d) - 0.5_wp - yxorigin(d) ) / yxdir(d)
          column(d) = dimnext(d) - 1
       ELSE
          dimnext(d) = CEILING( yxorigin(d) + 0.5_wp ) - 1
          dimdelta(d) = -1
          dimnextdist(d) = ( dimnext(d) - 0.5_wp - yxorigin(d) ) / yxdir(d)
          column(d) = dimnext(d)
       ENDIF
    ENDDO

    ntrack = 0
    DO
!
!--    Along what dimension will the next wall crossing be?
       seldim = MINLOC( dimnextdist, 1 )
       nextdist = dimnextdist(seldim)

       IF ( nextdist > lastdist )  THEN
          ntrack = ntrack + 1
!
!--       Calculate index of the grid with global indices (column(1),column(2)) in the array
!--       nzterrt/b and plantt and id of the coresponding processor
          CALL radiation_calc_global_offset( column(2), column(1), 0, 1, iproc = ip,               &
                                             offs_glob = ig )

          IF ( ip_last < 0 )  THEN
             horz_full_entry  = -HUGE( 1.0_wp )
             horz_mixed_entry = -HUGE( 1.0_wp )
          ELSE
             horz_full_entry  = ( REAL( nzterrb(ig), wp ) + 0.5_wp - origin(1) ) / lastdist
             horz_mixed_entry = ( REAL( nzterrt(ig), wp ) + 0.5_wp - origin(1) ) / lastdist
          ENDIF
          horz_full_exit  = ( REAL( nzterrb(ig), wp ) + 0.5_wp - origin(1) ) / nextdist
          horz_mixed_exit = ( REAL( nzterrt(ig), wp ) + 0.5_wp - origin(1) ) / nextdist
          is_mixed_col = ( nzterrt(ig) /= nzterrb(ig) )
!
!--       Identify vertical full obstacles hit by rays in current column, mixed rays need to be
!--       checked if they are already obstructed.
          DO WHILE ( lowest_mixed_ray > lowest_free_ray )
             IF ( zdirs(lowest_mixed_ray) > horz_full_entry )  EXIT
!
!--          This may only happen after 1st column, so lastdir and lastcolumn are valid
             IF ( target_surfl(lowest_mixed_ray) < 0 )  THEN
                CALL request_itarget( lastdir, CEILING( -0.5_wp + origin(1) +                   &
                                                        zdirs(lowest_mixed_ray) * lastdist ),   &
                                      lastcolumn(1), lastcolumn(2),                             &
                                      target_surfl(lowest_mixed_ray) )
                ray_ntrack(lowest_mixed_ray) = ntrack - 1
                zstop(lowest_mixed_ray) = -999.0_wp
             ENDIF
             lowest_mixed_ray = lowest_mixed_ray - 1
          ENDDO
!
!--       Identify vertical full obstacles hit by rays in current column, free rays need no
!--       individual checks.
          DO WHILE ( lowest_mixed_ray > 0 )
             IF ( zdirs(lowest_mixed_ray) > horz_full_entry )  EXIT
!
!--          This may only happen after 1st column, so lastdir and lastcolumn are valid
             CALL request_itarget( lastdir, CEILING( -0.5_wp + origin(1) +                      &
                                                     zdirs(lowest_mixed_ray) * lastdist ),      &
                                   lastcolumn(1), lastcolumn(2), target_surfl(lowest_mixed_ray) )
             ray_ntrack(lowest_mixed_ray) = ntrack - 1
             zstop(lowest_mixed_ray) = -999.0_wp
             lowest_mixed_ray = lowest_mixed_ray - 1
          ENDDO
          IF ( lowest_free_ray > lowest_mixed_ray )  lowest_free_ray = lowest_mixed_ray
!
!--       Identify targets for vertical mixed obstacles.
!--       lowest_mixed_ray now points to bottom of vertical mixed obstacles.
          IF ( is_mixed_col  .AND.  ip_last >= 0 )  THEN
!
!--          Load vertical surfaces belonging to previous column
             vert_col(:) = -999
             DO  isurf = surfg_col_start(ig_last), surfg_col_start(ig_last+1)-1
                IF ( surf(id, isurf) == lastdir )  THEN
                   vert_col(surf(iz, isurf)) = isurf
                ENDIF
             ENDDO
!
!--          Previously mixed rays need to be checked whether they are obstructed
             DO  iray = lowest_mixed_ray, lowest_free_ray+1, -1
                IF ( zdirs(iray) > horz_mixed_entry )  EXIT
                IF ( target_surfl(iray) >= 0 )  CYCLE
                target_surfl(iray) = vert_col( CEILING( -0.5_wp + origin(1) +                      &
                                                        zdirs(iray) * lastdist) )  ! Contains -999 if missing surface
                IF ( target_surfl(iray) >= 0 )  THEN
                   ray_ntrack(lowest_mixed_ray) = ntrack - 1
                   zstop(lowest_mixed_ray) = -999.0_wp
                ENDIF
             ENDDO
!
!--          Previously free rays cannot be obstructed yet
             iray = lowest_free_ray
             DO WHILE ( iray >= 1 )
                IF ( zdirs(iray) > horz_mixed_entry )  EXIT
                target_surfl(iray) = vert_col( CEILING( -0.5_wp + origin(1) + zdirs(iray) *     &
                                                        lastdist ) )  ! Contains -999 if missing surface
                IF ( target_surfl(iray) >= 0 )  THEN
                   ray_ntrack(lowest_mixed_ray) = ntrack - 1
                   zstop(lowest_mixed_ray) = -999.0_wp
                ENDIF
                iray = iray - 1
             ENDDO
!
!--          Extend mixed rays by raising the lowest_free ray (remains unchanged if the previous
!--          loop exited immediately, becomes 0 if all rays are full/mixed)
             lowest_free_ray = iray
          ENDIF  ! End of mixed horizon
!
!--       Identify horizontal full obstacles hit by rays in current column, mixed rays need to be
!--       checked if they are already obstructed.
          DO WHILE ( lowest_mixed_ray > lowest_free_ray )
             IF ( zdirs(lowest_mixed_ray) > horz_full_exit )  EXIT
             IF ( target_surfl(lowest_mixed_ray) < 0 )  THEN
                CALL request_itarget( iup, nzterrb(ig)+1, column(1), column(2),                 &
                                      target_surfl(lowest_mixed_ray) )
                ray_ntrack(lowest_mixed_ray) = ntrack
                zstop(lowest_mixed_ray) = REAL( nzterrb(ig), wp ) + 0.5_wp
             ENDIF
             lowest_mixed_ray = lowest_mixed_ray - 1
          ENDDO
!
!--       Identify horizontal full obstacles hit by rays in current column, free rays need no
!--       individual checks.
          DO WHILE ( lowest_mixed_ray > 0 )
             IF ( zdirs(lowest_mixed_ray) > horz_full_exit )  EXIT
             CALL request_itarget( iup, nzterrb(ig)+1, column(1), column(2),                    &
                                   target_surfl(lowest_mixed_ray) )
             ray_ntrack(lowest_mixed_ray) = ntrack
             zstop(lowest_mixed_ray) = REAL( nzterrb(ig), wp ) + 0.5_wp
             lowest_mixed_ray = lowest_mixed_ray - 1
          ENDDO
          IF ( lowest_free_ray > lowest_mixed_ray )  lowest_free_ray = lowest_mixed_ray
!
!--       Identify targets for horizontal mixed obstacles.
!--       lowest_mixed_ray now points _above_ horizontal full obstacles.
          IF ( is_mixed_col )  THEN
!
!--          Load horizontal surfaces corresponding to current column
             up_col(:) = - 999
             down_col(:) = - 999
             DO  isurf = surfg_col_start(ig), surfg_col_start(ig+1)-1
                SELECT CASE ( surf(id, isurf) )
                CASE ( iup )
                   up_col(surf(iz, isurf)) = isurf
                CASE ( idown )
                   down_col(surf(iz, isurf)) = isurf
                ENDSELECT
             ENDDO
!
!--          Previously mixed rays need to be checked whether they are obstructed
             DO  iray = lowest_mixed_ray, lowest_free_ray+1, -1
                IF ( zdirs(iray) > MAX( horz_mixed_entry, horz_mixed_exit ) )  EXIT
                IF ( target_surfl(iray) >= 0 )  CYCLE
                IF ( zdirs(iray) <= 0.0_wp )  THEN
!
!--                Downward pointed ray, cycle k down from entry to exit, search for upward
!--                oriented faces
                   DO  k = FLOOR( 0.5_wp + origin(1) + zdirs(iray) * lastdist ),                &
                           CEILING( -0.5_wp + origin(1) + zdirs(iray) * nextdist ) + 1, - 1
                      target_surfl(iray) = up_col(k)  ! Contains -999 if missing surface
                      IF ( target_surfl(iray) >= 0 )  THEN
                         ray_ntrack(lowest_mixed_ray) = ntrack
                         zstop(lowest_mixed_ray) = REAL( k, wp ) - 0.5_wp
                         EXIT
                      ENDIF
                   ENDDO
                ELSE
!
!--                Upward pointed ray, cycle k up from entry to exit, search for downward
!--                oriented faces
                   DO  k = CEILING( -0.5_wp + origin(1) + zdirs(iray) * lastdist ),             &
                           FLOOR( 0.5_wp + origin(1) + zdirs(iray) * nextdist ) - 1
                      target_surfl(iray) = down_col(k)  ! Contains -999 if missing surface
                      IF ( target_surfl(iray) >= 0 )  THEN
                         ray_ntrack(lowest_mixed_ray) = ntrack
                         zstop(lowest_mixed_ray) = REAL( k, wp ) + 0.5_wp
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
!
!--          Previously free rays cannot be obstructed yet
             iray = lowest_free_ray
             DO WHILE ( iray >= 1 )
                IF ( zdirs(iray) > MAX( horz_mixed_entry, horz_mixed_exit ) )  EXIT
                IF ( zdirs(iray) <= 0.0_wp )  THEN
!
!--                Downward pointed ray, cycle k down from entry to exit, search upward oriented
!--                faces
                   DO  k = FLOOR( 0.5_wp + origin(1) + zdirs(iray) * lastdist ),                &
                           CEILING( -0.5_wp + origin(1) + zdirs(iray) * nextdist) + 1, - 1
                      target_surfl(iray) = up_col(k)  ! Contains -999 if missing surface
                      IF ( target_surfl(iray) >= 0 )  THEN
                         ray_ntrack(lowest_mixed_ray) = ntrack
                         zstop(lowest_mixed_ray) = REAL( k, wp ) - 0.5_wp
                         EXIT
                      ENDIF
                   ENDDO
                ELSE
!
!--                Upward pointed ray, cycle k up from entry to exit, search downward oriented
!--                faces
                   DO  k = CEILING( -0.5_wp + origin(1) + zdirs(iray) * lastdist ),             &
                           FLOOR( 0.5_wp + origin(1) + zdirs(iray) * nextdist ) - 1
                      target_surfl(iray) = down_col(k)  ! Contains -999 if missing surface
                      IF ( target_surfl(iray) >= 0 )  THEN
                         ray_ntrack(lowest_mixed_ray) = ntrack
                         zstop(lowest_mixed_ray) = REAL( k, wp ) + 0.5_wp
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF
                iray = iray - 1
             ENDDO
!
!--          Extend mixed rays by raising the lowest_free ray (remains unchanged if the previous
!--          loop exited immediately, becomes 0 if all rays are full/mixed)
             lowest_free_ray = iray
          ENDIF  ! End of mixed horizon
          full_horizon = MAX( full_horizon, horz_full_entry, horz_full_exit )

          IF ( plant_canopy )  THEN
             rt2_track(:, ntrack) = column(:)
             rt2_track_dist(ntrack) = nextdist
          ENDIF
       ENDIF

       lastcolumn(:) = column(:)
       column(seldim) = column(seldim) + dimdelta(seldim)
       IF ( column(seldim) < 0  .OR.  column(seldim) > glob_max(seldim) )  EXIT
!
!--    If all (potentially) free rays exited the urban layer upwards, stop
       IF ( lowest_mixed_ray <= 0 )  EXIT
       IF ( origin(1) + zdirs(lowest_mixed_ray) * nextdist >= REAL(nz_urban_t, wp) + 0.5_wp )  EXIT
!
!--    Save wall direction of coming building column (= this air column)
       IF ( seldim == 1 )  THEN
          IF ( dimdelta(seldim) == 1 )  THEN
             lastdir = isouth
          ELSE
             lastdir = inorth
          ENDIF
       ELSE
          IF ( dimdelta(seldim) == 1 )  THEN
             lastdir = iwest
          ELSE
             lastdir = ieast
          ENDIF
       ENDIF

       ip_last = ip
       ig_last = ig
       lastdist = nextdist
       dimnext(seldim) = dimnext(seldim) + dimdelta(seldim)
       dimnextdist(seldim) = ( dimnext(seldim) - 0.5_wp - yxorigin(seldim) ) / yxdir(seldim)
    ENDDO

    IF ( plant_canopy )  THEN
!
!--    Request LAD WHERE applicable

#if defined( __parallel )
       IF ( raytrace_mpi_rma )  THEN
!
!--       Send requests for lad_s to appropriate processor
          DO  i = 1, ntrack
             CALL radiation_calc_global_offset( rt2_track(2,i), rt2_track(1,i), 0, 1,              &
                                                offs_glob = ig )

             IF ( calc_svf )  THEN
!
!--             For fixed view resolution, we need plant canopy even for rays to opposing surfaces
                lowest_lad = nzterrb(ig) + 1
             ELSE
!
!--             We only need LAD for rays directed above full horizon (to sky)
                lowest_lad = CEILING( -0.5_wp + origin(1) +                                        &
                                      MIN( full_horizon * rt2_track_dist(i-1),  &  ! Entry
                                           full_horizon * rt2_track_dist(i) ) )  ! Exit
             ENDIF
!
!--          Skip asking for LAD where all plant canopy is under requested level
             IF ( plantt(ig) < lowest_lad )  CYCLE

             CALL radiation_calc_global_offset( rt2_track(2,i), rt2_track(1,i),                    &
                                                lowest_lad - nz_urban_b, nz_plant, iproc = ip,     &
                                                offs_proc = wdisp )
             wcount = plantt(ig) - lowest_lad + 1
!
!--          TODO: send request ASAP - even during raytracing
             CALL MPI_GET( rt2_track_lad(lowest_lad:plantt(ig), i), wcount, MPI_REAL, ip, wdisp,   &
                           wcount, MPI_REAL, win_lad, ierr )
             IF ( ierr /= 0  .AND.  debug_output )  THEN
                WRITE( debug_string, * ) 'Error MPI_GET2:', ierr,                                  &
                                         rt2_track_lad(lowest_lad:plantt(ig), i), wcount, ip,      &
                                         wdisp, win_lad
                CALL debug_message( debug_string, 'info' )
             ENDIF
          ENDDO
!
!--       Wait for all pending local requests to complete
!--       TODO: Wait selectively for each column later when needed
          CALL MPI_WIN_FLUSH_LOCAL_ALL( win_lad, ierr )
          IF ( ierr /= 0  .AND.  debug_output )  THEN
             WRITE( debug_string, * ) 'Error MPI_WIN_FLUSH_LOCAL_ALL2:', ierr, win_lad
             CALL debug_message( debug_string, 'info' )
          ENDIF

       ELSE  ! raytrace_mpi_rma = .F.
          DO  i = 1, ntrack
             CALL radiation_calc_global_offset( rt2_track(2,i), rt2_track(1,i), 0, nz_plant,       &
                                                offs_glob = ig )
             rt2_track_lad(nz_urban_b:plantt_max, i) = sub_lad_g(ig:ig+nly-1)
          ENDDO
       ENDIF
#else
       DO  i = 1, ntrack
          rt2_track_lad(nz_urban_b:plantt_max, i) =                                                &
                                      sub_lad(rt2_track(1,i),rt2_track(2,i),nz_urban_b:plantt_max)
       ENDDO
#endif
    ENDIF  ! plant_canopy

#if defined( __parallel )
!
!-- Wait for all gridsurf requests to complete
    CALL MPI_WIN_FLUSH_LOCAL_ALL( win_gridsurf, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
       WRITE( debug_string, * ) 'Error MPI_WIN_FLUSH_LOCAL_ALL3:', ierr, win_gridsurf
       CALL debug_message( debug_string, 'info' )
    ENDIF
#endif
    itarget(:) = target_surfl(:)
    DEALLOCATE( target_surfl )

    IF ( plant_canopy )  THEN
!
!--    Skip the PCB around origin if requested (for MRT, the PCB might not be there)
       IF ( skip_1st_pcb  .AND.  NINT( origin(1) ) <= plantt_max )  THEN
          rt2_track_lad( NINT( origin(1), iwp ), 1 ) = 0.0_wp
       ENDIF
!
!--    Assert that we have space allocated for CSFs
       maxboxes = ( ntrack + MAX( CEILING( origin(1) - 0.5_wp ) - nz_urban_b, nz_urban_t -         &
                                  CEILING( origin(1) - 0.5_wp ) ) ) * nrays
       IF ( calc_svf )  maxboxes = maxboxes * 2 ! Two passes for each ray
       IF ( ncsfl + maxboxes > ncsfla )  THEN
!
!--       Use this code for growing by fixed exponential increments (equivalent to case where ncsfl
!--       always increases by 1)
!--       k = CEILING(grow_factor ** real(CEILING(log(real(ncsfl + maxboxes, KIND=wp)) &
!--                                              / log(grow_factor)), KIND=wp))
!--       Or use this code to simply always keep some extra space after growing
          k = CEILING( REAL( ncsfl + maxboxes, KIND = wp ) * grow_factor )
          CALL merge_and_grow_csf(k)
       ENDIF
!
!--    Calculate transparencies and store new CSFs
       zbottom = REAL( nz_urban_b, wp ) - 0.5_wp
       ztop = REAL( plantt_max, wp ) + 0.5_wp
!
!--    Reverse direction of radiation (face->sky), only when calc_svf
       IF ( calc_svf )  THEN
          DO  i = 1, ntrack  ! For each column
             dxxyy = ( ( dy * yxdir(1) )**2 + ( dx * yxdir(2) )**2 ) *                             &
                     ( rt2_track_dist(i) - rt2_track_dist(i-1) )**2
             CALL radiation_calc_global_offset( rt2_track(2,i), rt2_track(1,i), 0, 1, iproc = ip )

             DO  k = 1, nrays  ! For each ray
!
!--             NOTE 6778:
!--             With traditional svf discretization, CSFs under the horizon (i.e. for surface to
!--             surface radiation) were created in raytrace(). With angular discretization, we
!--             must create CSFs under horizon only for one direction, otherwise we would have
!--             duplicate amount of energy. Although we could choose either of the two directions
!--             (they differ only by discretization error with no bias), we choose the backward
!--             direction, because it tends to cumulate high canopy sink factors closer to raytrace
!--             origin, i.e. it should potentially cause less moiree.
                IF ( i > ray_ntrack(k) )  CYCLE
                zorig = origin(1) + zdirs(k) * rt2_track_dist(i-1)
!
!--             zorig is always nearer to origin than zexit
                IF ( zorig <= zbottom  .OR.  zorig >= ztop )  CYCLE

                zsgn = INT( SIGN( 1.0_wp, zdirs(k) ), iwp )
                rt2_dist(1) = 0.0_wp
                IF ( zdirs(k) == 0.0_wp )  THEN  ! Ray is exactly horizontal
                   nz = 2
                   rt2_dist(nz) = SQRT( dxxyy )
                   kz = CEILING( -0.5_wp + zorig, iwp )
                ELSE
                   IF ( i == ray_ntrack(k)  .AND.  zstop(k) >= 0.0_wp )  THEN
                      zexit = zstop(k)
                   ELSE
                      zexit = MIN( MAX( origin(1) + zdirs(k) * rt2_track_dist(i), zbottom ), ztop )
                   ENDIF
                   zb0 = FLOOR( zorig * zsgn - 0.5_wp ) + 1  ! Because it must be greater than orig
                   zb1 = CEILING( zexit * zsgn - 0.5_wp ) - 1  ! Because it must be smaller than exit
                   nz = MAX( zb1 - zb0 + 3, 2 )
                   rt2_dist(nz) = SQRT( ( (zexit - zorig ) * dz(1) )**2 + dxxyy )
                   qdist = rt2_dist(nz) / ( zexit - zorig )
                   rt2_dist(2:nz-1) = (/ ( ( ( REAL( l, wp ) + 0.5_wp ) * zsgn - zorig ) * qdist,  &
                                         l = zb0, zb1 ) /)
                   kz = zb0 * zsgn
                ENDIF

                DO  l = 2, nz
                   IF ( rt2_track_lad(kz, i) > 0.0_wp )  THEN
                      curtrans = EXP( - ext_coef * rt2_track_lad(kz, i) * ( rt2_dist(l) -          &
                                                                            rt2_dist(l-1) ) )

                      IF ( create_csf )  THEN
                         ncsfl = ncsfl + 1
                         acsf(ncsfl)%ip = ip
                         acsf(ncsfl)%itx = rt2_track(2,i)
                         acsf(ncsfl)%ity = rt2_track(1,i)
                         acsf(ncsfl)%itz = kz
                         acsf(ncsfl)%isurfs = iorig
                         acsf(ncsfl)%rcvf = ( 1.0_wp - curtrans ) * transparency(k) * vffrac(k)
                      ENDIF

                      transparency(k) = transparency(k) * curtrans
                   ENDIF
                   kz = kz + zsgn
                ENDDO  ! l = 1, nz - 1
             ENDDO  ! k = 1, nrays
          ENDDO  ! i = 1, ntrack
!
!--       Reset rays above horizon to transparent (see NOTE 6778)
          WHERE( itarget < 0 )  transparency = 1.0_wp
       ENDIF
!
!--    Forward direction of radiation (sky->face), always
       DO  i = ntrack, 1, -1  ! For each column backwards
          dxxyy = ( ( dy * yxdir(1) )**2 + ( dx * yxdir(2) )**2 ) *                                &
                  ( rt2_track_dist(i) - rt2_track_dist(i-1) )**2
          CALL radiation_calc_global_offset( rt2_track(2,i), rt2_track(1,i), 0, 1, iproc = ip )

          DO  k = 1, nrays  ! For each ray
!
!--          See NOTE 6778 above
             IF ( itarget(k) >= 0 )  CYCLE
             IF ( i > ray_ntrack(k) )  CYCLE

             zexit = origin(1) + zdirs(k) * rt2_track_dist(i-1)
             IF ( zexit <= zbottom  .OR.  zexit >= ztop )  CYCLE

             zsgn = - INT( SIGN( 1.0_wp, zdirs(k) ), iwp )
             rt2_dist(1) = 0.0_wp
             IF ( zdirs(k) == 0.0_wp )  THEN  ! Ray is exactly horizontal
                nz = 2
                rt2_dist(nz) = SQRT( dxxyy )
                kz = NINT( zexit, iwp )
             ELSE
                IF ( i == ray_ntrack(k)  .AND.  zstop(k) >= 0.0_wp )  THEN
                   zorig = zstop(k)
                ELSE
                   zorig = MIN( MAX( origin(1) + zdirs(k) * rt2_track_dist(i), zbottom ), ztop )
                ENDIF
                zb0 = FLOOR( zorig * zsgn - 0.5_wp ) + 1  ! Because it must be greater than orig
                zb1 = CEILING( zexit * zsgn - 0.5_wp ) - 1  ! Because it must be smaller than exit
                nz = MAX( zb1 - zb0 + 3, 2 )
                rt2_dist(nz) = SQRT( ( ( zexit - zorig ) * dz(1) )**2 + dxxyy )
                qdist = rt2_dist(nz) / ( zexit - zorig )
                rt2_dist(2:nz-1) = (/ ( ( ( REAL( l, wp ) + 0.5_wp ) * zsgn - zorig ) * qdist,     &
                                      l = zb0, zb1 ) /)
                kz = zb0 * zsgn
             ENDIF

             DO  l = 2, nz
                IF ( rt2_track_lad(kz, i) > 0.0_wp )  THEN
                   curtrans = EXP( - ext_coef * rt2_track_lad(kz, i) *                             &
                                  ( rt2_dist(l) - rt2_dist(l-1) ) )

                   IF ( create_csf )  THEN
                      ncsfl = ncsfl + 1
                      acsf(ncsfl)%ip = ip
                      acsf(ncsfl)%itx = rt2_track(2,i)
                      acsf(ncsfl)%ity = rt2_track(1,i)
                      acsf(ncsfl)%itz = kz
                      acsf(ncsfl)%isurfs = -1
                      acsf(ncsfl)%rcvf = ( 1.0_wp - curtrans ) * transparency(k) * aorig * vffrac(k)
                   ENDIF   ! create_csf

                   transparency(k) = transparency(k) * curtrans
                ENDIF
                kz = kz + zsgn
             ENDDO  ! l = 1, nz - 1
          ENDDO  ! k = 1, nrays
       ENDDO  ! i = 1, ntrack
    ENDIF  ! plant_canopy

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Use MPI one-sided operation Get to request index of target surface at the point of intersection
!> of a ray with a full obstacle (terrain, building). Although it does not assign the value
!> immediately (only after MPI_WIN_FLUSH_ALL), it is used to assign targets that are (become) below
!> the `lowest_mixed_ray` index, therefore its value is not tested again within `raytrace_2d`.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE request_itarget( d, z, y, x, isurfl )

    INTEGER(iwp), INTENT(IN) ::  d, z, y, x  !<

    INTEGER(iwp), TARGET, INTENT(OUT) ::  isurfl  !<

#if defined( __parallel )
    INTEGER(iwp) ::  iproc  !<

    INTEGER(KIND=MPI_ADDRESS_KIND) ::  target_displ  !< index of the grid in the local gridsurf array

!
!-- Calculate target processor and index in the remote local target gridsurf array
    CALL radiation_calc_global_offset( x, y, ( z - nz_urban_b ) * nsurf_type_u + d,                &
                                       nz_urban * nsurf_type_u, iproc = iproc,                     &
                                       offs_proc = target_displ )
!
!-- Send MPI_GET request to obtain index target_surfl(i)
    CALL MPI_GET( isurfl, 1, MPI_INTEGER, iproc, target_displ, 1, MPI_INTEGER, win_gridsurf, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
       WRITE( debug_string, * ) 'Error MPI_GET3:', ierr, isurfl, iproc, target_displ, win_gridsurf
       CALL debug_message( debug_string, 'info' )
    ENDIF
#else
!
!-- Set index target_surfl(i)
    isurfl = gridsurf(d,z,y,x)
#endif

 END SUBROUTINE request_itarget

 END SUBROUTINE raytrace_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Part of new localized raytracing. Initializes raytracing of an angular section and sends the
!> first raytracing request message (to itself, because raytracing always starts in the local
!> process).
!> For a brief description of the algorighm, see
!> https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/rtm#Localizedraytracingparallelizationscheme
!> and https://gitlab.palm-model.org/palm/model/-/merge_requests/247
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE raytrace_init( task, origin, iaz, nrays, zdirs, iorig, aorig, vffrac, s, r, dtype )

    IMPLICIT NONE

    CHARACTER, INTENT(IN) ::  task !< raytracing task (used also as message code)

    INTEGER(iwp)             ::  d      !< direction index
    INTEGER(iwp)             ::  iray   !< ray index

    INTEGER(iwp), INTENT(IN) ::  dtype  !< MPI datatype associated with the message
    INTEGER(iwp), INTENT(IN) ::  iaz    !< discretized azimuth index
    INTEGER(iwp), INTENT(IN) ::  iorig  !< index of origin face for csf
    INTEGER(iwp), INTENT(IN) ::  nrays  !< number of rays (z directions) to raytrace

    REAL(wp), INTENT(IN) ::  aorig   !< origin face area for csf

    REAL(wp), DIMENSION(2)                 ::  yxorigin  !< horizontal copy of origin (y,x)
    REAL(wp), DIMENSION(3), INTENT(IN)     ::  origin    !< z,y,x coordinates of ray origin
    REAL(wp), DIMENSION(nrays), INTENT(IN) ::  vffrac    !< view factor fractions of each ray for csf
    REAL(wp), DIMENSION(nrays), INTENT(IN) ::  zdirs     !< list of z directions to raytrace (z/hdist in grid,
                                                         !< zenith->nadir)

    TYPE(t_traced_section), INTENT(INOUT)           ::  s  !< raytraced section message structure
    TYPE(t_traced_ray), DIMENSION(:), INTENT(INOUT) ::  r  !< array of ray message structures


!
!-- Fill in raytracing section structure
    s%message = task
    s%iorig = iorig
    s%lowest_free_ray = nrays
    s%lowest_mixed_ray = nrays
    s%nrays = nrays
    s%aorig = aorig
    s%lastdist = 0.0_wp
    s%origin = origin
    s%iaz = iaz

    DO  iray = 1, nrays
       r(iray)%itarget = -1
       r(iray)%vffrac = vffrac(iray)
       r(iray)%transp = 1.0_wp
       r(iray)%zdir = zdirs(iray)
    ENDDO
!
!-- Since all face coordinates have values *.5 and we'd like to use integers, all these have
!-- 0.5 added
    yxorigin(:) = origin(2:3)
    DO  d = 1, 2
       IF ( discr_azim_yxdir(d,iaz) == 0.0_wp )  THEN
          s%dimnext(d) = HUGE( 1 )
          s%dimnextdist(d) = HUGE( 1.0_wp )
       ELSEIF ( discr_azim_yxdir(d,iaz) > 0.0_wp )  THEN
          s%dimnext(d) = FLOOR( yxorigin(d) + 0.5_wp ) + 1
          s%dimnextdist(d) = ( s%dimnext(d) - 0.5_wp - yxorigin(d) ) / discr_azim_yxdir(d,iaz)
       ELSE
          s%dimnext(d) = CEILING( yxorigin(d) + 0.5_wp ) - 1
          s%dimnextdist(d) = ( s%dimnext(d) - 0.5_wp - yxorigin(d) ) / discr_azim_yxdir(d,iaz)
       ENDIF
    ENDDO
!
!-- Raytracing starts in the local process, so we call the subroutine directly
    CALL raytrace_segment( s, r, dtype )

 END SUBROUTINE raytrace_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Part of new localized raytracing. Copies data from raytracing message after the last segment has
!> been processed.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE raytrace_finalize( s, r )

    IMPLICIT NONE

    INTEGER(iwp) ::  i        !< solar direction index
    INTEGER(iwp) ::  iorig_l  !< local index of source face
    INTEGER(iwp) ::  iray     !< ray index

    TYPE(t_traced_section), INTENT(INOUT)           ::  s  !< raytraced section message structure
    TYPE(t_traced_ray), DIMENSION(:), INTENT(INOUT) ::  r  !< array of ray message structures


!
!-- Save regular view-factor outputs for SVF and MRTF
    IF ( s%message == 'F'  .OR.  s%message == 'M' )  THEN
       DO  iray = 1, s%nrays
          nsvf_ins = nsvf_ins + 1
          itarget(nsvf_ins) = r(iray)%itarget
          vffrac(nsvf_ins) = r(iray)%vffrac
          ztransp(nsvf_ins) = r(iray)%transp
       ENDDO
    ENDIF
!
!-- Save sky-view factors and direct solar transparency except for downward oriented surfaces.
!-- Later for  for slanted surfaces, we will use generic zdir approach
    IF ( s%message == 'F'  .AND.  r(1)%zdir > 0.0_wp )  THEN
       iorig_l = s%iorig - surfstart(myid)
       DO  iray = 1, s%nrays
          IF ( r(iray)%itarget < 0 )  THEN
             skyvf(iorig_l) = skyvf(iorig_l) + r(iray)%vffrac
             skyvft(iorig_l) = skyvft(iorig_l) + r(iray)%vffrac * r(iray)%transp

             IF ( iray <= raytrace_discrete_elevs / 2 )  THEN
                i = dsidir_rev(iray-1, s%iaz-1)
                IF ( i /= -1 )  dsitrans(iorig_l, i) = r(iray)%transp
             ENDIF
          ENDIF
       ENDDO
    ENDIF
!
!-- Save direct solar transparency for PCGB
    IF ( s%message == 'P' )  THEN
       DO  iray = 1, s%nrays
          IF ( r(iray)%itarget < 0 )  THEN
             i = dsidir_rev(iray-1, s%iaz-1)
             IF ( i /= -1 )  dsitransc(s%iorig, i) = r(iray)%transp
          ENDIF
       ENDDO
       nsvf_ins = nsvf_ins + 1 ! counting only azimuths here
    ELSEIF ( s%message == 'M' )  THEN
       DO  iray = 1, s%nrays
          IF ( r(iray)%itarget < 0 )  THEN
             mrtsky(s%iorig) = mrtsky(s%iorig) + r(iray)%vffrac
             mrtskyt(s%iorig) = mrtskyt(s%iorig) + r(iray)%vffrac * r(iray)%transp

             IF ( iray <= raytrace_discrete_elevs / 2 )  THEN
                i = dsidir_rev(iray-1, s%iaz-1)
                IF ( i /= -1 )  mrtdsit(s%iorig, i) = r(iray)%transp
             ENDIF
          ENDIF
       ENDDO
    ENDIF

 END SUBROUTINE raytrace_finalize


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Part of new localized raytracing. Performs 2D raytracing of a segment of rays that passes the
!> current subdomain.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE raytrace_segment( s, r, dtype )

    IMPLICIT NONE

    INTEGER(iwp) ::  d          !< dimension index
    INTEGER(iwp) ::  iray       !< index into zdirs
    INTEGER(iwp) ::  k, l       !< running indices
    INTEGER(iwp) ::  kz         !< z running index
    INTEGER(iwp) ::  lastdir    !< wall direction before hitting this column
    INTEGER(iwp) ::  nz         !< number of layers for plant canopy
    INTEGER(iwp) ::  seldim     !< dimension to be incremented
    INTEGER(iwp) ::  topo_full  !< full topography height (below full-3d) at current column
    INTEGER(iwp) ::  topo_mix   !< mixed topography height (above full-3d) at current column
    INTEGER(iwp) ::  zb0        !< z boundary for plant canopy
    INTEGER(iwp) ::  zb1        !< z boundary for plant canopy
    INTEGER(iwp) ::  zsgn       !< sign for vertical direction for plant canopy

    INTEGER(iwp), INTENT(IN) ::  dtype  !< MPI datatype associated with the message

    INTEGER(iwp), DIMENSION(2) ::  column      !< grid column being crossed
    INTEGER(iwp), DIMENSION(2) ::  lastcolumn  !< previous grid column
    INTEGER(iwp), DIMENSION(2) ::  dimdelta    !< dimension direction = +- 1

    LOGICAL ::  backward      !< direction of raytracing
    LOGICAL ::  create_csf    !< whether to create canopy sink factors
    LOGICAL ::  is_mixed_col  !< whether current column contains full-3D geometry

    REAL(wp) ::  curtrans          !< transparency of current PC box crossing
    REAL(wp) ::  dxxyy             !< square of real horizontal distance
    REAL(wp) ::  horz_full_entry   !< full horizon at entry to column
    REAL(wp) ::  horz_full_exit    !< full horizon at exit from column
    REAL(wp) ::  horz_mixed_entry  !< mixed horizon at entry to column
    REAL(wp) ::  horz_mixed_exit   !< mixed horizon at exit from column
    REAL(wp) ::  nextdist          !< end of current crossing
    REAL(wp) ::  qdist             !< ratio of real distance to z coord difference
    REAL(wp) ::  zentry            !< z coordinate of ray column entry
    REAL(wp) ::  zexit             !< z coordinate of ray column exit
    REAL(wp) ::  ztop              !< plant canopy top

    REAL(wp), DIMENSION(2) ::  glob_max  !< global domain boundary (north y, east x)
    REAL(wp), DIMENSION(2) ::  loc_max   !< local domain boundary (north y, east x)
    REAL(wp), DIMENSION(2) ::  loc_min   !< local domain boundary (south y, west x)
    REAL(wp), DIMENSION(2) ::  yxorigin  !< horizontal copy of origin (y,x)

    TYPE(t_traced_section), INTENT(INOUT)           ::  s  !< raytraced section message structure
    TYPE(t_traced_ray), DIMENSION(:), INTENT(INOUT) ::  r  !< array of ray message structures

#if defined( __parallel )
    INTEGER(iwp) ::  next_proc  !< next process that will continue the raytracing
#endif


    backward = s%message == 'b'
    create_csf = ( s%message == 'f'  .OR.  s%message == 'b' )
    yxorigin(:) = s%origin(2:3)
    rt2_dist(1) = 0.0_wp
    glob_max(1) = ny
    loc_max(1) = nyn
    loc_min(1) = nys
    glob_max(2) = nx
    loc_max(2) = nxr
    loc_min(2) = nxl

    DO  d = 1, 2
       IF ( discr_azim_yxdir(d,s%iaz) == 0.0_wp )  THEN
          dimdelta(d) = HUGE( 1 )
          column(d) = NINT(yxorigin(d))
       ELSEIF ( discr_azim_yxdir(d,s%iaz) > 0.0_wp )  THEN
          dimdelta(d) = 1
          IF ( backward )  THEN
             column(d) = s%dimnext(d)
          ELSE
             column(d) = s%dimnext(d) - 1
          ENDIF
       ELSE
          dimdelta(d) = -1
          IF ( backward )  THEN
             column(d) = s%dimnext(d) - 1
          ELSE
             column(d) = s%dimnext(d)
          ENDIF
       ENDIF
    ENDDO

    topo_full = topo_top_ind(column(1), column(2), 0)
    topo_mix  = topo_top_ind(column(1), column(2), 5)
    is_mixed_col = ( topo_mix /= topo_full )
    IF ( s%lastdist > 0.0_wp )  THEN
       horz_full_entry  = ( REAL( topo_full, wp ) + 0.5_wp - s%origin(1) ) / s%lastdist
       horz_mixed_entry = ( REAL( topo_mix,  wp ) + 0.5_wp - s%origin(1) ) / s%lastdist
    ELSE
       horz_full_entry  = -HUGE( 1.0_wp )
       horz_mixed_entry = -HUGE( 1.0_wp )
    ENDIF

    IF ( .NOT. backward )  THEN
       DO
!
!--       Along what dimension will the next wall crossing be?
          seldim = MINLOC( s%dimnextdist, 1 )
          nextdist = s%dimnextdist(seldim)

          horz_full_exit  = ( REAL( topo_full, wp ) + 0.5_wp - s%origin(1) ) / nextdist
          horz_mixed_exit = ( REAL( topo_mix,  wp ) + 0.5_wp - s%origin(1) ) / nextdist
          dxxyy = ( ( dy * discr_azim_yxdir(1,s%iaz) )**2 +                                        &
                    ( dx * discr_azim_yxdir(2,s%iaz) )**2 ) * ( nextdist - s%lastdist )**2
!
!--       Identify horizontal full obstacles hit by rays in current column, mixed rays need to be
!--       checked if they are already obstructed.
          DO WHILE ( s%lowest_mixed_ray > s%lowest_free_ray )
             IF ( r(s%lowest_mixed_ray)%zdir > horz_full_exit )  EXIT
             IF ( r(s%lowest_mixed_ray)%itarget < 0 )  THEN
                r(s%lowest_mixed_ray)%itarget = gridsurf(iup, topo_full+1, column(1), column(2))
                IF ( plant_canopy)  CALL process_ray_canopy( s%lowest_mixed_ray, .FALSE.,          &
                                                             REAL( topo_full, wp ) + 0.5_wp)
             ENDIF
             s%lowest_mixed_ray = s%lowest_mixed_ray - 1
          ENDDO
!
!--       Identify horizontal full obstacles hit by rays in current column, free rays need no
!--       individual checks.
          DO WHILE ( s%lowest_mixed_ray > 0 )
             IF ( r(s%lowest_mixed_ray)%zdir > horz_full_exit )  EXIT
             r(s%lowest_mixed_ray)%itarget = gridsurf(iup, topo_full+1, column(1), column(2))
             IF ( plant_canopy)  CALL process_ray_canopy( s%lowest_mixed_ray,  .FALSE.,            &
                                                          REAL( topo_full, wp) + 0.5_wp )
             s%lowest_mixed_ray = s%lowest_mixed_ray - 1
          ENDDO
          IF ( s%lowest_free_ray > s%lowest_mixed_ray )  s%lowest_free_ray = s%lowest_mixed_ray
!
!--       Identify targets for horizontal mixed obstacles.
!--       lowest_mixed_ray now points _above_ horizontal full obstacles.
          IF ( is_mixed_col )  THEN
!
!--          Previously mixed rays need to be checked whether they are obstructed
             iray = s%lowest_mixed_ray
             DO WHILE ( iray > s%lowest_free_ray )
                IF ( r(iray)%zdir > MAX( horz_mixed_entry, horz_mixed_exit ) )  THEN
!
!--                Above mixed horizon, just process canopy for remaining mixed rays
                   IF ( plant_canopy)  THEN
                      DO WHILE ( iray > s%lowest_free_ray )
                         IF ( r(iray)%itarget < 0 )  CALL process_ray_canopy( iray, .FALSE. )
                         iray = iray - 1
                      ENDDO
                   ENDIF
                   EXIT
                ENDIF
                IF ( r(iray)%itarget >= 0 )  THEN
                   iray = iray - 1
                   CYCLE
                ENDIF
                IF ( r(iray)%zdir <= 0.0_wp )  THEN
!
!--                Downward pointed ray, cycle k down from entry to exit, search for upward
!--                oriented faces
                   DO  k = FLOOR( 0.5_wp + s%origin(1) + r(iray)%zdir * s%lastdist ),              &
                           CEILING( -0.5_wp + s%origin(1) + r(iray)%zdir * nextdist ) + 1, - 1
                      r(iray)%itarget = gridsurf(iup, k, column(1), column(2))
                      IF ( r(iray)%itarget >= 0 )  THEN
                         IF ( plant_canopy )  CALL process_ray_canopy( iray, .FALSE.,              &
                                                                       REAL(k, wp) - 0.5_wp )
                         EXIT
                      ENDIF
                   ENDDO
                ELSE
!
!--                Upward pointed ray, cycle k up from entry to exit, search for downward
!--                oriented faces
                   DO  k = CEILING( -0.5_wp + s%origin(1) + r(iray)%zdir * s%lastdist ),           &
                           FLOOR( 0.5_wp + s%origin(1) + r(iray)%zdir * nextdist ) - 1
                      r(iray)%itarget = gridsurf(idown, k, column(1), column(2))
                      IF ( r(iray)%itarget >= 0 )  THEN
                         IF ( plant_canopy )  CALL process_ray_canopy( iray, .FALSE.,              &
                                                                       REAL(k, wp) + 0.5_wp )
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF
!
!--             If unstopped, then the whole ray is processed for canopy
                IF ( plant_canopy )  THEN
                   IF ( r(iray)%itarget < 0 )  CALL process_ray_canopy( iray, .FALSE. )
                ENDIF
                iray = iray - 1
             ENDDO
!
!--          Previously free rays cannot be obstructed yet
             iray = s%lowest_free_ray
             DO WHILE ( iray >= 1 )
                IF ( r(iray)%zdir > MAX( horz_mixed_entry, horz_mixed_exit ) )  EXIT
                IF ( r(iray)%zdir <= 0.0_wp )  THEN
!
!--                Downward pointed ray, cycle k down from entry to exit, search upward oriented
!--                faces
                   DO  k = FLOOR( 0.5_wp + s%origin(1) + r(iray)%zdir * s%lastdist ),              &
                           CEILING( -0.5_wp + s%origin(1) + r(iray)%zdir * nextdist) + 1, - 1
                      r(iray)%itarget = gridsurf(iup, k, column(1), column(2))
                      IF ( r(iray)%itarget >= 0 )  THEN
                         IF ( plant_canopy )  CALL process_ray_canopy( iray, .FALSE.,              &
                                                                       REAL(k, wp) - 0.5_wp )
                         EXIT
                      ENDIF
                   ENDDO
                ELSE
!
!--                Upward pointed ray, cycle k up from entry to exit, search downward oriented
!--                faces
                   DO  k = CEILING( -0.5_wp + s%origin(1) + r(iray)%zdir * s%lastdist ),           &
                           FLOOR( 0.5_wp + s%origin(1) + r(iray)%zdir * nextdist ) - 1
                      r(iray)%itarget = gridsurf(idown, k, column(1), column(2))
                      IF ( r(iray)%itarget >= 0 )  THEN
                         IF ( plant_canopy )  CALL process_ray_canopy( iray, .FALSE.,              &
                                                                       REAL(k, wp) + 0.5_wp )
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF
!
!--             If unstopped, then the whole ray is processed for canopy
                IF ( plant_canopy )  THEN
                   IF ( r(iray)%itarget < 0 )  CALL process_ray_canopy( iray, .FALSE. )
                ENDIF
                iray = iray - 1
             ENDDO
!
!--          Extend mixed rays by raising the lowest_free ray (remains unchanged if the previous
!--          loop exited immediately, becomes 0 if all rays are full/mixed)
             s%lowest_free_ray = iray
!
!--          Above free horizon, just process canopy for remaining rays
             IF ( plant_canopy )  THEN
                DO WHILE ( iray >= 1 )
                   CALL process_ray_canopy( iray, .FALSE. )
                   iray = iray - 1
                ENDDO
             ENDIF
          ELSE ! is_mixed_col
!
!--          Process canopy for remaining mixed and free rays
             IF ( plant_canopy)  THEN
                iray = s%lowest_mixed_ray
                DO WHILE ( iray > s%lowest_free_ray )
                   IF ( r(iray)%itarget >= 0 )  THEN
                      CALL process_ray_canopy( iray, .FALSE. )
                   ENDIF
                   iray = iray - 1
                ENDDO
                DO WHILE (iray >= 1 )
                   CALL process_ray_canopy( iray, .FALSE. )
                   iray = iray - 1
                ENDDO
             ENDIF
          ENDIF  ! is_mixed_col

          lastcolumn(:) = column(:)
          column(seldim) = column(seldim) + dimdelta(seldim)
!
!--       Save wall direction of coming building column (= this air column)
          IF ( seldim == 1 )  THEN
             IF ( dimdelta(seldim) == 1 )  THEN
                lastdir = isouth
             ELSE
                lastdir = inorth
             ENDIF
          ELSE
             IF ( dimdelta(seldim) == 1 )  THEN
                lastdir = iwest
             ELSE
                lastdir = ieast
             ENDIF
          ENDIF
          s%lastdist = nextdist
          s%dimnext(seldim) = s%dimnext(seldim) + dimdelta(seldim)
          s%dimnextdist(seldim) = ( s%dimnext(seldim) - 0.5_wp - yxorigin(seldim) ) /              &
                                  discr_azim_yxdir(seldim, s%iaz)
!
!--       If we have left whole domain, exit now (no vertical surfaces at lateral boundaries).
          IF ( column(seldim) < 0  .OR.  column(seldim) > glob_max(seldim) )  EXIT
!
!--       If all (potentially) free rays exited the urban layer upwards, stop
          IF ( s%lowest_mixed_ray <= 0 )  EXIT
          IF ( s%origin(1) + r(s%lowest_mixed_ray)%zdir * s%lastdist >=                            &
               REAL(nz_urban_t, wp) + 0.5_wp )  EXIT

          topo_full = topo_top_ind(column(1), column(2), 0)
          topo_mix  = topo_top_ind(column(1), column(2), 5)
          is_mixed_col = ( topo_mix /= topo_full )
          horz_full_entry  = ( REAL( topo_full, wp ) + 0.5_wp - s%origin(1) ) / s%lastdist
          horz_mixed_entry = ( REAL( topo_mix,  wp ) + 0.5_wp - s%origin(1) ) / s%lastdist
!
!--       Identify vertical full obstacles hit by rays in current column, mixed rays need to be
!--       checked if they are already obstructed.
          DO WHILE ( s%lowest_mixed_ray > s%lowest_free_ray )
             IF ( r(s%lowest_mixed_ray)%zdir > horz_full_entry )  EXIT
             IF ( r(s%lowest_mixed_ray)%itarget < 0 )  THEN
                 r(s%lowest_mixed_ray)%itarget = gridsurf(lastdir,                                 &
                                                          CEILING( -0.5_wp + s%origin(1) +         &
                                                                   r(s%lowest_mixed_ray)%zdir *    &
                                                                   s%lastdist ),                   &
                                                          lastcolumn(1), lastcolumn(2))
             ENDIF
             s%lowest_mixed_ray = s%lowest_mixed_ray - 1
          ENDDO
!
!--       Identify vertical full obstacles hit by rays in current column, free rays need no
!--       individual checks.
          DO WHILE ( s%lowest_mixed_ray > 0 )
             IF ( r(s%lowest_mixed_ray)%zdir > horz_full_entry )  EXIT
             r(s%lowest_mixed_ray)%itarget = gridsurf(lastdir,                                     &
                                                      CEILING( -0.5_wp + s%origin(1) +             &
                                                               r(s%lowest_mixed_ray)%zdir *        &
                                                               s%lastdist ),                       &
                                                      lastcolumn(1), lastcolumn(2))
             s%lowest_mixed_ray = s%lowest_mixed_ray - 1
          ENDDO
          IF ( s%lowest_free_ray > s%lowest_mixed_ray )  s%lowest_free_ray = s%lowest_mixed_ray
!
!--       Identify targets for vertical mixed obstacles.
!--       lowest_mixed_ray now points to bottom of vertical mixed obstacles.
          IF ( is_mixed_col )  THEN
!
!--          Previously mixed rays need to be checked whether they are obstructed
             DO  iray = s%lowest_mixed_ray, s%lowest_free_ray+1, -1
                IF ( r(iray)%zdir > horz_mixed_entry )  EXIT
                IF ( r(iray)%itarget >= 0 )  CYCLE
                r(iray)%itarget = gridsurf(lastdir, CEILING( -0.5_wp + s%origin(1) +               &
                                                               r(iray)%zdir * s%lastdist ),        &
                                           lastcolumn(1), lastcolumn(2))
                                     ! contains -999 if no surface at given d,z,y,x
             ENDDO
!
!--          Previously free rays cannot be obstructed yet
             iray = s%lowest_free_ray
             DO WHILE ( iray >= 1 )
                IF ( r(iray)%zdir > horz_mixed_entry )  EXIT

                r(iray)%itarget = gridsurf(lastdir, CEILING( -0.5_wp + s%origin(1) +               &
                                                               r(iray)%zdir * s%lastdist ),        &
                                           lastcolumn(1), lastcolumn(2))
                iray = iray - 1
             ENDDO
!
!--          Extend mixed rays by raising the lowest_free ray (remains unchanged if the previous
!--          loop exited immediately, becomes 0 if all rays are full/mixed)
             s%lowest_free_ray = iray
          ENDIF  ! End of mixed horizon

#if defined( __parallel )
          IF ( column(seldim) < loc_min(seldim)  .OR.  column(seldim) > loc_max(seldim) )  THEN
!
!--          Continue raytracing at next process
             CALL radiation_calc_global_offset( column(2), column(1), 0, 1, iproc = next_proc )

             CALL MPI_BSEND( MPI_BOTTOM, 1, dtype, next_proc, lrt_msg_tag, comm2d, ierr )
             IF ( ierr /= 0  .AND.  debug_output )  THEN
                WRITE( debug_string, * ) 'MPI Bsend Error 2:', ierr
                CALL debug_message( debug_string, 'info' )
             ENDIF
             RETURN
          ENDIF
#endif
       ENDDO
!
!--    Finished forward raytracing. Prepare to reverse if applicable.
       IF ( s%message == 'f'  .AND.  plant_canopy )  THEN
!
!--       Reset transparency for rays to the sky
          iray = s%lowest_mixed_ray
          DO WHILE ( iray > s%lowest_free_ray )
             IF ( r(iray)%itarget < 0 )  r(iray)%transp = 1.0_wp
             iray = iray - 1
          ENDDO
          DO WHILE ( iray > 0 )
             r(iray)%transp = 1.0_wp
             iray = iray - 1
          ENDDO
!
!--       Reverse indices and continue with backward direction (in the same process). We have
!--       stopped one column *after* the boundary (or last relevant column), so we can go one
!--       step back immediately.
          column(seldim) = column(seldim) - dimdelta(seldim)
          s%dimnext(seldim) = s%dimnext(seldim) - dimdelta(seldim)
          DO  d = 1, 2
             IF ( discr_azim_yxdir(d,s%iaz) /= 0.0_wp )  THEN
                s%dimnext(d) = s%dimnext(d) - dimdelta(d)
                s%dimnextdist(d) = ( s%dimnext(d) - 0.5_wp - yxorigin(d) ) /                       &
                                   discr_azim_yxdir(d,s%iaz)
             ENDIF
          ENDDO
          s%message = 'b'
       ELSE
!
!--       Return to start column and call finalize directly at origin process
!--       (no backward raytracing)
          s%message = CHAR( ICHAR( s%message ) - 32 ) ! convert to uppercase
#if defined( __parallel )
          DO  d = 1, 2
             IF ( discr_azim_yxdir(d,s%iaz) == 0.0_wp )  THEN
                column(d) = NINT(yxorigin(d))
             ELSE IF ( discr_azim_yxdir(d,s%iaz) > 0.0_wp )  THEN
                column(d) = FLOOR( yxorigin(d) + 0.5_wp )
             ELSE
                column(d) = CEILING( yxorigin(d) + 0.5_wp ) - 1
             ENDIF
          ENDDO
          CALL radiation_calc_global_offset( column(2), column(1), 0, 1, iproc = next_proc )
          CALL MPI_BSEND( MPI_BOTTOM, 1, dtype, next_proc, lrt_msg_tag, comm2d, ierr )
          IF ( ierr /= 0  .AND.  debug_output )  THEN
             WRITE( debug_string, * ) 'MPI Bsend Error 3:', ierr
             CALL debug_message( debug_string, 'info' )
          ENDIF
#else
          CALL raytrace_finalize( s, r )
#endif
          RETURN
       ENDIF
    ENDIF ! .NOT. backward
!
!-- Backward raytracing remote request - for plant canopy only
    DO
!
!--    Along what dimension will the next wall crossing be?
       seldim = MAXLOC( s%dimnextdist, 1 )
       nextdist = s%dimnextdist(seldim)
       IF ( nextdist <= 0.0_wp )  THEN
          nextdist = 0.0_wp
       ENDIF

       dxxyy = ( ( dy * discr_azim_yxdir(1,s%iaz) )**2 +                                           &
                 ( dx * discr_azim_yxdir(2,s%iaz) )**2 ) * ( nextdist - s%lastdist )**2
       topo_full = topo_top_ind(column(1), column(2), 0)

       DO  iray = 1, s%nrays
!
!--       We only process rays from the sky, so they cannot start in the middle of a cell.
          IF ( r(iray)%itarget < 0 )  CALL process_ray_canopy( iray, .TRUE. )
       ENDDO

       IF ( nextdist <= 0.0_wp )  EXIT

       column(seldim) = column(seldim) - dimdelta(seldim)
       s%lastdist = nextdist
       s%dimnext(seldim) = s%dimnext(seldim) - dimdelta(seldim)
       s%dimnextdist(seldim) = ( s%dimnext(seldim) - 0.5_wp - yxorigin(seldim) ) /                 &
                               discr_azim_yxdir(seldim,s%iaz)

#if defined( __parallel )
       IF ( column(seldim) < loc_min(seldim)  .OR.  column(seldim) > loc_max(seldim) )  THEN
!
!--       Continue raytracing at next process
          CALL radiation_calc_global_offset( column(2), column(1), 0, 1, iproc = next_proc )

          CALL MPI_BSEND( MPI_BOTTOM, 1, dtype, next_proc, lrt_msg_tag, comm2d, ierr )
          IF ( ierr /= 0  .AND.  debug_output )  THEN
             WRITE( debug_string, * ) 'MPI Bsend Error 2:', ierr
             CALL debug_message( debug_string, 'info' )
          ENDIF
          RETURN
       ENDIF
#else
       IF ( dtype == 0 )  CONTINUE ! silence compiler waring about unused parameter
#endif
    ENDDO
!
!-- Save raytracing results
    s%message = 'F'
    CALL raytrace_finalize( s, r )

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Traces a section of ray through one column with potential plant canopy and creates canopy view
!> factors where applicable.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE process_ray_canopy( c_ray, backward, zstop )

       INTEGER(iwp), INTENT(IN)       ::  c_ray     !< ray index
       INTEGER(iwp)                   ::  ngrow     !< new size of acsf

       LOGICAL, INTENT(IN)            ::  backward  !< backward direction sky->face

       REAL(wp)                       ::  rtopo     !< local topo height
       REAL(wp), INTENT(IN), OPTIONAL ::  zstop     !< z-coordinate of the obstacle for a ray that is
                                                    !< partially blocked in the column

       zentry = s%origin(1) + r(c_ray)%zdir * s%lastdist
       zexit = s%origin(1) + r(c_ray)%zdir * nextdist
       ztop = REAL( pct(column(1), column(2)), wp ) + 0.5_wp
       zsgn = INT( SIGN( 1.0_wp, r(c_ray)%zdir ), iwp )

       IF ( r(c_ray)%zdir == 0.0_wp )  THEN  ! Ray is exactly horizontal, always full-length
          IF ( zentry > ztop )  RETURN
          nz = 2
          rt2_dist(nz) = SQRT( dxxyy )
          kz = CEILING( -0.5_wp + zentry, iwp )
       ELSE
          IF ( backward )  THEN
             IF ( zexit >= ztop )  RETURN ! zexit is always closer to origin
             zsgn = -zsgn
             IF ( PRESENT( zstop ) )  zentry = zstop
          ELSE
             IF ( zentry >= ztop )  RETURN ! zentry is always closer to origin
             IF ( PRESENT( zstop ) )  zexit = zstop
          ENDIF
!
!--       Cut by PC top on top and by terrain on bottom (bottom should not be necessary, but
!--       rounding errors sometimes place it below).
          rtopo = REAL( topo_full, wp ) + 0.5_wp
          zentry = MAX( rtopo, MIN( ztop, zentry ) )
          zexit =  MAX( rtopo, MIN( ztop, zexit  ) )

          zb0 = FLOOR( zentry * zsgn - 0.5_wp ) + 1  ! Because it must be greater than orig
          zb1 = CEILING( zexit * zsgn - 0.5_wp ) - 1  ! Because it must be smaller than exit
          nz = MAX( zb1 - zb0 + 3, 2 )
          rt2_dist(nz) = SQRT( ( ( zexit - zentry ) * dz(1) )**2 + dxxyy )
          qdist = rt2_dist(nz) / ( zexit - zentry )
          rt2_dist(2:nz-1) = (/ ( ( ( REAL( l, wp ) + 0.5_wp ) * zsgn - zentry ) * qdist,          &
                                l = zb0, zb1 ) /)
          kz = zb0 * zsgn
       ENDIF

       IF ( s%message == 'p'  .AND.  s%lastdist == 0.0_wp )  THEN
!
!--       Skip first plant canopy box
          l = 3
          kz = kz + zsgn
       ELSE
          l = 2
       ENDIF

       DO WHILE ( l <= nz )
          IF ( lad_s(kz-topo_full, column(1), column(2)) > 0.0_wp )  THEN
             curtrans = EXP( - ext_coef * lad_s(kz-topo_full, column(1), column(2)) *              &
                               ( rt2_dist(l) - rt2_dist(l-1) ) )

             IF ( create_csf )  THEN
                IF ( ncsfl + 1 > ncsfla )  THEN
                   ngrow = CEILING( REAL( ncsfl + 1, KIND = wp ) * grow_factor )
                   CALL merge_and_grow_csf(ngrow)
                ENDIF
                ncsfl = ncsfl + 1
                acsf(ncsfl)%ip = myid ! purely for usage with legacy-compatible quicksort_csf
                acsf(ncsfl)%itx = column(2)
                acsf(ncsfl)%ity = column(1)
                acsf(ncsfl)%itz = kz
                IF ( backward )  THEN
!
!--                Backward from sky (rays from faces are skipped in raytrace)
                   acsf(ncsfl)%isurfs = -1
                   acsf(ncsfl)%rcvf = ( 1.0_wp - curtrans ) * r(c_ray)%transp * s%aorig *          &
                                      r(c_ray)%vffrac
                ELSE
                   acsf(ncsfl)%isurfs = s%iorig
                   acsf(ncsfl)%rcvf = ( 1.0_wp - curtrans ) * r(c_ray)%transp * r(c_ray)%vffrac
                ENDIF
             ENDIF

             r(c_ray)%transp = r(c_ray)%transp * curtrans
          ENDIF
          l = l + 1
          kz = kz + zsgn
       ENDDO
    END SUBROUTINE process_ray_canopy

 END SUBROUTINE raytrace_segment


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks completion of local and previous process in localized raytracing and passes on the message
!> if applicable.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE lrt_check_completion( s, dtype )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  dtype  !< MPI datatype associated with the message
    INTEGER(iwp)             ::  i

    TYPE(t_traced_section), INTENT(INOUT) ::  s  !< raytraced section message structure


    IF ( lrt_prev_process_complete  .AND.  lrt_local_complete )  THEN
       IF ( myid == numprocs - 1 )  THEN
!
!--       Send stop info to all processes, incl. self to consume the last open MPI_REQUEST
          DO  i = 0, numprocs - 1
             s%message = 't'
             CALL MPI_BSEND( MPI_BOTTOM, 1, dtype, i, lrt_msg_tag, comm2d, ierr )
             IF ( ierr /= 0  .AND.  debug_output )  THEN
                WRITE( debug_string, * ) 'MPI Bsend Error 3:', ierr
                CALL debug_message( debug_string, 'info' )
             ENDIF
          ENDDO
       ELSE
!
!--       Send complete command to the next process
          s%message = 'c'
          CALL MPI_BSEND( MPI_BOTTOM, 1, dtype, myid+1, lrt_msg_tag, comm2d, ierr )
          IF ( ierr /= 0  .AND.  debug_output )  THEN
             WRITE( debug_string, * ) 'MPI Bsend Error 4:', ierr
             CALL debug_message( debug_string, 'info' )
          ENDIF
       ENDIF
    ENDIF

 END SUBROUTINE lrt_check_completion
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Processes all pending (queued) messages in localized raytracing.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE lrt_process_pending( wait_for )

    IMPLICIT NONE

    INTEGER(iwp) ::  clock1  !< current system clock
    INTEGER(iwp) ::  clock2  !< current system clock
    INTEGER(iwp) ::  dtype   !< MPI datatype associated with the message

    INTEGER(iwp), DIMENSION(MPI_STATUS_SIZE) ::  wait_status  !<

    LOGICAL             ::  ready
    LOGICAL, INTENT(IN) ::  wait_for  !< wait for at least one message

    TYPE(t_traced_section), POINTER ::  section  !< raytraced section message structure

    TYPE(t_traced_ray), DIMENSION(:), POINTER ::  rays  !< array of ray message structures


    ready = .NOT. wait_for
    DO
       IF ( ready )  THEN
          CALL MPI_TEST(lrt_req_incoming, ready, wait_status, ierr)
          IF ( ierr /= 0  .AND.  debug_output )  THEN
             WRITE( debug_string, * ) 'MPI Error MPI_WAIT:', ierr
             CALL debug_message( debug_string, 'info' )
          ENDIF
          IF ( .NOT. ready )  EXIT ! no pending messages
       ELSE
!
!--       We cannot use cpu_log, because that has too much management code and we are measuring
!--       very small times very often
          CALL SYSTEM_CLOCK( COUNT=clock1 )
          CALL MPI_WAIT(lrt_req_incoming, wait_status, ierr)
          CALL SYSTEM_CLOCK( COUNT=clock2 )
          lrt_wait_time = lrt_wait_time + ( clock2 - clock1 )
          IF ( ierr /= 0  .AND.  debug_output )  THEN
             WRITE( debug_string, * ) 'MPI Error MPI_WAIT:', ierr
             CALL debug_message( debug_string, 'info' )
          ENDIF
          ready = .TRUE.
       ENDIF

       IF ( lrt_msg_incoming_section%message == 't' )  THEN
!
!--       Termination (all procs finished)
          lrt_unfinished = .FALSE.
          EXIT
       ENDIF
!
!--    Swap incoming and processing messages and listen again
       section => lrt_msg_incoming_section
       rays => lrt_msg_incoming_rays
       dtype = lrt_msg_type_incoming

       lrt_msg_incoming_section => lrt_msg_processing_section
       lrt_msg_incoming_rays => lrt_msg_processing_rays
       lrt_msg_type_incoming = lrt_msg_type_processing

       lrt_msg_processing_section => section
       lrt_msg_processing_rays => rays
       lrt_msg_type_processing = dtype

       CALL MPI_IRECV( MPI_BOTTOM, 1, lrt_msg_type_incoming, MPI_ANY_SOURCE, lrt_msg_tag, comm2d,  &
                       lrt_req_incoming, ierr )
       IF ( ierr /= 0  .AND.  debug_output )  THEN
          WRITE( debug_string, * ) 'MPI Irecv Error:', ierr
          CALL debug_message( debug_string, 'info' )
       ENDIF
!
!--    Process incoming message
       SELECT CASE ( section%message )
          CASE ( 'f', 'b', 'p', 'm' )
!
!--          Raytrace forward, backward
             CALL raytrace_segment( section, rays, dtype )
          CASE ( 'F', 'P', 'M' )
!
!--          Direct finalization of raytracing
             CALL raytrace_finalize( section, rays )
          CASE ( 'c' )
!
!--          Complete raytracing from prev process
             lrt_prev_process_complete = .TRUE.
             CALL lrt_check_completion( section, dtype )
       END SELECT
    ENDDO

 END SUBROUTINE lrt_process_pending
#endif


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocates memory for incoming message in localized raytracing and initiates the respective MPI
!> datatypes.
!> NOTE: We have tried using MPI_ALLOC_MEM/MPI_FREE_MEM together with C_F_POINTER, which should use
!> a more optimized allocation where available (and fall back to standard allocation), but for some
!> reason it does not work at all on Intel compilers, therefore we allocate to standard Fortran
!> pointers.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE lrt_allocate_message( section, rays, dtype, tot_size )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(OUT)           ::  dtype     !< MPI derived datatype for the allocated message
    INTEGER(iwp), INTENT(OUT), OPTIONAL ::  tot_size  !< Total size of message in bytes

    INTEGER(iwp), DIMENSION(2) ::  lengths  !< buffer lengths

    TYPE(t_traced_section), POINTER, INTENT(OUT) ::  section  !<

    TYPE(t_traced_ray), DIMENSION(:), POINTER, INTENT(OUT) ::  rays  !<

#if defined( __parallel )
    INTEGER(MPI_ADDRESS_KIND), DIMENSION(2) ::  displ  !< buffer displacements


!
!-- Allocate section structure
    lengths(1) = STORAGE_SIZE(section) / 8
#endif
    ALLOCATE( section )
#if defined( __parallel )
    CALL MPI_GET_ADDRESS( section, displ(1), ierr )
!
!-- Allocate array of structures for rays
    lengths(2) = ( STORAGE_SIZE(rays) / 8 ) * raytrace_discrete_elevs
#endif
    ALLOCATE( rays(raytrace_discrete_elevs) )
#if defined( __parallel )
    CALL MPI_GET_ADDRESS( rays, displ(2), ierr )
!
!-- Initialize MPI datatype
    CALL MPI_TYPE_CREATE_HINDEXED( 2, lengths, displ, MPI_BYTE, dtype, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
        WRITE( debug_string, * ) 'Error MPI_TYPE_CREATE_HINDEXED:', ierr
        CALL debug_message( debug_string, 'info' )
    ENDIF
    CALL MPI_TYPE_COMMIT( dtype, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
        WRITE( debug_string, * ) 'Error MPI_TYPE_COMMIT:', ierr
        CALL debug_message( debug_string, 'info' )
    ENDIF
#else
    IF ( dtype == 0 )  CONTINUE ! silence compiler waring about unused parameter
#endif

    IF ( PRESENT( tot_size ) )  THEN
       tot_size = SUM( lengths )
    ENDIF

 END SUBROUTINE lrt_allocate_message


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Frees memory and MPI datatype associated with a localized raytracing message.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE lrt_free_message( section, rays, dtype )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(INOUT) ::  dtype  !<

    TYPE(t_traced_section), POINTER, INTENT(INOUT) ::  section  !<

    TYPE(t_traced_ray), DIMENSION(:), POINTER, INTENT(INOUT) ::  rays  !<


#if defined( __parallel )
    CALL MPI_TYPE_FREE( dtype, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
        WRITE ( debug_string, * ) 'Error MPI_TYPE_FREE:', ierr
        CALL debug_message( debug_string, 'info' )
    ENDIF
#else
    IF ( dtype == 0 )  CONTINUE ! silence compiler waring about unused parameter
#endif

    DEALLOCATE( rays )
    DEALLOCATE( section )

 END SUBROUTINE lrt_free_message


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Performs simplified 2-D raytracing in 2.5-D geometry with terrain/buildings and opaque part of
!> plant canopy. Returns horizon heights for each vertical level at (y,x) coordinates of origin.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE trace_horizons( iorig, jorig, yxdir, horizon )

    IMPLICIT NONE

    INTEGER(iwp) ::  d         !< index
    INTEGER(iwp) ::  k         !< index
    INTEGER(iwp) ::  k_bottom  !< lowest meaningful level at origin
    INTEGER(iwp) ::  ig        !< 1D index of grid column in global 2D array
    INTEGER(iwp) ::  seldim    !< dimension to be incremented

    INTEGER(iwp), INTENT(IN) ::  iorig  !< coordinates of origin
    INTEGER(iwp), INTENT(IN) ::  jorig  !< coordinates of origin

    INTEGER(iwp), DIMENSION(2) ::  column    !< grid column being crossed
    INTEGER(iwp), DIMENSION(2) ::  dimnext   !< next dimension increments along path
    INTEGER(iwp), DIMENSION(2) ::  dimdelta  !< dimension direction = +- 1

    REAL(wp) ::  lastdist  !< beginning of current crossing
    REAL(wp) ::  nextdist  !< end of current crossing

    REAL(wp), DIMENSION(2) ::  ddim         !<
    REAL(wp), DIMENSION(2) ::  glob_max     !< global domain boundary (north y, east x)
    REAL(wp), DIMENSION(2) ::  dimnextdist  !< distance for each dimension increments
    REAL(wp), DIMENSION(2) ::  yxorigin     !< coordinates of origin (y,x)
    REAL(wp), DIMENSION(2) ::  yxdir_d      !< yxdir / ddim

    REAL(wp), DIMENSION(2), INTENT(IN) ::  yxdir  !< y,x *unit* vector of ray direction (in physical coords)

    REAL(wp), DIMENSION(nz_urban_b:nz_urban_t), INTENT(OUT) ::  horizon  !< tan(horizon) in physical coords for each level
                                                                         !< at yxorigin

    REAL(wp), PARAMETER ::  eps = 1E-10_wp  !< epsilon for value comparison


    ddim = (/ dy, dx /)
    yxdir_d(:) = yxdir(:) / ddim(:)
    yxorigin = (/ REAL( jorig, wp ), REAL( iorig, wp ) /)
    k_bottom = opaque_top_l(jorig, iorig) + 1
    horizon(:) = -HUGE( 1.0_wp )
    glob_max(1) = ny
    glob_max(2) = nx
    lastdist = 0.0_wp

!-- Since all face coordinates have values *.5 and we'd like to use integers, all these have
!-- 0.5 added
    DO  d = 1, 2
       IF ( yxdir(d) == 0.0_wp )  THEN
          dimnext(d) = HUGE( 1 )
          dimdelta(d) = HUGE( 1 )
          dimnextdist(d) = HUGE( 1.0_wp )
          column(d) = NINT(yxorigin(d))
       ELSEIF ( yxdir(d) > 0.0_wp )  THEN
          dimnext(d) = FLOOR( yxorigin(d) + 0.5_wp ) + 1
          dimdelta(d) = 1
          dimnextdist(d) = ( dimnext(d) - 0.5_wp - yxorigin(d) ) / yxdir_d(d)
          column(d) = dimnext(d) - 1
       ELSE
          dimnext(d) = CEILING( yxorigin(d) + 0.5_wp ) - 1
          dimdelta(d) = -1
          dimnextdist(d) = ( dimnext(d) - 0.5_wp - yxorigin(d) ) / yxdir_d(d)
          column(d) = dimnext(d)
       ENDIF
    ENDDO

    DO
!
!--    Along what dimension will the next wall crossing be?
       seldim = MINLOC( dimnextdist, 1 )
       nextdist = dimnextdist(seldim)

       IF ( nextdist - lastdist >= eps )  THEN
!
!--       Calculate index of the grid with global indices (column(1),column(2)) in the array
!--       opaque_top
          CALL radiation_calc_global_offset( column(2), column(1), 0, 1, offs_glob = ig )

          IF ( lastdist == 0.0_wp )  THEN
             DO k = k_bottom, nz_urban_t
                horizon(k) = MAX( horizon(k),                                             &
                                  ( REAL( opaque_top(ig) - k, wp ) + 0.5_wp ) * dz(1) / nextdist )   ! exit
             ENDDO
          ELSE
             DO k = k_bottom, nz_urban_t
                horizon(k) = MAX( horizon(k),                                             &
                                  ( REAL( opaque_top(ig) - k, wp ) + 0.5_wp ) * dz(1) / lastdist,  & ! entry
                                  ( REAL( opaque_top(ig) - k, wp ) + 0.5_wp ) * dz(1) / nextdist )   ! exit
             ENDDO
          ENDIF

       ENDIF

       column(seldim) = column(seldim) + dimdelta(seldim)
       IF ( column(seldim) < 0  .OR.  column(seldim) > glob_max(seldim) )  EXIT

       lastdist = nextdist
       dimnext(seldim) = dimnext(seldim) + dimdelta(seldim)
       dimnextdist(seldim) = ( dimnext(seldim) - 0.5_wp - yxorigin(seldim) ) / yxdir_d(seldim)
    ENDDO

 END SUBROUTINE trace_horizons


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Calculates apparent solar positions for all timesteps and stores discretized positions for RTM.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_presimulate_solar_pos

    USE control_parameters,                                                                        &
        ONLY:  rotation_angle

    IMPLICIT NONE

    INTEGER(iwp) ::  it, i, j  !< loop indices

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  dsidir_tmp  !< dsidir_tmp[:,i] = unit vector of i-th
                                                          !< appreant solar direction

    ALLOCATE( dsidir_rev(0:raytrace_discrete_elevs/2-1,0:raytrace_discrete_azims-1) )
    dsidir_rev(:,:) = -1
    ALLOCATE( dsidir_tmp(3, raytrace_discrete_elevs/2*raytrace_discrete_azims) )
    ndsidir = 0
    sun_direction = .TRUE.

!
!-- Process spinup time if configured
    IF ( spinup_time > 0.0_wp )  THEN
       DO  it = 0, CEILING( spinup_time / dt_spinup )
          CALL simulate_pos( it * dt_spinup - spinup_time )
       ENDDO
    ENDIF
!
!-- Process simulation time
    DO  it = 0, CEILING( ( end_time - spinup_time ) / dt_radiation )
       CALL simulate_pos( it * dt_radiation )
    ENDDO
!
!-- Allocate global vars which depend on ndsidir
    ALLOCATE( dsidir ( 3, ndsidir ) )
    dsidir(:,:) = dsidir_tmp(:, 1:ndsidir)
    DEALLOCATE( dsidir_tmp )

    ALLOCATE( dsitrans(nsurfl, ndsidir) )
    ALLOCATE( dsitransc(npcbl, ndsidir) )
    IF ( nmrtbl > 0 )  ALLOCATE( mrtdsit(nmrtbl, ndsidir) )

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Simuates a single position
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE simulate_pos( time_since_reference_local )

    REAL(wp) ::  solar_azim  !< solar azimuth in rotated model coordinates

    REAL(wp), INTENT(IN) ::  time_since_reference_local  !< local time since reference

!
!-- Update apparent solar position based on modified t_s_r_p
    CALL get_date_time( time_since_reference_local, day_of_year = day_of_year,                     &
                        second_of_day = second_of_day )
    CALL calc_zenith( day_of_year, second_of_day )
    IF ( cos_zenith > 0 )  THEN
!
!--    Identify solar direction vector (discretized number) 1)
       solar_azim = ATAN2( sun_dir_lon, sun_dir_lat ) * ( 180.0_wp / pi ) - rotation_angle
       i = MODULO( NINT( solar_azim / 360.0_wp * REAL( raytrace_discrete_azims, KIND = wp )        &
                         - 0.5_wp, iwp ), raytrace_discrete_azims )
       j = FLOOR( ACOS( cos_zenith ) / pi * REAL( raytrace_discrete_elevs, KIND = wp ) )
       IF ( dsidir_rev(j, i) == -1 )  THEN
          ndsidir = ndsidir + 1
          dsidir_tmp(:, ndsidir) =                                                                 &
             (/ COS( (REAL( j, wp ) + 0.5_wp ) * pi        / REAL( raytrace_discrete_elevs, wp ) ),&
                SIN( (REAL( j, wp ) + 0.5_wp ) * pi        / REAL( raytrace_discrete_elevs, wp ) ) &
              * COS( (REAL( i, wp ) + 0.5_wp ) * 2.0_wp*pi / REAL( raytrace_discrete_azims, wp ) ),&
                SIN( (REAL( j, wp ) + 0.5_wp ) * pi        / REAL( raytrace_discrete_elevs, wp ) ) &
              * SIN( (REAL( i, wp ) + 0.5_wp ) * 2.0_wp*pi / REAL( raytrace_discrete_azims, wp ) ) /)
          dsidir_rev(j, i) = ndsidir
       ENDIF
    ENDIF
 END SUBROUTINE simulate_pos

 END SUBROUTINE radiation_presimulate_solar_pos


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Determines whether two faces are oriented towards each other in RTM. Since the surfaces follow
!> the gird box surfaces, it checks first whether the two surfaces are directed in the same
!> direction, then it checks if the two surfaces are located in confronted direction but facing away
!> from each other, e.g. <--| |-->
!--------------------------------------------------------------------------------------------------!
 PURE LOGICAL FUNCTION surface_facing( x, y, z, d, x2, y2, z2, d2 )

    IMPLICIT NONE

    INTEGER(iwp),INTENT(IN) ::  x, y, z, d, x2, y2, z2, d2  !<

    surface_facing = .FALSE.
!
!-- First check: are the two surfaces directed in the same direction
    IF ( d == iup     .AND.  d2 == iup    )  RETURN
    IF ( d == isouth  .AND.  d2 == isouth )  RETURN
    IF ( d == inorth  .AND.  d2 == inorth )  RETURN
    IF ( d == iwest   .AND.  d2 == iwest  )  RETURN
    IF ( d == ieast   .AND.  d2 == ieast  )  RETURN
!
!-- Second check: are surfaces facing away from each other
    SELECT CASE (d)
        CASE (iup)                   !< Upward facing surfaces
            IF ( z2 < z )  RETURN
        CASE (isouth)                !< Southward facing surfaces
            IF ( y2 > y )  RETURN
        CASE (inorth)                !< Northward facing surfaces
            IF ( y2 < y )  RETURN
        CASE (iwest)                 !< Westward facing surfaces
            IF ( x2 > x )  RETURN
        CASE (ieast)                 !< Eastward facing surfaces
            IF ( x2 < x )  RETURN
    END SELECT

    SELECT CASE (d2)
        CASE (iup)                   !< Ground, roof
            IF ( z < z2 )  RETURN
        CASE (isouth)                !< South facing
            IF ( y > y2 )  RETURN
        CASE (inorth)                !< North facing
            IF ( y < y2 )  RETURN
        CASE (iwest)                 !< West facing
            IF ( x > x2 )  RETURN
        CASE (ieast)                 !< East facing
            IF ( x < x2 )  RETURN
        CASE (-1)
            CONTINUE
    END SELECT

    surface_facing = .TRUE.

 END FUNCTION surface_facing


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads svf, svfsurf, csf, csfsurf and mrt factors data from saved file. This allows to skip their
!> calculation during of RTM init phase. SVF means sky view factors and CSF means canopy sink
!> factors.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_read_svf

    IMPLICIT NONE

    CHARACTER(rad_version_len) ::  rad_version_field  !<

    CHARACTER(LEN=64) ::  rad_version_in  !< rad version on file

    INTEGER(iwp) ::  i                      !<

    INTEGER(iwp) ::  ndsidir_bin_file = 0   !<
    INTEGER(iwp) ::  nmrtbl_bin_file  = 0   !<
    INTEGER(iwp) ::  npcbl_bin_file   = 0   !<
    INTEGER(iwp) ::  nsurfl_bin_file  = 0   !<

    INTEGER(idp) ::  ncsfl_tot              !<
    INTEGER(idp) ::  ndsidir_from_file = 0  !<
    INTEGER(idp) ::  ndsidir_tot            !<
    INTEGER(idp) ::  nmrtbl_from_file  = 0  !<
    INTEGER(idp) ::  nmrtbl_tot             !<
    INTEGER(idp) ::  nmrtf_tot              !<
    INTEGER(idp) ::  npcbl_from_file   = 0  !<
    INTEGER(idp) ::  npcbl_tot              !<
    INTEGER(idp) ::  nsurfl_from_file  = 0  !<
    INTEGER(idp) ::  nsurfl_tot             !<
    INTEGER(idp) ::  nsvfl_tot              !< total (sum across all PEs) counter for the different IO variables
#if defined( __parallel )
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  tmp_int !< temporary array to store 0d dimension sizes
#endif

    INTEGER(idp), DIMENSION(4) ::  global_sum  !<
    INTEGER(idp), DIMENSION(4) ::  local_sum   !< variables to compute total counter

    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_end    !< global end index     (I8)
    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_start  !< global start index   (I8)
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  end_index     !< local end index
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  start_index   !< local start index

    LOGICAL ::  data_to_read  !< flag indicating if data is available for current variable

#if defined( __parallel )
    REAL(wp), ALLOCATABLE, DIMENSION(:)     ::  tmp_1d    !< dummy array to read 1d arrays
#endif
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::  tmp_2d    !< dummy array to read 2d arrays
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  tmp_3d  !< dummy array used for reordering of 3d arrays


    IF ( TRIM( restart_data_format_input ) == 'fortran_binary' )  THEN

       CALL location_message( 'reading sky view factors in Fortran binary format', 'start' )

       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN

!
!--          numprocs_previous_run is only known in case of reading restart data. If a new initial
!--          run which reads svf data is started the following query will be skipped.
             IF ( initializing_actions == 'read_restart_data' )  THEN

                IF ( numprocs_previous_run /= numprocs )  THEN
                   WRITE( message_string, * ) 'A different number of processors between the run ', &
                                              'that has written the svf data and the one that ',   &
                                              'will read it is not allowed'
                   CALL message( 'radiation_read_svf', 'RAD0053', 1, 2, 0, 6, 0 )
                ENDIF

             ENDIF

!
!--          Open binary file
             CALL check_open( 88 )

!
!--          Read and check version
             READ ( 88 ) rad_version_field
             IF ( TRIM( rad_version_field ) /= TRIM( rad_version ) )  THEN
                 WRITE( message_string, * ) 'Version of binary SVF file "',                        &
                                            TRIM( rad_version_field ), '" does not match ',        &
                                            'the version of model "', TRIM( rad_version ), '"'
                 CALL message( 'radiation_read_svf', 'RAD0054', 1, 2, 0, 6, 0 )
             ENDIF

!
!--          Read nsvfl, ncsfl, nsurfl, nmrtf
             READ ( 88 ) nsvfl, ncsfl, nsurfl_bin_file, npcbl_bin_file, ndsidir_bin_file,       &
                         nmrtbl_bin_file, nmrtf

             IF ( nsvfl < 0  .OR.  ncsfl < 0 )  THEN
                 WRITE( message_string, * ) 'wrong number of SVF or CSF'
                 CALL message( 'radiation_read_svf', 'RAD0055', 1, 2, 0, 6, 0 )
             ELSE
                 WRITE( debug_string , * ) 'Number of SVF, CSF, and nsurfl to read', nsvfl, ncsfl,     &
                                       nsurfl_bin_file
                 IF ( debug_output )  CALL debug_message( debug_string, 'info' )
             ENDIF

             IF ( nsurfl_bin_file /= nsurfl )  THEN
                 WRITE( message_string, * ) 'nsurfl from SVF file does not match calculated ',     &
                                            'nsurfl from radiation_interaction_init'
                 CALL message( 'radiation_read_svf', 'RAD0056', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( npcbl_bin_file /= npcbl )  THEN
                 WRITE( message_string, * ) 'npcbl from SVF file does not match calculated npcbl', &
                                            ' from radiation_interaction_init'
                 CALL message( 'radiation_read_svf', 'RAD0057', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( ndsidir_bin_file /= ndsidir )  THEN
                 WRITE( message_string, * ) 'ndsidir from SVF file does not match calculated ',    &
                                            'ndsidir from radiation_presimulate_solar_pos'
                 CALL message( 'radiation_read_svf', 'RAD0058', 1, 2, 0, 6, 0 )
             ENDIF

             IF ( nmrtbl_bin_file /= nmrtbl )  THEN
                 WRITE( message_string, * ) 'nmrtbl from SVF file does not match calculated ',     &
                                            'nmrtbl from radiation_interaction_init'
                 CALL message( 'radiation_read_svf', 'RAD0059', 1, 2, 0, 6, 0 )
             ELSE
                 WRITE( debug_string, * ) 'Number of nmrtf to read ', nmrtf
                 IF ( debug_output )  CALL debug_message( debug_string, 'info' )
             ENDIF

!
!--          Arrays skyvf, skyvft, dsitrans and dsitransc are allready allocated in
!--          radiation_interaction_init and radiation_presimulate_solar_pos
             IF ( nsurfl > 0 )  THEN
                READ( 88 ) skyvf
                READ( 88 ) skyvft
                READ( 88 ) dsitrans
             ENDIF

             IF ( plant_canopy  .AND.  npcbl > 0 )  THEN
                READ( 88 )  dsitransc
             ENDIF

!
!--          The allocation of svf, svfsurf, csf, csfsurf, mrtf, mrtft, and mrtfsurf happens in
!--          routine radiation_calc_svf which is not called if the program enters
!--          radiation_read_svf. Therefore these arrays have to be allocated in the following.
             IF ( nsvfl > 0 )  THEN
                ALLOCATE( svf(ndsvf,nsvfl) )
                ALLOCATE( svfsurf(idsvf,nsvfl) )
                READ( 88 ) svf
                READ( 88 ) svfsurf
             ENDIF

             IF ( plant_canopy  .AND.  ncsfl > 0 )  THEN
                ALLOCATE( csf(ndcsf,ncsfl) )
                ALLOCATE( csfsurf(idcsf,ncsfl) )
                READ( 88 ) csf
                READ( 88 ) csfsurf
             ENDIF

             IF ( nmrtbl > 0 )  THEN
                READ( 88 ) mrtsky
                READ( 88 ) mrtskyt
                READ( 88 ) mrtdsit
             ENDIF

             IF ( nmrtf > 0 )  THEN
                ALLOCATE( mrtf(nmrtf) )
                ALLOCATE( mrtft(nmrtf) )
                ALLOCATE( mrtfsurf(2,nmrtf) )
                READ( 88 ) mrtf
                READ( 88 ) mrtft
                READ( 88 ) mrtfsurf
             ENDIF

#if defined( __parallel )
             READ( 88 ) nsend_radx, nrecv_radx, niters_radx, nmaxsend_radx
             ALLOCATE( isurf_send_radx(nsend_radx) )
             ALLOCATE( isurf_recv_radx(nrecv_radx) )
             ALLOCATE( disp_send_radx(0:numprocs) )
             ALLOCATE( disp_recv_radx(0:numprocs) )
             ALLOCATE( radx_send(nsend_radx) )
             ALLOCATE( surfoutl_recv(nrecv_radx) )
             ALLOCATE( surfouts_recv(nrecv_radx) )
             ALLOCATE( radx_send_surfinl(nrecv_radx) )
             ALLOCATE( surfinl_recv(nsend_radx) )
             ALLOCATE( disp_sendbuf_radx(numprocs) )
             ALLOCATE( disp_recvbuf_radx(numprocs) )
             ALLOCATE( num_sendbuf_radx(numprocs) )
             ALLOCATE( num_recvbuf_radx(numprocs) )
             READ( 88 ) isurf_send_radx
             READ( 88 ) isurf_recv_radx
             READ( 88 ) disp_send_radx
             READ( 88 ) disp_recv_radx
#endif

             IF ( radiation_volumetric_flux ) THEN
!
!--              skyvf_vol is already allocated in radiation_interaction_init
                 ALLOCATE( shadow_top(nys:nyn,nxl:nxr,ndsidir) )
                 READ( 88 ) shadow_top
                 READ( 88 ) skyvf_vol
             ENDIF
!
!--          Close binary file
             CALL close_file( 88 )

          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

       CALL location_message( 'reading sky view factors in Fortran binary format', 'finished' )


    ELSEIF ( restart_data_format_input(1:3) == 'mpi' )  THEN

!
!--    In case of MPI-IO data is treated like surface data and the respective routines from
!--    restart_data_mpi_io_mod are used for reading. Contrary to restart data, the virtual processor
!--    grid can not be changed between writing and reading svf data.
       CALL location_message( 'reading sky view factors in MPI-IO format', 'start' )
!
!--    Open MPI-IO SVF file for reading the global data.
       CALL rd_mpi_io_open( 'READ', 'SVFIN' // TRIM( coupling_char ),                              &
                            open_for_global_io_only = .TRUE. )
!
!--    Check general header.
       IF ( tgh%pes_along_x /= npex  .OR.  tgh%pes_along_y /= npey )  THEN
!
!--       Force re-calculation of svfs.
          read_svf = .FALSE.
          WRITE( message_string, '(A,I7,A,I7,A,I7,A,I7,A)' )                                       &
              'virtual PE grid has changed between previous and current run &npex_prev = ',        &
              tgh%pes_along_x, ' npey_prev = ', tgh%pes_along_y, ' npex_new = ', npex,             &
              ' npey_new = ', npey, '&svf will be re-calculated'
          CALL message( 'radiation_read_svf', 'RAD0060', 0, 0, 0, 6, 0 )
          RETURN
       ENDIF
!
!--    Read global variables
       CALL rrd_mpi_io( 'rad_version', rad_version_in )
       CALL rrd_mpi_io( 'nsvfl', nsvfl_tot )
       CALL rrd_mpi_io( 'ncsfl', ncsfl_tot )
       CALL rrd_mpi_io( 'nsurfl', nsurfl_from_file )
       CALL rrd_mpi_io( 'npcbl', npcbl_from_file )
       CALL rrd_mpi_io( 'ndsidir', ndsidir_from_file )
       CALL rrd_mpi_io( 'nmrtbl', nmrtbl_from_file )
       CALL rrd_mpi_io( 'nmrtf', nmrtf_tot )

       CALL rd_mpi_io_close
!
!--    Compute global values of local counters.
       local_sum(1) = nsurfl
       local_sum(2) = npcbl
       local_sum(3) = ndsidir
       local_sum(4) = nmrtbl
#if defined( __parallel )
       CALL MPI_ALLREDUCE( local_sum, global_sum, SIZE(local_sum), MPI_INTEGER8, MPI_SUM, comm2d,  &
                           ierr)
#else
       global_sum = local_sum
#endif
       nsurfl_tot  = global_sum(1)
       npcbl_tot   = global_sum(2)
       ndsidir_tot = global_sum(3)
       nmrtbl_tot  = global_sum(4)
!
!--    Check for errors.
       nx_on_file = tgh%total_nx-1
       ny_on_file = tgh%total_ny-1

       IF ( nx_on_file /= nx  .OR.  ny_on_file /= ny )  THEN
          WRITE( message_string, '(A,4(A,I7))' )                                                   &
               'total number of grid points along x and y in file SVFIN do not match current run', &
               '&nx_on_file = ', nx_on_file, ' ny_on_file = ', ny_on_file, ' nx = ', nx, ' ny = ', &
               ny
          CALL message( 'radiation_read_svf', 'RAD0061', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( TRIM( rad_version_in ) /= TRIM( rad_version ) )  THEN
          WRITE( message_string, * ) 'Version of binary SVF file "', TRIM( rad_version_field ),    &
                                     '" does not match the version of model "',                    &
                                     TRIM( rad_version ), '"'
          CALL message( 'radiation_read_svf', 'RAD0054', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( nsvfl_tot < 0  .OR.  ncsfl_tot < 0 )  THEN
          WRITE( message_string, * ) 'Wrong number of SVF or CSF'
          CALL message( 'radiation_read_svf', 'RAD0055', 1, 2, 0, 6, 0 )
       ELSE
          WRITE( debug_string , * ) 'Number of SVF, CSF, and nsurfl to read', nsvfl_tot, ncsfl_tot,  &
                                          nsurfl_tot
          IF ( debug_output )  CALL debug_message( debug_string, 'info' )
       ENDIF

       IF ( nsurfl_from_file /= nsurfl_tot )  THEN
          WRITE( message_string, * ) 'nsurfl from SVF file does not match calculated ',        &
             'nsurfl from radiation_interaction_init'
          CALL message( 'radiation_read_svf', 'RAD0056', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( npcbl_from_file /= npcbl_tot )  THEN
          WRITE( message_string, * ) 'npcbl from SVF file does not match calculated npcbl ',   &
             'from radiation_interaction_init'
          CALL message( 'radiation_read_svf', 'RAD0057', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( ndsidir_from_file /= ndsidir_tot )  THEN
          WRITE( message_string, * ) 'ndsidir from SVF file does not match calculated ',       &
             'ndsidir from radiation_presimulate_solar_pos'
          CALL message( 'radiation_read_svf', 'RAD0058', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( nmrtbl_from_file /= nmrtbl_tot )  THEN
          WRITE( message_string, * ) 'nmrtbl from SVF file does not match calculated nmrtbl ', &
             'from radiation_interaction_init'
          CALL message( 'radiation_read_svf', 'RAD0059', 1, 2, 0, 6, 0 )
       ELSE
          WRITE( debug_string, * ) 'Number of nmrtf to read ', nmrtf_tot
          IF ( debug_output )  CALL debug_message( debug_string, 'info' )
       ENDIF

!
!--    Open MPI-IO SVF file for read local data.
       CALL rd_mpi_io_open( 'READ', 'SVFIN' // TRIM( coupling_char ) )

#if defined( __parallel )
!
!--    Reading core-dependent view-factor data sizes and restore them on the respective variable.
       ALLOCATE( tmp_int(0:numprocs-1) )
       tmp_int = 0

       CALL rrd_mpi_io_global_array( 'nsend_radx', tmp_int )
       nsend_radx = tmp_int(myid)

       CALL rrd_mpi_io_global_array( 'nrecv_radx', tmp_int )
       nrecv_radx = tmp_int(myid)

       CALL rrd_mpi_io_global_array( 'niters_radx', tmp_int )
       niters_radx = tmp_int(myid)

       CALL rrd_mpi_io_global_array( 'nmaxsend_radx', tmp_int )
       nmaxsend_radx = tmp_int(myid)

       DEALLOCATE( tmp_int )
!
!--    Allocate core-dependent buffer arrays.
       ALLOCATE( isurf_send_radx(nsend_radx) )
       ALLOCATE( isurf_recv_radx(nrecv_radx) )
       ALLOCATE( disp_send_radx(0:numprocs) )
       ALLOCATE( disp_recv_radx(0:numprocs) )
       ALLOCATE( radx_send(nsend_radx) )
       ALLOCATE( surfoutl_recv(nrecv_radx) )
       ALLOCATE( surfouts_recv(nrecv_radx) )
       ALLOCATE( radx_send_surfinl(nrecv_radx) )
       ALLOCATE( surfinl_recv(nsend_radx) )
       ALLOCATE( disp_sendbuf_radx(numprocs) )
       ALLOCATE( disp_recvbuf_radx(numprocs) )
       ALLOCATE( num_sendbuf_radx(numprocs) )
       ALLOCATE( num_recvbuf_radx(numprocs) )
!
!--    Read isurf_send_radx. Even though the array is of type integer, it is treated as a real
!--    (MPI-IO mechanism is only targeted towards real-type arrays.
!--    This requires a type conversion using the Fortran function NINT.
       CALL rrd_mpi_io( 'int_send_buf_global_start', global_start )
       CALL rrd_mpi_io( 'int_send_buf_global_end', global_end )

       CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start,       &
                                         global_end )

       ALLOCATE( tmp_1d(nsend_radx) )
       CALL rrd_mpi_io_surface( 'isurf_send_radx', tmp_1d )
       isurf_send_radx = NINT( tmp_1d )
       DEALLOCATE( tmp_1d )
!
!--    Read isurf_recv_radx.
       CALL rrd_mpi_io( 'int_recv_buf_global_start', global_start )
       CALL rrd_mpi_io( 'int_recv_buf_global_end', global_end )

       CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start,       &
                                         global_end )

       ALLOCATE( tmp_1d(nrecv_radx) )
       CALL rrd_mpi_io_surface( 'isurf_recv_radx', tmp_1d )
       isurf_recv_radx = NINT( tmp_1d )
       DEALLOCATE( tmp_1d )
!
!--    Read disp_send_radx.
       CALL rrd_mpi_io( 'n_source_proc_disp_global_start', global_start )
       CALL rrd_mpi_io( 'n_source_proc_disp_global_end', global_end )

       CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start,       &
                                         global_end )

       ALLOCATE( tmp_1d(0:numprocs) )
       CALL rrd_mpi_io_surface( 'disp_send_radx', tmp_1d )
       disp_send_radx = NINT( tmp_1d )
       DEALLOCATE( tmp_1d )
!
!--    Read disp_recv_radx.
       CALL rrd_mpi_io( 'n_target_proc_disp_global_start', global_start )
       CALL rrd_mpi_io( 'n_target_proc_disp_global_end', global_end )

       CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start,       &
                                         global_end )

       ALLOCATE( tmp_1d(0:numprocs) )
       CALL rrd_mpi_io_surface( 'disp_recv_radx', tmp_1d )
       disp_recv_radx = NINT( tmp_1d )
       DEALLOCATE( tmp_1d )
#endif

       IF ( nsurfl_tot > 0 )  THEN
!
!--       Read global indices.
          CALL rrd_mpi_io( 'nsurfl_global_start', global_start )
          CALL rrd_mpi_io( 'nsurfl_global_end', global_end )
!
!--       Set file types of variables and compute local indices.
          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start,    &
                                            global_end )
          nsurfl = end_index(nyn,nxr)

          IF ( data_to_read )  THEN

             CALL rrd_mpi_io_surface( 'skyvf', skyvf )
             CALL rrd_mpi_io_surface( 'skyvft', skyvft )
!
!--          To avoid another overlay of rrd_mpi_io_surface, dsitrans is read as REAL tmp_2d array.
!--          The order of dimensions of dsitrans is different to the order expected by
!--          rrd_mpi_io_surface. Therefor a tranpose of tmp_2d is required.
             ALLOCATE( tmp_2d(SIZE(dsitrans,2),SIZE(dsitrans,1)) )
             CALL rrd_mpi_io_surface( 'dsitrans', tmp_2d )
             dsitrans = TRANSPOSE( tmp_2d )
             DEALLOCATE( tmp_2d )

          ENDIF

       ENDIF

       IF ( npcbl_tot > 0 )  THEN

          CALL rrd_mpi_io( 'npcbl_global_start', global_start )
          CALL rrd_mpi_io( 'npcbl_global_end', global_end )

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start,    &
                                            global_end )
          npcbl = end_index(nyn,nxr)

          IF ( data_to_read )  THEN

             ALLOCATE( tmp_2d(SIZE(dsitransc,2),SIZE(dsitransc,1)) )
             CALL rrd_mpi_io_surface( 'dsitransc', tmp_2d )
             dsitransc = TRANSPOSE( tmp_2d )
             DEALLOCATE( tmp_2d )

          ENDIF

       ENDIF

       IF ( nsvfl_tot > 0 )  THEN

          CALL rrd_mpi_io( 'nsvfl_global_start', global_start )
          CALL rrd_mpi_io( 'nsvfl_global_end', global_end )

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start,    &
                                            global_end )
          nsvfl = end_index(nyn,nxr)

          IF ( .NOT. ALLOCATED( svf )     )  ALLOCATE( svf(ndsvf,nsvfl) )
          IF ( .NOT. ALLOCATED( svfsurf ) )  ALLOCATE( svfsurf(idsvf,nsvfl) )

          IF ( data_to_read )  THEN

             CALL rrd_mpi_io_surface( 'svf', svf )
             ALLOCATE( tmp_2d(SIZE(svfsurf,1),SIZE(svfsurf,2)) )
             CALL rrd_mpi_io_surface( 'svfsurf', tmp_2d )
             svfsurf = tmp_2d
             DEALLOCATE( tmp_2d )

          ENDIF

       ENDIF

       IF ( plant_canopy  )  THEN

          IF ( ncsfl_tot > 0 )  THEN

             CALL rrd_mpi_io( 'ncsfl_global_start', global_start )
             CALL rrd_mpi_io( 'ncsfl_global_end', global_end )

             CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start, &
                                               global_end )
             ncsfl = end_index(nyn,nxr)

             IF ( .NOT. ALLOCATED( csf )     )  ALLOCATE( csf(ndcsf,ncsfl) )
             IF ( .NOT. ALLOCATED( csfsurf ) )  ALLOCATE( csfsurf(idcsf,ncsfl) )

             IF ( data_to_read )  THEN

                CALL rrd_mpi_io_surface( 'csf', csf )
                ALLOCATE( tmp_2d(SIZE(csfsurf,1),SIZE(csfsurf,2)) )
                CALL rrd_mpi_io_surface( 'csfsurf', tmp_2d )
                csfsurf = tmp_2d
                DEALLOCATE( tmp_2d )

             ENDIF

          ENDIF

       ENDIF

       IF ( nmrtbl_tot > 0 )  THEN

          CALL rrd_mpi_io( 'nmrtbl_global_start', global_start )
          CALL rrd_mpi_io( 'nmrtbl_global_end', global_end )

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start,    &
                                            global_end )
          nmrtbl = end_index(nyn,nxr)

          IF ( data_to_read )  THEN

             CALL rrd_mpi_io_surface( 'mrtsky', mrtsky )
             CALL rrd_mpi_io_surface( 'mrtskyt', mrtskyt )
             ALLOCATE( tmp_2d(SIZE(mrtdsit,2),SIZE(mrtdsit,1)) )
             CALL rrd_mpi_io_surface( 'mrtdsit', tmp_2d )
             mrtdsit = TRANSPOSE( tmp_2d )
             DEALLOCATE( tmp_2d )

          ENDIF

       ENDIF

       IF ( nmrtf_tot > 0 )  THEN

          CALL rrd_mpi_io( 'nmrtf_global_start', global_start )
          CALL rrd_mpi_io( 'nmrtf_global_end', global_end )

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_read, global_start,    &
                                            global_end )
          nmrtf = end_index(nyn,nxr)

          IF ( .NOT. ALLOCATED( mrtf )     )  ALLOCATE( mrtf(nmrtf) )
          IF ( .NOT. ALLOCATED( mrtft )    )  ALLOCATE( mrtft(nmrtf) )
          IF ( .NOT. ALLOCATED( mrtfsurf ) )  ALLOCATE( mrtfsurf(2,nmrtf) )

          IF ( data_to_read )  THEN

             CALL rrd_mpi_io_surface( 'mrtf', mrtf )
             CALL rrd_mpi_io_surface( 'mrtft', mrtft )
             ALLOCATE( tmp_2d(SIZE(mrtfsurf,1),SIZE(mrtfsurf,2)) )
             CALL rrd_mpi_io_surface( 'mrtfsurf', tmp_2d )
             mrtfsurf = tmp_2d
             DEALLOCATE( tmp_2d )

          ENDIF

       ENDIF

       IF ( radiation_volumetric_flux )  THEN

          IF ( .NOT. ALLOCATED(shadow_top) )  ALLOCATE( shadow_top(nys:nyn,nxl:nxr,ndsidir) )

          IF ( SIZE( shadow_top, DIM=3 ) > 0 )  THEN
             CALL rrd_mpi_io( 'shadow_top', shadow_top, SIZE( shadow_top, DIM=3 ) )
          ENDIF

          ALLOCATE( tmp_3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'skyvf_vol', tmp_3d )
          skyvf_vol = tmp_3d(nz_urban_b:nz_urban_t,nys:nyn,nxl:nxr)
          DEALLOCATE( tmp_3d )

       ENDIF

       CALL rd_mpi_io_close

       CALL location_message( 'reading sky view factors in MPI-IO format', 'finished' )

    ENDIF

 END SUBROUTINE radiation_read_svf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_rrd_global_ftn( found )

    LOGICAL, INTENT(OUT) ::  found !< control flag to indicate whether a variable has been found in the module or not


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'albedo_eff' )
          READ( 13 )  albedo_eff

       CASE ( 'dt_radiation' )
          READ( 13 )  dt_radiation

       CASE ( 'emissivity_eff' )
          READ( 13 )  emissivity_eff

       CASE ( 't_rad_eff' )
          READ( 13 )  t_rad_eff

       CASE ( 'time_radiation' )
          READ( 13 )  time_radiation

       CASE DEFAULT
          found = .FALSE.

    END SELECT

 END SUBROUTINE radiation_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_rrd_global_mpi

    CALL rrd_mpi_io( 'albedo_eff',     albedo_eff     )
    CALL rrd_mpi_io( 'dt_radiation',   dt_radiation   )
    CALL rrd_mpi_io( 'emissivity_eff', emissivity_eff )
    CALL rrd_mpi_io( 't_rad_eff',      t_rad_eff      )
    CALL rrd_mpi_io( 'time_radiation', time_radiation )

 END SUBROUTINE radiation_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine stores svf, svfsurf, csf, csfsurf and mrt data to a file. The stored factors can be
!> reused in future simulation with the same geometry structure of the surfaces and resolved plant
!> canopy.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_write_svf

    IMPLICIT NONE

    CHARACTER(LEN=16), DIMENSION(7) ::  counter_name  !< names of sky view factor counter

    INTEGER(iwp), PARAMETER ::  max_i4_value  = 2147483647  !< maximum positive INTEGER(4) value, 2**31 - 1

    INTEGER(iwp) ::  i            !<
    INTEGER(iwp) ::  ic           !<
#if defined( __parallel )
    INTEGER(iwp) ::  ierr         !<
#endif
    INTEGER(iwp) ::  ind          !<
    INTEGER(iwp) ::  ipcgb        !<
    INTEGER(iwp) ::  isurf        !<
    INTEGER(iwp) ::  j            !<
    INTEGER(iwp) ::  jc           !<
    INTEGER(idp) ::  ncsfl_tot    !<
    INTEGER(idp) ::  ndsidir_tot  !<
    INTEGER(idp) ::  nmrtbl_tot   !<
    INTEGER(idp) ::  nmrtf_tot    !<
    INTEGER(idp) ::  npcbl_tot    !<
    INTEGER(idp) ::  nsurfl_tot   !<
    INTEGER(idp) ::  nsvfl_tot    !< total (sum over all PEs) counter for the different IO variables
#if defined( __parallel )
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  tmp_int !< temporary array to store 0d dimension sizes
#endif

    INTEGER(idp), DIMENSION(7) ::  global_sum  !< idp to allow check, if total number of values > 2G
    INTEGER(idp), DIMENSION(7) ::  local_sum   !< variables to compute total counter

    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  end_index     !< local end index
    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_end    !< global end index
    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_start  !< global start index
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  lo_no         !< local number of values
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  start_index   !< local start index

    LOGICAL ::  data_to_write                         !< flag indicating if data is available for writing

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::  tmp  !<

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  tmp_3d  !< dummy array used for reordering of 3d arrays


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL location_message( 'writing sky view factors in Fortran binary format', 'start' )

       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN
!
!--          Open binary file
             CALL check_open( 89 )

             WRITE( 89 ) rad_version
             WRITE( 89 ) nsvfl, ncsfl, nsurfl, npcbl, ndsidir, nmrtbl, nmrtf
             IF ( nsurfl > 0 )  THEN
                WRITE( 89 ) skyvf
                WRITE( 89 ) skyvft
                WRITE( 89 ) dsitrans
             ENDIF
             IF ( npcbl > 0 )  THEN
                WRITE( 89 ) dsitransc
             ENDIF
             IF ( nsvfl > 0 )  THEN
                WRITE( 89 ) svf
                WRITE( 89 ) svfsurf
             ENDIF
             IF ( plant_canopy  .AND.  ncsfl > 0 )  THEN
                 WRITE( 89 ) csf
                 WRITE( 89 ) csfsurf
             ENDIF
             IF ( nmrtbl > 0 )  THEN
                WRITE( 89 ) mrtsky
                WRITE( 89 ) mrtskyt
                WRITE( 89 ) mrtdsit
             ENDIF
             IF ( nmrtf > 0 )  THEN
                 WRITE( 89 ) mrtf
                 WRITE( 89 ) mrtft
                 WRITE( 89 ) mrtfsurf
             ENDIF
#if defined( __parallel )
             WRITE( 89 ) nsend_radx, nrecv_radx, niters_radx, nmaxsend_radx
             WRITE( 89 ) isurf_send_radx
             WRITE( 89 ) isurf_recv_radx
             WRITE( 89 ) disp_send_radx
             WRITE( 89 ) disp_recv_radx
#endif
             IF ( radiation_volumetric_flux ) THEN
                 WRITE( 89 ) shadow_top
                 WRITE( 89 ) skyvf_vol
             ENDIF
!
!--          Close binary file
             CALL close_file( 89 )

          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

       CALL location_message( 'writing sky view factors in Fortran binary format', 'finished' )

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
!
!--    Sky view factor data is treated like surface data and the respective routines from
!--    restart_data_mpi_io_mod are used for writing. In contrary to restart data, the virtual
!--    processor grid can not be changed between writing and reading svf data.
       CALL location_message( 'writing sky view factors in MPI-IO format', 'start' )

!
!--    Open MPI-IO SVF file. It must not contain outer boundaries of the total domain.
       CALL rd_mpi_io_open( 'write', 'SVFOUT' // TRIM( coupling_char ),                            &
                            write_total_domain_boundaries = .FALSE. )
!
!--    Write global variables.
       CALL wrd_mpi_io( 'rad_version', rad_version )
!
!--    Set names of sky view counters
       counter_name(1) = 'nsvfl '
       counter_name(2) = 'ncsfl '
       counter_name(3) = 'nsurfl '
       counter_name(4) = 'npcbl '
       counter_name(5) = 'ndsidir '
       counter_name(6) = 'nmrtbl '
       counter_name(7) = 'nmrtf '
!
!--    Sum local number of skyview factor values on all PEs.
       local_sum(1) = nsvfl
       local_sum(2) = ncsfl
       local_sum(3) = nsurfl
       local_sum(4) = npcbl
       local_sum(5) = ndsidir
       local_sum(6) = nmrtbl
       local_sum(7) = nmrtf
#if defined( __parallel )
       CALL MPI_ALLREDUCE( local_sum, global_sum, SIZE( local_sum ), MPI_INTEGER8, MPI_SUM,        &
                           comm2d, ierr)
#else
       global_sum = local_sum
#endif
!
!--    Check, if total number of respective skyview values do not exceed 2**31-1
       DO  i = 1, 7
          IF ( global_sum(i) > max_i4_value )  THEN
             WRITE( message_string, '(A,A,I12,A)' ) 'number of sky view factor values for ',       &
                                    TRIM( counter_name(i) ) // ' = ',  global_sum(i),' is > 2**31-1'
             CALL message( 'radiation_write_svf', 'RAD0062', 0, 0, 0, 6, 0 )
          ENDIF
       ENDDO

       nsvfl_tot   = global_sum(1)
       ncsfl_tot   = global_sum(2)
       nsurfl_tot  = global_sum(3)
       npcbl_tot   = global_sum(4)
       ndsidir_tot = global_sum(5)
       nmrtbl_tot  = global_sum(6)
       nmrtf_tot   = global_sum(7)

#if defined( __parallel )
!
!--    Write core-dependent view-factor data sizes. Gather this data from all cores before
!--    writing to the output file.
       ALLOCATE( tmp_int(0:numprocs-1) )

       tmp_int = 0
       tmp_int(myid) = nsend_radx
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_int, SIZE( tmp_int ), MPI_INTEGER, MPI_SUM, comm2d, ierr)
       CALL wrd_mpi_io_global_array( 'nsend_radx', tmp_int )

       tmp_int = 0
       tmp_int(myid) = nrecv_radx
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_int, SIZE( tmp_int ), MPI_INTEGER, MPI_SUM, comm2d, ierr)
       CALL wrd_mpi_io_global_array( 'nrecv_radx', tmp_int )

       tmp_int = 0
       tmp_int(myid) = niters_radx
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_int, SIZE( tmp_int ), MPI_INTEGER, MPI_SUM, comm2d, ierr)
       CALL wrd_mpi_io_global_array( 'niters_radx', tmp_int )

       tmp_int = 0
       tmp_int(myid) = nmaxsend_radx
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_int, SIZE( tmp_int ), MPI_INTEGER, MPI_SUM, comm2d, ierr)
       CALL wrd_mpi_io_global_array( 'nmaxsend_radx', tmp_int )
       DEALLOCATE( tmp_int )
!
!--    To write core-dependent view-factor data, use the surface-data IO mechanism, even though
!--    this data is not really surface data. All the data arrays are referred to start- and end
!--    index (nys,nxl). As the restart MPI-IO mechanism is currently only targeted to read and
!--    write real-type surface data, convert the integer to real type.
       IF ( ALLOCATED( isurf_send_radx ) )  THEN
          start_index = nsend_radx + 1
          end_index   = nsend_radx
          start_index(nys,nxl) = 1
!
!--       Write the data. First, compute global start and end indices.
          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,   &
                                            global_end )
          CALL wrd_mpi_io( 'int_send_buf_global_start', global_start )
          CALL wrd_mpi_io( 'int_send_buf_global_end', global_end )
          IF ( data_to_write )  CALL wrd_mpi_io_surface( 'isurf_send_radx',                        &
                                                         REAL( isurf_send_radx, KIND = wp ) )
       ENDIF

       IF ( ALLOCATED( isurf_recv_radx ) )  THEN
          start_index = nrecv_radx + 1
          end_index   = nrecv_radx
          start_index(nys,nxl) = 1

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,   &
                                            global_end )
          CALL wrd_mpi_io( 'int_recv_buf_global_start', global_start )
          CALL wrd_mpi_io( 'int_recv_buf_global_end', global_end )
          IF ( data_to_write )  CALL wrd_mpi_io_surface( 'isurf_recv_radx',                        &
                                                         REAL( isurf_recv_radx, KIND = wp ) )
       ENDIF

       IF ( ALLOCATED( disp_send_radx ) )  THEN
          start_index = SIZE( disp_send_radx ) + 1
          end_index   = SIZE( disp_send_radx )
          start_index(nys,nxl) = 1

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,   &
                                            global_end )
          CALL wrd_mpi_io( 'n_source_proc_disp_global_start', global_start )
          CALL wrd_mpi_io( 'n_source_proc_disp_global_end', global_end )
          IF ( data_to_write )  CALL wrd_mpi_io_surface( 'disp_send_radx',                         &
                                                         REAL( disp_send_radx, KIND = wp ) )
       ENDIF

       IF ( ALLOCATED( disp_recv_radx ) )  THEN
          start_index = SIZE( disp_recv_radx ) + 1
          end_index   = SIZE( disp_recv_radx )
          start_index(nys,nxl) = 1

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,   &
                                            global_end )
          CALL wrd_mpi_io( 'n_target_proc_disp_global_start', global_start )
          CALL wrd_mpi_io( 'n_target_proc_disp_global_end', global_end )
          IF ( data_to_write )  CALL wrd_mpi_io_surface( 'disp_recv_radx',                         &
                                                         REAL( disp_recv_radx, KIND = wp ) )
       ENDIF
#endif
!
!--    Write total counters in header section of MPI-IO file.
       CALL wrd_mpi_io( 'nsvfl', nsvfl_tot )
       CALL wrd_mpi_io( 'ncsfl', ncsfl_tot )
       CALL wrd_mpi_io( 'nsurfl', nsurfl_tot )
       CALL wrd_mpi_io( 'npcbl', npcbl_tot )
       CALL wrd_mpi_io( 'ndsidir', ndsidir_tot )
       CALL wrd_mpi_io( 'nmrtbl', nmrtbl_tot )
       CALL wrd_mpi_io( 'nmrtf', nmrtf_tot )
!
!--    Write local data.
!--    All svf values are treated as surface values and use the respective routines from
!--    restart_data_mpi_io_mod.
       IF ( nsurfl > 0 )  THEN

          lo_no = 0
!
!--       Count surface values on individual grid cells
          DO  isurf = 1, nsurfl
             jc = surfl(iy, isurf)
             ic = surfl(ix, isurf)
             lo_no(jc,ic) = lo_no(jc,ic) + 1
          ENDDO
!
!--       Create local index array similar to surface routines.
          ind = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                start_index(j,i) = ind
                end_index(j,i)   = start_index(j,i) + lo_no(j,i) - 1
                ind = ind + lo_no(j,i)
             ENDDO
          ENDDO

       ELSE

          start_index = 1
          end_index = 0

       ENDIF
!
!--    Each PE has to call the next block, therefore nsurfl_tot is used.
!--    This is required to use MPI_FILE_WRITE_ALL for writing.
       IF ( nsurfl_tot > 0 )  THEN
!
!--       Set file types of variables for this block and compute global indices.
          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,   &
                                            global_end )

          CALL wrd_mpi_io( 'nsurfl_global_start', global_start )
          CALL wrd_mpi_io( 'nsurfl_global_end', global_end )

          IF ( data_to_write )  THEN

             CALL wrd_mpi_io_surface( 'skyvf', skyvf )
             CALL wrd_mpi_io_surface( 'skyvft', skyvft )
!
!--          To avoid another overlay of rrd_mpi_io_surface, dsitrans is written as REAL tmp array.
!--          The order of dimensions of dsitrans is different to the order expected by
!--          rrd_mpi_io_surface. Therefor a tranpose of tmp is required.
             ALLOCATE( tmp(SIZE(dsitrans,2),SIZE(dsitrans,1)) )
             tmp = TRANSPOSE( dsitrans )
             CALL wrd_mpi_io_surface( 'dsitrans', tmp )
             DEALLOCATE( tmp )

          ENDIF

       ENDIF

       IF ( npcbl > 0 )  THEN

          lo_no = 0
          DO  isurf = 1, npcbl
             jc = pcbl(iy,isurf)
             ic = pcbl(ix,isurf)
             lo_no(jc,ic) = lo_no(jc,ic) + 1
          ENDDO

          ind = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                start_index(j,i) = ind
                end_index(j,i)   = start_index(j,i) + lo_no(j,i) - 1
                ind = ind+lo_no(j,i)
             ENDDO
          ENDDO

       ELSE

          start_index = 1
          end_index = 0

       ENDIF

       IF ( npcbl_tot > 0 )  THEN

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,   &
                                            global_end )
          CALL wrd_mpi_io( 'npcbl_global_start', global_start )
          CALL wrd_mpi_io( 'npcbl_global_end', global_end )

          IF ( data_to_write )  THEN

             ALLOCATE( tmp(SIZE(dsitransc,2),SIZE(dsitransc,1)) )
             tmp = TRANSPOSE( dsitransc )
             CALL wrd_mpi_io_surface( 'dsitransc', tmp )
             DEALLOCATE( tmp )
          ENDIF
       ENDIF

       IF ( nsvfl > 0 )  THEN

          lo_no = 0
          DO  j = 1, SIZE( svfsurf, 2 )
             isurf = svfsurf(1,j)
             jc = surfl(iy,isurf)
             ic = surfl(ix,isurf)
             lo_no(jc,ic) = lo_no(jc,ic) + 1
          ENDDO

          ind = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                start_index(j,i) = ind
                end_index(j,i)   = start_index(j,i) + lo_no(j,i) - 1
                ind = ind+lo_no(j,i)
             ENDDO
          ENDDO

       ELSE

          start_index = 1
          end_index = 0
       ENDIF

       IF ( nsvfl_tot > 0 )  THEN

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,   &
                                            global_end )
          CALL wrd_mpi_io( 'nsvfl_global_start', global_start )
          CALL wrd_mpi_io( 'nsvfl_global_end', global_end )

          IF ( data_to_write )  THEN

             CALL wrd_mpi_io_surface( 'svf', svf )
             ALLOCATE( tmp(SIZE(svfsurf,1),SIZE(svfsurf,2)) )
             tmp(:,:) = svfsurf(:,:)
             CALL wrd_mpi_io_surface( 'svfsurf', tmp )
             DEALLOCATE( tmp )

          ENDIF

       ENDIF

       IF ( plant_canopy  )  THEN

          lo_no = 0
          IF ( ncsfl > 0 )  THEN

             DO  j = 1, ncsfl
                ipcgb = csfsurf(1, j)
                jc = pcbl(iy,ipcgb)
                ic = pcbl(ix,ipcgb)
                lo_no(jc,ic) = lo_no(jc,ic) + 1
             ENDDO

             ind = 1
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   start_index(j,i) = ind
                   end_index(j,i)   = start_index(j,i) + lo_no(j,i) - 1
                   ind = ind+lo_no(j,i)
                ENDDO
             ENDDO

          ELSE

             start_index = 1
             end_index = -1

          ENDIF

          IF ( ncsfl_tot > 0 )  THEN

             CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,&
                                               global_end )
             CALL wrd_mpi_io( 'ncsfl_global_start', global_start )
             CALL wrd_mpi_io( 'ncsfl_global_end', global_end )

             IF ( data_to_write )  THEN

                IF ( ALLOCATED( csf ) )  THEN
                   CALL wrd_mpi_io_surface( 'csf', csf )
                ELSE
                   ALLOCATE( tmp(ndcsf,0) )
                   CALL wrd_mpi_io_surface( 'csf', tmp )
                   DEALLOCATE( tmp )
                ENDIF

                IF ( ALLOCATED( csfsurf ) )  THEN
                   ALLOCATE( tmp(SIZE(csfsurf,1),SIZE(csfsurf,2)) )
                   tmp(:,:) = csfsurf(:,:)
                   CALL wrd_mpi_io_surface( 'csfsurf', tmp )
                   DEALLOCATE( tmp )
                ELSE
                   ALLOCATE( tmp(idcsf,ncsfl) )
                   CALL wrd_mpi_io_surface( 'csfsurf', tmp )
                   DEALLOCATE( tmp )
                ENDIF

             ENDIF

          ENDIF

       ENDIF

       IF ( nmrtbl > 0 )  THEN

          lo_no = 0
          DO  j = 1, nmrtbl
             jc = mrtbl(iy,j)
             ic = mrtbl(ix,j)
             lo_no(jc,ic) = lo_no(jc,ic) + 1
          ENDDO

          ind = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                start_index(j,i) = ind
                end_index(j,i)   = start_index(j,i) + lo_no(j,i) - 1
                ind = ind+lo_no(j,i)
             ENDDO
          ENDDO

       ELSE

          start_index = 1
          end_index = 0
       ENDIF

       IF ( nmrtbl_tot > 0 )  THEN

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,   &
                                            global_end )
          CALL wrd_mpi_io( 'nmrtbl_global_start', global_start )
          CALL wrd_mpi_io( 'nmrtbl_global_end', global_end )

          IF ( data_to_write )  THEN

             CALL wrd_mpi_io_surface( 'mrtsky', mrtsky )
             CALL wrd_mpi_io_surface( 'mrtskyt', mrtskyt )
             ALLOCATE( tmp(SIZE(mrtdsit,2),SIZE(mrtdsit,1)) )
             tmp = TRANSPOSE( mrtdsit )
             CALL wrd_mpi_io_surface( 'mrtdsit', tmp )
             DEALLOCATE( tmp )

          ENDIF

       ENDIF

       IF ( nmrtf > 0 )  THEN

          lo_no = 0
          DO  j = 1, nmrtf
             isurf = mrtfsurf(1,j)
             jc = surfl(iy,isurf)
             ic = surfl(ix,isurf)
             lo_no(jc,ic) = lo_no(jc,ic) + 1
          ENDDO

          ind = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                start_index(j,i) = ind
                end_index(j,i)   = start_index(j,i) + lo_no(j,i) - 1
                ind = ind+lo_no(j,i)
             ENDDO
          ENDDO

       ELSE
          start_index = 1
          end_index = 0
       ENDIF

       IF ( nmrtf_tot > 0 )  THEN

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, data_to_write, global_start,   &
                                            global_end)
          CALL wrd_mpi_io( 'nmrtf_global_start', global_start )
          CALL wrd_mpi_io( 'nmrtf_global_end', global_end )

          IF ( data_to_write )  THEN

             CALL wrd_mpi_io_surface ( 'mrtf', mrtf )
             CALL wrd_mpi_io_surface ( 'mrtft', mrtft )

             ALLOCATE( tmp(SIZE(mrtfsurf,1),SIZE(mrtfsurf,2)) )
             tmp = mrtfsurf
             CALL wrd_mpi_io_surface( 'mrtfsurf', tmp )
             DEALLOCATE( tmp )

          ENDIF

       ENDIF

       IF ( radiation_volumetric_flux )  THEN

          IF ( SIZE( shadow_top, DIM=3 ) > 0 )  THEN
             CALL wrd_mpi_io( 'shadow_top', shadow_top, SIZE( shadow_top, DIM=3 ) )
          ENDIF

          ALLOCATE( tmp_3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          tmp_3d = 0.0_wp
          tmp_3d(nz_urban_b:nz_urban_t,nys:nyn,nxl:nxr) = skyvf_vol
          CALL wrd_mpi_io( 'skyvf_vol', tmp_3d )
          DEALLOCATE( tmp_3d )

       ENDIF

       CALL rd_mpi_io_close

       CALL location_message( 'writing sky view factors in MPI-IO format', 'finished' )

    ENDIF

 END SUBROUTINE radiation_write_svf


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Block of auxiliary subroutines for RTM:
!> 1. quicksort and corresponding comparison
!> 2. merge_and_grow_csf for implementation of "dynamical growing" array for csf
!--------------------------------------------------------------------------------------------------!
!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_itarget( itarget, vffrac, ztransp, first, last )

    IMPLICIT NONE

    INTEGER(iwp) ::  i, j  !<
    INTEGER(iwp) ::  x, t  !<

    INTEGER(iwp), INTENT(IN) :: first, last  !<

    INTEGER(iwp), DIMENSION(:), INTENT(INOUT) ::  itarget  !<

    REAL(wp) ::  tr  !<

    REAL(wp), DIMENSION(:), INTENT(INOUT) ::  vffrac, ztransp  !<

    IF ( first >= last )  RETURN
    x = itarget(( first + last ) / 2)
    i = first
    j = last
    DO
        DO WHILE ( itarget(i) < x )
           i = i+1
        ENDDO
        DO WHILE ( x < itarget(j) )
            j = j-1
        ENDDO
        IF ( i >= j ) EXIT
        t = itarget(i);  itarget(i) = itarget(j);  itarget(j) = t
        tr = vffrac(i);  vffrac(i) = vffrac(j);  vffrac(j) = tr
        tr = ztransp(i);  ztransp(i) = ztransp(j);  ztransp(j) = tr
        i = i+1
        j = j-1
    ENDDO
    IF ( first < i-1 )  CALL quicksort_itarget( itarget, vffrac, ztransp, first, i - 1 )
    IF ( j+1 < last )  CALL quicksort_itarget( itarget, vffrac, ztransp, j + 1, last )

 END SUBROUTINE quicksort_itarget


 PURE FUNCTION csf_lt( csf1, csf2 ) RESULT( res )

    LOGICAL ::  res  !<

    TYPE(t_csf), INTENT(in) ::  csf1,csf2  !<

    IF ( csf1%ip < csf2%ip     .OR.                                                                &
       ( csf1%ip == csf2%ip    .AND.  csf1%itx < csf2%itx )  .OR.                                  &
       ( csf1%ip == csf2%ip    .AND.  csf1%itx == csf2%itx   .AND.  csf1%ity < csf2%ity )  .OR.    &
       ( csf1%ip == csf2%ip    .AND.  csf1%itx == csf2%itx   .AND.  csf1%ity == csf2%ity   .AND.   &
         csf1%itz < csf2%itz ) .OR.                                                                &
       ( csf1%ip == csf2%ip    .AND.  csf1%itx == csf2%itx   .AND.  csf1%ity == csf2%ity   .AND.   &
        csf1%itz == csf2%itz   .AND.  csf1%isurfs < csf2%isurfs ) )  THEN
        res = .TRUE.
    ELSE
        res = .FALSE.
    ENDIF

 END FUNCTION csf_lt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @Todo: Missing subroutine description!
!--------------------------------------------------------------------------------------------------!
!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_csf( csfl, first, last )

    IMPLICIT NONE

    INTEGER(iwp) ::  i, j  !<

    INTEGER(iwp), INTENT(IN) ::  first, last  !<

    TYPE(t_csf) ::  x, t  !<

    TYPE(t_csf), DIMENSION(:), INTENT(INOUT) ::  csfl  !<


    IF ( first >= last )  RETURN
    x = csfl(( first + last ) / 2)
    i = first
    j = last
    DO
        DO WHILE ( csf_lt(csfl(i),x) )
            i = i+1
        ENDDO
        DO WHILE ( csf_lt(x,csfl(j)) )
            j = j-1
        ENDDO
        IF ( i >= j )  EXIT
        t = csfl(i);  csfl(i) = csfl(j);  csfl(j) = t
        i = i+1
        j = j-1
    ENDDO
    IF ( first < i-1 )  CALL quicksort_csf( csfl, first, i-1 )
    IF ( j+1 < last )  CALL quicksort_csf( csfl, j+1, last )

 END SUBROUTINE quicksort_csf


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Grows the CSF array in RTM exponentially when it is full. During that, the ray canopy sink
!> factors with common source face and target plant canopy grid cell are merged together so that the
!> size doesn't grow out of control.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE merge_and_grow_csf( newsize )

   INTEGER(iwp) ::  iread, iwrite  !<

   INTEGER(iwp), INTENT(IN) ::  newsize  !< new array size after grow, must be >= ncsfl
                                         !< or -1 to shrink to minimum

    TYPE(t_csf), DIMENSION(:), POINTER ::  acsfnew  !<


    IF ( newsize == -1 )  THEN
!
!--     Merge in-place
        acsfnew => acsf
    ELSE
!
!--     Allocate new array
        IF ( mcsf == 0 )  THEN
            ALLOCATE( acsf1(newsize) )
            acsfnew => acsf1
        ELSE
            ALLOCATE( acsf2(newsize) )
            acsfnew => acsf2
        ENDIF
    ENDIF

    IF ( ncsfl >= 1 )  THEN
!
!--     Sort csf in place (quicksort)
        CALL quicksort_csf( acsf, 1, ncsfl )
!
!--     While moving to a new array, aggregate canopy sink factor records with identical box & source
        acsfnew(1) = acsf(1)
        iwrite = 1
        DO  iread = 2, ncsfl
!
!--         Here acsf(kcsf) already has values from acsf(icsf)
            IF ( acsfnew(iwrite)%itx == acsf(iread)%itx                                            &
                 .AND.  acsfnew(iwrite)%ity == acsf(iread)%ity                                     &
                 .AND.  acsfnew(iwrite)%itz == acsf(iread)%itz                                     &
                 .AND.  acsfnew(iwrite)%isurfs == acsf(iread)%isurfs )  THEN

                acsfnew(iwrite)%rcvf = acsfnew(iwrite)%rcvf + acsf(iread)%rcvf
!
!--         Advance reading index, keep writing index
            ELSE
!
!--             Not identical, just advance and copy
                iwrite = iwrite + 1
                acsfnew(iwrite) = acsf(iread)
            ENDIF
        ENDDO
        ncsfl = iwrite
    ENDIF

    IF ( newsize == -1 )  THEN
!
!--     Allocate new array and copy shrinked data
        IF ( mcsf == 0 )  THEN
            ALLOCATE( acsf1(ncsfl) )
            acsf1(1:ncsfl) = acsf2(1:ncsfl)
        ELSE
            ALLOCATE( acsf2(ncsfl) )
            acsf2(1:ncsfl) = acsf1(1:ncsfl)
        ENDIF
    ENDIF
!
!-- Deallocate old array
    IF ( mcsf == 0 )  THEN
        mcsf = 1
        acsf => acsf1
        DEALLOCATE( acsf2 )
    ELSE
        mcsf = 0
        acsf => acsf2
        DEALLOCATE( acsf1 )
    ENDIF
    ncsfla = newsize

    IF ( debug_output )  THEN
       WRITE( debug_string, '(A,2I12)' ) 'Grow acsf2:', ncsfl, ncsfla
       CALL debug_message( debug_string, 'info' )
    ENDIF

 END SUBROUTINE merge_and_grow_csf


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> @Todo: Missing subroutine description!
!--------------------------------------------------------------------------------------------------!
!-- quicksort.f -*-f90-*-
!-- Author: t-nissie, adaptation J.Resler
!-- License: GPLv3
!-- Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE quicksort_csf2( kpcsflt, pcsflt, first, last )

    IMPLICIT NONE

    INTEGER(iwp) ::  i, j  !<

    INTEGER(iwp), INTENT(IN) ::  first, last  !<

    INTEGER(iwp), DIMENSION(kdcsf) ::  x, t1  !<

    INTEGER(iwp), DIMENSION(:,:), INTENT(INOUT) ::  kpcsflt  !<

    REAL(wp), DIMENSION(ndcsf) ::  t2  !<

    REAL(wp), DIMENSION(:,:), INTENT(INOUT) ::  pcsflt  !<


    IF ( first >= last )  RETURN
    x = kpcsflt(:, ( first + last ) / 2 )
    i = first
    j = last
    DO
        DO WHILE ( csf_lt2(kpcsflt(:,i),x) )
            i = i+1
        ENDDO
        DO WHILE ( csf_lt2(x,kpcsflt(:,j)) )
            j = j-1
        ENDDO
        IF ( i >= j )  EXIT
        t1 = kpcsflt(:,i);  kpcsflt(:,i) = kpcsflt(:,j);  kpcsflt(:,j) = t1
        t2 = pcsflt(:,i);  pcsflt(:,i) = pcsflt(:,j);  pcsflt(:,j) = t2
        i=i+1
        j=j-1
    ENDDO
    IF ( first < i-1 )  CALL quicksort_csf2( kpcsflt, pcsflt, first, i-1 )
    IF ( j+1 < last )  CALL quicksort_csf2( kpcsflt, pcsflt, j+1, last )

 END SUBROUTINE quicksort_csf2


 PURE FUNCTION csf_lt2( item1, item2 ) result( res )

    INTEGER(iwp), DIMENSION(kdcsf), INTENT(IN)  :: item1, item2  !<

    LOGICAL :: res  !<

    res = ( ( item1(3) < item2(3) )  .OR.  ( item1(3) == item2(3)  .AND.  item1(2) < item2(2) )    &
            .OR.  ( item1(3) == item2(3)  .AND.  item1(2) == item2(2)  .AND.  item1(1) < item2(1) )&
            .OR.  ( item1(3) == item2(3)  .AND.  item1(2) == item2(2)  .AND.  item1(1) == item2(1) &
            .AND.  item1(4) < item2(4) ) )

 END FUNCTION csf_lt2


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> For given coordinates, calculates indices within a global 3D (or 2D if nlayers=1) field, e.g. an
!> MPI one-sided window or an array which has been created using e.g. MPI_ALLGATHER.
!--------------------------------------------------------------------------------------------------!
 PURE SUBROUTINE radiation_calc_global_offset( i, j, k, nlayers, iproc, offs_proc, offs_glob )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN)            ::  i        !< x-coordinate
    INTEGER(iwp), INTENT(IN)            ::  j        !< y-coordinate
    INTEGER(iwp), INTENT(IN)            ::  k        !< z-coordinate
    INTEGER(iwp), INTENT(IN)            ::  nlayers  !< number of z-layers
    INTEGER(iwp), INTENT(OUT), OPTIONAL ::  iproc    !< MPI process rank
#if defined( __parallel )
    INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(OUT), OPTIONAL ::  offs_proc   !< offset within MPI proc
#else
    INTEGER(iwp), INTENT(OUT), OPTIONAL ::  offs_proc   !(actually unused without __parallel)
#endif
    INTEGER(iwp), INTENT(OUT), OPTIONAL ::  offs_glob   !< global offset

    INTEGER(iwp) ::  iproc_l  !< local variable for iproc
    INTEGER(iwp) ::  oproc_l  !< local variable for offs_proc


    iproc_l = ipx(i) * npey + ipy(j)
    IF ( PRESENT( iproc ) )  iproc = iproc_l
    IF ( PRESENT( offs_proc )  .OR.  PRESENT( offs_glob ) )  THEN
       oproc_l = (i - nxl_pe(ipx(i))) * (nyn_pe(ipy(j)) - nys_pe(ipy(j)) + 1) * nlayers +          & ! columns before
                 (j - nys_pe(ipy(j))) * nlayers                                         +          & ! rows in column
                 k
       IF ( PRESENT( offs_proc ) )  offs_proc = oproc_l
       IF ( PRESENT( offs_glob ) )  THEN
          IF ( iproc_l == 0 )  THEN
              offs_glob = oproc_l
          ELSE
              offs_glob = nnxyd(iproc_l) * nlayers + oproc_l
          ENDIF
       ENDIF
    ENDIF

 END SUBROUTINE radiation_calc_global_offset


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Performs MPI ALLTOALL exchange for integer and floating-point data, optionally splitting the
!> exchange to multiple iterations with maximum number of items per iteration.
!--------------------------------------------------------------------------------------------------!
#if defined( __parallel )
 SUBROUTINE radiation_exchange_alltoall( ntosend, npint, npfloat, isendbuf, fsendbuf, nrecv,       &
                                         irecvbuf, frecvbuf )

    IMPLICIT NONE

    INTEGER(iwp) ::  i, j      !< iterators
    INTEGER(iwp) ::  iproc     !< process iterator
    INTEGER(iwp) ::  iter      !< current iteration
    INTEGER(iwp) ::  niters    !< local number of iterations needed
    INTEGER(iwp) ::  nitersg   !< global no. of iterations needed
    INTEGER(iwp) ::  nmaxsend  !< max no. of records sent to each process in each iteration

    INTEGER(iwp), INTENT(IN) ::  npint    !< no. of integers in a record
    INTEGER(iwp), INTENT(IN) ::  npfloat  !< no. of floats in a record

    INTEGER(iwp), INTENT(OUT) ::  nrecv  !< total no. of records received

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  drecv        !< received data displacements per proc
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  drecvnow     !< current receive displacements
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  drecvnow_np  !< drecvnow times npint or npfloat
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  dsend        !< sent data displacements per process
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  dsendnow     !< current send displacements
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  dsendnow_np  !< dsendnow times npint or npfloat
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nrecvnow     !< no. of items to receive in current iteration
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nrecvnow_np  !< nrecvnow times npint or npfloat
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nsendnow     !< no. of items to send in current iteration
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  nsendnow_np  !< nsendnow times npint or npfloat
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ntorecv      !< no. of records to receive from each process

    INTEGER(iwp), DIMENSION(0:), INTENT(IN) ::  isendbuf  !< send buffer with integers
    INTEGER(iwp), DIMENSION(0:), INTENT(IN) ::  ntosend   !< number of records to send to each process

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::  irecvbuf  !< int receive buffer (will be allocated to proper size)

    REAL(wp), DIMENSION(0:), INTENT(IN) ::  fsendbuf  !< send buffer with floats

    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) ::  frecvbuf  !< float receive buffer


    ALLOCATE( dsend(0:numprocs) )
    ALLOCATE( drecv(0:numprocs) )
    ALLOCATE( drecvnow(0:numprocs-1) )
    ALLOCATE( drecvnow_np(0:numprocs-1) )
    ALLOCATE( dsendnow(0:numprocs-1) )
    ALLOCATE( dsendnow_np(0:numprocs-1) )
    ALLOCATE( nrecvnow(0:numprocs-1) )
    ALLOCATE( nrecvnow_np(0:numprocs-1) )
    ALLOCATE( nsendnow(0:numprocs-1) )
    ALLOCATE( nsendnow_np(0:numprocs-1) )
    ALLOCATE( ntorecv(0:numprocs-1) )
!
!-- Exchange send and receive sizes
    CALL MPI_ALLTOALL( ntosend, 1, MPI_INTEGER, ntorecv, 1, MPI_INTEGER, comm2d, ierr )
    IF ( ierr /= 0  .AND.  debug_output )  THEN
       WRITE( debug_string, * ) 'Error at MPI_ALLTOALL1:', ierr, ntosend, ntorecv
       CALL debug_message( debug_string, 'info' )
    ENDIF
!
!-- Calculate initial displacements
    i = 0
    j = 0
    DO  iproc = 0, numprocs-1
       dsend(iproc) = i
       dsendnow(iproc) = i
       drecv(iproc) = j
       drecvnow(iproc) = j
       i = i + ntosend(iproc)
       j = j + ntorecv(iproc)
    ENDDO
    dsend(numprocs) = i  ! Behind last pos = sum of all to send
    drecv(numprocs) = j  ! Behind last pos = sum of all to receive
    nrecv = j
!
!-- Allocate receive buffers
    ALLOCATE( irecvbuf(0:nrecv*npint-1) )
    ALLOCATE( frecvbuf(0:nrecv*npfloat-1) )
!
!-- Determine number of iterations among all processes
!-- (e.g. this process may have nothing to send and receive, yet some other still might)
    IF ( bufsize_alltoall <= 0 )  THEN
       nitersg = 1
       nmaxsend = HUGE( nitersg )
    ELSE
       nmaxsend = bufsize_alltoall
       niters = ( MAXVAL( ntosend(:) ) + nmaxsend - 1 ) / nmaxsend
       CALL MPI_ALLREDUCE( niters, nitersg, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
       IF ( nitersg > 1 )  THEN
          WRITE( debug_string, '("The MPI_ALLTOALL call has been split to ",I8," iterations ' //   &
                               'of max. ",I12," records each.")' ) nitersg, bufsize_alltoall
          CALL debug_message( debug_string, 'info' )
       ENDIF
    ENDIF
!
!-- Iterate ALLTOALLV using max-sized buffers
    DO  iter = 1, nitersg
       nsendnow(:) = MIN( dsend(1:) - dsendnow(:), nmaxsend )
       nrecvnow(:) = MIN( drecv(1:) - drecvnow(:), nmaxsend )
!
!--    Send integer data
       nsendnow_np(:) = nsendnow(:) * npint
       dsendnow_np(:) = dsendnow(:) * npint
       nrecvnow_np(:) = nrecvnow(:) * npint
       drecvnow_np(:) = drecvnow(:) * npint
       CALL MPI_ALLTOALLV( isendbuf, nsendnow_np, dsendnow_np, MPI_INTEGER,                        &
                           irecvbuf, nrecvnow_np, drecvnow_np, MPI_INTEGER, comm2d, ierr )

       IF ( ierr /= 0  .AND.  debug_output )  THEN
          WRITE( debug_string, * ) 'Error at MPI_ALLTOALLV 1:', ierr, iter, nmaxsend, dsend,       &
                                   dsendnow, nsendnow, drecv, drecvnow, nrecvnow
          CALL debug_message( debug_string, 'info' )
       ENDIF
!
!--    Send floating point data
       nsendnow_np(:) = nsendnow(:) * npfloat
       dsendnow_np(:) = dsendnow(:) * npfloat
       nrecvnow_np(:) = nrecvnow(:) * npfloat
       drecvnow_np(:) = drecvnow(:) * npfloat
       CALL MPI_ALLTOALLV( fsendbuf, nsendnow_np, dsendnow_np, MPI_REAL,                           &
                           frecvbuf, nrecvnow_np, drecvnow_np, MPI_REAL, comm2d, ierr )

       IF ( ierr /= 0  .AND.  debug_output )  THEN
          WRITE( debug_string, * ) 'Error at MPI_ALLTOALLV 2:', ierr, iter, nmaxsend, dsend,       &
                                   dsendnow, nsendnow, drecv, drecvnow, nrecvnow
          CALL debug_message( debug_string, 'info' )
       ENDIF
!
!--    Shift displacements for next iteration
       dsendnow(:) = dsendnow(:) + nsendnow(:)
       drecvnow(:) = drecvnow(:) + nrecvnow(:)
    ENDDO

    DEALLOCATE( drecv, drecvnow, drecvnow_np, dsend, dsendnow, dsendnow_np, nrecvnow, nrecvnow_np, &
                nsendnow, nsendnow_np, ntorecv )

 END SUBROUTINE radiation_exchange_alltoall
#endif


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!--------------------------------------------------------------------------------------------------!
SUBROUTINE radiation_3d_data_averaging( mode, variable )


    USE control_parameters

    USE indices

    USE kinds

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  variable  !<
    CHARACTER(LEN=*) ::  mode      !<

    CHARACTER(LEN=varnamelength) ::  var  !<

    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  imrt                !< index of MRT
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  l, m                !< index of current surface element
    INTEGER(iwp) ::  ids, idsint, isurf  !<

!
!-- Find the real name of the variable
    ids = -1
    l = -1
    var = TRIM( variable )
    DO  i = 0, nd-1
       k = LEN( TRIM( var ) )
       j = LEN( TRIM( dirname(i) ) )
       IF ( k - j + 1 >= 1 )  THEN
          IF ( TRIM( var(k-j+1:k) ) == TRIM( dirname(i) ) )  THEN
              ids = i
              idsint = dirint(ids)
              var = var(:k-j)
              EXIT
          ENDIF
       ENDIF
    ENDDO
    IF ( ids == -1 )  THEN
        var = TRIM( variable )
    ENDIF

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( var ) )
!
!--          Block of large scale (e.g. RRTMG) radiation output variables
             CASE ( 'rad_net*' )
                IF ( .NOT. ALLOCATED( rad_net_av ) )  THEN
                   ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_net_av = 0.0_wp

             CASE ( 'rad_lw_in*' )
                IF ( .NOT. ALLOCATED( rad_lw_in_xy_av ) )  THEN
                   ALLOCATE( rad_lw_in_xy_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_in_xy_av = 0.0_wp

             CASE ( 'rad_lw_out*' )
                IF ( .NOT. ALLOCATED( rad_lw_out_xy_av ) )  THEN
                   ALLOCATE( rad_lw_out_xy_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_out_xy_av = 0.0_wp

             CASE ( 'rad_sw_in*' )
                IF ( .NOT. ALLOCATED( rad_sw_in_xy_av ) )  THEN
                   ALLOCATE( rad_sw_in_xy_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_in_xy_av = 0.0_wp

             CASE ( 'rad_sw_out*' )
                IF ( .NOT. ALLOCATED( rad_sw_out_xy_av ) )  THEN
                   ALLOCATE( rad_sw_out_xy_av(nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_out_xy_av = 0.0_wp

             CASE ( 'rad_lw_in' )
                IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
                   ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_in_av = 0.0_wp

             CASE ( 'rad_lw_out' )
                IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
                   ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_out_av = 0.0_wp

             CASE ( 'rad_lw_cs_hr' )
                IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
                   ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_cs_hr_av = 0.0_wp

             CASE ( 'rad_lw_hr' )
                IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
                   ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_lw_hr_av = 0.0_wp

             CASE ( 'rad_sw_in' )
                IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
                   ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_in_av = 0.0_wp

             CASE ( 'rad_sw_out' )
                IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
                   ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_out_av = 0.0_wp

             CASE ( 'rad_sw_cs_hr' )
                IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
                   ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_cs_hr_av = 0.0_wp

             CASE ( 'rad_sw_hr' )
                IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
                   ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rad_sw_hr_av = 0.0_wp
!
!--          Block of RTM output variables
             CASE ( 'rtm_rad_net' )
!
!--              Array of complete radiation balance
                 IF ( .NOT.  ALLOCATED( surfradnet_av ) )  THEN
                     ALLOCATE( surfradnet_av(nsurfl) )
                 ENDIF
                 surfradnet_av = 0.0_wp

             CASE ( 'rtm_rad_insw' )
!
!--              Array of sw radiation falling to surface after i-th reflection
                 IF ( .NOT.  ALLOCATED( surfinsw_av ) )  THEN
                     ALLOCATE( surfinsw_av(nsurfl) )
                 ENDIF
                 surfinsw_av = 0.0_wp

             CASE ( 'rtm_rad_inlw' )
!
!--              Array of lw radiation falling to surface after i-th reflection
                 IF ( .NOT.  ALLOCATED( surfinlw_av ) )  THEN
                     ALLOCATE( surfinlw_av(nsurfl) )
                 ENDIF
                 surfinlw_av = 0.0_wp

             CASE ( 'rtm_rad_inswdir' )
!
!--              Array of direct sw radiation falling to surface from sun
                 IF ( .NOT.  ALLOCATED( surfinswdir_av ) )  THEN
                     ALLOCATE( surfinswdir_av(nsurfl) )
                 ENDIF
                 surfinswdir_av = 0.0_wp

             CASE ( 'rtm_rad_inswdif' )
!
!--              Array of difusion sw radiation falling to surface from sky and borders of the domain
                 IF ( .NOT.  ALLOCATED( surfinswdif_av ) )  THEN
                     ALLOCATE( surfinswdif_av(nsurfl) )
                 ENDIF
                 surfinswdif_av = 0.0_wp

             CASE ( 'rtm_rad_inswref' )
!
!--              Array of sw radiation falling to surface from reflections
                 IF ( .NOT.  ALLOCATED( surfinswref_av ) )  THEN
                     ALLOCATE( surfinswref_av(nsurfl) )
                 ENDIF
                 surfinswref_av = 0.0_wp

             CASE ( 'rtm_rad_inlwdif' )
!
!--             Array of sw radiation falling to surface after i-th reflection
                IF ( .NOT.  ALLOCATED( surfinlwdif_av ) )  THEN
                     ALLOCATE( surfinlwdif_av(nsurfl) )
                 ENDIF
                 surfinlwdif_av = 0.0_wp

             CASE ( 'rtm_rad_inlwref' )
!
!--              Array of lw radiation falling to surface from reflections
                 IF ( .NOT.  ALLOCATED( surfinlwref_av ) )  THEN
                     ALLOCATE( surfinlwref_av(nsurfl) )
                 ENDIF
                 surfinlwref_av = 0.0_wp

             CASE ( 'rtm_rad_outsw' )
!
!--              Array of sw radiation emitted from surface after i-th reflection
                 IF ( .NOT.  ALLOCATED( surfoutsw_av ) )  THEN
                     ALLOCATE( surfoutsw_av(nsurfl) )
                 ENDIF
                 surfoutsw_av = 0.0_wp

             CASE ( 'rtm_rad_outlw' )
!
!--              Array of lw radiation emitted from surface after i-th reflection
                 IF ( .NOT.  ALLOCATED( surfoutlw_av ) )  THEN
                     ALLOCATE( surfoutlw_av(nsurfl) )
                     surfoutlw_av = 0.0_wp
                 ENDIF
             CASE ( 'rtm_rad_ressw' )
!
!--              Array of residua of sw radiation absorbed in surface after last reflection
                 IF ( .NOT.  ALLOCATED( surfins_av ) )  THEN
                     ALLOCATE( surfins_av(nsurfl) )
                 ENDIF
                 surfins_av = 0.0_wp

             CASE ( 'rtm_rad_reslw' )
!
!--              Array of residua of lw radiation absorbed in surface after last reflection
                 IF ( .NOT.  ALLOCATED( surfinl_av ) )  THEN
                     ALLOCATE( surfinl_av(nsurfl) )
                 ENDIF
                 surfinl_av = 0.0_wp

             CASE ( 'rtm_rad_pc_inlw' )
!
!--              Array of of lw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED( pcbinlw_av ) )  THEN
                     ALLOCATE( pcbinlw_av(1:npcbl) )
                     pcbinlw_av = 0.0_wp
                 ENDIF

             CASE ( 'rtm_rad_pc_insw' )
!
!--              Array of of sw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED( pcbinsw_av ) )  THEN
                     ALLOCATE( pcbinsw_av(1:npcbl) )
                 ENDIF
                 pcbinsw_av = 0.0_wp

             CASE ( 'rtm_rad_pc_inswdir' )
!
!--              Array of of direct sw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED( pcbinswdir_av ) )  THEN
                     ALLOCATE( pcbinswdir_av(1:npcbl) )
                 ENDIF
                 pcbinswdir_av = 0.0_wp

             CASE ( 'rtm_rad_pc_inswdif' )
!
!--              Array of of diffuse sw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED( pcbinswdif_av ) )  THEN
                     ALLOCATE( pcbinswdif_av(1:npcbl) )
                 ENDIF
                 pcbinswdif_av = 0.0_wp

             CASE ( 'rtm_rad_pc_inswref' )
!
!--              Array of of reflected sw radiation absorbed in plant canopy
                 IF ( .NOT.  ALLOCATED( pcbinswref_av ) )  THEN
                     ALLOCATE( pcbinswref_av(1:npcbl) )
                 ENDIF
                 pcbinswref_av = 0.0_wp

             CASE ( 'rtm_rad_pc_sw_in' )
!
!--              Array of incoming sw radiation in plant canopy
                 IF ( .NOT.  ALLOCATED( pcinsw_av ) )  THEN
                     ALLOCATE( pcinsw_av(1:npcbl) )
                 ENDIF
                 pcinsw_av = 0.0_wp

             CASE ( 'rtm_rad_pc_sw_dir' )
!
!--              Array of incoming direct sw radiation in plant canopy
                 IF ( .NOT.  ALLOCATED( pcinswdir_av ) )  THEN
                     ALLOCATE( pcinswdir_av(1:npcbl) )
                 ENDIF
                 pcinswdir_av = 0.0_wp

             CASE ( 'rtm_rad_pc_sw_dif' )
!
!--              Array of incoming diffuse sw radiation in plant canopy
                 IF ( .NOT.  ALLOCATED( pcinswdif_av ) )  THEN
                     ALLOCATE( pcinswdif_av(1:npcbl) )
                 ENDIF
                 pcinswdif_av = 0.0_wp

             CASE ( 'rtm_mrt_sw' )
                IF ( .NOT. ALLOCATED( mrtinsw_av ) )  THEN
                   ALLOCATE( mrtinsw_av(nmrtbl) )
                ENDIF
                mrtinsw_av = 0.0_wp

             CASE ( 'rtm_mrt_lw' )
                IF ( .NOT. ALLOCATED( mrtinlw_av ) )  THEN
                   ALLOCATE( mrtinlw_av(nmrtbl) )
                ENDIF
                mrtinlw_av = 0.0_wp

             CASE ( 'rtm_mrt' )
                IF ( .NOT. ALLOCATED( mrt_av ) )  THEN
                   ALLOCATE( mrt_av(nmrtbl) )
                ENDIF
                mrt_av = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( var ) )
!
!--       Sum-up surface-related radiation quantities. Only the flux from the uppermost upward-
!--       facing surface is taken.
          CASE ( 'rad_net*' )
             IF ( ALLOCATED( rad_net_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            rad_net_av(j,i) = rad_net_av(j,i) + surf_lsm%rad_net(m)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            rad_net_av(j,i) = rad_net_av(j,i) + surf_usm%rad_net(m)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in*' )
             IF ( ALLOCATED( rad_lw_in_xy_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            rad_lw_in_xy_av(j,i) = rad_lw_in_xy_av(j,i) + surf_lsm%rad_lw_in(m)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            rad_lw_in_xy_av(j,i) = rad_lw_in_xy_av(j,i) + surf_usm%rad_lw_in(m)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out*' )
             IF ( ALLOCATED( rad_lw_out_xy_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            rad_lw_out_xy_av(j,i) = rad_lw_out_xy_av(j,i) + surf_lsm%rad_lw_out(m)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            rad_lw_out_xy_av(j,i) = rad_lw_out_xy_av(j,i) + surf_usm%rad_lw_out(m)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in*' )
             IF ( ALLOCATED( rad_sw_in_xy_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            rad_sw_in_xy_av(j,i) = rad_sw_in_xy_av(j,i) + surf_lsm%rad_sw_in(m)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            rad_sw_in_xy_av(j,i) = rad_sw_in_xy_av(j,i) + surf_usm%rad_sw_in(m)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out*' )
             IF ( ALLOCATED( rad_sw_out_xy_av ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn

                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         IF ( surf_lsm%upward_top(m) )  THEN
                            rad_sw_out_xy_av(j,i) = rad_sw_out_xy_av(j,i) + surf_lsm%rad_sw_out(m)
                         ENDIF
                      ENDDO

                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         IF ( surf_usm%upward_top(m) )  THEN
                            rad_sw_out_xy_av(j,i) = rad_sw_out_xy_av(j,i) + surf_usm%rad_sw_out(m)
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in' )
             IF ( ALLOCATED( rad_lw_in_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_in_av(k,j,i) = rad_lw_in_av(k,j,i) + rad_lw_in(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out' )
             IF ( ALLOCATED( rad_lw_out_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_out_av(k,j,i) = rad_lw_out_av(k,j,i) + rad_lw_out(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_cs_hr' )
             IF ( ALLOCATED( rad_lw_cs_hr_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_cs_hr_av(k,j,i) = rad_lw_cs_hr_av(k,j,i) + rad_lw_cs_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_hr' )
             IF ( ALLOCATED( rad_lw_hr_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_hr_av(k,j,i) = rad_lw_hr_av(k,j,i) + rad_lw_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in' )
             IF ( ALLOCATED( rad_sw_in_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_in_av(k,j,i) = rad_sw_in_av(k,j,i) + rad_sw_in(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out' )
             IF ( ALLOCATED( rad_sw_out_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_out_av(k,j,i) = rad_sw_out_av(k,j,i) + rad_sw_out(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_cs_hr' )
             IF ( ALLOCATED( rad_sw_cs_hr_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_cs_hr_av(k,j,i) = rad_sw_cs_hr_av(k,j,i) + rad_sw_cs_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_hr' )
             IF ( ALLOCATED( rad_sw_hr_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_hr_av(k,j,i) = rad_sw_hr_av(k,j,i) + rad_sw_hr(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--       Block of RTM output variables.
          CASE ( 'rtm_rad_net' )
!
!--           Array of complete radiation balance
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                    surfradnet_av(isurf) = surfradnet_av(isurf) + surfinsw(isurf) -                &
                                           surfoutsw(isurf) + surfinlw(isurf) - surfoutlw(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_insw' )
!
!--           Array of sw radiation falling to surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinsw_av(isurf) = surfinsw_av(isurf) + surfinsw(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inlw' )
!
!--           Array of lw radiation falling to surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinlw_av(isurf) = surfinlw_av(isurf) + surfinlw(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswdir' )
!
!--           Array of direct sw radiation falling to surface from sun
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinswdir_av(isurf) = surfinswdir_av(isurf) + surfinswdir(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswdif' )
!
!--           Array of diffusion sw radiation falling to surface from sky and borders of the domain
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinswdif_av(isurf) = surfinswdif_av(isurf) + surfinswdif(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswref' )
!
!--           Array of sw radiation falling to surface from reflections
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinswref_av(isurf) = surfinswref_av(isurf) + surfinsw(isurf) -             &
                                             surfinswdir(isurf) - surfinswdif(isurf)
                 ENDIF
              ENDDO


          CASE ( 'rtm_rad_inlwdif' )
!
!--           Array of sw radiation falling to surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinlwdif_av(isurf) = surfinlwdif_av(isurf) + surfinlwdif(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inlwref' )
!
!--           Array of lw radiation falling to surface from reflections
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinlwref_av(isurf) = surfinlwref_av(isurf) +                               &
                                             surfinlw(isurf) - surfinlwdif(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_outsw' )
!
!--           Array of sw radiation emitted from surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfoutsw_av(isurf) = surfoutsw_av(isurf) + surfoutsw(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_outlw' )
!
!--           Array of lw radiation emitted from surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfoutlw_av(isurf) = surfoutlw_av(isurf) + surfoutlw(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_ressw' )
!
!--           Array of residua of sw radiation absorbed in surface after last reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfins_av(isurf) = surfins_av(isurf) + surfins(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_reslw' )
!
!--           Array of residua of lw radiation absorbed in surface after last reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinl_av(isurf) = surfinl_av(isurf) + surfinl(isurf)
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_pc_inlw' )
              DO  l = 1, npcbl
                 pcbinlw_av(l) = pcbinlw_av(l) + pcbinlw(l)
              ENDDO

          CASE ( 'rtm_rad_pc_insw' )
              DO  l = 1, npcbl
                 pcbinsw_av(l) = pcbinsw_av(l) + pcbinsw(l)
              ENDDO

          CASE ( 'rtm_rad_pc_inswdir' )
              DO  l = 1, npcbl
                 pcbinswdir_av(l) = pcbinswdir_av(l) + pcbinswdir(l)
              ENDDO

          CASE ( 'rtm_rad_pc_inswdif' )
              DO  l = 1, npcbl
                 pcbinswdif_av(l) = pcbinswdif_av(l) + pcbinswdif(l)
              ENDDO

          CASE ( 'rtm_rad_pc_inswref' )
              DO  l = 1, npcbl
                 pcbinswref_av(l) = pcbinswref_av(l) + pcbinsw(l) - pcbinswdir(l) - pcbinswdif(l)
              ENDDO

          CASE ( 'rtm_rad_pc_sw_in' )
              DO  l = 1, npcbl
                 pcinsw_av(l) = pcinsw_av(l) + pcinsw(l)
              ENDDO

          CASE ( 'rtm_rad_pc_sw_dir' )
              DO  l = 1, npcbl
                 pcinswdir_av(l) = pcinswdir_av(l) + pcinswdir(l)
              ENDDO

          CASE ( 'rtm_rad_pc_sw_dif' )
              DO  l = 1, npcbl
                 pcinswdif_av(l) = pcinswdif_av(l) + pcinswdif(l)
              ENDDO

          CASE ( 'rtm_mrt_sw' )
             IF ( ALLOCATED( mrtinsw_av ) )  THEN
                mrtinsw_av(:) = mrtinsw_av(:) + mrtinsw(:)
             ENDIF

          CASE ( 'rtm_mrt_lw' )
             IF ( ALLOCATED( mrtinlw_av ) )  THEN
                mrtinlw_av(:) = mrtinlw_av(:) + mrtinlw(:)
             ENDIF

          CASE ( 'rtm_mrt' )
             IF ( ALLOCATED( mrt_av ) )  THEN
                mrt_av(:) = mrt_av(:) + mrt(:)
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( var ) )
!
!--       Block of large scale (e.g. RRTMG) radiation output variables
          CASE ( 'rad_net*' )
             IF ( ALLOCATED( rad_net_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_net_av(j,i) = rad_net_av(j,i) / REAL( average_count_3d, KIND = wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in*' )
             IF ( ALLOCATED( rad_lw_in_xy_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_lw_in_xy_av(j,i) = rad_lw_in_xy_av(j,i) /                                &
                                             REAL( average_count_3d, KIND = wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out*' )
             IF ( ALLOCATED( rad_lw_out_xy_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_lw_out_xy_av(j,i) = rad_lw_out_xy_av(j,i) /                              &
                                              REAL( average_count_3d, KIND = wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in*' )
             IF ( ALLOCATED( rad_sw_in_xy_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_sw_in_xy_av(j,i) = rad_sw_in_xy_av(j,i) /                                &
                                             REAL( average_count_3d, KIND = wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out*' )
             IF ( ALLOCATED( rad_sw_out_xy_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      rad_sw_out_xy_av(j,i) = rad_sw_out_xy_av(j,i) /                              &
                                              REAL( average_count_3d, KIND = wp )
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_in' )
             IF ( ALLOCATED( rad_lw_in_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_in_av(k,j,i) = rad_lw_in_av(k,j,i) /                               &
                                               REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_out' )
             IF ( ALLOCATED( rad_lw_out_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_out_av(k,j,i) = rad_lw_out_av(k,j,i) /                             &
                                                REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_cs_hr' )
             IF ( ALLOCATED( rad_lw_cs_hr_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_cs_hr_av(k,j,i) = rad_lw_cs_hr_av(k,j,i) /                         &
                                                  REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_lw_hr' )
             IF ( ALLOCATED( rad_lw_hr_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_lw_hr_av(k,j,i) = rad_lw_hr_av(k,j,i) /                               &
                                               REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_in' )
             IF ( ALLOCATED( rad_sw_in_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_in_av(k,j,i) = rad_sw_in_av(k,j,i) /                               &
                                               REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_out' )
             IF ( ALLOCATED( rad_sw_out_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_out_av(k,j,i) = rad_sw_out_av(k,j,i) /                             &
                                                REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_cs_hr' )
             IF ( ALLOCATED( rad_sw_cs_hr_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_cs_hr_av(k,j,i) = rad_sw_cs_hr_av(k,j,i) /                         &
                                                  REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'rad_sw_hr' )
             IF ( ALLOCATED( rad_sw_hr_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rad_sw_hr_av(k,j,i) = rad_sw_hr_av(k,j,i) /                               &
                                               REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
!
!--       Block of RTM output variables
          CASE ( 'rtm_rad_net' )
!
!--           Array of complete radiation balance
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfradnet_av(isurf) = surfradnet_av(isurf) /                                 &
                                            REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_insw' )
!
!--           Array of sw radiation falling to surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinsw_av(isurf) = surfinsw_av(isurf) / REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inlw' )
!
!--           Array of lw radiation falling to surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinlw_av(isurf) = surfinlw_av(isurf) / REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswdir' )
!
!--           Array of direct sw radiation falling to surface from sun
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinswdir_av(isurf) = surfinswdir_av(isurf) /                               &
                                             REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswdif' )
!
!--           Array of diffusion sw radiation falling to surface from sky and borders of the domain
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinswdif_av(isurf) = surfinswdif_av(isurf) /                               &
                                             REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inswref' )
!
!--           Array of sw radiation falling to surface from reflections
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinswref_av(isurf) = surfinswref_av(isurf) /                               &
                                             REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inlwdif' )
!
!--           Array of sw radiation falling to surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinlwdif_av(isurf) = surfinlwdif_av(isurf) /                               &
                                             REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_inlwref' )
!
!--           Array of lw radiation falling to surface from reflections
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinlwref_av(isurf) = surfinlwref_av(isurf) /                               &
                                             REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_outsw' )
!
!--           Array of sw radiation emitted from surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfoutsw_av(isurf) = surfoutsw_av(isurf) / REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_outlw' )
!
!--           Array of lw radiation emitted from surface after i-th reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfoutlw_av(isurf) = surfoutlw_av(isurf) / REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_ressw' )
!
!--           Array of residua of sw radiation absorbed in surface after last reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfins_av(isurf) = surfins_av(isurf) / REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_reslw' )
!
!--           Array of residua of lw radiation absorbed in surface after last reflection
              DO  isurf = 1, nsurfl
                 IF ( surfl(id,isurf) == idsint )  THEN
                     surfinl_av(isurf) = surfinl_av(isurf) / REAL( average_count_3d, KIND = wp )
                 ENDIF
              ENDDO

          CASE ( 'rtm_rad_pc_inlw' )
              DO  l = 1, npcbl
                 pcbinlw_av(l) = pcbinlw_av(l) / REAL( average_count_3d, KIND = wp )
              ENDDO

          CASE ( 'rtm_rad_pc_insw' )
              DO  l = 1, npcbl
                 pcbinsw_av(l) = pcbinsw_av(l) / REAL( average_count_3d, KIND = wp )
              ENDDO

          CASE ( 'rtm_rad_pc_inswdir' )
              DO  l = 1, npcbl
                 pcbinswdir_av(l) = pcbinswdir_av(l) / REAL( average_count_3d, KIND = wp )
              ENDDO

          CASE ( 'rtm_rad_pc_inswdif' )
              DO  l = 1, npcbl
                 pcbinswdif_av(l) = pcbinswdif_av(l) / REAL( average_count_3d, KIND = wp )
              ENDDO

          CASE ( 'rtm_rad_pc_inswref' )
              DO  l = 1, npcbl
                 pcbinswref_av(l) = pcbinswref_av(l) / REAL( average_count_3d, KIND = wp )
              ENDDO

          CASE ( 'rtm_rad_pc_sw_in' )
              DO  l = 1, npcbl
                 pcinsw_av(l) = pcinsw_av(l) / REAL( average_count_3d, KIND = wp )
              ENDDO

          CASE ( 'rtm_rad_pc_sw_dir' )
              DO  l = 1, npcbl
                 pcinswdir_av(l) = pcinswdir_av(l) / REAL( average_count_3d, KIND = wp )
              ENDDO

          CASE ( 'rtm_rad_pc_sw_dif' )
              DO  l = 1, npcbl
                 pcinswdif_av(l) = pcinswdif_av(l) / REAL( average_count_3d, KIND = wp )
              ENDDO

          CASE ( 'rtm_mrt_sw' )
             IF ( ALLOCATED( mrtinsw_av ) )  THEN
                DO  imrt = 1, nmrtbl
                   mrtinsw_av(imrt) = mrtinsw_av(imrt) / REAL( average_count_3d, KIND = wp )
                ENDDO
             ENDIF

          CASE ( 'rtm_mrt_lw' )
             IF (  ALLOCATED( mrtinlw_av ) )  THEN
                DO  imrt = 1, nmrtbl
                   mrtinlw_av(imrt) = mrtinlw_av(imrt) / REAL( average_count_3d, KIND = wp )
                ENDDO
             ENDIF

          CASE ( 'rtm_mrt' )
             IF ( ALLOCATED( mrt_av ) )  THEN
                DO  imrt = 1, nmrtbl
                   mrt_av(imrt) = mrt_av(imrt) / REAL( average_count_3d, KIND = wp )
                ENDDO
             ENDIF

       END SELECT

    ENDIF

END SUBROUTINE radiation_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific averaging of surface data
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_surface_data_averaging( trimvar, n_out )

    CHARACTER(LEN=*), INTENT(IN) ::  trimvar  !< dummy variable for current output variable

    INTEGER(iwp), INTENT(IN) ::  n_out  !< counter variables for surface output


!
!-- So far we have no averaged variables, so just silence compiler waring about unused parameter
    IF ( trimvar == ''  .AND.  n_out == 0 )  CONTINUE

 END SUBROUTINE radiation_surface_data_averaging


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE radiation_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  variable  !<

    CHARACTER(LEN=*), INTENT(OUT) ::  grid_x  !<
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_y  !<
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_z  !<

    CHARACTER(LEN=varnamelength) ::  var  !<

    LOGICAL, INTENT(OUT) ::  found  !<

    found  = .TRUE.

!
!-- Check for the grid
    var = TRIM( variable )
!
!-- RTM directional variables
    IF ( var(1:12) == 'rtm_rad_net_'      .OR.  var(1:13) == 'rtm_rad_insw_'     .OR.              &
         var(1:13) == 'rtm_rad_inlw_'     .OR.  var(1:16) == 'rtm_rad_inswdir_'  .OR.              &
         var(1:16) == 'rtm_rad_inswdif_'  .OR.  var(1:16) == 'rtm_rad_inswref_'  .OR.              &
         var(1:16) == 'rtm_rad_inlwdif_'  .OR.  var(1:16) == 'rtm_rad_inlwref_'  .OR.              &
         var(1:14) == 'rtm_rad_outsw_'    .OR.  var(1:14) == 'rtm_rad_outlw_'    .OR.              &
         var(1:14) == 'rtm_rad_ressw_'    .OR.  var(1:14) == 'rtm_rad_reslw_'    .OR.              &
         var == 'rtm_rad_pc_inlw'         .OR.  var == 'rtm_rad_pc_insw'         .OR.              &
         var == 'rtm_rad_pc_inswdir'      .OR.  var == 'rtm_rad_pc_inswdif'      .OR.              &
         var == 'rtm_rad_pc_inswref'      .OR.  var(1:7) == 'rtm_svf'            .OR.              &
         var == 'rtm_rad_pc_sw_in'        .OR.  var == 'rtm_rad_pc_sw_dir'       .OR.              &
         var == 'rtm_rad_pc_sw_dif'       .OR.                                                     &
         var(1:7) == 'rtm_dif'            .OR.  var(1:12) == 'rtm_surfalb_'      .OR.              &
         var(1:13) == 'rtm_surfemis_'     .OR.  var == 'rtm_mrt'                 .OR.              &
         var == 'rtm_mrt_sw'              .OR.  var == 'rtm_mrt_lw'              .OR.              &
         var == 'rtm_rad_vol_sw' )  THEN

         found = .TRUE.
         grid_x = 'x'
         grid_y = 'y'
         grid_z = 'zu'
    ELSE

       SELECT CASE ( TRIM( var ) )

          CASE ( 'rad_lw_cs_hr', 'rad_lw_hr', 'rad_sw_cs_hr', 'rad_sw_hr', 'rad_lw_cs_hr_xy',      &
                 'rad_lw_hr_xy', 'rad_sw_cs_hr_xy', 'rad_sw_hr_xy', 'rad_lw_cs_hr_xz',             &
                 'rad_lw_hr_xz', 'rad_sw_cs_hr_xz', 'rad_sw_hr_xz', 'rad_lw_cs_hr_yz',             &
                 'rad_lw_hr_yz', 'rad_sw_cs_hr_yz', 'rad_sw_hr_yz' )
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zu'

          CASE ( 'rad_lw_in', 'rad_lw_out', 'rad_sw_in', 'rad_sw_out', 'rad_lw_in_xy',             &
                 'rad_lw_out_xy', 'rad_sw_in_xy','rad_sw_out_xy','rad_lw_in_xz','rad_lw_out_xz',   &
                 'rad_sw_in_xz','rad_sw_out_xz', 'rad_lw_in_yz', 'rad_lw_out_yz', 'rad_sw_in_yz',  &
                 'rad_sw_out_yz' )
             grid_x = 'x'
             grid_y = 'y'
             grid_z = 'zw'


          CASE DEFAULT
             found  = .FALSE.
             grid_x = 'none'
             grid_y = 'none'
             grid_z = 'none'

           END SELECT
       ENDIF

 END SUBROUTINE radiation_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 2D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do,    &
                                      nzt_do )

    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER(LEN=*) ::  grid      !<
    CHARACTER(LEN=*) ::  mode      !<
    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  av      !<
    INTEGER(iwp) ::  i       !<
    INTEGER(iwp) ::  j       !<
    INTEGER(iwp) ::  k       !<
    INTEGER(iwp) ::  m       !< index of surface element at grid point (j,i)
    INTEGER(iwp) ::  nzb_do  !<
    INTEGER(iwp) ::  nzt_do  !<

    LOGICAL ::  found  !<
    LOGICAL ::  two_d  !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !<

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'rad_net*_xy' )  ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type.
!--                Only upward faced horizontal outputs are considered here.
!--                Natural-type surfaces.
                   DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_lsm%rad_net(m), local_pf(i,j,nzb+1),       &
                                                   surf_lsm%upward(m) )
                   ENDDO
!
!--                Urban-type surfaces.
                   DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_usm%rad_net(m), local_pf(i,j,nzb+1),       &
                                                   surf_usm%upward(m) )
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_net_av ) )  THEN
                ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
                rad_net_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_net_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_lw_in*_xy' )  ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type.
!--                Natural-type surfaces.
                   DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_lsm%rad_lw_in(m), local_pf(i,j,nzb+1),     &
                                                   surf_lsm%upward(m) )
                   ENDDO
!
!--                Urban-type surfaces.
                   DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_usm%rad_lw_in(m), local_pf(i,j,nzb+1),     &
                                                   surf_usm%upward(m) )
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_in_xy_av ) )  THEN
                ALLOCATE( rad_lw_in_xy_av(nysg:nyng,nxlg:nxrg) )
                rad_lw_in_xy_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_lw_in_xy_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_lw_out*_xy' )  ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type.
!--                Natural-type surfaces.
                   DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_lsm%rad_lw_out(m), local_pf(i,j,nzb+1),    &
                                                   surf_lsm%upward(m) )
                   ENDDO
!
!--                Urban-type surfaces.
                   DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_usm%rad_lw_out(m), local_pf(i,j,nzb+1),    &
                                                   surf_usm%upward(m) )
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_out_xy_av ) )  THEN
                ALLOCATE( rad_lw_out_xy_av(nysg:nyng,nxlg:nxrg) )
                rad_lw_out_xy_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_lw_out_xy_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_sw_in*_xy' )  ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type.
!--                Natural-type surfaces.
                   DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_lsm%rad_sw_in(m), local_pf(i,j,nzb+1),     &
                                                   surf_lsm%upward(m) )
                   ENDDO
!
!--                Urban-type surfaces.
                   DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_usm%rad_sw_in(m), local_pf(i,j,nzb+1),     &
                                                   surf_usm%upward(m) )
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_in_xy_av ) )  THEN
                ALLOCATE( rad_sw_in_xy_av(nysg:nyng,nxlg:nxrg) )
                rad_sw_in_xy_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_sw_in_xy_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_sw_out*_xy' )  ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
!
!--                Obtain rad_net from its respective surface type.
!--                Natural-type surfaces.
                   DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_lsm%rad_sw_out(m), local_pf(i,j,nzb+1),     &
                                                   surf_lsm%upward(m) )
                   ENDDO
!
!--                Urban-type surfaces.
                   DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                      local_pf(i,j,nzb+1) = MERGE( surf_usm%rad_sw_out(m), local_pf(i,j,nzb+1),     &
                                                   surf_usm%upward(m) )
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_out_xy_av ) )  THEN
                ALLOCATE( rad_sw_out_xy_av(nysg:nyng,nxlg:nxrg) )
                rad_sw_out_xy_av = 0.0_wp
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = rad_sw_out_xy_av(j,i)
                ENDDO
             ENDDO
          ENDIF
          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'rad_lw_in_xy', 'rad_lw_in_xz', 'rad_lw_in_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_lw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
               ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_in_av = 0.0_wp
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_in_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_lw_out_xy', 'rad_lw_out_xz', 'rad_lw_out_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_lw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
               ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_out_av = 0.0_wp
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_out_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_lw_cs_hr_xy', 'rad_lw_cs_hr_xz', 'rad_lw_cs_hr_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_lw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
               ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_cs_hr_av = 0.0_wp
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_cs_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_lw_hr_xy', 'rad_lw_hr_xz', 'rad_lw_hr_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_lw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
               ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_hr_av= 0.0_wp
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_lw_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_sw_in_xy', 'rad_sw_in_xz', 'rad_sw_in_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_sw_in(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
               ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_in_av = 0.0_wp
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_in_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_sw_out_xy', 'rad_sw_out_xz', 'rad_sw_out_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_sw_out(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
               ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_out_av = 0.0_wp
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb, nzt+1
                      local_pf(i,j,k) = rad_sw_out_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zu'

       CASE ( 'rad_sw_cs_hr_xy', 'rad_sw_cs_hr_xz', 'rad_sw_cs_hr_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_sw_cs_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
               ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_cs_hr_av = 0.0_wp
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_cs_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE ( 'rad_sw_hr_xy', 'rad_sw_hr_xz', 'rad_sw_hr_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_sw_hr(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
            IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
               ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_hr_av = 0.0_wp
            ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = rad_sw_hr_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          IF ( mode == 'xy' )  grid = 'zw'

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT

 END SUBROUTINE radiation_data_output_2d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER(LEN=*) ::  variable !<

    CHARACTER(LEN=varnamelength) ::  var, surfid  !<

    INTEGER(iwp) ::  av                                          !<
    INTEGER(iwp) ::  i, j, k, l                                  !<
    INTEGER(iwp) ::  nzb_do                                      !<
    INTEGER(iwp) ::  nzt_do                                      !<
    INTEGER(iwp) ::  ids,idsint,isurf,isvf,isurfs,isurflt,ipcgb  !<
    INTEGER(iwp) ::  is, js, ks, istat                           !<

    LOGICAL ::  found  !<

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !<


    found = .TRUE.
    var = TRIM( variable )
!
!-- Check if variable belongs to radiation related variables (starts with rad or rtm)
    IF ( LEN( var ) < 3  )  THEN
       found = .FALSE.
       RETURN
    ENDIF

    IF ( var(1:3) /= 'rad'  .AND.  var(1:3) /= 'rtm' )  THEN
       found = .FALSE.
       RETURN
    ENDIF

    ids = -1
    DO  i = 0, nd-1
       k = LEN( TRIM( var ) )
       j = LEN( TRIM( dirname(i) ) )
       IF ( k - j + 1 >= 1 )  THEN
          IF ( TRIM( var(k-j+1:k) ) == TRIM( dirname(i) ) )  THEN
             ids = i
             idsint = dirint(ids)
             var = var(:k-j)
             EXIT
          ENDIF
       ENDIF
    ENDDO
    IF ( ids == -1 )  THEN
        var = TRIM( variable )
    ENDIF

    IF ( (var(1:8) == 'rtm_svf_'  .OR.  var(1:8) == 'rtm_dif_')  .AND.  LEN( TRIM( var ) ) >= 13 ) &
    THEN
!
!--     svf values to particular surface
        surfid = var(9:)
        i = INDEX( surfid, '_' )
        j = INDEX( surfid(i+1:), '_' )
        READ( surfid(1:i-1), *, IOSTAT = istat ) is
        IF ( istat == 0 )  THEN
            READ( surfid(i+1:i+j-1), *, IOSTAT = istat ) js
        ENDIF
        IF ( istat == 0 )  THEN
            READ( surfid(i+j+1:), *, IOSTAT = istat ) ks
        ENDIF
        IF ( istat == 0 )  THEN
            var = var(1:7)
        ENDIF
    ENDIF

    SELECT CASE ( TRIM( var ) )
!
!--   Block of large scale radiation model (e.g. RRTMG) output variables.
      CASE ( 'rad_sw_in' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_sw_in(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
               ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_in_av = 0.0_wp
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_in_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_out' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_sw_out(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
               ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_out_av = 0.0_wp
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_out_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_cs_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_sw_cs_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
               ALLOCATE( rad_sw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_cs_hr_av = 0.0_wp
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_cs_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_sw_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_sw_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
               ALLOCATE( rad_sw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_sw_hr_av = 0.0_wp
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_sw_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_in' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_lw_in(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
               ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_in_av = 0.0_wp
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_in_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_out' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_lw_out(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
               ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_out_av = 0.0_wp
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_out_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_cs_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_lw_cs_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
               ALLOCATE( rad_lw_cs_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
               rad_lw_cs_hr_av = 0.0_wp
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_cs_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rad_lw_hr' )
         IF ( av == 0 )  THEN
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     IF ( BTEST( topo_flags(k,j,i), 0 ) )  local_pf(i,j,k) = rad_lw_hr(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
               ALLOCATE( rad_lw_hr_av(nzb+1:nzt+1,nysg:nyng,nxlg:nxrg) )
              rad_lw_hr_av = 0.0_wp
            ENDIF
            DO  i = nxl, nxr
               DO  j = nys, nyn
                  DO  k = nzb_do, nzt_do
                     local_pf(i,j,k) = rad_lw_hr_av(k,j,i)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      CASE ( 'rtm_rad_net' )
!
!--      Array of complete radiation balance
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) =                      &
                  surfinsw(isurf) - surfoutsw(isurf) +  surfinlw(isurf) - surfoutlw(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfradnet_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_insw' )
!
!--      Array of sw radiation falling to surface after i-th reflection
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                 local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinsw(isurf)
               ELSE
                 local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinsw_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inlw' )
!
!--      Array of lw radiation falling to surface after i-th reflection
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlw(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlw_av(isurf)
               ENDIF
             ENDIF
         ENDDO

      CASE ( 'rtm_rad_inswdir' )
!
!--      Array of direct sw radiation falling to surface from sun
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswdir(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswdir_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inswdif' )
!
!--      Array of difusion sw radiation falling to surface from sky and borders of the domain
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswdif(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswdif_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inswref' )
!
!--      Array of sw radiation falling to surface from reflections
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) =                      &
                  surfinsw(isurf) - surfinswdir(isurf) - surfinswdif(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinswref_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inlwdif' )
!
!--      Array of difusion lw radiation falling to surface from sky and borders of the domain
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlwdif(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlwdif_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_inlwref' )
!
!--      Array of lw radiation falling to surface from reflections
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlw(isurf) -    &
                                                                              surfinlwdif(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinlwref_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_outsw' )
!
!--      Array of sw radiation emitted from surface after i-th reflection
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfoutsw(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfoutsw_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_outlw' )
!
!--      Array of lw radiation emitted from surface after i-th reflection
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfoutlw(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfoutlw_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_ressw' )
!
!--      Average of array of residua of sw radiation absorbed in surface after last reflection
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfins(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfins_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_reslw' )
!
!--      Average of array of residua of lw radiation absorbed in surface after last reflection
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               IF ( av == 0 )  THEN
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinl(isurf)
               ELSE
                  local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = surfinl_av(isurf)
               ENDIF
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_inlw' )
!
!--      Array of lw radiation absorbed by plant canopy
         DO  ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinlw(ipcgb)
            ELSE
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinlw_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_insw' )
!
!--      Array of sw radiation absorbed by plant canopy
         DO  ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinsw(ipcgb)
            ELSE
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinsw_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_inswdir' )
!
!--      Array of direct sw radiation absorbed by plant canopy
         DO  ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswdir(ipcgb)
            ELSE
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswdir_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_inswdif' )
!
!--      Array of diffuse sw radiation absorbed by plant canopy
         DO  ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswdif(ipcgb)
            ELSE
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswdif_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_inswref' )
!
!--      Array of reflected sw radiation absorbed by plant canopy
         DO  ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) =                            &
               pcbinsw(ipcgb) - pcbinswdir(ipcgb) - pcbinswdif(ipcgb)
            ELSE
               local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcbinswref_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_sw_in' )
!
!--      Array of incoming sw radiation to plant canopy
         DO  ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcinsw(ipcgb)
            ELSE
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcinsw_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_sw_dir' )
!
!--      Array of direct incoming sw radiation to plant canopy
         DO  ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcinswdir(ipcgb)
            ELSE
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcinswdir_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_pc_sw_dif' )
!
!--      Array of diffuse incoming sw radiation to plant canopy
         DO  ipcgb = 1, npcbl
            IF ( av == 0 )  THEN
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcinswdif(ipcgb)
            ELSE
              local_pf(pcbl(ix,ipcgb),pcbl(iy,ipcgb),pcbl(iz,ipcgb)) = pcinswdif_av(ipcgb)
            ENDIF
         ENDDO

      CASE ( 'rtm_rad_vol_sw' )
!
!--      3-d volumetric output of shortwave radiative flux density
          IF ( av == 0 .AND. radiation_volumetric_flux )  THEN
            local_pf(nxl:nxr,nys:nyn,nz_urban_b:nz_urban_t) =                    &
                 RESHAPE( swflux_vol(nz_urban_b:nz_urban_t,nys:nyn,nxl:nxr),     &
                          (/ nxr-nxl+1, nyn-nys+1, nz_urban_t-nz_urban_b+1 /),   &
                          ORDER = (/ 3,2,1 /) )
         ENDIF

      CASE ( 'rtm_mrt_sw' )
         local_pf = 0.0_wp
         IF ( av == 0 )  THEN
            DO  l = 1, nmrtbl
               local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrtinsw(l)
            ENDDO
         ELSE
            IF ( ALLOCATED( mrtinsw_av ) )  THEN
               DO  l = 1, nmrtbl
                  local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrtinsw_av(l)
               ENDDO
            ENDIF
         ENDIF

      CASE ( 'rtm_mrt_lw' )
         local_pf = 0.0_wp
         IF ( av == 0 )  THEN
            DO  l = 1, nmrtbl
               local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrtinlw(l)
            ENDDO
         ELSE
            IF ( ALLOCATED( mrtinlw_av ) )  THEN
               DO  l = 1, nmrtbl
                  local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrtinlw_av(l)
               ENDDO
            ENDIF
         ENDIF

      CASE ( 'rtm_mrt' )
         local_pf = 0.0_wp
         IF ( av == 0 )  THEN
            DO  l = 1, nmrtbl
               local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrt(l)
            ENDDO
         ELSE
            IF ( ALLOCATED( mrt_av ) )  THEN
               DO  l = 1, nmrtbl
                  local_pf(mrtbl(ix,l),mrtbl(iy,l),mrtbl(iz,l)) = mrt_av(l)
               ENDDO
            ENDIF
         ENDIF

      CASE ( 'rtm_svf', 'rtm_dif' )
!
!--      Shape view factors or iradiance factors to selected surface
         IF ( TRIM( var ) == 'rtm_svf' )  THEN
             k = 1
         ELSE
             k = 2
         ENDIF
         DO  isvf = 1, nsvfl
            isurflt = svfsurf(1, isvf)
            isurfs = svfsurf(2, isvf)

            IF ( surf(ix,isurfs) == is  .AND.  surf(iy,isurfs) == js  .AND. surf(iz,isurfs) == ks  &
                 .AND.  surfl(id,isurflt) == idsint )  THEN
!
!--            Correct source surface
               local_pf(surfl(ix,isurflt),surfl(iy,isurflt),surfl(iz,isurflt)) = svf(k,isvf)
            ENDIF
         ENDDO

      CASE ( 'rtm_surfalb' )
!
!--      Surface albedo
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = albedo_surf(isurf)
            ENDIF
         ENDDO

      CASE ( 'rtm_surfemis' )
!
!--      Surface emissivity, weighted average
         DO  isurf = 1, nsurfl
            IF ( surfl(id,isurf) == idsint )  THEN
               local_pf(surfl(ix,isurf),surfl(iy,isurf),surfl(iz,isurf)) = emiss_surf(isurf)
            ENDIF
         ENDDO

      CASE DEFAULT
         found = .FALSE.

    END SELECT

 END SUBROUTINE radiation_data_output_3d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining masked data output
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_mask( av, variable, found, local_pf, mid )

    USE control_parameters

    USE indices

    USE kinds


    IMPLICIT NONE

    CHARACTER (LEN=*) ::  variable  !<

    CHARACTER(LEN=5) ::  grid  !< flag to distinquish between staggered grids

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


       CASE ( 'rad_lw_in' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_lw_in
          ELSE
             to_be_resorted => rad_lw_in_av
          ENDIF

       CASE ( 'rad_lw_out' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_lw_out
          ELSE
             to_be_resorted => rad_lw_out_av
          ENDIF

       CASE ( 'rad_lw_cs_hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_lw_cs_hr
          ELSE
             to_be_resorted => rad_lw_cs_hr_av
          ENDIF

       CASE ( 'rad_lw_hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_lw_hr
          ELSE
             to_be_resorted => rad_lw_hr_av
          ENDIF

       CASE ( 'rad_sw_in' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_sw_in
          ELSE
             to_be_resorted => rad_sw_in_av
          ENDIF

       CASE ( 'rad_sw_out' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_sw_out
          ELSE
             to_be_resorted => rad_sw_out_av
          ENDIF

       CASE ( 'rad_sw_cs_hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_sw_cs_hr
          ELSE
             to_be_resorted => rad_sw_cs_hr_av
          ENDIF

       CASE ( 'rad_sw_hr' )
          IF ( av == 0 )  THEN
             to_be_resorted => rad_sw_hr
          ELSE
             to_be_resorted => rad_sw_hr_av
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
                topo_top_index = topo_top_ind(mask_j(mid,j), mask_i(mid,i), 0 )
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



 END SUBROUTINE radiation_data_output_mask


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Define radiation surface output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_data_output_surf( av, trimvar, found )

    CHARACTER (LEN=*), INTENT(IN)    ::  trimvar    !< variable name

    INTEGER(iwp), INTENT(IN) ::  av  !< flag for (non-)average output

    INTEGER(iwp) ::  offset_lsm  !< offset for LSM surfaces in surf_out
    INTEGER(iwp) ::  offset_usm  !< offset for USM surfaces in surf_out
    INTEGER(iwp) ::  i           !< horizontal coordinate
    INTEGER(iwp) ::  j           !< horizontal coordinate
    INTEGER(iwp) ::  iso         !< running index for surf_out elements
    INTEGER(iwp) ::  isurf_rtm   !< running index for RTM surface elements

    LOGICAL, INTENT(INOUT) ::  found  !< flag if output variable is found


!
!-- The code in the cycles depends on the order of the execution. Do not parallelize by OpenMP!
!-- Surfaces in surf_out are (by definition) ordered by simply stacking surf_def followed by
!-- surf_lsm and then surf_usm. Surfaces in RTM are ordered by i, j and within that as in surf_usm
!-- followed by surf_lsm.
    offset_lsm = surf_def%ns ! should be zero if RTM is enabled
    offset_usm = offset_lsm + surf_lsm%ns

    found = .TRUE.
    isurf_rtm = 1

    SELECT CASE ( TRIM( trimvar ) )

       CASE ( 'rtm_skyvf' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  iso = surf_usm%start_index(j,i)+offset_usm, surf_usm%end_index(j,i)+offset_usm
                   surf_out%var_out(iso) = skyvf(isurf_rtm)
                   isurf_rtm = isurf_rtm + 1
                ENDDO
                DO  iso = surf_lsm%start_index(j,i)+offset_lsm, surf_lsm%end_index(j,i)+offset_lsm
                   surf_out%var_out(iso) = skyvf(isurf_rtm)
                   isurf_rtm = isurf_rtm + 1
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'rtm_skyvft' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  iso = surf_usm%start_index(j,i)+offset_usm, surf_usm%end_index(j,i)+offset_usm
                   surf_out%var_out(iso) = skyvft(isurf_rtm)
                   isurf_rtm = isurf_rtm + 1
                ENDDO
                DO  iso = surf_lsm%start_index(j,i)+offset_lsm, surf_lsm%end_index(j,i)+offset_lsm
                   surf_out%var_out(iso) = skyvft(isurf_rtm)
                   isurf_rtm = isurf_rtm + 1
                ENDDO
             ENDDO
          ENDDO

       CASE DEFAULT
          found = .FALSE.

    END SELECT
!
!-- So far we have no averaged variables, so just silence compiler waring about unused parameter.
    IF ( av == 0 )  CONTINUE

 END SUBROUTINE radiation_data_output_surf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate module specific statistics for radiation model, i.e. timeseries. Profiles could be
!> added in the future.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_statistics( mode, sr )

    CHARACTER (LEN=*) ::  mode  !< output mode

    INTEGER(iwp) ::  sr    !< number of statistic region


    IF ( mode == 'time_series' )  THEN
!
!--    Store time series date.
       ts_value(dots_start_index_rtm,sr)   = hom(nzb,1,99,sr)    ! rad_net
       ts_value(dots_start_index_rtm+1,sr) = hom(nzb,1,100,sr)   ! rad_lw_in
       ts_value(dots_start_index_rtm+2,sr) = hom(nzb,1,101,sr)   ! rad_lw_out
       ts_value(dots_start_index_rtm+3,sr) = hom(nzb,1,102,sr)   ! rad_sw_in
       ts_value(dots_start_index_rtm+4,sr) = hom(nzb,1,103,sr)   ! rad_sw_out

       IF ( radiation_scheme /= 'tenstream' )  THEN

          IF ( average_radiation ) THEN
             ts_value(dots_start_index_rtm+5,sr) = t_rad_eff
             ts_value(dots_start_index_rtm+6,sr) = emissivity_eff
             ts_value(dots_start_index_rtm+7,sr) = albedo_eff
          ENDIF

          IF ( radiation_scheme == 'rrtmg' )  THEN
             ts_value(dots_start_index_rtm+8,sr)  = hom(nzb,1,108,sr)  ! rrtm_aldif
             ts_value(dots_start_index_rtm+9,sr)  = hom(nzb,1,109,sr)  ! rrtm_aldir
             ts_value(dots_start_index_rtm+10,sr) = hom(nzb,1,110,sr)  ! rrtm_asdif
             ts_value(dots_start_index_rtm+11,sr) = hom(nzb,1,111,sr)  ! rrtm_asdir
          ENDIF

       ENDIF

    ENDIF

 END SUBROUTINE radiation_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sampling of radiation variables along customized measurement coordinates.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_vm_sampling( variable, var_atmos, i_atmos, j_atmos, k_atmos, ns_atmos,       &
                                   var_soil, i_soil, j_soil, k_soil, ns_soil, sampled )

    CHARACTER(LEN=*) ::  variable  !< treated variable

    INTEGER(iwp) ::  i         !< grid index in x-direction
    INTEGER(iwp) ::  j         !< grid index in y-direction
    INTEGER(iwp) ::  m         !< running index over all virtual observation coordinates
    INTEGER(iwp) ::  ns_atmos  !< number of sampling points for atmosphere and surface variables
    INTEGER(iwp) ::  ns_soil   !< number of sampling points for soil variables

    INTEGER(iwp), DIMENSION(1:ns_atmos) ::  i_atmos  !< sampling index in x-direction for atmosphere variables
    INTEGER(iwp), DIMENSION(1:ns_atmos) ::  j_atmos  !< sampling index in y-direction for atmosphere variables
    INTEGER(iwp), DIMENSION(1:ns_atmos) ::  k_atmos  !< sampling index in z-direction for atmosphere variables

    INTEGER(iwp), DIMENSION(1:ns_soil) ::   i_soil   !< sampling index in x-direction for soil variables
    INTEGER(iwp), DIMENSION(1:ns_soil) ::   j_soil   !< sampling index in y-direction for soil variables
    INTEGER(iwp), DIMENSION(1:ns_soil) ::   k_soil   !< sampling index in z-direction for soil variables

    LOGICAL ::  sampled !< flag indicating whether a variable has been sampled

    REAL(wp), DIMENSION(1:ns_atmos) ::  var_atmos  !< array to store atmosphere variables

    REAL(wp), DIMENSION(1:ns_soil) ::  var_soil   !< array to store soil variables


    SELECT CASE ( TRIM( variable ) )
!
!--    Shortwave incoming radiation at surface - diffuse part.
       CASE ( 'rsddif' )
          DO  m = 1, ns_atmos
             j = j_atmos(m)
             i = i_atmos(m)
             var_atmos(m) = rad_sw_in_diff(j,i)
          ENDDO
          sampled = .TRUE.

       CASE DEFAULT

    END SELECT
!
!-- Avoid compiler warning for unused variables by constructing an if condition which is never
!-- fulfilled.
    IF ( .FALSE.  .AND.  ns_atmos < 0  .AND.  ns_soil < 0 )  THEN
       i_soil = i_soil
       j_soil = j_soil
       k_soil = k_soil
       k_atmos = k_atmos
       var_soil = var_soil
    ENDIF

 END SUBROUTINE radiation_vm_sampling


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes global restart data
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_wrd_global

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL wrd_write_string( 'albedo_eff' )
       WRITE( 14 )  albedo_eff

       CALL wrd_write_string( 'dt_radiation' )
       WRITE( 14 )  dt_radiation

       CALL wrd_write_string( 'emissivity_eff' )
       WRITE( 14 )  emissivity_eff

       CALL wrd_write_string( 't_rad_eff' )
       WRITE( 14 )  t_rad_eff

       CALL wrd_write_string( 'time_radiation' )
       WRITE( 14 )  time_radiation

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       CALL wrd_mpi_io( 'albedo_eff',     albedo_eff     )
       CALL wrd_mpi_io( 'dt_radiation',   dt_radiation   )
       CALL wrd_mpi_io( 'emissivity_eff', emissivity_eff )
       CALL wrd_mpi_io( 't_rad_eff',      t_rad_eff      )
       CALL wrd_mpi_io( 'time_radiation', time_radiation )

    ENDIF

 END SUBROUTINE radiation_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes local (subdomain) restart data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_wrd_local

    IMPLICIT NONE

    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_end_index    !< end index for surface data (MPI-IO)
    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_start_index  !< start index for surface data (MPI-IO)

    LOGICAL ::  surface_data_to_write  !< switch for MPI-I/O if PE has surface data to write

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  tmp  !< temporary array for reading from file


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       IF ( ALLOCATED( rad_net_av ) )  THEN
          CALL wrd_write_string( 'rad_net_av' )
          WRITE( 14 )  rad_net_av
       ENDIF

       IF ( ALLOCATED( rad_lw_in_xy_av ) )  THEN
          CALL wrd_write_string( 'rad_lw_in_xy_av' )
          WRITE( 14 )  rad_lw_in_xy_av
       ENDIF

       IF ( ALLOCATED( rad_lw_out_xy_av ) )  THEN
          CALL wrd_write_string( 'rad_lw_out_xy_av' )
          WRITE( 14 )  rad_lw_out_xy_av
       ENDIF

       IF ( ALLOCATED( rad_sw_in_xy_av ) )  THEN
          CALL wrd_write_string( 'rad_sw_in_xy_av' )
          WRITE( 14 )  rad_sw_in_xy_av
       ENDIF

       IF ( ALLOCATED( rad_sw_out_xy_av ) )  THEN
          CALL wrd_write_string( 'rad_sw_out_xy_av' )
          WRITE( 14 )  rad_sw_out_xy_av
       ENDIF

       IF ( ALLOCATED( rad_lw_in ) )  THEN
          CALL wrd_write_string( 'rad_lw_in' )
          WRITE( 14 )  rad_lw_in
       ENDIF

       IF ( ALLOCATED( rad_lw_in_av ) )  THEN
          CALL wrd_write_string( 'rad_lw_in_av' )
          WRITE( 14 )  rad_lw_in_av
       ENDIF

       IF ( ALLOCATED( rad_lw_in_diff) )  THEN
          CALL wrd_write_string( 'rad_lw_in_diff' )
          WRITE( 14 )  rad_lw_in_diff
       ENDIF

       IF ( ALLOCATED( rad_lw_out ) )  THEN
          CALL wrd_write_string( 'rad_lw_out' )
          WRITE( 14 )  rad_lw_out
       ENDIF

       IF ( ALLOCATED( rad_lw_out_av) )  THEN
          CALL wrd_write_string( 'rad_lw_out_av' )
          WRITE( 14 )  rad_lw_out_av
       ENDIF

       IF ( ALLOCATED( rad_lw_cs_hr) )  THEN
          CALL wrd_write_string( 'rad_lw_cs_hr' )
          WRITE( 14 )  rad_lw_cs_hr
       ENDIF

       IF ( ALLOCATED( rad_lw_cs_hr_av) )  THEN
          CALL wrd_write_string( 'rad_lw_cs_hr_av' )
          WRITE( 14 )  rad_lw_cs_hr_av
       ENDIF

       IF ( ALLOCATED( rad_lw_hr) )  THEN
          CALL wrd_write_string( 'rad_lw_hr' )
          WRITE( 14 )  rad_lw_hr
       ENDIF

       IF ( ALLOCATED( rad_lw_hr_av) )  THEN
          CALL wrd_write_string( 'rad_lw_hr_av' )
          WRITE( 14 )  rad_lw_hr_av
       ENDIF

       IF ( ALLOCATED( rad_sw_in) )  THEN
          CALL wrd_write_string( 'rad_sw_in' )
          WRITE( 14 )  rad_sw_in
       ENDIF

       IF ( ALLOCATED( rad_sw_in_diff) )  THEN
          CALL wrd_write_string( 'rad_sw_in_diff' )
          WRITE( 14 )  rad_sw_in_diff
       ENDIF

       IF ( ALLOCATED( rad_sw_in_dir) )  THEN
          CALL wrd_write_string( 'rad_sw_in_dir' )
          WRITE( 14 )  rad_sw_in_dir
       ENDIF

       IF ( ALLOCATED( rad_sw_in_av) )  THEN
          CALL wrd_write_string( 'rad_sw_in_av' )
          WRITE( 14 )  rad_sw_in_av
       ENDIF

       IF ( ALLOCATED( rad_sw_out) )  THEN
          CALL wrd_write_string( 'rad_sw_out' )
          WRITE( 14 )  rad_sw_out
       ENDIF

       IF ( ALLOCATED( rad_sw_out_av) )  THEN
          CALL wrd_write_string( 'rad_sw_out_av' )
          WRITE( 14 )  rad_sw_out_av
       ENDIF

       IF ( ALLOCATED( rad_sw_cs_hr) )  THEN
          CALL wrd_write_string( 'rad_sw_cs_hr' )
          WRITE( 14 )  rad_sw_cs_hr
       ENDIF

       IF ( ALLOCATED( rad_sw_cs_hr_av) )  THEN
          CALL wrd_write_string( 'rad_sw_cs_hr_av' )
          WRITE( 14 )  rad_sw_cs_hr_av
       ENDIF

       IF ( ALLOCATED( rad_sw_hr) )  THEN
          CALL wrd_write_string( 'rad_sw_hr' )
          WRITE( 14 )  rad_sw_hr
       ENDIF

       IF ( ALLOCATED( rad_sw_hr_av) )  THEN
          CALL wrd_write_string( 'rad_sw_hr_av' )
          WRITE( 14 )  rad_sw_hr_av
       ENDIF
!
!--    Write surface related data. Write for LSM and USM separately.
!--    Land-surface radiation.
       CALL wrd_write_string( 'ns_on_file_lsm_rad' )
       WRITE ( 14 )  surf_lsm%ns

       CALL wrd_write_string( 'lsm_start_index_rad' )
       WRITE ( 14 )  surf_lsm%start_index

       CALL wrd_write_string( 'lsm_end_index_rad' )
       WRITE ( 14 )  surf_lsm%end_index

       IF ( ALLOCATED( surf_lsm%rad_lw_in) )  THEN
          CALL wrd_write_string( 'surf_lsm%rad_lw_in' )
          WRITE( 14 )  surf_lsm%rad_lw_in
       ENDIF

       IF ( ALLOCATED( surf_lsm%rad_lw_out) )  THEN
          CALL wrd_write_string( 'surf_lsm%rad_lw_out' )
          WRITE( 14 )  surf_lsm%rad_lw_out
       ENDIF

       IF ( ALLOCATED( surf_lsm%rad_sw_in) )  THEN
          CALL wrd_write_string( 'surf_lsm%rad_sw_in' )
          WRITE( 14 )  surf_lsm%rad_sw_in
       ENDIF

       IF ( ALLOCATED( surf_lsm%rad_sw_out) )  THEN
          CALL wrd_write_string( 'surf_lsm%rad_sw_out' )
          WRITE( 14 )  surf_lsm%rad_sw_out
       ENDIF

       IF ( ALLOCATED( surf_lsm%rad_net) )  THEN
          CALL wrd_write_string( 'surf_lsm%rad_net' )
          WRITE( 14 )  surf_lsm%rad_net
       ENDIF

       IF ( ALLOCATED( surf_lsm%rad_lw_out_change_0) )  THEN
          CALL wrd_write_string( 'surf_lsm%rad_lw_out_change_0' )
          WRITE( 14 )  surf_lsm%rad_lw_out_change_0
       ENDIF
!
!--    Urban-surface radiation.
       CALL wrd_write_string( 'ns_on_file_usm_rad' )
       WRITE( 14 )  surf_usm%ns

       CALL wrd_write_string( 'usm_start_index_rad' )
       WRITE( 14 )  surf_usm%start_index

       CALL wrd_write_string( 'usm_end_index_rad' )
       WRITE( 14 )  surf_usm%end_index

       IF ( ALLOCATED( surf_usm%rad_lw_in) )  THEN
          CALL wrd_write_string( 'surf_usm%rad_lw_in' )
          WRITE( 14 )  surf_usm%rad_lw_in
       ENDIF

       IF ( ALLOCATED( surf_usm%rad_lw_out) )  THEN
          CALL wrd_write_string( 'surf_usm%rad_lw_out' )
          WRITE( 14 )  surf_usm%rad_lw_out
       ENDIF

       IF ( ALLOCATED( surf_usm%rad_sw_in) )  THEN
          CALL wrd_write_string( 'surf_usm%rad_sw_in' )
          WRITE( 14 )  surf_usm%rad_sw_in
       ENDIF

       IF ( ALLOCATED( surf_usm%rad_sw_out) )  THEN
          CALL wrd_write_string( 'surf_usm%rad_sw_out' )
          WRITE( 14 )  surf_usm%rad_sw_out
       ENDIF

       IF ( ALLOCATED( surf_usm%rad_net) )  THEN
          CALL wrd_write_string( 'surf_usm%rad_net' )
          WRITE( 14 )  surf_usm%rad_net
       ENDIF

       IF ( ALLOCATED( surf_usm%rad_lw_out_change_0) )  THEN
          CALL wrd_write_string( 'surf_usm%rad_lw_out_change_0' )
          WRITE( 14 )  surf_usm%rad_lw_out_change_0
       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       IF ( ALLOCATED( rad_net_av )       )  CALL wrd_mpi_io( 'rad_net_av', rad_net_av )
       IF ( ALLOCATED( rad_lw_in_xy_av )  )  CALL wrd_mpi_io( 'rad_lw_in_xy_av', rad_lw_in_xy_av )
       IF ( ALLOCATED( rad_lw_out_xy_av ) )  CALL wrd_mpi_io( 'rad_lw_out_xy_av', rad_lw_out_xy_av )
       IF ( ALLOCATED( rad_sw_in_xy_av )  )  CALL wrd_mpi_io( 'rad_sw_in_xy_av', rad_sw_in_xy_av )
       IF ( ALLOCATED( rad_sw_out_xy_av ) )  CALL wrd_mpi_io( 'rad_sw_out_xy_av', rad_sw_out_xy_av )

       IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.            &
            radiation_scheme == 'external' )                                                       &
       THEN
          IF ( ALLOCATED( rad_lw_in ) )  THEN
             tmp = rad_lw_in(0,:,:)
             CALL wrd_mpi_io( 'rad_lw_in', tmp )
          ENDIF
          IF ( ALLOCATED( rad_lw_in_av ) )  THEN
             tmp = rad_lw_in_av(0,:,:)
             CALL wrd_mpi_io( 'rad_lw_in_av', tmp )
          ENDIF
          IF ( ALLOCATED( rad_lw_out ) )  THEN
             tmp = rad_lw_out(0,:,:)
             CALL wrd_mpi_io( 'rad_lw_out', tmp )
          ENDIF
          IF ( ALLOCATED( rad_lw_out_av ) )  THEN
             tmp = rad_lw_out_av(0,:,:)
             CALL wrd_mpi_io( 'rad_lw_out_av', tmp )
          ENDIF
       ELSE
          IF ( ALLOCATED( rad_lw_in )    )  CALL wrd_mpi_io( 'rad_lw_in', rad_lw_in )
          IF ( ALLOCATED( rad_lw_in_av ) )  CALL wrd_mpi_io( 'rad_lw_in_av', rad_lw_in_av )
          IF ( ALLOCATED( rad_lw_out )   )  CALL wrd_mpi_io( 'rad_lw_out', rad_lw_out )
          IF ( ALLOCATED( rad_lw_out_av) )  CALL wrd_mpi_io( 'rad_lw_out_av', rad_lw_out_av )
       ENDIF

       IF ( ALLOCATED( rad_lw_in_diff) )  CALL wrd_mpi_io( 'rad_lw_in_diff', rad_lw_in_diff )
       IF ( ALLOCATED( rad_sw_in_diff) )  CALL wrd_mpi_io( 'rad_sw_in_diff', rad_sw_in_diff )
       IF ( ALLOCATED( rad_sw_in_dir)  )  CALL wrd_mpi_io( 'rad_sw_in_dir', rad_sw_in_dir   )

       IF ( ALLOCATED( rad_lw_cs_hr)    )  CALL wrd_mpi_io( 'rad_lw_cs_hr', rad_lw_cs_hr )
       IF ( ALLOCATED( rad_lw_cs_hr_av) )  CALL wrd_mpi_io( 'rad_lw_cs_hr_av', rad_lw_cs_hr_av )
       IF ( ALLOCATED( rad_lw_hr)       )  CALL wrd_mpi_io( 'rad_lw_hr', rad_lw_hr )
       IF ( ALLOCATED( rad_lw_hr_av)    )  CALL wrd_mpi_io( 'rad_lw_hr_av', rad_lw_hr_av )

       IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.            &
            radiation_scheme == 'external' )                                                       &
       THEN
          IF ( ALLOCATED( rad_sw_in ) )  THEN
             tmp = rad_sw_in(0,:,:)
             CALL wrd_mpi_io( 'rad_sw_in', tmp )
          ENDIF
          IF ( ALLOCATED( rad_sw_in_av ) )  THEN
             tmp = rad_sw_in_av(0,:,:)
             CALL wrd_mpi_io( 'rad_sw_in_av', tmp )
          ENDIF
          IF ( ALLOCATED( rad_sw_out ) )  THEN
             tmp = rad_sw_out(0,:,:)
             CALL wrd_mpi_io( 'rad_sw_out', tmp )
          ENDIF
          IF ( ALLOCATED( rad_sw_out_av ) )  THEN
             tmp = rad_sw_out_av(0,:,:)
             CALL wrd_mpi_io( 'rad_sw_out_av', tmp )
          ENDIF
       ELSE
          IF ( ALLOCATED( rad_sw_in)     )  CALL wrd_mpi_io( 'rad_sw_in', rad_sw_in )
          IF ( ALLOCATED( rad_sw_in_av)  )  CALL wrd_mpi_io( 'rad_sw_in_av', rad_sw_in_av )
          IF ( ALLOCATED( rad_sw_out)    )  CALL wrd_mpi_io( 'rad_sw_out', rad_sw_out )
          IF ( ALLOCATED( rad_sw_out_av) )  CALL wrd_mpi_io( 'rad_sw_out_av', rad_sw_out_av )
       ENDIF
       IF ( ALLOCATED( rad_sw_cs_hr)    )  CALL wrd_mpi_io( 'rad_sw_cs_hr', rad_sw_cs_hr )
       IF ( ALLOCATED( rad_sw_cs_hr_av) )  CALL wrd_mpi_io( 'rad_sw_cs_hr_av', rad_sw_cs_hr_av )
       IF ( ALLOCATED( rad_sw_hr)       )  CALL wrd_mpi_io( 'rad_sw_hr', rad_sw_hr )
       IF ( ALLOCATED( rad_sw_hr_av)    )  CALL wrd_mpi_io( 'rad_sw_hr_av', rad_sw_hr_av )

!
!--    Write local surface data. Distinguish between LSM and USM surfaces. This way, usage of
!--    lengthy surface-restore routines can be avoided. Start with LSM surface data.
       CALL rd_mpi_io_surface_filetypes( surf_lsm%start_index, surf_lsm%end_index,                 &
                                         surface_data_to_write, global_start_index,                &
                                         global_end_index )

       CALL wrd_mpi_io( 'rad_lsm_global_start', global_start_index )
       CALL wrd_mpi_io( 'rad_lsm_global_end',   global_end_index   )

       IF ( surface_data_to_write )  THEN
          CALL wrd_mpi_io_surface( 'surf_lsm%rad_lw_in',           surf_lsm%rad_lw_in  )
          CALL wrd_mpi_io_surface( 'surf_lsm%rad_lw_out',          surf_lsm%rad_lw_out )
          CALL wrd_mpi_io_surface( 'surf_lsm%rad_sw_in',           surf_lsm%rad_sw_in  )
          CALL wrd_mpi_io_surface( 'surf_lsm%rad_sw_out',          surf_lsm%rad_sw_out )
          CALL wrd_mpi_io_surface( 'surf_lsm%rad_net',             surf_lsm%rad_net    )
          CALL wrd_mpi_io_surface( 'surf_lsm%rad_lw_out_change_0', surf_lsm%rad_lw_out_change_0 )
       ENDIF
!
!--    USM surface data.
       CALL rd_mpi_io_surface_filetypes( surf_usm%start_index, surf_usm%end_index,                 &
                                         surface_data_to_write, global_start_index,                &
                                         global_end_index )

       CALL wrd_mpi_io( 'rad_usm_global_start', global_start_index )
       CALL wrd_mpi_io( 'rad_usm_global_end',   global_end_index   )

       IF ( surface_data_to_write )  THEN
          CALL wrd_mpi_io_surface( 'surf_usm%rad_lw_in',           surf_usm%rad_lw_in  )
          CALL wrd_mpi_io_surface( 'surf_usm%rad_lw_out',          surf_usm%rad_lw_out )
          CALL wrd_mpi_io_surface( 'surf_usm%rad_sw_in',           surf_usm%rad_sw_in  )
          CALL wrd_mpi_io_surface( 'surf_usm%rad_sw_out',          surf_usm%rad_sw_out )
          CALL wrd_mpi_io_surface( 'surf_usm%rad_net',             surf_usm%rad_net    )
          CALL wrd_mpi_io_surface( 'surf_usm%rad_lw_out_change_0', surf_usm%rad_lw_out_change_0 )
       ENDIF
!
    ENDIF

 END SUBROUTINE radiation_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf,    &
                                     nync, nyn_on_file, nysf, nysc, nys_on_file, tmp_2d, tmp_3d,   &
                                     found )


    USE control_parameters

    USE kinds

    USE pegrid


    IMPLICIT NONE

    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  ns_on_file_lsm_rad  !< number of land-surface elements on file
    INTEGER(iwp) ::  ns_on_file_usm_rad  !< number of urban-surface elements on file
    INTEGER(iwp) ::  nxlc                !<
    INTEGER(iwp) ::  nxlf                !<
    INTEGER(iwp) ::  nxl_on_file         !<
    INTEGER(iwp) ::  nxrc                !<
    INTEGER(iwp) ::  nxrf                !<
    INTEGER(iwp) ::  nxr_on_file         !<
    INTEGER(iwp) ::  nync                !<
    INTEGER(iwp) ::  nynf                !<
    INTEGER(iwp) ::  nyn_on_file         !<
    INTEGER(iwp) ::  nysc                !<
    INTEGER(iwp) ::  nysf                !<
    INTEGER(iwp) ::  nys_on_file         !<

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  end_index_lsm_on_file    !< end index of LSM surface elements at (j,i)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  end_index_usm_on_file    !< end index of USM surface elements at (j,i)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  start_index_lsm_on_file  !< start index of LSM surface elements at (j,i)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  start_index_usm_on_file  !< start index of LSM surface elements at (j,i)

    LOGICAL, INTENT(OUT) ::  found  !<

    REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::  tmp_2d  !<

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::  tmp_3d  !<

    REAL(wp), DIMENSION(0:0,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) ::  tmp_3d2  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  tmp_surf_lsm  !< temporary variable to read surface data (LSM)
    REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::  tmp_surf_usm  !< temporary variable to read surface data (LSM)

    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'rad_net_av' )
          IF ( .NOT. ALLOCATED( rad_net_av ) )  THEN
             ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_2d
          rad_net_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_in_xy_av' )
          IF ( .NOT. ALLOCATED( rad_lw_in_xy_av ) )  THEN
             ALLOCATE( rad_lw_in_xy_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_2d
          rad_lw_in_xy_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_out_xy_av' )
          IF ( .NOT. ALLOCATED( rad_lw_out_xy_av ) )  THEN
             ALLOCATE( rad_lw_out_xy_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_2d
          rad_lw_out_xy_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_in_xy_av' )
          IF ( .NOT. ALLOCATED( rad_sw_in_xy_av ) )  THEN
             ALLOCATE( rad_sw_in_xy_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_2d
          rad_sw_in_xy_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                               &
          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_out_xy_av' )
          IF ( .NOT. ALLOCATED( rad_sw_out_xy_av ) )  THEN
             ALLOCATE( rad_sw_out_xy_av(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_2d
          rad_sw_out_xy_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
          tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_in' )
          IF ( .NOT. ALLOCATED( rad_lw_in ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_lw_in(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                READ( 13 )  tmp_3d2
                rad_lw_in(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
                tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ( 13 )  tmp_3d
                rad_lw_in(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                             &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_in_av' )
          IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_lw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                READ( 13 )  tmp_3d2
                rad_lw_in_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                        &
                tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ( 13 )  tmp_3d
                rad_lw_in_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_in_diff' )
          IF ( .NOT. ALLOCATED( rad_lw_in_diff ) )  THEN
             ALLOCATE( rad_lw_in_diff(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_2d
          rad_lw_in_diff(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
                                                    tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_out' )
          IF ( .NOT. ALLOCATED( rad_lw_out ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_lw_out(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant' .OR.       &
                  radiation_scheme == 'external' )  THEN
                READ( 13 )  tmp_3d2
                rad_lw_out(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
                tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ( 13 )  tmp_3d
                rad_lw_out(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_out_av' )
          IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_lw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                READ( 13 )  tmp_3d2
                rad_lw_out_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                       &
                tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ( 13 )  tmp_3d
                rad_lw_out_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                         &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_lw_cs_hr' )
          IF ( .NOT. ALLOCATED( rad_lw_cs_hr ) )  THEN
             ALLOCATE( rad_lw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_3d
          rad_lw_cs_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
          tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_cs_hr_av' )
          IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )  THEN
             ALLOCATE( rad_lw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_3d
          rad_lw_cs_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                             &
          tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_hr' )
          IF ( .NOT. ALLOCATED( rad_lw_hr ) )  THEN
             ALLOCATE( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_3d
          rad_lw_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                   &
          tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_lw_hr_av' )
          IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )  THEN
             ALLOCATE( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_3d
          rad_lw_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
          tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_in' )
          IF ( .NOT. ALLOCATED( rad_sw_in ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_sw_in(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                READ( 13 )  tmp_3d2
                rad_sw_in(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                           &
                tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ( 13 )  tmp_3d
                rad_sw_in(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                             &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_in_av' )
          IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_sw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                READ( 13 )  tmp_3d2
                rad_sw_in_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                        &
                tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ( 13 )  tmp_3d
                rad_sw_in_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_in_diff' )
          IF ( .NOT. ALLOCATED( rad_sw_in_diff ) )  THEN
             ALLOCATE( rad_sw_in_diff(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_2d
          rad_sw_in_diff(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
                                                    tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_in_dir' )
          IF ( .NOT. ALLOCATED( rad_sw_in_dir ) )  THEN
             ALLOCATE( rad_sw_in_dir(nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_2d
          rad_sw_in_dir(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                 &
                                                    tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_out' )
          IF ( .NOT. ALLOCATED( rad_sw_out ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_sw_out(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                READ( 13 )  tmp_3d2
                rad_sw_out(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                          &
                tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ( 13 )  tmp_3d
                rad_sw_out(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                            &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_out_av' )
          IF ( .NOT. ALLOCATED( rad_sw_out_av ) )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                ALLOCATE( rad_sw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
             ELSE
                ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
          ENDIF
          IF ( k == 1 )  THEN
             IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.      &
                  radiation_scheme == 'external' )  THEN
                READ( 13 )  tmp_3d2
                rad_sw_out_av(0:0,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                       &
                tmp_3d2(0:0,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ELSE
                READ( 13 )  tmp_3d
                rad_sw_out_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                         &
                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             ENDIF
          ENDIF

       CASE ( 'rad_sw_cs_hr' )
          IF ( .NOT. ALLOCATED( rad_sw_cs_hr ) )  THEN
             ALLOCATE( rad_sw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_3d
          rad_sw_cs_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
          tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_cs_hr_av' )
          IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )  THEN
             ALLOCATE( rad_sw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_3d
          rad_sw_cs_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                             &
          tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_hr' )
          IF ( .NOT. ALLOCATED( rad_sw_hr ) )  THEN
             ALLOCATE( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_3d
          rad_sw_hr(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                   &
          tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'rad_sw_hr_av' )
          IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )  THEN
             ALLOCATE( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ( 13 )  tmp_3d
          rad_lw_hr_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
          tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'ns_on_file_lsm_rad')
          IF ( k == 1 )  THEN
             READ ( 13 ) ns_on_file_lsm_rad
!
!--          In case of changing mpi topology, this routine could be called more than once.
!--          Hence, arrays need to be deallocated before allocated again.
             IF ( ALLOCATED( tmp_surf_lsm ) )  DEALLOCATE( tmp_surf_lsm )
!
!--          Allocate temporary arrays for reading data on file. Note, the size of allocated surface
!--          elements do not necessarily need to match the size of present surface elements on
!--          current processor, as the number of processors between restarts can change.
             ALLOCATE( tmp_surf_lsm(1:ns_on_file_lsm_rad) )
          ENDIF

       CASE ( 'lsm_start_index_rad' )
          IF ( k == 1 )  THEN
             IF ( ALLOCATED( start_index_lsm_on_file ) )  DEALLOCATE( start_index_lsm_on_file )
             ALLOCATE ( start_index_lsm_on_file(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) )
             READ ( 13 )  start_index_lsm_on_file
          ENDIF

       CASE ( 'lsm_end_index_rad' )
          IF ( k == 1 )  THEN
             IF ( ALLOCATED( end_index_lsm_on_file ) )  DEALLOCATE( end_index_lsm_on_file )
             ALLOCATE ( end_index_lsm_on_file(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) )
             READ ( 13 )  end_index_lsm_on_file
          ENDIF

       CASE ( 'surf_lsm%rad_lw_in' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_lsm%rad_lw_in ) )                                         &
                ALLOCATE( surf_lsm%rad_lw_in(1:surf_lsm%ns) )
             READ ( 13 )  tmp_surf_lsm
          ENDIF
          CALL surface_restore_elements( surf_lsm%rad_lw_in, tmp_surf_lsm,                         &
                                         surf_lsm%start_index, start_index_lsm_on_file,            &
                                         end_index_lsm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_lsm%rad_lw_out' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_lsm%rad_lw_out ) )                                        &
                ALLOCATE( surf_lsm%rad_lw_out(1:surf_lsm%ns) )
             READ ( 13 )  tmp_surf_lsm
          ENDIF
          CALL surface_restore_elements( surf_lsm%rad_lw_out, tmp_surf_lsm,                        &
                                         surf_lsm%start_index, start_index_lsm_on_file,            &
                                         end_index_lsm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_lsm%rad_sw_in' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_lsm%rad_sw_in ) )                                         &
                ALLOCATE( surf_lsm%rad_sw_in(1:surf_lsm%ns) )
             READ ( 13 )  tmp_surf_lsm
          ENDIF
          CALL surface_restore_elements( surf_lsm%rad_sw_in, tmp_surf_lsm,                         &
                                         surf_lsm%start_index, start_index_lsm_on_file,            &
                                         end_index_lsm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_lsm%rad_sw_out' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_lsm%rad_sw_out ) )                                        &
                ALLOCATE( surf_lsm%rad_sw_out(1:surf_lsm%ns) )
             READ ( 13 )  tmp_surf_lsm
          ENDIF
          CALL surface_restore_elements( surf_lsm%rad_sw_out, tmp_surf_lsm,                        &
                                         surf_lsm%start_index, start_index_lsm_on_file,            &
                                         end_index_lsm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_lsm%rad_net' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_lsm%rad_net ) )                                           &
                ALLOCATE( surf_lsm%rad_net(1:surf_lsm%ns) )
             READ ( 13 )  tmp_surf_lsm
          ENDIF
          CALL surface_restore_elements( surf_lsm%rad_net, tmp_surf_lsm,                           &
                                         surf_lsm%start_index, start_index_lsm_on_file,            &
                                         end_index_lsm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_lsm%rad_lw_out_change_0' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_lsm%rad_lw_out_change_0 ) )                               &
                ALLOCATE( surf_lsm%rad_lw_out_change_0(1:surf_lsm%ns) )
             READ ( 13 )  tmp_surf_lsm
          ENDIF
          CALL surface_restore_elements( surf_lsm%rad_lw_out_change_0, tmp_surf_lsm,               &
                                         surf_lsm%start_index, start_index_lsm_on_file,            &
                                         end_index_lsm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'ns_on_file_usm_rad')
          IF ( k == 1 )  THEN
             READ ( 13 ) ns_on_file_usm_rad
!
!--          In case of changing mpi topology, this routine could be called more than once.
!--          Hence, arrays need to be deallocated before allocated again.
             IF ( ALLOCATED( tmp_surf_usm ) )  DEALLOCATE( tmp_surf_usm )
!
!--          Allocate temporary arrays for reading data on file. Note, the size of allocated surface
!--          elements do not necessarily need to match the size of present surface elements on
!--          current processor, as the number of processors between restarts can change.
             ALLOCATE( tmp_surf_usm(1:ns_on_file_usm_rad) )
          ENDIF

       CASE ( 'usm_start_index_rad' )
          IF ( k == 1 )  THEN
             IF ( ALLOCATED( start_index_usm_on_file ) )  DEALLOCATE( start_index_usm_on_file )
             ALLOCATE ( start_index_usm_on_file(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) )
             READ ( 13 )  start_index_usm_on_file
          ENDIF

       CASE ( 'usm_end_index_rad' )
          IF ( k == 1 )  THEN
             IF ( ALLOCATED( end_index_usm_on_file ) )  DEALLOCATE( end_index_usm_on_file )
             ALLOCATE ( end_index_usm_on_file(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) )
             READ ( 13 )  end_index_usm_on_file
          ENDIF
       CASE ( 'surf_usm%rad_lw_in' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm%rad_lw_in ) )                                         &
                ALLOCATE( surf_usm%rad_lw_in(1:surf_usm%ns) )
             READ ( 13 )  tmp_surf_usm
          ENDIF
          CALL surface_restore_elements( surf_usm%rad_lw_in, tmp_surf_usm,                         &
                                         surf_usm%start_index, start_index_usm_on_file,            &
                                         end_index_usm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_usm%rad_lw_out' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm%rad_lw_out ) )                                        &
                ALLOCATE( surf_usm%rad_lw_out(1:surf_usm%ns) )
             READ ( 13 )  tmp_surf_usm
          ENDIF
          CALL surface_restore_elements( surf_usm%rad_lw_out, tmp_surf_usm,                        &
                                         surf_usm%start_index, start_index_usm_on_file,            &
                                         end_index_usm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_usm%rad_sw_in' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm%rad_sw_in ) )                                         &
                ALLOCATE( surf_usm%rad_sw_in(1:surf_usm%ns) )
             READ ( 13 )  tmp_surf_usm
          ENDIF
          CALL surface_restore_elements( surf_usm%rad_sw_in, tmp_surf_usm,                         &
                                         surf_usm%start_index, start_index_usm_on_file,            &
                                         end_index_usm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_usm%rad_sw_out' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm%rad_sw_out ) )                                        &
                ALLOCATE( surf_usm%rad_sw_out(1:surf_usm%ns) )
             READ ( 13 )  tmp_surf_usm
          ENDIF
          CALL surface_restore_elements( surf_usm%rad_sw_out, tmp_surf_usm,                        &
                                         surf_usm%start_index, start_index_usm_on_file,            &
                                         end_index_usm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_usm%rad_net' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm%rad_net ) )                                           &
                ALLOCATE( surf_usm%rad_net(1:surf_usm%ns) )
             READ ( 13 )  tmp_surf_usm
          ENDIF
          CALL surface_restore_elements( surf_usm%rad_net, tmp_surf_usm,                           &
                                         surf_usm%start_index, start_index_usm_on_file,            &
                                         end_index_usm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 'surf_usm%rad_lw_out_change_0' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm%rad_lw_out_change_0 ) )                               &
                ALLOCATE( surf_usm%rad_lw_out_change_0(1:surf_usm%ns) )
             READ ( 13 )  tmp_surf_usm
          ENDIF
          CALL surface_restore_elements( surf_usm%rad_lw_out_change_0, tmp_surf_usm,               &
                                         surf_usm%start_index, start_index_usm_on_file,            &
                                         end_index_usm_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,&
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE radiation_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE radiation_rrd_local_mpi

    USE control_parameters

    USE indices

    USE kinds


    IMPLICIT NONE

    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_end_index
    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_start_index

    LOGICAL ::  array_found  !< flag indicating that enquired array was found
    LOGICAL ::  data_to_read !< dummy variable

    REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  tmp  !< temporary array for reading from file

!
!-- Read arrays needed for time-averaging. Note, this is skipped in case of cyclic_fill runs.
    IF ( .NOT. cyclic_fill_initialization )  THEN
       CALL rd_mpi_io_check_array( 'rad_net_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rad_net_av ) )  ALLOCATE( rad_net_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_net_av', rad_net_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_lw_in_xy_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rad_lw_in_xy_av ) )  ALLOCATE( rad_lw_in_xy_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_lw_in_xy_av', rad_lw_in_xy_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_lw_out_xy_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rad_lw_out_xy_av ) )  ALLOCATE( rad_lw_out_xy_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_lw_out_xy_av', rad_lw_out_xy_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_sw_in_xy_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rad_sw_in_xy_av ) )  ALLOCATE( rad_sw_in_xy_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_sw_in_xy_av', rad_sw_in_xy_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_sw_out_xy_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rad_sw_out_xy_av ) )  ALLOCATE( rad_sw_out_xy_av(nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_sw_out_xy_av', rad_sw_out_xy_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_lw_in_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.         &
               radiation_scheme == 'external' )  THEN
             IF ( .NOT. ALLOCATED( rad_lw_in_av ) )  ALLOCATE( rad_lw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'rad_lw_in_av', tmp )
             rad_lw_in_av(0,:,:) = tmp
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_in_av ) )                                                &
             ALLOCATE( rad_lw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'rad_lw_in_av', rad_lw_in_av )
          ENDIF
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_lw_out_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.         &
               radiation_scheme == 'external' )  THEN
             IF ( .NOT. ALLOCATED( rad_lw_out_av ) )  ALLOCATE( rad_lw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'rad_lw_out_av', tmp )
             rad_lw_out_av(0,:,:) = tmp
          ELSE
             IF ( .NOT. ALLOCATED( rad_lw_out_av ) )                                               &
             ALLOCATE( rad_lw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'rad_lw_out_av', rad_lw_out_av )
          ENDIF
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_lw_cs_hr_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rad_lw_cs_hr_av ) )                                                &
          ALLOCATE( rad_lw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_lw_cs_hr_av', rad_lw_cs_hr_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_lw_hr_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rad_lw_hr_av ) )                                                   &
          ALLOCATE( rad_lw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_lw_hr_av', rad_lw_hr_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_sw_in_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.         &
               radiation_scheme == 'external' )  THEN
             IF ( .NOT. ALLOCATED( rad_sw_in_av ) )  ALLOCATE( rad_sw_in_av(0:0,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'rad_sw_in_av', tmp )
             rad_sw_in_av(0,:,:) = tmp
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_in_av ) )                                                &
             ALLOCATE( rad_sw_in_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'rad_sw_in_av', rad_sw_in_av )
          ENDIF
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_sw_out_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.         &
               radiation_scheme == 'external' )  THEN
             IF ( .NOT. ALLOCATED( rad_sw_out_av ) )                                               &
             ALLOCATE( rad_sw_out_av(0:0,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'rad_sw_out_av', tmp )
             rad_sw_out_av(0,:,:) = tmp
          ELSE
             IF ( .NOT. ALLOCATED( rad_sw_out_av ) )                                               &
             ALLOCATE( rad_sw_out_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             CALL rrd_mpi_io( 'rad_sw_out_av', rad_sw_out_av )
          ENDIF
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_sw_cs_hr_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rad_sw_cs_hr_av ) )                                                &
          ALLOCATE( rad_sw_cs_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_sw_cs_hr_av', rad_sw_cs_hr_av )
       ENDIF

       CALL rd_mpi_io_check_array( 'rad_sw_hr_av' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( rad_sw_hr_av ) )                                                   &
          ALLOCATE( rad_sw_hr_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_sw_hr_av', rad_sw_hr_av )
       ENDIF
    ENDIF

!
!-- Arrays required for radiation forcing terms or the RTM in a restart run.
    CALL rd_mpi_io_check_array( 'rad_lw_in' , found = array_found )
    IF ( array_found )  THEN
       IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.            &
            radiation_scheme == 'external' )  THEN
          IF ( .NOT. ALLOCATED( rad_lw_in ) )  ALLOCATE( rad_lw_in(0:0,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_lw_in', tmp )
          rad_lw_in(0,:,:) = tmp
       ELSE
          IF ( .NOT. ALLOCATED( rad_lw_in ) )  ALLOCATE( rad_lw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_lw_in', rad_lw_in )
       ENDIF
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_lw_in_diff' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( rad_lw_in_diff ) )  ALLOCATE( rad_lw_in_diff(nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'rad_lw_in_diff', rad_lw_in_diff )
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_lw_out' , found = array_found )
    IF ( array_found )  THEN
       IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.            &
            radiation_scheme == 'external' )  THEN
          IF ( .NOT. ALLOCATED( rad_lw_out ) )  ALLOCATE( rad_lw_out(0:0,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_lw_out', tmp )
          rad_lw_out(0,:,:) = tmp
       ELSE
          IF ( .NOT. ALLOCATED( rad_lw_out ) )  ALLOCATE( rad_lw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_lw_out', rad_lw_out )
       ENDIF
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_lw_cs_hr' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( rad_lw_cs_hr ) )                                                      &
       ALLOCATE( rad_lw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'rad_lw_cs_hr', rad_lw_cs_hr )
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_lw_hr' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( rad_lw_hr ) )  ALLOCATE( rad_lw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'rad_lw_hr', rad_lw_hr )
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_sw_in' , found = array_found )
    IF ( array_found )  THEN
       IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.            &
            radiation_scheme == 'external' )  THEN
          IF ( .NOT. ALLOCATED( rad_sw_in ) )  ALLOCATE( rad_sw_in(0:0,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_sw_in', tmp )
          rad_sw_in(0,:,:) = tmp
       ELSE
          IF ( .NOT. ALLOCATED( rad_sw_in ) )  ALLOCATE( rad_sw_in(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_sw_in', rad_sw_in )
       ENDIF
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_sw_in_diff' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( rad_sw_in_diff ) )  ALLOCATE( rad_sw_in_diff(nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'rad_sw_in_diff', rad_sw_in_diff )
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_sw_in_dir' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( rad_sw_in_dir ) )  ALLOCATE( rad_sw_in_dir(nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'rad_sw_in_dir', rad_sw_in_dir )
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_sw_out' , found = array_found )
    IF ( array_found )  THEN
       IF ( radiation_scheme == 'clear-sky'  .OR.  radiation_scheme == 'constant'  .OR.            &
            radiation_scheme == 'external' )  THEN
          IF ( .NOT. ALLOCATED( rad_sw_out ) )  ALLOCATE( rad_sw_out(0:0,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_sw_out',  tmp)
          rad_sw_out(0,:,:) = tmp
       ELSE
          IF ( .NOT. ALLOCATED( rad_sw_out ) )                                                     &
          ALLOCATE( rad_sw_out(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'rad_sw_out',  rad_sw_out )
       ENDIF
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_sw_cs_hr' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( rad_sw_cs_hr ) )                                                      &
       ALLOCATE( rad_sw_cs_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'rad_sw_cs_hr', rad_sw_cs_hr )
    ENDIF

    CALL rd_mpi_io_check_array( 'rad_sw_hr' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( rad_sw_hr ) )  ALLOCATE( rad_sw_hr(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'rad_sw_hr', rad_sw_hr )
    ENDIF

!
!-- Now, read surface data. Start with LSM data.
!-- At the moment reading of surface data in combination with cyclic fill is not realized,
!-- so that this is skipped for the moment.
    IF ( cyclic_fill_initialization )  RETURN

    CALL rd_mpi_io_check_array( 'rad_lsm_global_start' , found = array_found )
    IF ( array_found )  CALL rrd_mpi_io( 'rad_lsm_global_start', global_start_index )

    CALL rd_mpi_io_check_array( 'rad_lsm_global_end' , found = array_found )
    IF ( array_found )  CALL rrd_mpi_io( 'rad_lsm_global_end',   global_end_index   )
!
!-- Check if data input for surface-type variables is required. Note, only invoke routine if
!-- surface restart data is on file. In case of cyclic fill initialization this is not necessarily
!-- guaranteed. To check this use the array_found control flag.
    IF ( array_found )  THEN
       CALL rd_mpi_io_surface_filetypes( surf_lsm%start_index, surf_lsm%end_index,                 &
                                         data_to_read, global_start_index, global_end_index )
    ELSE
       data_to_read = .FALSE.
    ENDIF

    IF ( data_to_read )  THEN
       CALL rd_mpi_io_check_array( 'surf_lsm%rad_lw_in' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_lsm%rad_lw_in ) )                                             &
             ALLOCATE( surf_lsm%rad_lw_in(1:surf_lsm%ns) )
          CALL rrd_mpi_io_surface( 'surf_lsm%rad_lw_in', surf_lsm%rad_lw_in )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_lsm%rad_lw_out' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_lsm%rad_lw_out ) )                                            &
             ALLOCATE( surf_lsm%rad_lw_out(1:surf_lsm%ns) )
          CALL rrd_mpi_io_surface( 'surf_lsm%rad_lw_out', surf_lsm%rad_lw_out )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_lsm%rad_sw_in' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_lsm%rad_sw_in ) )                                             &
             ALLOCATE( surf_lsm%rad_sw_in(1:surf_lsm%ns) )
          CALL rrd_mpi_io_surface( 'surf_lsm%rad_sw_in', surf_lsm%rad_sw_in )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_lsm%rad_sw_out' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_lsm%rad_sw_out ) )                                            &
             ALLOCATE( surf_lsm%rad_sw_out(1:surf_lsm%ns) )
          CALL rrd_mpi_io_surface( 'surf_lsm%rad_sw_out', surf_lsm%rad_sw_out )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_lsm%rad_net' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_lsm%rad_net ) )                                               &
             ALLOCATE( surf_lsm%rad_net(1:surf_lsm%ns) )
          CALL rrd_mpi_io_surface( 'surf_lsm%rad_net', surf_lsm%rad_net )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_lsm%rad_lw_out_change_0' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_lsm%rad_lw_out_change_0 ) )                                   &
             ALLOCATE( surf_lsm%rad_lw_out_change_0(1:surf_lsm%ns) )
          CALL rrd_mpi_io_surface( 'surf_lsm%rad_lw_out_change_0', surf_lsm%rad_lw_out_change_0 )
       ENDIF
    ENDIF
!
!-- USM data.
    CALL rd_mpi_io_check_array( 'rad_usm_global_start' , found = array_found )
    IF ( array_found )  CALL rrd_mpi_io( 'rad_usm_global_start', global_start_index )

    CALL rd_mpi_io_check_array( 'rad_usm_global_end' , found = array_found )
    IF ( array_found )  CALL rrd_mpi_io( 'rad_usm_global_end',   global_end_index   )

    IF ( array_found )  THEN
       CALL rd_mpi_io_surface_filetypes( surf_usm%start_index, surf_usm%end_index,                 &
                                         data_to_read, global_start_index, global_end_index )
    ELSE
       data_to_read = .FALSE.
    ENDIF

    IF ( data_to_read )  THEN
       CALL rd_mpi_io_check_array( 'surf_usm%rad_lw_in' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_usm%rad_lw_in ) )                                             &
             ALLOCATE( surf_usm%rad_lw_in(1:surf_usm%ns) )
          CALL rrd_mpi_io_surface( 'surf_usm%rad_lw_in', surf_usm%rad_lw_in )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_usm%rad_lw_out' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_usm%rad_lw_out ) )                                            &
             ALLOCATE( surf_usm%rad_lw_out(1:surf_usm%ns) )
          CALL rrd_mpi_io_surface( 'surf_usm%rad_lw_out', surf_usm%rad_lw_out )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_usm%rad_sw_in' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_usm%rad_sw_in ) )                                             &
             ALLOCATE( surf_usm%rad_sw_in(1:surf_usm%ns) )
          CALL rrd_mpi_io_surface( 'surf_usm%rad_sw_in', surf_usm%rad_sw_in )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_usm%rad_sw_out' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_usm%rad_sw_out ) )                                            &
             ALLOCATE( surf_usm%rad_sw_out(1:surf_usm%ns) )
          CALL rrd_mpi_io_surface( 'surf_usm%rad_sw_out', surf_usm%rad_sw_out )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_usm%rad_net' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_usm%rad_net ) )                                               &
             ALLOCATE( surf_usm%rad_net(1:surf_usm%ns) )
          CALL rrd_mpi_io_surface( 'surf_usm%rad_net', surf_usm%rad_net )
       ENDIF

       CALL rd_mpi_io_check_array( 'surf_usm%rad_lw_out_change_0' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_usm%rad_lw_out_change_0 ) )                                   &
             ALLOCATE( surf_usm%rad_lw_out_change_0(1:surf_usm%ns) )
          CALL rrd_mpi_io_surface( 'surf_usm%rad_lw_out_change_0', surf_usm%rad_lw_out_change_0 )
       ENDIF
    ENDIF

 END SUBROUTINE radiation_rrd_local_mpi


 END MODULE radiation_model_mod
