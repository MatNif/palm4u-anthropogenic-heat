!> @file chemistry_model_mod.f90
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
! Copyright 2017-2021 Leibniz Universitaet Hannover
! Copyright 2017-2021 Karlsruhe Institute of Technology
! Copyright 2017-2021 Freie Universitaet Berlin
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Renate Forkel
! @author Farah Kanani-Suehring
! @author Klaus Ketelsen
! @author Basit Khan
! @author Sabine Banzhaf
! @author Edward C. Chan
! @author Ilona Ilona Jäkel
! @author Matthias Sühring
! @author Pecanode GmbH
!
!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Chemistry model for PALM-4U
!> @todo Adjust chem_rrd_local to CASE structure of others modules. It is not allowed to use the
!>       chemistry model in a precursor run and additionally not using it in a main run
!> @todo Implement turbulent inflow of chem spcs in inflow_turbulence. Do we need this? Not done for salsa either.
!
!--------------------------------------------------------------------------------------------------!

 MODULE chemistry_model_mod

    USE advec_s_pw_mod,                                                                            &
        ONLY:  advec_s_pw

    USE advec_s_up_mod,                                                                            &
        ONLY:  advec_s_up

    USE advec_ws,                                                                                  &
        ONLY:  advec_s_ws,                                                                         &
               ws_init_flags_scalar

    USE arrays_3d,                                                                                 &
        ONLY:  dzw,                                                                                &
               exner,                                                                              &
               hyp,                                                                                &
               pt,                                                                                 &
               q,                                                                                  &
               ql,                                                                                 &
               rdf_sc,                                                                             &
               rho_air_zw,                                                                         &
               tend,                                                                               &
               zu

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model

    USE chem_gasphase_mod,                                                                         &
        ONLY:  atol,                                                                               &
               chem_gasphase_integrate,                                                            &
               cs_mech,                                                                            &
               get_mechanism_name,                                                                 &
               nkppctrl,                                                                           &
               nmaxfixsteps,                                                                       &
               nphot,                                                                              &
               nreact,                                                                             &
               nspec,                                                                              &
               nvar,                                                                               &
               phot_names,                                                                         &
               rtol,                                                                               &
               spc_names,                                                                          &
               t_steps,                                                                            &
               vl_dim

    USE chem_modules

    USE chem_photolysis_mod,                                                                       &
        ONLY:  photolysis_control

    USE control_parameters,                                                                        &
        ONLY:  advanced_div_correction,                                                            &
               allow_negative_scalar_values,                                                       &
               air_chemistry,                                                                      &
               bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_radiation_l,                                                                     &
               bc_lr_cyc,                                                                          &
               bc_ns_cyc,                                                                          &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               cyclic_fill_initialization,                                                         &
               debug_output,                                                                       &
               dt_3d,                                                                              &
               humidity,                                                                           &
               initializing_actions,                                                               &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               max_pr_user,                                                                        &
               message_string,                                                                     &
               monotonic_limiter_z,                                                                &
               nesting_offline,                                                                    &
               omega,                                                                              &
               restart_data_format_output,                                                         &
               scalar_advec,                                                                       &
               time_since_reference_point,                                                         &
               timestep_scheme,                                                                    &
               tsc,                                                                                &
               use_prescribed_profile_data,                                                        &
               ws_scheme_sca

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE diffusion_s_mod,                                                                           &
        ONLY:  diffusion_s

    USE indices,                                                                                   &
        ONLY:  advc_flags_s,                                                                       &
               nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_flags

    USE kinds

    USE pegrid,                                                                                    &
        ONLY: myid,                                                                                &
              threads_per_task

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    USE radiation_model_mod,                                                                       &
        ONLY:  cos_zenith

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rrd_mpi_io,                                                                         &
               rd_mpi_io_check_array,                                                              &
               wrd_mpi_io

    USE statistics

    USE surface_mod,                                                                               &
        ONLY:  bc_hv,                                                                              &
               ind_pav_green,                                                                      &
               ind_veg_wall,                                                                       &
               ind_wat_win,                                                                        &
               surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_type,                                                                          &
               surf_top,                                                                           &
               surf_usm


    IMPLICIT NONE

    PRIVATE
    SAVE

    INTEGER, DIMENSION(nkppctrl)                           ::  icntrl       !< 20 integer parameters for fine tuning KPP code

    REAL(KIND=wp), PUBLIC ::  cs_time_step = 0.0_wp

    REAL(KIND=wp), DIMENSION(nkppctrl)                     ::  rcntrl       !< 20 real parameters for fine tuning of KPP code
                                                                            !< (e.g starting internal timestep of solver)

    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  spec_conc_1  !< pointer for swapping of timelevels for conc
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  spec_conc_2  !< pointer for swapping of timelevels for conc
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  spec_conc_3  !< pointer for swapping of timelevels for conc
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  freq_1       !< pointer for phtolysis frequncies
                                                                            !< (only 1 timelevel required) (e.g. solver type)

!
!-- Parameter needed for Deposition calculation using DEPAC model (van Zanten et al., 2010)
    INTEGER(iwp), PARAMETER ::  nlu_dep = 15            !< Number of DEPAC landuse classes (lu's)
    INTEGER(iwp), PARAMETER ::  ncmp = 10               !< Number of DEPAC gas components
    INTEGER(iwp), PARAMETER ::  nposp = 69              !< Number of possible species for deposition
!
!-- DEPAC landuse classes as defined in LOTOS-EUROS model v2.1
    INTEGER(iwp) ::  ilu_grass              = 1
    INTEGER(iwp) ::  ilu_arable             = 2
    INTEGER(iwp) ::  ilu_permanent_crops    = 3
    INTEGER(iwp) ::  ilu_coniferous_forest  = 4
    INTEGER(iwp) ::  ilu_deciduous_forest   = 5
    INTEGER(iwp) ::  ilu_water_sea          = 6
    INTEGER(iwp) ::  ilu_urban              = 7
    INTEGER(iwp) ::  ilu_other              = 8
    INTEGER(iwp) ::  ilu_desert             = 9
    INTEGER(iwp) ::  ilu_ice                = 10
    INTEGER(iwp) ::  ilu_savanna            = 11
    INTEGER(iwp) ::  ilu_tropical_forest    = 12
    INTEGER(iwp) ::  ilu_water_inland       = 13
    INTEGER(iwp) ::  ilu_mediterrean_scrub  = 14
    INTEGER(iwp) ::  ilu_semi_natural_veg   = 15

!
!-- NH3/SO2 ratio regimes:
    INTEGER(iwp), PARAMETER ::  iratns_low      = 1       !< low ratio NH3/SO2
    INTEGER(iwp), PARAMETER ::  iratns_high     = 2       !< high ratio NH3/SO2
    INTEGER(iwp), PARAMETER ::  iratns_very_low = 3       !< very low ratio NH3/SO2
!
!-- Default:
    INTEGER, PARAMETER ::  iratns_default = iratns_low
!
!-- Set alpha for f_light (4.57 is conversion factor from 1./(mumol m-2 s-1) to W m-2
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  alpha = (/                                         &
                     0.009_wp, 0.009_wp, 0.009_wp, 0.006_wp, 0.006_wp, -999.0_wp, -999.0_wp,       &
                     0.009_wp, -999.0_wp, -999.0_wp, 0.009_wp, 0.006_wp, -999.0_wp, 0.009_wp,      &
                     0.008_wp                          /) * 4.57_wp
!
!-- Set temperatures per land use for f_temp
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  tmin = (/                                          &
                     12.0_wp,  12.0_wp,  12.0_wp,  0.0_wp,  0.0_wp, -999.0_wp, -999.0_wp, 12.0_wp, &
                     -999.0_wp, -999.0_wp, 12.0_wp,  0.0_wp, -999.0_wp, 12.0_wp,  8.0_wp           &
                                                      /)
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  topt = (/                                          &
                     26.0_wp, 26.0_wp,  26.0_wp, 18.0_wp, 20.0_wp, -999.0_wp, -999.0_wp, 26.0_wp,  &
                     -999.0_wp, -999.0_wp, 26.0_wp, 20.0_wp, -999.0_wp, 26.0_wp, 24.0_wp           &
                                                      /)
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  tmax = (/                                          &
                     40.0_wp, 40.0_wp,  40.0_wp, 36.0_wp,   35.0_wp, -999.0_wp, -999.0_wp, 40.0_wp,&
                     -999.0_wp, -999.0_wp,  40.0_wp, 35.0_wp, -999.0_wp, 40.0_wp, 39.0_wp          &
                                                      /)
!
!-- Set f_min:
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  f_min = (/                                         &
                     0.01_wp, 0.01_wp, 0.01_wp, 0.1_wp, 0.1_wp, -999.0_wp, -999.0_wp, 0.01_wp,     &
                     -999.0_wp, -999.0_wp, 0.01_wp, 0.1_wp, -999.0_wp, 0.01_wp, 0.04_wp            &
                                                       /)

!
!-- Set maximal conductance (m/s)
!-- (R T/P) = 1/41000 mmol/m3 is given for 20 deg C to go from  mmol O3/m2/s to m/s
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  g_max = (/                                         &
                     270.0_wp, 300.0_wp, 300.0_wp, 140.0_wp, 150.0_wp, -999.0_wp, -999.0_wp,       &
                     270.0_wp, -999.0_wp, -999.0_wp, 270.0_wp, 150.0_wp, -999.0_wp, 300.0_wp,      &
                     422.0_wp                          /) / 41000.0_wp
!
!-- Set max, min for vapour pressure deficit vpd
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  vpd_max = (/                                       &
                     1.3_wp, 0.9_wp, 0.9_wp, 0.5_wp, 1.0_wp, -999.0_wp, -999.0_wp, 1.3_wp,         &
                     -999.0_wp, -999.0_wp, 1.3_wp, 1.0_wp, -999.0_wp, 0.9_wp, 2.8_wp               &
                                                         /)
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  vpd_min = (/                                       &
                     3.0_wp, 2.8_wp, 2.8_wp, 3.0_wp, 3.25_wp, -999.0_wp, -999.0_wp, 3.0_wp,        &
                     -999.0_wp, -999.0_wp, 3.0_wp,  3.25_wp, -999.0_wp, 2.8_wp, 4.5_wp             &
                                                         /)

!
!-- Soil resistance (numbers matched with lu_classes and component numbers)
    !     grs          ara        crp       cnf         dec         wat        urb      oth
    !     des          ice        sav       trf         wai         med        sem
    REAL(wp), PARAMETER ::  rsoil(nlu_dep,ncmp) = RESHAPE( (/                                      &
         1000.0_wp,  200.0_wp,  200.0_wp,  200.0_wp,  200.0_wp, 2000.0_wp,  400.0_wp, 1000.0_wp,   &
         2000.0_wp, 2000.0_wp, 1000.0_wp,  200.0_wp, 2000.0_wp,  200.0_wp,  400.0_wp,              &    !< O3
         1000.0_wp, 1000.0_wp, 1000.0_wp, 1000.0_wp, 1000.0_wp,   10.0_wp, 1000.0_wp, 1000.0_wp,   &
         1000.0_wp,  500.0_wp, 1000.0_wp, 1000.0_wp,   10.0_wp, 1000.0_wp, 1000.0_wp,              &    !< SO2
         1000.0_wp, 1000.0_wp, 1000.0_wp, 1000.0_wp, 1000.0_wp, 2000.0_wp, 1000.0_wp, 1000.0_wp,   &
         1000.0_wp, 2000.0_wp, 1000.0_wp, 1000.0_wp, 2000.0_wp, 1000.0_wp, 1000.0_wp,              &    !< NO2
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, 2000.0_wp, 1000.0_wp, -999.0_wp,   &
         2000.0_wp, 2000.0_wp, -999.0_wp, -999.0_wp, 2000.0_wp, -999.0_wp, -999.0_wp,              &    !< NO
          100.0_wp,  100.0_wp,  100.0_wp,  100.0_wp,  100.0_wp,   10.0_wp,  100.0_wp,  100.0_wp,   &
          100.0_wp, 1000.0_wp,  100.0_wp,  100.0_wp,   10.0_wp,  100.0_wp,  100.0_wp,              &    !< NH3
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, 2000.0_wp, 1000.0_wp, -999.0_wp,   &
         2000.0_wp, 2000.0_wp, -999.0_wp, -999.0_wp, 2000.0_wp, -999.0_wp, -999.0_wp,              &    !< CO
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp,   &
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp,              &    !< NO3
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp,   &
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp,              &    !< HNO3
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp,   &
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp,              &    !< N2O5
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp,   &
         -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp, -999.0_wp /),           &    !< H2O2
                                                                                 (/nlu_dep,ncmp/) )

    REAL(wp), PARAMETER ::  rsoil_frozen(ncmp) = (/ 2000.0_wp,                                     & !o3
                                                    500.0_wp,                                      & !so2
                                                    2000.0_wp,                                     & !no2
                                                    -999.0_wp,                                     & !no
                                                    1000.0_wp,                                     & !nh3
                                                    -999.0_wp,                                     & !co
                                                    -999.0_wp,                                     & !no3
                                                    -999.0_wp,                                     & !hno3
                                                    -999.0_wp,                                     & !n2o5
                                                    -999.0_wp /)                                     !h2o2
    REAL(wp), PARAMETER ::  rsoil_wet(ncmp)    = (/ 2000.0_wp,                                     & !o3
                                                    10.0_wp ,                                      & !so2
                                                    2000.0_wp,                                     & !no2
                                                    -999.0_wp,                                     & !no
                                                    10.0_wp  ,                                     & !nh3
                                                    -999.0_wp,                                     & !co
                                                    -999.0_wp,                                     & !no3
                                                    -999.0_wp,                                     & !hno3
                                                    -999.0_wp,                                     & !n2o5
                                                    -999.0_wp /)                                     !h2o2

    PUBLIC nreact
    PUBLIC nspec               !< number of gas phase chemical species including constant compound (e.g. N2)
    PUBLIC nvar                !< number of variable gas phase chemical species (nvar <= nspec)
    PUBLIC spc_names           !< names of gas phase chemical species (come from KPP) (come from KPP)
    PUBLIC spec_conc_2
!
!-- Interface section
    INTERFACE chem_actions
       MODULE PROCEDURE chem_actions
       MODULE PROCEDURE chem_actions_ij
    END INTERFACE chem_actions

    INTERFACE chem_3d_data_averaging
       MODULE PROCEDURE chem_3d_data_averaging
    END INTERFACE chem_3d_data_averaging

    INTERFACE chem_boundary_conditions
       MODULE PROCEDURE chem_boundary_conditions
    END INTERFACE chem_boundary_conditions

    INTERFACE chem_check_data_output
       MODULE PROCEDURE chem_check_data_output
    END INTERFACE chem_check_data_output

    INTERFACE chem_data_output_2d
       MODULE PROCEDURE chem_data_output_2d
    END INTERFACE chem_data_output_2d

    INTERFACE chem_data_output_3d
       MODULE PROCEDURE chem_data_output_3d
    END INTERFACE chem_data_output_3d

    INTERFACE chem_data_output_mask
       MODULE PROCEDURE chem_data_output_mask
    END INTERFACE chem_data_output_mask

    INTERFACE chem_check_data_output_pr
       MODULE PROCEDURE chem_check_data_output_pr
    END INTERFACE chem_check_data_output_pr

    INTERFACE chem_check_parameters
       MODULE PROCEDURE chem_check_parameters
    END INTERFACE chem_check_parameters

    INTERFACE chem_define_netcdf_grid
       MODULE PROCEDURE chem_define_netcdf_grid
    END INTERFACE chem_define_netcdf_grid

    INTERFACE chem_header
       MODULE PROCEDURE chem_header
    END INTERFACE chem_header

    INTERFACE chem_init_arrays
       MODULE PROCEDURE chem_init_arrays
    END INTERFACE chem_init_arrays

    INTERFACE chem_init
       MODULE PROCEDURE chem_init
    END INTERFACE chem_init

    INTERFACE chem_init_profiles
       MODULE PROCEDURE chem_init_profiles
    END INTERFACE chem_init_profiles

    INTERFACE chem_integrate
       MODULE PROCEDURE chem_integrate
       MODULE PROCEDURE chem_integrate_ij
    END INTERFACE chem_integrate

    INTERFACE chem_parin
       MODULE PROCEDURE chem_parin
    END INTERFACE chem_parin

    INTERFACE chem_non_advective_processes
       MODULE PROCEDURE chem_non_advective_processes
       MODULE PROCEDURE chem_non_advective_processes_ij
    END INTERFACE chem_non_advective_processes

    INTERFACE chem_exchange_horiz_bounds
       MODULE PROCEDURE chem_exchange_horiz_bounds
    END INTERFACE chem_exchange_horiz_bounds

    INTERFACE chem_prognostic_equations
       MODULE PROCEDURE chem_prognostic_equations
       MODULE PROCEDURE chem_prognostic_equations_ij
    END INTERFACE chem_prognostic_equations

    INTERFACE chem_rrd_local
       MODULE PROCEDURE chem_rrd_local_ftn
       MODULE PROCEDURE chem_rrd_local_mpi
    END INTERFACE chem_rrd_local

    INTERFACE chem_statistics
       MODULE PROCEDURE chem_statistics
    END INTERFACE chem_statistics

    INTERFACE chem_swap_timelevel
       MODULE PROCEDURE chem_swap_timelevel
    END INTERFACE chem_swap_timelevel

    INTERFACE chem_wrd_local
       MODULE PROCEDURE chem_wrd_local
    END INTERFACE chem_wrd_local

    INTERFACE chem_depo
       MODULE PROCEDURE chem_depo
       MODULE PROCEDURE chem_depo_ij
    END INTERFACE chem_depo

    PUBLIC chem_3d_data_averaging,                                                                 &
           chem_actions,                                                                           &
           chem_boundary_conditions,                                                               &
           chem_check_data_output,                                                                 &
           chem_check_data_output_pr,                                                              &
           chem_check_parameters,                                                                  &
           chem_data_output_2d,                                                                    &
           chem_data_output_3d,                                                                    &
           chem_data_output_mask,                                                                  &
           chem_define_netcdf_grid,                                                                &
           chem_header,                                                                            &
           chem_init,                                                                              &
           chem_init_arrays,                                                                       &
           chem_init_profiles,                                                                     &
           chem_integrate,                                                                         &
           chem_parin,                                                                             &
           chem_prognostic_equations,                                                              &
           chem_rrd_local,                                                                         &
           chem_statistics,                                                                        &
           chem_swap_timelevel,                                                                    &
           chem_wrd_local,                                                                         &
           chem_depo,                                                                              &
           chem_non_advective_processes,                                                           &
           chem_exchange_horiz_bounds

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for averaging 3D data of chemical species. Due to the fact that the averaged chem
!> arrays are allocated in chem_init, no if-query concerning the allocation is required (in any
!> mode). Attention: If you just specify an averaged output quantity in the _p3dr file during
!> restarts the first output includes the time between the beginning of the restart run and the
!> first output time (not necessarily the whole averaging_interval you have specified in your
!> _p3d/_p3dr file ).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_3d_data_averaging( mode, variable )

    USE control_parameters

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz_2d


    CHARACTER (LEN=*) ::  mode     !<
    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  i                  !< grid index x direction
    INTEGER(iwp) ::  j                  !< grid index y direction
    INTEGER(iwp) ::  k                  !< grid index z direction
    INTEGER(iwp) ::  lsp                !< running index for chem spcs
    INTEGER(iwp) ::  m                  !< running index surface type

    IF ( ( variable(1:3) == 'kc_'  .OR.  variable(1:3) == 'em_' )  )  THEN

       IF ( mode == 'allocate' )  THEN

          DO  lsp = 1, nspec
             IF ( TRIM( variable(1:3) ) == 'kc_'  .AND.                                            &
                  TRIM( variable(4:) )  == TRIM( chem_species(lsp)%name ) )  THEN
                IF ( .NOT. ALLOCATED( chem_species(lsp)%conc_av ) )  THEN
                   ALLOCATE( chem_species(lsp)%conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   chem_species(lsp)%conc_av = 0.0_wp
                ENDIF
             ENDIF
          ENDDO

       ELSEIF ( mode == 'sum' )  THEN

          DO  lsp = 1, nspec
             IF ( TRIM( variable(1:3) ) == 'kc_'  .AND.                                            &
                  TRIM( variable(4:) )  == TRIM( chem_species(lsp)%name ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         chem_species(lsp)%conc_av(k,j,i) = chem_species(lsp)%conc_av(k,j,i) +     &
                                                            chem_species(lsp)%conc(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ELSEIF ( TRIM( variable(4:) ) == TRIM( 'cssws*' ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  m = surf_def%start_index(j,i), surf_def%end_index(j,i)
                         chem_species(lsp)%cssws_av(j,i) = chem_species(lsp)%cssws_av(j,i) +       &
                                                           MERGE( surf_def%cssws(lsp,m), 0.0_wp,   &
                                                                  surf_def%upward(m) )
                      ENDDO
                      DO  m = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
                         chem_species(lsp)%cssws_av(j,i) = chem_species(lsp)%cssws_av(j,i) +       &
                                                           MERGE( surf_lsm%cssws(lsp,m), 0.0_wp,   &
                                                                  surf_lsm%upward(m) )
                      ENDDO
                      DO  m = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
                         chem_species(lsp)%cssws_av(j,i) = chem_species(lsp)%cssws_av(j,i) +       &
                                                           MERGE( surf_usm%cssws(lsp,m), 0.0_wp,   &
                                                                  surf_usm%upward(m) )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDDO

       ELSEIF ( mode == 'average' )  THEN

          DO  lsp = 1, nspec
             IF ( TRIM( variable(1:3) ) == 'kc_'  .AND.                                            &
                  TRIM( variable(4:) )  == TRIM( chem_species(lsp)%name ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         chem_species(lsp)%conc_av(k,j,i) = chem_species(lsp)%conc_av(k,j,i) /     &
                                                            REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO

             ELSEIF ( TRIM( variable(4:) ) == TRIM( 'cssws*' ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      chem_species(lsp)%cssws_av(j,i) = chem_species(lsp)%cssws_av(j,i) /          &
                                                        REAL( average_count_3d, KIND = wp )
                   ENDDO
                ENDDO
                CALL exchange_horiz_2d( chem_species(lsp)%cssws_av )
             ENDIF
          ENDDO
       ENDIF

    ENDIF

 END SUBROUTINE chem_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to initialize and set all boundary conditions for chemical species
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_boundary_conditions( horizontal_conditions_only )

    USE arrays_3d,                                                                                 &
        ONLY:  dzu

    INTEGER(iwp) ::  i                            !< grid index x direction.
    INTEGER(iwp) ::  j                            !< grid index y direction.
    INTEGER(iwp) ::  k                            !< grid index z direction.
    INTEGER(iwp) ::  lsp                          !< running index for chem spcs.
    INTEGER(iwp) ::  m                            !< running index surface elements.

    LOGICAL, OPTIONAL ::  horizontal_conditions_only  !< switch to set horizontal bc only


    IF ( .NOT. PRESENT( horizontal_conditions_only ) )  THEN
!
!--    Boundary condtions for chemical species at horizontal walls
       DO  lsp = 1, nspec
!
!--       Surface conditions:
          IF ( ibc_cs_b == 0 )  THEN
!
!--          Dirichlet:
!--          Run loop over all non-natural and natural walls. Note, in wall-datatype the k,j,i
!--          coordinate belong to the atmospheric grid point, therefore, set s_p at k+koff,
!--          j+joff, i+ioff, respectively.
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_hv%ns
                i = bc_hv%i(m)
                j = bc_hv%j(m)
                k = bc_hv%k(m)
                chem_species(lsp)%conc_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) =        &
                            chem_species(lsp)%conc(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m))
             ENDDO

          ELSEIF ( ibc_cs_b == 1 )  THEN
!
!--          Neumann:
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_hv%ns
                i = bc_hv%i(m)
                j = bc_hv%j(m)
                k = bc_hv%k(m)
                chem_species(lsp)%conc_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) =        &
                                                                   chem_species(lsp)%conc_p(k,j,i)
             ENDDO

          ENDIF

       ENDDO

!
!--    Top boundary conditions for chemical species
       DO  lsp = 1, nspec
          IF ( ibc_cs_t == 0 )  THEN
             chem_species(lsp)%conc_p(nzt+1,:,:) = chem_species(lsp)%conc(nzt+1,:,:)
          ELSEIF ( ibc_cs_t == 1 )  THEN
             chem_species(lsp)%conc_p(nzt+1,:,:) = chem_species(lsp)%conc_p(nzt,:,:)
          ELSEIF ( ibc_cs_t == 2 )  THEN
             chem_species(lsp)%conc_p(nzt+1,:,:) = chem_species(lsp)%conc_p(nzt,:,:) +             &
                                                   bc_cs_t_val(lsp) * dzu(nzt+1)
          ENDIF
       ENDDO

!
!--    Lateral boundary conditions.
!--    Dirichlet conditions have been already set when chem_species concentration is initialized.
!--    The initially set value is not touched during time integration, hence, this boundary value
!--    remains at a constant value.
!--    Neumann conditions:
       IF ( bc_radiation_cs_s )  THEN
          DO  lsp = 1, nspec
             chem_species(lsp)%conc_p(:,nys-1,:) = chem_species(lsp)%conc_p(:,nys,:)
          ENDDO
       ENDIF
       IF ( bc_radiation_cs_n )  THEN
          DO  lsp = 1, nspec
             chem_species(lsp)%conc_p(:,nyn+1,:) = chem_species(lsp)%conc_p(:,nyn,:)
          ENDDO
       ENDIF
       IF ( bc_radiation_cs_l )  THEN
          DO  lsp = 1, nspec
             chem_species(lsp)%conc_p(:,:,nxl-1) = chem_species(lsp)%conc_p(:,:,nxl)
          ENDDO
       ENDIF
       IF ( bc_radiation_cs_r )  THEN
          DO  lsp = 1, nspec
             chem_species(lsp)%conc_p(:,:,nxr+1) = chem_species(lsp)%conc_p(:,:,nxr)
          ENDDO
       ENDIF

    ELSE
!
!--    Lateral Neumann booundary conditions for timelevel t.
!--    This branch is executed when routine is called after the non-advective processes / before the
!--    prognostic equations.
       IF ( bc_radiation_cs_s )  THEN
          DO  lsp = 1, nspec
             chem_species(lsp)%conc(:,nys-1,:) = chem_species(lsp)%conc(:,nys,:)
          ENDDO
       ENDIF
       IF ( bc_radiation_cs_n )  THEN
          DO  lsp = 1, nspec
             chem_species(lsp)%conc(:,nyn+1,:) = chem_species(lsp)%conc(:,nyn,:)
          ENDDO
       ENDIF
       IF ( bc_radiation_cs_l )  THEN
          DO  lsp = 1, nspec
             chem_species(lsp)%conc(:,:,nxl-1) = chem_species(lsp)%conc(:,:,nxl)
          ENDDO
       ENDIF
       IF ( bc_radiation_cs_r )  THEN
          DO  lsp = 1, nspec
             chem_species(lsp)%conc(:,:,nxr+1) = chem_species(lsp)%conc(:,:,nxr)
          ENDDO
       ENDIF

    ENDIF

 END SUBROUTINE chem_boundary_conditions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for checking data output for chemical species
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_check_data_output( var, unit, i, ilen, k )


    CHARACTER (LEN=*) ::  unit     !<
    CHARACTER (LEN=*) ::  var      !<

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  ilen
    INTEGER(iwp) ::  lsp
    INTEGER(iwp) ::  k

    CHARACTER(LEN=16)    ::  spec_name

!
!-- Next statement is to avoid compiler warnings about unused variables
    IF ( ( i + ilen + k ) > 0  .OR.  var(1:1) == ' ' )  CONTINUE

    unit = 'illegal'

    spec_name = TRIM( var(4:) )             !< var 1:3 is 'kc_' or 'em_'.

    IF ( TRIM( var(1:3) ) == 'em_' )  THEN
       DO  lsp=1,nspec
          IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name ) )  THEN
             unit = 'mol m-2 s-1'
          ENDIF
!
!--       It is possible to plant PM10 and PM25 into the gasphase chemistry code as passive species
!--       (e.g. 'passive' in GASPHASE_PREPROC/mechanisms): set unit to micrograms per m**3 for PM10
!--       and PM25 (PM2.5)
          IF (spec_name(1:2) == 'PM')  THEN
             unit = 'kg m-2 s-1'
          ENDIF
       ENDDO

    ELSE

       DO  lsp=1,nspec
          IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name ) )  THEN
             unit = 'ppm'
          ENDIF
!
!--       It is possible to plant PM10 and PM25 into the gasphase chemistry code as passive species
!--       (e.g. 'passive' in GASPHASE_PREPROC/mechanisms): set unit to kilograms per m**3 for PM10
!--       and PM25 (PM2.5)
          IF (spec_name(1:2) == 'PM')  THEN
            unit = 'kg m-3'
          ENDIF

          IF ( spec_name(1:3) == 'POL' ) unit = 'nConc m-3'

       ENDDO

       DO  lsp=1,nphot
          IF ( TRIM( spec_name ) == TRIM( phot_frequen(lsp)%name ) )  THEN
             unit = 'sec-1'
          ENDIF
       ENDDO
    ENDIF


    RETURN
 END SUBROUTINE chem_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for checking data output of profiles for chemistry model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_check_data_output_pr( variable, var_count, unit, dopr_unit )

    USE arrays_3d

    USE control_parameters,                                                                        &
        ONLY:  data_output_pr,                                                                     &
               message_string

    USE profil_parameter

    USE statistics


    CHARACTER (LEN=*)  ::  unit      !< unit
    CHARACTER (LEN=*)  ::  variable  !< variable name
    CHARACTER (LEN=*)  ::  dopr_unit !< unit
    CHARACTER (LEN=16) ::  spec_name !< species name extracted from output string

    INTEGER(iwp) ::  index_start     !< start index of the species name in data-output string
    INTEGER(iwp) ::  index_end       !< end index of the species name in data-output string
    INTEGER(iwp) ::  var_count       !< number of data-output quantity
    INTEGER(iwp) ::  lsp             !< running index over species

    SELECT CASE ( TRIM( variable(1:3 ) ) )

       CASE ( 'kc_' )

          IF ( .NOT. air_chemistry )  THEN
             message_string = 'data_output_pr = "' // TRIM( data_output_pr(var_count) ) //          &
                              '" not implemented for air_chemistry = .FALSE.'
             CALL message( 'chem_check_data_output_pr', 'CHM0001', 1, 2, 0, 6, 0 )

          ENDIF
!
!--       Output of total fluxes is not allowed to date.
          IF ( TRIM( variable(1:4) ) == 'kc_w' )  THEN
             IF ( TRIM( variable(1:5) ) /= 'kc_w*'  .AND.  TRIM( variable(1:5) ) /= 'kc_w"' )  THEN
                message_string = 'data_output_pr = "' // TRIM( data_output_pr(var_count) ) //      &
                                 '" currently not implemented'

                CALL message( 'chem_check_data_output_pr', 'CHM0002', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
!
!--       Check for profile output of first-order moments, i.e. variable(4:) equals a species name.
          spec_name = TRIM( variable(4:) )
          DO  lsp = 1, nspec
             IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name ) )  THEN
                cs_pr_count_sp                 = cs_pr_count_sp + 1
                cs_pr_index_sp(cs_pr_count_sp) = lsp
                dopr_index(var_count)          = pr_palm + cs_pr_count_sp +                        &
                                                 cs_pr_count_fl_sgs + cs_pr_count_fl_res
                dopr_unit                      = 'ppm'
                IF ( spec_name(1:2) == 'PM')  THEN
                   dopr_unit = 'kg m-3'
                ENDIF
                hom(:,2,dopr_index(var_count),:) = SPREAD( zu, 2, statistic_regions+1 )
                unit                              = dopr_unit

                hom_index_spec(cs_pr_count_sp)    = dopr_index(var_count)
             ENDIF
          ENDDO
!
!--       Check for profile output of fluxes. variable(index_start:index_end) equals a species name.
!--       Start with SGS components.
          IF ( TRIM( variable(1:5) ) == 'kc_w"' )  THEN
             DO  lsp = 1, nspec
                index_end   = LEN( TRIM( variable ) ) - 1
                index_start = 6
                spec_name = TRIM( variable(index_start:index_end) )
                IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name ) )  THEN
                   cs_pr_count_fl_sgs                     = cs_pr_count_fl_sgs + 1
                   cs_pr_index_fl_sgs(cs_pr_count_fl_sgs) = lsp
                   dopr_index(var_count)                  = pr_palm + cs_pr_count_sp +          &
                                                            cs_pr_count_fl_sgs +                &
                                                            cs_pr_count_fl_res
                   dopr_unit                              = 'm ppm s-1'
                   IF ( spec_name(1:2) == 'PM')  THEN
                      dopr_unit = 'kg m-2 s-1'
                   ENDIF
                   hom(:,2,dopr_index(var_count),:)     = SPREAD( zu, 2, statistic_regions+1 )
                   unit                                 = dopr_unit

                   hom_index_fl_sgs(cs_pr_count_fl_sgs) = dopr_index(var_count)
                ENDIF
             ENDDO
          ENDIF
!
!--       Proceed with resolved-scale fluxes.
          IF ( TRIM( variable(1:5) ) == 'kc_w*' )  THEN
             spec_name = TRIM( variable(6:) )
             DO  lsp = 1, nspec
                index_start = 6
                index_end   = LEN( TRIM( variable ) ) - 1
                spec_name = TRIM( variable(index_start:index_end) )
                IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name ) )  THEN
                   cs_pr_count_fl_res                     = cs_pr_count_fl_res + 1
                   cs_pr_index_fl_res(cs_pr_count_fl_res) = lsp
                   dopr_index(var_count)                  = pr_palm + cs_pr_count_sp +             &
                                                            cs_pr_count_fl_sgs +                   &
                                                            cs_pr_count_fl_res
                   dopr_unit                              = 'm ppm s-1'
                   IF ( spec_name(1:2) == 'PM')  THEN
                      dopr_unit = 'kg m-2 s-1'
                   ENDIF
                   hom(:,2, dopr_index(var_count),:)    = SPREAD( zu, 2, statistic_regions+1 )
                   unit                                 = dopr_unit

                   hom_index_fl_res(cs_pr_count_fl_res) = dopr_index(var_count)
                ENDIF
             ENDDO
          ENDIF
       CASE DEFAULT
          unit = 'illegal'

    END SELECT

 END SUBROUTINE chem_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for chemistry_model_mod
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_check_parameters

    USE control_parameters,                                                                        &
        ONLY:  bc_lr,                                                                              &
               bc_ns,                                                                              &
               data_output_pr

    USE radiation_model_mod,                                                                       &
        ONLY:  radiation_volumetric_flux


    INTEGER(iwp) ::  i             !< loop variable
    INTEGER(iwp) ::  lsp           !< running index for chem spcs.
    INTEGER(iwp) ::  lsp_usr       !< running index for user defined chem spcs

    LOGICAL  ::  found


    IF ( .NOT. radiation_volumetric_flux  .AND.  photolysis_shading )  THEN
       message_string = 'photolyis shading requires 3D volumetric radiation flux set to true'
       CALL message( 'chem_check_parameters', 'CHM0003', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check for chemistry time-step
    IF ( call_chem_at_all_substeps )  THEN
       message_string = 'call_chem_at_all_substeps should only be used for test purposes'
       CALL message( 'chem_check_parameters', 'CHM0004', 0, 1, 0, 6, 0 )
    ENDIF
!
!-- Check for photolysis scheme
    IF ( ( photolysis_scheme /= 'simple' ) .AND. ( photolysis_scheme /= 'constant' )  )  THEN
       message_string = 'unknown photolysis scheme = "' // TRIM( photolysis_scheme ) // '"'
       CALL message( 'chem_check_parameters', 'CHM0005', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check for chemical mechanism used
    CALL get_mechanism_name
    IF ( chem_mechanism /= TRIM( cs_mech ) )  THEN
       message_string = 'unknown chem_mechanism = "' // TRIM( chem_mechanism ) // '"'
       CALL message( 'chem_check_parameters', 'CHM0006', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check bottom boundary condition and set internal steering parameter
    IF ( bc_cs_b == 'dirichlet' )  THEN
       ibc_cs_b = 0
    ELSEIF ( bc_cs_b == 'neumann' )  THEN
       ibc_cs_b = 1
    ELSE
       message_string = 'unknown boundary condition: bc_cs_b ="' // TRIM( bc_cs_b ) // '"'
       CALL message( 'chem_check_parameters', 'CHM0007', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check top boundary condition and set internal steering parameter
    IF ( bc_cs_t == 'dirichlet' )  THEN
       ibc_cs_t = 0
    ELSEIF ( bc_cs_t == 'neumann' )  THEN
       ibc_cs_t = 1
    ELSEIF ( bc_cs_t == 'initial_gradient' )  THEN
       ibc_cs_t = 2
    ELSEIF ( bc_cs_t == 'nested' )  THEN
       ibc_cs_t = 3
    ELSE
       message_string = 'unknown boundary condition: bc_c_t ="' // TRIM( bc_cs_t ) // '"'
       CALL message( 'chem_check_parameters', 'CHM0007', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- If nesting_chem = .F., set top boundary condition to its default value
    IF ( .NOT. nesting_chem  .AND.  ibc_cs_t == 3  )  THEN
       ibc_cs_t = 2
       bc_cs_t = 'initial_gradient'
    ENDIF

!
!-- Check left and right boundary conditions. First set default value if not set by user.
    IF ( bc_cs_l == 'undefined' )  THEN
       IF ( bc_lr == 'cyclic' )  THEN
          bc_cs_l = 'cyclic'
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          bc_cs_l = 'dirichlet'
       ELSEIF ( bc_lr == 'radiation/dirichlet' )  THEN
          bc_cs_l = 'neumann'
       ENDIF
    ENDIF
    IF ( bc_cs_r == 'undefined' )  THEN
       IF ( bc_lr == 'cyclic' )  THEN
          bc_cs_r = 'cyclic'
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          bc_cs_r = 'neumann'
       ELSEIF ( bc_lr == 'radiation/dirichlet' )  THEN
          bc_cs_r = 'dirichlet'
       ENDIF
    ENDIF
    IF ( bc_cs_l /= 'dirichlet'  .AND.  bc_cs_l /= 'neumann'  .AND.  bc_cs_l /= 'cyclic' )  THEN
       message_string = 'unknown boundary condition: bc_cs_l = "' // TRIM( bc_cs_l ) // '"'
       CALL message( 'chem_check_parameters','CHM0007', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( bc_cs_r /= 'dirichlet'  .AND.  bc_cs_r /= 'neumann'  .AND.  bc_cs_r /= 'cyclic' )  THEN
       message_string = 'unknown boundary condition: bc_cs_r = "' // TRIM( bc_cs_r ) // '"'
       CALL message( 'chem_check_parameters','CHM0007', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check north and south boundary conditions. First set default value if not set by user.
    IF ( bc_cs_n == 'undefined' )  THEN
       IF ( bc_ns == 'cyclic' )  THEN
          bc_cs_n = 'cyclic'
       ELSEIF ( bc_ns == 'dirichlet/radiation' )  THEN
          bc_cs_n = 'dirichlet'
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          bc_cs_n = 'neumann'
       ENDIF
    ENDIF
    IF ( bc_cs_s == 'undefined' )  THEN
       IF ( bc_ns == 'cyclic' )  THEN
          bc_cs_s = 'cyclic'
       ELSEIF ( bc_ns == 'dirichlet/radiation' )  THEN
          bc_cs_s = 'neumann'
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          bc_cs_s = 'dirichlet'
       ENDIF
    ENDIF
    IF ( bc_cs_n /= 'dirichlet'  .AND.  bc_cs_n /= 'neumann'  .AND.  bc_cs_n /= 'cyclic' )  THEN
       message_string = 'unknown boundary condition: bc_cs_n = "' // TRIM( bc_cs_n ) // '"'
       CALL message( 'chem_check_parameters','CHM0007', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( bc_cs_s /= 'dirichlet'  .AND.  bc_cs_s /= 'neumann'  .AND.  bc_cs_s /= 'cyclic' )  THEN
       message_string = 'unknown boundary condition: bc_cs_s = "' // TRIM( bc_cs_s ) // '"'
       CALL message( 'chem_check_parameters','CHM0007', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Cyclic conditions must be set identically at opposing boundaries
    IF ( ( bc_cs_l == 'cyclic' .AND. bc_cs_r /= 'cyclic' )  .OR.                                   &
         ( bc_cs_r == 'cyclic' .AND. bc_cs_l /= 'cyclic' ) )  THEN
       message_string = 'boundary conditions bc_cs_l and bc_cs_r must both be cyclic or non-cyclic'
       CALL message( 'chem_check_parameters','CHM0008', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ( bc_cs_n == 'cyclic' .AND. bc_cs_s /= 'cyclic' )  .OR.                                   &
         ( bc_cs_s == 'cyclic' .AND. bc_cs_n /= 'cyclic' ) )  THEN
       message_string = 'boundary conditions bc_cs_n and bc_cs_s must both be cyclic or non-cyclic'
       CALL message( 'chem_check_parameters','CHM0008', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Set the switches that control application of horizontal boundary conditions at the boundaries
!-- of the total domain
    IF ( bc_cs_n == 'dirichlet'  .AND.  nyn == ny )  bc_dirichlet_cs_n = .TRUE.
    IF ( bc_cs_n == 'neumann'    .AND.  nyn == ny )  bc_radiation_cs_n = .TRUE.
    IF ( bc_cs_s == 'dirichlet'  .AND.  nys ==  0 )  bc_dirichlet_cs_s = .TRUE.
    IF ( bc_cs_s == 'neumann'    .AND.  nys ==  0 )  bc_radiation_cs_s = .TRUE.
    IF ( bc_cs_l == 'dirichlet'  .AND.  nxl ==  0 )  bc_dirichlet_cs_l = .TRUE.
    IF ( bc_cs_l == 'neumann'    .AND.  nxl ==  0 )  bc_radiation_cs_l = .TRUE.
    IF ( bc_cs_r == 'dirichlet'  .AND.  nxr == nx )  bc_dirichlet_cs_r = .TRUE.
    IF ( bc_cs_r == 'neumann'    .AND.  nxr == nx )  bc_radiation_cs_r = .TRUE.

!
!-- Set the communicator to be used for ghost layer data exchange
!-- 1: cyclic, 2: cyclic along x, 3: cyclic along y, 4: non-cyclic
    IF ( bc_cs_l == 'cyclic' )  THEN
       IF ( bc_cs_s == 'cyclic' )  THEN
          communicator_chem = 1
       ELSE
          communicator_chem = 2
       ENDIF
    ELSE
       IF ( bc_cs_s == 'cyclic' )  THEN
          communicator_chem = 3
       ELSE
          communicator_chem = 4
       ENDIF
    ENDIF

!
!-- chem_check_parameters is called before the array chem_species is allocated!
!-- temporary switch of this part of the check
!> TODO: this workaround definitely needs to be removed from here!!!
    CALL chem_init_internal
!
!-- Check for initial chem species input
    lsp_usr = 1
    lsp     = 1
    DO WHILE ( cs_name (lsp_usr) /= 'novalue')
       found = .FALSE.
       DO  lsp = 1, nvar
          IF ( TRIM( cs_name (lsp_usr) ) == TRIM( chem_species(lsp)%name) )  THEN
             found = .TRUE.
             EXIT
          ENDIF
       ENDDO
       IF ( .NOT.  found )  THEN
          message_string = 'unused/incorrect input for initial surface value = "' //               &
                            TRIM( cs_name(lsp_usr) ) // '"'
          CALL message( 'chem_check_parameters', 'CHM0009', 1, 2, 0, 6, 0 )
       ENDIF
       lsp_usr = lsp_usr + 1
    ENDDO
!
!-- Check for surface  emission flux chem species
    lsp_usr = 1
    lsp     = 1
    DO WHILE ( surface_csflux_name (lsp_usr) /= 'novalue')
       found = .FALSE.
       DO  lsp = 1, nvar
          IF ( TRIM( surface_csflux_name (lsp_usr) ) == TRIM( chem_species(lsp)%name ) )  THEN
             found = .TRUE.
             EXIT
          ENDIF
       ENDDO
       IF ( .NOT.  found )  THEN
          message_string = 'unused/incorrect input of chemical species for surface emission ' //   &
                           'fluxes = "' // TRIM( surface_csflux_name(lsp_usr) ) // '"'
          CALL message( 'chem_check_parameters', 'CHM0010', 1, 2, 0, 6, 0 )
       ENDIF
       lsp_usr = lsp_usr + 1
    ENDDO
!
!-- Determine the number of chemistry profiles and append them to the standard data output.
    i = 1

    DO  WHILE ( data_output_pr(i)  /= ' '  .AND.  i <= SIZE( data_output_pr ) )
       IF ( TRIM( data_output_pr(i)(1:3) ) == 'kc_' )  THEN
          max_pr_cs = max_pr_cs + 1
       ENDIF
       i = i + 1
    ENDDO

 END SUBROUTINE chem_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 2D output variables for chemical species
!> @todo: Remove "mode" from argument list, not used.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )


    CHARACTER (LEN=*) ::  grid       !<
    CHARACTER (LEN=*) ::  mode       !<
    CHARACTER (LEN=*) ::  variable   !<

    INTEGER(iwp) ::  av              !< flag to control data output of instantaneous or
                                     !< time-averaged data
    INTEGER(iwp) ::  nzb_do          !< lower limit of the domain (usually nzb)
    INTEGER(iwp) ::  nzt_do          !< upper limit of the domain (usually nzt+1)

    LOGICAL      ::  found           !<
    LOGICAL      ::  two_d           !< flag parameter that indicates 2D variables (horizontal cross
                                     !< sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb:nzt+1) ::  local_pf

!
!-- Local variables.
    CHARACTER(LEN=16)    ::  spec_name
    INTEGER(iwp) ::  lsp
    INTEGER(iwp) ::  i               !< grid index along x-direction
    INTEGER(iwp) ::  j               !< grid index along y-direction
    INTEGER(iwp) ::  k               !< grid index along z-direction
    INTEGER(iwp) ::  m               !< running indices for surfaces
    INTEGER(iwp) ::  char_len        !< length of a character string


!
!-- Next statement is to avoid compiler warnings about unused variables
    IF ( mode(1:1) == ' '  .OR.  two_d )  CONTINUE

    found = .FALSE.
    char_len  = LEN_TRIM( variable )

    spec_name = TRIM( variable(4:char_len-3) )
!
!-- Output of emission values, i.e. surface fluxes cssws.
    IF ( variable(1:3) == 'em_' )  THEN

       local_pf = 0.0_wp

       DO  lsp = 1, nvar
          IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name) )  THEN
!
!--          No average output for now.
             DO  m = 1, surf_lsm%ns
                local_pf(surf_lsm%i(m),surf_lsm%j(m),nzb+1) =                                  &
                                           MERGE( surf_lsm%cssws(lsp,m),                       &
                                                  local_pf(surf_lsm%i(m),surf_lsm%j(m),nzb+1), &
                                                  surf_lsm%upward(m) )
             ENDDO
             DO  m = 1, surf_usm%ns
                local_pf(surf_usm%i(m),surf_usm%j(m),nzb+1) =                                  &
                                           MERGE( surf_usm%cssws(lsp,m),                       &
                                                  local_pf(surf_usm%i(m),surf_usm%j(m),nzb+1), &
                                                  surf_usm%upward(m) )
             ENDDO
             grid = 'zu'
             found = .TRUE.
          ENDIF
       ENDDO

    ELSE

       DO  lsp=1,nspec
          IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name )  .AND.                          &
                ( (variable(char_len-2:) == '_xy')  .OR.                                           &
                  (variable(char_len-2:) == '_xz')  .OR.                                           &
                  (variable(char_len-2:) == '_yz') ) )  THEN
             IF (av == 0)  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                            local_pf(i,j,k) = chem_species(lsp)%conc(k,j,i)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO

             ELSE
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb_do, nzt_do
                         IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                            local_pf(i,j,k) = chem_species(lsp)%conc_av(k,j,i)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
             grid = 'zu'
             found = .TRUE.
          ENDIF
       ENDDO
    ENDIF

    RETURN

 END SUBROUTINE chem_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 3D output variables for chemical species
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    CHARACTER (LEN=*)    ::  variable     !<

    INTEGER(iwp)         ::  av         !<
    INTEGER(iwp) ::  nzb_do             !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do             !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found                !<

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf
!
!-- Local variables
    CHARACTER(LEN=16)    ::  spec_name

    INTEGER(iwp)         ::  i
    INTEGER(iwp)         ::  j
    INTEGER(iwp)         ::  k
    INTEGER(iwp)         ::  m       !< running indices for surfaces
    INTEGER(iwp)         ::  lsp     !< running index for chem spcs


    found = .FALSE.
    IF ( .NOT. (variable(1:3) == 'kc_' .OR. variable(1:3) == 'em_' ) )  THEN
       RETURN
    ENDIF

    spec_name = TRIM( variable(4:) )

    IF ( variable(1:3) == 'em_' )  THEN

       DO  lsp = 1, nvar   !!! cssws - nvar species, chem_species - nspec species !!!
          IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name) )  THEN

             local_pf = 0.0_wp
!
!--          no average for now
             DO  m = 1, surf_usm%ns
                local_pf(surf_usm%i(m),surf_usm%j(m),surf_usm%k(m)) =                              &
                        local_pf(surf_usm%i(m),surf_usm%j(m),surf_usm%k(m)) + surf_usm%cssws(lsp,m)
             ENDDO
             DO  m = 1, surf_lsm%ns
                local_pf(surf_lsm%i(m),surf_lsm%j(m),surf_lsm%k(m)) =                              &
                        local_pf(surf_lsm%i(m),surf_lsm%j(m),surf_lsm%k(m)) + surf_lsm%cssws(lsp,m)
             ENDDO
             found = .TRUE.
          ENDIF
       ENDDO
    ELSE
      DO  lsp = 1, nspec
         IF ( TRIM( spec_name ) == TRIM( chem_species(lsp)%name) )  THEN
            IF (av == 0)  THEN
               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     DO  k = nzb_do, nzt_do
                        IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                           local_pf(i,j,k) = chem_species(lsp)%conc(k,j,i)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ELSE

               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     DO  k = nzb_do, nzt_do
                        IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                           local_pf(i,j,k) = chem_species(lsp)%conc_av(k,j,i)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
            found = .TRUE.
         ENDIF
      ENDDO

       DO  lsp=1,nphot
          IF ( TRIM( spec_name ) == TRIM( phot_frequen(lsp)%name ) )  THEN
            IF (av == 0)  THEN
               DO  i = nxl, nxr
                  DO  j = nys, nyn
                     DO  k = nzb_do, nzt_do
                        IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                           local_pf(i,j,k) = phot_frequen(lsp)%freq(k,j,i)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
            found = .TRUE.
          ENDIF
       ENDDO
    ENDIF

    RETURN

 END SUBROUTINE chem_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining mask output variables for chemical species
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_data_output_mask( av, variable, found, local_pf, mid )

    USE control_parameters

    CHARACTER(LEN=16) ::  spec_name
    CHARACTER(LEN=*)  ::  variable    !<

    INTEGER(iwp) ::  av              !< flag to control data output of instantaneous or
                                     !< time-averaged data
    INTEGER(iwp) ::  i               !< grid index along x-direction
    INTEGER(iwp) ::  im              !< loop index for masked variables
    INTEGER(iwp) ::  j               !< grid index along y-direction
    INTEGER(iwp) ::  jm              !< loop index for masked variables
    INTEGER(iwp) ::  k               !< grid index along z-direction
    INTEGER(iwp) ::  kk              !< masked output index along z-direction
    INTEGER(iwp) ::  ktt             !< k index of lowest non-terrain grid point
    INTEGER(iwp) ::  lsp
    INTEGER(iwp) ::  mid             !< masked output running index

    LOGICAL ::  found

    REAL(wp), DIMENSION(mask_size_l(mid,1),mask_size_l(mid,2),mask_size_l(mid,3)) ::  local_pf   !<


!
!-- Local variables.

    spec_name = TRIM( variable(4:) )
    found = .FALSE.

    DO  lsp=1,nspec
       IF (TRIM( spec_name ) == TRIM( chem_species(lsp)%name) )  THEN
          IF (av == 0)  THEN
             IF ( .NOT. mask_surface(mid) )  THEN

                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k = 1, mask_size(mid,3)
                          local_pf(i,j,k) = chem_species(lsp)%conc( mask_k(mid,k),                 &
                                                                    mask_j(mid,j),                 &
                                                                    mask_i(mid,i)      )
                      ENDDO
                   ENDDO
                ENDDO

             ELSE
!
!--             Terrain-following masked output
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
!
!--                   Get k index of the lowest non-terrain grid point
                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)
                      ktt = MINLOC( MERGE( 1, 0, BTEST( topo_flags(:,jm,im), 5 ) ),                &
                                    DIM = 1 ) - 1
                      DO  k = 1, mask_size_l(mid,3)
                         kk = MIN( ktt + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
                            local_pf(i,j,k) = chem_species(lsp)%conc(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO

             ENDIF
          ELSE
             IF ( .NOT. mask_surface(mid) )  THEN

                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
                      DO  k =  1, mask_size_l(mid,3)
                          local_pf(i,j,k) = chem_species(lsp)%conc_av(  mask_k(mid,k),             &
                                                                        mask_j(mid,j),             &
                                                                        mask_i(mid,i)         )
                      ENDDO
                   ENDDO
                ENDDO

             ELSE
!
!--             Terrain-following masked output
                DO  i = 1, mask_size_l(mid,1)
                   DO  j = 1, mask_size_l(mid,2)
!
!--                   Get k index of the lowest non-terrain grid point
                      im = mask_i(mid,i)
                      jm = mask_j(mid,j)
                      ktt = MINLOC( MERGE( 1, 0, BTEST( topo_flags(:,jm,im), 5 )),                 &
                                    DIM = 1 ) - 1
                      DO  k = 1, mask_size_l(mid,3)
                         kk = MIN( ktt + mask_k(mid,k) - 1, nzt+1 )
!
!--                      Set value if not in building.
                         IF ( .NOT. BTEST( topo_flags(kk,jm,im), 6 ) )  THEN
                            local_pf(i,j,k) = chem_species(lsp)%conc_av(kk,jm,im)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO

             ENDIF

          ENDIF
          found = .TRUE.
          EXIT
       ENDIF
    ENDDO

    RETURN

 END SUBROUTINE chem_data_output_mask


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )


    CHARACTER (LEN=*), INTENT(IN)  ::  var          !<

    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x       !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y       !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z       !<

    LOGICAL, INTENT(OUT)           ::  found        !<

    found  = .TRUE.

    IF ( var(1:3) == 'kc_' .OR. var(1:3) == 'em_' )  THEN                !< always the same grid for
                                                                         !< chemistry variables
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'zu'
    ELSE
       found  = .FALSE.
       grid_x = 'none'
       grid_y = 'none'
       grid_z = 'none'
    ENDIF

 END SUBROUTINE chem_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining header output for chemistry model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_header( io )

    USE radiation_model_mod,                                                                       &
        ONLY:  radiation_volumetric_flux

    CHARACTER (LEN=80)  :: docsflux_chr
    CHARACTER (LEN=80)  :: docsinit_chr

    INTEGER(iwp), INTENT(IN) ::  io            !< Unit of the output file

    INTEGER(iwp)  :: cs_fixed
    INTEGER(iwp)  :: lsp                       !< running index for chem spcs

!
!-- Get name of chemical mechanism from chem_gasphase_mod
    CALL get_mechanism_name
!
!-- Write chemistry model  header
    WRITE( io, 1 )
!
!-- Gasphase reaction status
    IF ( chem_gasphase_on )  THEN
       WRITE( io, 2 )
    ELSE
       WRITE( io, 3 )
    ENDIF
!
!-- Emission mode info
    WRITE( io, 4 ) emiss_read_legacy_mode
!
!-- At the moment the evaluation is done with both emiss_lod and mode_emis but once salsa has been
!-- migrated to emiss_lod the .OR. mode_emis conditions can be removed
    IF     ( ( emiss_lod == 1 )  .OR.  ( mode_emis == 'DEFAULT' ) )        THEN
       WRITE( io, 5 )
    ELSEIF ( ( emiss_lod == 0 )  .OR.  ( mode_emis == 'PARAMETERIZED' ) )  THEN
       WRITE( io, 6 )
    ELSEIF ( ( emiss_lod == 2 )  .OR.  ( mode_emis == 'PRE-PROCESSED' ) )  THEN
       WRITE( io, 7 )
    ENDIF

    IF     (  emis_biogenic_lod == 0 )  THEN
       WRITE( io, 15 )
    ELSEIF (  emis_biogenic_lod == 1 )  THEN
       WRITE( io, 16 )
    ELSEIF (  emis_biogenic_lod == 2 )  THEN
       WRITE( io, 17 )
    ENDIF
!
!-- Photolysis scheme info
    IF ( photolysis_scheme == "simple" )  THEN
       WRITE( io, 8 )
    ELSEIF (photolysis_scheme == "constant" )  THEN
       WRITE( io, 9 )
    ENDIF

    IF  ( radiation_volumetric_flux .AND. photolysis_shading )  THEN
       WRITE( io, 20 )
    ELSEIF ( .NOT.  photolysis_shading ) THEN
       WRITE( io, 21 )
    ENDIF



!
!-- ISORROPIA
    IF ( chem_isorropia )  THEN
       WRITE( io, * ) '   --> ISORROPIA coupling activated'
       WRITE( io, * ) '       Problem type                                CNTRL(1) : ',            &
                      chem_isorropia_problem_type
       WRITE( io, * ) '       Aerosol state                               CNTRL(2) : ',            &
                      chem_isorropia_aerosol_state
       WRITE( io, * ) '       MDR weighting method                        WFTYPI   : ',            &
                      chem_isorropia_mdr_weight_method
       WRITE( io, * ) '       Activity coefficient algorithm              IACALCI  : ',            &
                      chem_isorropia_activity_coefficient_method
       WRITE( io, * ) '       ISORROPIA solver convergence                EPSI     : ',            &
                      chem_isorropia_solver_tolerance
       WRITE( io, * ) '       ISORROPIA solver maximum iterations         MAXTI    : ',            &
                      chem_isorropia_max_iteration
       WRITE( io, * ) '       Activity coefficient solver sweeps          NSWEEPI  : ',            &
                      chem_isorropia_max_activity_sweep
       WRITE( io, * ) '       Activity coefficient solver tolerance       EPSACTI  : ',            &
                      chem_isorropia_activity_tolerance
       WRITE( io, * ) '       Subdivisions for root tracking              NDIV     : ',            &
                      chem_isorropia_root_subdivisions
       WRITE( io, * ) '       Mass conservation mode (ISORROPIA II only)  NADJI    : ',            &
                      chem_isorropia_mass_conservation_mode
    ENDIF
!
!-- Emission flux info
    lsp = 1
    docsflux_chr ='Chemical species for surface emission flux: '
    DO WHILE ( surface_csflux_name(lsp) /= 'novalue' )
       docsflux_chr = TRIM( docsflux_chr ) // ' ' // TRIM( surface_csflux_name(lsp) ) // ','
       IF ( LEN_TRIM( docsflux_chr ) >= 75 )  THEN
          WRITE( io, 10 ) docsflux_chr
          docsflux_chr = '       '
       ENDIF
       lsp = lsp + 1
    ENDDO

    IF ( docsflux_chr /= '' )  THEN
       WRITE( io, 10 ) docsflux_chr
    ENDIF
!
!-- Initialization of Surface and profile chemical species
    lsp = 1
    docsinit_chr ='Chemical species for initial surface and profile emissions: '
    DO WHILE ( cs_name(lsp) /= 'novalue' )
       docsinit_chr = TRIM( docsinit_chr ) // ' ' // TRIM( cs_name(lsp) ) // ','
       IF ( LEN_TRIM( docsinit_chr ) >= 75 )  THEN
          WRITE( io, 11 ) docsinit_chr
          docsinit_chr = '       '
       ENDIF
       lsp = lsp + 1
    ENDDO

    IF ( docsinit_chr /= '' )  THEN
       WRITE( io, 11 ) docsinit_chr
    ENDIF

    IF ( nesting_chem )  WRITE( io, 12 ) nesting_chem
    IF ( nesting_offline_chem .AND. nesting_offline )  WRITE( io, 13 ) nesting_offline_chem

    WRITE( io, 14 ) TRIM( bc_cs_b ), TRIM( bc_cs_t ), TRIM( bc_cs_s ), TRIM( bc_cs_n ),            &
                    TRIM( bc_cs_l ), TRIM( bc_cs_r )

!
!-- Number of variable and fix chemical species and number of reactions
    cs_fixed = nspec - nvar
    WRITE( io, * ) '   --> Chemical Mechanism          : ', cs_mech
    WRITE( io, * ) '   --> Chemical species, variable  : ', nvar
    WRITE( io, * ) '   --> Chemical species, fixed     : ', cs_fixed
    WRITE( io, * ) '   --> Total number of reactions   : ', nreact
    WRITE( io, * ) '   --> Gas phase chemistry solver  : ', icntrl(3)
    WRITE( io, * ) '   --> Vector length (vector mode if > 1): ', vl_dim

    IF ( allow_negative_scalar_values )  THEN
       WRITE( io, 18 )
    ELSE
       WRITE( io, 19 )
    ENDIF


1   FORMAT (//' Chemistry model information:'/' ----------------------------'/)
2   FORMAT ('    --> Chemical reactions are turned on')
3   FORMAT ('    --> Chemical reactions are turned off')
4   FORMAT ('    --> Legacy emission read mode: ',L3,/,                                            &
            '        All emissions data will be loaded prior to start of simulation')
5   FORMAT ('    --> Emission mode = DEFAULT ')
6   FORMAT ('    --> Emission mode = PARAMETERIZED (LOD 0)')
7   FORMAT ('    --> Emission mode = PRE-PROCESSED (LOD 2)')
8   FORMAT ('    --> Photolysis scheme used =  simple ')
9   FORMAT ('    --> Photolysis scheme used =  constant ')
10  FORMAT (/'    ',A)
11  FORMAT (/'    ',A)
12  FORMAT (/'    Self nesting for chemistry variables (if nested_run): ', L1 )
13  FORMAT (/'    Offline nesting for chemistry variables : ', L1 )
14  FORMAT (/'    Boundary conditions for chemical species:', /                                     &
             '       bottom/top:   ',A10,' / ',A10, /                                               &
             '       north/south:  ',A10,' / ',A10, /                                               &
             '       left/right:   ',A10,' / ',A10)
15  FORMAT ('    --> Biogenic emission mode = DEFAULT ')
16  FORMAT ('    --> Biogenic emission mode = MEGAN ')
17  FORMAT ('    --> Biogenic emission mode = BAUME ')
18  FORMAT (/'     Negative values of chemical species due to dispersion errors are permitted!')
19  FORMAT (/'     Negative values of chemical species due to dispersion errors are cut, which',/  &
             '     which may appear in results as (small) artificial source of chemical species.')
20  FORMAT ('        3D radiation flux and photolyis shading are turned on')
21  FORMAT ('        Photolyis shading turned off')

 END SUBROUTINE chem_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine initializating chemistry_model_mod specific arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_init_arrays
!
!-- Please use this place to allocate required arrays

 END SUBROUTINE chem_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine initializating chemistry_model_mod.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_init

    USE chem_emis_biogenic_mod,                                                                    &
        ONLY:  chem_emis_biogenic_init

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_generic_init

    USE chem_emis_domestic_mod,                                                                    &
        ONLY:  chem_emis_domestic_init

    USE chem_emis_nonstationary_mod,                                                               &
        ONLY:  chem_emis_nonstationary_init

    USE chem_emis_pollen_mod,                                                                      &
        ONLY:  chem_emis_pollen_init

    USE chem_emis_pt_source_mod,                                                                   &
        ONLY:  chem_emis_pt_source_init

    USE chem_emis_traffic_mod,                                                                     &
        ONLY:  chem_emis_traffic_init

    USE chem_emissions_mod,                                                                        &
        ONLY:  chem_emissions_header_init,                                                         &
               chem_emissions_init

    USE chem_isorropia_mod,                                                                        &
        ONLY:  chem_isorropia_init

    USE chem_modules,                                                                              &
        ONLY:  chem_isorropia,                                                                     &
               chem_wet_deposition,                                                                &
               emis_biogenic,                                                                      &
               emis_domestic,                                                                      &
               emis_generic,                                                                       &
               emis_nonstationary,                                                                 &
               emis_pollen,                                                                        &
               emis_pt_source,                                                                     &
               emis_traffic

    USE chem_wet_deposition_mod,                                                                   &
        ONLY:  chem_wet_deposition_init

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  init_3d


    INTEGER(iwp) ::  i !< running index x dimension
    INTEGER(iwp) ::  j !< running index y dimension
    INTEGER(iwp) ::  n !< running index for chemical species


    IF ( debug_output )  CALL debug_message( 'chem_init', 'start' )
!
!-- Next statement is to avoid compiler warning about unused variables.
    IF ( ( ilu_arable + ilu_coniferous_forest + ilu_deciduous_forest + ilu_mediterrean_scrub +     &
           ilu_permanent_crops + ilu_savanna + ilu_semi_natural_veg + ilu_tropical_forest +        &
           ilu_urban ) == 0 )  CONTINUE

!
!-- NB Calls specific emisisons initialization subroutines for legacy mode and on-demand mode.
    IF  ( emissions_anthropogenic )  THEN

       IF  ( emiss_read_legacy_mode )  THEN
          CALL chem_emissions_init
       ELSE
          CALL chem_emissions_header_init
       ENDIF

    ENDIF

!
!-- Initiate ISORROPIA.
    IF ( chem_isorropia )  CALL chem_isorropia_init( )
!
!-- Initiate wet deposition.
    IF ( chem_wet_deposition )  CALL chem_wet_deposition_init( )

!
!-- Initiate activated emission modes.
    IF ( emis_biogenic      )  CALL chem_emis_biogenic_init
    IF ( emis_generic       )  CALL chem_emis_generic_init
    IF ( emis_domestic      )  CALL chem_emis_domestic_init
    IF ( emis_nonstationary )  CALL chem_emis_nonstationary_init
    IF ( emis_pollen        )  CALL chem_emis_pollen_init
    IF ( emis_pt_source     )  CALL chem_emis_pt_source_init
    IF ( emis_traffic       )  CALL chem_emis_traffic_init

!
!-- Chemistry variables will be initialized if availabe from dynamic input file. Note, it is
!-- possible to initialize only part of the chemistry variables from dynamic input.
    IF ( INDEX( initializing_actions, 'read_from_file' ) /= 0 )  THEN
       DO  n = 1, nspec
          IF ( init_3d%from_file_chem(n) )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   chem_species(n)%conc(:,j,i) = init_3d%chem_init(:,n)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDIF
!
!-- Initialize also the new time level. This is especially required in restart runs to properly
!-- maintain the lateral boundary conditions.
    DO  n = 1, nspec
       chem_species(n)%conc_p = chem_species(n)%conc
    ENDDO

    IF ( debug_output )  CALL debug_message( 'chem_init', 'end' )

 END SUBROUTINE chem_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine initializating chemistry_model_mod
!> internal workaround for chem_species dependency in chem_check_parameters
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_init_internal

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  chem_emis,                                                                          &
               chem_emis_att,                                                                      &
               input_pids_dynamic,                                                                 &
               init_3d,                                                                            &
               netcdf_data_input_chemistry_data

    USE pegrid

!
!-- Local variables
    INTEGER(iwp) ::  i                 !< running index in x-direction
    INTEGER(iwp) ::  j                 !< running index in y-direction
    INTEGER(iwp) ::  k                 !< running index in z-direction
    INTEGER(iwp) ::  lsp               !< running index for chem spcs

    REAL(wp)     ::  flag              !< flag for masking topography/building grid points
!
!-- NB reads netcdf data only under legacy mode

!    IF ( emissions_anthropogenic )  THEN
!       CALL netcdf_data_input_chemistry_data( chem_emis_att, chem_emis )
!    ENDIF

    IF ( emissions_anthropogenic )  THEN
       IF ( emiss_read_legacy_mode )  THEN
          CALL netcdf_data_input_chemistry_data( chem_emis_att, chem_emis )
       ENDIF
    ENDIF

!
!-- Allocate memory for chemical species
    ALLOCATE( chem_species(nspec) )
    ALLOCATE( spec_conc_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) )
    ALLOCATE( spec_conc_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) )
    ALLOCATE( spec_conc_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nspec) )
    ALLOCATE( phot_frequen(nphot) )
    ALLOCATE( freq_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nphot) )
    ALLOCATE( bc_cs_t_val(nspec) )
!
!-- Initialize arrays
    spec_conc_1(:,:,:,:) = 0.0_wp
    spec_conc_2(:,:,:,:) = 0.0_wp
    spec_conc_3(:,:,:,:) = 0.0_wp
    freq_1(:,:,:,:) = 0.0_wp

!
!-- Allocate array to store locally summed-up resolved-scale vertical fluxes.
    IF ( scalar_advec == 'ws-scheme' )  THEN
       ALLOCATE( sums_ws_l(nzb:nzt+1,0:threads_per_task-1,nspec) )
       sums_ws_l = 0.0_wp
    ENDIF

    DO  lsp = 1, nspec
       chem_species(lsp)%name    = spc_names(lsp)

       chem_species(lsp)%conc   (nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_1 (:,:,:,lsp)
       chem_species(lsp)%conc_p (nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_2 (:,:,:,lsp)
       chem_species(lsp)%tconc_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg)       => spec_conc_3 (:,:,:,lsp)

       ALLOCATE (chem_species(lsp)%cssws_av(nysg:nyng,nxlg:nxrg))
       chem_species(lsp)%cssws_av    = 0.0_wp
!
!--    The following block can be useful when emission module is not applied. &
!--    If emission module is applied the following block will be overwritten.
       ALLOCATE (chem_species(lsp)%flux_s_cs(nzb+1:nzt,0:threads_per_task-1))
       ALLOCATE (chem_species(lsp)%diss_s_cs(nzb+1:nzt,0:threads_per_task-1))
       ALLOCATE (chem_species(lsp)%flux_l_cs(nzb+1:nzt,nys:nyn,0:threads_per_task-1))
       ALLOCATE (chem_species(lsp)%diss_l_cs(nzb+1:nzt,nys:nyn,0:threads_per_task-1))
       chem_species(lsp)%flux_s_cs = 0.0_wp
       chem_species(lsp)%flux_l_cs = 0.0_wp
       chem_species(lsp)%diss_s_cs = 0.0_wp
       chem_species(lsp)%diss_l_cs = 0.0_wp
!
!--   Allocate memory for initial concentration profiles (concentration values come from namelist)
!--   (@todo (FK): Because of this, chem_init is called in palm before check_parameters, since
!--                conc_pr_init is used there.
!--                We have to find another solution since chem_init should eventually be called from
!--                init_3d_model!!)
       ALLOCATE ( chem_species(lsp)%conc_pr_init(0:nz+1) )
       chem_species(lsp)%conc_pr_init(:) = 0.0_wp

    ENDDO

!
!-- For chemistry variables lateral boundary conditions can be set non-cyclic while
!-- the other scalars may have cyclic boundary conditions.
!-- However, large gradients near the boundaries may produce stationary numerical
!-- oscillations near the lateral boundaries when a higher-order scheme is
!-- applied near these boundaries.
!-- To get rid-off this, set-up additional flags that control the order of the scalar advection
!-- scheme near the lateral boundaries for passive scalars with non-cyclic bcs
    IF ( scalar_advec == 'ws-scheme' )  THEN
       ALLOCATE( cs_advc_flags_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!--    In case of non-cyclic boundary conditions set topo_flags bit 31
!--    Bit 31 is used to identify extended degradation zones.
!--    Note, since several also other modules like Salsa or other future one may access this bit but
!--    may have other boundary conditions, the original value of topo_flags bit 31 must not
!--    be modified. Hence, store the boundary conditions directly on cs_advc_flags_s.
!--    cs_advc_flags_s will be later overwritten in ws_init_flags_scalar and bit 31 won't be used to
!--    control the numerical order.
!--    Initialize with flag 31 only.
       cs_advc_flags_s = 0
       cs_advc_flags_s = MERGE( IBSET( cs_advc_flags_s, 31 ), 0, BTEST( topo_flags, 31 ) )

       IF ( bc_dirichlet_cs_n .OR. bc_dirichlet_cs_s )  THEN
          IF ( nys == 0  )  THEN
             DO  i = 1, nbgp
                cs_advc_flags_s(:,nys-i,:) = cs_advc_flags_s(:,nys,:)
             ENDDO
          ENDIF
          IF ( nyn == ny )  THEN
             DO  i = 1, nbgp
                cs_advc_flags_s(:,nyn+i,:) = cs_advc_flags_s(:,nyn,:)
             ENDDO
          ENDIF
       ENDIF
       IF ( bc_dirichlet_cs_l .OR. bc_dirichlet_cs_r )  THEN
          IF ( nxl == 0  )  THEN
             DO  i = 1, nbgp
                cs_advc_flags_s(:,:,nxl-i) = cs_advc_flags_s(:,:,nxl)
             ENDDO
          ENDIF
          IF ( nxr == nx )  THEN
             DO  i = 1, nbgp
                cs_advc_flags_s(:,:,nxr+i) = cs_advc_flags_s(:,:,nxr)
             ENDDO
          ENDIF

       ENDIF
!
!--    To initialize advection flags appropriately, pass the boundary flags.
!--    The extensive_degrad argument indicates that a passive scalar is treated, where the
!--    horizontal advection terms are degraded already 2 grid points before the lateral boundary
!--    to avoid stationary oscillations at large-gradients.
!--    Also, extended degradation zones are applied, where horizontal advection of passive scalars
!--    is discretized by first-order scheme at all grid points that in the vicinity of buildings
!--    (<= 3 grid points), even if no building is within the numerical stencil, first-order
!--    scheme is used.
!--    At the fourth and fifth grid point apart from the building, the order of the horizontal
!--    advection scheme is successively increased.
!--    These extended degradation zones are used to avoid stationary numerical oscillations, which
!--    are responsible for high concentration maxima that may appear under shear-free stable
!--    conditions.
       CALL ws_init_flags_scalar( bc_dirichlet_cs_l  .OR.  bc_radiation_cs_l,                      &
                                  bc_dirichlet_cs_n  .OR.  bc_radiation_cs_n,                      &
                                  bc_dirichlet_cs_r  .OR.  bc_radiation_cs_r,                      &
                                  bc_dirichlet_cs_s  .OR.  bc_radiation_cs_s,                      &
                                  cs_advc_flags_s, extensive_degrad = .TRUE.,                      &
                                  alternative_communicator = communicator_chem )
    ENDIF
!
!-- Initial concentration of profiles is prescribed by parameters cs_profile and cs_heights in the
!-- namelist &chemistry_parameters.
    CALL chem_init_profiles
!
!-- In case there is dynamic input file, create a list of names for chemistry initial input files.
!-- Also, initialize array that indicates whether the respective variable is on file or not.
    IF ( input_pids_dynamic )  THEN
       ALLOCATE( init_3d%var_names_chem(1:nspec) )
       ALLOCATE( init_3d%from_file_chem(1:nspec) )
       init_3d%from_file_chem(:) = .FALSE.

       DO  lsp = 1, nspec
          init_3d%var_names_chem(lsp) = init_3d%init_char // TRIM( chem_species(lsp)%name )
       ENDDO
    ENDIF
!
!-- Initialize model variables
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.                                &
         .NOT. cyclic_fill_initialization )                                                        &
    THEN
!
!--    First model run of a possible job queue.
!--    Initial profiles of the variables must be computed.
       IF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
!
!--       Transfer initial profiles to the arrays of the 3D model
!--       Concentrations within buildings are set to zero.
          DO  lsp = 1, nspec
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = 1, nzt+1
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      chem_species(lsp)%conc(k,j,i) = chem_species(lsp)%conc_pr_init(k) * flag
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       ELSEIF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .OR.                  &
                INDEX( initializing_actions, 'interpolate_from_parent' ) /= 0 )  THEN

          DO  lsp = 1, nspec
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt+1
                      flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                      chem_species(lsp)%conc(k,j,i) = chem_species(lsp)%conc_pr_init(k) * flag
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       ENDIF
!
    ENDIF
!
!-- Initial old and new time levels. Note, this has to be done also in restart runs
    DO  lsp = 1, nvar
       chem_species(lsp)%tconc_m = 0.0_wp
       chem_species(lsp)%conc_p  = chem_species(lsp)%conc
    ENDDO

    DO  lsp = 1, nphot
       phot_frequen(lsp)%name = phot_names(lsp)
       phot_frequen(lsp)%freq(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  =>  freq_1(:,:,:,lsp)
    ENDDO

!    CALL photolysis_init   ! probably also required for restart


    RETURN

 END SUBROUTINE chem_init_internal


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining initial vertical profiles of chemical species (given by namelist parameters
!> chem_profiles and chem_heights)  --> which should work analogically to parameters u_profile,
!> v_profile and uv_heights)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_init_profiles

    USE chem_modules

!
!-- Local variables
    INTEGER ::  k          !< running index in z-direction
    INTEGER ::  lsp        !< running index for number of species in derived data type species_def
    INTEGER ::  lsp_usr    !< running index for number of species (user defined)  in cs_names,
                           !< cs_profiles etc
    INTEGER ::  npr_lev    !< the next available profile lev


!
!-- Parameter "cs_profile" and "cs_heights" are used to prescribe user defined initial profiles
!-- and heights. If parameter "cs_profile" is not prescribed then initial surface values
!-- "cs_surface" are used as constant initial profiles for each species. If "cs_profile" and
!-- "cs_heights" are prescribed, their values will!override the constant profile given by
!-- "cs_surface".
!     IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
    lsp_usr = 1
    DO  WHILE ( TRIM( cs_name( lsp_usr ) ) /= 'novalue' )   !'novalue' is the default
       DO  lsp = 1, nspec                                !
!
!--       Create initial profile (conc_pr_init) for each chemical species
          IF ( TRIM( chem_species(lsp)%name ) == TRIM( cs_name(lsp_usr) ) )  THEN
             IF ( cs_profile(lsp_usr,1) == 9999999.9_wp )  THEN
!
!--            Set a vertically constant profile based on the surface conc (cs_surface(lsp_usr)) of
!--            each species
                DO  k = 0, nzt+1
                   chem_species(lsp)%conc_pr_init(k) = cs_surface(lsp_usr)
                ENDDO
             ELSE
                IF ( cs_heights(1,1) /= 0.0_wp )  THEN
                   WRITE( message_string, * ) 'illegal cs_heights(1,1) = ', cs_heights(1,1)
                   CALL message( 'chem_check_parameters', 'CHM0011', 1, 2, 0, 6, 0 )
                ENDIF

                use_prescribed_profile_data = .TRUE.

                npr_lev = 1
!                chem_species(lsp)%conc_pr_init(0) = 0.0_wp
                DO  k = 1, nzt+1
                   IF ( npr_lev < 100 )  THEN
                      DO  WHILE ( cs_heights(lsp_usr, npr_lev+1) <= zu(k) )
                         npr_lev = npr_lev + 1
                         IF ( npr_lev == 100 )  THEN
                            message_string = 'number of chem spcs exceeding the limit'
                            CALL message( 'chem_check_parameters', 'CHM0012', 1, 2, 0, 6, 0 )
                            EXIT
                         ENDIF
                      ENDDO
                   ENDIF
                   IF ( npr_lev < 100  .AND.  cs_heights(lsp_usr,npr_lev+1) /= 9999999.9_wp )  THEN
                      chem_species(lsp)%conc_pr_init(k) = cs_profile(lsp_usr,npr_lev) +            &
                           ( zu(k) - cs_heights(lsp_usr, npr_lev) ) /                              &
                           ( cs_heights(lsp_usr,npr_lev+1) - cs_heights(lsp_usr,npr_lev) ) *       &
                           ( cs_profile(lsp_usr,npr_lev+1) - cs_profile(lsp_usr,npr_lev) )
                   ELSE
                      chem_species(lsp)%conc_pr_init(k) = cs_profile(lsp_usr, npr_lev)
                   ENDIF
                ENDDO
             ENDIF
!
!--       If a profile is prescribed explicity using cs_profiles and cs_heights, then
!--       chem_species(lsp)%conc_pr_init is populated with the specific "lsp" based on the
!--       cs_profiles(lsp_usr,:)  and cs_heights(lsp_usr,:).
          ENDIF

       ENDDO

       lsp_usr = lsp_usr + 1
    ENDDO
!     ENDIF

 END SUBROUTINE chem_init_profiles


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to integrate chemical species in the given chemical mechanism. Vector-optimized
!> version.
!--------------------------------------------------------------------------------------------------!
  SUBROUTINE chem_integrate

    REAL(wp), PARAMETER ::  fr2ppm  = 1.0E6_wp     !< Conversion factor fraction to ppm
    REAL(wp), PARAMETER ::  p_std   = 101325.0_wp  !< standard pressure (Pa)
    REAL(wp), PARAMETER ::  ppm2fr  = 1.0E-6_wp    !< Conversion factor ppm to fraction
    REAL(wp), PARAMETER ::  t_std   = 273.15_wp    !< standard pressure (Pa)
    REAL(wp), PARAMETER ::  vmolcm  = 22.414E3_wp  !< Mole volume (22.414 l) in cm^3
    REAL(wp), PARAMETER ::  xna     = 6.022E23_wp  !< Avogadro number (molecules/mol)

    INTEGER(iwp) ::  i      !< running index for x-direction
    INTEGER(iwp) ::  j      !< running index for y-direction
    INTEGER(iwp) ::  ks     !< start index of treated index space
    INTEGER(iwp) ::  ke     !< end index of treated index space
    INTEGER(iwp) ::  lph    !< running index for photolysis frequencies
    INTEGER(iwp) ::  lsp    !< running index for chem species
    INTEGER(iwp) ::  m      !< dummy variable used to determine the start and end index
    INTEGER(iwp) ::  nr_1d  !< number of prognostic grid points per subdomain

    INTEGER, DIMENSION(20) ::  istatus  !< return value of gasphase integration

    REAL(wp)      ::  conv     !< conversion factor
    REAL(kind=wp) ::  dt_chem  !< chemistry timestep

    REAL(wp),DIMENSION(size(rcntrl))   :: rcntrl_local

    REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   ::  tmp_fact
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   ::  tmp_fact_i  !< conversion factor between molecules cm^{-3} and ppm
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   ::  tmp_qvap
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:)   ::  tmp_temp

    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:) ::  tmp_conc
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:) ::  tmp_phot


!
!-- Set chem_gasphase_on to .FALSE. if you want to skip computation of gas phase chemistry.
    IF ( chem_gasphase_on )  THEN
!
!--    Compute length of time step
       IF ( call_chem_at_all_substeps )  THEN
          dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_chem = dt_3d
       ENDIF

       nr_1d = (nxr-nxl+1)*(nyn-nys+1)*(nzt-nzb+2)

       ALLOCATE( tmp_fact(nr_1d) )
       ALLOCATE( tmp_fact_i(nr_1d) )
       ALLOCATE( tmp_qvap(nr_1d) )
       ALLOCATE( tmp_temp(nr_1d) )
       ALLOCATE( tmp_conc(nr_1d,nspec) )
       ALLOCATE( tmp_phot(nr_1d,nphot) )

       IF ( MAXVAL( rcntrl ) > 0.0 )  THEN
          IF( time_since_reference_point <= 2.0_wp * dt_3d )  THEN
             rcntrl_local = 0
          ELSE
             rcntrl_local = rcntrl
          ENDIF
       ELSE
          rcntrl_local = 0
       ENDIF

       m = 1

       cs_time_step = dt_chem
!
!--    Pre-calculate arrays for concentration, water vapor, photolysis frequency, etc.,
!--    for the entire subdomain and pass them as a whole to chem_integrate.
!--    In the following loop, the ks and ke indices represent the index space of a vertical
!--    column (at j,i) with respect to a 1D-array for all prognostic grid point of the
!--    subdomain.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             ks = m
             ke = m + ( nzt - nzb + 1 )
             m  = m + ( nzt - nzb + 2 )

             tmp_temp(ks:ke) = pt(nzb+1:nzt,j,i) * exner(nzb+1:nzt)
!
!--          Convert ppm to molecules/cm**3
!--          tmp_fact = 1.e-6_wp*6.022e23_wp/(22.414_wp*1000._wp) * 273.15_wp *
!--                     hyp(nzb+1:nzt)/( 101300.0_wp * tmp_temp )
             conv = ppm2fr * xna / vmolcm
             tmp_fact(ks:ke) = conv * t_std * hyp(nzb+1:nzt) / ( tmp_temp(ks:ke) * p_std )
             tmp_fact_i(ks:ke) = 1.0_wp / tmp_fact(ks:ke)

             IF ( humidity )  THEN
                IF ( bulk_cloud_model )  THEN
                   tmp_qvap(ks:ke) = ( q(nzb+1:nzt,j,i) - ql(nzb+1:nzt,j,i) ) *                    &
                                     xm_air / xm_h2o * fr2ppm * tmp_fact(ks:ke)
                ELSE
                   tmp_qvap(ks:ke) = q(nzb+1:nzt,j,i) * xm_air / xm_h2o * fr2ppm * tmp_fact(ks:ke)
                ENDIF
             ELSE
!
!--             Constant value for q if water vapor is not computed.
                tmp_qvap(ks:ke) = 0.01 * xm_air / xm_h2o * fr2ppm * tmp_fact(ks:ke)
             ENDIF

             DO  lsp = 1, nspec
                tmp_conc(ks:ke,lsp) = chem_species(lsp)%conc(nzb+1:nzt,j,i) * tmp_fact(ks:ke)
             ENDDO

             DO lph = 1, nphot
                tmp_phot(ks:ke,lph) = phot_frequen(lph)%freq(nzb+1:nzt,j,i)
             ENDDO

          ENDDO
       ENDDO

       CALL chem_gasphase_integrate( dt_chem, tmp_conc, tmp_temp, tmp_qvap, tmp_fact, tmp_phot,    &
                                     icntrl_i = icntrl, rcntrl_i = rcntrl_local, istatus=istatus )

       m = 1
       DO  i = nxl, nxr
          DO  j = nys, nyn
             ks = m
             ke = m + ( nzt - nzb + 1 )
             m  = m + ( nzt - nzb + 2 )
             DO  lsp = 1, nspec
                chem_species(lsp)%conc(nzb+1:nzt,j,i) = tmp_conc(ks:ke,lsp) * tmp_fact_i(ks:ke)
             ENDDO
          ENDDO
       ENDDO

    ENDIF

 END SUBROUTINE chem_integrate


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to integrate chemical species in the given chemical mechanism. Cache-optimized
!> version.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_integrate_ij( i, j )

    REAL(wp), PARAMETER ::  fr2ppm  = 1.0e6_wp              !< Conversion factor fraction to ppm
    REAL(wp), PARAMETER ::  p_std   = 101325.0_wp           !< standard pressure (Pa)
    REAL(wp), PARAMETER ::  ppm2fr  = 1.0e-6_wp             !< Conversion factor ppm to fraction
    REAL(wp), PARAMETER ::  t_std   = 273.15_wp             !< standard pressure (Pa)
    REAL(wp), PARAMETER ::  vmolcm  = 22.414e3_wp           !< Mole volume (22.414 l) in cm^3
    REAL(wp), PARAMETER ::  xna     = 6.022e23_wp           !< Avogadro number (molecules/mol)

    INTEGER,INTENT(IN) ::  i
    INTEGER,INTENT(IN) ::  j
!
!-- Local variables
    INTEGER(iwp) ::  lph  !< running index for photolysis frequencies
    INTEGER(iwp) ::  lsp  !< running index for chem spcs.

    INTEGER, DIMENSION(20) ::  istatus

    INTEGER,DIMENSION(nzb+1:nzt) ::  nacc  !< Number of accepted steps
    INTEGER,DIMENSION(nzb+1:nzt) ::  nrej  !< Number of rejected steps

    REAL(wp)      ::  conv     !< conversion factor
    REAL(KIND=wp) ::  dt_chem  !<

    REAL(wp),DIMENSION(size(rcntrl)) :: rcntrl_local  !<

    REAL(KIND=wp), DIMENSION(nzb+1:nzt) :: tmp_fact    !<
    REAL(KIND=wp), DIMENSION(nzb+1:nzt) :: tmp_fact_i  !< conversion factor between  molecules cm^{-3} and ppm
    REAL(KIND=wp), DIMENSION(nzb+1:nzt) :: tmp_qvap    !<
    REAL(KIND=wp), DIMENSION(nzb+1:nzt) :: tmp_temp    !<

    REAL(KIND=wp), DIMENSION(nzb+1:nzt,nspec) :: tmp_conc  !<
    REAL(KIND=wp), DIMENSION(nzb+1:nzt,nphot) :: tmp_phot  !<


!
!-- Set chem_gasphase_on to .FALSE. if you want to skip computation of gas phase chemistry.
    IF ( chem_gasphase_on )  THEN
       nacc = 0
       nrej = 0

       tmp_temp(:) = pt(nzb+1:nzt,j,i) * exner(nzb+1:nzt)
!
!--    Convert ppm to molecules/cm**3
!--    tmp_fact = 1.e-6_wp*6.022e23_wp/(22.414_wp*1000._wp) * 273.15_wp *
!--               hyp(nzb+1:nzt)/( 101300.0_wp * tmp_temp )
       conv = ppm2fr * xna / vmolcm
       tmp_fact(:) = conv * t_std * hyp(nzb+1:nzt) / ( tmp_temp(:) * p_std )
       tmp_fact_i = 1.0_wp / tmp_fact

       IF ( humidity )  THEN
          IF ( bulk_cloud_model )  THEN
             tmp_qvap(:) = ( q(nzb+1:nzt,j,i) - ql(nzb+1:nzt,j,i) ) *                              &
                             xm_air / xm_h2o * fr2ppm * tmp_fact(:)
          ELSE
             tmp_qvap(:) = q(nzb+1:nzt,j,i) * xm_air / xm_h2o * fr2ppm * tmp_fact(:)
          ENDIF
       ELSE
!
!--       Constant value for q if water vapor is not computed.
          tmp_qvap(:) = 0.01_wp * xm_air / xm_h2o * fr2ppm * tmp_fact(:)
       ENDIF

       DO  lsp = 1, nspec
          tmp_conc(:,lsp) = chem_species(lsp)%conc(nzb+1:nzt,j,i) * tmp_fact(:)
       ENDDO

       DO lph = 1, nphot
          tmp_phot(:,lph) = phot_frequen(lph)%freq(nzb+1:nzt,j,i)
       ENDDO
!
!--    Compute length of time step.
       IF ( call_chem_at_all_substeps )  THEN
          dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_chem = dt_3d
       ENDIF

       cs_time_step = dt_chem

       IF ( MAXVAL( rcntrl ) > 0.0 )  THEN    ! Only if rcntrl is set
          IF( time_since_reference_point <= 2.0_wp * dt_3d)  THEN
             rcntrl_local = 0
          ELSE
             rcntrl_local = rcntrl
          ENDIF
       ELSE
          rcntrl_local = 0
       ENDIF

       CALL chem_gasphase_integrate ( dt_chem, tmp_conc, tmp_temp, tmp_qvap, tmp_fact, tmp_phot,   &
                                      icntrl_i = icntrl, rcntrl_i = rcntrl_local, xnacc = nacc,    &
                                      xnrej = nrej, istatus=istatus )

       DO  lsp = 1,nspec
          chem_species(lsp)%conc(nzb+1:nzt,j,i) = tmp_conc(:,lsp) * tmp_fact_i(:)
       ENDDO

    ENDIF

 END SUBROUTINE chem_integrate_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining parin for &chemistry_parameters for chemistry model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_parin

    USE chem_modules
    USE control_parameters
    USE pegrid
    USE statistics

    CHARACTER(LEN=8)   ::  solver_type  !<
    CHARACTER(LEN=100) ::  line         !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status  !< Status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    REAL(wp), DIMENSION(nmaxfixsteps) ::  my_steps  !< List of fixed timesteps  my_step(1) = 0.0 automatic stepping


    NAMELIST /chemistry_parameters/                                                                &
         bc_cs_b,                                                                                  &
         bc_cs_l,                                                                                  &
         bc_cs_n,                                                                                  &
         bc_cs_r,                                                                                  &
         bc_cs_s,                                                                                  &
         bc_cs_t,                                                                                  &
         call_chem_at_all_substeps,                                                                &
         chem_gasphase_on,                                                                         &
         chem_isorropia,                                                                           &
         chem_isorropia_activity_coefficient_method,                                               &
         chem_isorropia_activity_tolerance,                                                        &
         chem_isorropia_aerosol_state,                                                             &
         chem_isorropia_mass_conservation_mode,                                                    &
         chem_isorropia_max_iteration,                                                             &
         chem_isorropia_max_activity_sweep,                                                        &
         chem_isorropia_mdr_weight_method,                                                         &
         chem_isorropia_root_subdivisions,                                                         &
         chem_isorropia_problem_type,                                                              &
         chem_isorropia_solver_tolerance,                                                          &
         chem_isorropia_update_interval,                                                           &
         chem_mechanism,                                                                           &
         chem_wet_deposition,                                                                      &
         chem_wet_deposition_model_override,                                                       &
         chem_wet_deposition_cloud_level_lower,                                                    &
         chem_wet_deposition_cloud_level_upper,                                                    &
         chem_wet_deposition_rain_rate,                                                            &
         chem_wet_deposition_update_interval,                                                      &
         cs_heights,                                                                               &
         cs_name,                                                                                  &
         cs_profile,                                                                               &
         cs_surface,                                                                               &
         deposition_dry,                                                                           &
         daytype_mdh,                                                                              &
         ebio_dt,                                                                                  &
         ebio_ef_pft,                                                                              &
         ebio_ef_tree,                                                                             &
         ebio_emis_name,                                                                           &
         ebio_max_emis_day,                                                                        &
         ebio_pft,                                                                                 &
         ebio_ppfd_factor,                                                                         &
         ebio_rad_method,                                                                          &
         ebio_soilm_method,                                                                        &
         ebio_tree,                                                                                &
         emissions_anthropogenic,                                                                  &
         emiss_factor_main,                                                                        &
         emiss_factor_side,                                                                        &
         emiss_interpolate,                                                                        &
         emiss_lod,                                                                                &
         emiss_read_legacy_mode,                                                                   &
         emis_biogenic,                                                                            &
         emis_biogenic_lod,                                                                        &
         emis_domestic,                                                                            &
         emis_domestic_base_temperature,                                                           &
         emis_domestic_compact_factors,                                                            &
         emis_domestic_energy_demands,                                                             &
         emis_domestic_heating_degree,                                                             &
         emis_domestic_lod,                                                                        &
         emis_domestic_sampling_k,                                                                 &
         emis_domestic_species_emission_factors,                                                   &
         emis_domestic_species_names,                                                              &
         emis_domestic_update_interval,                                                            &
         emis_generic,                                                                             &
         emis_nonstationary,                                                                       &
         emis_pollen,                                                                              &
         epol_ignore_precip,                                                                       &
         epol_ignore_solar,                                                                        &
         epol_model,                                                                               &
         epol_pool_reset_hour,                                                                     &
         epol_seasonal_factors,                                                                    &
         epol_specs_names,                                                                         &
         epol_tke_scheme,                                                                          &
         epol_tke_sgs_fraction,                                                                    &
         epol_tree_specs,                                                                          &
         epol_tuning_factors,                                                                      &
         epol_update_interval,                                                                     &
         epol_vegetation_specs,                                                                    &
         emis_pt_source,                                                                           &
         emis_pt_source_annual_values,                                                             &
         emis_pt_source_k_spread,                                                                  &
         emis_pt_source_k_weights,                                                                 &
         emis_pt_source_leap_year,                                                                 &
         emis_pt_source_locations_ijk,                                                             &
         emis_pt_source_species_names,                                                             &
         emis_traffic,                                                                             &
         emis_traffic_lod,                                                                         &
         icntrl,                                                                                   &
         main_street_id,                                                                           &
         max_street_id,                                                                            &
         mode_emis,                                                                                &
         my_steps,                                                                                 &
         nesting_chem,                                                                             &
         nesting_offline_chem,                                                                     &
         photolysis_scheme,                                                                        &
         photolysis_shading,                                                                       &
         rcntrl,                                                                                   &
         side_street_id,                                                                           &
         surface_csflux,                                                                           &
         surface_csflux_name,                                                                      &
         switch_off_module,                                                                        &
         time_fac_type,                                                                            &
         wall_csflux
!
!-- Analogically to chem_names(nspj) we could invent chem_surfaceflux(nspj) and chem_topflux(nspj)
!-- so this way we could prescribe a specific flux value for each species
    !>  chemistry_parameters for initial profiles
    !>  cs_names = 'O3', 'NO2', 'NO', ...   to set initial profiles)
    !>  cs_heights(1,:) = 0.0, 100.0, 500.0, 2000.0, .... (height levels where concs will be prescribed for O3)
    !>  cs_heights(2,:) = 0.0, 200.0, 400.0, 1000.0, .... (same for NO2 etc.)
    !>  cs_profiles(1,:) = 10.0, 20.0, 20.0, 30.0, .....  (chem spcs conc at height lvls chem_heights(1,:)) etc.
    !>  If the respective concentration profile should be constant with height, then use "cs_surface( number of spcs)"
    !>  then write these cs_surface values to chem_species(lsp)%conc_pr_init(:)

!
!-- Read chem namelist.
    icntrl    = 0
    rcntrl    = 0.0_wp
    my_steps  = 0.0_wp
    photolysis_scheme = 'simple'
    atol = 1.0_wp
    rtol = 0.01_wp
!
!-- Move to the beginning of the namelist file and try to find and read the namelist named.
!-- chemistry_parameters.
    REWIND( 11 )
    READ( 11, chemistry_parameters, IOSTAT=io_status )
!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!      chemistry_parameters namelist was found and read correctly. Switch on chemistry model.
       IF ( .NOT. switch_off_module )  air_chemistry = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    chemistry_parameters namelist was found, but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'chemistry_parameters', line )

    ENDIF

!
!-- Synchronize emiss_lod and mod_emis only if emissions_anthropogenic is activated in the namelist.
!-- Otherwise their values are "don't care"
    IF ( emissions_anthropogenic )  THEN

!
!--    Check for emission mode for chem species
       IF ( emiss_lod < 0 )  THEN   !- if LOD not defined in namelist
          IF ( ( mode_emis /= 'PARAMETERIZED'  )    .AND.                                          &
               ( mode_emis /= 'DEFAULT'        )    .AND.                                          &
               ( mode_emis /= 'PRE-PROCESSED'  ) )  THEN
             message_string = 'unknown mode_emiss = "' // TRIM( mode_emis ) // '"'
             CALL message( 'chem_parin', 'CHM0013', 1, 2, 0, 6, 0 )
          ENDIF
       ELSE
          IF ( ( emiss_lod /= 0 )    .AND.                                                         &
               ( emiss_lod /= 1 )    .AND.                                                         &
               ( emiss_lod /= 2 ) )  THEN
             WRITE( message_string, * ) 'illegal emiss_lod = ', emiss_lod
             CALL message( 'chem_parin', 'CHM0014', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

!
!--    Conflict resolution for emiss_lod and mode_emis.
!--    1) if emiss_lod is defined, have mode_emis assume same setting as emiss_lod
!--    2) if emiss_lod it not defined, have emiss_lod assuem same setting as mode_emis
!--    This check is in place to retain backward compatibility with salsa until the code is
!--    migrated completely to emiss_lod.
       IF  ( emiss_lod >= 0 ) THEN

          SELECT CASE  ( emiss_lod )
!
!--          Synchronize mode_emis to defined emiss_lod (mode_emis will be depreciated in
!--          future releases)
             CASE (0)  !- parameterized mode
                mode_emis = 'PARAMETERIZED'
             CASE (1)  !- default mode
                mode_emis = 'DEFAULT'
             CASE (2)  !- preprocessed mode
                mode_emis = 'PRE-PROCESSED'
          END SELECT

       ELSE ! if emiss_lod is not set

          SELECT CASE ( mode_emis )
             CASE ('PARAMETERIZED')
                emiss_lod = 0
             CASE ('DEFAULT')
                emiss_lod = 1
             CASE ('PRE-PROCESSED')
                emiss_lod = 2
          END SELECT

          message_string = 'emiss_lod undefined.  Using existing mode_emis setting&'    //         &
                           'NOTE - mode_emis will be depreciated in future releases.&'   //        &
                           'Please use emiss_lod to define emission mode.'
          CALL message( 'chem_parin', 'CHM0015', 0, 1, 0, 6, 0 )
       ENDIF

!
!-- NB input check for emission read mode.
!--    legacy : business as usual (everything read / set up at start of run)
!--    new    : emission based on timestamp, and for lod2 data is loaded on an hourly basis

!
!-- NB handler for emiss_read_legacy_mode
!-- * emiss_read_legacy_mode is defaulted to TRUE
!-- * if emiss_read_legacy_mode is TRUE and LOD is 0 or 1,
!--       force emission_read_legacy_mode to TRUE (not yet implemented)
       IF ( .NOT. emiss_read_legacy_mode )  THEN    !< if new read mode selected

          IF ( emiss_lod < 2 )  THEN            !< check LOD compatibility

             message_string = 'New emission read mode currently unavailable for LODs 0 and 1.&' // &
                              'Reverting to legacy emission read mode.'
             CALL message( 'chem_parin', 'CHM0016', 0, 0, 0, 6, 0 )

             emiss_read_legacy_mode = .TRUE.

          ELSE                                  !< notify new read mode

             message_string = 'New emission read mode activated.& LOD 2 emissions will be ' //     &
                              'updated on-demand according to indicated timestamps.'
             CALL message( 'chem_parin', 'CHM0017', 0, 0, 0, 6, 0 )

          ENDIF

       ENDIF ! if emiss_read_legacy_mode

    ENDIF  ! if emissions_anthropengic

    t_steps = my_steps

!
!-- Set Solver Type
    IF ( icntrl(3) == 0 )  THEN
       solver_type = 'rodas3'           !Default
    ELSEIF ( icntrl(3) == 1 )  THEN
       solver_type = 'ros2'
    ELSEIF ( icntrl(3) == 2 )  THEN
       solver_type = 'ros3'
    ELSEIF ( icntrl(3) == 3 )  THEN
       solver_type = 'ro4'
    ELSEIF ( icntrl(3) == 4 )  THEN
       solver_type = 'rodas3'
    ELSEIF ( icntrl(3) == 5 )  THEN
       solver_type = 'rodas4'
    ELSEIF ( icntrl(3) == 6 )  THEN
       solver_type = 'Rang3'
    ELSE
       WRITE( message_string, * ) 'illegal Rosenbrock-solver type icntrl(3) = ', icntrl(3)
       CALL message( 'chem_parin', 'CHM0018', 1, 2, 0, 6, 0 )
    END IF

    RETURN

 END SUBROUTINE chem_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_actions( location )

    USE chem_emis_biogenic_mod,                                                                    &
        ONLY:  chem_emis_biogenic_update

    USE chem_emis_domestic_mod,                                                                    &
        ONLY:  chem_emis_domestic_update,                                                          &
               chem_emis_domestic_cleanup

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_generic_update,                                                           &
               chem_emis_generic_cleanup

    USE chem_emis_nonstationary_mod,                                                               &
        ONLY:  chem_emis_nonstationary_update

    USE chem_emis_pollen_mod,                                                                      &
        ONLY:  chem_emis_pollen_update

    USE chem_emis_pt_source_mod,                                                                   &
        ONLY:  chem_emis_pt_source_update,                                                         &
               chem_emis_pt_source_cleanup

    USE chem_emis_traffic_mod,                                                                     &
        ONLY:  chem_emis_traffic_update,                                                           &
               chem_emis_traffic_cleanup

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_reset_source,                                                        &
               chem_emis_vsrc_cleanup

    USE chem_emissions_mod,                                                                        &
        ONLY:  chem_emissions_setup,                                                               &
               chem_emissions_update_on_demand

    USE chem_isorropia_mod,                                                                        &
        ONLY:  chem_isorropia_update,                                                              &
               chem_isorropia_cleanup

    USE chem_wet_deposition_mod,                                                                   &
        ONLY:  chem_wet_deposition_cleanup

    USE chem_modules,                                                                              &
        ONLY:  chem_isorropia,                                                                     &
               chem_wet_deposition,                                                                &
               emis_biogenic,                                                                      &
               emis_generic,                                                                       &
               emis_domestic,                                                                      &
               emis_nonstationary,                                                                 &
               emis_pollen,                                                                        &
               emis_pt_source,                                                                     &
               emis_traffic

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  chem_emis,                                                                          &
               chem_emis_att

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

    INTEGER(iwp) ::  hour                !< hour of current time
    INTEGER(iwp) ::  hour_call_emis = -1 !< last hour where emission was called

    SELECT CASE ( location )

       CASE ( 'before_prognostic_equations' )
!
!--       Chemical reactions and deposition.
          IF ( chem_gasphase_on )  THEN
!
!--          If required, calculate photolysis frequencies -
!--          UNFINISHED: Why not before the intermediate timestep loop?
             IF ( intermediate_timestep_count ==  1 )  THEN
                CALL photolysis_control
             ENDIF

          ENDIF

       CASE ( 'before_timestep' )
!
!--       Set array used to sum-up resolved scale fluxes to zero.
          IF ( ws_scheme_sca )  THEN
             sums_ws_l = 0.0_wp
          ENDIF

       CASE ( 'update_emission_sources' )
!
!--       If required, consider chemical emissions.
!--       Allows for emission update mode in legacy mode as well as on-demand mode.
          IF ( emissions_anthropogenic )  THEN
             IF ( emiss_read_legacy_mode )  THEN
!
!--             Get hourly index and updates emission data when the hour is passed.
                CALL get_date_time( time_since_reference_point, hour = hour )

                IF ( hour_call_emis /= hour )  THEN
                   CALL chem_emissions_setup( chem_emis_att, chem_emis, n_matched_vars )
                   hour_call_emis = hour
                ENDIF
             ELSE
                CALL chem_emissions_update_on_demand
             ENDIF
          ENDIF
!
!--       Updates emission sources for all activated modes.
!--       Attention: The way it is done here uses an Euler timestep after all
!--       the prognostics have been done and the ghost-point exchange has been carried out.
!--       Also, it need to be checked if all boundary-conditions are set correctly!
          CALL chem_emis_vsrc_reset_source( )
          IF ( emis_biogenic  )  CALL chem_emis_biogenic_update( )
          IF ( emis_generic   )  CALL chem_emis_generic_update( )
          IF ( emis_domestic  )  CALL chem_emis_domestic_update( )
          IF ( emis_pollen    )  CALL chem_emis_pollen_update( )
          IF ( emis_pt_source )  CALL chem_emis_pt_source_update( )
          IF ( emis_traffic   )  CALL chem_emis_traffic_update( )
!
!--       Impose nonstationary emissions sources.
          IF ( emis_nonstationary )  CALL chem_emis_nonstationary_update
!
!--       SIA update (w/ ISORROPIA).
          IF ( chem_isorropia )  CALL chem_isorropia_update( )

       CASE ( 'after_time_integration' )
!
!--       Deallocates all dynamic memory.
          IF ( chem_isorropia      )  CALL chem_isorropia_cleanup( )
          IF ( chem_wet_deposition )  CALL chem_wet_deposition_cleanup( )
          IF ( emis_domestic       )  CALL chem_emis_domestic_cleanup( )
          IF ( emis_generic        )  CALL chem_emis_generic_cleanup( )
          IF ( emis_pt_source      )  CALL chem_emis_pt_source_cleanup( )
          IF ( emis_traffic        )  CALL chem_emis_traffic_cleanup( )
          CALL chem_emis_vsrc_cleanup ( )

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE chem_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid points i,j
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE chem_actions_ij( i, j, location )

    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string

    INTEGER(iwp)  ::  dummy                     !< call location string

    INTEGER(iwp),      INTENT(IN) ::  i         !< grid index in x-direction
    INTEGER(iwp),      INTENT(IN) ::  j         !< grid index in y-direction

    IF ( air_chemistry    )   dummy = i + j

    SELECT CASE ( location )

       CASE DEFAULT
          CONTINUE

    END SELECT


 END SUBROUTINE chem_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_non_advective_processes()

    USE chem_modules,                                                                              &
        ONLY:  emis_pollen

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_assign_source

    USE chem_emis_pollen_mod,                                                                      &
        ONLY:  chem_emis_pollen_update_settling

    USE chem_wet_deposition_mod,                                                                   &
        ONLY:  chem_wet_deposition_update_ij

    INTEGER(iwp) ::  i  !< grid index in x-direction
    INTEGER(iwp) ::  j  !< grid index in y-direction


!
!-- Calculation of chemical reactions and deposition.
    IF ( intermediate_timestep_count == 1  .OR.  call_chem_at_all_substeps )  THEN
!
!--    For scalars or short vectors use one z column as base for spitting into vector chunks.
!--    vl_dim is created by kpp.
       IF ( vl_dim <= 32 )  THEN
          !$OMP PARALLEL PRIVATE (i,j)
          !$OMP DO schedule(static,1)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                CALL chem_emis_vsrc_assign_source( i, j )
             ENDDO
          ENDDO
          !$OMP END PARALLEL
!
!--    For long vectors use the whole local 3D array column as base for spitting into vector
!--    chunks. This method uses more memory as above
       ELSE
          !$OMP PARALLEL PRIVATE (i,j)
          !$OMP DO schedule(static,1)
          DO  i = nxl, nxr
             DO  j = nys, nyn
                CALL chem_emis_vsrc_assign_source( i, j )
             ENDDO
          ENDDO
          !$OMP END PARALLEL
       ENDIF

       IF ( chem_gasphase_on )  THEN
          CALL cpu_log( log_point_s(19), 'chem.reactions', 'start' )
!
!--       For scalars or short vectors use one z column as base for spitting into vector chunks.
!--       vl_dim is created by kpp.
          IF ( vl_dim <= 32 )  THEN
             !$OMP PARALLEL PRIVATE (i,j)
             !$OMP DO schedule(static,1)
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   CALL chem_integrate( i, j )
                ENDDO
             ENDDO
             !$OMP END PARALLEL
!
!--       For long vectors use the whole local 3D array column as base for spitting into vector
!--       chunks. This method uses more memory as above
          ELSE

             CALL chem_integrate

          ENDIF

          CALL cpu_log( log_point_s(19), 'chem.reactions', 'stop' )

       ENDIF

       IF ( deposition_dry )  THEN
          CALL cpu_log( log_point_s(24), 'chem.deposition', 'start' )
          CALL chem_depo
          CALL cpu_log( log_point_s(24), 'chem.deposition', 'stop' )
       ENDIF

       IF ( emis_pollen )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                CALL chem_emis_pollen_update_settling( i, j )
             ENDDO
          ENDDO
       ENDIF

       IF ( chem_wet_deposition )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                CALL chem_wet_deposition_update_ij( i, j )
             ENDDO
          ENDDO
       ENDIF

    ENDIF

 END SUBROUTINE chem_non_advective_processes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid points i,j.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_non_advective_processes_ij( i, j )

    USE chem_modules,                                                                              &
        ONLY:  emis_pollen

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_assign_source

    USE chem_emis_pollen_mod,                                                                      &
        ONLY:  chem_emis_pollen_update_settling

    USE chem_wet_deposition_mod,                                                                   &
        ONLY:  chem_wet_deposition_update_ij

    INTEGER(iwp), INTENT(IN) ::  i  !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  j  !< grid index in y-direction


!
!-- Calculation of chemical reactions and deposition.
!-- It would have been nice to have time measurements for chemistry and deposition here.
!-- Unfortunately measurements within i,j loops degrade performance ince they are calles so often
!-- and the counter for this measurement gets extremely huge values. Therefore, no measurements
!-- here.
    IF ( intermediate_timestep_count == 1  .OR.  call_chem_at_all_substeps )  THEN

       IF ( chem_gasphase_on )  THEN
          CALL chem_emis_vsrc_assign_source ( i, j )
          CALL chem_integrate( i, j )
       ENDIF

       IF ( deposition_dry )  THEN
          CALL chem_depo( i, j )
       ENDIF

       IF ( emis_pollen )  THEN
          CALL chem_emis_pollen_update_settling( i , j )
       ENDIF

       IF ( chem_wet_deposition )  THEN
           CALL chem_wet_deposition_update_ij( i, j )
       ENDIF

    ENDIF

 END SUBROUTINE chem_non_advective_processes_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for exchange horiz of chemical quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_exchange_horiz_bounds( location )

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string

    INTEGER(iwp) ::  lsp  !<
    INTEGER(iwp) ::  n    !<


    SELECT CASE ( location )

       CASE ( 'before_prognostic_equation' )
!
!--       Loop over chemical species.
          CALL cpu_log( log_point_s(84), 'chem.exch-horiz', 'start' )
          DO  lsp = 1, nvar
             CALL exchange_horiz( chem_species(lsp)%conc, nbgp,                                    &
                                  alternative_communicator = communicator_chem )
          ENDDO

          CALL chem_boundary_conditions( horizontal_conditions_only = .TRUE. )

          CALL cpu_log( log_point_s(84), 'chem.exch-horiz', 'stop' )

       CASE ( 'after_prognostic_equation' )

          IF ( air_chemistry )  THEN
             DO  n = 1, nvar
                CALL exchange_horiz( chem_species(n)%conc_p, nbgp,                                 &
                                     alternative_communicator = communicator_chem )
             ENDDO
          ENDIF

       CASE ( 'after_anterpolation' )

          IF ( air_chemistry )  THEN
             DO  n = 1, nvar
                CALL exchange_horiz( chem_species(n)%conc, nbgp,                                   &
                                     alternative_communicator = communicator_chem )
             ENDDO
          ENDIF

    END SELECT

 END SUBROUTINE chem_exchange_horiz_bounds


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine calculating prognostic equations for chemical species (vector-optimized).
!> Routine is called separately for each chemical species over a loop from prognostic_equations.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_prognostic_equations()

    INTEGER ::  i   !< running index
    INTEGER ::  j   !< running index
    INTEGER ::  k   !< running index

    INTEGER(iwp) ::  ilsp   !<


    CALL cpu_log( log_point_s(25), 'chem.advec+diff+prog', 'start' )

    DO  ilsp = 1, nvar
!
!--    Tendency terms for chemical species
       tend = 0.0_wp
!
!--    Advection terms
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( ws_scheme_sca )  THEN
             sums_wschs_ws_l(nzb:,0:) => sums_ws_l(:,:,ilsp)

             IF ( .NOT. advanced_div_correction )  THEN
                CALL advec_s_ws( cs_advc_flags_s, chem_species(ilsp)%conc, 'kc',                   &
                                 bc_dirichlet_cs_l  .OR.  bc_radiation_cs_l,                       &
                                 bc_dirichlet_cs_n  .OR.  bc_radiation_cs_n,                       &
                                 bc_dirichlet_cs_r  .OR.  bc_radiation_cs_r,                       &
                                 bc_dirichlet_cs_s  .OR.  bc_radiation_cs_s )
             ELSE
                CALL advec_s_ws( cs_advc_flags_s, chem_species(ilsp)%conc, 'kc',                   &
                                 bc_dirichlet_cs_l  .OR.  bc_radiation_cs_l,                       &
                                 bc_dirichlet_cs_n  .OR.  bc_radiation_cs_n,                       &
                                 bc_dirichlet_cs_r  .OR.  bc_radiation_cs_r,                       &
                                 bc_dirichlet_cs_s  .OR.  bc_radiation_cs_s,                       &
                                 advanced_div_correction )
             ENDIF
          ELSE
             CALL advec_s_pw( chem_species(ilsp)%conc )
          ENDIF
       ELSE
          CALL advec_s_up( chem_species(ilsp)%conc )
       ENDIF
!
!--    Diffusion terms  (the last three arguments are zero)
       CALL diffusion_s( chem_species(ilsp)%conc, surf_top%cssws(ilsp,:), surf_def%cssws(ilsp,:),  &
                         surf_lsm%cssws(ilsp,:), surf_usm%cssws(ilsp,:) )
!
!--    Prognostic equation for chemical species
       DO  i = nxl, nxr
          DO  j = nys, nyn
             !following directive is required to vectorize on Intel19
             !DIR$ IVDEP
             DO  k = nzb+1, nzt
                chem_species(ilsp)%conc_p(k,j,i) =   chem_species(ilsp)%conc(k,j,i)                &
                     + ( dt_3d  *                                                                  &
                     (   tsc(2) * tend(k,j,i)                                                      &
                     + tsc(3) * chem_species(ilsp)%tconc_m(k,j,i)                                  &
                     )                                                                             &
                     - tsc(5) * rdf_sc(k)                                                          &
                     * ( chem_species(ilsp)%conc(k,j,i) - chem_species(ilsp)%conc_pr_init(k) )     &
                     )                                                                             &
                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                IF ( chem_species(ilsp)%conc_p(k,j,i) < 0.0_wp  .AND.                              &
                     .NOT. allow_negative_scalar_values )                                          &
                THEN
                   chem_species(ilsp)%conc_p(k,j,i) = 0.1_wp * chem_species(ilsp)%conc(k,j,i)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      chem_species(ilsp)%tconc_m(k,j,i) = tend(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
               intermediate_timestep_count_max )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
                      chem_species(ilsp)%tconc_m(k,j,i) = - 9.5625_wp * tend(k,j,i)                &
                                                     + 5.3125_wp * chem_species(ilsp)%tconc_m(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF

    ENDDO

    CALL cpu_log( log_point_s(25), 'chem.advec+diff+prog', 'stop' )

 END SUBROUTINE chem_prognostic_equations


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine calculating prognostic equations for chemical species (cache-optimized).
!> Routine is called separately for each chemical species over a loop from prognostic_equations.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_prognostic_equations_ij( i, j, i_omp_start, tn )

    INTEGER(iwp),INTENT(IN) :: i, j, i_omp_start, tn

    INTEGER(iwp) :: ilsp
!
!-- local variables

    INTEGER :: k

    DO  ilsp = 1, nvar
!
!--    Tendency-terms for chem spcs.
       tend(:,j,i) = 0.0_wp
!
!--    Advection terms
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( ws_scheme_sca )  THEN
             sums_wschs_ws_l(nzb:,0:) => sums_ws_l(nzb:nzt+1,0:threads_per_task-1,ilsp)
             CALL advec_s_ws( cs_advc_flags_s,                                                     &
                              i,                                                                   &
                              j,                                                                   &
                              chem_species(ilsp)%conc,                                             &
                              'kc',                                                                &
                              chem_species(ilsp)%flux_s_cs,                                        &
                              chem_species(ilsp)%diss_s_cs,                                        &
                              chem_species(ilsp)%flux_l_cs,                                        &
                              chem_species(ilsp)%diss_l_cs,                                        &
                              i_omp_start,                                                         &
                              tn,                                                                  &
                              bc_dirichlet_cs_l  .OR.  bc_radiation_cs_l,                          &
                              bc_dirichlet_cs_n  .OR.  bc_radiation_cs_n,                          &
                              bc_dirichlet_cs_r  .OR.  bc_radiation_cs_r,                          &
                              bc_dirichlet_cs_s  .OR.  bc_radiation_cs_s,                          &
                              monotonic_limiter_z )
          ELSE
             CALL advec_s_pw( i, j, chem_species(ilsp)%conc )
          ENDIF
       ELSE
          CALL advec_s_up( i, j, chem_species(ilsp)%conc )
       ENDIF
!
!--    Diffusion terms (the last three arguments are zero)
       CALL diffusion_s( i, j, chem_species(ilsp)%conc, surf_top%cssws(ilsp,:),                    &
                         surf_def%cssws(ilsp,:), surf_lsm%cssws(ilsp,:), surf_usm%cssws(ilsp,:) )
!
!--    Prognostic equation for chem spcs
       DO  k = nzb+1, nzt
          chem_species(ilsp)%conc_p(k,j,i) = chem_species(ilsp)%conc(k,j,i) + ( dt_3d  *           &
               ( tsc(2) * tend(k,j,i) +                                                            &
               tsc(3) * chem_species(ilsp)%tconc_m(k,j,i) )                                        &
               - tsc(5) * rdf_sc(k)                                                                &
               * ( chem_species(ilsp)%conc(k,j,i) - chem_species(ilsp)%conc_pr_init(k) )           &
               )                                                                                   &
               * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

          IF ( chem_species(ilsp)%conc_p(k,j,i) < 0.0_wp  .AND.                                    &
               .NOT. allow_negative_scalar_values )                                                &
          THEN
             chem_species(ilsp)%conc_p(k,j,i) = 0.1_wp * chem_species(ilsp)%conc(k,j,i)    !FKS6
          ENDIF
       ENDDO
!
!--    Calculate tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             DO  k = nzb+1, nzt
                chem_species(ilsp)%tconc_m(k,j,i) = tend(k,j,i)
             ENDDO
          ELSEIF ( intermediate_timestep_count < &
               intermediate_timestep_count_max )  THEN
             DO  k = nzb+1, nzt
                chem_species(ilsp)%tconc_m(k,j,i) = -9.5625_wp * tend(k,j,i) +                     &
                                                     5.3125_wp * chem_species(ilsp)%tconc_m(k,j,i)
             ENDDO
          ENDIF
       ENDIF

    ENDDO

 END SUBROUTINE chem_prognostic_equations_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync,   &
                                nyn_on_file, nysf, nysc, nys_on_file, tmp_3d, found )

    USE control_parameters


    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  lsp             !<
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

    LOGICAL, INTENT(OUT) :: found

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) &
                 :: tmp_3d   !< 3D array to temp store data


    found = .FALSE.


    IF ( ALLOCATED( chem_species ) )  THEN

       DO  lsp = 1, nspec

          IF ( restart_string(1:length) == TRIM( chem_species(lsp)%name) )  THEN

             IF ( k == 1 )  READ ( 13 )  tmp_3d
             chem_species(lsp)%conc(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                   &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          ELSEIF (restart_string(1:length) == TRIM( chem_species(lsp)%name ) // '_av' )  THEN

             IF ( .NOT. ALLOCATED( chem_species(lsp)%conc_av ) )  THEN
                ALLOCATE( chem_species(lsp)%conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ) )
             ENDIF
             IF ( k == 1 )  READ ( 13 )  tmp_3d
             chem_species(lsp)%conc_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                &
                                                   tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
             found = .TRUE.

          ENDIF

       ENDDO

    ENDIF

 END SUBROUTINE chem_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_rrd_local_mpi

    IMPLICIT NONE

    INTEGER(iwp) ::  lsp  !<

    LOGICAL      ::  array_found  !<

!
!-- Restart input of time-averaged quantities is skipped in case of cyclic-fill initialization.
!-- This case, input of time-averaged data is useless and can lead to faulty averaging.
    IF ( .NOT. cyclic_fill_initialization )  THEN

       DO  lsp = 1, nspec

          CALL rrd_mpi_io( TRIM( chem_species(lsp)%name ), chem_species(lsp)%conc )

          CALL rd_mpi_io_check_array( TRIM( chem_species(lsp)%name )//'_av' , found = array_found )
          IF ( array_found )  THEN
             IF ( .NOT. ALLOCATED( chem_species(lsp)%conc_av ) )  THEN
                ALLOCATE( chem_species(lsp)%conc_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             ENDIF
             CALL rrd_mpi_io( TRIM( chem_species(lsp)%name )//'_av', chem_species(lsp)%conc_av )
          ENDIF

       ENDDO

    ENDIF

 END SUBROUTINE chem_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
!> Description:
!> Calculation of horizontally averaged profiles
!> This routine is called for every statistic region (sr) defined by the user,
!> but at least for the region "total domain" (sr=0).
!> quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_statistics( mode, sr, tn )


    USE arrays_3d

    USE statistics

    CHARACTER (LEN=*) ::  mode   !<

    INTEGER(iwp) ::  i                   !< running index on x-axis
    INTEGER(iwp) ::  j                   !< running index on y-axis
    INTEGER(iwp) ::  k                   !< vertical index counter
    INTEGER(iwp) ::  lpr                 !< running index chem spcs
    INTEGER(iwp) ::  m                   !< running index for surface elements
!$  INTEGER(iwp) ::  omp_get_thread_num  !< intrinsic OMP function
    INTEGER(iwp) ::  sr                  !< statistical region
    INTEGER(iwp) ::  surf_e              !< end surface index
    INTEGER(iwp) ::  surf_s              !< start surface index
    INTEGER(iwp) ::  tn                  !< thread number

    REAL(wp)                                            ::  flag     !< topography masking flag
    REAL(wp), DIMENSION(nzb:nzt+1,0:threads_per_task-1) ::  sums_tmp !< temporary array used to sum-up profiles

    IF ( mode == 'profiles' )  THEN
!
!--    Sum-up profiles for the species
       tn = 0
       !$OMP PARALLEL PRIVATE( i, j, k, tn, lpr, sums_tmp )
       !$ tn = omp_get_thread_num()
       !$OMP DO
       DO  lpr = 1, cs_pr_count_sp
          sums_tmp(:,tn) = 0.0_wp
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt+1
                   sums_tmp(k,tn) = sums_tmp(k,tn) +                                               &
                        chem_species(cs_pr_index_sp(lpr))%conc(k,j,i) *                            &
                        rmask(j,i,sr)  *                                                           &
                        MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 22 ) )
                ENDDO
             ENDDO
          ENDDO
          sums_l(nzb:nzt+1,hom_index_spec(lpr),tn) = sums_tmp(nzb:nzt+1,tn)
       ENDDO
       !$OMP END PARALLEL
!
!--    Sum-up profiles for vertical fluxes of the the species. Note, in case of WS5 scheme the
!--    profiles of resolved-scale fluxes have been already summed-up, while resolved-scale fluxes
!--    need to be calculated in case of PW scheme.
!--    For summation employ a temporary array.
       !$OMP PARALLEL PRIVATE( i, j, k, tn, lpr, sums_tmp, flag )
       !$ tn = omp_get_thread_num()
       !$OMP DO
       DO  lpr = 1, cs_pr_count_fl_sgs
          sums_tmp(:,tn) = 0.0_wp
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = nzb, nzt
                   flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 23 ) ) *                &
                          MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 9  ) )
                   sums_tmp(k,tn) = sums_tmp(k,tn) -                                               &
                        0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )                                       &
                               * ( chem_species(cs_pr_index_fl_sgs(lpr))%conc(k+1,j,i) -           &
                                   chem_species(cs_pr_index_fl_sgs(lpr))%conc(k,j,i) ) *           &
                        ddzu(k+1) * rmask(j,i,sr)  * flag
                ENDDO
!
!--             Add surface fluxes (?Is the order mandatory or could it be done in one cycle?)
!--             Sum up only from upward-facing surfaces.
                surf_s = surf_def%start_index(j,i)
                surf_e = surf_def%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_def%k(m) + surf_def%koff(m)
                   sums_tmp(k,tn) = sums_tmp(k,tn) +                                               &
                                    MERGE( surf_def%cssws(cs_pr_index_fl_sgs(lpr),m), 0.0_wp,      &
                                           surf_def%upward(m) )
                ENDDO
                surf_s = surf_lsm%start_index(j,i)
                surf_e = surf_lsm%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_lsm%k(m) + surf_lsm%koff(m)
                   sums_tmp(k,tn) = sums_tmp(k,tn) +                                               &
                                    MERGE( surf_lsm%cssws(cs_pr_index_fl_sgs(lpr),m), 0.0_wp,      &
                                           surf_lsm%upward(m) )
                ENDDO
                surf_s = surf_usm%start_index(j,i)
                surf_e = surf_usm%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_usm%k(m) + surf_usm%koff(m)
                   sums_tmp(k,tn) = sums_tmp(k,tn) +                                               &
                                    MERGE( surf_usm%cssws(cs_pr_index_fl_sgs(lpr),m), 0.0_wp,      &
                                           surf_usm%upward(m) )
                ENDDO
             ENDDO
          ENDDO
          sums_l(nzb:nzt+1,hom_index_fl_sgs(lpr),tn) = sums_tmp(nzb:nzt+1,tn)
       ENDDO
       !$OMP END PARALLEL
!
!--    Resolved-scale fluxes from the WS5 scheme
       IF ( ws_scheme_sca )  THEN
          !$OMP PARALLEL PRIVATE( tn, lpr )
          !$ tn = omp_get_thread_num()
          !$OMP DO
          DO  lpr = 1, cs_pr_count_fl_res
             sums_l(nzb:nzt+1,hom_index_fl_res(lpr),tn) =                                          &
                                                    scalarflux_output_conversion(nzb:nzt+1) *      &
                                                    sums_ws_l(nzb:nzt+1,tn,cs_pr_index_fl_res(lpr))
          ENDDO
          !$OMP END PARALLEL
       ENDIF

    ELSEIF ( mode == 'time_series' )  THEN
!      @todo
    ENDIF

 END SUBROUTINE chem_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine for swapping of timelevels for chemical species called out from subroutine
!> swap_timelevel
!--------------------------------------------------------------------------------------------------!


 SUBROUTINE chem_swap_timelevel( level )


    INTEGER(iwp), INTENT(IN) ::  level
!
!-- Local variables
    INTEGER(iwp)             ::  lsp


    IF ( level == 0 )  THEN
       DO  lsp=1, nvar
          chem_species(lsp)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_1(:,:,:,lsp)
          chem_species(lsp)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_2(:,:,:,lsp)
       ENDDO
    ELSE
       DO  lsp=1, nvar
          chem_species(lsp)%conc(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    => spec_conc_2(:,:,:,lsp)
          chem_species(lsp)%conc_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  => spec_conc_1(:,:,:,lsp)
       ENDDO
    ENDIF

    RETURN
 END SUBROUTINE chem_swap_timelevel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to write restart data for chemistry model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_wrd_local


    INTEGER(iwp) ::  lsp  !< running index for chem spcs.

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       DO  lsp = 1, nspec
          CALL wrd_write_string( TRIM( chem_species(lsp)%name ) )
          WRITE ( 14 )  chem_species(lsp)%conc
          IF ( ALLOCATED( chem_species(lsp)%conc_av ) )  THEN
             CALL wrd_write_string( TRIM( chem_species(lsp)%name )//'_av' )
             WRITE ( 14 )  chem_species(lsp)%conc_av
          ENDIF
       ENDDO

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       DO  lsp = 1, nspec
          CALL wrd_mpi_io( TRIM( chem_species(lsp)%name ), chem_species(lsp)%conc )
          IF ( ALLOCATED( chem_species(lsp)%conc_av ) )  THEN
             CALL wrd_mpi_io( TRIM( chem_species(lsp)%name ) // '_av', chem_species(lsp)%conc_av )
          ENDIF
       ENDDO

    ENDIF

 END SUBROUTINE chem_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to calculate the deposition of gases and PMs. For now deposition only takes place on
!> lsm and usm horizontal upward faceing surfaces. Default surfaces are NOT considered.
!> The deposition of particlesis derived following Zhang et al., 2001, gases are deposited using
!> the DEPAC module (van Zanten et al., 2010).
!>
!> @TODO: Consider deposition on vertical surfaces
!> @TODO: Consider overlaying horizontal surfaces
!> @TODO: Consider resolved vegetation
!> @TODO: Check error messages
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_depo

!
!-- List of names of possible tracers.
    CHARACTER(LEN=*), PARAMETER ::  pspecnames(nposp) =                                            &
                                    (/ 'NO2           ',                       &    !< NO2
                                       'NO            ',                       &    !< NO
                                       'O3            ',                       &    !< O3
                                       'CO            ',                       &    !< CO
                                       'form          ',                       &    !< FORM
                                       'ald           ',                       &    !< ALD
                                       'pan           ',                       &    !< PAN
                                       'mgly          ',                       &    !< MGLY
                                       'par           ',                       &    !< PAR
                                       'ole           ',                       &    !< OLE
                                       'eth           ',                       &    !< ETH
                                       'tol           ',                       &    !< TOL
                                       'cres          ',                       &    !< CRES
                                       'xyl           ',                       &    !< XYL
                                       'SO4a_f        ',                       &    !< SO4a_f
                                       'SO2           ',                       &    !< SO2
                                       'HNO2          ',                       &    !< HNO2
                                       'CH4           ',                       &    !< CH4
                                       'NH3           ',                       &    !< NH3
                                       'NO3           ',                       &    !< NO3
                                       'OH            ',                       &    !< OH
                                       'HO2           ',                       &    !< HO2
                                       'N2O5          ',                       &    !< N2O5
                                       'SO4a_c        ',                       &    !< SO4a_c
                                       'NH4a_f        ',                       &    !< NH4a_f
                                       'NO3a_f        ',                       &    !< NO3a_f
                                       'NO3a_c        ',                       &    !< NO3a_c
                                       'C2O3          ',                       &    !< C2O3
                                       'XO2           ',                       &    !< XO2
                                       'XO2N          ',                       &    !< XO2N
                                       'cro           ',                       &    !< CRO
                                       'HNO3          ',                       &    !< HNO3
                                       'H2O2          ',                       &    !< H2O2
                                       'iso           ',                       &    !< ISO
                                       'ispd          ',                       &    !< ISPD
                                       'to2           ',                       &    !< TO2
                                       'open          ',                       &    !< OPEN
                                       'terp          ',                       &    !< TERP
                                       'ec_f          ',                       &    !< EC_f
                                       'ec_c          ',                       &    !< EC_c
                                       'pom_f         ',                       &    !< POM_f
                                       'pom_c         ',                       &    !< POM_c
                                       'ppm_f         ',                       &    !< PPM_f
                                       'ppm_c         ',                       &    !< PPM_c
                                       'na_ff         ',                       &    !< Na_ff
                                       'na_f          ',                       &    !< Na_f
                                       'na_c          ',                       &    !< Na_c
                                       'na_cc         ',                       &    !< Na_cc
                                       'na_ccc        ',                       &    !< Na_ccc
                                       'dust_ff       ',                       &    !< dust_ff
                                       'dust_f        ',                       &    !< dust_f
                                       'dust_c        ',                       &    !< dust_c
                                       'dust_cc       ',                       &    !< dust_cc
                                       'dust_ccc      ',                       &    !< dust_ccc
                                       'tpm10         ',                       &    !< tpm10
                                       'tpm25         ',                       &    !< tpm25
                                       'tss           ',                       &    !< tss
                                       'tdust         ',                       &    !< tdust
                                       'tc            ',                       &    !< tc
                                       'tcg           ',                       &    !< tcg
                                       'tsoa          ',                       &    !< tsoa
                                       'tnmvoc        ',                       &    !< tnmvoc
                                       'SOxa          ',                       &    !< SOxa
                                       'NOya          ',                       &    !< NOya
                                       'NHxa          ',                       &    !< NHxa
                                       'NO2_obs       ',                       &    !< NO2_obs
                                       'tpm10_biascorr',                       &    !< tpm10_biascorr
                                       'tpm25_biascorr',                       &    !< tpm25_biascorr
                                       'O3_biascorr   ' /)                          !< O3_biascorr

    INTEGER(iwp) ::  day_of_year  !< current day of the year
    INTEGER(iwp) ::  i            !< grid index in x-direction
    INTEGER(iwp) ::  i_pspec      !< index for matching depac gas component
    INTEGER(iwp) ::  j            !< grid index in y-direction
    INTEGER(iwp) ::  k            !< matching k to surface m at i,j
    INTEGER(iwp) ::  lsp          !< running index for chem spcs.
    INTEGER(iwp) ::  luv_palm     !< index of PALM LSM vegetation_type at current surface element
    INTEGER(iwp) ::  lup_palm     !< index of PALM LSM pavement_type at current surface element
    INTEGER(iwp) ::  luw_palm     !< index of PALM LSM water_type at current surface element
    INTEGER(iwp) ::  luu_palm     !< index of PALM USM walls/roofs at current surface element
    INTEGER(iwp) ::  lug_palm     !< index of PALM USM green walls/roofs at current surface element
    INTEGER(iwp) ::  lud_palm     !< index of PALM USM windows at current surface element
    INTEGER(iwp) ::  luv_dep      !< matching DEPAC LU to luv_palm
    INTEGER(iwp) ::  lup_dep      !< matching DEPAC LU to lup_palm
    INTEGER(iwp) ::  luw_dep      !< matching DEPAC LU to luw_palm
    INTEGER(iwp) ::  luu_dep      !< matching DEPAC LU to luu_palm
    INTEGER(iwp) ::  lug_dep      !< matching DEPAC LU to lug_palm
    INTEGER(iwp) ::  lud_dep      !< matching DEPAC LU to lud_palm
    INTEGER(iwp) ::  m            !< index for horizontal surfaces
    INTEGER(iwp) ::  mm           !< running index for horizontal surfaces
    INTEGER(iwp) ::  pspec        !< running index
!
!-- Vegetation (assign PALM classes to DEPAC land use classes).
    INTEGER(iwp) ::  ind_luv_user = 0             !<  ERROR as no class given in PALM
    INTEGER(iwp) ::  ind_luv_b_soil = 1           !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_luv_mixed_crops = 2      !<  assigned to ilu_arable
    INTEGER(iwp) ::  ind_luv_s_grass = 3          !<  assigned to ilu_grass
    INTEGER(iwp) ::  ind_luv_ev_needle_trees = 4  !<  assigned to ilu_coniferous_forest
    INTEGER(iwp) ::  ind_luv_de_needle_trees = 5  !<  assigned to ilu_coniferous_forest
    INTEGER(iwp) ::  ind_luv_ev_broad_trees = 6   !<  assigned to ilu_tropical_forest
    INTEGER(iwp) ::  ind_luv_de_broad_trees = 7   !<  assigned to ilu_deciduous_forest
    INTEGER(iwp) ::  ind_luv_t_grass = 8          !<  assigned to ilu_grass
    INTEGER(iwp) ::  ind_luv_desert = 9           !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_luv_tundra = 10          !<  assigned to ilu_other
    INTEGER(iwp) ::  ind_luv_irr_crops = 11       !<  assigned to ilu_arable
    INTEGER(iwp) ::  ind_luv_semidesert = 12      !<  assigned to ilu_other
    INTEGER(iwp) ::  ind_luv_ice = 13             !<  assigned to ilu_ice
    INTEGER(iwp) ::  ind_luv_marsh = 14           !<  assigned to ilu_other
    INTEGER(iwp) ::  ind_luv_ev_shrubs = 15       !<  assigned to ilu_mediterrean_scrub
    INTEGER(iwp) ::  ind_luv_de_shrubs = 16       !<  assigned to ilu_mediterrean_scrub
    INTEGER(iwp) ::  ind_luv_mixed_forest = 17    !<  assigned to ilu_coniferous_forest(ave(decid+conif))
    INTEGER(iwp) ::  ind_luv_intrup_forest = 18   !<  assigned to ilu_other (ave(other+decid))
!
!-- Water.
    INTEGER(iwp) ::  ind_luw_user = 0             !<  ERROR as no class given in PALM
    INTEGER(iwp) ::  ind_luw_lake = 1             !<  assigned to ilu_water_inland
    INTEGER(iwp) ::  ind_luw_river = 2            !<  assigned to ilu_water_inland
    INTEGER(iwp) ::  ind_luw_ocean = 3            !<  assigned to ilu_water_sea
    INTEGER(iwp) ::  ind_luw_pond = 4             !<  assigned to ilu_water_inland
    INTEGER(iwp) ::  ind_luw_fountain = 5         !<  assigned to ilu_water_inland
!
!-- Pavement
    INTEGER(iwp) ::  ind_lup_user = 0             !<  ERROR as no class given in PALM
    INTEGER(iwp) ::  ind_lup_asph_conc = 1        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_asph = 2             !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_conc = 3             !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_sett = 4             !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_pav_stones = 5       !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_cobblest = 6         !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_metal = 7            !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_wood = 8             !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_gravel = 9           !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_f_gravel = 10        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_pebblest = 11        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_woodchips = 12       !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_tartan = 13          !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_art_turf = 14        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_clay = 15            !<  assigned to ilu_desert
!
!-- Particle parameters according to the respective aerosol classes (PM25, PM10).
    INTEGER(iwp) ::  ind_p_size = 1  !< index for partsize in particle_pars
    INTEGER(iwp) ::  ind_p_dens = 2  !< index for rhopart in particle_pars
    INTEGER(iwp) ::  ind_p_slip = 3  !< index for slipcor in particle_pars
    INTEGER(iwp) ::  nwet            !< wetness indicator dor DEPAC; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
    INTEGER(iwp) ::  part_type       !< index for particle type (PM10 or PM25) in particle_pars

    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  lup_dep_v   !< lsm
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  luv_dep_v   !<
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  luw_dep_v   !<
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  luu_dep_v   !< usm
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  lug_dep_v   !<
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  lud_dep_v   !<
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  mvec        !< pre-computed index of upward-facing surface
    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  nwet_v      !<

    LOGICAL ::  match_lsm  !< flag indicating natural-type surface
    LOGICAL ::  match_usm  !< flag indicating urban-type surface

    LOGICAL, DIMENSION(nys:nyn,nxl:nxr) ::  do_pav_green    !< flag array indicating pavement or greened walls
    LOGICAL, DIMENSION(nys:nyn,nxl:nxr) ::  do_veg_wall     !< flag array indicating vegetation or walls
    LOGICAL, DIMENSION(nys:nyn,nxl:nxr) ::  do_wat_win      !< flag array indicating water or windows
    LOGICAL, DIMENSION(nys:nyn,nxl:nxr) ::  sai_present     !<

    REAL(wp) ::  bud          !< overall budget at current surface element
    REAL(wp) ::  dens         !< density at layer k at i,j
    REAL(wp) ::  dh           !< vertical grid size
    REAL(wp) ::  diffusivity  !< diffusivity
    REAL(wp) ::  dt_chem      !< length of chem time step
    REAL(wp) ::  dt_dh        !< dt_chem/dh
    REAL(wp) ::  inv_dh       !< inverse of vertical grid size
    REAL(wp) ::  lai          !< leaf area index at current surface element
    REAL(wp) ::  qv_tmp       !< surface mixing ratio at current surface element
    REAL(wp) ::  r_aero_surf  !< aerodynamic resistance (s/m) at current surface element
    REAL(wp) ::  rb           !< quasi-laminar boundary layer resistance (s/m)
    REAL(wp) ::  rc_tot       !< total canopy resistance (s/m)
    REAL(wp) ::  rh_surf      !< relative humidity at current surface element
    REAL(wp) ::  rs           !< Sedimentaion resistance (s/m)
    REAL(wp) ::  slinnfac     !<
    REAL(wp) ::  solar_rad    !< solar radiation, direct and diffuse, at current surface element
    REAL(wp) ::  temp_tmp     !< temperatur at i,j,k
    REAL(wp) ::  ts           !< surface temperatur in degrees celsius
    REAL(wp) ::  ustar_surf   !< ustar at current surface element
    REAL(wp) ::  vd_lu        !< deposition velocity (m/s)
    REAL(wp) ::  visc         !< Viscosity
    REAL(wp) ::  vs           !< Sedimentation velocity
    REAL(wp) ::  z0h_surf     !< roughness length for heat at current surface element

    REAL(wp), DIMENSION(nspec) ::  ccomp_tot  !< total compensation point (ug/m3), for now kept to zero for all species!

    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  dens_v         !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  inv_dh_v       !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  lai_v          !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  rb_v           !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  rc_tot_v       !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  rh_surf_v      !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  r_aero_surf_v  !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  sai_v          !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  solar_rad_v    !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  temp_tmp_v     !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  ts_v           !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  ustar_surf_v   !<
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::  visc_v         !<

    REAL(wp), DIMENSION(nys:nyn,nxl:nxr,nspec) ::  bud_lud_v  !< budget for USM windows at current surface element
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr,nspec) ::  bud_lug_v  !< budget for USM green surfaces at current surface element
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr,nspec) ::  bud_lup_v  !< budget for LSM pavement type at current surface element
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr,nspec) ::  bud_luu_v  !< budget for USM walls/roofs at current surface element
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr,nspec) ::  bud_luv_v  !< budget for LSM vegetation type at current surface element
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr,nspec) ::  bud_luw_v  !< budget for LSM water type at current surface element
!
!-- Particle parameters (PM10 (1), PM25 (2)) partsize (diameter in m), rhopart (density in kg/m3),
!-- slipcor (slip correction factor dimensionless, Seinfeld and Pandis 2006, Table 9.3).
    REAL(wp), DIMENSION(1:3,1:2), PARAMETER ::  particle_pars = RESHAPE( (/                        &
                                                                   8.0E-6_wp, 1.14E3_wp, 1.016_wp, &
                                                                   0.7E-6_wp, 1.14E3_wp, 1.082_wp  &
                                                                         /), (/ 3, 2 /) )

!
!-- Get current day of the year.
    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year )
!
!-- LSM or USM horizintal upward facing surface present at i,j. Note, default surfaces are not
!-- considered for deposition.
!-- First, check if upward-facing LSM surface exist and find its index. If more than one
!-- upward-facing surface exist, take the last one.
    mvec = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
             IF ( surf_lsm%upward(mm) )  mvec(j,i) = mm
          ENDDO
       ENDDO
    ENDDO

    do_veg_wall  = .FALSE.
    luv_dep_v    = 0
    do_pav_green = .FALSE.
    lup_dep_v    = 0
    do_wat_win   = .FALSE.
    luw_dep_v    = 0
    sai_present = .FALSE.

!
!-- Initialize budgets.
    bud_luv_v = 0.0_wp
    bud_lup_v = 0.0_wp
    bud_luw_v = 0.0_wp

    DO  i = nxl, nxr
       DO  j = nys, nyn
          m = mvec(j,i)

          match_lsm = ( m /= 0 )
!
!--       For LSM surfaces
          IF ( match_lsm )  THEN

             k = surf_lsm%k(m)
!
!--          Get needed variables for surface element m.
             ustar_surf  = surf_lsm%us(m)
             z0h_surf    = surf_lsm%z0h(m)
             r_aero_surf = surf_lsm%r_a(m)
             solar_rad   = surf_lsm%rad_sw_dir(m) + surf_lsm%rad_sw_dif(m)
             solar_rad_v(j,i) = solar_rad

             lai_v(j,i) = surf_lsm%lai(m)
             sai_v(j,i) = lai + 1
             sai_present(j,i) = .TRUE.
!
!--          For small grid spacing neglect R_a.
             IF ( dzw(k) <= 1.0 )  THEN
                r_aero_surf = 0.0_wp
             ENDIF
             r_aero_surf_v(j,i) = r_aero_surf
!
!--          Initialize lu's
             luv_palm = 0
             luv_dep = 0
             lup_palm = 0
             lup_dep = 0
             luw_palm = 0
             luw_dep = 0
!
!--          Get land use for i,j and assign to DEPAC lu.
             IF ( surf_lsm%frac(m,ind_veg_wall) > 0 )  THEN
                luv_palm = surf_lsm%vegetation_type(m)
                IF ( luv_palm == ind_luv_user )  THEN
                   message_string = 'no lsm-vegetation type defined'
                   CALL message( 'chem_depo', 'CHM0019', 1, 2, 0, 6, 0 )
                ELSEIF ( luv_palm == ind_luv_b_soil )  THEN
                   luv_dep = 9
                ELSEIF ( luv_palm == ind_luv_mixed_crops )  THEN
                   luv_dep = 2
                ELSEIF ( luv_palm == ind_luv_s_grass )  THEN
                   luv_dep = 1
                ELSEIF ( luv_palm == ind_luv_ev_needle_trees )  THEN
                   luv_dep = 4
                ELSEIF ( luv_palm == ind_luv_de_needle_trees )  THEN
                   luv_dep = 4
                ELSEIF ( luv_palm == ind_luv_ev_broad_trees )  THEN
                   luv_dep = 12
                ELSEIF ( luv_palm == ind_luv_de_broad_trees )  THEN
                   luv_dep = 5
                ELSEIF ( luv_palm == ind_luv_t_grass )  THEN
                   luv_dep = 1
                ELSEIF ( luv_palm == ind_luv_desert )  THEN
                   luv_dep = 9
                ELSEIF ( luv_palm == ind_luv_tundra )  THEN
                   luv_dep = 8
                ELSEIF ( luv_palm == ind_luv_irr_crops )  THEN
                   luv_dep = 2
                ELSEIF ( luv_palm == ind_luv_semidesert )  THEN
                   luv_dep = 8
                ELSEIF ( luv_palm == ind_luv_ice )  THEN
                   luv_dep = 10
                ELSEIF ( luv_palm == ind_luv_marsh )  THEN
                   luv_dep = 8
                ELSEIF ( luv_palm == ind_luv_ev_shrubs )  THEN
                   luv_dep = 14
                ELSEIF ( luv_palm == ind_luv_de_shrubs )  THEN
                   luv_dep = 14
                ELSEIF ( luv_palm == ind_luv_mixed_forest )  THEN
                   luv_dep = 4
                ELSEIF ( luv_palm == ind_luv_intrup_forest )  THEN
                   luv_dep = 8
                ENDIF
                do_veg_wall(j,i) = .TRUE.
             ENDIF
             luv_dep_v(j,i) = luv_dep

             IF ( surf_lsm%frac(m,ind_pav_green) > 0 )  THEN
                lup_palm = surf_lsm%pavement_type(m)
                IF ( lup_palm == ind_lup_user )  THEN
                   message_string = 'no lsm-pavement type defined'
                   CALL message( 'chem_depo', 'CHM0019', 1, 2, 0, 6, 0 )
!
!--             On Pavement, lup_dep is always 9.
                ELSEIF ( lup_palm == ind_lup_asph_conc )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_asph )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_conc )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_sett )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_pav_stones )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_cobblest )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_metal )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_wood )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_gravel )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_f_gravel )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_pebblest )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_woodchips )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_tartan )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_art_turf )  THEN
                   lup_dep = 9
                ELSEIF ( lup_palm == ind_lup_clay )  THEN
                   lup_dep = 9
                ENDIF
                do_pav_green(j,i) = .TRUE.
             ENDIF
             lup_dep_v(j,i) = lup_dep

             IF ( surf_lsm%frac(m,ind_wat_win) > 0 )  THEN
                luw_palm = surf_lsm%water_type(m)
                IF ( luw_palm == ind_luw_user )  THEN
                   message_string = 'no lsm-water type defined'
                   CALL message( 'chem_depo', 'CHM0019', 1, 2, 0, 6, 0 )
                ELSEIF ( luw_palm ==  ind_luw_lake )  THEN
                   luw_dep = 13
                ELSEIF ( luw_palm == ind_luw_river )  THEN
                   luw_dep = 13
                ELSEIF ( luw_palm == ind_luw_ocean )  THEN
                   luw_dep = 6
                ELSEIF ( luw_palm == ind_luw_pond )  THEN
                   luw_dep = 13
                ELSEIF ( luw_palm == ind_luw_fountain )  THEN
                   luw_dep = 13
                ENDIF
                do_wat_win(j,i) = .TRUE.
             ENDIF
             luw_dep_v(j,i) = luw_dep

!
!--          Set wetness indicator to dry or wet for lsm vegetation or pavement.
             IF ( surf_lsm%c_liq(m) > 0 )  THEN
                nwet = 1
             ELSE
                nwet = 0
             ENDIF
             nwet_v(j,i) = nwet
!
!--          Compute length of time step
             IF ( call_chem_at_all_substeps )  THEN
                dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
             ELSE
                dt_chem = dt_3d
             ENDIF

             dh = dzw(k)
             inv_dh = 1.0_wp / dh
             dt_dh = dt_chem / dh

!--          Temperature at i,j,k.
             temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp
             temp_tmp_v(j,i) = temp_tmp

             ts        = temp_tmp - 273.15_wp  !< in degrees celcius
             ts_v(j,i) = temp_tmp - 273.15_wp  !< in degrees celcius
!
!--          Viscosity of air.
             visc = 1.496E-6_wp * temp_tmp**1.5_wp / (temp_tmp + 120.0_wp )
             visc_v(j,i) = visc
!
!--          Air density at k.
             dens = rho_air_zw(k)
             dens_v(j,i) = dens
!
!--          Calculate relative humidity from specific humidity for DEPAC.
             qv_tmp = MAX( q(k,j,i), 0.0_wp )
             rh_surf = relativehumidity_from_specifichumidity(qv_tmp, temp_tmp, hyp(k) )
             rh_surf_v(j,i) = rh_surf
          ENDIF
          inv_dh_v(j,i)    = inv_dh
       ENDDO
    ENDDO
!
!-- Check if surface fraction (vegetation, pavement or water) > 0 and calculate vd and budget
!-- for each surface fraction. Then derive overall budget taking into account the surface fractions.
!-- Vegetation treatment.
    IF ( ANY( do_veg_wall ) )  THEN

       slinnfac = 1.0_wp
!
!--    Get deposition velocity vd.
       DO  lsp = 1, nvar
!
!--       Initialize values.
          vs     = 0.0_wp
          vd_lu  = 0.0_wp
          rs     = 0.0_wp
          rb     = 0.0_wp
          rc_tot = 0.0_wp
!
!--       Particulate matter.
          IF ( spc_names(lsp) == 'PM10'  .OR.  spc_names(lsp) == 'PM25' )  THEN
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
             ENDIF

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_lsm = ( m /= 0 )

                   IF ( match_lsm  .AND.  do_veg_wall(j,i) )  THEN
!
!--                   Compute sedimentation velocity.
                      vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),&
                                                              particle_pars(ind_p_size, part_type),&
                                                              particle_pars(ind_p_slip, part_type),&
                                                              visc_v(j,i))

                      ustar_surf  = surf_lsm%us(m)
                      r_aero_surf = surf_lsm%r_a(m)

                      CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                   &
                                                  particle_pars(ind_p_size, part_type),            &
                                                  particle_pars(ind_p_slip, part_type),            &
                                                  nwet_v(j,i), temp_tmp_v(j,i), dens_v(j,i),       &
                                                  visc_v(j,i), luv_dep_v(j,i), r_aero_surf,        &
                                                  ustar_surf )

                      k  = surf_lsm%k(m)
                      dh = dzw(k)
                      dt_dh = dt_chem / dh
                      bud_luv_v(j,i,lsp) = -chem_species(lsp)%conc(k,j,i) *                        &
                                           ( 1.0_wp - EXP( - vd_lu * dt_dh ) ) * dh

                   ENDIF
                ENDDO
             ENDDO
!
!--       Gases.
          ELSE
!
!--          Read spc_name of current species for gas parameter.
             IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                i_pspec = 0
                DO  pspec = 1, nposp
                   IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                      i_pspec = pspec
                   ENDIF
                ENDDO

             ELSE
!
!--             For now, other species are not deposited.
                CYCLE
             ENDIF
!
!--          Diffusivity for DEPAC relevant gases. Use default value.
             diffusivity = 0.11E-4_wp
!
!--          Overwrite with known coefficients of diffusivity from Massman (1998).
             IF ( spc_names(lsp) == 'NO2' )  diffusivity = 0.136E-4_wp
             IF ( spc_names(lsp) == 'NO'  )  diffusivity = 0.199E-4_wp
             IF ( spc_names(lsp) == 'O3'  )  diffusivity = 0.144E-4_wp
             IF ( spc_names(lsp) == 'CO'  )  diffusivity = 0.176E-4_wp
             IF ( spc_names(lsp) == 'SO2' )  diffusivity = 0.112E-4_wp
             IF ( spc_names(lsp) == 'CH4' )  diffusivity = 0.191E-4_wp
             IF ( spc_names(lsp) == 'NH3' )  diffusivity = 0.191E-4_wp

!
!--          Get quasi-laminar boundary layer resistance rb:
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_lsm = ( m /= 0 )
                   IF ( match_lsm  .AND.  do_veg_wall(j,i) )  THEN
!
!--                   No vegetation on bare soil, desert or ice.
                      luv_palm = surf_lsm%vegetation_type(m)
                      IF ( ( luv_palm == ind_luv_b_soil )  .OR.                                    &
                           ( luv_palm == ind_luv_desert )  .OR.                                    &
                           ( luv_palm == ind_luv_ice    ) )  THEN

                         lai_v(j,i) = 0.0_wp
                         sai_v(j,i) = 0.0_wp
                         sai_present(j,i) = .FALSE.

                      ENDIF
                      ustar_surf  = surf_lsm%us(m)
                      ustar_surf_v(j,i)  = surf_lsm%us(m)
                      z0h_surf    = surf_lsm%z0h(m)
                      solar_rad   = surf_lsm%rad_sw_dir(m) + surf_lsm%rad_sw_dif(m)
                      solar_rad_v(j,i) = solar_rad

!--                   Temperature at i,j,k.
                      temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp
                      ts       = temp_tmp - 273.15_wp  ! in degrees celcius

                      luv_dep = luv_dep_v(j,i)
!

!--                   Get quasi-laminar boundary layer resistance rb.
                      CALL get_rb_cell( ( luv_dep == ilu_water_sea )  .OR.                         &
                                        ( luv_dep == ilu_water_inland ), z0h_surf, ustar_surf,     &
                                        diffusivity, rb_v(j,i) )
                   ENDIF
                ENDDO
             ENDDO
!
!--          Get rc_tot.
             CALL drydepos_gas_depac( spc_names(lsp), mvec, ts_v, ustar_surf_v, sai_present,       &
                                      solar_rad_v, cos_zenith, rh_surf_v, lai_v, sai_v, nwet_v,    &
                                      luv_dep_v, 2, rc_tot_v, ccomp_tot(lsp), hyp(nzb),            &
                                      diffusivity )

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_lsm = ( m /= 0 )
                   IF ( match_lsm  .AND.  do_veg_wall(j,i) )  THEN
!
!--                   Calculate budget
                      IF ( rc_tot_v(j,i) <= 0.0 )  THEN
                         bud_luv_v(j,i,lsp) = 0.0_wp
                      ELSE
                         k  = surf_lsm%k(m)
                         dh = dzw(k)
                         dt_dh = dt_chem / dh
                         vd_lu = 1.0_wp / ( r_aero_surf_v(j,i) + rb_v(j,i) + rc_tot_v(j,i) )
                         bud_luv_v(j,i,lsp) = -( chem_species(lsp)%conc(k,j,i) - ccomp_tot(lsp) ) *&
                                               ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh
                      ENDIF

                   ENDIF
                ENDDO
             ENDDO

          ENDIF
       ENDDO
    ENDIF
!
!-- Pavement.
    IF ( ANY( do_pav_green ) )  THEN
!
!--    No vegetation on pavements:
       lai_v = 0.0_wp
       sai_v = 0.0_wp
       sai_present = .FALSE.

       slinnfac = 1.0_wp
!
!--    Get vd.
       DO  lsp = 1, nvar
!
!--       Initialize.
          vs     = 0.0_wp
          vd_lu  = 0.0_wp
          rs     = 0.0_wp
          rb     = 0.0_wp
          rc_tot = 0.0_wp
!
!--       Particulate matter.
          IF ( spc_names(lsp) == 'PM10'  .OR.  spc_names(lsp) == 'PM25' )  THEN
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
             ENDIF

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_lsm = ( m /= 0 )
                   IF ( match_lsm  .AND.  do_pav_green(j,i) )  THEN
!
!--                   Sedimentation velocity.
                      vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),&
                                                              particle_pars(ind_p_size, part_type),&
                                                              particle_pars(ind_p_slip, part_type),&
                                                              visc_v(j,i))

                      ustar_surf  = surf_lsm%us(m)

                      CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                   &
                                                  particle_pars(ind_p_size, part_type),            &
                                                  particle_pars(ind_p_slip, part_type),            &
                                                  nwet_v(j,i), temp_tmp_v(j,i), dens_v(j,i),       &
                                                  visc_v(j,i), lup_dep_v(j,i),                     &
                                                  r_aero_surf_v(j,i), ustar_surf )

                      k  = surf_lsm%k(m)
                      dh = dzw(k)
                      dt_dh = dt_chem / dh
                      bud_lup_v(j,i,lsp) = -chem_species(lsp)%conc(k,j,i) *                        &
                                           ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

                   ENDIF
                ENDDO
             ENDDO
!
!--       Gases.
          ELSE
!
!--          Read spc_name of current species for gas parameter.
             IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                i_pspec = 0
                DO  pspec = 1, nposp
                   IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                      i_pspec = pspec
                   ENDIF
                ENDDO

             ELSE
!
!--             For now, other species are not deposited.
                CYCLE
             ENDIF
!
!--          Diffusivity for DEPAC relevant gases. Use default value.
             diffusivity = 0.11E-4_wp
!
!--          Overwrite with known coefficients of diffusivity from Massman (1998)
             IF ( spc_names(lsp) == 'NO2' )  diffusivity = 0.136E-4_wp
             IF ( spc_names(lsp) == 'NO'  )  diffusivity = 0.199E-4_wp
             IF ( spc_names(lsp) == 'O3'  )  diffusivity = 0.144E-4_wp
             IF ( spc_names(lsp) == 'CO'  )  diffusivity = 0.176E-4_wp
             IF ( spc_names(lsp) == 'SO2' )  diffusivity = 0.112E-4_wp
             IF ( spc_names(lsp) == 'CH4' )  diffusivity = 0.191E-4_wp
             IF ( spc_names(lsp) == 'NH3' )  diffusivity = 0.191E-4_wp

             lai_v = 0.0_wp
             sai_v = 0.0_wp
             sai_present = .FALSE.

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_lsm = ( m /= 0 )
                   IF ( match_lsm .AND. do_pav_green(j,i) )  THEN

                      ustar_surf  = surf_lsm%us(m)
                      ustar_surf_v(j,i)  = surf_lsm%us(m)
                      z0h_surf    = surf_lsm%z0h(m)
                      solar_rad   = surf_lsm%rad_sw_dir(m) + surf_lsm%rad_sw_dif(m)
                      solar_rad_v(j,i) = solar_rad

!--                   Temperature at i,j,k.
                      temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp
                      ts       = temp_tmp - 273.15_wp  ! in degrees celcius

                      lup_dep = lup_dep_v(j,i)
!
!--                   Get quasi-laminar boundary layer resistance rb.
                      CALL get_rb_cell( ( lup_dep == ilu_water_sea )  .OR.                         &
                                        ( lup_dep == ilu_water_inland ), z0h_surf, ustar_surf,     &
                                          diffusivity, rb_v(j,i) )
                   ENDIF
                ENDDO
             ENDDO

!
!--          Get rc_tot.
             CALL drydepos_gas_depac( spc_names(lsp), mvec, ts_v, ustar_surf_v, sai_present,       &
                                      solar_rad_v, cos_zenith, rh_surf_v, lai_v, sai_v, nwet_v,    &
                                      lup_dep_v, 2, rc_tot_v, ccomp_tot(lsp), hyp(nzb),            &
                                      diffusivity )

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_lsm = ( m /= 0 )

                   IF ( match_lsm  .AND.  do_pav_green(j,i) )  THEN
!
!--                   Calculate budget.
                      IF ( rc_tot_v(j,i) <= 0.0 )  THEN
                         bud_lup_v(j,i,lsp) = 0.0_wp
                      ELSE
                         k  = surf_lsm%k(m)
                         dh = dzw(k)
                         dt_dh = dt_chem / dh

                         vd_lu = 1.0_wp / ( r_aero_surf_v(j,i) + rb_v(j,i) + rc_tot_v(j,i) )
                         bud_lup_v(j,i,lsp) = -( chem_species(lsp)%conc(k,j,i) - ccomp_tot(lsp) ) *&
                                               ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh
                      ENDIF

                   ENDIF
                ENDDO
             ENDDO

          ENDIF
       ENDDO
    ENDIF
!
!-- Water.
    bud_luw_v = 0.0
    IF ( ANY( do_wat_win ) )  THEN
!
!--    No vegetation on water:
       lai_v = 0.0_wp
       sai_v = 0.0_wp
       sai_present = .FALSE.

       slinnfac    = 1.0_wp
!
!--    Get vd
       DO  lsp = 1, nvar
!
!--       Initialize
          vs     = 0.0_wp
          vd_lu  = 0.0_wp
          rs     = 0.0_wp
          rb     = 0.0_wp
          rc_tot = 0.0_wp
!
!--       Particulate matter.
          IF ( spc_names(lsp) == 'PM10' .OR.  spc_names(lsp) == 'PM25' )  THEN
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
             ENDIF

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_lsm = ( m /= 0 )
                   IF ( match_lsm  .AND.  do_wat_win(j,i) )  THEN
!
!--                   Sedimentation velocity:
                      vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),&
                                                              particle_pars(ind_p_size, part_type),&
                                                              particle_pars(ind_p_slip, part_type),&
                                                              visc_v(j,i))

                      ustar_surf  = surf_lsm%us(m)

                      CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                   &
                                                  particle_pars(ind_p_size, part_type),            &
                                                  particle_pars(ind_p_slip, part_type),            &
                                                  nwet_v(j,i), temp_tmp_v(j,i), dens_v(j,i),       &
                                                  visc_v(j,i), luw_dep_v(j,i),                     &
                                                  r_aero_surf_v(j,i), ustar_surf )

                      k  = surf_lsm%k(m)
                      dh = dzw(k)
                      dt_dh = dt_chem / dh
                      bud_luw_v(j,i,lsp) = -chem_species(lsp)%conc(k,j,i) *                        &
                                           ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

                   ENDIF
                ENDDO
             ENDDO
!
!--       Gases.
          ELSE
!
!--          Read spc_name of current species for gas PARAMETER.
             IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                i_pspec = 0
                DO  pspec = 1, nposp
                   IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                      i_pspec = pspec
                   ENDIF
                ENDDO

             ELSE
!
!--             For now, other species are not deposited.
                CYCLE
             ENDIF
!
!--          Diffusivity for DEPAC relevant gases. Use default value.
             diffusivity = 0.11E-4_wp
!
!--          Overwrite with known coefficients of diffusivity from Massman (1998).
             IF ( spc_names(lsp) == 'NO2' )  diffusivity = 0.136E-4_wp
             IF ( spc_names(lsp) == 'NO'  )  diffusivity = 0.199E-4_wp
             IF ( spc_names(lsp) == 'O3'  )  diffusivity = 0.144E-4_wp
             IF ( spc_names(lsp) == 'CO'  )  diffusivity = 0.176E-4_wp
             IF ( spc_names(lsp) == 'SO2' )  diffusivity = 0.112E-4_wp
             IF ( spc_names(lsp) == 'CH4' )  diffusivity = 0.191E-4_wp
             IF ( spc_names(lsp) == 'NH3' )  diffusivity = 0.191E-4_wp


             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_lsm = ( m /= 0 )
                   IF ( match_lsm .AND. do_pav_green(j,i) )  THEN

                      ustar_surf  = surf_lsm%us(m)
                      ustar_surf_v(j,i)  = surf_lsm%us(m)
                      z0h_surf    = surf_lsm%z0h(m)
                      solar_rad   = surf_lsm%rad_sw_dir(m) + surf_lsm%rad_sw_dif(m)
                      solar_rad_v(j,i) = solar_rad

!--                   Temperature at i,j,k.
                      temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp
                      ts       = temp_tmp - 273.15_wp  ! in degrees celcius

                      luw_dep = luw_dep_v(j,i)
!
!--                   Get quasi-laminar boundary layer resistance rb:
                      CALL get_rb_cell( ( luw_dep == ilu_water_sea )  .OR.                         &
                                        ( luw_dep == ilu_water_inland ), z0h_surf, ustar_surf,     &
                                          diffusivity, rb_v(j,i) )
                   ENDIF
                ENDDO
             ENDDO
!
!--          Get rc_tot.
             CALL drydepos_gas_depac( spc_names(lsp), mvec, ts_v, ustar_surf_v, sai_present,       &
                                      solar_rad_v, cos_zenith, rh_surf_v, lai_v, sai_v, nwet_v,    &
                                      luw_dep_v, 2, rc_tot_v, ccomp_tot(lsp), hyp(nzb),            &
                                      diffusivity )

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_lsm = ( m /= 0 )
                   IF ( match_lsm  .AND.  do_wat_win(j,i) )  THEN
!
!--                   Calculate budget.
                      IF ( rc_tot_v(j,i) <= 0.0 )  THEN
                         bud_luw_v(j,i,lsp) = 0.0_wp
                      ELSE
                         k  = surf_lsm%k(m)
                         dh = dzw(k)
                         dt_dh = dt_chem / dh

                         vd_lu = 1.0_wp / ( r_aero_surf_v(j,i) + rb_v(j,i) + rc_tot_v(j,i) )
                         bud_luw_v(j,i,lsp) = -( chem_species(lsp)%conc(k,j,i) - ccomp_tot(lsp) ) *&
                                               ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh
                      ENDIF

                   ENDIF
                ENDDO
             ENDDO

          ENDIF
       ENDDO
    ENDIF
!
!-- Calculate overall budget for surface m and adapt concentration.
    DO  i = nxl, nxr
!NEC$ ivdep
       DO  j = nys, nyn
          m = mvec(j,i)

          match_lsm = ( m /= 0 )
!
!--       For LSM surfaces.
          IF ( match_lsm )  THEN
             DO  lsp = 1, nspec

                k  = surf_lsm%k(m)

                bud = surf_lsm%frac(m,ind_veg_wall)  * bud_luv_v(j,i,lsp) +                        &
                      surf_lsm%frac(m,ind_pav_green) * bud_lup_v(j,i,lsp) +                        &
                      surf_lsm%frac(m,ind_wat_win)   * bud_luw_v(j,i,lsp)
!
!--             Compute new concentration.
                chem_species(lsp)%conc(k,j,i) = MAX( 0.0_wp, chem_species(lsp)%conc(k,j,i) +       &
                                                             bud * inv_dh_v(j,i) )

             ENDDO

          ENDIF
       ENDDO
    ENDDO
!
!-- For USM surfaces.
!-- First, check if upward-facing USM surface exist and find its index. If more than one
!-- upward-facing surface exist, take the last one.
    mvec = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
             IF ( surf_usm%upward(mm) )  mvec(j,i) = mm
          ENDDO
       ENDDO
    ENDDO

    do_veg_wall  = .FALSE.
    luu_dep_v    = 0
    do_pav_green = .FALSE.
    lug_dep_v    = 0
    do_wat_win   = .FALSE.
    lud_dep_v    = 0
    sai_present  = .FALSE.

!
!-- Initialize budgets.
    bud_luu_v = 0.0_wp
    bud_lug_v = 0.0_wp
    bud_lud_v = 0.0_wp

    DO  i = nxl, nxr
       DO  j = nys, nyn
          m = mvec(j,i)

          match_usm = ( m /= 0 )

          IF ( match_usm )  THEN
             k = surf_usm%k(m)
!
!--          Get needed variables for surface element m.
             ustar_surf  = surf_usm%us(m)
             z0h_surf    = surf_usm%z0h(m)
             r_aero_surf = surf_usm%r_a(m)
             solar_rad   = surf_usm%rad_sw_dir(m) + surf_usm%rad_sw_dif(m)
             solar_rad_v(j,i) = solar_rad
             lai_v(j,i) = surf_usm%lai(m)
             sai_v(j,i) = lai + 1
             sai_present(j,i) = .TRUE.

!--          For small grid spacing neglect r_a.
             IF ( dzw(k) <= 1.0_wp )  THEN
                r_aero_surf = 0.0_wp
             ENDIF
             r_aero_surf_v(j,i) = r_aero_surf
!
!--          Initialize lu's.
             luu_palm = 0
             luu_dep  = 0
             lug_palm = 0
             lug_dep  = 0
             lud_palm = 0
             lud_dep  = 0
!
!--          Assign DEPAC lu for green wall surfaces.
             IF ( surf_usm%frac(m,ind_pav_green) > 0 )  THEN
!
!--             For green urban surfaces, i.e. green roofs assume LU short grass.
                lug_palm = ind_luv_s_grass
                lug_dep = 1
                do_pav_green(j,i) = .TRUE.
             ENDIF
             lug_dep_v(j,i) = lug_dep
!
!--          Wall surfaces.
             IF ( surf_usm%frac(m,ind_veg_wall) > 0 )  THEN
!
!--             For walls in USM assume concrete walls/roofs,
!--             assumed LU class desert as also assumed for pavements in LSM.
                luu_palm = ind_lup_conc
                luu_dep = 9
                do_veg_wall(j,i) = .TRUE.
             ENDIF
             luu_dep_v(j,i) = luu_dep
!
!--          Window surfaces.
             IF ( surf_usm%frac(m,ind_wat_win) > 0 )  THEN
!
!--             For windows in USM assume metal as this is as close as we get,
!--             assumed LU class desert as also assumed for pavements in LSM.
                lud_palm = ind_lup_metal
                lud_dep = 9
                do_wat_win(j,i) = .TRUE.
             ENDIF
             lud_dep_v(j,i) = lud_dep

!
!--          @TODO: Activate these lines as soon as new ebsolver branch is merged:
!--          Set wetness indicator to dry or wet for usm vegetation or pavement.
             !IF ( surf_usm%c_liq(m) > 0 )  THEN
             !   nwet = 1
             !ELSE
             nwet = 0
             !ENDIF
             nwet_v(j,i) = nwet
!
!--          Compute length of time step.
             IF ( call_chem_at_all_substeps )  THEN
                dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
             ELSE
                dt_chem = dt_3d
             ENDIF

             dh = dzw(k)
             inv_dh = 1.0_wp / dh
             dt_dh = dt_chem / dh
!
!--          Temperature at i,j,k
             temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp
             temp_tmp_v(j,i) = temp_tmp

             ts       = temp_tmp - 273.15_wp  !< in degrees celcius
             ts_v(j,i) = temp_tmp - 273.15_wp  !< in degrees celcius
!
!--          Viscosity of air
             visc = 1.496E-6_wp * temp_tmp**1.5_wp / ( temp_tmp + 120.0_wp )
             visc_v(j,i) = visc
!
!--          Air density at k
             dens = rho_air_zw(k)
             dens_v(j,i) = dens
!
!--          Calculate relative humidity from specific humidity for DEPAC
             qv_tmp = MAX( q(k,j,i), 0.0_wp )
             rh_surf = relativehumidity_from_specifichumidity( qv_tmp, temp_tmp, hyp(k) )
             rh_surf_v(j,i) = rh_surf
          ENDIF
          inv_dh_v(j,i) = inv_dh
       ENDDO
    ENDDO
!
!-- Check if surface fraction (vegetation, pavement or water) > 0 and calculate vd and budget
!-- for each surface fraction. Then derive overall budget taking into account the surface
!-- fractions.
!-- Walls/roofs
    IF ( ANY( do_veg_wall ) )  THEN

       slinnfac = 1.0_wp
!
!--    Get deposition velocity vd.
       DO  lsp = 1, nvar
!
!--       Initialize
          vs     = 0.0_wp
          vd_lu  = 0.0_wp
          rs     = 0.0_wp
          rb     = 0.0_wp
          rc_tot = 0.0_wp
!
!--       Particulate matter.
          IF ( spc_names(lsp) == 'PM10'  .OR.  spc_names(lsp) == 'PM25' )  THEN
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
             ENDIF

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_usm = ( m /= 0 )
                   IF ( match_usm  .AND.  do_veg_wall(j,i) )  THEN
!
!--                   Sedimentation velocity:
                      vs = slinnfac *                                                              &
                           sedimentation_velocity( particle_pars(ind_p_dens, part_type),           &
                                                   particle_pars(ind_p_size, part_type),           &
                                                   particle_pars(ind_p_slip, part_type),           &
                                                   visc_v(j,i))

                      ustar_surf  = surf_usm%us(m)
                      r_aero_surf = surf_usm%r_a(m)

                      CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                   &
                                                  particle_pars(ind_p_size, part_type),            &
                                                  particle_pars(ind_p_slip, part_type),            &
                                                  nwet_v(j,i), temp_tmp_v(j,i), dens_v(j,i),       &
                                                  visc_v(j,i), luu_dep_v(j,i),                     &
                                                  r_aero_surf, ustar_surf )

                      k  = surf_usm%k(m)
                      dh = dzw(k)
                      dt_dh = dt_chem / dh
                      bud_luu_v(j,i,lsp) = -chem_species(lsp)%conc(k,j,i) *                        &
                                           ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

                   ENDIF
                ENDDO
             ENDDO
!
!--       Gases.
          ELSE
!
!--          Read spc_name of current species for gas parameter.
             IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                i_pspec = 0
                DO  pspec = 1, nposp
                   IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                      i_pspec = pspec
                   ENDIF
                ENDDO
             ELSE
!
!--             For now, other species are not deposited.
                CYCLE
             ENDIF
!
!--          Diffusivity for DEPAC relevant gases. Use default value.
             diffusivity = 0.11E-4_wp
!
!--          Overwrite with known coefficients of diffusivity from Massman (1998).
             IF ( spc_names(lsp) == 'NO2' )  diffusivity = 0.136E-4_wp
             IF ( spc_names(lsp) == 'NO'  )  diffusivity = 0.199E-4_wp
             IF ( spc_names(lsp) == 'O3'  )  diffusivity = 0.144E-4_wp
             IF ( spc_names(lsp) == 'CO'  )  diffusivity = 0.176E-4_wp
             IF ( spc_names(lsp) == 'SO2' )  diffusivity = 0.112E-4_wp
             IF ( spc_names(lsp) == 'CH4' )  diffusivity = 0.191E-4_wp
             IF ( spc_names(lsp) == 'NH3' )  diffusivity = 0.191E-4_wp
!
!--          Get quasi-laminar boundary layer resistance rb.
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_usm = ( m /= 0 )
                   IF ( match_usm  .AND.  do_veg_wall(j,i) )  THEN

                      lai_v(j,i) = 0.0_wp
                      sai_v(j,i) = 0.0_wp
                      sai_present(j,i) = .FALSE.

                      ustar_surf  = surf_usm%us(m)
                      ustar_surf_v(j,i)  = surf_usm%us(m)
                      z0h_surf    = surf_usm%z0h(m)
                      solar_rad   = surf_usm%rad_sw_dir(m) + surf_usm%rad_sw_dif(m)
                      solar_rad_v(j,i) = solar_rad
!
!--                   Temperature at i,j,k.
                      temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp
                      ts       = temp_tmp - 273.15_wp  ! in degrees celcius

                      luu_dep = luu_dep_v(j,i)
!
!--                   Get quasi-laminar boundary layer resistance rb.
                      CALL get_rb_cell( ( luu_dep == ilu_water_sea )  .OR.                         &
                                        ( luu_dep == ilu_water_inland ), z0h_surf, ustar_surf,     &
                                          diffusivity, rb_v(j,i) )
                   ENDIF
                ENDDO
             ENDDO
!
!--          Get rc_tot.
             CALL drydepos_gas_depac( spc_names(lsp), mvec, ts_v, ustar_surf_v, sai_present,       &
                                      solar_rad_v, cos_zenith, rh_surf_v, lai_v, sai_v, nwet_v,    &
                                      luu_dep_v, 2, rc_tot_v, ccomp_tot(lsp), hyp(nzb),            &
                                      diffusivity )

!
!--          Calculate budget.
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_usm = ( m /= 0 )
                   IF ( match_usm  .AND.  do_veg_wall(j,i) )  THEN
!
!--                   Calculate budget
                      IF ( rc_tot_v(j,i) <= 0.0_wp )  THEN
                         bud_luu_v(j,i,lsp) = 0.0_wp
                      ELSE
                         k  = surf_usm%k(m)
                         dh = dzw(k)
                         dt_dh = dt_chem / dh
                         vd_lu = 1.0_wp / ( r_aero_surf_v(j,i) + rb_v(j,i) + rc_tot_v(j,i) )
                         bud_luu_v(j,i,lsp) = -( chem_species(lsp)%conc(k,j,i) - ccomp_tot(lsp) ) *&
                                               ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh
                      ENDIF

                   ENDIF
                ENDDO
             ENDDO

          ENDIF
       ENDDO
    ENDIF

!
!-- Green usm surfaces.
    IF ( ANY( do_pav_green ) )  THEN

       slinnfac = 1.0_wp
!
!--    Get vd.
       DO  lsp = 1, nvar
!
!--       Initialize.
          vs     = 0.0_wp
          vd_lu  = 0.0_wp
          rs     = 0.0_wp
          rb     = 0.0_wp
          rc_tot = 0.0_wp
!
!--       Particulate matter.
          IF ( spc_names(lsp) == 'PM10'  .OR.  spc_names(lsp) == 'PM25')  THEN
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
             ENDIF

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_usm = ( m /= 0 )
                   IF ( match_usm  .AND.  do_pav_green(j,i) )  THEN
!
!--                   Sedimentation velocity:
                      vs = slinnfac *                                                              &
                           sedimentation_velocity( particle_pars(ind_p_dens, part_type),           &
                                                   particle_pars(ind_p_size, part_type),           &
                                                   particle_pars(ind_p_slip, part_type),           &
                                                   visc_v(j,i))

                      ustar_surf  = surf_usm%us(m)

                      CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                   &
                                                  particle_pars(ind_p_size, part_type),            &
                                                  particle_pars(ind_p_slip, part_type),            &
                                                  nwet_v(j,i), temp_tmp_v(j,i), dens_v(j,i),       &
                                                  visc_v(j,i), lug_dep_v(j,i),                     &
                                                  r_aero_surf_v(j,i), ustar_surf )

                      k  = surf_usm%k(m)
                      dh = dzw(k)
                      dt_dh = dt_chem / dh
                      bud_lug_v(j,i,lsp) = -chem_species(lsp)%conc(k,j,i) *                        &
                                           ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

                   ENDIF
                ENDDO
             ENDDO
!
!--       Gases.
          ELSE
!
!--          Read spc_name of current species for gas parameter.
             IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                i_pspec = 0
                DO  pspec = 1, nposp
                   IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                      i_pspec = pspec
                   ENDIF
                ENDDO
             ELSE
!
!--             For now, other species are not deposited.
                CYCLE
             ENDIF
!
!--          Diffusivity for DEPAC relevant gases. Use default value.
             diffusivity = 0.11E-4_wp
!
!--          Overwrite with known coefficients of diffusivity from Massman (1998).
             IF ( spc_names(lsp) == 'NO2' )  diffusivity = 0.136E-4_wp
             IF ( spc_names(lsp) == 'NO'  )  diffusivity = 0.199E-4_wp
             IF ( spc_names(lsp) == 'O3'  )  diffusivity = 0.144E-4_wp
             IF ( spc_names(lsp) == 'CO'  )  diffusivity = 0.176E-4_wp
             IF ( spc_names(lsp) == 'SO2' )  diffusivity = 0.112E-4_wp
             IF ( spc_names(lsp) == 'CH4' )  diffusivity = 0.191E-4_wp
             IF ( spc_names(lsp) == 'NH3' )  diffusivity = 0.191E-4_wp
!
!--          Get quasi-laminar boundary layer resistance rb.
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_usm = ( m /= 0 )
                   IF ( match_usm  .AND.  do_pav_green(j,i) )  THEN

                      lug_palm = ind_luv_s_grass
                      IF ( ( lug_palm == ind_luv_b_soil )  .OR.                                    &
                           ( lug_palm == ind_luv_desert )  .OR.                                    &
                           ( lug_palm == ind_luv_ice ) ) THEN

                         lai_v(j,i) = 0.0_wp
                         sai_v(j,i) = 0.0_wp
                         sai_present(j,i) = .FALSE.

                      ENDIF

                      ustar_surf  = surf_usm%us(m)
                      ustar_surf_v(j,i)  = surf_usm%us(m)
                      z0h_surf    = surf_usm%z0h(m)
                      solar_rad   = surf_usm%rad_sw_dir(m) + surf_usm%rad_sw_dif(m)
                      solar_rad_v(j,i) = solar_rad

!--                   Temperature at i,j,k.
                      temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp
                      ts       = temp_tmp - 273.15_wp  ! in degrees celcius

                      lug_dep = lug_dep_v(j,i)
!
!--                   Get quasi-laminar boundary layer resistance rb.
                      CALL get_rb_cell( ( lug_dep == ilu_water_sea )  .OR.                         &
                                        ( lug_dep == ilu_water_inland ), z0h_surf, ustar_surf,     &
                                          diffusivity, rb_v(j,i) )
                   ENDIF
                ENDDO
             ENDDO
!
!--          Get rc_tot
             CALL drydepos_gas_depac( spc_names(lsp), mvec, ts_v, ustar_surf_v, sai_present,       &
                                      solar_rad_v, cos_zenith, rh_surf_v, lai_v, sai_v, nwet_v,    &
                                      lug_dep_v, 2, rc_tot_v, ccomp_tot(lsp), hyp(nzb),            &
                                      diffusivity )

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_usm = ( m /= 0 )
                   IF ( match_usm  .AND.  do_pav_green(j,i) )  THEN
!
!--                   Calculate budget.
                      IF ( rc_tot_v(j,i) <= 0.0_wp )  THEN
                         bud_lug_v(j,i,lsp) = 0.0_wp
                      ELSE
                         k  = surf_usm%k(m)
                         dh = dzw(k)
                         dt_dh = dt_chem / dh

                         vd_lu = 1.0_wp / ( r_aero_surf_v(j,i) + rb_v(j,i) + rc_tot_v(j,i) )
                         bud_lug_v(j,i,lsp) = -( chem_species(lsp)%conc(k,j,i) - ccomp_tot(lsp) ) *&
                                               ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh
                      ENDIF

                   ENDIF
                ENDDO
             ENDDO

          ENDIF
       ENDDO
    ENDIF
!
!-- Windows.
    IF ( ANY( do_wat_win) )  THEN

       slinnfac = 1.0_wp
!
!--    Get vd.
       DO  lsp = 1, nvar
!
!--       Initialize.
          vs     = 0.0_wp
          vd_lu  = 0.0_wp
          rs     = 0.0_wp
          rb     = 0.0_wp
          rc_tot = 0.0_wp
!
!--       Particulate matter.
          IF ( spc_names(lsp) == 'PM10'  .OR.  spc_names(lsp) == 'PM25')  THEN
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
             ENDIF

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_usm = ( m /= 0 )
                   IF ( match_usm  .AND.  do_wat_win(j,i) )  THEN
!
!--                   Sedimentation velocity.
                      vs = slinnfac *                                                              &
                           sedimentation_velocity( particle_pars(ind_p_dens, part_type),           &
                                                   particle_pars(ind_p_size, part_type),           &
                                                   particle_pars(ind_p_slip, part_type),           &
                                                   visc_v(j,i))

                      ustar_surf  = surf_usm%us(m)

                      CALL drydepo_aero_zhang_vd( vd_lu, rs,                                       &
                                                  vs,                                              &
                                                  particle_pars(ind_p_size, part_type),            &
                                                  particle_pars(ind_p_slip, part_type),            &
                                                  nwet_v(j,i), temp_tmp_v(j,i), dens_v(j,i),       &
                                                  visc_v(j,i), lud_dep_v(j,i),                     &
                                                  r_aero_surf_v(j,i), ustar_surf )

                      k  = surf_usm%k(m)
                      dh = dzw(k)
                      dt_dh = dt_chem / dh
                      bud_lud_v(j,i,lsp) = -chem_species(lsp)%conc(k,j,i) *                        &
                                           ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

                   ENDIF
                ENDDO
             ENDDO
!
!--       Gases.
          ELSE
!
!--          Read spc_name of current species for gas PARAMETER.
             IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                i_pspec = 0
                DO  pspec = 1, nposp
                   IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                      i_pspec = pspec
                   END IF
                ENDDO
             ELSE
!
!--             For now, other species are not deposited.
                CYCLE
             ENDIF
!
!--          Diffusivity for DEPAC relevant gases. Use default value.
             diffusivity = 0.11E-4_wp
!
!--          Overwrite with known coefficients of diffusivity from Massman (1998).
             IF ( spc_names(lsp) == 'NO2' )  diffusivity = 0.136E-4_wp
             IF ( spc_names(lsp) == 'NO'  )  diffusivity = 0.199E-4_wp
             IF ( spc_names(lsp) == 'O3'  )  diffusivity = 0.144E-4_wp
             IF ( spc_names(lsp) == 'CO'  )  diffusivity = 0.176E-4_wp
             IF ( spc_names(lsp) == 'SO2' )  diffusivity = 0.112E-4_wp
             IF ( spc_names(lsp) == 'CH4' )  diffusivity = 0.191E-4_wp
             IF ( spc_names(lsp) == 'NH3' )  diffusivity = 0.191E-4_wp
!
!--          Get quasi-laminar boundary layer resistance rb,
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_usm = ( m /= 0 )
                   IF ( match_usm  .AND.  do_wat_win(j,i) )  THEN
!
!--                   No vegetation on water.
                      lai_v(j,i) = 0.0_wp
                      sai_v(j,i) = 0.0_wp
                      sai_present(j,i) = .FALSE.

                      ustar_surf  = surf_usm%us(m)
                      ustar_surf_v(j,i)  = surf_usm%us(m)
                      z0h_surf    = surf_usm%z0h(m)
                      solar_rad   = surf_usm%rad_sw_dir(m) + surf_usm%rad_sw_dif(m)
                      solar_rad_v(j,i) = solar_rad

!--                   Temperature at i,j,k
                      temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp
                      ts       = temp_tmp - 273.15_wp  ! in degrees celcius

                      lud_dep = lud_dep_v(j,i)
!
!--                   Get quasi-laminar boundary layer resistance rb.
                      CALL get_rb_cell( ( lud_dep == ilu_water_sea )  .OR.                         &
                                        ( lud_dep == ilu_water_inland ), z0h_surf, ustar_surf,     &
                                          diffusivity, rb_v(j,i) )
                   ENDIF
                ENDDO
             ENDDO
!
!--          Get rc_tot.
             CALL drydepos_gas_depac( spc_names(lsp), mvec, ts_v, ustar_surf_v, sai_present,       &
                                      solar_rad_v, cos_zenith, rh_surf_v, lai_v, sai_v, nwet_v,    &
                                      lud_dep_v, 2, rc_tot_v, ccomp_tot(lsp), hyp(nzb),            &
                                      diffusivity )

             DO  i = nxl, nxr
                DO  j = nys, nyn
                   m = mvec(j,i)

                   match_usm = ( m /= 0 )
                   IF ( match_usm .AND. do_wat_win(j,i) )  THEN
!
!--                   Calculate budget.
                      IF ( rc_tot_v(j,i) <= 0.0 )  THEN
                         bud_lud_v(j,i,lsp) = 0.0_wp
                      ELSE
                         k  = surf_usm%k(m)
                         dh = dzw(k)
                         dt_dh = dt_chem/dh

                         vd_lu = 1.0_wp / ( r_aero_surf_v(j,i) + rb_v(j,i) + rc_tot_v(j,i) )
                         bud_lud_v(j,i,lsp) = -( chem_species(lsp)%conc(k,j,i) - ccomp_tot(lsp) ) *&
                                               ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh
                      ENDIF

                   ENDIF
                ENDDO
             ENDDO

          ENDIF
       ENDDO
    ENDIF
!
!-- Calculate overall budget for surface m and adapt concentration.
    DO  i = nxl, nxr
!NEC$ ivdep
       DO  j = nys, nyn
          m = mvec(j,i)

          match_usm = ( m /= 0 )
!
!--       For USM surfaces.
          IF ( match_usm )  THEN

             DO  lsp = 1, nspec
                k  = surf_usm%k(m)

                bud = surf_usm%frac(m,ind_veg_wall)  * bud_luu_v(j,i,lsp) +                        &
                      surf_usm%frac(m,ind_pav_green) * bud_lug_v(j,i,lsp) +                        &
                      surf_usm%frac(m,ind_wat_win)   * bud_lud_v(j,i,lsp)
!
!--             Compute new concentration.
                chem_species(lsp)%conc(k,j,i) = MAX( 0.0_wp, chem_species(lsp)%conc(k,j,i) +       &
                                                             bud * inv_dh_v(j,i) )

             ENDDO
          ENDIF
       ENDDO
    ENDDO

 END SUBROUTINE chem_depo


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to calculate the deposition of gases and PMs. For now deposition only takes place on
!> lsm and usm horizontal upward faceing surfaces. Default surfaces are NOT considered.
!> The deposition of particlesis derived following Zhang et al., 2001, gases are deposited using
!>  the DEPAC module (van Zanten et al., 2010).
!>
!> @TODO: Consider deposition on vertical surfaces
!> @TODO: Consider overlaying horizontal surfaces
!> @TODO: Consider resolved vegetation
!> @TODO: Check error messages
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_depo_ij( i, j )

!
!-- List of names of possible tracers
    CHARACTER(LEN=*), PARAMETER ::  pspecnames(nposp) =                                            &
                                    (/ 'NO2           ',                       &    !< NO2
                                       'NO            ',                       &    !< NO
                                       'O3            ',                       &    !< O3
                                       'CO            ',                       &    !< CO
                                       'form          ',                       &    !< FORM
                                       'ald           ',                       &    !< ALD
                                       'pan           ',                       &    !< PAN
                                       'mgly          ',                       &    !< MGLY
                                       'par           ',                       &    !< PAR
                                       'ole           ',                       &    !< OLE
                                       'eth           ',                       &    !< ETH
                                       'tol           ',                       &    !< TOL
                                       'cres          ',                       &    !< CRES
                                       'xyl           ',                       &    !< XYL
                                       'SO4a_f        ',                       &    !< SO4a_f
                                       'SO2           ',                       &    !< SO2
                                       'HNO2          ',                       &    !< HNO2
                                       'CH4           ',                       &    !< CH4
                                       'NH3           ',                       &    !< NH3
                                       'NO3           ',                       &    !< NO3
                                       'OH            ',                       &    !< OH
                                       'HO2           ',                       &    !< HO2
                                       'N2O5          ',                       &    !< N2O5
                                       'SO4a_c        ',                       &    !< SO4a_c
                                       'NH4a_f        ',                       &    !< NH4a_f
                                       'NO3a_f        ',                       &    !< NO3a_f
                                       'NO3a_c        ',                       &    !< NO3a_c
                                       'C2O3          ',                       &    !< C2O3
                                       'XO2           ',                       &    !< XO2
                                       'XO2N          ',                       &    !< XO2N
                                       'cro           ',                       &    !< CRO
                                       'HNO3          ',                       &    !< HNO3
                                       'H2O2          ',                       &    !< H2O2
                                       'iso           ',                       &    !< ISO
                                       'ispd          ',                       &    !< ISPD
                                       'to2           ',                       &    !< TO2
                                       'open          ',                       &    !< OPEN
                                       'terp          ',                       &    !< TERP
                                       'ec_f          ',                       &    !< EC_f
                                       'ec_c          ',                       &    !< EC_c
                                       'pom_f         ',                       &    !< POM_f
                                       'pom_c         ',                       &    !< POM_c
                                       'ppm_f         ',                       &    !< PPM_f
                                       'ppm_c         ',                       &    !< PPM_c
                                       'na_ff         ',                       &    !< Na_ff
                                       'na_f          ',                       &    !< Na_f
                                       'na_c          ',                       &    !< Na_c
                                       'na_cc         ',                       &    !< Na_cc
                                       'na_ccc        ',                       &    !< Na_ccc
                                       'dust_ff       ',                       &    !< dust_ff
                                       'dust_f        ',                       &    !< dust_f
                                       'dust_c        ',                       &    !< dust_c
                                       'dust_cc       ',                       &    !< dust_cc
                                       'dust_ccc      ',                       &    !< dust_ccc
                                       'tpm10         ',                       &    !< tpm10
                                       'tpm25         ',                       &    !< tpm25
                                       'tss           ',                       &    !< tss
                                       'tdust         ',                       &    !< tdust
                                       'tc            ',                       &    !< tc
                                       'tcg           ',                       &    !< tcg
                                       'tsoa          ',                       &    !< tsoa
                                       'tnmvoc        ',                       &    !< tnmvoc
                                       'SOxa          ',                       &    !< SOxa
                                       'NOya          ',                       &    !< NOya
                                       'NHxa          ',                       &    !< NHxa
                                       'NO2_obs       ',                       &    !< NO2_obs
                                       'tpm10_biascorr',                       &    !< tpm10_biascorr
                                       'tpm25_biascorr',                       &    !< tpm25_biascorr
                                       'O3_biascorr   ' /)                          !< O3_biascorr

    INTEGER(iwp) ::  i            !< grid index in x-direction
    INTEGER(iwp) ::  i_pspec      !< index for matching depac gas component
    INTEGER(iwp) ::  j            !< grid index in y-direction
    INTEGER(iwp) ::  k            !< matching k to surface m at i,j
    INTEGER(iwp) ::  lsp          !< running index for chem spcs.
    INTEGER(iwp) ::  luv_palm     !< index of PALM LSM vegetation_type at current surface element
    INTEGER(iwp) ::  lup_palm     !< index of PALM LSM pavement_type at current surface element
    INTEGER(iwp) ::  luw_palm     !< index of PALM LSM water_type at current surface element
    INTEGER(iwp) ::  luu_palm     !< index of PALM USM walls/roofs at current surface element
    INTEGER(iwp) ::  lug_palm     !< index of PALM USM green walls/roofs at current surface element
    INTEGER(iwp) ::  lud_palm     !< index of PALM USM windows at current surface element
    INTEGER(iwp) ::  luv_dep      !< matching DEPAC LU to luv_palm
    INTEGER(iwp) ::  lup_dep      !< matching DEPAC LU to lup_palm
    INTEGER(iwp) ::  luw_dep      !< matching DEPAC LU to luw_palm
    INTEGER(iwp) ::  luu_dep      !< matching DEPAC LU to luu_palm
    INTEGER(iwp) ::  lug_dep      !< matching DEPAC LU to lug_palm
    INTEGER(iwp) ::  lud_dep      !< matching DEPAC LU to lud_palm
    INTEGER(iwp) ::  m            !< index for horizontal surfaces
    INTEGER(iwp) ::  mm           !< running index for horizontal surfaces
    INTEGER(iwp) ::  pspec        !< running index
!
!-- Vegetation.
!.. Assign PALM classes to DEPAC land use classes.
    INTEGER(iwp) ::  ind_luv_user = 0             !<  ERROR as no class given in PALM
    INTEGER(iwp) ::  ind_luv_b_soil = 1           !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_luv_mixed_crops = 2      !<  assigned to ilu_arable
    INTEGER(iwp) ::  ind_luv_s_grass = 3          !<  assigned to ilu_grass
    INTEGER(iwp) ::  ind_luv_ev_needle_trees = 4  !<  assigned to ilu_coniferous_forest
    INTEGER(iwp) ::  ind_luv_de_needle_trees = 5  !<  assigned to ilu_coniferous_forest
    INTEGER(iwp) ::  ind_luv_ev_broad_trees = 6   !<  assigned to ilu_tropical_forest
    INTEGER(iwp) ::  ind_luv_de_broad_trees = 7   !<  assigned to ilu_deciduous_forest
    INTEGER(iwp) ::  ind_luv_t_grass = 8          !<  assigned to ilu_grass
    INTEGER(iwp) ::  ind_luv_desert = 9           !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_luv_tundra = 10          !<  assigned to ilu_other
    INTEGER(iwp) ::  ind_luv_irr_crops = 11       !<  assigned to ilu_arable
    INTEGER(iwp) ::  ind_luv_semidesert = 12      !<  assigned to ilu_other
    INTEGER(iwp) ::  ind_luv_ice = 13             !<  assigned to ilu_ice
    INTEGER(iwp) ::  ind_luv_marsh = 14           !<  assigned to ilu_other
    INTEGER(iwp) ::  ind_luv_ev_shrubs = 15       !<  assigned to ilu_mediterrean_scrub
    INTEGER(iwp) ::  ind_luv_de_shrubs = 16       !<  assigned to ilu_mediterrean_scrub
    INTEGER(iwp) ::  ind_luv_mixed_forest = 17    !<  assigned to ilu_coniferous_forest(ave(decid+conif))
    INTEGER(iwp) ::  ind_luv_intrup_forest = 18   !<  assigned to ilu_other (ave(other+decid))
!
!-- Water.
    INTEGER(iwp) ::  ind_luw_user = 0             !<  ERROR as no class given in PALM
    INTEGER(iwp) ::  ind_luw_lake = 1             !<  assigned to ilu_water_inland
    INTEGER(iwp) ::  ind_luw_river = 2            !<  assigned to ilu_water_inland
    INTEGER(iwp) ::  ind_luw_ocean = 3            !<  assigned to ilu_water_sea
    INTEGER(iwp) ::  ind_luw_pond = 4             !<  assigned to ilu_water_inland
    INTEGER(iwp) ::  ind_luw_fountain = 5         !<  assigned to ilu_water_inland
!
!-- Pavement.
    INTEGER(iwp) ::  ind_lup_user = 0             !<  ERROR as no class given in PALM
    INTEGER(iwp) ::  ind_lup_asph_conc = 1        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_asph = 2             !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_conc = 3             !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_sett = 4             !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_pav_stones = 5       !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_cobblest = 6         !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_metal = 7            !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_wood = 8             !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_gravel = 9           !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_f_gravel = 10        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_pebblest = 11        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_woodchips = 12       !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_tartan = 13          !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_art_turf = 14        !<  assigned to ilu_desert
    INTEGER(iwp) ::  ind_lup_clay = 15            !<  assigned to ilu_desert
!
!-- Particle parameters according to the respective aerosol classes (PM25, PM10).
    INTEGER(iwp) ::  ind_p_size = 1  !< index for partsize in particle_pars
    INTEGER(iwp) ::  ind_p_dens = 2  !< index for rhopart in particle_pars
    INTEGER(iwp) ::  ind_p_slip = 3  !< index for slipcor in particle_pars
    INTEGER(iwp) ::  nwet            !< wetness indicator dor DEPAC; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow
    INTEGER(iwp) ::  part_type       !< index for particle type (PM10 or PM25) in particle_pars

    LOGICAL ::  match_lsm  !< flag indicating natural-type surface
    LOGICAL ::  match_usm  !< flag indicating urban-type surface

    REAL(wp) ::  dens           !< density at layer k at i,j
    REAL(wp) ::  dh             !< vertical grid size
    REAL(wp) ::  diffusivity    !< diffusivity
    REAL(wp) ::  dt_chem        !< length of chem time step
    REAL(wp) ::  dt_dh          !< dt_chem/dh
    REAL(wp) ::  inv_dh         !< inverse of vertical grid size
    REAL(wp) ::  lai            !< leaf area index at current surface element
    REAL(wp) ::  ppm2ugm3       !< conversion factor from ppm to ug/m3
    REAL(wp) ::  qv_tmp         !< surface mixing ratio at current surface element
    REAL(wp) ::  r_aero_surf    !< aerodynamic resistance (s/m) at current surface element
    REAL(wp) ::  rb             !< quasi-laminar boundary layer resistance (s/m)
    REAL(wp) ::  rc_tot         !< total canopy resistance (s/m)
    REAL(wp) ::  rh_surf        !< relative humidity at current surface element
    REAL(wp) ::  rs             !< Sedimentaion resistance (s/m)
    REAL(wp) ::  sai            !< surface area index at current surface element assumed to be lai + 1
    REAL(wp) ::  slinnfac       !<
    REAL(wp) ::  solar_rad      !< solar radiation, direct and diffuse, at current surface element
    REAL(wp) ::  temp_tmp       !< temperatur at i,j,k
    REAL(wp) ::  ts             !< surface temperatur in degrees celsius
    REAL(wp) ::  ustar_surf     !< ustar at current surface element
    REAL(wp) ::  vd_lu          !< deposition velocity (m/s)
    REAL(wp) ::  visc           !< Viscosity
    REAL(wp) ::  vs             !< Sedimentation velocity
    REAL(wp) ::  z0h_surf       !< roughness length for heat at current surface element

    REAL(wp), DIMENSION(nspec) ::  bud       !< overall budget at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_lud   !< budget for USM windows at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_lug   !< budget for USM green surfaces at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_lup   !< budget for LSM pavement type at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_luu   !< budget for USM walls/roofs at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_luv   !< budget for LSM vegetation type at current surface element
    REAL(wp), DIMENSION(nspec) ::  bud_luw   !< budget for LSM water type at current surface element
    REAL(wp), DIMENSION(nspec) ::  ccomp_tot !< total compensation point (ug/m3), for now kept to zero for all species!
    REAL(wp), DIMENSION(nspec) ::  conc_ijk  !< concentration at i,j,k

!
!-- Particle parameters (PM10 (1), PM25 (2)) partsize (diameter in m), rhopart (density in kg/m3),
!-- slipcor (slip correction factor dimensionless, Seinfeld and Pandis 2006, Table 9.3).
    REAL(wp), DIMENSION(1:3,1:2), PARAMETER ::  particle_pars = RESHAPE( (/                        &
                                                         8.0e-6_wp, 1.14e3_wp, 1.016_wp,           &   !<  1
                                                         0.7e-6_wp, 1.14e3_wp, 1.082_wp            &   !<  2
                                                                        /), (/ 3, 2 /) )!

!
!-- LSM or USM horizintal upward facing surface present at i,j. Note, default surfaces are not
!-- considered for deposition.
!-- First, check if upward-facing LSM surface exist and find its index. If more than one
!-- upward-facing surface exist, take the last one.
    m = 0
    DO  mm = surf_lsm%start_index(j,i), surf_lsm%end_index(j,i)
       IF ( surf_lsm%upward(mm) )  m = mm
    ENDDO
    match_lsm = ( m /= 0 )
!
!-- For LSM surfaces
    IF ( match_lsm )  THEN

       k = surf_lsm%k(m)
!
!--    Get needed variables for surface element m
       ustar_surf  = surf_lsm%us(m)
       z0h_surf    = surf_lsm%z0h(m)
       r_aero_surf = surf_lsm%r_a(m)
       solar_rad   = surf_lsm%rad_sw_dir(m) + surf_lsm%rad_sw_dif(m)

       lai = surf_lsm%lai(m)
       sai = lai + 1
!
!--    For small grid spacing neglect R_a
       IF ( dzw(k) <= 1.0 )  THEN
          r_aero_surf = 0.0_wp
       ENDIF
!
!--    Initialize lu's
       luv_palm = 0
       luv_dep = 0
       lup_palm = 0
       lup_dep = 0
       luw_palm = 0
       luw_dep = 0
!
!--    Initialize budgets
       bud_luv = 0.0_wp
       bud_lup = 0.0_wp
       bud_luw = 0.0_wp
!
!--    Get land use for i,j and assign to DEPAC lu
       IF ( surf_lsm%frac(m,ind_veg_wall) > 0 )  THEN
          luv_palm = surf_lsm%vegetation_type(m)
          IF ( luv_palm == ind_luv_user )  THEN
             message_string = 'no lsm-vegetation type defined'
             CALL message( 'chem_depo_ij', 'CHM0019', 1, 2, 0, 6, 0 )
          ELSEIF ( luv_palm == ind_luv_b_soil )  THEN
             luv_dep = 9
          ELSEIF ( luv_palm == ind_luv_mixed_crops )  THEN
             luv_dep = 2
          ELSEIF ( luv_palm == ind_luv_s_grass )  THEN
             luv_dep = 1
          ELSEIF ( luv_palm == ind_luv_ev_needle_trees )  THEN
             luv_dep = 4
          ELSEIF ( luv_palm == ind_luv_de_needle_trees )  THEN
             luv_dep = 4
          ELSEIF ( luv_palm == ind_luv_ev_broad_trees )  THEN
             luv_dep = 12
          ELSEIF ( luv_palm == ind_luv_de_broad_trees )  THEN
             luv_dep = 5
          ELSEIF ( luv_palm == ind_luv_t_grass )  THEN
             luv_dep = 1
          ELSEIF ( luv_palm == ind_luv_desert )  THEN
             luv_dep = 9
          ELSEIF ( luv_palm == ind_luv_tundra )  THEN
             luv_dep = 8
          ELSEIF ( luv_palm == ind_luv_irr_crops )  THEN
             luv_dep = 2
          ELSEIF ( luv_palm == ind_luv_semidesert )  THEN
             luv_dep = 8
          ELSEIF ( luv_palm == ind_luv_ice )  THEN
             luv_dep = 10
          ELSEIF ( luv_palm == ind_luv_marsh )  THEN
             luv_dep = 8
          ELSEIF ( luv_palm == ind_luv_ev_shrubs )  THEN
             luv_dep = 14
          ELSEIF ( luv_palm == ind_luv_de_shrubs )  THEN
             luv_dep = 14
          ELSEIF ( luv_palm == ind_luv_mixed_forest )  THEN
             luv_dep = 4
          ELSEIF ( luv_palm == ind_luv_intrup_forest )  THEN
             luv_dep = 8
          ENDIF
       ENDIF

       IF ( surf_lsm%frac(m,ind_pav_green) > 0 )  THEN
          lup_palm = surf_lsm%pavement_type(m)
          IF ( lup_palm == ind_lup_user )  THEN
             message_string = 'no lsm-pavement type defined'
             CALL message( 'chem_depo_ij', 'CHM0019', 1, 2, 0, 6, 0 )
          ELSEIF ( lup_palm == ind_lup_asph_conc )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_asph )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_conc )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_sett )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_pav_stones )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_cobblest )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_metal )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_wood )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_gravel )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_f_gravel )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_pebblest )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_woodchips )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_tartan )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_art_turf )  THEN
             lup_dep = 9
          ELSEIF ( lup_palm == ind_lup_clay )  THEN
             lup_dep = 9
          ENDIF
       ENDIF

       IF ( surf_lsm%frac(m,ind_wat_win) > 0 )  THEN
          luw_palm = surf_lsm%water_type(m)
          IF ( luw_palm == ind_luw_user )  THEN
             message_string = 'no lsm-water type defined'
             CALL message( 'chem_depo_ij', 'CHM0019', 1, 2, 0, 6, 0 )
          ELSEIF ( luw_palm ==  ind_luw_lake )  THEN
             luw_dep = 13
          ELSEIF ( luw_palm == ind_luw_river )  THEN
             luw_dep = 13
          ELSEIF ( luw_palm == ind_luw_ocean )  THEN
             luw_dep = 6
          ELSEIF ( luw_palm == ind_luw_pond )  THEN
             luw_dep = 13
          ELSEIF ( luw_palm == ind_luw_fountain )  THEN
             luw_dep = 13
          ENDIF
       ENDIF
!
!--    Set wetness indicator to dry or wet for lsm vegetation or pavement
       IF ( surf_lsm%c_liq(m) > 0 )  THEN
          nwet = 1
       ELSE
          nwet = 0
       ENDIF
!
!--    Compute length of time step
       IF ( call_chem_at_all_substeps )  THEN
          dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_chem = dt_3d
       ENDIF

       dh = dzw(k)
       inv_dh = 1.0_wp / dh
       dt_dh = dt_chem/dh
!
!--    Concentration at i,j,k
       DO  lsp = 1, nspec
          conc_ijk(lsp) = chem_species(lsp)%conc(k,j,i)
       ENDDO

!--    Temperature at i,j,k
       temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp

       ts       = temp_tmp - 273.15  !< in degrees celcius
!
!--    Viscosity of air
       visc = 1.496e-6 * temp_tmp**1.5 / (temp_tmp + 120.0)
!
!--    Air density at k
       dens = rho_air_zw(k)
!
!--    Calculate relative humidity from specific humidity for DEPAC
       qv_tmp = MAX( q(k,j,i), 0.0_wp)
       rh_surf = relativehumidity_from_specifichumidity( qv_tmp, temp_tmp, hyp(k) )
!
!-- Check if surface fraction (vegetation, pavement or water) > 0 and calculate vd and budget
!-- for each surface fraction. Then derive overall budget taking into account the surface fractions.
!
!--    Vegetation
       IF ( surf_lsm%frac(m,ind_veg_wall) > 0 )  THEN

!
!--       No vegetation on bare soil, desert or ice:
          IF ( ( luv_palm == ind_luv_b_soil )  .OR.  ( luv_palm == ind_luv_desert )  .OR.          &
               ( luv_palm == ind_luv_ice ) )  THEN

             lai = 0.0_wp
             sai = 0.0_wp

          ENDIF

          slinnfac = 1.0_wp
!
!--       Get deposition velocity vd
          DO  lsp = 1, nvar
!
!--          Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs,                                             &
                                            vs,                                                    &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc,                            &
                                            luv_dep,                                               &
                                            r_aero_surf, ustar_surf )

                bud_luv(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh


             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc,                            &
                                            luv_dep ,                                              &
                                            r_aero_surf, ustar_surf )

                bud_luv(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

             ELSE !< GASES
!
!--             Read spc_name of current species for gas parameter
                IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO

                ELSE
!
!--             For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--             ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9              xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                 c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at k:
                ppm2ugm3 =  (dens/xm_air) * 0.001_wp  !< (mole air)/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value
                diffusivity            = 0.11e-4
!
!--             Overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( ( luv_dep == ilu_water_sea )  .OR.                               &
                                  ( luv_dep == ilu_water_inland ), z0h_surf, ustar_surf,           &
                                  diffusivity, rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac_ij( spc_names(lsp), ts, ustar_surf, solar_rad, cos_zenith, &
                                            rh_surf, lai, sai, nwet, luv_dep, 2,                   &
                                            rc_tot, ccomp_tot(lsp), hyp(nzb), diffusivity )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN

                   bud_luv(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / ( r_aero_surf + rb + rc_tot )

                   bud_luv(lsp) = - ( conc_ijk(lsp) - ccomp_tot(lsp) ) *                           &
                                    ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF
!
!--    Pavement
       IF ( surf_lsm%frac(m,ind_pav_green) > 0 )  THEN
!
!--       No vegetation on pavements:
          lai = 0.0_wp
          sai = 0.0_wp

          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--       Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity:
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc, lup_dep,                   &
                                            r_aero_surf, ustar_surf )

                bud_lup(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh


             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity:
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc, lup_dep,                   &
                                            r_aero_surf, ustar_surf )

                bud_lup(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

             ELSE  !<GASES
!
!--             Read spc_name of current species for gas parameter
                IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO

                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9               xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                 c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at lowest layer:
                ppm2ugm3 = ( dens / xm_air ) * 0.001_wp  !< (mole air)/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value
                diffusivity = 0.11E-4
!
!--             Overwrite with known coefficients of diffusivity from Massman (1998).
                IF ( spc_names(lsp) == 'NO2' )  diffusivity = 0.136E-4
                IF ( spc_names(lsp) == 'NO'  )  diffusivity = 0.199E-4
                IF ( spc_names(lsp) == 'O3'  )  diffusivity = 0.144E-4
                IF ( spc_names(lsp) == 'CO'  )  diffusivity = 0.176E-4
                IF ( spc_names(lsp) == 'SO2' )  diffusivity = 0.112E-4
                IF ( spc_names(lsp) == 'CH4' )  diffusivity = 0.191E-4
                IF ( spc_names(lsp) == 'NH3' )  diffusivity = 0.191E-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( ( luv_dep == ilu_water_sea )  .OR.                               &
                                  ( luv_dep == ilu_water_inland ), z0h_surf, ustar_surf,           &
                                  diffusivity, rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac_ij( spc_names(lsp), ts, ustar_surf, solar_rad, cos_zenith, &
                                            rh_surf, lai, sai, nwet, lup_dep,2, rc_tot,            &
                                            ccomp_tot(lsp), hyp(nzb), diffusivity )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN
                   bud_lup(lsp) = 0.0_wp
                ELSE
                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )
                   bud_lup(lsp) = - ( conc_ijk(lsp) - ccomp_tot(lsp) ) *                           &
                                  ( 1.0_wp - EXP( - vd_lu * dt_dh ) ) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF
!
!--    Water
       IF ( surf_lsm%frac(m,ind_wat_win) > 0 )  THEN
!
!--       No vegetation on water:
          lai = 0.0_wp
          sai = 0.0_wp
          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--          Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity:
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc, luw_dep,                   &
                                            r_aero_surf, ustar_surf )

                bud_luw(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity:
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type), &
                                                        particle_pars(ind_p_slip, part_type), &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc, luw_dep,                   &
                                            r_aero_surf, ustar_surf )

                bud_luw(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

             ELSE  !<GASES
!
!--             Read spc_name of current species for gas PARAMETER
                IF ( ANY(pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO

                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9               xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                 c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at lowest layer:
                ppm2ugm3 = (dens/xm_air) * 0.001_wp  !< (mole air)/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value
                diffusivity            = 0.11e-4
!
!--             Overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( ( luw_dep == ilu_water_sea )  .OR.                               &
                                  ( luw_dep == ilu_water_inland ), z0h_surf, ustar_surf,           &
                                    diffusivity, rb )

!--             Get rc_tot
                CALL drydepos_gas_depac_ij( spc_names(lsp), ts, ustar_surf, solar_rad, cos_zenith, &
                                            rh_surf, lai, sai, nwet, luw_dep,                      &
                                            2, rc_tot, ccomp_tot(lsp), hyp(nzb), diffusivity )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN

                   bud_luw(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )

                   bud_luw(lsp) = - ( conc_ijk(lsp) - ccomp_tot(lsp) ) *                           &
                                  ( 1.0_wp - EXP(-vd_lu * dt_dh ) ) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF


       bud = 0.0_wp
!
!--    Calculate overall budget for surface m and adapt concentration
       DO  lsp = 1, nspec

          bud(lsp) = surf_lsm%frac(m,ind_veg_wall)  * bud_luv(lsp) +                               &
                     surf_lsm%frac(m,ind_pav_green) * bud_lup(lsp) +                               &
                     surf_lsm%frac(m,ind_wat_win)   * bud_luw(lsp)
!
!--       Compute new concentration:
          conc_ijk(lsp) = conc_ijk(lsp) + bud(lsp) * inv_dh

          chem_species(lsp)%conc(k,j,i) = MAX( 0.0_wp, conc_ijk(lsp) )

       ENDDO

    ENDIF
!
!-- For USM surfaces
!-- First, check if upward-facing USM surface exist and find its index. If more than one
!-- upward-facing surface exist, take the last one.
    m = 0
    DO  mm = surf_usm%start_index(j,i), surf_usm%end_index(j,i)
       IF ( surf_usm%upward(mm) )  m = mm
    ENDDO
    match_usm = ( m /= 0 )

    IF ( match_usm )  THEN
       k = surf_usm%k(m)
!
!--    Get needed variables for surface element m
       ustar_surf  = surf_usm%us(m)
       z0h_surf    = surf_usm%z0h(m)
       r_aero_surf = surf_usm%r_a(m)
       solar_rad   = surf_usm%rad_sw_dir(m) + surf_usm%rad_sw_dif(m)
       lai = surf_usm%lai(m)
       sai = lai + 1
!
!--    For small grid spacing neglect R_a
       IF ( dzw(k) <= 1.0 )  THEN
          r_aero_surf = 0.0_wp
       ENDIF
!
!--    Initialize lu's
       luu_palm = 0
       luu_dep = 0
       lug_palm = 0
       lug_dep = 0
       lud_palm = 0
       lud_dep = 0
!
!--    Initialize budgets
       bud_luu  = 0.0_wp
       bud_lug = 0.0_wp
       bud_lud = 0.0_wp
!
!--    Get land use for i,j and assign to DEPAC lu
       IF ( surf_usm%frac(m,ind_pav_green) > 0 )  THEN
!
!--       For green urban surfaces (e.g. green roofs assume LU short grass
          lug_palm = ind_luv_s_grass
          lug_dep = 1
       ENDIF

       IF ( surf_usm%frac(m,ind_veg_wall) > 0 )  THEN
!
!--       For walls in USM assume concrete walls/roofs,
!--       assumed LU class desert as also assumed for pavements in LSM
          luu_palm = ind_lup_conc
          luu_dep = 9
       ENDIF

       IF ( surf_usm%frac(m,ind_wat_win) > 0 )  THEN
!
!--       For windows in USM assume metal as this is as close as we get, assumed LU class desert as
!--       also assumed for pavements in LSM.
          lud_palm = ind_lup_metal
          lud_dep = 9
       ENDIF
!
!--    @TODO: Activate these lines as soon as new ebsolver branch is merged:
!--    Set wetness indicator to dry or wet for usm vegetation or pavement
       !IF ( surf_usm%c_liq(m) > 0 )  THEN
       !   nwet = 1
       !ELSE
       nwet = 0
       !ENDIF
!
!--    Compute length of time step
       IF ( call_chem_at_all_substeps )  THEN
          dt_chem = dt_3d * weight_pres(intermediate_timestep_count)
       ELSE
          dt_chem = dt_3d
       ENDIF

       dh = dzw(k)
       inv_dh = 1.0_wp / dh
       dt_dh = dt_chem / dh
!
!--    Concentration at i,j,k
       DO  lsp = 1, nspec
          conc_ijk(lsp) = chem_species(lsp)%conc(k,j,i)
       ENDDO
!
!--    Temperature at i,j,k
       temp_tmp = pt(k,j,i) * ( hyp(k) / 100000.0_wp )**0.286_wp

       ts       = temp_tmp - 273.15  !< in degrees celcius
!
!--    Viscosity of air
       visc = 1.496e-6 * temp_tmp**1.5 / ( temp_tmp + 120.0 )
!
!--    Air density at k
       dens = rho_air_zw(k)
!
!--    Calculate relative humidity from specific humidity for DEPAC
       qv_tmp = MAX( q(k,j,i), 0.0_wp )
       rh_surf = relativehumidity_from_specifichumidity( qv_tmp, temp_tmp, hyp(k) )
!
!--    Check if surface fraction (vegetation, pavement or water) > 0 and calculate vd and budget for
!--    each surface fraction. Then derive overall budget taking into account the surface fractions.

!--    Walls/roofs
       IF ( surf_usm%frac(m,ind_veg_wall) > 0 )  THEN
!
!--       No vegetation on non-green walls:
          lai = 0.0_wp
          sai = 0.0_wp

          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--          Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF (spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc, luu_dep,                   &
                                            r_aero_surf, ustar_surf )

                bud_luu(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( - vd_lu * dt_dh ) ) * dh

             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc, luu_dep ,                  &
                                            r_aero_surf, ustar_surf )

                bud_luu(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( - vd_lu * dt_dh ) ) * dh

             ELSE  !< GASES
!
!--             Read spc_name of current species for gas parameter
                IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO
                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9              xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                 c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at k:
                ppm2ugm3 =  (dens/xm_air) * 0.001_wp  !< (mole air)/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value
                diffusivity            = 0.11e-4
!
!--             Overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( ( luu_dep == ilu_water_sea )  .OR.                               &
                                  ( luu_dep == ilu_water_inland ), z0h_surf, ustar_surf,           &
                                    diffusivity, rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac_ij( spc_names(lsp), ts, ustar_surf, solar_rad, cos_zenith, &
                                            rh_surf, lai, sai, nwet, luu_dep,                      &
                                            2, rc_tot, ccomp_tot(lsp), hyp(nzb), diffusivity )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN

                   bud_luu(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )

                   bud_luu(lsp) = - ( conc_ijk(lsp) - ccomp_tot(lsp) ) *                           &
                                  ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF
!
!--    Green usm surfaces
       IF ( surf_usm%frac(m,ind_pav_green) > 0 )  THEN

!
!--       No vegetation on bare soil, desert or ice:
          IF ( ( lug_palm == ind_luv_b_soil )  .OR.  ( lug_palm == ind_luv_desert ) .OR.           &
               ( lug_palm == ind_luv_ice ) ) THEN

             lai = 0.0_wp
             sai = 0.0_wp

          ENDIF


          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--          Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc, lug_dep,                   &
                                            r_aero_surf, ustar_surf )

                bud_lug(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc, lug_dep,                   &
                                            r_aero_surf, ustar_surf )

                bud_lug(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( -vd_lu * dt_dh ) ) * dh

             ELSE  !< GASES
!
!--             Read spc_name of current species for gas parameter
                IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO
                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9               xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                  c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at k:
                ppm2ugm3 =  (dens/xm_air) * 0.001_wp  ! (mole air)/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value
                diffusivity            = 0.11e-4
!
!--             Overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( ( lug_dep == ilu_water_sea )  .OR.                               &
                                  ( lug_dep == ilu_water_inland ), z0h_surf, ustar_surf,           &
                                    diffusivity, rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac_ij( spc_names(lsp), ts, ustar_surf, solar_rad, cos_zenith, &
                                            rh_surf, lai, sai, nwet, lug_dep,                      &
                                            2, rc_tot, ccomp_tot(lsp), hyp(nzb), diffusivity )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN

                   bud_lug(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / ( r_aero_surf + rb + rc_tot )

                   bud_lug(lsp) = - ( conc_ijk(lsp) - ccomp_tot(lsp) ) *                           &
                                  ( 1.0_wp - EXP( -vd_lu * dt_dh ) )  * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF
!
!--    Windows
       IF ( surf_usm%frac(m,ind_wat_win) > 0 )  THEN
!
!--       No vegetation on windows:
          lai = 0.0_wp
          sai = 0.0_wp

          slinnfac = 1.0_wp
!
!--       Get vd
          DO  lsp = 1, nvar
!
!--          Initialize
             vs = 0.0_wp
             vd_lu = 0.0_wp
             rs = 0.0_wp
             rb = 0.0_wp
             rc_tot = 0.0_wp
             IF ( spc_names(lsp) == 'PM10' )  THEN
                part_type = 1
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc,                            &
                                            lud_dep, r_aero_surf, ustar_surf )

                bud_lud(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( - vd_lu * dt_dh ) ) * dh

             ELSEIF ( spc_names(lsp) == 'PM25' )  THEN
                part_type = 2
!
!--             Sedimentation velocity
                vs = slinnfac * sedimentation_velocity( particle_pars(ind_p_dens, part_type),      &
                                                        particle_pars(ind_p_size, part_type),      &
                                                        particle_pars(ind_p_slip, part_type),      &
                                                        visc)

                CALL drydepo_aero_zhang_vd( vd_lu, rs, vs,                                         &
                                            particle_pars(ind_p_size, part_type),                  &
                                            particle_pars(ind_p_slip, part_type),                  &
                                            nwet, temp_tmp, dens, visc, lud_dep,                   &
                                            r_aero_surf, ustar_surf )

                bud_lud(lsp) = - conc_ijk(lsp) * ( 1.0_wp - EXP( - vd_lu * dt_dh ) ) * dh

             ELSE  !< GASES
!
!--             Read spc_name of current species for gas PARAMETER
                IF ( ANY( pspecnames(:) == spc_names(lsp) ) )  THEN
                   i_pspec = 0
                   DO  pspec = 1, nposp
                      IF ( pspecnames(pspec) == spc_names(lsp) )  THEN
                         i_pspec = pspec
                      END IF
                   ENDDO
                ELSE
!
!--                For now species not deposited
                   CYCLE
                ENDIF
!
!--             Factor used for conversion from ppb to ug/m3 :
!--                 ppb (mole tr)/(mole air)/ppb (kg tr)/(mole tr) (ug tr)/(kg tr) &
!--                 (mole air)/(kg air) (kg air)/(m3 air) (kg air(ug/m3)/ppb/(kg/mole) = / (kg/mole)
!--                 c           1e-9               xm_tracer         1e9       /       xm_air            dens
!--             thus:
!--                  c_in_ppb * xm_tracer * [ dens / xm_air ] = c_in_ugm3
!--             Use density at k:

                ppm2ugm3 =  (dens/xm_air) * 0.001_wp  ! (mole air)/m3
!
!--             Diffusivity for DEPAC relevant gases
!--             Use default value
                diffusivity = 0.11e-4
!
!--             Overwrite with known coefficients of diffusivity from Massman (1998)
                IF ( spc_names(lsp) == 'NO2' ) diffusivity = 0.136e-4
                IF ( spc_names(lsp) == 'NO'  ) diffusivity = 0.199e-4
                IF ( spc_names(lsp) == 'O3'  ) diffusivity = 0.144e-4
                IF ( spc_names(lsp) == 'CO'  ) diffusivity = 0.176e-4
                IF ( spc_names(lsp) == 'SO2' ) diffusivity = 0.112e-4
                IF ( spc_names(lsp) == 'CH4' ) diffusivity = 0.191e-4
                IF ( spc_names(lsp) == 'NH3' ) diffusivity = 0.191e-4
!
!--             Get quasi-laminar boundary layer resistance rb:
                CALL get_rb_cell( ( lud_dep == ilu_water_sea )  .OR.                               &
                                  ( lud_dep == ilu_water_inland ), z0h_surf, ustar_surf,           &
                                    diffusivity, rb )
!
!--             Get rc_tot
                CALL drydepos_gas_depac_ij( spc_names(lsp), ts, ustar_surf, solar_rad, cos_zenith, &
                                            rh_surf, lai, sai, nwet, lud_dep,                      &
                                            2, rc_tot, ccomp_tot(lsp), hyp(nzb), diffusivity )
!
!--             Calculate budget
                IF ( rc_tot <= 0.0 )  THEN

                   bud_lud(lsp) = 0.0_wp

                ELSE

                   vd_lu = 1.0_wp / (r_aero_surf + rb + rc_tot )

                   bud_lud(lsp) = - ( conc_ijk(lsp) - ccomp_tot(lsp) ) *                           &
                                  ( 1.0_wp - EXP( -vd_lu * dt_dh ) )  * dh
                ENDIF

             ENDIF
          ENDDO
       ENDIF


       bud = 0.0_wp
!
!--    Calculate overall budget for surface m and adapt concentration
       DO  lsp = 1, nspec

          bud(lsp) = surf_usm%frac(m,ind_veg_wall)  * bud_luu(lsp) +                             &
                     surf_usm%frac(m,ind_pav_green) * bud_lug(lsp) +                             &
                     surf_usm%frac(m,ind_wat_win)   * bud_lud(lsp)
!
!--       Compute new concentration
          conc_ijk(lsp) = conc_ijk(lsp) + bud(lsp) * inv_dh

          chem_species(lsp)%conc(k,j,i) = MAX( 0.0_wp, conc_ijk(lsp) )

       ENDDO

    ENDIF


 END SUBROUTINE chem_depo_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute total canopy (or surface) resistance Rc for gases.
!> Vector-optimized version. (For further documentation, please see drydepos_gas_depac_ij.
!>
!> DEPAC:
!> Code of the DEPAC routine and corresponding subroutines below from the DEPAC module of the
!> LOTOS-EUROS model (Manders et al., 2017).
!>
!> Original DEPAC routines by RIVM and TNO (2015), for Documentation see van Zanten et al., 2010.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE drydepos_gas_depac( compnam, mvec, t, ust, sai_present, solar_rad, sinphi, rh, lai,    &
                                sai, nwet, lu, iratns, rc_tot, ccomp_tot, p, diffusivity )

    CHARACTER(LEN=*), INTENT(IN) ::  compnam  !< component name: 'HNO3','NO','NO2','O3','SO2','NH3'

    INTEGER                  ::  i              !< grid index in x-direction
    INTEGER(iwp)             ::  icmp           !< component number taken from component name, paramteres matched with include files
    INTEGER(iwp), PARAMETER  ::  icmp_o3   = 1  !< o3
    INTEGER(iwp), PARAMETER  ::  icmp_so2  = 2  !< so2
    INTEGER(iwp), PARAMETER  ::  icmp_no2  = 3  !< no2
    INTEGER(iwp), PARAMETER  ::  icmp_no   = 4  !< no
    INTEGER(iwp), PARAMETER  ::  icmp_nh3  = 5  !< nh3
    INTEGER(iwp), PARAMETER  ::  icmp_co   = 6  !< co
    INTEGER(iwp), PARAMETER  ::  icmp_no3  = 7  !< no3
    INTEGER(iwp), PARAMETER  ::  icmp_hno3 = 8  !< hno3
    INTEGER(iwp), PARAMETER  ::  icmp_n2o5 = 9  !< n2o5
    INTEGER(iwp), PARAMETER  ::  icmp_h2o2 = 10 !< h2o2
    INTEGER(iwp), INTENT(IN) ::  iratns         !< index for NH3/SO2 ratio used for SO2:
                                                !< iratns = 1: low NH3/SO2
                                                !< iratns = 2: high NH3/SO2
                                                !< iratns = 3: very low NH3/SO2
    INTEGER                  ::  j              !< grid index in y-direction
    INTEGER                  ::  m              !< surface index


    INTEGER(iwp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr) ::  mvec    !< pre-calculated upward-facing surface index
    INTEGER(iwp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr) ::  lu      !< land use type, lu = 1,...,nlu
    INTEGER(iwp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr) ::  nwet    !< wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow

    LOGICAL      :: match_lsm      !< flag indicating LSM surface at given j,i
    LOGICAL      :: do_rw_constant

    LOGICAL, DIMENSION(nys:nyn,nxl:nxr)              ::  ready       !< Rc has been set:
                                                                     !< = 1 -> constant Rc
                                                                     !< = 2 -> temperature dependent Rc
    LOGICAL, INTENT(IN), DIMENSION(nys:nyn,nxl:nxr)  ::  sai_present !< vegetation is present for current land use

    REAL(wp), INTENT(OUT) ::  ccomp_tot        !< total compensation point (ug/m3) (= 0 for species that don't have a compensation)
    REAL(wp)              ::  csoil            !< soil compensation point (ug/m3)
    REAL(wp)              ::  cstom            !< stomatal compensation point (ug/m3)
    REAL(wp)              ::  cw               !< external leaf surface compensation point
    REAL(wp)              ::  const_val        !< Input value for rw_constant (ug/m3)
    REAL(wp), PARAMETER   ::  dO3 = 0.13E-4_wp !< diffusion coefficient of ozon (m2/s)
    REAL(wp), INTENT(IN)  ::  diffusivity      !< diffusivity
    REAL(wp), INTENT(IN)  ::  p                !< pressure (Pa)
    REAL(wp)              ::  rsoil_eff        !< effective soil resistance
    REAL(wp), INTENT(IN)  ::  sinphi           !< sin of solar elevation angle
    REAL(wp)              ::  vpd
    REAL(wp)              ::  wet_soil         !< wet soil resistance


    REAL(wp), DIMENSION(nlu_dep) ::  dry_soil1
    REAL(wp), DIMENSION(nys:nyn) ::  dry_soil2
    REAL(wp), DIMENSION(nys:nyn) ::  rinc

    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::               gc_tot     !< total canopy conductance (m/s)
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::               gsoil_eff  !< effective soil conductance (m/s)
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::               gstom      !< stomatal conductance (m/s)
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) ::               gw         !< external leaf conductance (m/s)

    REAL(wp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr)  ::  lai        !< one-sided leaf area index (-)
                                                                     !< hemisphere not possible)
    REAL(wp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr)  ::  rh         !< relative humidity (%)


    REAL(wp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr)  ::  sai        !< surface area index (-) (lai + branches and stems)

    REAL(wp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr)  ::  solar_rad  !< solar radiation, dirict+diffuse (W/m2)
    REAL(wp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr)  ::  t          !< temperature (C)
    REAL(wp), INTENT(IN), DIMENSION(nys:nyn,nxl:nxr)  ::  ust        !< friction velocity (m/s)
    REAL(wp), INTENT(OUT), DIMENSION(nys:nyn,nxl:nxr) ::  rc_tot     !< total canopy resistance Rc (s/m)

!
!-- Define component number.
    SELECT CASE ( TRIM( compnam ) )

    CASE ( 'O3', 'o3' )
       icmp = icmp_o3

    CASE ( 'SO2', 'so2' )
       icmp = icmp_so2

    CASE ( 'NO2', 'no2' )
       icmp = icmp_no2

    CASE ( 'NO', 'no' )
       icmp = icmp_no

    CASE ( 'NH3', 'nh3' )
       icmp = icmp_nh3

    CASE ( 'CO', 'co' )
       icmp = icmp_co

    CASE ( 'NO3', 'no3' )
       icmp = icmp_no3

    CASE ( 'HNO3', 'hno3' )
       icmp = icmp_hno3

    CASE ( 'N2O5', 'n2o5' )
       icmp = icmp_n2o5

    CASE ( 'H2O2', 'h2o2' )
       icmp = icmp_h2o2

    CASE DEFAULT
!
!--    Component not part of DEPAC --> not deposited.
       RETURN

    END SELECT

!
!-- Inititalize
    gw        = 0.0_wp
    gstom     = 0.0_wp
    gsoil_eff = 0.0_wp
    gc_tot    = 0.0_wp
    cw        = 0.0_wp
    cstom     = 0.0_wp
    csoil     = 0.0_wp
    rc_tot    = 0.0_wp
    ccomp_tot = 0.0_wp
    ready = .FALSE.
!
!-- Set Rc (i.e. rc_tot) in special cases.
!-- Next statements substitute CALL rc_special( icmp, compnam, lu, t, nwet, rc_tot, ready,
!-- ccomp_tot ).
    SELECT CASE ( TRIM( compnam ) )
       CASE( 'HNO3', 'N2O5', 'NO3', 'H2O2' )
!
!--       No separate resistances for HNO3; just one total canopy resistance.
!--       No snow (nwet == 9)
          rc_tot = 10.0_wp
          ready = .TRUE.

       CASE( 'NO', 'CO' )
          WHERE ( lu == ilu_water_sea  .OR.  lu == ilu_water_inland )      ! water
             rc_tot = 2000.0_wp
             ready = .TRUE.
          ELSE WHERE ( nwet == 1 )   ! wet
             rc_tot = 2000.0_wp
             ready = .TRUE.
          END WHERE

       CASE( 'NO2', 'O3', 'SO2', 'NH3' )
!         No Snow

       CASE DEFAULT
          message_string = 'component "'// TRIM( compnam ) // '" not supported'
          CALL message( 'drydepos_gas_depac', 'CHM0020', 1, 2, 0, 6, 0 )
    END SELECT
!
!--    External conductance.
!       CALL rc_gw( compnam, iratns, t, rh, nwet, sai_present, sai,gw )

    do_rw_constant = .FALSE.

    IF ( TRIM( compnam ) == 'NO2' )  THEN
       const_val = 2000.0
       do_rw_constant = .TRUE.
    ENDIF

    IF ( TRIM( compnam ) == 'NO'  .OR.  TRIM( compnam ) == 'CO' ) THEN
       const_val = -9999.0
       do_rw_constant = .TRUE.
    ENDIF

    IF ( TRIM( compnam ) == 'O3' )  THEN
       const_val = 2500.0
       do_rw_constant = .TRUE.
    ENDIF

    IF ( do_rw_constant )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             m = mvec(j,i)

             match_lsm = ( m /= 0 )

             IF ( match_lsm  .AND.  .NOT. ready(j,i) )  THEN
                CALL rw_constant( const_val, sai_present(j,i), gw(j,i) )
             ENDIF
          ENDDO
       ENDDO
    ENDIF
!
!-- SO2 and NH3 are not part of chemical phsatp setup.
!-- The following IF block is not tested yet.
    IF ( TRIM( compnam ) == 'SO2' ) THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             m = mvec(j,i)

             match_lsm = ( m /= 0 )

             IF ( match_lsm .AND. .NOT. ready(j,i) )  THEN
                CALL rw_so2( t(j,i), nwet(j,i), rh(j,i), iratns, sai_present(j,i), gw(j,i) )
             ENDIF
          ENDDO
       ENDDO
    ELSEIF ( TRIM(compnam) == 'NH3' )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             m = mvec(j,i)

             match_lsm = ( m /= 0 )

             IF ( match_lsm  .AND.  .NOT. ready(j,i) )  THEN
                CALL rw_nh3_sutton( t(j,i), rh(j,i), sai_present(j,i), gw(j,i) )
                gw(j,i) = sai(j,i) * gw(j,i)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
!
!--    Stomatal conductance.
!       CALL rc_gstom( icmp, compnam, lu, lai_present, lai, solar_rad, sinphi, t, rh, diffusivity,  &
!                      gstom, p )

    SELECT CASE ( TRIM( compnam ) )

       CASE( 'NO', 'CO' )
!
!--       For no stomatal uptake is neglected:
          gstom = 0.0_wp

       CASE( 'NO2', 'O3', 'SO2', 'NH3' )
!
!--       If vegetation present:
          IF ( ANY(sai_present) )  THEN

             IF ( ANY(solar_rad > 0.0_wp) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      m = mvec(j,i)

                      match_lsm = ( m /= 0 )

                      IF ( match_lsm  .AND.  .NOT. ready(j,i) )  THEN
                         CALL rc_get_vpd( t(j,i), rh(j,i), vpd )
                         CALL rc_gstom_emb( lu(j,i), solar_rad(j,i), t(j,i), vpd, sai_present(j,i),&
                                            lai(j,i), sinphi, gstom(j,i), p )
                      ENDIF
                   ENDDO
                ENDDO
                gstom = gstom * diffusivity / dO3       !< Gstom of Emberson is derived for ozone
             ELSE
                gstom = 0.0_wp
             ENDIF
          ELSE
!
!--          No vegetation; zero conductance (infinite resistance):
             gstom = 0.0_wp
          ENDIF

       CASE DEFAULT
          message_string = 'component "'// TRIM( compnam ) // '" not supported'
          CALL message( 'drydeo_gas_depac', 'CHM0020', 1, 2, 0, 6, 0 )

    END SELECT
!
!-- Effective soil conductance
!    CALL rc_gsoil_eff( icmp, lu, sai, ust, nwet, t, gsoil_eff )
!
!-- rc_gsoil_eff contains a two-fold indirect addressing. This does not vectorize.
!-- Introducing temporaray arrays and precompute values avoids indirect addressing.
!-- The following block vectorizes completely-
    DO  i = nxl, nxr
       dry_soil1(:) = rsoil( :, icmp )

       DO  j = nys, nyn
!
!--       Compute in canopy (in crop) resistance.
          CALL rc_rinc( lu(j,i), sai(j,i), ust(j,i), rinc(j) )
!
!--       Check for missing deposition path.
          IF ( missing( rinc(j) ) )   THEN
             gsoil_eff(j,i) = 0.0_wp
             ready(j,i) = .TRUE.
          ENDIF
       ENDDO


       DO  j = nys, nyn
!
!--       Non-frozen soil - dry.
          IF ( missing( dry_soil1( lu(j,i)) ) )  THEN
             dry_soil2( j) = -9999.0_wp
          ELSE
             dry_soil2( j) = dry_soil1( lu(j,i) ) + rinc(j)
          ENDIF
       ENDDO

       IF ( missing( rsoil_wet( icmp ) ) )  THEN
          wet_soil = -9999.0_wp
       ELSE
          wet_soil = rsoil_wet( icmp ) + rinc(j)
       ENDIF

       DO  j = nys, nyn
          m = mvec(j,i)

          match_lsm = ( m /= 0 )

          IF ( match_lsm  .AND.  .NOT. ready(j,i) )  THEN
!
!--          Frozen soil (temperature below 0).
             IF ( t(j,i) < 0.0_wp )  THEN
                IF ( missing( rsoil_frozen( icmp ) ) )  THEN
                   rsoil_eff = -9999.0_wp
                ELSE
                   rsoil_eff = rsoil_frozen( icmp ) + rinc(j)
                ENDIF
             ELSE
!
!--             Non-frozen soil - wet.
                IF ( nwet(j,i) == 0 )  THEN

                   rsoil_eff = dry_soil2( j )
                ELSE
                   rsoil_eff = wet_soil
                ENDIF
             ENDIF
!
!--          Compute conductance.
             IF ( rsoil_eff > 0.0_wp )  THEN
                gsoil_eff(j,i) = 1.0_wp / rsoil_eff
             ELSE
                gsoil_eff(j,i) = 0.0_wp
             ENDIF
          ENDIF
       ENDDO
    ENDDO

!
!-- Total canopy conductance (gc_tot) and resistance Rc (rc_tot).
    DO  i = nxl, nxr
       DO  j = nys, nyn
          m = mvec(j,i)

          match_lsm = ( m /= 0 )

          IF ( match_lsm  .AND.  .NOT. ready(j,i) )  THEN
              CALL rc_rctot( gstom(j,i), gsoil_eff(j,i), gw(j,i), gc_tot(j,i), rc_tot(j,i) )
          ENDIF
       ENDDO
    ENDDO

 END SUBROUTINE drydepos_gas_depac


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute total canopy (or surface) resistance Rc for gases.
!> Call for each i,j (cache-optimized version).
!>
!> DEPAC:
!> Code of the DEPAC routine and corresponding subroutines below from the DEPAC module of the
!> LOTOS-EUROS model (Manders et al., 2017).
!>
!> Original DEPAC routines by RIVM and TNO (2015), for Documentation see van Zanten et al., 2010.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE drydepos_gas_depac_ij( compnam, t, ust, solar_rad, sinphi, rh, lai, sai,               &
                                   nwet, lu, iratns, rc_tot, ccomp_tot, p, diffusivity )
!
!--   Some of depac arguments are OPTIONAL:
!--    A. compute Rc_tot without compensation points (ccomp_tot will be zero):
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot,&
!--                    ccomp_tot, [smi])
!--    B. compute Rc_tot with compensation points (used for LOTOS-EUROS):
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot,&
!--                    ccomp_tot, [smi], c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water)
!--
!--    C. compute effective Rc based on compensation points (used for OPS):
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot,&
!--                    ccomp_tot, [smi], c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, ra,  &
!--                    rb, rc_eff)
!--    X1. Extra (OPTIONAL) output variables:
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot,&
!--                    ccomp_tot, [smi], c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, ra,  &
!--                    rb, rc_eff, gw_out, gstom_out, gsoil_eff_out, cw_out, cstom_out, csoil_out,    &
!--                    lai_out, sai_out)
!--    X2. Extra (OPTIONAL) needed for stomatal ozone flux calculation (only sunlit leaves):
!--        CALL depac (compnam, day_of_year, lat, t, ust, glrad, sinphi, rh, nwet, lu, iratns, rc_tot,&
!--                    ccomp_tot, [smi], c_ave_prev_nh3, c_ave_prev_so2, catm, gamma_soil_water, ra,  &
!--                    rb, rc_eff, gw_out, gstom_out, gsoil_eff_out, cw_out, cstom_out, csoil_out,    &
!--                    lai_out, sai_out, calc_stom_o3flux, frac_sto_o3_lu, fac_surface_area_2_PLA)


    CHARACTER(LEN=*), INTENT(IN) ::  compnam         !< component name
                                                     !< 'HNO3','NO','NO2','O3','SO2','NH3'
    INTEGER(iwp), INTENT(IN) ::  iratns              !< index for NH3/SO2 ratio used for SO2:
                                                     !< iratns = 1: low NH3/SO2
                                                     !< iratns = 2: high NH3/SO2
                                                     !< iratns = 3: very low NH3/SO2
    INTEGER(iwp), INTENT(IN) ::  lu                  !< land use type, lu = 1,...,nlu
    INTEGER(iwp), INTENT(IN) ::  nwet                !< wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow

    REAL(wp), INTENT(IN) ::  diffusivity             !< diffusivity
    REAL(wp), INTENT(IN) ::  lai                     !< one-sidedleaf area index (-)
    REAL(wp), INTENT(IN) ::  p                       !< pressure (Pa)
    REAL(wp), INTENT(IN) ::  rh                      !< relative humidity (%)
    REAL(wp), INTENT(IN) ::  sai                     !< surface area index (-) (lai + branches and
                                                     !< stems)
    REAL(wp), INTENT(IN) ::  sinphi                  !< sin of solar elevation angle
    REAL(wp), INTENT(IN) ::  solar_rad               !< solar radiation, dirict+diffuse (W/m2)
    REAL(wp), INTENT(IN) ::  t                       !< temperature (C)
    REAL(wp), INTENT(IN) ::  ust                     !< friction velocity (m/s)

    REAL(wp), INTENT(OUT) ::  ccomp_tot              !< total compensation point (ug/m3)
                                                     !< [= 0 for species that don't have a compensation
    REAL(wp), INTENT(OUT) ::  rc_tot                 !< total canopy resistance Rc (s/m)

!-- Local variables:
!
!-- Component number taken from component name, paramteres matched with include files
    INTEGER(iwp) ::  icmp
!
!-- Component numbers:
    INTEGER(iwp), PARAMETER ::  icmp_o3   = 1
    INTEGER(iwp), PARAMETER ::  icmp_so2  = 2
    INTEGER(iwp), PARAMETER ::  icmp_no2  = 3
    INTEGER(iwp), PARAMETER ::  icmp_no   = 4
    INTEGER(iwp), PARAMETER ::  icmp_nh3  = 5
    INTEGER(iwp), PARAMETER ::  icmp_co   = 6
    INTEGER(iwp), PARAMETER ::  icmp_no3  = 7
    INTEGER(iwp), PARAMETER ::  icmp_hno3 = 8
    INTEGER(iwp), PARAMETER ::  icmp_n2o5 = 9
    INTEGER(iwp), PARAMETER ::  icmp_h2o2 = 10

    LOGICAL ::  ready                                !< Rc has been set:
                                                     !< = 1 -> constant Rc
                                                     !< = 2 -> temperature dependent Rc
!
!-- Vegetation indicators:
    LOGICAL ::  lai_present                          !< leaves are present for current land use type
    LOGICAL ::  sai_present                          !< vegetation is present for current land use
                                                     !< type

    REAL(wp) ::  csoil                               !< soil compensation point (ug/m3)
    REAL(wp) ::  cstom                               !< stomatal compensation point (ug/m3)
    REAL(wp) ::  cw                                  !< external leaf surface compensation point
                                                     !< (ug/m3)
    REAL(wp) ::  gc_tot                              !< total canopy conductance (m/s)
    REAL(wp) ::  gsoil_eff                           !< effective soil conductance (m/s)
    REAL(wp) ::  gstom                               !< stomatal conductance (m/s)
    REAL(wp) ::  gw                                  !< external leaf conductance (m/s)
!    REAL(wp) ::  laimax                              !< maximum leaf area index (-)
!
!-- Define component number
    SELECT CASE ( TRIM( compnam ) )

    CASE ( 'O3', 'o3' )
       icmp = icmp_o3

    CASE ( 'SO2', 'so2' )
       icmp = icmp_so2

    CASE ( 'NO2', 'no2' )
       icmp = icmp_no2

    CASE ( 'NO', 'no' )
       icmp = icmp_no

    CASE ( 'NH3', 'nh3' )
       icmp = icmp_nh3

    CASE ( 'CO', 'co' )
       icmp = icmp_co

    CASE ( 'NO3', 'no3' )
       icmp = icmp_no3

    CASE ( 'HNO3', 'hno3' )
       icmp = icmp_hno3

    CASE ( 'N2O5', 'n2o5' )
       icmp = icmp_n2o5

    CASE ( 'H2O2', 'h2o2' )
       icmp = icmp_h2o2

    CASE default
!
!--    Component not part of DEPAC --> not deposited
       RETURN

    END SELECT

!
!-- Inititalize
    gw        = 0.0_wp
    gstom     = 0.0_wp
    gsoil_eff = 0.0_wp
    gc_tot    = 0.0_wp
    cw        = 0.0_wp
    cstom     = 0.0_wp
    csoil     = 0.0_wp
!
!-- Check whether vegetation is present:
    lai_present = ( lai > 0.0 )
    sai_present = ( sai > 0.0 )
!
!-- Set Rc (i.e. rc_tot) in special cases:
    CALL rc_special( icmp, compnam, lu, t, nwet, rc_tot, ready, ccomp_tot )
!
!-- If Rc is not set:
    IF ( .NOT. ready ) then
!
!--    External conductance:
       CALL rc_gw( compnam, iratns, t, rh, nwet, sai_present, sai,gw )
!
!--    Stomatal conductance:
       CALL rc_gstom( icmp, compnam, lu, lai_present, lai, solar_rad, sinphi, t, rh, diffusivity,  &
                      gstom, p )
!
!--    Effective soil conductance:
       CALL rc_gsoil_eff( icmp, lu, sai, ust, nwet, t, gsoil_eff )
!
!--    Total canopy conductance (gc_tot) and resistance Rc (rc_tot):
       CALL rc_rctot( gstom, gsoil_eff, gw, gc_tot, rc_tot )

    ENDIF

 END SUBROUTINE drydepos_gas_depac_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute total canopy resistance in special cases
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rc_special( icmp, compnam, lu, t, nwet, rc_tot, ready, ccomp_tot )


    CHARACTER(LEN=*), INTENT(IN)  ::  compnam     !< component name

    INTEGER(iwp), INTENT(IN)  ::  icmp            !< component index
    INTEGER(iwp), INTENT(IN)  ::  lu              !< land use type, lu = 1,...,nlu
    INTEGER(iwp), INTENT(IN)  ::  nwet            !< wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow

    LOGICAL, INTENT(OUT) ::  ready               !< Rc has been set
                                                 !< = 1 -> constant Rc

    REAL(wp), INTENT(IN)  ::  t                   !< temperature (C)

    REAL(wp), INTENT(OUT) ::  ccomp_tot          !< total compensation point (ug/m3)
    REAL(wp), INTENT(OUT) ::  rc_tot             !< total canopy resistance Rc (s/m)

!
!-- Next line is to avoid compiler warning about unused variable
    IF ( icmp == 0 )  CONTINUE
!
!-- rc_tot is not yet set:
    ready = .FALSE.
!
!-- Default compensation point in special CASEs = 0:
    ccomp_tot = 0.0_wp

    SELECT CASE( TRIM( compnam ) )
    CASE( 'HNO3', 'N2O5', 'NO3', 'H2O2' )
!
!--    No separate resistances for HNO3; just one total canopy resistance:
       IF ( t < -5.0_wp .AND. nwet == 9 )  THEN
!
!--       T < 5 C and snow:
          rc_tot = 50.0_wp
       ELSE
!
!--       All other circumstances:
          rc_tot = 10.0_wp
       ENDIF
       ready = .TRUE.

    CASE( 'NO', 'CO' )
       IF ( lu == ilu_water_sea .OR. lu == ilu_water_inland )  THEN       ! water
          rc_tot = 2000.0_wp
          ready = .TRUE.
       ELSEIF ( nwet == 1 )  THEN  !< wet
          rc_tot = 2000.0_wp
          ready = .TRUE.
       ENDIF
    CASE( 'NO2', 'O3', 'SO2', 'NH3' )
!
!--    snow surface:
       IF ( nwet == 9 )  THEN
!
!--       To be activated when snow is implemented
          !CALL rc_snow(ipar_snow(icmp),t,rc_tot)
          ready = .TRUE.
       ENDIF
    CASE default
       message_string = 'component "'// TRIM( compnam ) // '" not supported'
       CALL message( 'rc_special', 'CHM0020', 1, 2, 0, 6, 0 )
    END SELECT

 END SUBROUTINE rc_special


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute external conductance
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rc_gw( compnam, iratns, t, rh, nwet, sai_present, sai, gw )

!
!-- Input/output variables:
    CHARACTER(LEN=*), INTENT(IN) ::  compnam      !< component name
    INTEGER(iwp), INTENT(IN) ::  iratns           !< index for NH3/SO2 ratio;
                                                  !< iratns = 1: low NH3/SO2
                                                  !< iratns = 2: high NH3/SO2
                                                  !< iratns = 3: very low NH3/SO2

    INTEGER(iwp), INTENT(IN) ::  nwet             !< wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow

    LOGICAL, INTENT(IN) ::  sai_present

    REAL(wp), INTENT(IN) ::  rh                   !< relative humidity (%)
    REAL(wp), INTENT(IN) ::  sai                  !< one-sided leaf area index (-)
    REAL(wp), INTENT(IN) ::  t                    !< temperature (C)

    REAL(wp), INTENT(OUT) ::  gw                  !< external leaf conductance (m/s)

    SELECT CASE( TRIM( compnam ) )

    CASE( 'NO2' )
       CALL rw_constant( 2000.0_wp, sai_present, gw )

    CASE( 'NO', 'CO' )
       CALL rw_constant( -9999.0_wp, sai_present, gw )   !< see Erisman et al, 1994 section 3.2.3

    CASE( 'O3' )
       CALL rw_constant( 2500.0_wp, sai_present, gw )

    CASE( 'SO2' )
       CALL rw_so2( t, nwet, rh, iratns, sai_present, gw )

    CASE( 'NH3' )
       CALL rw_nh3_sutton( t, rh, sai_present, gw )
!
!--    Conversion from leaf resistance to canopy resistance by multiplying with sai:
       gw = sai * gw

    CASE default
       message_string = 'component "'// TRIM( compnam ) // '" not supported'
       CALL message( 'rc_gw', 'CHM0020', 1, 2, 0, 6, 0 )
    END SELECT

 END SUBROUTINE rc_gw


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute external leaf conductance for SO2
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rw_so2( t, nwet, rh, iratns, sai_present, gw )

!
!-- Input/output variables:
    INTEGER(iwp), INTENT(IN) ::  iratns      !< index for NH3/SO2 ratio:
                                             !< iratns = 1: low NH3/SO2
                                             !< iratns = 2: high NH3/SO2
                                             !< iratns = 3: very low NH3/SO2
    INTEGER(iwp), INTENT(IN) ::  nwet        !< wetness indicator; nwet=0 -> dry; nwet=1 -> wet; nwet=9 -> snow

    LOGICAL, INTENT(IN) ::  sai_present

    REAL(wp), INTENT(IN) ::  rh              !< relative humidity (%)
    REAL(wp), INTENT(IN) ::  t               !< temperature (C)

    REAL(wp), INTENT(OUT) ::  gw             !< external leaf conductance (m/s)
!
!-- Local variables:
    REAL(wp) ::  rw                          !< external leaf resistance (s/m)
!
!-- Check if vegetation present:
    IF ( sai_present )  THEN

       IF ( nwet == 0 )  THEN
!
!--   ------------------------
!--         dry surface
!--   ------------------------
!--         T > -1 C
          IF ( t > -1.0_wp )  THEN
             IF ( rh < 81.3_wp )  THEN
                rw = 25000.0_wp * EXP( -0.0693_wp * rh )
             ELSE
                rw = 0.58e12 * EXP( -0.278_wp * rh ) + 10.0_wp
             ENDIF
          ELSE
             ! -5 C < T <= -1 C
             IF ( t > -5.0_wp )  THEN
                rw = 200.0_wp
             ELSE
                ! T <= -5 C
                rw = 500.0_wp
             ENDIF
          ENDIF
       ELSE
!
!--   ------------------------
!--         wet surface
!--   ------------------------
          rw = 10.0_wp !see Table 5, Erisman et al, 1994 Atm. Environment, 0 is impl. as 10
       ENDIF
!
!--    Very low NH3/SO2 ratio:
       IF ( iratns == iratns_very_low ) rw = rw + 50.0_wp
!
!--      Conductance:
       gw = 1.0_wp / rw
    ELSE
!
!--      No vegetation:
       gw = 0.0_wp
    ENDIF

 END SUBROUTINE rw_so2


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute external leaf conductance for NH3, following Sutton & Fowler, 1993
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rw_nh3_sutton( tsurf, rh,sai_present, gw )

!
!-- Input/output variables:
    LOGICAL, INTENT(IN) ::  sai_present

    REAL(wp), INTENT(IN) ::  rh             !< relative humidity (%)
    REAL(wp), INTENT(IN) ::  tsurf          !< surface temperature (C)

    REAL(wp), INTENT(OUT) ::  gw            !< external leaf conductance (m/s)
!
!-- Local variables:
    REAL(wp) ::  rw                         !< external leaf resistance (s/m)
    REAL(wp) ::  sai_grass_haarweg          !< surface area index at experimental site Haarweg
!
!-- Fix sai_grass at value valid for Haarweg data for which gamma_w parametrization is derived
    sai_grass_haarweg = 3.5_wp
!
!-- Calculation rw:
!--                    100 - rh
!--    rw = 2.0 * exp(----------)
!--                      12

    IF ( sai_present )  THEN
!
!--    External resistance according to Sutton & Fowler, 1993
       rw = 2.0_wp * EXP( ( 100.0_wp - rh ) / 12.0_wp )
       rw = sai_grass_haarweg * rw
!
!--    Frozen soil (from Depac v1):
       IF ( tsurf < 0.0_wp ) rw = 200.0_wp
!
!--    Conductance:
       gw = 1.0_wp / rw
    ELSE
       ! no vegetation:
       gw = 0.0_wp
    ENDIF

 END SUBROUTINE rw_nh3_sutton


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute external leaf conductance
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rw_constant( rw_val, sai_present, gw )

!
!-- Input/output variables:
    LOGICAL, INTENT(IN) ::  sai_present

    REAL(wp), INTENT(IN) ::  rw_val       !< constant value of Rw

    REAL(wp), INTENT(OUT) ::  gw          !< wernal leaf conductance (m/s)
!
!-- Compute conductance:
    IF ( sai_present  .AND.  .NOT. missing(rw_val) )  THEN
       gw = 1.0_wp / rw_val
    ELSE
       gw = 0.0_wp
    ENDIF

 END SUBROUTINE rw_constant


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute stomatal conductance
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rc_gstom( icmp, compnam, lu, lai_present, lai, solar_rad, sinphi, t, rh, diffusivity,  &
                      gstom, p )

!
!-- input/output variables:
    CHARACTER(LEN=*), INTENT(IN) ::  compnam       !< component name

    INTEGER(iwp), INTENT(IN) ::  icmp              !< component index
    INTEGER(iwp), INTENT(IN) ::  lu                !< land use type , lu = 1,...,nlu

    LOGICAL, INTENT(IN) ::  lai_present

    REAL(wp), INTENT(IN) ::  diffusivity           !< diffusion coefficient of the gas involved
    REAL(wp), INTENT(IN) ::  lai                   !< one-sided leaf area index
    REAL(wp), INTENT(IN) ::  rh                    !< relative humidity (%)
    REAL(wp), INTENT(IN) ::  sinphi                !< sin of solar elevation angle
    REAL(wp), INTENT(IN) ::  solar_rad             !< solar radiation, dirict+diffuse (W/m2)
    REAL(wp), INTENT(IN) ::  t                     !< temperature (C)

    REAL(wp), OPTIONAL,INTENT(IN) :: p             !< pressure (Pa)

    REAL(wp), INTENT(OUT) ::  gstom                !< stomatal conductance (m/s)
!
!-- Local variables
    REAL(wp), PARAMETER ::  dO3 = 0.13e-4          !< diffusion coefficient of ozon (m2/s)

    REAL(wp) ::  vpd                               !< vapour pressure deficit (kPa)

!
!-- Next line is to avoid compiler warning about unused variables
    IF ( icmp == 0 )  CONTINUE

    SELECT CASE( TRIM( compnam ) )

    CASE( 'NO', 'CO' )
!
!--    For no stomatal uptake is neglected:
       gstom = 0.0_wp

    CASE( 'NO2', 'O3', 'SO2', 'NH3' )
!
!--    If vegetation present:
       IF ( lai_present )  THEN

          IF ( solar_rad > 0.0_wp )  THEN
             CALL rc_get_vpd( t, rh, vpd )
             CALL rc_gstom_emb( lu, solar_rad, t, vpd, lai_present, lai, sinphi, gstom, p )
             gstom = gstom * diffusivity / dO3       !< Gstom of Emberson is derived for ozone
          ELSE
             gstom = 0.0_wp
          ENDIF
       ELSE
!
!--       No vegetation; zero conductance (infinite resistance):
          gstom = 0.0_wp
       ENDIF

    CASE default
       message_string = 'component "'// TRIM( compnam ) // '" not supported'
       CALL message( 'rc_gstom', 'CHM0020', 1, 2, 0, 6, 0 )
    END SELECT

 END SUBROUTINE rc_gstom


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to compute stomatal conductance according to Emberson
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rc_gstom_emb( lu, solar_rad, T, vpd, lai_present, lai, sinp, Gsto, p )
!
!>  History
!>   Original code from Lotos-Euros, TNO, M. Schaap
!>   2009-08, M.C. van Zanten, Rivm
!>     Updated and extended.
!>   2009-09, Arjo Segers, TNO
!>     Limitted temperature influence to range to avoid floating point exceptions.

!> Method

!>   Code based on Emberson et al, 2000, Env. Poll., 403-413
!>   Notation conform Unified EMEP Model Description Part 1, ch 8
!
!>   In the calculation of f_light the modification of L. Zhang 2001, AE to the PARshade and PARsun
!>   parametrizations of Norman 1982 are applied
!>   f_phen and f_SWP are set to 1
!
!>   Land use types DEPAC versus Emberson (Table 5.1, EMEP model description)
!>   DEPAC                     Emberson
!>     1 = grass                 GR = grassland
!>     2 = arable land           TC = temperate crops ( lai according to RC = rootcrops)
!>     3 = permanent crops       TC = temperate crops ( lai according to RC = rootcrops)
!>     4 = coniferous forest     CF = tempareate/boREAL(wp) coniferous forest
!>     5 = deciduous forest      DF = temperate/boREAL(wp) deciduous forest
!>     6 = water                 W  = water
!>     7 = urban                 U  = urban
!>     8 = other                 GR = grassland
!>     9 = desert                DE = desert
!
!-- Emberson specific declarations
!
!-- Input/output variables:
    INTEGER(iwp), INTENT(IN) ::  lu             !< land use type, lu = 1,...,nlu

    LOGICAL, INTENT(IN) ::  lai_present

    REAL(wp), INTENT(IN) ::  lai                !< one-sided leaf area index
    REAL(wp), INTENT(IN) ::  sinp               !< sin of solar elevation angle
    REAL(wp), INTENT(IN) ::  solar_rad          !< solar radiation, dirict+diffuse (W/m2)
    REAL(wp), INTENT(IN) ::  t                  !< temperature (C)
    REAL(wp), INTENT(IN) ::  vpd                !< vapour pressure deficit (kPa)


    REAL(wp), OPTIONAL, INTENT(IN) ::  p        !< pressure (Pa)

    REAL(wp), INTENT(OUT) ::  gsto              !< stomatal conductance (m/s)
!
!-- Local variables:
    REAL(wp), PARAMETER ::  p_sealevel = 1.01325e05    !< Pa

    REAL(wp) ::  bt
    REAL(wp) ::  f_env
    REAL(wp) ::  f_light
    REAL(wp) ::  f_phen
    REAL(wp) ::  f_temp
    REAL(wp) ::  f_swp
    REAL(wp) ::  f_vpd
    REAL(wp) ::  laishade
    REAL(wp) ::  laisun
    REAL(wp) ::  pardiff
    REAL(wp) ::  pardir
    REAL(wp) ::  parshade
    REAL(wp) ::  parsun
    REAL(wp) ::  pres
    REAL(wp) ::  sinphi

!
!-- Check whether vegetation is present:
    IF ( lai_present )  THEN

       ! Calculation of correction factors for stomatal conductance
       IF ( sinp <= 0.0_wp )  THEN
          sinphi = 0.0001_wp
       ELSE
          sinphi = sinp
       END IF
!
!--    Ratio between actual and sea-level pressure is used to correct for height in the computation
!--    of par; should not exceed sea-level pressure therefore ...
       IF (  PRESENT( p ) )  THEN
          pres = MIN( p, p_sealevel )
       ELSE
          pres = p_sealevel
       ENDIF
!
!--    Direct and diffuse par, Photoactive (=visible) radiation:
       CALL par_dir_diff( solar_rad, sinphi, pres, p_sealevel, pardir, pardiff )
!
!--    Par for shaded leaves (canopy averaged):
       parshade = pardiff * EXP( -0.5 * lai**0.7 ) + 0.07 * pardir * ( 1.1 - 0.1 * lai ) *         &
                  EXP( -sinphi )     !< Norman,1982
       IF ( solar_rad > 200.0_wp  .AND.  lai > 2.5_wp )  THEN
          parshade = pardiff * EXP( -0.5 * lai**0.8 ) + 0.07 * pardir * ( 1.1 - 0.1 * lai ) *      &
                     EXP( -sinphi )  !< Zhang et al., 2001
       END IF
!
!--    Par for sunlit leaves (canopy averaged):
!--    alpha -> mean angle between leaves and the sun is fixed at 60 deg -> i.e. cos alpha = 0.5
       parsun = pardir * 0.5/sinphi + parshade             !< Norman, 1982
       IF ( solar_rad > 200.0_wp  .AND.  lai > 2.5_wp )  THEN
          parsun = pardir**0.8 * 0.5 / sinphi + parshade   !< Zhang et al., 2001
       END IF
!
!--    Leaf area index for sunlit and shaded leaves:
       IF ( sinphi > 0 )  THEN
          laisun = 2 * sinphi * ( 1 - EXP( -0.5 * lai / sinphi ) )
          laishade = lai - laisun
       ELSE
          laisun = 0
          laishade = lai
       END IF

       f_light = ( laisun * ( 1 - EXP( -1.0_wp * alpha(lu) * parsun ) ) +                          &
                   laishade * ( 1 - EXP( -1.0_wp * alpha(lu) * parshade ) ) ) / lai

       f_light = MAX(f_light,f_min(lu))
!
!--    Temperature influence; only non-zero within range [tmin,tmax]:
       IF ( ( tmin(lu) < t )  .AND.  ( t < tmax(lu) ) )  THEN
          bt = ( tmax(lu) - topt(lu) ) / ( topt(lu) - tmin(lu) )
          f_temp = ( ( t - tmin(lu) ) / ( topt(lu) - tmin(lu) ) ) *                                &
                   ( ( tmax(lu) - t ) / ( tmax(lu) - topt(lu) ) )**bt
       ELSE
          f_temp = 0.0_wp
       END IF
       f_temp = MAX( f_temp, f_min(lu) )
!
!--    Vapour pressure deficit influence
       f_vpd = MIN( 1.0_wp, ( ( 1.0_wp - f_min(lu) ) * ( vpd_min(lu) - vpd ) /                     &
                              ( vpd_min(lu) - vpd_max(lu) ) + f_min(lu) ) )
       f_vpd = MAX( f_vpd, f_min(lu) )

       f_swp = 1.0_wp
!
!--    Influence of phenology on stom. conductance
!--    Ignored for now in DEPAC since influence of f_phen on lu classes in use is negligible.
!--    When other EMEP classes (e.g. med. broadleaf) are used f_phen might be too important to
!--    ignore.
       f_phen = 1.0_wp
!
!--    Evaluate total stomatal conductance
       f_env = f_temp * f_vpd * f_swp
       f_env = MAX( f_env,f_min(lu) )
       gsto = g_max(lu) * f_light * f_phen * f_env
!
!--    gstom expressed per m2 leafarea;
!--    This is converted with lai to m2 surface.
       gsto = lai * gsto    ! in m/s

    ELSE
       gsto = 0.0_wp
    ENDIF

 END SUBROUTINE rc_gstom_emb


!--------------------------------------------------------------------------------------------------!
 !> par_dir_diff
 !>     Weiss, A., Norman, J.M. (1985) Partitioning solar radiation into direct and diffuse, visible
 !>     and near-infrared components. Agric. Forest Meteorol. 34, 205-213.
 !>     From a SUBROUTINE obtained from Leiming Zhang,
 !>     Meteorological Service of Canada
 !>     Leiming uses solar irradiance. This should be equal to global radiation and
 !>     Willem Asman set it to global radiation (here defined as solar radiation, dirict+diffuse)
 !>
 !>     @todo Check/connect/replace with radiation_model_mod variables
 !-------------------------------------------------------------------------------------------------!
 SUBROUTINE par_dir_diff( solar_rad, sinphi, pres, pres_0, par_dir, par_diff )


    REAL(wp), INTENT(IN) ::  pres            !< actual pressure (to correct for height) (Pa)
    REAL(wp), INTENT(IN) ::  pres_0          !< pressure at sea level (Pa)
    REAL(wp), INTENT(IN) ::  sinphi          !< sine of the solar elevation
    REAL(wp), INTENT(IN) ::  solar_rad       !< solar radiation, dirict+diffuse (W m-2)

    REAL(wp), INTENT(OUT) ::  par_diff       !< par diffuse: visible (photoactive) diffuse radiation
                                             !< (W m-2)
    REAL(wp), INTENT(OUT) ::  par_dir        !< par direct : visible (photoactive) direct beam
                                             !< radiation (W m-2)

    REAL(wp) ::  fv                          !< par direct beam fraction (dimensionless)
    REAL(wp) ::  ratio                       !< ratio measured to potential solar radiation
                                             !< (dimensionless)
    REAL(wp) ::  rdm                         !< potential direct beam near-infrared radiation
                                             !< (W m-2); "potential" means clear-sky
    REAL(wp) ::  rdn                         !< potential diffuse near-infrared radiation (W m-2)
    REAL(wp) ::  rdu                         !< visible (par) direct beam radiation (W m-2)
    REAL(wp) ::  rdv                         !< potential visible (par) diffuse radiation (W m-2)
    REAL(wp) ::  rn                          !< near-infrared radiation (W m-2)
    REAL(wp) ::  rv                          !< visible radiation (W m-2)
    REAL(wp) ::  sv                          !< total visible radiation
    REAL(wp) ::  ww                          !< water absorption in the near infrared for 10 mm of
                                             !< precipitable water

!
!-- Calculate visible (PAR) direct beam radiation
!-- 600 W m-2 represents average amount of par (400-700 nm wavelength) at the top of the atmosphere;
!-- this is roughly 0.45*solar constant (solar constant=1320 Wm-2)
    rdu = 600.0_wp* EXP( -0.185_wp * ( pres / pres_0 ) / sinphi ) * sinphi
!
!-- Calculate potential visible diffuse radiation
    rdv = 0.4_wp * ( 600.0_wp - rdu ) * sinphi
!
!-- Calculate the water absorption in the-near infrared
    ww = 1320 * 10**( -1.195_wp + 0.4459_wp * LOG10( 1.0_wp / sinphi ) - 0.0345_wp *               &
                      ( LOG10( 1.0_wp / sinphi ) )**2 )
!
!-- Calculate potential direct beam near-infrared radiation
    rdm = ( 720.0_wp * EXP( -0.06_wp * ( pres / pres_0) / sinphi ) - ww ) * sinphi  !< 720 = solar
                                                                                    !< constant - 600
!
!-- Calculate potential diffuse near-infrared radiation
    rdn = 0.6_wp * ( 720 - rdm - ww ) * sinphi
!
!-- Compute visible and near-infrared radiation
    rv = MAX( 0.1_wp, rdu + rdv )
    rn = MAX( 0.01_wp, rdm + rdn )
!
!-- Compute ratio between input global radiation (here defined as solar radiation, dirict+diffuse)
!-- and total radiation computed here.
    ratio = MIN( 0.89_wp, solar_rad / ( rv + rn ) )
!
!-- Calculate total visible radiation
    sv = ratio * rv
!
!-- Calculate fraction of par in the direct beam
    fv = MIN( 0.99_wp, ( 0.9_wp - ratio ) / 0.7_wp )              !< help variable
    fv = MAX( 0.01_wp, rdu / rv * ( 1.0_wp - fv**0.6667_wp ) )    !< fraction of par in the direct
                                                                  !< beam
!
!-- Compute direct and diffuse parts of par
    par_dir = fv * sv
    par_diff = sv - par_dir

 END SUBROUTINE par_dir_diff


!--------------------------------------------------------------------------------------------------!
!> rc_get_vpd: get vapour pressure deficit (kPa)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rc_get_vpd( temp, rh, vpd )

!
!-- Input/output variables:
    REAL(wp), INTENT(IN) ::  rh    !< relative humidity (%)
    REAL(wp), INTENT(IN) ::  temp    !< temperature (C)

    REAL(wp), INTENT(OUT) ::  vpd    !< vapour pressure deficit (kPa)
!
!-- Local variables:
    REAL(wp) ::  esat
!
!-- Fit parameters:
    REAL(wp), PARAMETER ::  a1 = 6.113718e-01
    REAL(wp), PARAMETER ::  a2 = 4.43839e-02
    REAL(wp), PARAMETER ::  a3 = 1.39817e-03
    REAL(wp), PARAMETER ::  a4 = 2.9295e-05
    REAL(wp), PARAMETER ::  a5 = 2.16e-07
    REAL(wp), PARAMETER ::  a6 = 3.0e-09
!
!-- esat is saturation vapour pressure (kPa) at temp(C) following Monteith (1973)
    esat = a1 + a2 * temp + a3 * temp**2 + a4 * temp**3 + a5 * temp**4 + a6 * temp**5
    vpd  = esat * ( 1 - rh / 100 )

 END SUBROUTINE rc_get_vpd


!--------------------------------------------------------------------------------------------------!
!> rc_gsoil_eff: compute effective soil conductance
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rc_gsoil_eff( icmp, lu, sai, ust, nwet, t, gsoil_eff )

!
!-- Input/output variables:
    INTEGER(iwp), INTENT(IN) ::  icmp          !< component index
    INTEGER(iwp), INTENT(IN) ::  nwet          !< index for wetness
                                               !< nwet = 0 -> dry; nwet = 1 -> wet; nwet = 9 -> snow
                                               !< N.B. this routine cannot be called with nwet = 9,
                                               !< nwet = 9 should be handled outside this routine.
    INTEGER(iwp), INTENT(IN) ::  lu            !< land use type, lu = 1,..., nlu

    REAL(wp), INTENT(IN) ::  sai               !< surface area index
    REAL(wp), INTENT(IN) ::  t                 !< temperature (C)
    REAL(wp), INTENT(IN) ::  ust               !< friction velocity (m/s)

    REAL(wp), INTENT(OUT) ::  gsoil_eff        !< effective soil conductance (m/s)
!
!-- local variables:
    REAL(wp) ::  rinc                          !< in canopy resistance  (s/m)
    REAL(wp) ::  rsoil_eff                     !< effective soil resistance (s/m)
!
!-- Compute in canopy (in crop) resistance:
    CALL rc_rinc( lu, sai, ust, rinc )
!
!-- Check for missing deposition path:
    IF ( missing(rinc) )  THEN
       rsoil_eff = -9999.0_wp
    ELSE
!
!--    Frozen soil (temperature below 0):
       IF ( t < 0.0_wp )  THEN
          IF ( missing( rsoil_frozen( icmp ) ) )  THEN
             rsoil_eff = -9999.0_wp
          ELSE
             rsoil_eff = rsoil_frozen( icmp ) + rinc
          ENDIF
       ELSE
!
!--       Non-frozen soil; dry:
          IF ( nwet == 0 )  THEN
             IF ( missing( rsoil( lu, icmp ) ) )  THEN
                rsoil_eff = -9999.0_wp
             ELSE
                rsoil_eff = rsoil( lu, icmp ) + rinc
             ENDIF
!
!--       Non-frozen soil; wet:
          ELSEIF ( nwet == 1 )  THEN
             IF ( missing( rsoil_wet( icmp ) ) )  THEN
                rsoil_eff = -9999.0_wp
             ELSE
                rsoil_eff = rsoil_wet( icmp ) + rinc
             ENDIF
          ELSE
             WRITE( message_string, * ) 'illegal nwet = ', nwet
             CALL message( 'rc_gsoil_eff', 'CHM0021', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF
!
!-- Compute conductance:
    IF ( rsoil_eff > 0.0_wp )  THEN
       gsoil_eff = 1.0_wp / rsoil_eff
    ELSE
       gsoil_eff = 0.0_wp
    ENDIF

 END SUBROUTINE rc_gsoil_eff


!--------------------------------------------------------------------------------------------------!
!> rc_rinc: compute in canopy (or in crop) resistance van Pul and Jacobs, 1993, BLM
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rc_rinc( lu, sai, ust, rinc )

!
!-- Input/output variables:
    INTEGER(iwp), INTENT(IN) ::  lu          !< land use class, lu = 1, ..., nlu

    REAL(wp), INTENT(IN) ::  sai             !< surface area index
    REAL(wp), INTENT(IN) ::  ust             !< friction velocity (m/s)

    REAL(wp), INTENT(OUT) ::  rinc           !< in canopy resistance (s/m)
!
!-- b = empirical constant for computation of rinc (in canopy resistance) (= 14 m-1 or -999 if not
!--     applicable)
!-- h = vegetation height (m)                     gra  ara crop con dec wat   urb   oth   des   ice   sav   trf  wai  med semi
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  b = (/ -999, 14, 14, 14, 14, -999, -999, -999, -999, -999, -999, 14, -999,  &
         14, 14 /)
    REAL(wp), DIMENSION(nlu_dep), PARAMETER ::  h = (/ -999, 1,  1,  20, 20, -999, -999, -999, -999, -999, -999, 20, -999,  &
         1 ,  1 /)
!
!-- Compute Rinc only for arable land, perm. crops, forest; otherwise Rinc = 0:
    IF ( b(lu) > 0.0_wp )  THEN
!
!--    Check for u* > 0 (otherwise denominator = 0):
       IF ( ust > 0.0_wp )  THEN
          rinc = b(lu) * h(lu) * sai/ust
       ELSE
          rinc = 1000.0_wp
       ENDIF
    ELSE
       IF ( lu == ilu_grass  .OR.  lu == ilu_other )  THEN
          rinc = -999.0_wp     !< no deposition path for grass, other, and semi-natural
       ELSE
          rinc = 0.0_wp        !< no in-canopy resistance
       ENDIF
    ENDIF

 END SUBROUTINE rc_rinc


!--------------------------------------------------------------------------------------------------!
!> rc_rctot: compute total canopy (or surface) resistance Rc
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE rc_rctot( gstom, gsoil_eff, gw, gc_tot, rc_tot )

!
!-- Input/output variables:
    REAL(wp), INTENT(IN) ::  gsoil_eff     !< effective soil conductance (s/m)
    REAL(wp), INTENT(IN) ::  gstom         !< stomatal conductance (s/m)
    REAL(wp), INTENT(IN) ::  gw            !< external leaf conductance (s/m)

    REAL(wp), INTENT(OUT) ::  gc_tot       !< total canopy conductance (m/s)
    REAL(wp), INTENT(OUT) ::  rc_tot       !< total canopy resistance Rc (s/m)
!
!-- Total conductance:
    gc_tot = gstom + gsoil_eff + gw
!
!-- Total resistance (note: gw can be negative, but no total emission allowed here):
    IF ( gc_tot <= 0.0_wp .OR. gw < 0.0_wp )  THEN
       rc_tot = -9999.0_wp
    ELSE
       rc_tot = 1.0_wp / gc_tot
    ENDIF

 END SUBROUTINE rc_rctot


!--------------------------------------------------------------------------------------------------!
!> missing: check for data that correspond with a missing deposition path this data is represented
!>          by -999
!--------------------------------------------------------------------------------------------------!
 LOGICAL FUNCTION missing( x )

    REAL(wp), INTENT(IN) ::  x

!
!-- Bandwidth for checking (in)equalities of floats
    REAL(wp), PARAMETER :: eps = 1.0e-5

    missing = ( ABS( x + 999.0_wp ) <= eps )

 END FUNCTION missing


 ELEMENTAL FUNCTION sedimentation_velocity( rhopart, partsize, slipcor, visc ) RESULT( vs )

!
!-- in/out
    REAL(wp), INTENT(IN) ::  rhopart                 !< particle density (kg/m3)
    REAL(wp), INTENT(IN) ::  partsize                !< particle size (m)
    REAL(wp), INTENT(IN) ::  slipcor                 !< slip correction factor (m)
    REAL(wp), INTENT(IN) ::  visc                    !< viscosity

    REAL(wp) ::  vs
!
!-- acceleration of gravity:
    REAL(wp), PARAMETER         ::  grav = 9.80665_wp   !< m/s2

!-- sedimentation velocity
    vs = rhopart * ( partsize**2 ) * grav * slipcor / ( 18.0_wp * visc )

 END FUNCTION sedimentation_velocity


 !-------------------------------------------------------------------------------------------------!
 !> Boundary-layer deposition resistance following Zhang (2001)
 !-------------------------------------------------------------------------------------------------!
 SUBROUTINE drydepo_aero_zhang_vd( vd, rs, vs1, partsize, slipcor, nwet, tsurf, dens1, viscos1,    &
                                   luc, ftop_lu, ustar )

!
!-- in/out
    INTEGER(iwp), INTENT(IN) ::  luc         !< DEPAC LU
    INTEGER(iwp), INTENT(IN) ::  nwet        !< 1=rain, 9=snowcover

    REAL(wp), INTENT(IN) ::  dens1           !< air density (kg/m3) in lowest layer
    REAL(wp), INTENT(IN) ::  ftop_lu         !< atmospheric resistnace Ra
    REAL(wp), INTENT(IN) ::  partsize        !< particle diameter (m)
    REAL(wp), INTENT(IN) ::  slipcor         !< slip correction factor
    REAL(wp), INTENT(IN) ::  tsurf           !< surface temperature (K)
    REAL(wp), INTENT(IN) ::  ustar           !< friction velocity u*
    REAL(wp), INTENT(IN) ::  viscos1         !< air viscosity in lowest layer
    REAL(wp), INTENT(IN) ::  vs1             !< sedimentation velocity in lowest layer

    REAL(wp), INTENT(OUT) ::  rs             !< sedimentaion resistance (s/m)
    REAL(wp), INTENT(OUT) ::  vd             !< deposition velocity (m/s)
!
!-- Constants
    REAL(wp), PARAMETER ::  grav     = 9.80665_wp             !< acceleration of gravity (m/s2)
    REAL(wp), PARAMETER ::  epsilon0 = 3.0_wp
    REAL(wp), PARAMETER ::  kb       = 1.38066E-23_wp
    REAL(wp), PARAMETER ::  pi       = 3.141592654_wp      !< pi

    REAL(wp), PARAMETER :: alfa_lu(nlu_dep) = &
         (/ 1.2_wp,  1.2_wp,   1.2_wp,  1.0_wp,  1.0_wp,   100.0_wp, 1.5_wp,  1.2_wp, 50.0_wp, 100.0_wp, &
              1.2_wp, 1.0_wp, 100.0_wp, 1.2_wp, 50.0_wp /)
    REAL(wp), PARAMETER :: gamma_lu(nlu_dep) = &
         (/ 0.54_wp, 0.54_wp,  0.54_wp, 0.56_wp, 0.56_wp,  0.50_wp,  0.56_wp, 0.54_wp, 0.58_wp, 0.50_wp, &
              0.54_wp, 0.56_wp, 0.50_wp, 0.54_wp, 0.54_wp /)
    REAL(wp), PARAMETER ::A_lu(nlu_dep) = &
         (/ 3.0_wp,  3.0_wp,   2.0_wp,  2.0_wp,  7.0_wp, -99.0_wp, 10.0_wp, 3.0_wp, -99.0_wp, -99.0_wp,  &
              3.0_wp, 7.0_wp, -99.0_wp, 2.0_wp, -99.0_wp /)
!
!--   grass  arabl crops conif decid  water  urba  othr  desr  ice   sav  trf   wai  med   sem
!
!-- local
    REAL(wp) ::  diff_part
    REAL(wp) ::  ebrown
    REAL(wp) ::  eimpac
    REAL(wp) ::  einterc
    REAL(wp) ::  kinvisc
    REAL(wp) ::  reffic
    REAL(wp) ::  schmidt
    REAL(wp) ::  stokes
!
!-- Kinetic viscosity & diffusivity
    kinvisc = viscos1 / dens1    !< only needed at surface

    diff_part = kb * tsurf * slipcor / ( 3.0_wp * pi * viscos1 * partsize )
!
!-- Schmidt number
    schmidt = kinvisc / diff_part
!
!-- Calculate collection efficiencie E
    Ebrown = Schmidt**( -gamma_lu(luc) )    !< Brownian diffusion
!
!-- Determine Stokes number, interception efficiency and sticking efficiency R (1 = no rebound)
    IF ( luc == ilu_ice  .OR.  nwet == 9  .OR.  luc == ilu_water_sea  .OR.                         &
         luc == ilu_water_inland )  THEN
       stokes = vs1 * ustar**2 / ( grav * kinvisc )
       einterc = 0.0_wp
       reffic = 1.0_wp
    ELSE IF ( luc == ilu_other  .OR.  luc == ilu_desert )  THEN     !<tundra of desert
       stokes = vs1 * ustar**2 / ( grav * kinvisc )
       einterc = 0.0_wp
       reffic = EXP( - stokes**0.5_wp )
    ELSE
       stokes = vs1 * ustar / ( grav * a_lu(luc) * 1.0E-3_wp )
       einterc = 0.5_wp * ( partsize / ( a_lu(luc) * 1.0E-3_wp ) )**2
       reffic = EXP( - stokes**0.5_wp )
    END IF
!
!-- When surface is wet all particles do not rebound:
    IF ( nwet == 1 )  reffic = 1.0_wp
!
!-- Determine impaction efficiency:
    eimpac = ( stokes / ( alfa_lu(luc) + stokes ) )**2
!
!-- Sedimentation resistance:
    rs = 1.0_wp / ( epsilon0 * MAX( 1.0E-5_wp, ustar ) * ( ebrown + eimpac + einterc ) * reffic )

!-- Deposition velocity according to Seinfeld and Pandis (2006; eq 19.7):
!--
!--              1
!--      vd = ------------------ + vs
!--           Ra + Rs + Ra*Rs*vs
!--
!-- where: Rs = Rb (in Seinfeld and Pandis, 2006)

    vd = 1.0_wp / ( ftop_lu + rs + ftop_lu * rs * vs1) + vs1


 END SUBROUTINE drydepo_aero_zhang_vd


 !-------------------------------------------------------------------------------------------------!
 !> Compute quasi-laminar boundary layer resistance as a function of landuse and tracer
 !> Original EMEP formulation by (Simpson et al, 2003) is used.
 !-------------------------------------------------------------------------------------------------!
 SUBROUTINE get_rb_cell( is_water, z0h, ustar, diffusivity, rb )

!
!-- in/out
    LOGICAL , INTENT(IN) ::  is_water

    REAL(wp), INTENT(IN) ::  diffusivity          !< coefficient of diffusivity
    REAL(wp), INTENT(IN) ::  ustar                !< friction velocity
    REAL(wp), INTENT(IN) ::  z0h                  !< roughness length for heat

    REAL(wp), INTENT(OUT) ::  rb                  !< boundary layer resistance
!
!-- const
    REAL(wp), PARAMETER ::  kappa_stab = 0.35     !< von Karman constant
    REAL(wp), PARAMETER ::  thk = 0.19e-4         !< thermal diffusivity of dry air 20 C
!
!-- Next line is to avoid compiler warning about unused variable
    IF ( is_water  .OR.  ( z0h + kappa_stab ) > 0.0_wp )  CONTINUE
!
!-- Use Simpson et al. (2003)
!-- @TODO: Check rb over water calculation, until then leave commented lines
!--  IF ( is_water )  THEN
!--   org: rb = 1.0_wp / ( kappa_stab * MAX( 0.01_wp, ustar ) ) * LOG( z0h / diffusivity * kappa_stab * MAX( 0.01_wp, ustar ) )
!--        rb = 1.0_wp / ( kappa_stab * MAX( 0.1_wp, ustar ) ) *  LOG( z0h / diffusivity * kappa_stab * MAX( 0.1_wp, ustar ) )
!--  ELSE
    rb = 5.0_wp / MAX( 0.01_wp, ustar ) * ( thk / diffusivity )**0.67_wp
!--  END IF

 END SUBROUTINE get_rb_cell


!--------------------------------------------------------------------------------------------------!
!>  Compute water vapor partial pressure (e_w) given specific humidity Q [(kg water)/(kg air)].
!>
!>  Use that gas law for volume V with temperature T holds for the total mixture as well as the
!>  water part:
!>
!>    R T / V = p_air / n_air = p_water / n_water
!>
!>  thus:
!>
!>    p_water = p_air n_water / n_air
!>
!>  Use:
!>    n_air =   m_air   /        xm_air
!>            [kg air]  /  [(kg air)/(mole air)]
!>  and:
!>    n_water =  m_air * Q  /     xm_water
!>              [kg water]  /  [(kg water)/(mole water)]
!>  thus:
!>    p_water = p_air Q / (xm_water/xm_air)
!--------------------------------------------------------------------------------------------------!

 ELEMENTAL FUNCTION watervaporpartialpressure( q, p ) RESULT( p_w )

!
!-- in/out
    REAL(wp), INTENT(IN) ::  p                      !< air pressure [Pa]
    REAL(wp), INTENT(IN) ::  q                      !< specific humidity [(kg water)/(kg air)]

    REAL(wp) ::  p_w                                !< water vapor partial pressure [Pa]
!
!-- Const
    REAL(wp), PARAMETER  ::  eps = xm_h2o / xm_air  !< mole mass ratio ~ 0.622
!
!-- Partial pressure of water vapor:
    p_w = p * q / eps

 END FUNCTION watervaporpartialpressure


!--------------------------------------------------------------------------------------------------!
!>  Saturation vapor pressure.
!>  From (Stull 1988, eq. 7.5.2d):
!>
!>      e_sat = p0 exp( 17.67 * (T-273.16) / (T-29.66) )     [Pa]
!>
!>  where:
!>      p0 = 611.2 [Pa]   : reference pressure
!>
!>  Arguments:
!>      T  [K]  : air temperature
!>  Result:
!>      e_sat_w  [Pa]  : saturation vapor pressure
!>
!>  References:
!>      Roland B. Stull, 1988
!>      An introduction to boundary layer meteorology.
!--------------------------------------------------------------------------------------------------!
 ELEMENTAL FUNCTION saturationvaporpressure( t ) RESULT( e_sat_w )

!
!-- in/out
    REAL(wp), INTENT(IN) ::  t            !< temperature [K]

    REAL(wp) ::  e_sat_w                  !< saturation vapor pressure  [Pa]
!
!-- Const
    REAL(wp), PARAMETER ::  p0 = 611.2   !< base pressure [Pa]
!
!-- Saturation vapor pressure:
    e_sat_w = p0 * EXP( 17.67_wp * ( t - 273.16_wp ) / ( t - 29.66_wp ) )    !< [Pa]

 END FUNCTION saturationvaporpressure


!--------------------------------------------------------------------------------------------------!
!>  Relative humidity RH [%] is by definition:
!>
!>           e_w             water vapor partial pressure
!>    Rh = -------- * 100
!>         e_sat_w           saturation vapor pressure
!--------------------------------------------------------------------------------------------------!

 ELEMENTAL FUNCTION relativehumidity_from_specifichumidity( q, t, p ) RESULT( rh )

!
!-- in/out
    REAL(wp), INTENT(IN) ::  p    !< air pressure [Pa]
    REAL(wp), INTENT(IN) ::  q    !< specific humidity [(kg water)/(kg air)]
    REAL(wp), INTENT(IN) ::  t    !< temperature [K]

    REAL(wp) ::  rh               !< relative humidity [%]
!
!-- Relative humidity:
    rh = watervaporpartialpressure( q, p ) / saturationvaporpressure( t ) * 100.0_wp

 END FUNCTION relativehumidity_from_specifichumidity


 END MODULE chemistry_model_mod
