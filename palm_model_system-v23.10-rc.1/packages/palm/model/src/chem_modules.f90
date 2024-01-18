!> @file chem_modules.f90
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
! Copyright 2017-2021 Karlsruhe Institute of Technology
! Copyright 2017-2021 Freie Universitaet Berlin
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Farah Kanani-Suehring
! @author Basit Khan
! @author Sabine Banzhaf
! @author Emmanuele Russo
! @author Edward C. Chan
! @author Ilona Jäkel
!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of global PALM-4U chemistry variables.
!--------------------------------------------------------------------------------------------------!
 MODULE chem_modules

    USE kinds

    IMPLICIT NONE

    REAL, PARAMETER ::  xm_air   =   28.964E-3             !< air      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_c     =   12.01115E-3           !< C        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_ca    =   40.07800E-3           !< Ca       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_cd    =  112.41000E-3           !< Cd       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_cl    =   35.45300E-3           !< Cl       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_dummy = 1000.0E-3               !< dummy    molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_f     =   18.99840E-3           !< F        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_h     =    1.00790E-3           !< H        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_k     =   39.09800E-3           !< K        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_mg    =   24.30500E-3           !< Mg       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_n     =   14.00670E-3           !< N        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_na    =   22.98977E-3           !< Na       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_o     =   15.99940E-3           !< O        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_pb    =  207.20000E-3           !< Pb       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_pb210 =  210.00000E-3           !< Pb (210) molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_rn222 =  222.00000E-3           !< Rn (222) molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_s     =   32.06400E-3           !< S        molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_co2   = xm_C + xm_O * 2         !< CO2      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_h2o   = xm_H * 2 + xm_O         !< H2O      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_hno3  = xm_H + xm_N + xm_O * 3  !< HNO3     molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_o3    = xm_O * 3                !< O3       molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_n2o5  = xm_N * 2 + xm_O * 5     !< N2O5     molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_nh4   = xm_N + xm_H * 4         !< NH4      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_no3   = xm_N + xm_O * 3         !< NO3      molecular weight (kg/mol)
    REAL, PARAMETER ::  xm_so4   = xm_S + xm_O * 4         !< SO4      molecular weight (kg/mol)

    CHARACTER(LEN=20) ::  bc_cs_b        = 'dirichlet'         !< namelist parameter: surface
                                                               !< boundary condition for concentration
    CHARACTER(LEN=20) ::  bc_cs_l        = 'undefined'         !< left boundary condition
    CHARACTER(LEN=20) ::  bc_cs_n        = 'undefined'         !< north boundary condition
    CHARACTER(LEN=20) ::  bc_cs_r        = 'undefined'         !< right boundary condition
    CHARACTER(LEN=20) ::  bc_cs_s        = 'undefined'         !< south boundary condition
    CHARACTER(LEN=20) ::  bc_cs_t        = 'initial_gradient'  !< namelist parameter: top boudary
                                                               !< condition for concentration
    CHARACTER(LEN=30) ::  chem_mechanism = 'phstatp'           !< namelist parameter: chemistry
                                                               !< mechanism
    CHARACTER(LEN=80) ::  daytype_mdh    = 'workday'           !< namelist parameter: type of day
                                                               !< - workday, weekend, holiday

    CHARACTER(LEN=80) ::  mode_emis      = 'PARAMETERIZED'     !< namelist parameter: mode of
                                                               !< chemistry emissions -
                                                               !< DEFAULT, EXPERT, PARAMETERIZED
    CHARACTER(LEN=10) ::  photolysis_scheme                    !< 'constant',
                                                               !< 'simple' (Simple parameterisation from MCM,
                                                               !< Saunders et al., 2003, Atmos. Chem. Phys., 3, 161-180
                                                               !< 'fastj'  (Wild et al., 2000, J. Atmos. Chem., 37, 245-282)
                                                               !< STILL NOT IMPLEMENTED
    CHARACTER(LEN=80) ::  time_fac_type  = 'MDH'               !< namelist parameter: type of time treatment in the mode_emis
                                                               !< DEFAULT - HOUR, MDH

    CHARACTER(LEN=11), DIMENSION(99) ::  cs_name             = 'novalue'  !< namelist parameter:
                                                                          !<names of species with given fluxes (see csflux)
    CHARACTER(LEN=11), DIMENSION(99) ::  cs_profile_name     = 'novalue'  !< namelist parameter:
    CHARACTER(LEN=11), DIMENSION(99) ::  data_output_pr_cs   = 'novalue'  !< namelist parameter:
    CHARACTER(LEN=11), DIMENSION(99) ::  matched_spc_name    = 'novalue'    !< namelist with mech
    CHARACTER(LEN=11), DIMENSION(99) ::  surface_csflux_name = 'novalue'  !< namelist parameter:

    INTEGER(iwp) ::  communicator_chem  !< stores the number of the MPI communicator to be used for ghost layer data exchange
                                        !< 1: cyclic, 2: cyclic along x, 3: cyclic along y,  4: non-cyclic
    INTEGER(iwp) ::  cs_pr_count_fl_res                    = 0      !< counter for vertical flux profiles of chem specs
                                                                    !< (resolved-scale).
    INTEGER(iwp) ::  cs_pr_count_fl_sgs                    = 0      !< counter for vertical flux profiles of chem specs (SGS)
    INTEGER(iwp) ::  cs_pr_count_sp                        = 0      !< counter for chemical species profiles
    INTEGER(iwp) ::  cs_vertical_gradient_level_ind(99,10) = -9999  !< grid index values of
                                                                    !< cs_vertical_gradient_level
    INTEGER(iwp) ::  emiss_lod                             = -1     !< namelist parameter: chem emission LOD (same as mode_emis)
                                                                    !< -1 = unassigned, 0 = parameterized, 1 = default,
                                                                    !< 2 = pre-processed
    INTEGER(iwp) ::  ibc_cs_b                                       !< integer flag for bc_cs_b
    INTEGER(iwp) ::  ibc_cs_t                                       !< integer flag for bc_cs_t
    INTEGER(iwp) ::  n_matched_bc                                   !< number of matched biogenic emissions
    INTEGER(iwp) ::  main_street_id                        = 0      !< namelist parameter: lower bound of main street IDs
                                                                    !< (OpenStreetMaps) for PARAMETERIZED mode
    INTEGER(iwp) ::  max_pr_cs                             = 0      !< number of chemistry profiles in output
    INTEGER(iwp) ::  max_street_id                         = 0      !< namelist parameter: upper bound of main street IDs
                                                                    !< (OpenStreetMaps) for PARAMETERIZED mode
    INTEGER(iwp) ::  n_matched_vars                                 !< number of matched emissions variables.
    INTEGER(iwp) ::  side_street_id                        = 0      !< namelist parameter: lower bound of side street IDs
                                                                    !< (OpenStreetMaps) for PARAMETERIZED mode

    INTEGER(iwp), DIMENSION(99) ::  cs_pr_index_fl_res = 0  !< index for chemical species sgs-flux profiles
    INTEGER(iwp), DIMENSION(99) ::  cs_pr_index_fl_sgs = 0  !< index for chemical species resolved-scale flux profiles
    INTEGER(iwp), DIMENSION(99) ::  cs_pr_index_sp     = 0  !< index for chemical species mean profiles
    INTEGER(iwp), DIMENSION(99) ::  hom_index_fl_res   = 0  !< index of the resolved-scale flux profile w.r.the hom array
    INTEGER(iwp), DIMENSION(99) ::  hom_index_fl_sgs   = 0  !< index of the SGS flux profile w.r.t the hom array
    INTEGER(iwp), DIMENSION(99) ::  hom_index_spec     = 0  !< index of the profile with respect to the hom array

    INTEGER(iwp), DIMENSION(:) ::  match_spec_nox(1:2)      !< results of matching the input and model's NOx
    INTEGER(iwp), DIMENSION(:) ::  match_spec_pm(1:3)       !< results of matching the input and model's PMs
    INTEGER(iwp), DIMENSION(:) ::  match_spec_sox(1:2)      !< results of matching the input and model's SOx

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_spec_input      !< index of input chem species for matching routine
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_spec_model      !< index of model chem species for matching routine
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_spec_voc_input  !< index of VOC input components matching model's VOCs
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  match_spec_voc_model  !< index of VOC model species matching the input VOCs comp.

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE  ::  cs_advc_flags_s  !< flags used to degrade order of advection scheme for
                                                                      !< chemical species near walls and lateral boundaries

    LOGICAL ::  bc_dirichlet_cs_l         = .FALSE.  !< flag for indicating a dirichlet condition at
                                                     !< the left boundary
    LOGICAL ::  bc_dirichlet_cs_n         = .FALSE.  !< flag for indicating a dirichlet condition at
                                                     !< the north boundary
    LOGICAL ::  bc_dirichlet_cs_r         = .FALSE.  !< flag for indicating a dirichlet condition at
                                                     !< the right boundary
    LOGICAL ::  bc_dirichlet_cs_s         = .FALSE.  !< flag for indicating a dirichlet condition at
                                                     !< the south boundary
    LOGICAL ::  bc_radiation_cs_l         = .FALSE.  !< flag for indicating a radiation/neumann
                                                     !< condition at the left boundary
    LOGICAL ::  bc_radiation_cs_n         = .FALSE.  !< flag for indicating a radiation/neumann
                                                     !< condition at the north boundary
    LOGICAL ::  bc_radiation_cs_r         = .FALSE.  !< flag for indicating a radiation/neumann
                                                     !< condition at the right boundary
    LOGICAL ::  bc_radiation_cs_s         = .FALSE.  !< flag for indicating a radiation/neumann
                                                     !< condition at the south boundary
    LOGICAL ::  constant_csflux(99)       = .TRUE.   !< internal flag, set to .FALSE. if no
                                                     !< surface_csflux is prescribed
    LOGICAL ::  call_chem_at_all_substeps = .FALSE.  !< namelist parameter: Never set true
                                                     !< except for tests
    LOGICAL ::  chem_gasphase_on          = .TRUE.   !< namelist parameter: flag to switch off
                                                     !< chemical reactions
    LOGICAL ::  deposition_dry            = .FALSE.  !< namelist parameter: flag for activation of
                                                     !< deposition calculation
    LOGICAL ::  emiss_interpolate         = .FALSE.  !< namelist parameter: flag to switch-on interpolation of the emission
                                                     !< between emission timestamps
    LOGICAL ::  emissions_anthropogenic   = .FALSE.  !< namelist parameter: flag for turning on
                                                     !< anthropogenic emissions
    LOGICAL ::  emission_output_required  = .TRUE.   !< internal flag for requiring emission outputs
    LOGICAL ::  emiss_read_legacy_mode    = .FALSE.  !< namelist parameter: flag to read emission
                                                     !< data using legacy mode
    LOGICAL ::  nesting_chem              = .TRUE.   !< apply self-nesting for the chemistry model
    LOGICAL ::  nesting_offline_chem      = .TRUE.   !< apply offline nesting for the chemistry model
    LOGICAL ::  photolysis_shading        = .FALSE.  !< namelist parameter to turn on shading in the potolysis scheme

    REAL(wp) ::  emiss_factor_main(99)             = -9999.0_wp  !< namelist parameter: weighting  factor
                                                                 !< for LOD 0 parameterized emissions
    REAL(wp) ::  emiss_factor_side(99)             = -9999.0_wp  !< namelist parameter: weighting  factor
                                                                 !< for LOD 0 parameterized emissions
    REAL(wp) ::  surface_csflux(99)                = 0.0_wp      !< namelist parameter: surface emission
                                                                 !< flux or LOD 0 parameterized emissions
    REAL(wp) ::  wall_csflux(99,0:5)               = 0.0_wp      !< namelist parameter: ...???

    REAL(wp), DIMENSION(99)     ::  cs_surface     = 0.0_wp        !< namelist parameter: chem species
                                                                   !< concentration at surface
    REAL(wp), DIMENSION(99,100) ::  cs_heights     = 9999999.9_wp  !< namelist parameter: height levels
                                                                   !< for initial chem species concentrations
    REAL(wp), DIMENSION(99,100) ::  cs_profile     = 9999999.9_wp  !< namelist parameter: chem species
                                                                   !< concentration values at cs_heights levels

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  bc_cs_t_val            !< vertical grad of chem spcs, near domain top
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  css                    !< scaling parameter for chem species

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  ebio_distribution  !< emissions values
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  emis_distribution  !< emissions final values (main module output)


    REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE, TARGET ::  sums_ws_l  !< subdomain sum of vertical chemistry flux w'ch'

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  cs_1       !< pointer for swapping of
                                                                     !< timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  cs_2       !< pointer for swapping of
                                                                     !< timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  cs_3       !< pointer for swapping of
                                                                     !< timelevels for respective quantity

    REAL(wp), DIMENSION(:,:,:), POINTER ::  cs     !< pointer: sgs chem spcs
    REAL(wp), DIMENSION(:,:,:), POINTER ::  cs_p   !< pointer: prognostic value of sgs chem spcs
    REAL(wp), DIMENSION(:,:,:), POINTER ::  tcs_m  !< pointer: to tcs array (temp)

!
!-  Define chemical variables within chem_species.
    TYPE species_def

       CHARACTER(LEN=15)                            ::  name         !< name of chemical species
       CHARACTER(LEN=15)                            ::  unit         !< unit (ppm for gases, kg m^-3
                                                                     !< for aerosol tracers)
       REAL(KIND=wp), ALLOCATABLE, DIMENSION(:)     ::  conc_pr_init !< initial profile of chemical species
       REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:)   ::  cssws_av     !< averaged fluxes of trace gases at surface
       REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:)   ::  flux_s_cs    !< 6th-order advective flux at
                                                                     !< south face of grid box of chemical species (='cs')
       REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:)   ::  diss_s_cs    !< artificial numerical dissipation
                                                                     !< flux at south face of grid box of chemical species
       REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:) ::  conc_av      !< averaged concentrations
       REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:) ::  flux_l_cs    !< 6th-order advective flux at
                                                                     !< left face of grid box of chemical species (='cs')
       REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:) ::  diss_l_cs    !< artificial numerical dissipation
                                                                     !< flux at left face of grid box of chemical species
       REAL(KIND=wp), POINTER, DIMENSION(:,:,:)     ::  conc         !< concentrations of trace gases
       REAL(KIND=wp), POINTER, DIMENSION(:,:,:)     ::  conc_p       !< conc at prognostic time level
       REAL(KIND=wp), POINTER, DIMENSION(:,:,:)     ::  tconc_m      !< weighted tendency of conc
                                                                     !< for previous sub-timestep (Runge-Kutta)
    END TYPE species_def

!
!-- Define photolysis frequencies in phot_frequen.
    TYPE photols_def
       CHARACTER(LEN=15)                        :: name  !< name of pgotolysis frequency
       CHARACTER(LEN=15)                        :: unit  !< unit (1/s)
       REAL(KIND=wp), POINTER, DIMENSION(:,:,:) :: freq  !< photolysis frequency
    END TYPE photols_def

    TYPE(species_def), ALLOCATABLE, DIMENSION(:), TARGET ::  chem_species
    TYPE(photols_def), ALLOCATABLE, DIMENSION(:), TARGET ::  phot_frequen

!
!-- NetCDF file default and fill values
    CHARACTER(LEN=*), PARAMETER ::  nc_str_novalue   = 'novalue'

    INTEGER, PARAMETER ::  nc_field_length    =   64
    INTEGER, PARAMETER ::  nc_integer_novalue = -127

!
!-- NetCDF file attributes for emission data files.
    CHARACTER(LEN=*), PARAMETER ::  nc_att_lod       = 'lod'
    CHARACTER(LEN=*), PARAMETER ::  nc_dim_ntime     = 'ntime'
    CHARACTER(LEN=*), PARAMETER ::  nc_dim_nspecies  = 'nspecies'
    CHARACTER(LEN=*), PARAMETER ::  nc_var_timestamp = 'timestamp'
    CHARACTER(LEN=*), PARAMETER ::  nc_var_species   = 'species'

!
!-- Activation indicators for emission modes (sectors).
    LOGICAL ::  emis_biogenic      = .FALSE.   !< biogenic emission mode
    LOGICAL ::  emis_generic       = .FALSE.   !< generic emission mode (LOD 2 only)
    LOGICAL ::  emis_domestic      = .FALSE.   !< household emission mode
    LOGICAL ::  emis_nonstationary = .FALSE.   !< nonstationary (moving) emission mode
    LOGICAL ::  emis_pollen        = .FALSE.   !< pollen emission ode (LOD 0 only)
    LOGICAL ::  emis_pt_source     = .FALSE.   !< emission mode activtation (LOD 0 only)
    LOGICAL ::  emis_traffic       = .FALSE.   !< traffic emission mode

!
!-- Emission mode levels of detail (LODs).
    INTEGER(iwp) ::  emis_biogenic_lod      = 0  !< biogenic emission mode LOD
    INTEGER(iwp) ::  emis_domestic_lod      = 2  !< household emission LOD
    INTEGER(iwp) ::  emis_nonstationary_lod = 2  !< LOD of nonstationary emission sources
    INTEGER(iwp) ::  emis_traffic_lod       = 2  !< traffic emission LOD

!
!-- Options for biogenic emission factors (prefix ebio).
    CHARACTER(LEN=11), DIMENSION(99)    ::  ebio_emis_name    = 'novalue'     !< namelist parameter
    CHARACTER(LEN=11)                   ::  ebio_soilm_method = 'bulk'        !< namelist parameter

    INTEGER(iwp) ::  ebio_max_emis_day = 200  !< namelist parameter: day of max emis of bvoc spcs
    INTEGER(iwp) ::  ebio_rad_method   =   0  !< radiation method used for biogenic emis calc

    INTEGER(KIND=ibp), DIMENSION(19) ::  ebio_pft  = nc_integer_novalue  !< plant functional type (pft) for bio spcs
    INTEGER(KIND=ibp), DIMENSION(99) ::  ebio_tree = nc_integer_novalue  !< namelist parameter: street tree spcs.

    REAL(wp) ::  ebio_dt = 0.0_wp            !< timestep for biogenic emission update
    REAL(wp) ::  ebio_ppfd_factor = 2.02_wp  !< namelist parameter: photosynthetic photoflux density

    REAL(wp), DIMENSION(99,19) ::  ebio_ef_pft  = 9999999.9_wp  !< namelist parameter (bio-compound, pft)
    REAL(wp), DIMENSION(99,99) ::  ebio_ef_tree = 9999999.9_wp  !< namelist parameter (bio-compound, tree_type) in um m-2 s-1

!
!-- Options for domestic emission LOD 0 (prefix emis_domestic).
    INTEGER(iwp) ::  emis_domestic_sampling_k       = 10         !< # vertical layers for temperature sampling

    REAL(wp) ::      emis_domestic_base_temperature = 15.0_wp    !< target base temperature [C]
    REAL(wp) ::      emis_domestic_heating_degree   = 2100.0_wp  !< heating degree count [K/year]
    REAL(wp) ::      emis_domestic_update_interval  = 300.0_wp   !< interval between emission souce updates [s]

!
!-- Doemestic emission building parameters (LOD!=2).
    INTEGER(iwp), PARAMETER ::  emis_domestic_max_bld_types = 6     !< building types (in static file)

    REAL(wp), DIMENSION(emis_domestic_max_bld_types) ::  emis_domestic_compact_factors = -1.0_wp  !< compactness factors [1/m]
    REAL(wp), DIMENSION(emis_domestic_max_bld_types) ::  emis_domestic_energy_demands  = -1.0_wp  !< annual NRG demands [kWh/m2]
!
!-- Doemestic emission species parameters (LOD!=2).
    INTEGER(iwp), PARAMETER ::  emis_domestic_max_species = 20 !< maximum number of species supported

    CHARACTER(LEN=nc_field_length), DIMENSION(emis_domestic_max_species) ::                        &
                                emis_domestic_species_names = ''                   !< species names

    REAL(wp), DIMENSION(emis_domestic_max_species) ::                                              &
                                emis_domestic_species_emission_factors = -1.0_wp   !< species annual emission factor [mol/TJ or kg/TJ]

!
!-- Pollen emission mode : General options (prefix emis_pollen).
    CHARACTER(LEN=nc_field_length) ::  epol_model      = 'zink'     !< pollen model used (Zink only at the moment)
    CHARACTER(LEN=nc_field_length) ::  epol_tke_scheme = 'default'  !< pollen tke schemes (default, dynamic, adhoc)

    LOGICAL ::  epol_ignore_precip = .TRUE.  !< influence on preceipation
    LOGICAL ::  epol_ignore_solar  = .FALSE.  !< influence on preceipation

    INTEGER(iwp) ::  epol_pool_reset_hour = 0  !< time of day (UTC hour) when pollen pool is reset

    REAL(wp) ::  epol_update_interval  = 300.0_wp  !< interval [s] between update of pollen source
    REAL(wp) ::  epol_tke_sgs_fraction =   0.1_wp  !< use for estimate local TKE in resolved scale

!
!-- Pollen emission mode : Species specific options (prefix emis_pollen).
    INTEGER(iwp), PARAMETER ::  emis_pollen_max_species = 10  !< maxium number of pollen species supported

    CHARACTER(LEN=nc_field_length), DIMENSION(emis_pollen_max_species) ::                          &
                                    epol_specs_names = ''     !< names of pollen species (viz KPP entries)

    INTEGER(iwp), DIMENSION(emis_pollen_max_species) ::                                            &
                  epol_tree_specs    = 999     !< tree species type
    INTEGER(iwp), DIMENSION(emis_pollen_max_species) ::                                            &
                  epol_vegetation_specs = 999  !< vegetation species type

    REAL(wp), DIMENSION(emis_pollen_max_species) ::                                                &
              epol_seasonal_factors = -1.0_wp                  !< yield factor due to seasonal conditions

    REAL(wp), DIMENSION(emis_pollen_max_species) ::                                                &
              epol_tuning_factors   = 1.0_wp                  !< tuning factors
!
!-- Options for emission mode for E-PRTR or GRETA point sources (prefix emis_pt_source).
    CHARACTER(LEN=nc_field_length), DIMENSION(99) ::  emis_pt_source_species_names = 'novalue'  !< emission species

    INTEGER(iwp) ::  emis_pt_source_k_spread = 3  !< vertical distribution

    INTEGER(iwp), DIMENSION(199,3) ::  emis_pt_source_locations_ijk = -1  !< grid cell locations

    LOGICAL ::  emis_pt_source_leap_year = .FALSE.  !< leap year correction

    REAL(wp), DIMENSION(19)     ::  emis_pt_source_k_weights     = 1.0_wp  !< vertical dist. weights

    REAL(wp), DIMENSION(199,99) ::  emis_pt_source_annual_values = 0.0_wp  !< annual output

!
!-- Options for wet deposition (prefix chem_wet_deposition)
    INTEGER(iwp) ::  chem_wet_deposition_cloud_level_lower = 10  !< vertical level of lower cloud layer
    INTEGER(iwp) ::  chem_wet_deposition_cloud_level_upper = 15  !< vertical level of lower cloud layer

    LOGICAL ::  chem_wet_deposition                = .FALSE.  !< main switch for wet deposition
    LOGICAL ::  chem_wet_deposition_model_override = .FALSE.  !< override cloud and precipitation

    REAL(wp) ::  chem_wet_deposition_rain_rate = 1.0_wp         !< precipitation override
    REAL(wp) ::  chem_wet_deposition_update_interval = 300.0_wp !< update interval

!
!-- Options for ISORROPIA.
    LOGICAL ::  chem_isorropia = .FALSE.  !< main switch for ISORROPIA

    INTEGER(iwp) ::  chem_isorropia_activity_coefficient_method = 1    !< IACALCI (0,1)
    INTEGER(iwp) ::  chem_isorropia_mass_conservation_mode      = 0    !< NADJI (ISORROPIA II only)
    INTEGER(iwp) ::  chem_isorropia_max_activity_sweep          = 4    !< NSWEEPI
    INTEGER(iwp) ::  chem_isorropia_max_iteration               = 100  !< MAXITI
    INTEGER(iwp) ::  chem_isorropia_mdr_weight_method           = 0    !< WFTYPI (0,1,2)
    INTEGER(iwp) ::  chem_isorropia_root_subdivisions           = 5    !< NDIVI

    REAL(wp) ::  chem_isorropia_activity_tolerance =   0.05_wp  !< EPSACTI
    REAL(wp) ::  chem_isorropia_aerosol_state      =   0.0_wp   !< CNTRL[2] (0,1)
    REAL(wp) ::  chem_isorropia_problem_type       =   0.0_wp   !< CNTRL[1] (0,1)
    REAL(wp) ::  chem_isorropia_solver_tolerance   =   1.0E-6   !< EPSI
    REAL(wp) ::  chem_isorropia_update_interval    = 300.0_wp   !< update interval

 END MODULE chem_modules
