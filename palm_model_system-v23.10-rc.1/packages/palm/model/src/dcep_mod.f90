!> @file dcep_mod.f90
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
! Copyright 2015-2021 Humboldt-Universitaet zu Berlin
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> DCEP model and interfaces.
!
!> DCEP is a multi-layer urban canopy parametrization scheme (Schubert and Grossman 2013) based on
!> the Building Effect Parametrization (BEB). DCEP scheme calculates the incoming and outgoing
!> longwave and shortwave radiation for the surfaces of an urban street canyon, i.e. the roof, wall
!> and ground surfaces. The building morphology of the urban area is used to characterize the urban
!> street canyon using its street and building width as well as its canyon length, and the height
!> distribution of buildings.
!--------------------------------------------------------------------------------------------------!
 MODULE dcep_mod

    USE arrays_3d,                                                                                 &
        ONLY:  dzw,                                                                                &
               pt,                                                                                 &
               rho_air,                                                                            &
               tend,                                                                               &
               u,                                                                                  &
               v,                                                                                  &
               zu,                                                                                 &
               zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  c_p,                                                                                &
               g,                                                                                  &
               kappa,                                                                              &
               pi,                                                                                 &
               sigma_sb

    USE control_parameters,                                                                        &
        ONLY:  dcep,                                                                               &
               debug_output,                                                                       &
               debug_output_timestep,                                                              &
               dt_3d,                                                                              &
               latitude,                                                                           &
               longitude,                                                                          &
               message_string,                                                                     &
               roughness_length,                                                                   &
               time_since_reference_point,                                                         &
               varnamelength

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

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
               nzt,                                                                                &
               topo_top_ind

    USE kinds

    USE palm_date_time_mod,                                                                        &
        ONLY:  date_time_str_len,                                                                  &
               get_date_time,                                                                      &
               hours_per_day,                                                                      &
               seconds_per_hour

    USE pegrid,                                                                                    &
        ONLY: myid

    USE radiation_model_mod,                                                                       &
        ONLY:  calc_zenith,                                                                        &
               cos_zenith,                                                                         &
               dcep_average_radiation,                                                             &
               radiation_scheme,                                                                   &
               rad_lw_in_diff,                                                                     &
               rad_sw_in,                                                                          &
               rad_sw_in_diff,                                                                     &
               rad_sw_in_dir,                                                                      &
               sun_direction,                                                                      &
               sun_dir_lat,                                                                        &
               sun_dir_lon

    USE surface_mod,                                                                               &
        ONLY:  albedo_urb,                                                                         &
               albedop_urb,                                                                        &
               emiss_urb,                                                                          &
               fr_urb,                                                                             &
               surf_lsm,                                                                           &
               t_grad_urb


    IMPLICIT NONE
    
    INTEGER(iwp), PARAMETER ::  n_num_integ = 5 !< constant for numerical integration

    INTEGER(iwp), ALLOCATABLE ::  nz(:,:,:,:)   !< number of wall levels; roof levels are nz+1
    INTEGER(iwp), ALLOCATABLE ::  nz_meso(:,:)  !< number of mesoscale height levels in urban layer

    INTEGER(iwp) ::  iurb_cdrag =  1  !< type of drag coefficient for walls
    INTEGER(iwp) ::  iurb_ls    =  2  !< type of urban length scale parametrization (namelist)
    INTEGER(iwp) ::  ke_ground  = 10  !< number of grid levels in the ground (namelist parameter)
    INTEGER(iwp) ::  ke_roof    = 10  !< number of grid levels in roofs (namelist parameter)
    INTEGER(iwp) ::  ke_uhl     = 13  !< number of urban height levels (namelist parameter)
    INTEGER(iwp) ::  ke_wall    = 10  !< number of grid levels in walls (namelist parameter)
    INTEGER(iwp) ::  n_uclass   =  1  !< number of urban classes in DCEP model (namelist parameter)
    INTEGER(iwp) ::  n_udir     =  4  !< maximum number of street directions (namelist parameter)
    INTEGER(iwp) ::  nz_mesom         !< maximum number of vertical levels in the mesoscale grid
    INTEGER(iwp) ::  raddimfull       !< number of elements in radiation matrix

    LOGICAL ::  limpufl    = .TRUE.    !< treat urban tendencied implicitly? (namelist parameter)
    LOGICAL ::  lrroofs    = .FALSE.   !< flag to use roof radiation for budget (namelist parameter)
    LOGICAL ::  ltintfix   = .FALSE.   !< fixed innermost temperature in urban surfaces (namelist parameter)
    LOGICAL ::  lurbradcor = .TRUE.    !< use correction factor for radiation reduction factor for radiation from canyon sides (namelist)
    LOGICAL ::  lurbvel    = .TRUE.    !< urban modification of wind velocity

    REAL(wp), PARAMETER ::  azero      = 100 * EPSILON( 1.0_wp )  !< small number (almost zero)
    REAL(wp), PARAMETER ::  eps_urb    = 1.0E-13_wp               !< threshold for urban grid cell
    REAL(wp), PARAMETER ::  fill_value = -9999.0_wp               !< value for the _FillValue attribute

    REAL(wp), ALLOCATABLE ::  alb_ground(:)          !< albedo of the ground
    REAL(wp), ALLOCATABLE ::  alb_roof(:)            !< albedo of the roof
    REAL(wp), ALLOCATABLE ::  alb_wall(:)            !< albedo of the wall
    REAL(wp), ALLOCATABLE ::  cdrag_wall(:,:,:,:)    !< drag constant for wall elements
    REAL(wp), ALLOCATABLE ::  cs_ground(:)           !< specific heat of the ground material [J m^3 K^-1]
    REAL(wp), ALLOCATABLE ::  cs_roof(:)             !< specific heat of the roof material [J m^3 K^-1]
    REAL(wp), ALLOCATABLE ::  cs_wall(:)             !< specific heat of the wall material [J m^3 K^-1]
    REAL(wp), ALLOCATABLE ::  dz_ground(:,:)         !< layer thickness of ground
    REAL(wp), ALLOCATABLE ::  dz_meso(:,:,:)         !< size of cells
    REAL(wp), ALLOCATABLE ::  dz_roof(:,:)           !< layer thickness of roof
    REAL(wp), ALLOCATABLE ::  dz_uhl(:)              !< height of wall element, value for ke_urb+1 is for air above  hightest roof level
    REAL(wp), ALLOCATABLE ::  dz_wall(:,:)           !< layer thickness of wall
    REAL(wp), ALLOCATABLE ::  emiss_ground(:)        !< emissivity of ground
    REAL(wp), ALLOCATABLE ::  emiss_roof(:)          !< emissivity of roof
    REAL(wp), ALLOCATABLE ::  emiss_wall(:)          !< emissivity of wall
    REAL(wp), ALLOCATABLE ::  fgos(:,:,:,:,:)        !< groud to sky of double canyon
    REAL(wp), ALLOCATABLE ::  fgow(:,:,:,:,:,:)      !< ground to wall of other canyon
    REAL(wp), ALLOCATABLE ::  fgs(:,:,:,:)           !< from sky to ground
    REAL(wp), ALLOCATABLE ::  fgw(:,:,:,:,:)         !< from ground to wall
    REAL(wp), ALLOCATABLE ::  fr_build(:,:,:,:)
    REAL(wp), ALLOCATABLE ::  fr_roof(:,:,:,:,:)     !< probability that a building has an height equal to z
    REAL(wp), ALLOCATABLE ::  fr_rur(:,:)            !< fraction of non-urban parts in a grid element
    REAL(wp), ALLOCATABLE ::  frs(:,:,:,:,:)         !< from sky to roof
    REAL(wp), ALLOCATABLE ::  fr_street(:,:,:,:)
    REAL(wp), ALLOCATABLE ::  fr_uclass(:,:,:)       !< fraction of urban class
    REAL(wp), ALLOCATABLE ::  fr_udir(:,:,:,:)       !< fraction of street direction
    REAL(wp), ALLOCATABLE ::  fr_wall(:,:,:,:,:)     !< Probability that a building has an height greater or equal to z, corresponds to roof heigths
    REAL(wp), ALLOCATABLE ::  frw(:,:,:,:,:,:)       !< from roof to wall (ke_uhl+1,ke_uhl,n_udir,n_uclass)
    REAL(wp), ALLOCATABLE ::  fsg(:,:,:,:)           !< from sky to ground
    REAL(wp), ALLOCATABLE ::  fsog(:,:,:,:,:)        !< double-canyon sky to ground
    REAL(wp), ALLOCATABLE ::  fsow(:,:,:,:,:,:)      !< double-canyon sky to wall
    REAL(wp), ALLOCATABLE ::  fsr(:,:,:,:,:)         !< from sky to roof
    REAL(wp), ALLOCATABLE ::  fsw(:,:,:,:,:)         !< from sky to wall
    REAL(wp), ALLOCATABLE ::  fwg(:,:,:,:,:)         !< from wall to ground
    REAL(wp), ALLOCATABLE ::  fwog(:,:,:,:,:,:)      !< wall to ground of other canyon
    REAL(wp), ALLOCATABLE ::  fwos(:,:,:,:,:,:)      !< wall to double-canyon sky
    REAL(wp), ALLOCATABLE ::  fwow(:,:,:,:,:,:,:)    !< wall to wall of other canyon
    REAL(wp), ALLOCATABLE ::  fwr(:,:,:,:,:,:)       !< from wall to roof (ke_uhl,ke_uhl+1,n_udir,n_uclass)
    REAL(wp), ALLOCATABLE ::  fws(:,:,:,:,:)         !< from wall to sky
    REAL(wp), ALLOCATABLE ::  fww(:,:,:,:,:,:)       !< from wall to wall (ke_uhl,ke_uhl,n_udir,n_uclass)
    REAL(wp), ALLOCATABLE ::  radcor(:,:,:,:)        !< radiation correction factor
    REAL(wp), ALLOCATABLE ::  rl_ground(:,:,:,:)     !< incoming longwave radiation on grounds
    REAL(wp), ALLOCATABLE ::  rlipiv(:,:,:,:,:)      !< pivot elements for LU decomposed longwave radiation matrix
    REAL(wp), ALLOCATABLE ::  rlmatrix(:,:,:,:,:,:)  !< coefficients for longwave radiation budget
    REAL(wp), ALLOCATABLE ::  rl_roof(:,:,:,:,:)     !< incoming longwave radiation on roofs
    REAL(wp), ALLOCATABLE ::  rl_wall(:,:,:,:,:)     !< incoming longwave radiation on walls
    REAL(wp), ALLOCATABLE ::  rs_ground(:,:,:,:)     !< incoming shortwave radiation on grounds
    REAL(wp), ALLOCATABLE ::  rsipiv(:,:,:,:,:)      !< pivot elements for LU decomposed shortwave radiation matrix
    REAL(wp), ALLOCATABLE ::  rsmatrix(:,:,:,:,:,:)  !< coefficients for longwave radiation budget
    REAL(wp), ALLOCATABLE ::  rs_roof(:,:,:,:,:)     !< incoming shortwave radiation on roofs
    REAL(wp), ALLOCATABLE ::  rs_wall(:,:,:,:,:)     !< incoming shortwave radiation on walls
    REAL(wp), ALLOCATABLE ::  shfl_ground(:,:,:,:)   !< sensible heat flux from grounds
    REAL(wp), ALLOCATABLE ::  shfl_roof(:,:,:,:,:)   !< sensible heat flux from roofs
    REAL(wp), ALLOCATABLE ::  shfl_urb(:,:)          !< total sensible heat flux from urban surfaces
    REAL(wp), ALLOCATABLE ::  shfl_wall(:,:,:,:,:)   !< sensible heat flux from walls
    REAL(wp), ALLOCATABLE ::  strfl_ground(:,:,:,:)  !< storage heat flux from ground
    REAL(wp), ALLOCATABLE ::  strfl_roof(:,:,:,:,:)  !< storage heat flux from roofs
    REAL(wp), ALLOCATABLE ::  strfl_urb(:,:)         !< total ground flux from urban surfaces
    REAL(wp), ALLOCATABLE ::  strfl_wall(:,:,:,:,:)  !< storage heat flux from walls

    REAL(wp), ALLOCATABLE, TARGET ::  angrotx_udir(:,:,:)    !< street directions x component in rotated system
    REAL(wp), ALLOCATABLE, TARGET ::  angroty_udir(:,:,:)    !< street directions y component in rotated system
    REAL(wp), ALLOCATABLE, TARGET ::  l_urb(:,:,:)           !< urban length scale
    REAL(wp), ALLOCATABLE, TARGET ::  tkeb_urb(:,:,:)        !< explicit component for TKE due to buoyancy
    REAL(wp), ALLOCATABLE, TARGET ::  tkes_urb(:,:,:)        !< explicit component for TKE due to shear stress
    REAL(wp), ALLOCATABLE, TARGET ::  tt_urb(:,:,:)          !< explicit component for t
    REAL(wp), ALLOCATABLE, TARGET ::  umfl_ground(:,:,:,:)   !< u momentum flux from grounds
    REAL(wp), ALLOCATABLE, TARGET ::  umfl_roof(:,:,:,:,:)   !< u momentum flux from roofs
    REAL(wp), ALLOCATABLE, TARGET ::  umfl_wall(:,:,:,:,:)   !< u momentum flux from walls
    REAL(wp), ALLOCATABLE, TARGET ::  ut_urb(:,:,:)          !< explicit component for u
    REAL(wp), ALLOCATABLE, TARGET ::  vmfl_ground(:,:,:,:)   !< v momentum flux from grounds
    REAL(wp), ALLOCATABLE, TARGET ::  vmfl_roof(:,:,:,:,:)   !< v momentum flux from roofs
    REAL(wp), ALLOCATABLE, TARGET ::  vmfl_wall(:,:,:,:,:)   !< v momentum flux from walls
    REAL(wp), ALLOCATABLE, TARGET ::  volairhlred_urb(:,:,:) !< volume of air (without building volume) of half levels reduced by 
                                                             !< the volume of air for which shfl_urb and u/vmfl_urb is accounted for
    REAL(wp), ALLOCATABLE, TARGET ::  volairhl_urb(:,:,:)    !< volume of air (without building volume) of half levels
    REAL(wp), ALLOCATABLE, TARGET ::  volair_urb(:,:,:)      !< volume of air (without building volume) of main levels
    REAL(wp), ALLOCATABLE, TARGET ::  vt_urb(:,:,:)          !< explicit component for v
    REAL(wp), ALLOCATABLE, TARGET ::  zeff_urb(:,:,:)        !< effective height over ground

    REAL(wp), ALLOCATABLE ::  tddz_ground(:,:)      !< thermal diffusicity divided by length of grounds
    REAL(wp), ALLOCATABLE ::  tddz_roof(:,:)        !< thermal diffusicity divided by length of roofs
    REAL(wp), ALLOCATABLE ::  tddz_wall(:,:)        !< thermal diffusicity divided by length of walls
    REAL(wp), ALLOCATABLE ::  td_ground (:,:)       !< ground thermal diffusivity [m^2 s^-1]
    REAL(wp), ALLOCATABLE ::  td_roof (:,:)         !< roof thermal diffusivity [m^2 s^-1]
    REAL(wp), ALLOCATABLE ::  td_wall (:,:)         !< wall thermal diffusivity [m^2 s^-1]
    REAL(wp), ALLOCATABLE ::  t_ground(:,:,:,:,:)   !< temperature in each layer of the ground
    REAL(wp), ALLOCATABLE ::  t_g_urb(:,:)          !< effective urban ground temperature
    REAL(wp), ALLOCATABLE ::  tint_ground(:)        !< fixed inner or initial temperature of the ground
    REAL(wp), ALLOCATABLE ::  tint_roof(:)          !< fixed inner or initial temperature of the roof
    REAL(wp), ALLOCATABLE ::  tint_wall(:)          !< fixed inner or initial temperature of the wall
    REAL(wp), ALLOCATABLE ::  t_roof(:,:,:,:,:,:)   !< temperature in each layer of the roof
    REAL(wp), ALLOCATABLE ::  t_wall(:,:,:,:,:,:)   !< temperature in each layer of the wall
    REAL(wp), ALLOCATABLE ::  umfl_urb(:,:)         !< total u momentum flux from urban surfaces
    REAL(wp), ALLOCATABLE ::  vmfl_urb(:,:)         !< total v momentum flux from urban surfaces
    REAL(wp), ALLOCATABLE ::  volhl(:,:,:,:,:)      !< volume of air defined on half levels
    REAL(wp), ALLOCATABLE ::  vol(:,:,:,:,:)        !< volume of air defined at main levels
    REAL(wp), ALLOCATABLE ::  width_build(:,:,:,:)  !< building width [m]
    REAL(wp), ALLOCATABLE ::  width_dblcan(:,:,:,:) !< width of canyon including building [m]
    REAL(wp), ALLOCATABLE ::  width_sglcan(:,:,:,:) !< width of single canyon
    REAL(wp), ALLOCATABLE ::  width_street(:,:,:,:) !< street width [m]
    REAL(wp), ALLOCATABLE ::  z0_ground(:)          !< ground's roughness length
    REAL(wp), ALLOCATABLE ::  z0_roof(:)            !< roof's roughness length
    REAL(wp), ALLOCATABLE ::  z0_urb(:,:)           !< average roughness length of urban ground surfaces
    REAL(wp), ALLOCATABLE ::  z_ground(:)           !< layer thickness of ground
    REAL(wp), ALLOCATABLE ::  z_meso(:,:,:)         !< height levels of the PALM starting at the bottom with 0
    REAL(wp), ALLOCATABLE ::  z_midmeso(:,:,:)      !< height levels of the PALM starting at the  bottom with 0
    REAL(wp), ALLOCATABLE ::  z_midu(:,:)           !< heigth of ith face of level in urban class, z_uhl(1,:)=0
    REAL(wp), ALLOCATABLE ::  z_roof(:)             !< layer thickness of roof
    REAL(wp), ALLOCATABLE ::  z_wall(:)             !< layer thickness of wall

    REAL(wp), DIMENSION(100) ::  albedo_ground     = 0.2_wp      !< albedo of the ground (read in)
    REAL(wp), DIMENSION(100) ::  albedo_roof       = 0.2_wp      !< albedo of the roof (read in)
    REAL(wp), DIMENSION(100) ::  albedo_wall       = 0.2_wp      !< albedo of the wall (read in)
    REAL(wp), DIMENSION(100) ::  ang_udir          = -9999.9_wp  !< street direction [degree], ang_udir(n_udir)
    REAL(wp), DIMENSION(100) ::  dzlayer_ground    = -1.0_wp     !< thickness of ground's layers (read in)
    REAL(wp), DIMENSION(100) ::  dzlayer_roof      = -1.0_wp     !< thickness of roof's layers (read in)
    REAL(wp), DIMENSION(100) ::  dzlayer_wall      = -1.0_wp     !< thickness of wall's layers (read in)
    REAL(wp), DIMENSION(100) ::  emissivity_ground = 0.95_wp     !< emissivity of ground (read in)
    REAL(wp), DIMENSION(100) ::  emissivity_roof   = 0.9_wp      !< emissivity of roof (read in)
    REAL(wp), DIMENSION(100) ::  emissivity_wall   = 0.9_wp      !< emissivity of wall (read in)
    REAL(wp), DIMENSION(100) ::  heatcap_ground    = 1.4E6_wp    !< specific heat of ground material [J m^3 K^-1] (read in)
    REAL(wp), DIMENSION(100) ::  heatcap_roof      = 1.E6_wp     !< specific heat of roof material [J m^3 K^-1] (read in)
    REAL(wp), DIMENSION(100) ::  heatcap_wall      = 1.E6_wp     !< specific heat of wall material [J m^3 K^-1] (read in)
    REAL(wp), DIMENSION(100) ::  rlength_ground    = 0.01_wp     !< ground's roughness length (namelist)
    REAL(wp), DIMENSION(100) ::  rlength_roof      = 0.01_wp     !< roof's roughness length (namelist)
    REAL(wp), DIMENSION(100) ::  thermdiff_ground  = 0.67E-6_wp  !< thermal diffusivity of ground's layers
    REAL(wp), DIMENSION(100) ::  thermdiff_roof    = 0.29E-6_wp  !< thermal diffusivity of roof's layers
    REAL(wp), DIMENSION(100) ::  thermdiff_wall    = 0.67E-6_wp  !< thermal diffusivity of wall's layers
    REAL(wp), DIMENSION(100) ::  tinterior_ground  = 292.0_wp    !< ground's initial and interior temperature (read in)
    REAL(wp), DIMENSION(100) ::  tinterior_roof    = 292.0_wp    !< roof's initial and interior temperature (read in)
    REAL(wp), DIMENSION(100) ::  tinterior_wall    = 292.0_wp    !< walls's initial and interior temperature (read in)
    REAL(wp), DIMENSION(100) ::  z_uhl(100)        = -99999.9_wp !< heigth of ith face of level in urban class (namelist)

    REAL(wp), DIMENSION(n_num_integ)  ::  legendrepoly_zero
    REAL(wp), DIMENSION(n_num_integ)  ::  weight

!
!-- PALM interfaces
    INTERFACE dcep_check_data_output
       MODULE PROCEDURE dcep_check_data_output
    END INTERFACE dcep_check_data_output

    INTERFACE dcep_check_parameters
       MODULE PROCEDURE dcep_check_parameters
    END INTERFACE dcep_check_parameters

    INTERFACE dcep_data_output_2d
       MODULE PROCEDURE dcep_data_output_2d
    END INTERFACE dcep_data_output_2d

    INTERFACE dcep_data_output_3d
       MODULE PROCEDURE dcep_data_output_3d
    END INTERFACE dcep_data_output_3d

    INTERFACE dcep_define_netcdf_grid
       MODULE PROCEDURE dcep_define_netcdf_grid
    END INTERFACE dcep_define_netcdf_grid

    INTERFACE dcep_init
       MODULE PROCEDURE dcep_init
    END INTERFACE dcep_init

    INTERFACE dcep_init_arrays
       MODULE PROCEDURE dcep_init_arrays
    END INTERFACE dcep_init_arrays

    INTERFACE dcep_main
       MODULE PROCEDURE dcep_main
    END INTERFACE dcep_main

    INTERFACE dcep_netcdf_input
       MODULE PROCEDURE dcep_netcdf_input
    END INTERFACE dcep_netcdf_input

    INTERFACE dcep_parin
       MODULE PROCEDURE dcep_parin
    END INTERFACE dcep_parin

    INTERFACE dcep_tendency
       MODULE PROCEDURE dcep_tendency
       MODULE PROCEDURE dcep_tendency_ji
    END INTERFACE dcep_tendency

    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC dcep_check_data_output,                                                                 &
           dcep_check_parameters,                                                                  &
           dcep_data_output_2d,                                                                    &
           dcep_data_output_3d,                                                                    &
           dcep_define_netcdf_grid,                                                                &
           dcep_init,                                                                              &
           dcep_init_arrays,                                                                       &
           dcep_main,                                                                              &
           dcep_netcdf_input,                                                                      &
           dcep_parin,                                                                             &
           dcep_tendency
!
!-- Public variables
    PUBLIC eps_urb,                                                                                &
           fr_urb,                                                                                 &
           ke_ground,                                                                              &
           ke_roof,                                                                                &
           ke_uhl,                                                                                 &
           ke_wall,                                                                                &
           nz_mesom,                                                                               &
           shfl_urb,                                                                               &
           z_ground,                                                                               &
           z_roof,                                                                                 &
           z_uhl,                                                                                  &
           z_wall

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> The main DCEP routine which performs the DCEP model calculations.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_main

    USE radiation_model_mod,                                                                       &
        ONLY:  radiation_calc_diffusion_radiation

    IMPLICIT NONE

    INTEGER(iwp) ::  day_of_year            !< day of the current year
    INTEGER(iwp) ::  id                     !< index for street direction
    INTEGER(iwp) ::  ii                     !< index
    INTEGER(iwp) ::  i                      !< index for spatial x direction
    INTEGER(iwp) ::  iu                     !< index for urban height
    INTEGER(iwp) ::  iurb                   !< index for urban class
    INTEGER(iwp) ::  iz                     !< index for mesoscale height
    INTEGER(iwp) ::  izz                    !< index for mesoscale height above topo
    INTEGER(iwp) ::  jj                     !< index
    INTEGER(iwp) ::  j                      !< index for spatial y direction
    INTEGER(iwp) ::  k_topo                 !< topography top index
    INTEGER(iwp) ::  nz_id                  !< index for max nz

    REAL(wp) ::  dzz                        !< overlap of mesoscale and urban height levels
    REAL(wp) ::  fact                       !< factor for single urban flux to yield mesoscale flux
    REAL(wp) ::  ptint(1:nz_mesom)          !< input of air temperature for interpolate air temp
    REAL(wp) ::  rho_airint(1:nz_mesom)     !< input of air density for interpolate air density
    REAL(wp) ::  rhoint(ke_uhl+1)           !< interpolated air density
    REAL(wp) ::  rtg(n_udir)                !< total radiation at the ground surface
    REAL(wp) ::  rtr(n_udir,ke_uhl+1)       !< total radiation at the roof surfaces
    REAL(wp) ::  rtw(2*n_udir,ke_uhl)       !< total radiation at the wall surfaces
    REAL(wp) ::  seb                        !< sum of horizontal TKE production due to buoyancy
    REAL(wp) ::  second_of_day              !< second of the current day
    REAL(wp) ::  ses                        !< sum of horizontal TKE production due to shear
    REAL(wp) ::  shflc                      !< total urban sensible flux of urban class
    REAL(wp) ::  shfld                      !< total urban sensible flux of street direction
    REAL(wp) ::  strflc                     !< total urban storage flux of urban class
    REAL(wp) ::  strfld                     !< total urban storage flux of street direction
    REAL(wp) ::  st                         !< sum of horizontal urban temperature fluxes
    REAL(wp) ::  sun_azi                    !< solar azimuth
    REAL(wp) ::  sun_el                     !< solar azimuth
    REAL(wp) ::  su                         !< sum of horizontal u fluxes
    REAL(wp) ::  sv                         !< sum of horizontal v fluxes
    REAL(wp) ::  t0int(ke_uhl+1)            !< interpolated reference temperature
    REAL(wp) ::  tfh_exp(n_udir,ke_uhl+1)   !< temperature        horizontal surfaces, B (explicit) term
    REAL(wp) ::  tfhg_exp(n_udir)           !< temperature        horizontal surfaces, B (explicit) term
    REAL(wp) ::  tfh_imp(2*n_udir,ke_uhl)   !< temperature          vertical surfaces, A (implicit) term
    REAL(wp) ::  tfv_exp(2*n_udir,ke_uhl)   !< temperature          vertical surfaces, B (explicit) term
    REAL(wp) ::  tint(ke_uhl+1)             !< interpolated temperature
    REAL(wp) ::  tkebc(nz_mesom+1)          !< TKE tendency buoyancy for urban class
    REAL(wp) ::  tkebh_exp(n_udir,ke_uhl+1) !< energy (TKE) buoyancy horizontal surfaces, B (explicit) term
    REAL(wp) ::  tkebhg_exp(n_udir)         !< energy (TKE) buoyancy horizontal surfaces, B (explicit) term
    REAL(wp) ::  tkesc(nz_mesom+1)          !< TKE tendency shear for urban class
    REAL(wp) ::  tkesh_exp(n_udir,ke_uhl+1) !< energy (TKE) shear    horizontal surfaces, B (explicit) term
    REAL(wp) ::  tkeshg_exp(n_udir)         !< energy (TKE) shear    horizontal surfaces, B (explicit) term
    REAL(wp) ::  tkesv_exp(2*n_udir,ke_uhl) !< energy (TKE) shear    vertical surfaces, B (explicit) term
    REAL(wp) ::  ttc(nz_mesom+1)            !< temperature tendency for urban class
    REAL(wp) ::  ufh_exp(n_udir,ke_uhl+1)   !< u (wind component) horizontal surfaces, B (explicit) term
    REAL(wp) ::  ufhg_exp(n_udir)           !< u (wind component) horizontal surfaces, B (explicit) term
    REAL(wp) ::  ufv_exp(2*n_udir,ke_uhl)   !< u (wind component)   vertical surfaces, B (explicit) term
    REAL(wp) ::  ufv_imp(2*n_udir,ke_uhl)   !< u (wind component)   vertical surfaces, A (implicit) term
    REAL(wp) ::  uint(ke_uhl+1)             !< interpolated u wind speed
    REAL(wp) ::  umflc                      !< total urban u momentum flux of urban class
    REAL(wp) ::  umfld                      !< total urban u momentum flux of street direction
    REAL(wp) ::  umid(1:nz_mesom)           !< u velocity defined at centre of grid cell
    REAL(wp) ::  utc(nz_mesom)              !< u component tendency for urban class
    REAL(wp) ::  veb
    REAL(wp) ::  vfh_exp(n_udir,ke_uhl+1)   !< v (wind component) horizontal surfaces, B (explicit) term
    REAL(wp) ::  vfhg_exp(n_udir)           !< v (wind component) horizontal surfaces, B (explicit) term
    REAL(wp) ::  vfv_exp(2*n_udir,ke_uhl)   !< v (wind component)   vertical surfaces, B (explicit) term
    REAL(wp) ::  vfv_imp(2*n_udir,ke_uhl)   !< v (wind component)   vertical surfaces, A (implicit) term
    REAL(wp) ::  vint(ke_uhl+1)             !< interpolated v wind speed
    REAL(wp) ::  vmflc                      !< total urban v momentum flux of urban class
    REAL(wp) ::  vmfld                      !< total urban v momentum flux of street direction
    REAL(wp) ::  vmid(1:nz_mesom)           !< v velocity defined at centre of grid cell
    REAL(wp) ::  vtb                        !< sum of vertical urban temperature fluxes
    REAL(wp) ::  vtc(nz_mesom)              !< v component tendency for urban class
    REAL(wp) ::  vub                        !< sum of vertical u fluxes
    REAL(wp) ::  vvb                        !< sum of vertical v fluxes

    
    IF ( debug_output )  CALL debug_message( 'dcep_main', 'start' )

    CALL cpu_log( log_point_s(22), 'DCEP', 'start' )
!
!-- Get the current sun angle.
    sun_direction = .TRUE.
    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year,                     &
                        second_of_day = second_of_day )
    CALL calc_zenith( day_of_year, second_of_day )
    sun_azi = ATAN2( sun_dir_lon, sun_dir_lat )
    sun_el  = ACOS( cos_zenith )
!
!--  Split downwelling shortwave radiation into a diffuse and a direct part. Note, if radiation
!--  scheme is RRTMG or diffuse radiation is externally prescribed, this is not required. Please
!--  note, in case of external radiation, the clear-sky model is applied during spinup, so that
!--  radiation needs to be split also in this case.
     IF ( radiation_scheme == 'constant'  .OR.  radiation_scheme == 'clear-sky'  .OR.              &
          radiation_scheme == 'external'  )  THEN
        CALL radiation_calc_diffusion_radiation
     ENDIF
!
!-- Call radiation routines
    CALL dcep_modify_lw( rad_lw_in_diff )
    CALL dcep_modify_sw( rad_sw_in_dir, rad_sw_in_diff, sun_el, sun_azi )
!
!-- Apply DCEP algorithm call radiation routines
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Initialize tendencies, fluxes and average temperature.
          ut_urb(:,j,i)   = 0.0_wp
          vt_urb(:,j,i)   = 0.0_wp
          tt_urb(:,j,i)   = 0.0_wp
          tkes_urb(:,j,i) = 0.0_wp
          tkeb_urb(:,j,i) = 0.0_wp

          strfl_urb(j,i) = 0.0_wp
          shfl_urb(j,i)  = 0.0_wp
          umfl_urb(j,i)  = 0.0_wp
          vmfl_urb(j,i)  = 0.0_wp

          t_g_urb(j,i) = 0.0_wp

          k_topo = topo_top_ind(j,i,0)
!
!--       Calculate grid cell centre wind velocity (velocity is 1/2 level right
!--       of centre), air temmperature and density.
          ii = MAX( i-1, nxl )
          jj = MAX( j-1, nys )
          DO  iz = 1, nz_meso(j,i)
             izz = iz + k_topo
             umid(iz) = 0.5_wp * ( u(izz,j,ii) + u(izz,j,i) )
             vmid(iz) = 0.5_wp * ( v(izz,jj,i) + v(izz,j,i) )
          ENDDO
          ptint(1:nz_meso(j,i)) = pt(k_topo+1:nz_meso(j,i),j,i)
          rho_airint(1:nz_meso(j,i)) = rho_air(k_topo+1:nz_meso(j,i))

          DO  iurb = 1, n_uclass
             IF ( fr_uclass(iurb,j,i) >= eps_urb )  THEN
!
!--             Initialize temporary urban class values.
                utc(:)   = 0.0_wp
                vtc(:)   = 0.0_wp
                ttc(:)   = 0.0_wp
                tkesc(:) = 0.0_wp
                tkebc(:) = 0.0_wp

                strflc = 0.0_wp
                shflc  = 0.0_wp
                umflc  = 0.0_wp
                vmflc  = 0.0_wp
!
!--             Heighest urban height level of all street directions.
                nz_id = MAXVAL( nz(iurb,:,j,i) )
!
!--             Interpolation on the "urban grid", main levels.
                uint(1:nz_id+1)   = interpol( nz_meso(j,i), nz_id, z_meso(1:nz_meso(j,i)+1,j,i),   &
                                              z_uhl(1:nz_id+2), umid(1:nz_meso(j,i)) )

                vint(1:nz_id+1)   = interpol( nz_meso(j,i), nz_id, z_meso(1:nz_meso(j,i)+1,j,i),   &
                                              z_uhl(1:nz_id+2), vmid(1:nz_meso(j,i)) )

                tint(1:nz_id+1)   = interpol( nz_meso(j,i), nz_id, z_meso(1:nz_meso(j,i)+1,j,i),   &
                                              z_uhl(1:nz_id+2), ptint(1:nz_meso(j,i)) )

                rhoint(1:nz_id+1) = interpol( nz_meso(j,i), nz_id, z_meso(1:nz_meso(j,i)+1,j,i),   &
                                              z_uhl(1:nz_id+2), rho_airint(1:nz_meso(j,i)) )
!
!--             Reference temperature = air temperature.
                t0int = tint
!
!--             Compute the implicit and explicit components of the sources or sinks on the
!--             "urban grid".
                DO  id = 1, n_udir

                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
!
!--                   Calculation at the ground surfaces.
                      CALL flux_flat( z_uhl(2), z0_ground(iurb), uint(1), vint(1), tint(1),        &
                                      t0int(1), t_ground(iurb,id,ke_ground,j,i), ufhg_exp(id),     &
                                      vfhg_exp(id), tfhg_exp(id), tkeshg_exp(id), tkebhg_exp(id) )
!
!--                   Calculation at the roof surfaces.
                      DO  iu = 1, nz(iurb,id,j,i) + 1
                         CALL flux_flat( dz_uhl(iu), z0_roof(iurb), uint(iu), vint(iu), tint(iu),  &
                                         t0int(iu), t_roof(iurb,id,ke_roof,iu,j,i), ufh_exp(id,iu),&
                                         vfh_exp(id,iu), tfh_exp(id,iu), tkesh_exp(id,iu),         &
                                         tkebh_exp(id,iu) )

                      ENDDO
!
!--                   Calculation at the wall surfaces.
                      DO  iu = 1, nz(iurb,id,j,i)

                         CALL flux_wall( uint(iu), vint(iu), tint(iu), rhoint(iu),                 &
                                         t_wall(iurb,2*id-1,ke_wall,iu,j,i), ufv_imp(2*id-1,iu),   &
                                         vfv_imp(2*id-1,iu), ufv_exp(2*id-1,iu),                   &
                                         vfv_exp(2*id-1,iu), tfh_imp(2*id-1,iu),                   &
                                         tfv_exp(2*id-1,iu), tkesv_exp(2*id-1,iu),                 &
                                         angrotx_udir(id,j,i), angroty_udir(id,j,i),               &
                                         cdrag_wall(iurb,id,j,i), dt_3d )

                         CALL flux_wall( uint(iu), vint(iu), tint(iu), rhoint(iu),                 &
                                         t_wall(iurb,2*id,ke_wall,iu,j,i), ufv_imp(2*id,iu),       &
                                         vfv_imp(2*id,iu), ufv_exp(2*id,iu), vfv_exp(2*id,iu),     &
                                         tfh_imp(2*id,iu), tfv_exp(2*id,iu), tkesv_exp(2*id,iu),   &
                                         angrotx_udir(id,j,i), angroty_udir(id,j,i),               &
                                         cdrag_wall(iurb,id,j,i), dt_3d )

                      ENDDO

                   ENDIF

                ENDDO
!
!--             Add tendencies explicitly.
                DO  id = 1, n_udir

                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN

                      DO  iu = 1, nz(iurb,id,j,i)
                         vfv_exp(2*id-1,iu) = vfv_exp(2*id-1,iu) + vfv_imp(2*id-1,iu) * vint(iu)
                         vfv_exp(2*id,iu)   = vfv_exp(2*id,iu)   + vfv_imp(2*id,iu)   * vint(iu)
                         ufv_exp(2*id-1,iu) = ufv_exp(2*id-1,iu) + ufv_imp(2*id-1,iu) * uint(iu)
                         ufv_exp(2*id,iu)   = ufv_exp(2*id,iu)   + ufv_imp(2*id,iu)   * uint(iu)
                         tfv_exp(2*id-1,iu) = tfv_exp(2*id-1,iu) + tfh_imp(2*id-1,iu) * tint(iu)
                         tfv_exp(2*id,iu)   = tfv_exp(2*id,iu)   + tfh_imp(2*id,iu)   * tint(iu)
                      ENDDO

                   ENDIF

                ENDDO
!
!--             Calculate effective explicit tendencies.
                IF ( limpufl )  THEN

                   DO  id = 1, n_udir
                      IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN

                         DO  iu = 1, nz(iurb,id,j,i)
                            vfv_exp(2*id-1, iu) = vfv_exp(2*id-1, iu) /                            &
                                                  ( 1.0_wp - vfv_imp(2*id-1,iu) * dt_3d )
                            vfv_exp(2*id,   iu) = vfv_exp(2*id, iu) /                              &
                                                  ( 1.0_wp - vfv_imp(2*id,iu)   * dt_3d )

                            ufv_exp(2*id-1, iu) = ufv_exp(2*id-1, iu) /                            &
                                                  ( 1.0_wp - ufv_imp(2*id-1,iu) * dt_3d )
                            ufv_exp(2*id, iu)   = ufv_exp(2*id,   iu) /                            &
                                                  ( 1.0_wp - ufv_imp(2*id,iu)   * dt_3d )

                            tfv_exp(2*id-1, iu) = tfv_exp(2*id-1, iu) /                            &
                                                  ( 1.0_wp - tfh_imp(2*id-1,iu) * dt_3d )
                            tfv_exp(2*id,   iu) = tfv_exp(2*id  , iu) /                            &
                                                  ( 1.0_wp - tfh_imp(2*id,iu)   * dt_3d )
                         ENDDO

                      ENDIF
                   ENDDO

                ENDIF

                DO  id = 1, n_udir

                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
!
!--                   Compute the surface temperatures.
!--                   1) ground temperature
                      CALL soil_temp( ke_ground, dz_ground(:,iurb), t_ground(iurb,id,:,j,i),       &
                                      tddz_ground(:,iurb), cs_ground(iurb), rs_ground(iurb,id,j,i),&
                                      rl_ground(iurb,id,j,i), dt_3d, emiss_ground(iurb),           &
                                      alb_ground(iurb), rtg(id), shfl_ground(iurb,id,j,i),         &
                                      strfl_ground(iurb,id,j,i) )
!
!--                   Average ground temperature.
                      t_g_urb(j,i) = t_g_urb(j,i) + t_ground(iurb,id,ke_ground,j,i) *              &
                                                    fr_udir(iurb,id,j,i) * fr_uclass(iurb,j,i)
!
!--                   2) roof surfaces
                      DO  iu = 1, nz(iurb,id,j,i)+1
                         CALL soil_temp( ke_roof, dz_roof(:,iurb), t_roof(iurb,id,:,iu,j,i),       &
                                         tddz_roof(:,iurb), cs_roof(iurb), rs_roof(iurb,id,iu,j,i),&
                                         rl_roof(iurb,id,iu,j,i), dt_3d, emiss_roof(iurb),         &
                                         alb_roof(iurb), rtr(id,iu), shfl_roof(iurb,id,iu,j,i),    &
                                         strfl_roof(iurb,id,iu,j,i) )
                      ENDDO
!
!--                   3) wall surfaces
                      DO  iu = 1, nz(iurb,id,j,i)

                         CALL soil_temp( ke_wall, dz_wall(:,iurb), t_wall(iurb,2*id-1,:,iu,j,i),   &
                                         tddz_wall(:,iurb), cs_wall(iurb),                         &
                                         rs_wall(iurb,2*id-1,iu,j,i), rl_wall(iurb,2*id-1,iu,j,i), &
                                         dt_3d, emiss_wall(iurb),  alb_wall(iurb), rtw(2*id-1,iu), &
                                         shfl_wall(iurb,2*id-1,iu,j,i),                            &
                                         strfl_wall(iurb,2*id-1,iu,j,i) )

                         CALL soil_temp( ke_wall, dz_wall(:,iurb), t_wall(iurb,2*id,:,iu,j,i),     &
                                         tddz_wall(:,iurb), cs_wall(iurb),                         &
                                         rs_wall(iurb,2*id,iu,j,i), rl_wall(iurb,2*id,iu,j,i),     &
                                         dt_3d, emiss_wall(iurb), alb_wall(iurb), rtw(2*id,iu),    &
                                         shfl_wall(iurb,2*id,iu,j,i), strfl_wall(iurb,2*id,iu,j,i) )
                      ENDDO
!
!--                   Sum up storage fluxes.
                      strfld = strfl_ground(iurb,id,j,i) * width_street(iurb,id,j,i)
                      DO  iu = 1, nz(iurb,id,j,i)
                         strfld = strfld + (                                                       &
                                     strfl_wall(iurb,2*id-1,iu,j,i) + strfl_wall(iurb,2*id,iu,j,i) &
                                           ) * dz_uhl(iu) * fr_wall(iurb,id,iu,j,i)
                      ENDDO
                      DO  iu = 1, nz(iurb,id,j,i)+1
                         strfld = strfld + strfl_roof(iurb,id,iu,j,i) * width_build(iurb,id,j,i)   &
                                           * fr_roof(iurb,id,iu,j,i)
                      ENDDO

                      strflc = strflc + strfld * fr_udir(iurb,id,j,i) / width_sglcan(iurb,id,j,i)

                   ENDIF

                ENDDO !id

                strfl_urb(j,i) = strfl_urb(j,i) + strflc * fr_uclass(iurb,j,i)
!
!--             Calculation of sensible heat and momentum fluxes.
                DO  id = 1, n_udir

                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN

                      shfl_ground(iurb,id,j,i) = -rhoint(1) * tfhg_exp(id) * c_p
                      umfl_ground(iurb,id,j,i) = -rhoint(1) * ufhg_exp(id)
                      vmfl_ground(iurb,id,j,i) = -rhoint(1) * vfhg_exp(id)

                      DO  iu = 1, nz(iurb,id,j,i)+1
                         shfl_roof(iurb,id,iu,j,i) = -rhoint(iu) * tfh_exp(id,iu) * c_p
                         umfl_roof(iurb,id,iu,j,i) = -rhoint(iu) * ufh_exp(id,iu)
                         vmfl_roof(iurb,id,iu,j,i) = -rhoint(iu) * vfh_exp(id,iu)
                      ENDDO
!
!--                   Implicit tendency already included in tvb
                      DO  iu = 1, nz(iurb,id,j,i)
                         shfl_wall(iurb,2*id-1,iu,j,i) = -rhoint(iu) * tfv_exp(2*id-1,iu) * c_p
                         shfl_wall(iurb,2*id,iu,j,i)   = -rhoint(iu) * tfv_exp(2*id,iu)   * c_p
                         umfl_wall(iurb,2*id-1,iu,j,i) = -rhoint(iu) * ufv_exp(2*id-1,iu)
                         umfl_wall(iurb,2*id,iu,j,i)   = -rhoint(iu) * ufv_exp(2*id,iu)
                         vmfl_wall(iurb,2*id-1,iu,j,i) = -rhoint(iu) * vfv_exp(2*id-1,iu)
                         vmfl_wall(iurb,2*id,iu,j,i)   = -rhoint(iu) * vfv_exp(2*id,iu)
                      ENDDO
!
!--                   Sum up fluxes for grid cell total.
                      shfld = shfl_ground(iurb,id,j,i) * width_street(iurb,id,j,i)
                      umfld = umfl_ground(iurb,id,j,i) * width_street(iurb,id,j,i)
                      vmfld = vmfl_ground(iurb,id,j,i) * width_street(iurb,id,j,i)

                      DO  iu = 1, nz(iurb,id,j,i)

                         shfld = shfld + (                                                         &
                                       shfl_wall(iurb,2*id-1,iu,j,i) + shfl_wall(iurb,2*id,iu,j,i) &
                                         ) * dz_uhl(iu) * fr_wall(iurb,id,iu,j,i)

                         umfld = umfld + (                                                         &
                                       umfl_wall(iurb,2*id-1,iu,j,i) + umfl_wall(iurb,2*id,iu,j,i) &
                                         ) * dz_uhl(iu) * fr_wall(iurb,id,iu,j,i)

                         vmfld = vmfld + (                                                         &
                                       vmfl_wall(iurb,2*id-1,iu,j,i) + vmfl_wall(iurb,2*id,iu,j,i) &
                                         ) * dz_uhl(iu) * fr_wall(iurb,id,iu,j,i)

                      ENDDO

                      DO  iu = 1, nz(iurb,id,j,i) + 1
                         shfld = shfld + shfl_roof(iurb,id,iu,j,i) * width_build(iurb,id,j,i)      &
                                         * fr_roof(iurb,id,iu,j,i)
                         umfld = umfld + umfl_roof(iurb,id,iu,j,i) * width_build(iurb,id,j,i)      &
                                         * fr_roof(iurb,id,iu,j,i)
                         vmfld = vmfld + vmfl_roof(iurb,id,iu,j,i) * width_build(iurb,id,j,i)      &
                                         * fr_roof(iurb,id,iu,j,i)
                      ENDDO

                      shflc = shflc + shfld * fr_udir(iurb,id,j,i) / width_sglcan(iurb,id,j,i)
                      umflc = umflc + umfld * fr_udir(iurb,id,j,i) / width_sglcan(iurb,id,j,i)
                      vmflc = vmflc + vmfld * fr_udir(iurb,id,j,i) / width_sglcan(iurb,id,j,i)

                   ENDIF
                ENDDO !id

                shfl_urb(j,i) = shfl_urb(j,i) + shflc * fr_uclass(iurb,j,i)
                umfl_urb(j,i) = umfl_urb(j,i) + umflc * fr_uclass(iurb,j,i)
                vmfl_urb(j,i) = vmfl_urb(j,i) + vmflc * fr_uclass(iurb,j,i)

!
!--             Interpolation on the "PALM grid".
!--             Horizontal surface impact except TKE.
                DO  id = 1, n_udir
                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
!
!--                   Street.
!--                   Surface/volume (Schubert et al 2013, Eq. 20) =
!--                   1/(volume correction factor * height)
                      fact = 1.0_wp / vol(iurb,id,1,j,i) / dz_meso(1,j,i) * fr_street(iurb,id,j,i) &
                             * fr_udir(iurb,id,j,i)
                      ttc(1) = ttc(1) + tfhg_exp(id) * fact
                      utc(1) = utc(1) + ufhg_exp(id) * fact
                      vtc(1) = vtc(1) + vfhg_exp(id) * fact

                      DO  iz = 1, nz_meso(j,i)
                         st  = 0.0_wp
                         su  = 0.0_wp
                         sv  = 0.0_wp
                         ses = 0.0_wp
                         seb = 0.0_wp
                         DO  iu = 1, nz(iurb,id,j,i) + 1
!
!--                         Only iu levels that are in iz PALM level.
                            IF ( z_meso(iz,j,i) < z_uhl(iu)  .AND.  z_meso(iz+1,j,i) > z_uhl(iu) ) &
                            THEN
                               st = st + fr_roof(iurb,id,iu,j,i) * tfh_exp(id,iu)
                               su = su + fr_roof(iurb,id,iu,j,i) * ufh_exp(id,iu)
                               sv = sv + fr_roof(iurb,id,iu,j,i) * vfh_exp(id,iu)
                            ENDIF
                         ENDDO
!
!--                      Flux/volume
                         fact = fr_build(iurb,id,j,i) / vol(iurb,id,iz,j,i) / dz_meso(iz,j,i) *    &
                                fr_udir(iurb,id,j,i)
                         ttc(iz) = ttc(iz) + st * fact
                         utc(iz) = utc(iz) + su * fact
                         vtc(iz) = vtc(iz) + sv * fact
                      ENDDO

                   ENDIF
                ENDDO
!
!--             Vertical surface impact except TKE.
                DO  iz = 1, nz_meso(j,i)
                   DO  id = 1, n_udir
                      IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
                         vtb = 0.0_wp
                         vub = 0.0_wp
                         vvb = 0.0_wp
                         DO  iu = 1, nz(iurb,id,j,i)
                            dzz = MAX( MIN( z_uhl(iu+1), z_meso(iz+1,j,i) ) -                      &
                                       MAX( z_uhl(iu), z_meso(iz,j,i) ) , 0.0_wp )
                            fact = dzz / width_sglcan(iurb,id,j,i)
                            vtb = vtb + fr_wall(iurb,id,iu,j,i)                                    &
                                        * ( tfv_exp(2*id-1,iu) + tfv_exp(2*id,iu) ) * fact
                            vub = vub + fr_wall(iurb,id,iu,j,i)                                    &
                                        * ( ufv_exp(2*id-1,iu) + ufv_exp(2*id,iu) ) * fact
                            vvb = vvb + fr_wall(iurb,id,iu,j,i)                                    &
                                        * ( vfv_exp(2*id-1,iu) + vfv_exp(2*id,iu) ) * fact
                         ENDDO

                         fact = 1.0_wp / vol(iurb,id,iz,j,i) / dz_meso(iz,j,i) *                   &
                                fr_udir(iurb,id,j,i)
                         ttc(iz) = ttc(iz) + vtb * fact
                         utc(iz) = utc(iz) + vub * fact
                         vtc(iz) = vtc(iz) + vvb * fact
                      ENDIF
                   ENDDO
                ENDDO
!
!--             TKE is defined on half levels.
!--             Horizontal surface impact.
                DO  id = 1, n_udir
                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
!
!--                   Street.
!--                   Surface/volume (Schubert et al 2013, Eq. 20) =
!--                   1/(volume correction factor * height).
                      fact = 1.0_wp / volhl(iurb,id,1,j,i) / z_midmeso(2,j,i)                      &
                             * fr_street(iurb,id,j,i) * fr_udir(iurb,id,j,i)
                      tkesc(1) = tkesc(1) + tkeshg_exp(id) * fact * dz_uhl(1) * 0.5_wp
                      tkebc(1) = tkebc(1) + tkebhg_exp(id) * fact * dz_uhl(1) * 0.5_wp
                      DO  iz = 1, nz_meso(j,i)
                         ses = 0.0_wp
                         seb = 0.0_wp
                         DO  iu = 1, nz(iurb,id,j,i)+1
                            dzz = MAX( MIN( 0.5_wp*(z_uhl(iu+1)+z_uhl(iu) ), z_midmeso(iz+1,j,i) ) &
                                     - MAX( z_uhl(iu) , z_midmeso(iz,j,i) ) , 0.0_wp )
                            ses = ses + fr_roof(iurb,id,iu,j,i) * tkesh_exp(id,iu) * dzz
                            seb = seb + fr_roof(iurb,id,iu,j,i) * tkebh_exp(id,iu) * dzz
                         ENDDO

!
!--                      Flux / volume.
                         fact = fr_build(iurb,id,j,i) / volhl(iurb,id,iz,j,i) /                    &
                                ( z_midmeso(iz+1,j,i) - z_midmeso(iz,j,i) ) * fr_udir(iurb,id,j,i)
                         tkesc(iz+1) = tkesc(iz+1) + ses * fact
                         tkebc(iz+1) = tkebc(iz+1) + seb * fact
                      ENDDO
                   ENDIF
                ENDDO
!
!--             Vertical surface impact.
                DO  iz = 1, nz_meso(j,i)
                   DO  id = 1, n_udir
                      IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
                         veb = 0.0_wp
                         DO  iu = 1, nz(iurb,id,j,i)
                            dzz = MAX( MIN( z_uhl(iu+1), z_midmeso(iz+1,j,i) ) -                   &
                                       MAX( z_uhl(iu  ), z_midmeso(iz,j,i) ), 0.0_wp )
                            fact = dzz / width_sglcan(iurb,id,j,i)
                            veb = veb + fr_wall(iurb,id,iu,j,i) *                                  &
                                        ( tkesv_exp(2*id-1,iu) + tkesv_exp(2*id,iu) ) * fact
                         ENDDO

                         fact = 1.0_wp / volhl(iurb,id,iz,j,i)                                     &
                                / ( z_midmeso(iz+1,j,i) - z_midmeso(iz,j,i) ) * fr_udir(iurb,id,j,i)
                         tkesc(iz+1) = tkesc(iz+1) + veb * fact
                      ENDIF
                   ENDDO
                ENDDO

                DO  iz = 1, nz_meso(j,i)
                   ut_urb(iz,j,i) = ut_urb(iz,j,i) + utc(iz) * fr_uclass(iurb,j,i)
                   vt_urb(iz,j,i) = vt_urb(iz,j,i) + vtc(iz) * fr_uclass(iurb,j,i)
                   tt_urb(iz,j,i) = tt_urb(iz,j,i) + ttc(iz) * fr_uclass(iurb,j,i)
                ENDDO

                DO  iz = 1, nz_meso(j,i)
                   tkes_urb(iz,j,i) = tkes_urb(iz,j,i) + tkesc(iz) * fr_uclass(iurb,j,i)
                   tkeb_urb(iz,j,i) = tkeb_urb(iz,j,i) + tkebc(iz) * fr_uclass(iurb,j,i)
                ENDDO

             ENDIF

          ENDDO ! iurb

       ENDDO  ! iy
    ENDDO  ! ix

    CALL cpu_log( log_point_s(22), 'DCEP', 'stop' )

    IF ( debug_output )  CALL debug_message( 'dcep_main', 'end' )

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of the flux at the ground, formulation of Louis (Louis, 1979)
!--------------------------------------------------------------------------------------------------!
 PURE SUBROUTINE flux_flat( dz, z0, ua, va, pt, pt0, ptg, uhb, vhb, thb, ehs, ehb )

    IMPLICIT NONE

    REAL(wp), INTENT(in) ::  dz   !< first vertical level
    REAL(wp), INTENT(in) ::  pt0  !< reference (potential) temperature
    REAL(wp), INTENT(in) ::  ptg  !< ground (potential) temperature
    REAL(wp), INTENT(in) ::  pt   !< (potential) temperature
    REAL(wp), INTENT(in) ::  ua   !< wind speed
    REAL(wp), INTENT(in) ::  va   !< wind speed
    REAL(wp), INTENT(in) ::  z0   !< roughness length
!
!-- Explicit component of the momentum, temperature and TKE sources
!-- or sinks on horizontal surfaces (roofs and street).
!-- The fluxes can be computed as follow: Fluxes of X = B
!-- Example: Momentum fluxes on horizontal surfaces =  ufh_exp
    REAL(wp), INTENT(out) ::  ehb   !< energy (TKE) buoyancy horizontal surfaces, B (explicit) term
    REAL(wp), INTENT(out) ::  ehs   !< energy (TKE) shear horizontal surfaces, B (explicit) term
    REAL(wp), INTENT(out) ::  thb   !< temperature        horizontal surfaces, B (explicit) term
    REAL(wp), INTENT(out) ::  uhb   !< u (wind component) horizontal surfaces, B (explicit) term
    REAL(wp), INTENT(out) ::  vhb   !< v (wind component) horizontal surfaces, B (explicit) term
!
!-- Local variables.
    REAL(wp) ::  aa      !<
    REAL(wp) ::  c       !<
    REAL(wp) ::  fbuw    !<
    REAL(wp) ::  fh      !<
    REAL(wp) ::  fm      !<
    REAL(wp) ::  ric     !<
    REAL(wp) ::  tchvel  !<
    REAL(wp) ::  tcmvel  !<
    REAL(wp) ::  utot    !<
    REAL(wp) ::  utotsq  !<
    REAL(wp) ::  zz      !<
!
!-- Emperical constants.
    REAL(wp), PARAMETER ::  b   = 9.4_wp
    REAL(wp), PARAMETER ::  cm  = 7.4_wp
    REAL(wp), PARAMETER ::  ch  = 5.3_wp
    REAL(wp), PARAMETER ::  rr  = 0.74_wp


!
!-- Size of difference of horizontal windspeed at level IU (1/2 level above flat surface iu)
!-- and at level iu (at the surface so zero wind speed).
    utotsq = ua**2 + va**2
    utotsq = MAX( utotsq, 0.0001_wp )
    utot   = SQRT( utotsq )
!
!-- Louis formulation: compute the bulk Richardson Number.
    zz = dz * 0.5_wp
!
!-- Here also usage of non-potential temperature is possible because no (or just small) pressure
!-- difference between both height levels.
!-- ric= 2.0 * g * zz * (pt-ptg) / ((pt+ptg) * (utotsq))
!-- dz/2 * 2 cancels (latter from pt average)
    ric = g * dz * ( pt - ptg ) / ( ( pt + ptg ) * utotsq )

    aa = kappa / LOG( zz / z0 )
!
!-- Determine the parameters fm and fh for stable, neutral and unstable conditions.
    IF ( ric > 0.0_wp )  THEN
!
!--    Neutral
       fm = 1.0_wp / ( 1.0_wp + 0.5_wp * b * ric )**2
       fh = fm
    ELSE
!
!--    Unstable
       c  = b * cm * aa * aa * ( zz / z0 )**0.5_wp
       fm = 1.0_wp - b * ric / ( 1.0_wp + c * SQRT( -ric ) )
       c  = c * ch / cm
       fh = 1.0_wp - b * ric / ( 1.0_wp + c * SQRT( -ric ) )
    ENDIF
!
!-- Transfer coefficitent * wind velocity.
    tchvel = aa * aa *utot
    tcmvel = tchvel * fm
    tchvel = tchvel * fh / rr

    fbuw = -tcmvel * utot
!
!-- If interested in temperature tendencies, non-potential is ok.
    thb  = -tchvel * ( pt - ptg )
!
!-- Direction
    uhb = fbuw * ua / utot
    vhb = fbuw * va / utot

    ehs = ( ABS( fbuw ) )**1.5_wp / ( aa * zz * SQRT( fm ) )
    ehb = g / pt0 * thb

 END SUBROUTINE flux_flat


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine computes the surface sources or sinks of momentum, tke,
!> and heat from vertical surfaces (walls).
!--------------------------------------------------------------------------------------------------!
 PURE SUBROUTINE flux_wall( ua, va, pt, da, ptw, uva, vva, uvb, vvb, tva, tvb,evs, drstx, drsty,   &
                            cdrag, dt )

    IMPLICIT NONE

    REAL(wp), INTENT(in) ::  cdrag  !< drag coefficient
    REAL(wp), INTENT(in) ::  da     !< air density
    REAL(wp), INTENT(in) ::  drstx  !< street directions for the current urban class, x component
    REAL(wp), INTENT(in) ::  drsty  !< street directions for the current urban class, y component
    REAL(wp), INTENT(in) ::  dt     !< time step
    REAL(wp), INTENT(in) ::  pt     !< potential temperature
    REAL(wp), INTENT(in) ::  ptw    !< Walls potential temperatures
    REAL(wp), INTENT(in) ::  ua     !< wind speed
    REAL(wp), INTENT(in) ::  va     !< wind speed
!
!-- Explicit and implicit component of the momentum, temperature and TKE sources or sinks on
!   vertical surfaces (walls). The fluxes can be computed as follow:
!   Fluxes of X = A*X + B Example:
!   Momentum fluxes on vertical surfaces = ufv_imp * uint + ufv_exp
    REAL(wp), INTENT(out) ::  evs    !< Energy (TKE) shear   Vertical surfaces, B (explicit) term
    REAL(wp), INTENT(out) ::  tva    !< Temperature          Vertical surfaces, A (implicit) term
    REAL(wp), INTENT(out) ::  tvb    !< Temperature          Vertical surfaces, B (explicit) term
    REAL(wp), INTENT(out) ::  uva    !< U (wind component)   Vertical surfaces, A (implicit) term
    REAL(wp), INTENT(out) ::  uvb    !< U (wind component)   Vertical surfaces, B (explicit) term
    REAL(wp), INTENT(out) ::  vva    !< V (wind component)   Vertical surfaces, A (implicit) term
    REAL(wp), INTENT(out) ::  vvb    !< V (wind component)   Vertical surfaces, B (explicit) term
!
!-- Local variables
    REAL(wp) ::  hc    !<
    REAL(wp) :: u_ort  !<
    REAL(wp) :: vett   !<


    vett = SQRT( ua**2 + va**2 )
!
!-- Street direction: cos(street angle) = drstx
!--                   sin(street angle) = drsty
!-- Orthogonal: cos(street angle + 90) = -sin(street angle) = -drsty
!--             sin(street angle + 90) =  cos(street angle) =  drstx
!--
!-- Length of orthogonal wind:
    u_ort = ABS( -drsty * ua + drstx * va )
!
!-- Orthogonal wind vector: (-drsty*ua + drstx*va) * (-drsty drstx)
    uva = -cdrag * u_ort / 2.0_wp * drsty * drsty
    vva = -cdrag * u_ort / 2.0_wp * drstx * drstx

    uvb =  cdrag * u_ort / 2.0_wp * drstx * drsty * va
    vvb =  cdrag * u_ort / 2.0_wp * drstx * drsty * ua

    hc = 5.678_wp * (1.09_wp + 0.23_wp *(vett / 0.3048_wp))
    hc = MIN( hc, da * c_p / dt )

    tvb = hc * ptw / da / c_p - hc / da / c_p * pt
    tva = 0.0_wp

    evs = cdrag * u_ort**3 * 0.5_wp

 END SUBROUTINE flux_wall


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  This routine solves the Fourier diffusion equation for heat in the material (wall, roof,
!>  or ground). Solution is done implicitely.
!>  Boundary conditions are:
!>  - fixed temperature at the interior
!>  - energy budget at the surface
!>
!>           | K1  | K2  | K3  | K4  |
!>  fixed T  |  |  |  |  |  |  |  |  |  energy budget
!>  or         T1    T2    T3    T4
!>  zero flux \___/ \___/ \___/ \___/
!>             dz1   dz2   dz3   dz4
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE soil_temp( nz, dz, temp, tddz, cs, rs, rl, dt, em, alb, rt, sf, gf )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in):: nz         !< Number of layers

    REAL(wp), INTENT(in) ::  alb          !< albedo of the surface
    REAL(wp), INTENT(in) ::  cs           !< specific volumetric heat of the material [J m^-3 K^-1]
    REAL(wp), INTENT(in) ::  dt           !< time step
    REAL(wp), INTENT(in) ::  dz(nz)       !< layer sizes [m]
    REAL(wp), INTENT(in) ::  em           !< emissivity of the surface
    REAL(wp), INTENT(in) ::  rl           !< downward flux of the longwave radiation
    REAL(wp), INTENT(in) ::  rs           !< solar radiation
    REAL(wp), INTENT(in) ::  tddz(nz-1)   !< average thermal diffusivity divided by length [m s^-1]

    REAL(wp), INTENT(inout) ::  sf        !< sensible heat flux at the surface; can be modified due to temperature limit
    REAL(wp), INTENT(inout) ::  temp(nz)  !< temperature in each layer [K]

    REAL(wp), INTENT(out) ::  gf          !< heat flux transferred from the surface toward the interior
    REAL(wp), INTENT(out) ::  rt          !< total radiation at the surface (solar+incoming long+outgoing long)
!
!-- Local variables
    INTEGER(iwp) ::  iz

    REAL(wp), PARAMETER ::  deltatmax = 7.0_wp  !< maximum allowed temperature change of the outermost surface layer in one time step

    REAL(wp) ::  alpha
    REAL(wp) ::  alphatemp
    REAL(wp) ::  a(nz,3)
    REAL(wp) ::  c(nz)
    REAL(wp) ::  deltat
    REAL(wp) ::  dtiadz
    REAL(wp) ::  sftemp


    alpha = ( 1.0_wp - alb ) * rs + em * rl - em * sigma_sb * temp(nz)**4 + sf
!
!-- Limit the temperature difference of the uppermost soil layer by
!-- reducing the sensible heat flux estimation is done explicitely.
    dtiadz = dt / dz(nz)
    deltat = dtiadz * ( alpha / cs - tddz(nz-1) * ( temp(nz) - temp(nz-1) ) )
    IF ( ABS( deltat ) > deltatmax )  THEN
       alphatemp = cs *                                                                            &
                   ( SIGN( deltatmax,deltat ) / dtiadz + tddz(nz-1) * ( temp(nz) - temp(nz-1) ) )
       sftemp = alphatemp - ( ( 1.0_wp-alb ) * rs + em * rl - em * sigma_sb * temp(nz)**4 )
       alpha  = alphatemp
       sf     = sftemp
    ENDIF

    IF ( ltintfix )  THEN
!
!--    Fixed inner temperature T(1,t+1) = T(1,t).
       a(1,1) = 0.0_wp
       a(1,2) = 1.0_wp
       a(1,3) = 0.0_wp
       c(1)   = temp(1)
    ELSE
       dtiadz = dt / dz(1)
!
!--    No flux inside.
       a(1,1) = 0.0_wp
       a(1,2) = 1.0_wp + tddz(1) * dtiadz
       a(1,3) = -tddz(1) * dtiadz
       c(1)   = temp(1)
    ENDIF

    DO  iz = 2, nz-1
       dtiadz  =  dt / dz(iz)
       a(iz,1) = -tddz(iz-1) * dtiadz
       a(iz,2) = 1.0_wp + ( tddz(iz-1) + tddz(iz) ) * dtiadz
       a(iz,3) = -tddz(iz) * dtiadz
       c(iz)   = temp(iz)
    ENDDO

!
!-- Outer boundary.
    dtiadz  = dt / dz(nz)
    a(nz,1) = -tddz(nz-1) * dtiadz
    a(nz,2) = 1.0_wp + tddz(nz-1) * dtiadz
    a(nz,3) = 0.0_wp
    c(nz)   = temp(nz) + alpha / cs * dtiadz

!
!-- Crank Nicholson.
    CALL invert( nz, a, c, temp )
    rt = ( 1.0_wp - alb ) * rs + em * rl - em * sigma_sb * temp(nz)**4
    gf = ( 1.0_wp - alb ) * rs + em * rl - em * sigma_sb * temp(nz)**4 + sf

 END SUBROUTINE soil_temp


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Inversion and resolution of a tridiagonal matrix
!>               A X = C
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE invert( n, a, c, x )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) :: n                 !< dimension
    REAL(wp), DIMENSION(n,3), INTENT(inout) :: a  !<  a(*,1) lower diagonal (Ai,i-1)
                                                  !<  a(*,2) principal diagonal (Ai,i)
                                                  !<  a(*,3) upper diagonal (Ai,i+1)
    REAL(wp), INTENT(inout) :: c(n)  !<

    REAL(wp), INTENT(out), DIMENSION(n) :: x  !<

    INTEGER :: i  !<


    DO  i = n-1, 1, -1
       c(i)   = c(i)   - a(i,3) * c(i+1)   / a(i+1,2)
       a(i,2) = a(i,2) - a(i,3) * a(i+1,1) / a(i+1,2)
    ENDDO

    DO  i = 2, n
       c(i) = c(i) - a(i,1) * c(i-1) / a(i-1,2)
    ENDDO

    DO  i = 1, n
       x(i) = c(i)/a(i,2)
    ENDDO

 END SUBROUTINE invert


 END SUBROUTINE dcep_main


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for radiation model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_check_parameters

    USE control_parameters,                                                                        &
        ONLY:  urban_surface

    USE surface_mod,                                                                               &
        ONLY:  vertical_surfaces_exist

    IMPLICIT NONE

!
!-- Stop in case of urban-surface because DCEP is not implemented for USM
    IF ( urban_surface )  THEN
       message_string = 'DCEP is not implemented when urban-surface model is applied'
       CALL message( 'dcep_check_parameters', 'DCP0001', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Stop in case of vertical land-surface because DCEP is not implemented for vertical LSM surfaces
    IF ( vertical_surfaces_exist )  THEN
       message_string = 'DCEP is not implemented when vertical land surfaces exist'
       CALL message( 'dcep_check_parameters', 'DCP0002', 1, 2, 0, 6, 0 )
    ENDIF

 END SUBROUTINE dcep_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Main radiation routine for longwave radiation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_modify_lw( rl )

    IMPLICIT NONE

    INTEGER(iwp) ::  id              !<
    INTEGER(iwp) ::  i               !< loop index
    INTEGER(iwp) ::  iu              !<
    INTEGER(iwp) ::  iurb            !<
    INTEGER(iwp) ::  j               !< loop index
    INTEGER(iwp) ::  raddim          !< dimension of radiation matrices

    REAL(wp), INTENT(in) ::  rl(nysg:nyng,nxlg:nxrg)  !< diffuse longwave radiation

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  bbb   !< matrix and inhomogeneity

    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) :: rlup

    REAL(wp) ::  rldown1dg            !<
    REAL(wp) ::  rldown1dr(ke_uhl+1)  !<
    REAL(wp) ::  rldown1dwe(ke_uhl)   !< down east wall
    REAL(wp) ::  rldown1dww(ke_uhl)   !< down west wall


    ALLOCATE( bbb(raddimfull) )

    DO  j = nys, nyn
       DO  i = nxl, nxr
          IF ( fr_urb(j,i) > eps_urb )  THEN
             DO  iurb = 1, n_uclass

                IF ( fr_uclass(iurb,j,i) >= eps_urb )  THEN

                   DO  id = 1, n_udir
                      IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
!
!--                      Assume sky radiation is the same everywhere.
!--                      Down west wall:
                         rldown1dww = rl(j,i)

!--                      Down east wall:
                         rldown1dwe = rl(j,i)
                         rldown1dr  = rl(j,i)
                         rldown1dg  = rl(j,i)

                         IF ( lrroofs )  THEN
                            raddim = 3 * nz(iurb,id,j,i) + 2
                         ELSE
                            raddim = 2 * nz(iurb,id,j,i) + 1
                         ENDIF
!
!--                      Calculate what's actually arriving on every surface.
                         bbb(:) = long_rad_inhom( i, j, id, iurb, raddimfull, raddim,              &
                                                  rldown1dww(1:nz(iurb,id,j,i)),                   &
                                                  rldown1dwe(1:nz(iurb,id,j,i)),                   &
                                                  rldown1dr(1:nz(iurb,id,j,i)+1), rldown1dg,       &
                                                  radcor(iurb,id,j,i) )
!
!--                      Solve radiation system.
                         bbb(1:raddim) = MATMUL( rlmatrix(iurb,id,1:raddim,1:raddim,j,i),          &
                                                 bbb(1:raddim) )
!
!--                      West wall.
                         DO  iu = 1, nz(iurb,id,j,i)
                            rl_wall(iurb,2*id-1,iu,j,i) = bbb(iu)
                         ENDDO

!--                      East wall.
                         DO  iu = nz(iurb,id,j,i) + 1, 2 * nz(iurb,id,j,i)
                            rl_wall(iurb,2*id,iu-nz(iurb,id,j,i),j,i) = bbb(iu)
                         ENDDO

                         rl_ground(iurb,id,j,i) = bbb(2*nz(iurb,id,j,i)+1)

                         IF ( lrroofs )  THEN
                            DO  iu = 2*nz(iurb,id,j,i)+2, 3*nz(iurb,id,j,i)+2
                               rl_roof(iurb,id,iu-(2*nz(iurb,id,j,i)+1),j,i) = bbb(iu)
                            ENDDO
                         ELSE
!
!--                         Simple approach: roofs get full radiation.
                            DO  iu = 1, nz(iurb,id,j,i)+1
                               rl_roof(iurb,id,iu,j,i) = rl(j,i)
                            ENDDO
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO

          ENDIF
       ENDDO
    ENDDO
!
!-- Calculate total emitted radiation.
    rlup = uprad( rl, rl_ground, rl_wall, rl_roof,  1.0_wp-emiss_ground(:),                        &
                  1.0_wp-emiss_wall(:), 1.0_wp-emiss_roof(:), .TRUE. )

    DO  j = nys, nyn
       DO  i = nxl, nxr

          IF ( fr_urb(j,i) > eps_urb )  THEN
             t_grad_urb(j,i) = ( ( rlup(j,i) - ( 1.0_wp - emiss_urb(j,i) ) * rl(j,i) ) /           &
                               sigma_sb/emiss_urb(j,i) )**0.25_wp
          ELSE
             emiss_urb(j,i)  = 0.0_wp
             t_grad_urb(j,i) = 0.0_wp
          ENDIF

       ENDDO
    ENDDO

    DEALLOCATE( bbb )

 END SUBROUTINE dcep_modify_lw


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Main radiation routine for shortwave radiation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_modify_sw( rs, rsdiff,sel, saz )

    IMPLICIT NONE

    REAL(wp), INTENT(in) ::  rsdiff(nys:nyn,nxl:nxr) !< diffuse solar
    REAL(wp), INTENT(in) ::  rs(nys:nyn,nxl:nxr)     !< direkt solar radiation
    REAL(wp), INTENT(in) ::  saz                     !< azimuth from radiation
    REAL(wp), INTENT(in) ::  sel                     !< cos zenith angle
!
!-- Local variables
    INTEGER(iwp) ::  i,j,iurb,id,iu,ju,ku            !< index
    INTEGER(iwp) ::  raddim                          !< radiation matrices

    REAL(wp) ::  aae     !<
    REAL(wp) ::  bsd     !<
    REAL(wp) ::  rd      !<
    REAL(wp) ::  rd2     !<
    REAL(wp) ::  rd3     !<
    REAL(wp) ::  temp    !<
    REAL(wp) ::  tzr1d   !<
    REAL(wp) ::  wsd     !<
    REAL(wp) ::  wsbsd   !<
    REAL(wp) ::  ws2bsd  !<

    REAL(wp), DIMENSION(:),  ALLOCATABLE :: bbb      !< matrix and inhomogeneity
    REAL(wp), DIMENSION(nys:nyn,nxl:nxr) :: rsup

    REAL(wp) ::  rsdiff1dg(n_udir)                   !<
    REAL(wp) ::  rsdiff1dr(ke_uhl+1,n_udir)          !<
    REAL(wp) ::  rsdiff1dwe(ke_uhl,n_udir)           !< down east wall
    REAL(wp) ::  rsdiff1dww(ke_uhl,n_udir)           !< down west wall
    REAL(wp) ::  rsdown1dg(n_udir)                   !<
    REAL(wp) ::  rsdown1dr(ke_uhl+1,n_udir)          !<
    REAL(wp) ::  rsdown1dwe(ke_uhl,n_udir)           !< down east wall
    REAL(wp) ::  rsdown1dww(ke_uhl,n_udir)           !< down west wall


!
!-- Inhomogenity is time dependent, need it anyways.
    ALLOCATE( bbb(raddimfull) )

    DO  i = nxl, nxr
       DO  j = nys, nyn
          IF ( fr_urb(j,i) > eps_urb )  THEN
!
!--          If no insolation, skip this site.
             IF ( rs(j,i) + rsdiff(j,i) < eps_urb )  THEN
                DO  iurb = 1, n_uclass
                   DO  id = 1, n_udir
                      DO  iu = 1, nz(iurb,id,j,i)
                         rs_wall(iurb,2*id-1,iu,j,i) = 0.0_wp
                         rs_wall(iurb,2*id,iu,j,i)   = 0.0_wp
                      ENDDO
                      rs_ground(iurb,id,j,i) = 0.0_wp
                      DO  iu = 1, nz(iurb,id,j,i)+1
                         rs_roof(iurb,id,iu,j,i) = 0.0_wp
                      ENDDO
                   ENDDO
                ENDDO
!
!--             Keep albedo from step before.
                CYCLE
             ENDIF

!
!--          The principle sky radiation is rsdiff everywhere.
             DO  iu = 1, ke_uhl
                rsdiff1dww(iu,:) = rsdiff(j,i)
                rsdiff1dwe(iu,:) = rsdiff(j,i)
                rsdiff1dr(iu,:)  = rsdiff(j,i)
                rsdiff1dg(:)     = rsdiff(j,i)
             ENDDO

             rsdiff1dr(ke_uhl+1,:) = rsdiff(j,i)
!
!--          COTAN of sun elevation.
             tzr1d = 1.0_wp / TAN( sel )

             DO  iurb = 1, n_uclass

                IF ( fr_uclass(iurb,j,i) >= eps_urb )  THEN
!
!--                Shadow calculation (direct radiation).
                   IF ( .NOT. lrroofs )  THEN

!--                   Every building has the same height.
                      DO  id = 1, n_udir
                         IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
                            rsdown1dg(id) = 0.0_wp
!
!--                         Use unrotated as in radiation.
                            aae = saz - ang_udir(id)

                            temp = SIN( aae )

                            DO  iu = 1, nz(iurb,id,j,i)
                               rsdown1dww(iu,id) = 0.0_wp
                               rsdown1dwe(iu,id) = 0.0_wp
!
!--                            Modify incoming solar direct.
!--                            Radiation at the wall.
                               IF ( ABS( temp ) > 1.0E-10 )  THEN
                                  IF ( fr_wall(iurb,id,iu,j,i) > 0.0_wp ) THEN
!
!--                                  Opposite building has to be at least as large as receiving
!--                                  element iu.
                                     DO  ju = iu, nz(iurb,id,j,i)
                                        rd = shade_wall( z_uhl(iu), z_uhl(iu+1), z_uhl(ju+1),      &
                                                         tzr1d, temp, width_street(iurb,id,j,i) )
                                        rsdown1dww(iu,id) = rsdown1dww(iu,id) +                    &
                                                            rs(j,i) * rd * fr_roof(iurb,id,ju+1,j,i)

                                        rd = shade_wall( z_uhl(iu), z_uhl(iu+1), z_uhl(ju+1),      &
                                                         tzr1d, -temp, width_street(iurb,id,j,i) )
                                        rsdown1dwe(iu,id) = rsdown1dwe(iu,id) +                    &
                                                            rs(j,i) * rd * fr_roof(iurb,id,ju+1,j,i)
                                     ENDDO
!
!-                                   Normalize (total sum of added probablilty above).
                                     rsdown1dww(iu,id) = rsdown1dww(iu,id) / fr_wall(iurb,id,iu,j,i)
                                     rsdown1dwe(iu,id) = rsdown1dwe(iu,id) / fr_wall(iurb,id,iu,j,i)
                                  ENDIF
                               ENDIF
                            ENDDO

                            IF ( ABS( temp ) > 1.0E-10 )  THEN
!
!--                            street width / sin chi
                               wsd = ABS( width_street(iurb,id,j,i) / temp )
!
!--                            We have to take the 0 height buildings into account, so ju is roof
!--                            height here.
                               DO  ju = 1, nz(iurb,id,j,i) + 1
                                  rd = MAX( 0.0_wp, wsd-z_uhl(ju) * tzr1d )
                                  rsdown1dg(id) = rsdown1dg(id) + rs(j,i) * rd *                   &
                                                                  fr_roof(iurb,id,ju,j,i) / wsd
                               ENDDO
                            ELSE
                               rsdown1dg(id) = rs(j,i)
                            ENDIF

                            DO  iu = 1, nz(iurb,id,j,i) + 1
                               rsdown1dr(iu,id) = rs(j,i)
                            ENDDO
                         ENDIF
                      ENDDO

                   ELSE

                      DO  id = 1, n_udir
                         IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
                            rsdown1dg(id) = 0.0_wp
!
!--                         Here we use unrotated as in radiation.
                            aae = saz - ang_udir(id)

                            temp = SIN( aae )

                            DO  iu = 1, nz(iurb,id,j,i)
                               rsdown1dww(iu,id) = 0.0_wp
                               rsdown1dwe(iu,id) = 0.0_wp
!
!--                            Modify incoming solar direct.
!--                            Radiation at the wall.
                               IF ( ABS( temp ) > 1.0E-10 )  THEN
                                  IF( fr_wall(iurb,id,iu,j,i) > 0.0_wp )  THEN
                                     DO  ku = iu+1, nz(iurb,id,j,i)+1
                                        rd = shade_wall( z_uhl(iu), z_uhl(iu+1), z_uhl(ku),        &
                                                         tzr1d, temp, width_dblcan(iurb,id,j,i) )
                                        DO  ju = 1, nz(iurb,id,j,i)+1
                                           IF ( ju >= iu+1 )  THEN
                                              rd2 = shade_wall( z_uhl(iu), z_uhl(iu+1), z_uhl(ju), &
                                                             tzr1d,temp, width_street(iurb,id,j,i) )
                                              rd3 = MIN( rd, rd2 )
                                           ELSE
                                              rd3 = rd
                                           ENDIF
                                           rsdown1dww(iu,id) = rsdown1dww(iu,id) + rs(j,i) * rd3 * &
                                                                         fr_roof(iurb,id,ju,j,i) * &
                                                                         fr_roof(iurb,id,ku,j,i) / &
                                                                         fr_wall(iurb,id,iu,j,i)
                                        ENDDO

                                        rd = shade_wall( z_uhl(iu), z_uhl(iu+1), z_uhl(ku),        &
                                                         tzr1d, -temp, width_dblcan(iurb,id,j,i) )
                                        DO  ju = 1, nz(iurb,id,j,i)+1
                                           IF ( ju >= iu+1 )  THEN
                                              rd2 = shade_wall( z_uhl(iu), z_uhl(iu+1), z_uhl(ju), &
                                                           tzr1d, -temp, width_street(iurb,id,j,i) )
                                              rd3 = MIN( rd, rd2 )
                                           ELSE
                                              rd3 = rd
                                           ENDIF
                                           rsdown1dwe(iu,id) = rsdown1dwe(iu,id) + rs(j,i) * rd3 * &
                                                                         fr_roof(iurb,id,ju,j,i) * &
                                                                         fr_roof(iurb,id,ku,j,i) / &
                                                                         fr_wall(iurb,id,iu,j,i)
                                        ENDDO
                                     ENDDO
                                  ENDIF
                               ENDIF
                            ENDDO

                            IF ( ABS( temp ) > 1.0E-10 )  THEN
!
!--                            Street width / sin chi.
                               wsd    = ABS( width_street(iurb,id,j,i) / temp )
                               bsd    = ABS( width_build(iurb,id,j,i)  / temp )
                               ws2bsd = ABS( width_dblcan(iurb,id,j,i) / temp )
                               wsbsd  = ABS( width_sglcan(iurb,id,j,i) / temp )
!
!--                            We have to take the with 0 height buildings into account, so ju is
!--                            roof height here.
                               DO  ju = 1, nz(iurb,id,j,i)+1
                                  rsdown1dr(ju,id) = 0.0_wp

                                  rd = wsd - z_uhl(ju) * tzr1d

                                  DO  ku = 1, nz(iurb,id,j,i)+1

                                     IF ( ju >= ku )  THEN
                                        rsdown1dr(ju,id) = rsdown1dr(ju,id) +                      &
                                                                       fr_roof(iurb,id,ku,j,i) * bsd
                                     ELSE
                                        rsdown1dr(ju,id) = rsdown1dr(ju,id) +                      &
                                                               fr_roof(iurb,id,ku,j,i) *           &
                                                               MIN( bsd, MAX( 0.0_wp, wsbsd -      &
                                                               ( z_uhl(ku) - z_uhl(ju) ) * tzr1d ) )
                                     ENDIF

                                     rd2 = ws2bsd - z_uhl(ku) * tzr1d
                                     rd = MAX( 0.0_wp, MIN( rd,rd2 ) )
                                     rsdown1dg(id) = rsdown1dg(id) + rd * fr_roof(iurb,id,ju,j,i)  &
                                                                        * fr_roof(iurb,id,ku,j,i)
                                  ENDDO
                                  rsdown1dr(ju,id) = rsdown1dr(ju,id) * rs(j,i) / bsd
                               ENDDO
                               rsdown1dg(id) = rsdown1dg(id) * rs(j,i) / wsd
                            ELSE
                               rsdown1dg(id) = rs(j,i)
                               rsdown1dr(1:nz(iurb,id,j,i)+1,id) = rs(j,i)
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDIF  ! end shadow calculation

                   DO  id = 1, n_udir

                      IF ( fr_udir(iurb,id,j,i) > eps_urb)  THEN

                         IF ( lrroofs )  THEN
                            raddim = 3 * nz(iurb,id,j,i) + 2
                         ELSE
                            raddim = 2 * nz(iurb,id,j,i) + 1
                         ENDIF
!
!--                      Calculate the actual total radiation from the sky for every element.
                         bbb(:) = short_rad_inhom( i, j, id, iurb, raddimfull, raddim,             &
                                                   rsdown1dww(1:nz(iurb,id,j,i),id),               &
                                                   rsdown1dwe(1:nz(iurb,id,j,i),id),               &
                                                   rsdown1dr(1:nz(iurb,id,j,i)+1,id),              &
                                                   rsdown1dg(id), rsdiff1dww(1:nz(iurb,id,j,i),id),&
                                                   rsdiff1dwe(1:nz(iurb,id,j,i),id),               &
                                                   rsdiff1dr(1:nz(iurb,id,j,i)+1,id),              &
                                                   rsdiff1dg(id), radcor(iurb,id,j,i) )

                         bbb(1:raddim) = MATMUL( rsmatrix(iurb,id,1:raddim,1:raddim,j,i),          &
                                                 bbb(1:raddim) )

                         DO  iu = 1, nz(iurb,id,j,i)
                            rs_wall(iurb,2*id-1,iu,j,i) = bbb(iu)
                         ENDDO

                         DO  iu = nz(iurb,id,j,i)+1, 2*nz(iurb,id,j,i)
                            rs_wall(iurb,2*id,iu-nz(iurb,id,j,i),j,i) = bbb(iu)
                         ENDDO

                         rs_ground(iurb,id,j,i) = bbb(2*nz(iurb,id,j,i)+1)

                         IF ( lrroofs )  THEN
                            DO  iu = 2*nz(iurb,id,j,i)+2, 3*nz(iurb,id,j,i)+2
                               rs_roof(iurb,id,iu-(2*nz(iurb,id,j,i)+1),j,i) = bbb(iu)
                            ENDDO
                         ELSE
                            DO  iu = 1, nz(iurb,id,j,i)+1
                               rs_roof(iurb,id,iu,j,i) = rs(j,i) + rsdiff(j,i)
                            ENDDO
                         ENDIF
                      ENDIF
                   ENDDO

                ENDIF
             ENDDO

          ENDIF
       ENDDO
    ENDDO

!
!-- How much radiation is emitted into the sky?
    rsup = uprad( rs+rsdiff, rs_ground, rs_wall, rs_roof, alb_ground(:),                           &
                  alb_wall(:), alb_roof(:), .FALSE. )

    DO  i = nxl, nxr
       DO  j = nys, nyn

          IF ( fr_urb(j,i) > eps_urb  .AND.  rs(j,i) > eps_urb )  THEN
!
!--          Get effective albedo for parallel radiation.
             albedop_urb(j,i) = ( rsup(j,i) - albedo_urb(j,i) * rsdiff(j,i) ) / rs(j,i)
          ENDIF

       ENDDO
    ENDDO

    DEALLOCATE( bbb )

 END SUBROUTINE dcep_modify_sw


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Inhomogenities of linear system of equations for longwave
!> radiation: includes radiation from the sky and heat radiation from other surfaces.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION long_rad_inhom( i, j, id, iurb, ndimfull, ndim, rlww, rlwe, rlr, rlg, rc )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) ::  i
    INTEGER(iwp), INTENT(in) ::  id
    INTEGER(iwp), INTENT(in) ::  iurb
    INTEGER(iwp), INTENT(in) ::  j
    INTEGER(iwp), INTENT(in) ::  ndim
    INTEGER(iwp), INTENT(in) ::  ndimfull

    REAL(wp), INTENT(in) ::  rc          !< west east roof radiation
    REAL(wp), INTENT(in) ::  rlg         !< west east roof radiation

    REAL(wp), DIMENSION(nz(iurb,id,j,i)+1), INTENT(in) ::  rlr   !< incoming
    REAL(wp), DIMENSION(nz(iurb,id,j,i)),   INTENT(in) ::  rlwe  !< incoming
    REAL(wp), DIMENSION(nz(iurb,id,j,i)),   INTENT(in) ::  rlww  !< incoming

    REAL(wp) ::  long_rad_inhom(ndimfull)  !<

!
!-- Local variables
    INTEGER(iwp) ::  iu  !<
    INTEGER(iwp) ::  ju  !<
    INTEGER(iwp) ::  ou  !<


!
!-- West wall
    DO  iu = 1, nz(iurb,id,j,i)
!
!--    1) Incoming east wall + ground radiation due to temperature
       long_rad_inhom(iu) = emiss_ground(iurb) * fgw(iurb,id,iu,j,i) * sigma_sb *                  &
                            t_ground(iurb,id,ke_ground,j,i)**4
       DO  ju = 1, nz(iurb,id,j,i)
          long_rad_inhom(iu) = long_rad_inhom(iu) + fr_wall(iurb,id,ju,j,i) * emiss_wall(iurb) *   &
                                                    sigma_sb * fww(iurb,id,ju,iu,j,i) *            &
                                                    t_wall(iurb,2*id,ke_wall,ju,j,i)**4
       ENDDO

       IF ( lrroofs )  THEN
!
!--       2) from roof
          DO  ju = 1, iu
             long_rad_inhom(iu) = long_rad_inhom(iu) + fr_roof(iurb,id,ju,j,i) * emiss_roof(iurb) *&
                                                       sigma_sb * frw(iurb,id,ju,iu,j,i) *         &
                                                       t_roof(iurb,id,ke_roof,ju,j,i)**4
          ENDDO
!
!--       3) from other canyon
          DO  ou = 1, nz(iurb,id,j,i)+1
!
!--          From sky and other ground.
             long_rad_inhom(iu) = long_rad_inhom(iu) + fsow(iurb,id,ou,iu,j,i) * rlww(iu) *        &
                                                       fr_roof(iurb,id,ou,j,i)                     &
                                                     + emiss_ground(iurb) * fgow(iurb,id,ou,iu,j,i)&
                                                       * sigma_sb                                  &
                                                       * t_ground(iurb,id,ke_ground,j,i)**4        &
                                                       * fr_roof(iurb,id,ou,j,i)
!
!--          From east wall other canyon and sky rad from side.
             DO  ju = 1, nz(iurb,id,j,i)
                long_rad_inhom(iu) = long_rad_inhom(iu) + fr_wall(iurb,id,ju,j,i) *                &
                                                          emiss_wall(iurb) * sigma_sb              &
                                                          * fwow(iurb,id,ju,ou,iu,j,i) *           &
                                                          t_wall(iurb,2*id,ke_wall,ju,j,i)**4 *    &
                                                          fr_roof(iurb,id,ou,j,i) &
                                                        + fwow(iurb,id,ju,ou,iu,j,i) * rlwe(ju) *  &
                                                          rc * (1.0_wp - fr_wall(iurb,id,ju,j,i))  &
                                                          * fr_roof(iurb,id,ou,j,i)
             ENDDO
          ENDDO
       ELSE
!
!--       4) from sky
          long_rad_inhom(iu) = long_rad_inhom(iu) + fsw(iurb,id,iu,j,i) * rlww(iu)
!
!--       Sky rad from side.
          DO  ju = 1, nz(iurb,id,j,i)
             long_rad_inhom(iu) = long_rad_inhom(iu) + fww(iurb,id,ju,iu,j,i) * rlwe(ju) * rc *    &
                                                       ( 1.0_wp - fr_wall(iurb,id,ju,j,i) )
          ENDDO
       ENDIF

    ENDDO
!
!-- East wall
    DO  iu = 1+nz(iurb,id,j,i), 2*nz(iurb,id,j,i)
!
!--    1) incoming west wall + ground radiation
       long_rad_inhom(iu) = emiss_ground(iurb) * fgw(iurb,id,iu-nz(iurb,id,j,i),j,i) * sigma_sb *  &
                            t_ground(iurb,id,ke_ground,j,i)**4
       DO  ju = 1, nz(iurb,id,j,i)
          long_rad_inhom(iu) = long_rad_inhom(iu) + fr_wall(iurb,id,ju,j,i) * emiss_wall(iurb) *   &
                                                    sigma_sb *                                     &
                                                    fww(iurb,id,ju,iu-nz(iurb,id,j,i),j,i) *       &
                                                    t_wall(iurb,2*id-1,ke_wall,ju,j,i)**4
       ENDDO

       IF ( lrroofs )  THEN
!
!--       2) from roofs
          DO  ju = 1, iu-nz(iurb,id,j,i)
             long_rad_inhom(iu) = long_rad_inhom(iu) + fr_roof(iurb,id,ju,j,i) * emiss_roof(iurb) * &
                                                       sigma_sb *                                   &
                                                       frw(iurb,id,ju,iu-nz(iurb,id,j,i),j,i) *     &
                                                       t_roof(iurb,id,ke_roof,ju,j,i)**4
          ENDDO
!
!--       3) from other canyon
          DO  ou = 1, nz(iurb,id,j,i)+1
!
!--          From sky and other ground.
             long_rad_inhom(iu) = long_rad_inhom(iu) + fsow(iurb,id,ou,iu-nz(iurb,id,j,i),j,i) *   &
                                                       rlwe(iu-nz(iurb,id,j,i)) *                  &
                                                       fr_roof(iurb,id,ou,j,i)                     &
                                                     + emiss_ground(iurb) *                        &
                                                       fgow(iurb,id,ou,iu-nz(iurb,id,j,i),j,i) *   &
                                                       sigma_sb *                                  &
                                                       t_ground(iurb,id,ke_ground,j,i)**4 *        &
                                                       fr_roof(iurb,id,ou,j,i)
             DO  ju = 1, nz(iurb,id,j,i)
!
!--             From west wall other canyon and rad from side.
                long_rad_inhom(iu) = long_rad_inhom(iu) + fr_wall(iurb,id,ju,j,i) *                &
                                                         emiss_wall(iurb) * sigma_sb *             &
                                                         fwow(iurb,id,ju,ou,iu-nz(iurb,id,j,i),j,i)&
                                                         * t_wall(iurb,2*id-1,ke_wall,ju,j,i)**4   &
                                                         * fr_roof(iurb,id,ou,j,i)                 &
                                                       + fwow(iurb,id,ju,ou,iu-nz(iurb,id,j,i),j,i)&
                                                         * rlwe(ju) * rc * ( 1.0_wp -              &
                                                                         fr_wall(iurb,id,ju,j,i) ) &
                                                         * fr_roof(iurb,id,ou,j,i)
             ENDDO
          ENDDO
       ELSE
!
!--       3) from sky
          long_rad_inhom(iu) = long_rad_inhom(iu) + fsw(iurb,id,iu-nz(iurb,id,j,i),j,i) *          &
                                                    rlwe(iu-nz(iurb,id,j,i))
!
!--       4) rad from side
          DO  ju = 1, nz(iurb,id,j,i)
             long_rad_inhom(iu) = long_rad_inhom(iu) + fww(iurb,id,ju,iu-nz(iurb,id,j,i),j,i) *    &
                                                       rlww(ju) * rc *                             &
                                                       ( 1.0_wp - fr_wall(iurb,id,ju,j,i) )
          ENDDO
       ENDIF

    ENDDO
!
!-- Ground
    long_rad_inhom(2*nz(iurb,id,j,i)+1) = 0.0_wp
!
!-- Thermal radiation walls.
    DO  iu = 1, nz(iurb,id,j,i)
       long_rad_inhom(2*nz(iurb,id,j,i)+1) = long_rad_inhom(2*nz(iurb,id,j,i)+1) +                 &
                                             emiss_wall(iurb) * sigma_sb * fwg(iurb,id,iu,j,i) *   &
                                             fr_wall(iurb,id,iu,j,i) * (                           &
                                                             t_wall(iurb,2*id,ke_wall,iu,j,i)**4 + &
                                                             t_wall(iurb,2*id-1,ke_wall,iu,j,i)**4 &
                                                                       )
    ENDDO

    IF ( lrroofs )  THEN
!
!--    Radiation from nearer side.
       DO  iu = 1, nz(iurb,id,j,i)
          long_rad_inhom(2*nz(iurb,id,j,i)+1) = long_rad_inhom(2*nz(iurb,id,j,i)+1) +              &
                                                fwg(iurb,id,iu,j,i) *                              &
                                                ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) * rlg * rc
       ENDDO
!
!--    Radiation from other canyon.
       DO  ou = 1, nz(iurb,id,j,i)+1
!
!--       From sky
          long_rad_inhom(2*nz(iurb,id,j,i)+1) = long_rad_inhom(2*nz(iurb,id,j,i)+1) +              &
                                                fsog(iurb,id,ou,j,i) * rlg * fr_roof(iurb,id,ou,j,i)
!
!--       West and east wall other canyon and sky rad from farer side.
          DO  iu = 1, nz(iurb,id,j,i)
             long_rad_inhom(2*nz(iurb,id,j,i)+1) = long_rad_inhom(2*nz(iurb,id,j,i)+1) +           &
                                                   emiss_wall(iurb) * sigma_sb *                   &
                                                   fwog(iurb,id,iu,ou,j,i) *                       &
                                                   fr_wall(iurb,id,iu,j,i) * (                     &
                                                           t_wall(iurb,2*id  ,ke_wall,iu,j,i)**4 + &
                                                           t_wall(iurb,2*id-1,ke_wall,iu,j,i)**4   &
                                                                             )                     &
                                                   * fr_roof(iurb,id,ou,j,i) +                     &
                                                   fwog(iurb,id,iu,ou,j,i) *                       &
                                                   ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) * rlg * rc &
                                                   * fr_roof(iurb,id,ou,j,i)
          ENDDO
       ENDDO
    ELSE
!
!--    From sky.
       long_rad_inhom(2*nz(iurb,id,j,i)+1) = long_rad_inhom(2*nz(iurb,id,j,i)+1) +                 &
                                             fsg(iurb,id,j,i) * rlg
!--    Twice from side.
       DO  iu = 1, nz(iurb,id,j,i)
          long_rad_inhom(2*nz(iurb,id,j,i)+1) = long_rad_inhom(2*nz(iurb,id,j,i)+1) +              &
                                                2.0_wp * fwg(iurb,id,iu,j,i) *                     &
                                                ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) * rlg * rc
      ENDDO
    ENDIF
!
!-- Roofs.
    IF ( lrroofs )  THEN
       DO  iu = 2*nz(iurb,id,j,i)+2, ndim
!
!--       Radiation from top urban canopy.
          long_rad_inhom(iu) = fsr(iurb,id,iu-(2*nz(iurb,id,j,i)+1),j,i) *                         &
                               rlr(iu-(2*nz(iurb,id,j,i)+1))
!
!--       From walls and radiation where no walls.
          DO  ju = iu-2*nz(iurb,id,j,i)-1, nz(iurb,id,j,i)
             long_rad_inhom(iu) = long_rad_inhom(iu) + emiss_wall(iurb) * sigma_sb *               &
                                                       fwr(iurb,id,ju,iu-2*nz(iurb,id,j,i)-1,j,i) *&
                                                       fr_wall(iurb,id,ju,j,i) * (                 &
                                                       t_wall(iurb,2*id,ke_wall,ju,j,i)**4 + &
                                                       t_wall(iurb,2*id-1,ke_wall,ju,j,i)**4 )     &
                                                     + 2.0_wp *                                    &
                                                       fwr(iurb,id,ju,iu-2*nz(iurb,id,j,i)-1,j,i) *&
                                                       ( 1.0_wp - fr_wall(iurb,id,ju,j,i) ) *      &
                                                       rlr(iu-(2*nz(iurb,id,j,i)+1)) * rc
          ENDDO
       ENDDO
    ENDIF

 END FUNCTION long_rad_inhom


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &dcep_par for new modules
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_parin

    IMPLICIT NONE

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /dcep_parameters/  albedo_ground,                                                     &
                                albedo_roof,                                                       &
                                albedo_wall,                                                       &
                                dcep_average_radiation,                                            &
                                dzlayer_ground,                                                    &
                                dzlayer_roof,                                                      &
                                dzlayer_wall,                                                      &
                                emissivity_ground,                                                 &
                                emissivity_roof,                                                   &
                                emissivity_wall,                                                   &
                                iurb_cdrag,                                                        &
                                iurb_ls,                                                           &
                                ke_ground,                                                         &
                                ke_roof,                                                           &
                                ke_uhl,                                                            &
                                ke_wall,                                                           &
                                limpufl,                                                           &
                                lrroofs,                                                           &
                                ltintfix,                                                          &
                                lurbradcor,                                                        &
                                lurbvel,                                                           &
                                n_uclass,                                                          &
                                n_udir,                                                            &
                                rlength_ground,                                                    &
                                rlength_roof,                                                      &
                                thermdiff_ground,                                                  &
                                thermdiff_roof,                                                    &
                                thermdiff_wall,                                                    &
                                tinterior_ground,                                                  &
                                tinterior_roof,                                                    &
                                tinterior_wall,                                                    &
                                z_uhl

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, dcep_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    dcep_parameters namelist was found and read correctly. Set flag that indicates that the
!--    dcep model is switched on.
       IF ( .NOT. switch_off_module )  dcep = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    dcep_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'dcep_parameters', line )

    ENDIF

 END SUBROUTINE dcep_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Allocate DCEP fields.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_init_arrays

    IMPLICIT NONE

!
!-- Urban characteristics
    ALLOCATE( fr_uclass(n_uclass,nys:nyn,nxl:nxr) )
    ALLOCATE( fr_udir(n_uclass,n_udir,nys:nyn,nxl:nxr) )
    ALLOCATE( fr_urb(nys:nyn,nxl:nxr) )
    ALLOCATE( fr_rur(nys:nyn,nxl:nxr) )
!
!-- Radiative characteristics of urban domain
    ALLOCATE( albedop_urb(nys:nyn,nxl:nxr) )
    ALLOCATE( emiss_urb(nys:nyn,nxl:nxr) )
    ALLOCATE( t_grad_urb(nys:nyn,nxl:nxr) )
!
!-- Street and urban street canyon characteristics
    ALLOCATE( width_street(n_uclass,n_udir,nys:nyn,nxl:nxr) )   ! Street width [m]
    ALLOCATE( width_build(n_uclass,n_udir,nys:nyn,nxl:nxr) )    ! Building width [m]
    ALLOCATE( width_sglcan(n_uclass,n_udir,nys:nyn,nxl:nxr) )   ! Street width [m]
    ALLOCATE( width_dblcan(n_uclass,n_udir,nys:nyn,nxl:nxr) )   ! Building width [m]
    ALLOCATE( fr_street(n_uclass,n_udir,nys:nyn,nxl:nxr) )      ! Street width [m]
    ALLOCATE( fr_build(n_uclass,n_udir,nys:nyn,nxl:nxr) )       ! Building width [m]
    ALLOCATE( angrotx_udir(n_udir,nys:nyn,nxl:nxr) )            ! Street direction, x component
    ALLOCATE( angroty_udir(n_udir,nys:nyn,nxl:nxr) )            ! Street direction, y component
!
!-- Probability that a building has a height equal to z, corresponds to roof heigths
    ALLOCATE( fr_roof(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )
!
!-- Probability that a building has a height greater or equal to z, corresponds to roof heigths
    ALLOCATE( fr_wall(n_uclass,n_udir,ke_uhl,nys:nyn,nxl:nxr) )
!
!-- Height of wall element, value for ke_urb+1 is for air above hightest roof level
    ALLOCATE( dz_uhl(ke_uhl+1) )

    ALLOCATE( fww(n_uclass,n_udir,ke_uhl,ke_uhl,nys:nyn,nxl:nxr) )   !  from wall to wall
    ALLOCATE( fwg(n_uclass,n_udir,ke_uhl,nys:nyn,nxl:nxr) )          !  from wall to ground
    ALLOCATE( fgw(n_uclass,n_udir,ke_uhl,nys:nyn,nxl:nxr) )          !  from ground to wall
    ALLOCATE( fsw(n_uclass,n_udir,ke_uhl,nys:nyn,nxl:nxr) )          !  from sky to wall
    ALLOCATE( fws(n_uclass,n_udir,ke_uhl,nys:nyn,nxl:nxr) )          !  from wall to sky
    ALLOCATE( fsg(n_uclass,n_udir,nys:nyn,nxl:nxr) )
    ALLOCATE( fgs(n_uclass,n_udir,nys:nyn,nxl:nxr) )

    ALLOCATE( umfl_ground(n_uclass,n_udir,nys:nyn,nxl:nxr) )
    ALLOCATE( umfl_roof(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )
    ALLOCATE( umfl_wall(n_uclass,2*n_udir,ke_uhl,nys:nyn,nxl:nxr) )
    ALLOCATE( vmfl_ground(n_uclass,n_udir,nys:nyn,nxl:nxr) )
    ALLOCATE( vmfl_roof(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )
    ALLOCATE( vmfl_wall(n_uclass,2*n_udir,ke_uhl,nys:nyn,nxl:nxr) )

    IF ( lrroofs )  THEN
       ALLOCATE( frw(n_uclass,n_udir,ke_uhl+1,ke_uhl,nys:nyn,nxl:nxr) )
       ALLOCATE( fwr(n_uclass,n_udir,ke_uhl,ke_uhl+1,nys:nyn,nxl:nxr) )
       ALLOCATE( fsr(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )
       ALLOCATE( frs(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )

       ALLOCATE( fwow(n_uclass,n_udir,ke_uhl,ke_uhl+1,ke_uhl,nys:nyn,nxl:nxr) )
       ALLOCATE( fgos(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )
       ALLOCATE( fsog(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )
       ALLOCATE( fgow(n_uclass,n_udir,ke_uhl+1,ke_uhl,nys:nyn,nxl:nxr) )
       ALLOCATE( fwog(n_uclass,n_udir,ke_uhl,ke_uhl+1,nys:nyn,nxl:nxr) )
       ALLOCATE( fsow(n_uclass,n_udir,ke_uhl+1,ke_uhl,nys:nyn,nxl:nxr) )
       ALLOCATE( fwos(n_uclass,n_udir,ke_uhl,ke_uhl+1,nys:nyn,nxl:nxr) )
    ENDIF

    ALLOCATE( z_midu(ke_uhl+2,n_uclass) )

    ALLOCATE( nz_meso(nys:nyn,nxl:nxr) )
    ALLOCATE( nz(n_uclass,n_udir,nys:nyn,nxl:nxr) )
!
!-- 2 levels more then nz_u, +1 for the highest roof which is on nz_u+1, +1 for calulation of
!-- Richardson number.
!-- Temperature in each layer of the wall:
    ALLOCATE( t_wall(n_uclass,2*n_udir,ke_wall,ke_uhl,nys:nyn,nxl:nxr) )
!
!-- Temperature in each layer of the ground:
    ALLOCATE( t_ground(n_uclass,n_udir,ke_ground,nys:nyn,nxl:nxr) )
!
!-- Temperature in each layer of the roof:
    ALLOCATE( t_roof(n_uclass,n_udir,ke_roof,ke_uhl+1,nys:nyn,nxl:nxr) )

    ALLOCATE( rl_ground(n_uclass,n_udir,nys:nyn,nxl:nxr) )
    ALLOCATE( rl_wall(n_uclass,2*n_udir,ke_uhl,nys:nyn,nxl:nxr) )
    ALLOCATE( rl_roof(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )

    ALLOCATE( rs_ground(n_uclass,n_udir,nys:nyn,nxl:nxr) )
    ALLOCATE( rs_wall(n_uclass,2*n_udir,ke_uhl,nys:nyn,nxl:nxr) )
    ALLOCATE( rs_roof(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )

    IF ( lrroofs )  THEN
!
!--   Radiation matrix: west wall, east wall, ground, roofs for every site, street direction,
!--   and urban class.
      ALLOCATE( rlmatrix(n_uclass,n_udir,2*ke_uhl+1+(ke_uhl+1),                                    &
                         2*ke_uhl+1+(ke_uhl+1),nys:nyn,nxl:nxr) )
      ALLOCATE( rsmatrix(n_uclass,n_udir,2*ke_uhl+1+(ke_uhl+1),                                    &
                         2*ke_uhl+1+(ke_uhl+1),nys:nyn,nxl:nxr) )
      ALLOCATE( rlipiv(n_uclass,n_udir,2*ke_uhl+1+(ke_uhl+1),nys:nyn,nxl:nxr) )
      ALLOCATE( rsipiv(n_uclass,n_udir,2*ke_uhl+1+(ke_uhl+1),nys:nyn,nxl:nxr) )
    ELSE
!
!--   Radiation matrix: west wall, east wall, ground for every site, street direction,
!--   and urban class.
      ALLOCATE( rlmatrix(n_uclass,n_udir,2*ke_uhl+1,2*ke_uhl+1,nys:nyn,nxl:nxr) )
      ALLOCATE( rsmatrix(n_uclass,n_udir,2*ke_uhl+1,2*ke_uhl+1,nys:nyn,nxl:nxr) )
      ALLOCATE( rlipiv(n_uclass,n_udir,2*ke_uhl+1,nys:nyn,nxl:nxr) )
      ALLOCATE( rsipiv(n_uclass,n_udir,2*ke_uhl+1,nys:nyn,nxl:nxr) )
    ENDIF

    ALLOCATE( radcor(n_uclass,n_udir,nys:nyn,nxl:nxr) )

    ALLOCATE( strfl_ground(n_uclass,n_udir,nys:nyn,nxl:nxr) )
    ALLOCATE( strfl_roof(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )
    ALLOCATE( strfl_wall(n_uclass,2*n_udir,ke_uhl,nys:nyn,nxl:nxr) )

    ALLOCATE( shfl_ground(n_uclass,n_udir,nys:nyn,nxl:nxr) )
    ALLOCATE( shfl_roof(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr) )
    ALLOCATE( shfl_wall(n_uclass,2*n_udir,ke_uhl,nys:nyn,nxl:nxr) )

    ALLOCATE( tddz_roof(ke_roof-1,n_uclass) )
    ALLOCATE( tddz_wall(ke_wall-1,n_uclass) )
    ALLOCATE( tddz_ground(ke_ground-1,n_uclass) )

    ALLOCATE( cdrag_wall(n_uclass,n_udir,nys:nyn,nxl:nxr) )
!
!-- Allocation of other arrays (were not originally in alloc_urban subroutine in DCEP source
!-- but in other locations).
    ALLOCATE( dz_ground(ke_ground,n_uclass) )
    ALLOCATE( dz_roof(ke_roof,n_uclass) )
    ALLOCATE( dz_wall(ke_wall,n_uclass) )

    ALLOCATE( z_ground(1:ke_ground) )
    ALLOCATE( z_roof(1:ke_roof) )
    ALLOCATE( z_wall(1:ke_wall) )

    ALLOCATE( alb_ground(n_uclass) )
    ALLOCATE( alb_wall(n_uclass) )
    ALLOCATE( alb_roof(n_uclass) )

    ALLOCATE( albedo_urb(nys:nyn,nxl:nxr) )
    ALLOCATE( z0_urb(nys:nyn,nxl:nxr) )

    ALLOCATE( emiss_ground(n_uclass) )
    ALLOCATE( emiss_wall(n_uclass) )
    ALLOCATE( emiss_roof(n_uclass) )

    ALLOCATE( strfl_urb(nys:nyn,nxl:nxr) )
    ALLOCATE( shfl_urb(nys:nyn,nxl:nxr) )
    ALLOCATE( umfl_urb(nys:nyn,nxl:nxr) )
    ALLOCATE( vmfl_urb(nys:nyn,nxl:nxr) )
    ALLOCATE( t_g_urb(nys:nyn,nxl:nxr) )

    ALLOCATE( z0_ground(n_uclass) )
    ALLOCATE( z0_roof(n_uclass) )

    ALLOCATE( cs_ground(n_uclass) )
    ALLOCATE( cs_wall(n_uclass) )
    ALLOCATE( cs_roof(n_uclass) )

    ALLOCATE( volair_urb(nzb:nzt+1,nys:nyn,nxl:nxr) )
    ALLOCATE( volairhl_urb(nzb:nzt+1,nys:nyn,nxl:nxr) )
    ALLOCATE( volairhlred_urb(nzb:nzt+1,nys:nyn,nxl:nxr) )

    ALLOCATE( tint_ground(n_uclass) )
    ALLOCATE( tint_roof(n_uclass) )
    ALLOCATE( tint_wall(n_uclass) )

    ALLOCATE( l_urb(nzb:nzt+1,nys:nyn,nxl:nxr) )
    ALLOCATE( zeff_urb(nzb:nzt+1,nys:nyn,nxl:nxr) )

    ALLOCATE( td_ground(ke_ground,n_uclass) )
    ALLOCATE( td_roof(ke_roof,n_uclass) )
    ALLOCATE( td_wall(ke_wall,n_uclass) )


 END SUBROUTINE dcep_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> dcep_init initiates all fields related to DCEP model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_init


    IMPLICIT NONE

    INTEGER(iwp) ::  id      !<
    INTEGER(iwp) ::  i       !< loop running index
    INTEGER(iwp) ::  imax    !<
    INTEGER(iwp) ::  iu      !<
    INTEGER(iwp) ::  iurb    !<
    INTEGER(iwp) ::  iz      !<
    INTEGER(iwp) ::  j       !<
    INTEGER(iwp) ::  ju      !<
    INTEGER(iwp) ::  k       !<
    INTEGER(iwp) ::  k_topo  !< topography top index
    INTEGER(iwp) ::  m       !<
    INTEGER(iwp) ::  maxz    !<
    INTEGER(iwp) ::  ou      !<
    INTEGER(iwp) ::  raddim  !< dimension of radiation matrices

    REAL(wp) ::  b_avg       !<
    REAL(wp) ::  d           !<
    REAL(wp) ::  dh          !<
    REAL(wp) ::  dlgtmp      !<
    REAL(wp) ::  dzz         !< average grid size in longitude and latitude
    REAL(wp) ::  lc          !<
    REAL(wp) ::  sftot       !<
    REAL(wp) ::  ssl         !<
    REAL(wp) ::  temp        !<
    REAL(wp) ::  temp2       !<
    REAL(wp) ::  ulu         !<
    REAL(wp) ::  vtot2       !< average grid size in longitude and latitude
    REAL(wp) ::  vtot        !< average grid size in longitude and latitude
    REAL(wp) ::  xang        !<
    REAL(wp) ::  yang        !<
    REAL(wp) ::  z0_avg      !<
    REAL(wp) ::  z1          !< average grid size in longitude and latitude
    REAL(wp) ::  z2          !< average grid size in longitude and latitude

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::  eye        !<
    REAL(wp), ALLOCATABLE, DIMENSION(:)   ::  one        !<
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::  tmpmatrix  !<
    REAL(wp), ALLOCATABLE, DIMENSION(:)   ::  vec        !<
    REAL(wp), ALLOCATABLE, DIMENSION(:)   ::  zero       !<

    IF ( debug_output )  CALL debug_message( 'dcep_init', 'start' )

!
!-- Input DCEP data via netCDF file.
    CALL dcep_netcdf_input
!
!-- Height of wall element, value for ke_urb+1 is for air above hightest roof level.
    dz_uhl = z_uhl(2:(ke_uhl+2)) - z_uhl(1:(ke_uhl+1))

    DO  iurb = 1, n_uclass
       z_midu(1,iurb) = 0.0_wp
       DO  iu = 1, ke_uhl+1
          z_midu(iu+1,iurb) = 0.5_wp * (z_uhl(iu+1) + z_uhl(iu))
       ENDDO
    ENDDO
!
!-- Adjust urban fraction (fr_urb = -9999. as a default) and calculate fr_rur.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          IF ( fr_urb(j,i) < eps_urb )  fr_urb(j,i) = 0.0_wp
          fr_rur(j,i) = 1.0_wp - fr_urb(j,i)
       ENDDO
    ENDDO
!
!-- Find mesoscale layers in urban layer, hhl is absolute height (not height above ground).
    maxz = MAXVAL( z_uhl )
    nz_mesom = -1

    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Starting from the bottom (above orography).
          k_topo = topo_top_ind(j,i,0)
          imax = 1

          DO WHILE ( ( zw(k_topo+imax) - zw(k_topo) ) < maxz )
             imax = imax + 1
          ENDDO

          nz_meso(j,i) = imax
       ENDDO
    ENDDO

    nz_mesom = MAXVAL(nz_meso)
!
!-- Save the z values into z_meso.
    ALLOCATE( dz_meso(1:nz_mesom,nys:nyn,nxl:nxr) )
    ALLOCATE( z_meso(1:nz_mesom+1,nys:nyn,nxl:nxr) )
    ALLOCATE( z_midmeso(1:nz_mesom+2,nys:nyn,nxl:nxr) )

    DO  i = nxl, nxr
       DO  j = nys, nyn
          k_topo = topo_top_ind(j,i,0)
          DO  iz = 1, nz_meso(j,i)
             dz_meso(iz,j,i) = dzw(iz+k_topo)
          ENDDO

          DO  iz = 1, nz_meso(j,i)+1
             z_meso(iz,j,i) = zw(iz+k_topo-1)
          ENDDO

          z_midmeso(1,j,i) = 0.0_wp
          DO  iz = 1, nz_meso(j,i)+1
             z_midmeso(iz+1,j,i) = z_meso(iz,j,i) + 0.5_wp * ( zw(iz+k_topo+1) - zw(iz+k_topo) )
          ENDDO
       ENDDO
    ENDDO
!
!-- Confirm that fields are correctly normalized to 1 and set minimum building width to 1.5.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          IF ( fr_urb(j,i) > eps_urb )  THEN
             DO  iurb = 1, n_uclass
                temp = SUM( fr_udir(iurb,:,j,i) )
                fr_udir(iurb,:,j,i) = fr_udir(iurb,:,j,i) / temp

                DO  id = 1, n_udir
                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
                      IF ( width_build(iurb,id,j,i) > 1.5_wp )  THEN
                         temp = SUM( fr_roof(iurb,id,1:ke_uhl+1,j,i) )
                         fr_roof(iurb,id,1:ke_uhl+1,j,i) = fr_roof(iurb,id,1:ke_uhl+1,j,i) / temp
                      ELSE
!
!--                      Remove buildings with building width of 1.5m.
                         fr_roof(iurb,id,:,j,i) = 0.0_wp
                         fr_roof(iurb,id,1,j,i) = 1.0_wp
                      ENDIF
                   ENDIF
                ENDDO

             ENDDO
          ENDIF
       ENDDO
    ENDDO
!
!-  Number of wall levels (nz).
    DO  i = nxl,nxr
       DO  j = nys, nyn
          DO  id = 1, n_udir
             DO  iurb = 1, n_uclass

                IF ( fr_urb(j,i) > eps_urb  .AND.  fr_udir(iurb,id,j,i) > eps_urb )  THEN
                   maxz = ke_uhl
                   DO WHILE ( fr_roof(iurb,id,maxz+1,j,i) < eps_urb )
                      maxz = maxz - 1
                   ENDDO
                   nz(iurb,id,j,i) = maxz
                ELSE
                   nz(iurb,id,j,i) = 0
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculate the surface characteristics (albedo, emissivity, heat capacity, interior temp.,
!-- roughness length).
!-- Albedo:
    alb_ground(:) = albedo_ground(1:n_uclass)
    alb_roof(:)   = albedo_roof(1:n_uclass)
    alb_wall(:)   = albedo_wall(1:n_uclass)
!
!-- Emissivity:
    emiss_ground(:) = emissivity_ground(1:n_uclass)
    emiss_roof(:)   = emissivity_roof(1:n_uclass)
    emiss_wall(:)   = emissivity_wall(1:n_uclass)
!
!-- Heat capacity:
    cs_ground(:) = heatcap_ground(1:n_uclass)
    cs_roof(:)   = heatcap_roof(1:n_uclass)
    cs_wall(:)   = heatcap_wall(1:n_uclass)
!
!-- Interior temperature:
    tint_ground(:) = tinterior_ground(1:n_uclass)
    tint_roof(:)   = tinterior_roof(1:n_uclass)
    tint_wall(:)   = tinterior_wall(1:n_uclass)
!
!-- Roughness length:
    z0_ground(:) = rlength_ground(1:n_uclass)
    z0_roof(:)   = rlength_roof(1:n_uclass)
!
!-- Thermal diffusivity of wall's layers:
    DO id = 1, ke_ground
       td_ground(id,:) = thermdiff_ground(1:n_uclass)
    ENDDO
    DO id = 1, ke_roof
       td_roof(id,:)   = thermdiff_roof(1:n_uclass)
    ENDDO
    DO id = 1, ke_wall
       td_wall(id,:)   = thermdiff_wall(1:n_uclass)
    ENDDO
!
!-- Initlialize urban surface temperatures (Note: we use fixed values for all classes for now).
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  id = 1, n_udir
             DO  iurb = 1, n_uclass

                t_ground(iurb,id,:,j,i) = tint_ground(iurb)

                DO  iu = 1, nz(iurb,id,j,i)
                   t_wall(iurb,2*id-1,:,iu,j,i) = tint_wall(iurb)
                   t_wall(iurb,2*id  ,:,iu,j,i) = tint_wall(iurb)
                ENDDO

                DO  iu = 1, nz(iurb,id,j,i) + 1
                   t_roof(iurb,id,:,iu,j,i) = tint_roof(iurb)
                ENDDO

             ENDDO
          ENDDO
       ENDDO
    ENDDO

    t_grad_urb  = 300.0_wp
    albedop_urb = albedo_urb
!
!-- Calculate size of canyons and fraction of street and buildings.
    fr_street = 0.0_wp
    fr_build  = 0.0_wp
    DO  i = nxl,nxr
       DO  j = nys, nyn
          DO  id = 1, n_udir
             DO  iurb = 1, n_uclass

                width_sglcan(iurb,id,j,i) = width_street(iurb,id,j,i) + width_build(iurb,id,j,i)
                width_dblcan(iurb,id,j,i) = 2.0_wp * width_street(iurb,id,j,i) +                   &
                                            width_build(iurb,id,j,i)

                IF ( width_sglcan(iurb,id,j,i) > eps_urb )  THEN
                   fr_street(iurb,id,j,i) = width_street(iurb,id,j,i) / width_sglcan(iurb,id,j,i)
                   fr_build(iurb,id,j,i)  = width_build(iurb,id,j,i)  / width_sglcan(iurb,id,j,i)
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculation of wall probablity.
    fr_wall = fill_value
    DO  i = nxl,nxr
       DO  j = nys, nyn
          DO  id = 1, n_udir
             DO  iurb = 1, n_uclass

                IF ( fr_urb(j,i) > eps_urb  .AND.  fr_udir(iurb,id,j,i) > eps_urb )  THEN
                   fr_wall(iurb,id,1,j,i) = 1.0_wp - fr_roof(iurb,id,1,j,i)
                   DO  iu = 2, nz(iurb,id,j,i)
                      fr_wall(iurb,id,iu,j,i) = MAX( 0.0_wp, fr_wall(iurb,id,iu-1,j,i) -           &
                                                             fr_roof(iurb,id,iu,j,i) )
                   ENDDO
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

!
!-- Thickness of ground, wall, roof layers.
    IF ( ANY( dzlayer_ground(1:ke_ground) < 0.0_wp ) )  THEN
       dzlayer_ground(1:10) = (/  0.2_wp, 0.12_wp, 0.08_wp,  0.05_wp, 0.03_wp,                     &
                                 0.02_wp, 0.02_wp, 0.01_wp, 0.005_wp, 0.0025_wp /)
    ENDIF

    IF ( ANY( dzlayer_roof(1:ke_roof) < 0.0_wp ) )  THEN
       dzlayer_roof(1:10) = (/ 0.02_wp, 0.02_wp, 0.02_wp, 0.02_wp,  0.02_wp,                        &
                               0.02_wp, 0.02_wp, 0.01_wp, 0.005_wp, 0.0025_wp /)
    ENDIF

    IF ( ANY( dzlayer_wall(1:ke_wall) < 0.0_wp ) )  THEN
       dzlayer_wall(1:10) = (/ 0.02_wp, 0.02_wp, 0.02_wp, 0.02_wp,  0.02_wp,                        &
                               0.02_wp, 0.02_wp, 0.01_wp, 0.005_wp, 0.0025_wp /)
    ENDIF

    DO  iurb = 1, n_uclass
       dz_ground(1:ke_ground,iurb) = dzlayer_ground(1:ke_ground)
       dz_roof(1:ke_roof,iurb) = dzlayer_roof(1:ke_roof)
       dz_wall(1:ke_wall,iurb) = dzlayer_wall(1:ke_wall)
    ENDDO

    DO  iurb = 1, n_uclass

       z_ground(1) = - dz_ground(1,iurb)
       z_roof(1)   = - dz_roof(1,iurb)
       z_wall(1)   = - dz_wall(1,iurb)

       DO  iz = 2, ke_ground
          z_ground(iz) = z_ground(iz-1) - dz_ground(iz,iurb)
       ENDDO

       DO  iz = 2, ke_roof
          z_roof(iz) = z_roof(iz-1) - dz_roof(iz,iurb)
       ENDDO

       DO  iz = 2, ke_wall
          z_wall(iz) = z_wall(iz-1) - dz_wall(iz,iurb)
       ENDDO

    ENDDO
!
!-- Mean of ground roughness lenghts.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          z0_urb(j,i) = 0.0_wp
          DO  iurb = 1, n_uclass
             z0_urb(j,i) = z0_urb(j,i) + z0_ground(iurb) * fr_uclass(iurb,j,i)
          ENDDO
       ENDDO
    ENDDO
!
!-- Update the roughness length.
    DO  m = 1, surf_lsm%ns
       i = surf_lsm%i(m)
       j = surf_lsm%j(m)
       surf_lsm%z0(m) = surf_lsm%z0(m) * fr_rur(j,i) + z0_urb(j,i) * fr_urb(j,i)
    ENDDO
!
!-- Calculate urban length scales.
!-- Note: they are not used yet, however, we calculated them for future development to consider
!--       the turbulence effect.
    zeff_urb = 0.0_wp
    l_urb    = 0.0_wp
    IF ( iurb_ls == 1 )  THEN
!
!--    Alberto's method for urban length scale.
!--    Calculation of the length scale taking into account the buildings effects.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             k_topo = topo_top_ind(j,i,0)
             DO  iz = k_topo+1, nz_meso(j,i)+1
                ulu = 0.0_wp
                ssl = 0.0_wp
                DO  id = 1, n_udir
                   DO  iurb = 1, n_uclass
                      DO  iu = 2, nz(iurb,id,j,i)+1
                         IF ( z_uhl(iu) > z_meso(iz,j,i) )  THEN
                            ulu = ulu + fr_roof(iurb,id,iu,j,i) * fr_udir(iurb,id,j,i) *           &
                                        fr_uclass(iurb,j,i) / z_uhl(iu)
                            ssl = ssl + fr_roof(iurb,id,iu,j,i) * fr_udir(iurb,id,j,i) *           &
                                        fr_uclass(iurb,j,i)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO

                IF ( ulu > azero )  THEN
                   l_urb(iz,j,i) = ssl / ulu
                ELSE
                   l_urb(iz,j,i) = 0.0_wp
                ENDIF
             ENDDO
             DO  iz = nz_meso(j,i)+2, nzt
                l_urb(iz,j,i) = 0.0_wp
             ENDDO
          ENDDO
       ENDDO
!
!--    Calculation of effective height over the surface.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             IF ( fr_urb(j,i) > eps_urb )  THEN
                k_topo = topo_top_ind(j,i,0)
                zeff_urb(k_topo,j,i) = 0.0_wp
                DO  iz = k_topo+1, nzt
                   sftot  = 0.0_wp
                   dlgtmp = 0.0_wp
                   DO  id = 1, n_udir
                      DO  iurb = 1, n_uclass
                         temp = zw(iz) - zw(k_topo)
                         sftot  = sftot  + width_street(iurb,id,j,i) * fr_udir(iurb,id,j,i) *      &
                                           fr_uclass(iurb,j,i)
                         dlgtmp = dlgtmp + width_street(iurb,id,j,i) * fr_udir(iurb,id,j,i) *      &
                                           fr_uclass(iurb,j,i) / temp
                         DO  iu = 1, nz(iurb,id,j,i) + 1
                            IF ( temp > z_uhl(iu) )  THEN
                               sftot = sftot + fr_roof(iurb,id,iu,j,i) * width_build(iurb,id,j,i)  &
                                               * fr_udir(iurb,id,j,i) * fr_uclass(iurb,j,i)
                               dlgtmp = dlgtmp + fr_roof(iurb,id,iu,j,i) * width_build(iurb,id,j,i)&
                                                 * fr_udir(iurb,id,j,i) * fr_uclass(iurb,j,i)      &
                                                 / ( temp - z_uhl(iu) )
                            ENDIF
                         ENDDO
                      ENDDO
                   ENDDO
                   zeff_urb(iz,j,i) = sftot / dlgtmp
                ENDDO
             ELSE
                zeff_urb(:,j,i) = 0.0_wp
             ENDIF
          ENDDO
       ENDDO

    ELSEIF ( iurb_ls == 2 )  THEN
!
!--    Coceal and Belcher 2004 method
!--    Calculate average building height and building fraction for displacement height calculation
       l_urb = 0.0_wp
       DO  i = nxl, nxr
          DO  j = nys, nyn
             IF ( fr_urb(j,i) > eps_urb )  THEN
                k_topo = topo_top_ind(j,i,0)

                DO  id = 1, n_udir
                   DO  iurb = 1, n_uclass
                      IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
                         b_avg = fr_build(iurb,id,j,i) * fr_urb(j,i)
                         dh = 1.0_wp + 4.0_wp**(-b_avg) * ( b_avg - 1)
                         z0_avg = fr_rur(j,i) * roughness_length / g + fr_urb(j,i) * z0_ground(iurb)

                         DO  iu = 1, nz(iurb,id,j,i)+1
                            lc = kappa * ( 1.0_wp - dh ) / dh * z_uhl(iu)
                            d  = dh * z_uhl(iu)
                            DO  iz = k_topo+1, nzt
                               temp = zw(iz) - zw(k_topo)
                               IF ( temp + z0_avg < z_uhl(iu) )  THEN
                                  l_urb(iz,j,i) = l_urb(iz,j,i) + fr_roof(iurb,id,iu,j,i) *        &
                                                                  fr_udir(iurb,id,j,i) *           &
                                                                  fr_uclass(iurb,j,i) / ( 1.0_wp / &
                                                                   ( temp + z0_avg ) + 1.0_wp / lc )
                               ELSE
                                  l_urb(iz,j,i) = l_urb(iz,j,i) + ( temp + z0_avg - d ) *          &
                                                                  fr_roof(iurb,id,iu,j,i) *        &
                                                                  fr_udir(iurb,id,j,i) *           &
                                                                  fr_uclass(iurb,j,i)
                               ENDIF
                            ENDDO
                         ENDDO

                      ENDIF
                   ENDDO
                ENDDO

             ENDIF
          ENDDO
       ENDDO

       zeff_urb = 0.0_wp

    ELSE ! IF ( iurb_ls == 1 )

       WRITE( message_string, * ) 'wrong value of type of urban length scale parametrization'
       CALL message( 'dcep_init', 'DCP0003', 1, 2, 0, 6, 0 )

    ENDIF
!
!-- Calculate drag coefficient for walls.
    IF ( iurb_cdrag == 1 )  THEN
!
!--    Use default from Martilli et al. (2002).
       cdrag_wall = 0.4_wp
    ELSEIF ( iurb_cdrag == 2 )  THEN
!
!--   Use values in Krayenhoff et al. (2015) depending of plan area density.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  id = 1, n_udir
                DO  iurb = 1, n_uclass
                   IF ( fr_urb(j,i) > eps_urb  .AND.  fr_udir(iurb,id,j,i) > eps_urb )  THEN
                      IF ( fr_build(iurb,id,j,i) > 0.33_wp )  THEN
                         cdrag_wall(iurb,id,j,i) = 3.67_wp
                      ELSE
                         cdrag_wall(iurb,id,j,i) = 7.30_wp * fr_build(iurb,id,j,i)**0.62_wp
                      ENDIF
                   ELSE
                      cdrag_wall(iurb,id,j,i) = 0.0_wp
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO

    ELSE

       WRITE( message_string, * ) 'wrong value of type of drag coefficient for walls (iurb_cdrag)'
       CALL message( 'dcep_init', 'DCP0004', 1, 2, 0, 6, 0 )

    ENDIF
!
!-- Cell volume correction factor.
    ALLOCATE( vol(n_uclass,n_udir,nz_mesom,nys:nyn,nxl:nxr) )
    volair_urb = 0.0_wp
    DO  i = nxl, nxr
       DO  j = nys, nyn

          IF ( fr_urb(j,i) > eps_urb )  THEN
             volair_urb(nzt-nz_mesom:nzt,j,i) = 1.0_wp
             k_topo = topo_top_ind(j,i,0)
             DO  iz = k_topo, nz_mesom-1
                volair_urb(iz,j,i) = 1.0_wp
                z1 = zw(iz  ) - zw(k_topo)
                z2 = zw(iz+1) - zw(k_topo)
                DO  id = 1, n_udir
                   DO  iurb = 1, n_uclass
                      IF ( fr_udir(iurb,id,j,i) > eps_urb)  THEN
                         vtot = 0.0_wp
                         DO  iu = 1, nz(iurb,id,j,i)
!
!--                         Height in the PALM cell.
                            dzz = MAX( MIN( z_uhl(iu+1), z2 ) - MAX( z_uhl(iu), z1 ), 0.0_wp )
                            vtot = vtot + fr_wall(iurb,id,iu,j,i) * dzz
                         ENDDO
                         vtot = vtot / ( zw(iz+1) - zw(iz) )
                         vol(iurb,id,iz+1,j,i) = 1.0_wp - vtot * fr_build(iurb,id,j,i)
                         volair_urb(iz+1,j,i) = volair_urb(iz+1,j,i) - vtot * fr_build(iurb,id,j,i)&
                                                                       * fr_udir(iurb,id,j,i) *    &
                                                                       fr_uclass(iurb,j,i)
                      ELSE
                         vol(iurb,id,iz+1,j,i) = 1.0_wp
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             volair_urb(:,j,i) = 1.0_wp
             vol(:,:,:,j,i) = 1.0_wp
          ENDIF

       ENDDO
    ENDDO
!
!-- Cell volume correction factor for volume between main levels.
!-- Please note that volume correction factor is not yet used but it will be used for TKE tendency.
    volairhl_urb               = 0.0_wp
    volairhlred_urb            = 0.0_wp
!
!-- There's no building in the top layer.
    volairhl_urb(nzt+1,:,:)    = 1.0_wp
    volairhlred_urb(nzt+1,:,:) = 1.0_wp

    ALLOCATE( volhl(n_uclass,n_udir,nz_mesom+1,nys:nyn,nxl:nxr) )

    DO  i = nxl, nxr
       DO  j = nys, nyn

          IF ( fr_urb(j,i) > eps_urb )  THEN

             volairhl_urb(nz_mesom+1:nzt+1,j,i) = 1.0_wp
             volairhlred_urb(nz_mesom+1:nzt+1,j,i) = 1.0_wp
             k_topo = topo_top_ind(j,i,0)
             z1 =  0.5_wp * ( zw(nz_mesom) + zw(nz_mesom+1) ) - zw(k_topo)
             DO  iz = k_topo, nz_mesom-1
                z2 = z1
                IF ( iz == k_topo )  THEN
                   z1 = 0.0_wp
                ELSE
                   z1 = 0.5_wp * ( zw(iz) + zw(iz+1) ) - zw(k_topo)
                ENDIF
                volairhl_urb(iz,j,i)    = 1.0_wp
                volairhlred_urb(iz,j,i) = 0.0_wp
                vtot2 = 0.0_wp
                DO  id = 1, n_udir
                   DO  iurb = 1, n_uclass
                      IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
                         vtot = 0.0_wp
                         DO  iu = 1, nz(iurb,id,j,i)
!
!--                         Height in the PALM cell.
                            dzz = MAX( MIN( z_uhl(iu+1),z2 ) - MAX( z_uhl(iu),z1 ), 0.0_wp )
                            vtot = vtot + fr_wall(iurb,id,iu,j,i) * dzz
!
!--                         Height in the PALM cell, MOST is assumed
!--                         only in half of the cell (similar to PALM lowest layer).
                            dzz = MAX( MIN( 0.5_wp * ( z_uhl(iu+1) + z_uhl(iu) ), z2 ) -           &
                                  MAX( z_uhl(iu), z1 ), 0.0_wp )
                            vtot2 = vtot2 + fr_roof(iurb,id,iu,j,i) * dzz
                         ENDDO
!
!--                      The layer over the uppermost roof only for calculation of tke flux volume
!--                      since no air volume reduction here; MOST is assumed only in half of the
!--                      cell (similar to PALM lowest layer).
                         dzz = MAX( MIN( 0.5_wp * ( z_uhl(nz(iurb,id,j,i)+2) +                     &
                                         z_uhl(nz(iurb,id,j,i)+1) ), z2 ) -                        &
                                    MAX( z_uhl(nz(iurb,id,j,i)+1), z1 ), 0.0_wp )
                         vtot2 = vtot2 + fr_roof(iurb,id,nz(iurb,id,j,i)+1,j,i) * dzz

                         vtot = vtot / (z2-z1)
                         vtot2 = vtot2 / (z2-z1)

                         volhl(iurb,id,iz-k_topo+1,j,i) = 1.0_wp - vtot * fr_build(iurb,id,j,i)

                         volairhl_urb(iz,j,i) = volairhl_urb(iz,j,i) - vtot *                      &
                                                                       fr_build(iurb,id,j,i) *     &
                                                                       fr_udir(iurb,id,j,i) *      &
                                                                       fr_uclass(iurb,j,i)
!
!--                      Sum up the volume which has to be substracted at the end.
                         volairhlred_urb(iz,j,i) = volairhlred_urb(iz,j,i) + vtot2 *               &
                                                                            fr_build(iurb,id,j,i) *&
                                                                            fr_udir(iurb,id,j,i) * &
                                                                            fr_uclass(iurb,j,i)
                      ELSE
                         volhl(iurb,id,iz-k_topo+1,j,i) = 1.0_wp
                      ENDIF
                   ENDDO
                ENDDO
!
!--             Substract the volume of "urban air" above surfaces.
                volairhlred_urb(iz,j,i) = volairhl_urb(iz,j,i) - volairhlred_urb(iz,j,i)
             ENDDO

          ELSE

             volairhl_urb(:,j,i)    = 1.0_wp
             volairhlred_urb(:,j,i) = 1.0_wp
             volhl(:,:,:,j,i)       = 1.0_wp

          ENDIF

       ENDDO
    ENDDO
!
!-- Calculate street direction (ang_udir).
    IF ( ANY( ang_udir(1:n_udir) < -9999.0_wp ) )  THEN
       ang_udir(1:n_udir) = (/ -45.0_wp, 0.0_wp, 45.0_wp, 90.0_wp /)
    ENDIF
!
!-- Convert from degree to radian.
    ang_udir(1:n_udir) = ang_udir(1:n_udir) * pi / 180.0_wp
!
!-- Calculate street direction in rotated system (angrot[x,y]_udir1d).
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  id = 1, n_udir
!
!--          ang_udir relative to N-S axis, clockwise, N is 0 degree;
!--          calculate components in rotated lon (angrotx_udir(id,j,i)) and lat
!--          (angrotx_udir(id,j,i)) direction.
!--          xang = COS(0.5_wp*pi - ang_udir(id))
             xang = SIN( ang_udir(id) )
!
!--          xang = SIN(0.5_wp*pi - ang_udir(id))
             yang = COS( ang_udir(id) )
             DO  iurb = 1, n_uclass
                angrotx_udir(id,j,i) = xang
                angroty_udir(id,j,i) = yang
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculate mean emmissivity.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          emiss_urb(j,i) = 0.0_wp
          IF ( fr_urb(j,i) > eps_urb )  THEN
             DO  id = 1, n_udir
                DO  iurb = 1, n_uclass

                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
                      DO  iu = 1, nz(iurb,id,j,i)+1
                         emiss_urb(j,i) = emiss_urb(j,i) + fr_roof(iurb,id,iu,j,i) *               &
                                                           ( emiss_roof(iurb) *                    &
                                                             width_build(iurb,id,j,i) +            &
                                                             emiss_ground(iurb) *                  &
                                                             width_street(iurb,id,j,i) +           &
                                                             emiss_wall(iurb) * z_uhl(iu) )        &
                                                           / ( width_sglcan(iurb,id,j,i) +         &
                                                               z_uhl(iu) ) * fr_uclass(iurb,j,i) * &
                                                               fr_udir(iurb,id,j,i)
                      ENDDO
                   ENDIF

                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
!
!-- Nodes of evaluation.
    legendrepoly_zero(1) = -1.0_wp / 3.0_wp * SQRT( 5.0_wp + 2.0_wp * SQRT( 10.0_wp / 7.0_wp ) )
    legendrepoly_zero(2) = -1.0_wp / 3.0_wp * SQRT( 5.0_wp - 2.0_wp * SQRT( 10.0_wp / 7.0_wp ) )
    legendrepoly_zero(3) =  0.0_wp
    legendrepoly_zero(4) =  1.0_wp / 3.0_wp * SQRT( 5.0_wp - 2.0_wp * SQRT( 10.0_wp / 7.0_wp ) )
    legendrepoly_zero(5) =  1.0_wp / 3.0_wp * SQRT( 5.0_wp + 2.0_wp * SQRT( 10.0_wp / 7.0_wp ) )
!
!-- Integration weights.
    weight(1) = ( 322.0_wp - 13.0_wp * SQRT( 70.0_wp ) ) / 900.0_wp
    weight(2) = ( 322.0_wp + 13.0_wp * SQRT( 70.0_wp ) ) / 900.0_wp
    weight(3) = 128.0_wp / 225.0_wp
    weight(4) = ( 322.0_wp + 13.0_wp * SQRT( 70.0_wp ) ) / 900.0_wp
    weight(5) = ( 322.0_wp - 13.0_wp * SQRT( 70.0_wp ) ) / 900.0_wp
!
!-- If lrroofs then calculate additional svf using the optional arguments of view_factors.
    DO  i = nxl,nxr
       DO  j = nys, nyn
          DO  id = 1, n_udir
             DO  iurb = 1, n_uclass
                IF ( fr_urb(j,i) > eps_urb  .AND.  fr_udir(iurb,id,j,i) > eps_urb )  THEN
                   IF ( lrroofs )  THEN
                      CALL view_factors( nz(iurb,id,j,i), z_uhl(:), width_street(iurb,id,j,i),     &
                                         fww(iurb,id,:,:,j,i), fwg(iurb,id,:,j,i),                 &
                                         fgw(iurb,id,:,j,i), fsg(iurb,id,j,i), fgs(iurb,id,j,i),   &
                                         fsw(iurb,id,:,j,i), fws(iurb,id,:,j,i),                   &
                                         frw(iurb,id,:,:,j,i), fwr(iurb,id,:,:,j,i),               &
                                         fsr(iurb,id,:,j,i), frs(iurb,id,:,j,i),                   &
                                         width_build(iurb,id,j,i) )
                   ELSE
                      CALL view_factors( nz(iurb,id,j,i), z_uhl(:), width_street(iurb,id,j,i),     &
                                         fww(iurb,id,:,:,j,i), fwg(iurb,id,:,j,i),                 &
                                         fgw(iurb,id,:,j,i), fsg(iurb,id,j,i), fgs(iurb,id,j,i),   &
                                         fsw(iurb,id,:,j,i), fws(iurb,id,:,j,i) )
                   ENDIF
                ELSE
                   fww(iurb,id,:,:,j,i) = 0.0_wp
                   fwg(iurb,id,:,j,i)   = 0.0_wp
                   fgw(iurb,id,:,j,i)   = 0.0_wp
                   fsg(iurb,id,j,i)     = 0.0_wp
                   fgs(iurb,id,j,i)     = 0.0_wp
                   fsw(iurb,id,:,j,i)   = 0.0_wp
                   fws(iurb,id,:,j,i)   = 0.0_wp
                   IF ( lrroofs )  THEN
                      frw(iurb,id,:,:,j,i) = 0.0_wp
                      fwr(iurb,id,:,:,j,i) = 0.0_wp
                      fsr(iurb,id,:,j,i)   = 0.0_wp
                      frs(iurb,id,:,j,i)   = 0.0_wp
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
!-- Calculate double canyon svf.
    IF ( lrroofs )  THEN
       DO  i = nxl, nxr
          DO j = nys, nyn
             DO id = 1, n_udir
                DO iurb = 1, n_uclass

                   IF ( fr_urb(j,i) > eps_urb  .AND.  fr_udir(iurb,id,j,i) > eps_urb)  THEN
!
!--                   Use full array to not create temporary array, loops in function use
!--                   nz(iurb,id,j,i).
                      fgow(iurb,id,:,:,j,i)   = ground_otherwall_vf( i, j, id, iurb )
                      fgos(iurb,id,:,j,i)     = ground_othersky_vf( i, j, id, iurb )
                      fwow(iurb,id,:,:,:,j,i) = wall_otherwall_vf( i, j, id, iurb )
                      fwos(iurb,id,:,:,j,i)   = wall_othersky_vf( i, j, id, iurb )
                   ENDIF

                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Calculate inverse svf for two canyons.
    IF ( lrroofs )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  id = 1, n_udir
                DO  iurb = 1, n_uclass

                   IF ( fr_urb(j,i) > eps_urb  .AND.  fr_udir(iurb,id,j,i) > eps_urb )  THEN
                      DO  ou = 1, nz(iurb,id,j,i)+1
                         fsog(iurb,id,ou,j,i) = width_dblcan(iurb,id,j,i) /                        &
                                                width_street(iurb,id,j,i) * fgos(iurb,id,ou,j,i)
                         DO  iu = 1, nz(iurb,id,j,i)
                            fwog(iurb,id,iu,ou,j,i) = dz_uhl(iu) / width_street(iurb,id,j,i) *     &
                                                      fgow(iurb,id,ou,iu,j,i)
                            fsow(iurb,id,ou,iu,j,i) = width_dblcan(iurb,id,j,i)/dz_uhl(iu) *       &
                                                      fwos(iurb,id,iu,ou,j,i)
                         ENDDO
                      ENDDO
                   ENDIF

                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF ( lrroofs )  THEN
       CALL test_svf_roof()
    ELSE
       CALL test_svf()
    ENDIF

    IF ( lrroofs )  THEN
       raddimfull = 3 * ke_uhl + 2
    ELSE
       raddimfull = 2 * ke_uhl + 1
    ENDIF

    ALLOCATE( tmpmatrix(raddimfull, raddimfull) )
    ALLOCATE( eye(raddimfull, raddimfull) )
!
!-- Calculate radiation matrices and their inverse.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  id = 1, n_udir
             DO  iurb = 1, n_uclass
                IF ( fr_urb(j,i) > eps_urb  .AND.  fr_udir(iurb,id,j,i) > eps_urb )  THEN

                   IF ( lrroofs )  THEN
                      raddim = 3 * nz(iurb,id,j,i) + 2
                   ELSE
                      raddim = 2 * nz(iurb,id,j,i) + 1
                   ENDIF
!
!--                Longwave diffuse radiation.
                   rlmatrix(iurb,id,:,:,j,i) = rad_matrix( i, j, id, iurb,                         &
                                                           raddimfull, raddim,                     &
                                                           1.0_wp - emiss_ground(iurb),            &
                                                           1.0_wp - emiss_wall(iurb),              &
                                                           1.0_wp - emiss_roof(iurb) )
                   eye = 0.0_wp

                   DO  k = 1, raddim
                      eye(k,k) = 1.0_wp
                   ENDDO
!
!--                Save original matrix.
                   tmpmatrix(1:raddim, 1:raddim) = rlmatrix(iurb,id,1:raddim,1:raddim,j,i)

                   CALL gaussj( rlmatrix(iurb,id,:,:,j,i), raddimfull, raddim, eye, raddimfull,    &
                                raddim )
!
!--                Save inverse.
                   rlmatrix(iurb,id,1:raddim,1:raddim,j,i) = eye(1:raddim,1:raddim)
!
!--                Test inverse: this should be identity.
                   tmpmatrix(1:raddim,1:raddim) =                                                  &
                   MATMUL( tmpmatrix(1:raddim,1:raddim), rlmatrix(iurb,id,1:raddim,1:raddim,j,i) )

                   eye = 0.0_wp
                   DO  k = 1, raddim
                      eye(k,k) = 1.0_wp
                   ENDDO
!
!-                 Now just zeros.
                   tmpmatrix(1:raddim,1:raddim) =                                                  &
                                              tmpmatrix(1:raddim,1:raddim) - eye(1:raddim,1:raddim)

                   IF ( MAXVAL( ABS( tmpmatrix(1:raddim,1:raddim) ) ) > 1.0E-12 )  THEN
                      message_string = 'inverse of Longwave radiation matrix not correct'
                      CALL message( 'dcep_mod: dcep_init', 'DCP0005', 1, 2, 0, 6, 0 )
                   ENDIF
!
!--                Shortwave diffuse radiation.
                   rsmatrix(iurb,id,:,:,j,i) = rad_matrix( i, j, id, iurb, raddimfull, raddim,     &
                                                           alb_ground(iurb), alb_wall(iurb),       &
                                                           alb_roof(iurb) )
!
!--                Save original matrix.
                   tmpmatrix(1:raddim,1:raddim) = rsmatrix(iurb,id,1:raddim,1:raddim,j,i)

                   CALL gaussj( rsmatrix(iurb,id,:,:,j,i), raddimfull, raddim, eye, raddimfull,    &
                                raddim )
!
!--                Save inverse.
                   rsmatrix(iurb,id,1:raddim,1:raddim,j,i) = eye(1:raddim,1:raddim)
!
!--                matrix * inverse matrix: this should be identity
                   tmpmatrix(1:raddim, 1:raddim) =                                                 &
                   MATMUL( tmpmatrix(1:raddim, 1:raddim), rsmatrix(iurb,id,1:raddim,1:raddim,j,i) )

                   eye = 0.0_wp
                   DO  k = 1, raddim
                      eye(k,k) = 1.0_wp
                   ENDDO
                   tmpmatrix(1:raddim,1:raddim) =                                                  &
                                              tmpmatrix(1:raddim,1:raddim) - eye(1:raddim,1:raddim)

                   IF ( MAXVAL( ABS( tmpmatrix(1:raddim,1:raddim) ) ) > 1.0E-12 )  THEN
                      message_string = 'inverse of shortwave radiation matrix not correct'
                      CALL message( 'dcep_mod: dcep_init', 'DCP0006', 1, 2, 0, 6, 0 )
                   ENDIF

                ELSE
                   rlmatrix(iurb,id,:,:,j,i) = 0.0_wp
                   rsmatrix(iurb,id,:,:,j,i) = 0.0_wp
                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE( tmpmatrix )
    DEALLOCATE( eye )
!
!-- Radiation reduction factor for radiation from canyon sides.
    IF ( lurbradcor )  THEN

       IF ( lrroofs)  THEN

          DO  i = nxl, nxr
             DO  j = nys, nyn
                IF ( fr_urb(j,i) > eps_urb )  THEN
                   DO  id = 1, n_udir
                      DO  iurb = 1, n_uclass

                         IF ( fr_udir(iurb,id,j,i) > eps_urb)  THEN
!
!--                         Sum up radiation from top for one canyon element.
                            temp = 0.0_wp
                            DO  ou = 1, nz(iurb,id,j,i)+1
!
!--                            Energy on ground.
                               temp = temp + fgos(iurb,id,ou,j,i) * fr_roof(iurb,id,ou,j,i)
!
!--                            On the walls.
                               DO  iu = 1, nz(iurb,id,j,i)
                                  temp = temp + 2.0_wp*fwos(iurb,id,iu,ou,j,i) *                   &
                                                fr_wall(iurb,id,iu,j,i) * fr_roof(iurb,id,ou,j,i)
                               ENDDO
!
!--                            On the walls.
                               temp = temp + frs(iurb,id,ou,j,i) * fr_roof(iurb,id,ou,j,i)
                            ENDDO
!
!--                         *sending/reference
                            temp = temp * width_dblcan(iurb,id,j,i) / width_sglcan(iurb,id,j,i)

                            temp2 = 0.0_wp
                            DO  iu = 1, nz(iurb,id,j,i)
!
!--                            From nearer side.
                               temp2 = temp2 + fgw(iurb,id,iu,j,i) *                               &
                                               ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) * dz_uhl(iu)

                               DO  ou = 1, nz(iurb,id,j,i)+1
!
!--                               From farer side.
                                  temp2 = temp2 + fgow(iurb,id,ou,iu,j,i) *                        &
                                                  ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) * dz_uhl(iu)&
                                                  * fr_roof(iurb,id,ou,j,i)
                                  DO  ju = 1, nz(iurb,id,j,i)
                                     temp2 = temp2 + 2.0_wp * fwow(iurb,id,ju,ou,iu,j,i) *         &
                                                     ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) *        &
                                                     fr_wall(iurb,id,ju,j,i) * dz_uhl(iu)  *       &
                                                     fr_roof(iurb,id,ou,j,i)
                                  ENDDO
!
!--                               On roof.
                                  temp2 = temp2 + 2.0_wp * frw(iurb,id,ou,iu,j,i) *                 &
                                                  ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) * dz_uhl(iu) &
                                                  * fr_roof(iurb,id,ou,j,i)
                               ENDDO
                            ENDDO

                            IF ( temp2 > eps_urb )  THEN
!
!--                            /reference.
                               temp2 =  temp2 / width_sglcan(iurb,id,j,i)
                               radcor(iurb,id,j,i) = ( 1.0_wp - temp ) / temp2
                            ELSE
                               radcor(iurb,id,j,i) = 1.0_wp
                            ENDIF

                         ENDIF

                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO

       ELSE !- lrroofs

          DO  i = nxl, nxr
             DO  j = nys, nyn

                IF ( fr_urb(j,i) > eps_urb )  THEN
                   DO  id = 1, n_udir
                      DO  iurb = 1, n_uclass

                         IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
!
!--                         Sum up svf from top of canyon.
                            temp = fgs(iurb,id,j,i)
                            DO  iu = 1, nz(iurb,id,j,i)
                               temp = temp + 2.0_wp * fws(iurb,id,iu,j,i) * fr_wall(iurb,id,iu,j,i)
                            ENDDO
!
!--                         *sending/reference (width_build part is absorbed by roofs)
!--                         temp = temp*ws1d/ws1d.
                            temp2 = 0.0_wp
                            DO  iu = 1, nz(iurb,id,j,i)
                               temp2 = temp2 + fgw(iurb,id,iu,j,i) *                               &
                                               ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) * dz_uhl(iu)
                               DO  ju = 1, nz(iurb,id,j,i)
                                  temp2 = temp2 + fww(iurb,id,ju,iu,j,i) *                         &
                                                  ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) *           &
                                                  fr_wall(iurb,id,ju,j,i) * dz_uhl(iu)
                               ENDDO
                            ENDDO

                            IF ( temp2 > eps_urb )  THEN
                               temp2 = 2.0_wp * temp2 / width_street(iurb,id,j,i)
                               radcor(iurb,id,j,i) = (1.0_wp - temp) / temp2
                            ELSE
                               radcor(iurb,id,j,i) = 1.0_wp
                            ENDIF

                         ENDIF

                      ENDDO
                   ENDDO
                ENDIF

             ENDDO
          ENDDO

       ENDIF

    ELSE !- lurbradcor
!
!--    No urban radiation correction means factor=1.
       radcor = 1.0_wp
    ENDIF

!
!-- Calculation of albedo for diffuse radiation by solving shortwave radiation matrix with
!-- diffuse insolation=1 and direct insolation=0.
    ALLOCATE( eye(nys:nyn,nxl:nxr) )
    ALLOCATE( tmpmatrix(raddimfull,raddimfull) )
    ALLOCATE( vec(raddimfull) )
    ALLOCATE( one(raddimfull) )
    ALLOCATE( zero(raddimfull) )

    eye  = 1.0_wp
    one  = 1.0_wp
    zero = 0.0_wp

    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  id = 1, n_udir
             DO  iurb = 1, n_uclass

                IF ( fr_urb(j,i) > eps_urb .AND. fr_udir(iurb,id,j,i) > eps_urb )  THEN

                   IF ( lrroofs )  THEN
                      raddim = 3 * nz(iurb,id,j,i) + 2
                   ELSE
                      raddim = 2 * nz(iurb,id,j,i) + 1
                   ENDIF

                   vec(:) = short_rad_inhom( i, j, id, iurb, raddimfull, raddim,                   &
                                             zero(1:nz(iurb,id,j,i)), zero(1:nz(iurb,id,j,i)),     &
                                             zero(1:nz(iurb,id,j,i)+1), zero(1),                   &
                                             one(1:nz(iurb,id,j,i)), one(1:nz(iurb,id,j,i)),       &
                                             one(1:nz(iurb,id,j,i)+1), one(1),                     &
                                             radcor(iurb,id,j,i) )

                   vec(1:raddim) = MATMUL( rsmatrix(iurb,id,1:raddim,1:raddim,j,i), vec(1:raddim) )
!
!--                Temporarily store in output fields (will be overwritten in first time step).
                   DO  iu = 1, nz(iurb,id,j,i)
                      rs_wall(iurb,2*id-1,iu,j,i) = vec(iu)
                   ENDDO

                   DO  iu = nz(iurb,id,j,i)+1, 2*nz(iurb,id,j,i)
                      rs_wall(iurb,2*id,iu-nz(iurb,id,j,i),j,i) = vec(iu)
                   ENDDO

                   rs_ground(iurb,id,j,i) = vec(2*nz(iurb,id,j,i)+1)

                   IF ( lrroofs )  THEN
                      DO  iu = 2*nz(iurb,id,j,i)+2, 3*nz(iurb,id,j,i)+2
                         rs_roof(iurb,id,iu-(2*nz(iurb,id,j,i)+1),j,i) = vec(iu)
                      ENDDO
                   ELSE
                      DO  iu = 1, nz(iurb,id,j,i)+1
                         rs_roof(iurb,id,iu,j,i) = 1.0_wp
                      ENDDO
                   ENDIF

                ENDIF

             ENDDO
          ENDDO
       ENDDO
    ENDDO

    albedo_urb(nys:nyn,nxl:nxr) = uprad( eye, rs_ground, rs_wall, rs_roof, alb_ground(:),          &
                                         alb_wall(:), alb_roof(:), .FALSE.)
    DEALLOCATE( eye )
    DEALLOCATE( tmpmatrix )
    DEALLOCATE( vec )
    DEALLOCATE( one )
    DEALLOCATE( zero )
!
!-- Calculate effective heat transport values.
    DO  iurb = 1, n_uclass
       tddz_ground(:,iurb) = calc_tddz( dz_ground(:,iurb), td_ground(:,iurb), ke_ground )
       tddz_roof(:,iurb)   = calc_tddz( dz_roof(:,iurb),   td_roof(:,iurb),   ke_roof   )
       tddz_wall(:,iurb)   = calc_tddz( dz_wall(:,iurb),   td_wall(:,iurb),   ke_wall   )
    ENDDO
!
!-- Allocate tendancy arrays (not possible in dcep_allocate).
    ALLOCATE( ut_urb(1:nz_mesom,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( vt_urb(1:nz_mesom,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( tt_urb(1:nz_mesom,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( tkes_urb(1:nz_mesom,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( tkeb_urb(1:nz_mesom,nysg:nyng,nxlg:nxrg) )

!
!-- Call the main DCEP routine after initialization.
    CALL dcep_main

    IF ( debug_output )  CALL debug_message( 'dcep_init', 'end' )

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  view factors should add up to 1
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE test_svf()

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  id    !<
    INTEGER(iwp) ::  iu    !<
    INTEGER(iwp) ::  iurb  !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  ju    !<

    REAL(wp) ::  svfsum  !<


!
!-- real svf(i->j)  =  eff. svf(j->i).
    DO  j = nys, nyn
       DO  i = nxl, nxr
          IF ( fr_urb(j,i) > eps_urb )  THEN
             DO  id = 1, n_udir
                DO  iurb = 1, n_uclass
                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
!
!--                   From wall.
                      DO  iu = 1, nz(iurb,id,j,i)
                         svfsum = fgw(iurb,id,iu,j,i) + fsw(iurb,id,iu,j,i)
                         DO  ju = 1, nz(iurb,id,j,i)
                            svfsum = svfsum + fww(iurb,id,ju,iu,j,i)
                         ENDDO
                         IF ( ABS( svfsum - 1.0_wp ) > 1.0E-13 )  THEN
                            WRITE( message_string, * ) 'sum of SVFs at wall element at i= ', i,    &
                                                      ' and j= ', j, ' should be 1 and not ', svfsum
                            CALL message( 'test_svf', 'DCP0007', 1, 2, 0, 6, 0 )
                         ENDIF
!
!--                      Normalization already good, applying it here does not help.
                      ENDDO
!
!--                   From ground.
                      svfsum = fsg(iurb,id,j,i)
                      DO iu = 1, nz(iurb,id,j,i)
                         svfsum = svfsum + 2.0_wp * fwg(iurb,id,iu,j,i)
                      ENDDO

                      IF ( ABS( svfsum - 1.0_wp ) > 1.0E-13 )  THEN
                         WRITE( message_string, * ) 'sum of SVFs at ground element at i= ', i,     &
                                                    ' and j= ', j, ' should be 1 and not ', svfsum
                         CALL message( 'test_svf', 'DCP0007', 1, 2, 0, 6, 0 )
                      ENDIF
!--                   Normalization already good, applying it here does not help

!
!--                   From sky.
                      svfsum = fgs(iurb,id,j,i)
                      DO  iu = 1, nz(iurb,id,j,i)
                         svfsum = svfsum + 2.0_wp * fws(iurb,id,iu,j,i)
                      ENDDO

                      IF ( ABS( svfsum - 1.0_wp ) > 1.d-13 )  THEN
                         WRITE( message_string, * ) 'sum of SVFs at sky element at i= ', i,        &
                                                    ' and j= ', j, ' should be 1 and not ', svfsum
                         CALL message( 'test_svf', 'DCP0007', 1, 2, 0, 6, 0 )
                   ENDIF
!--                   Normalization already good, applying it here does not help

                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    
 END SUBROUTINE test_svf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  test svf roof
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE test_svf_roof()

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  id    !<
    INTEGER(iwp) ::  iu    !<
    INTEGER(iwp) ::  iurb  !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  ju    !<
    INTEGER(iwp) ::  ou    !<

    REAL(wp) ::  svfsum

!
!-- real svf(i->j)  =  eff. svf(j->i).
    DO  j = nys, nyn
       DO  i = nxl, nxr
          IF ( fr_urb(j,i) > eps_urb )  THEN
             DO  id = 1, n_udir
                DO  iurb = 1, n_uclass
                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
!
!--                   From wall.
                      DO  ou = 1, nz(iurb,id,j,i)+1
                         DO  iu = 1, nz(iurb,id,j,i)
                            svfsum = fgw(iurb,id,iu,j,i)
                            svfsum = svfsum + fsow(iurb,id,ou,iu,j,i) + fgow(iurb,id,ou,iu,j,i)  &
                                            + frw(iurb,id,ou,iu,j,i)
                            DO  ju = 1, ou-1
                               svfsum = svfsum + fww(iurb,id,ju,iu,j,i)
                            ENDDO
                            DO  ju = 1, nz(iurb,id,j,i)
                               svfsum = svfsum + fwow(iurb,id,ju,ou,iu,j,i)
                            ENDDO
                            IF ( ABS( svfsum - 1.0_wp ) > 1.0E-14 )  THEN
!
!--                            Normalize to increase energy closure accuracy.
                               svfsum = 1.0_wp / svfsum
                               fgw(iurb,id,iu,j,i)     = svfsum * fgw(iurb,id,iu,j,i)
                               fsow(iurb,id,ou,iu,j,i) = svfsum * fsow(iurb,id,ou,iu,j,i)
                               fgow(iurb,id,ou,iu,j,i) = svfsum * fgow(iurb,id,ou,iu,j,i)
                               frw(iurb,id,ou,iu,j,i)  = svfsum * frw(iurb,id,ou,iu,j,i)
                               DO  ju = 1, ou-1
                                  fww(iurb,id,ju,iu,j,i) = svfsum * fww(iurb,id,ju,iu,j,i)
                               ENDDO
                               DO  ju = 1, nz(iurb,id,j,i)
                                  fwow(iurb,id,ju,ou,iu,j,i) = svfsum * fwow(iurb,id,ju,ou,iu,j,i)
                               ENDDO
                            ENDIF
                         ENDDO
!
!--                      From ground.
                         svfsum = fsog(iurb,id,ou,j,i)
                         DO  iu = 1, nz(iurb,id,j,i)
                            svfsum = svfsum + fwg(iurb,id,iu,j,i)
                         ENDDO
                         DO  iu = 1, ou-1
                            svfsum = svfsum + fwg(iurb,id,iu,j,i)
                         ENDDO
                         DO  iu = 1, nz(iurb,id,j,i)
                            svfsum = svfsum + fwog(iurb,id,iu,ou,j,i)
                         ENDDO
                         IF ( ABS( svfsum - 1.0_wp ) > 1.0E-14 )  THEN
!
!--                         Normalize to increase energy closure accuracy.
                            svfsum = 1.0_wp / svfsum
                            fsog(iurb,id,ou,j,i) = svfsum * fsog(iurb,id,ou,j,i)
                            DO  iu = 1, nz(iurb,id,j,i)
                               fwg(iurb,id,iu,j,i) = svfsum * fwg(iurb,id,iu,j,i)
                               fwog(iurb,id,iu,ou,j,i) = svfsum * fwog(iurb,id,iu,ou,j,i)
                            ENDDO
                         ENDIF
!
!--                      From sky.
                         svfsum = 2.0_wp * fgos(iurb,id,ou,j,i) + frs(iurb,id,ou,j,i)
                         DO  iu = 1, nz(iurb,id,j,i)
                            svfsum = svfsum + 2.0_wp * fwos(iurb,id,iu,ou,j,i)
                         ENDDO
                         DO iu = 1, ou-1
                            svfsum = svfsum + 2.0_wp * fws(iurb,id,iu,j,i) *                      &
                                              width_street(iurb,id,j,i) / width_dblcan(iurb,id,j,i)
                         ENDDO
                         IF ( ABS( svfsum - 1.0_wp) > 1.0E-14 )  THEN
!
!--                         Normalize to increase energy closure accuracy.
                            svfsum = 1.0_wp / svfsum
                            fgos(iurb,id,ou,j,i) = svfsum * fgos(iurb,id,ou,j,i)
                            frs(iurb,id,ou,j,i) = svfsum * frs(iurb,id,ou,j,i)
                            DO  iu = 1, nz(iurb,id,j,i)
                               fwos(iurb,id,iu,ou,j,i) = svfsum * fwos(iurb,id,iu,ou,j,i)
                            ENDDO
                            DO  iu = 1, ou-1
                               fws(iurb,id,iu,j,i) = svfsum * fws(iurb,id,iu,j,i)
                            ENDDO
                         ENDIF
!
!--                      From roof.
                         svfsum = fsr(iurb,id,ou,j,i)
                         DO  iu = 1, nz(iurb,id,j,i)
                            svfsum = svfsum + 2.0_wp * fwr(iurb,id,iu,ou,j,i)
                         ENDDO
                         IF ( ABS( svfsum - 1.0_wp) > 1.0E-14 )  THEN
!
!--                         Normalize to increase energy closure accuracy.
                            svfsum = 1.0_wp / svfsum
                            fsr(iurb,id,ou,j,i) = svfsum * fsr(iurb,id,ou,j,i)
                            DO iu = 1, nz(iurb,id,j,i)
                               fwr(iurb,id,iu,ou,j,i) = svfsum * fwr(iurb,id,iu,ou,j,i)
                            ENDDO
                         ENDIF

                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    
 END SUBROUTINE test_svf_roof


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Gauss-Jordan algorithm
!> This routine solve a linear system of n equations of the form
!>              A X = B
!> where  A is a matrix a(i,j)
!>        B a vector and X the solution
!> In output b is replaced by the solution
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gaussj( a, nfull, n, b, mfull, m )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) ::  n      !<
    INTEGER(iwp), INTENT(in) ::  nfull  !<
    INTEGER(iwp), INTENT(in) ::  m      !<
    INTEGER(iwp), INTENT(in) ::  mfull  !<

    REAL(wp), INTENT(inout), DIMENSION(nfull,nfull) ::  a  !<
    REAL(wp), INTENT(inout), DIMENSION(nfull,mfull) ::  b  !<

!
!-- Local variables
    INTEGER(iwp) ::  i      !<
    INTEGER(iwp) ::  imax   !<
    INTEGER(iwp) ::  itemp  !<
    INTEGER(iwp) ::  j      !<
    INTEGER(iwp) ::  k      !<

    INTEGER(iwp), DIMENSION(n) ::  p  !<

    REAL(wp) ::  fac   !<
    REAL(wp) ::  rmax  !<
    REAL(wp) ::  temp  !<

    REAL(wp), DIMENSION(n,m) ::  x  !<


    IF ( debug_output )  CALL debug_message( 'gaussj', 'start' )

!
!-- Index array so no matrix modification needed when pivoting.
    DO  i = 1, n
       p(i) = i
    ENDDO

    DO  k = 1, n-1

       rmax = ABS( a(p(k),k) )
       imax = k
!
!--    Find max.
       DO  i = k+1, n
          temp = ABS( a(p(i),k) )
          IF ( temp > rmax )  THEN
             rmax = temp
             imax = i
          ENDIF
       ENDDO
!
!--    Exchange.
       itemp   = p(k)
       p(k)    = p(imax)
       p(imax) = itemp

       temp = 1.0_wp / a(p(k),k)
       DO  i = k+1, n

          fac = temp * a(p(i),k)
          a(p(i),k) = 0.0_wp

          DO  j = k+1, n
             a(p(i),j) = a(p(i),j) - fac * a(p(k),j)
          ENDDO
          b(p(i),:) = b(p(i),:) - fac * b(p(k),:)
       ENDDO

    ENDDO

    DO  k = n, 1, -1
       x(k,1:m) = b(p(k),1:m)
       DO  j = k+1, n
          x(k,:) = x(k,:) - a(p(k),j) * x(j,:)
       ENDDO
       x(k,:) = x(k,:) / a(p(k),k)
    ENDDO

    b(1:n,1:m) = x(:,:)


 END SUBROUTINE gaussj


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> average K/dz = 1/2(dz1 K1 + dz2 K2) / 1/2(dz1+dz2) /
!>                    1/2(dz1+dz2)
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION calc_tddz( dz, td, nz )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) ::  nz  !<

    REAL(wp), INTENT(in), DIMENSION(nz) ::  dz  !<
    REAL(wp), INTENT(in), DIMENSION(nz) ::  td  !<

    INTEGER(iwp) ::  i  !<

    REAL(wp) ::  calc_tddz(nz-1)  !<  function type


    DO  i = 1, nz-1
       calc_tddz(i) = 2.0_wp * ( dz(i) * td(i) + dz(i+1) * td(i+1) ) / ( dz(i) + dz(i+1) )**2
    ENDDO

 END FUNCTION calc_tddz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> matrix for radiation exchange (for diffuse exchange)
!> albg, albw, albr give the reflected radiation fraction for ground,
!> wall, roof, resp.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION rad_matrix( i, j, id, iurb, ndimfull, ndim, albg, albw, albr )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) ::  i         !<
    INTEGER(iwp), INTENT(in) ::  id        !<
    INTEGER(iwp), INTENT(in) ::  iurb      !<
    INTEGER(iwp), INTENT(in) ::  j         !<
    INTEGER(iwp), INTENT(in) ::  ndim      !<
    INTEGER(iwp), INTENT(in) ::  ndimfull  !<

    REAL(wp), INTENT(in) ::  albg  !<
    REAL(wp), INTENT(in) ::  albr  !<
    REAL(wp), INTENT(in) ::  albw  !<

    REAL(wp) ::  rad_matrix(ndimfull,ndimfull)  !<  function type
!
!-- Local variables
    INTEGER(iwp) ::  iu  !< loop index
    INTEGER(iwp) ::  ju  !< loop index
    INTEGER(iwp) ::  ou  !< loop index


    rad_matrix = 12345.0_wp
!
!-- West wall (A7).
    DO  iu = 1, nz(iurb,id,j,i)
!
!--    No terms from other parts of western wall.
       DO  ju = 1, nz(iurb,id,j,i)
          rad_matrix(iu,ju) = 0.0_wp
       ENDDO

       rad_matrix(iu,iu) = 1.0_wp
!
!--    Reflection from eastern wall:
       DO  ju = nz(iurb,id,j,i)+1, 2*nz(iurb,id,j,i)
          rad_matrix(iu,ju) = -albw * fww(iurb,id,ju-nz(iurb,id,j,i),iu,j,i) *                   &
                                      fr_wall(iurb,id,ju-nz(iurb,id,j,i),j,i)
       ENDDO
!
!--    From ground:
       rad_matrix(iu,2*nz(iurb,id,j,i)+1) = -albg * fgw(iurb,id,iu,j,i)

       IF ( lrroofs )  THEN
          DO  ou = 1, nz(iurb,id,j,i)+1
!
!--          Reflection from eastern wall of other canyon:
             DO  ju = nz(iurb,id,j,i)+1, 2*nz(iurb,id,j,i)
                rad_matrix(iu,ju) = rad_matrix(iu,ju) - albw *                                     &
                                                        fwow(iurb,id,ju-nz(iurb,id,j,i),ou,iu,j,i) &
                                                        * fr_wall(iurb,id,ju-nz(iurb,id,j,i),j,i)  &
                                                        * fr_roof(iurb,id,ou,j,i)
             ENDDO
!
!--          Also radiation from ground of other canyon:
             rad_matrix(iu,2*nz(iurb,id,j,i)+1) =                                                  &
             rad_matrix(iu,2*nz(iurb,id,j,i)+1) - albg * fgow(iurb,id,ou,iu,j,i) *                 &
                                                  fr_roof(iurb,id,ou,j,i)
          ENDDO
!
!--       From roof:
          DO  ju = 2*nz(iurb,id,j,i)+2, 2*nz(iurb,id,j,i)+1+iu
             rad_matrix(iu,ju) = -albr * frw(iurb,id,ju-(2*nz(iurb,id,j,i)+1),iu,j,i) *          &
                                 fr_roof(iurb,id,ju-(2*nz(iurb,id,j,i)+1),j,i)
          ENDDO

          DO  ju = 2*nz(iurb,id,j,i)+2+iu, ndim
             rad_matrix(iu,ju) = 0.0_wp
          ENDDO
       ENDIF
    ENDDO
!
!-- East wall.
    DO  iu = 1+nz(iurb,id,j,i), 2*nz(iurb,id,j,i)
!
!--    From west wall:
       DO  ju = 1, nz(iurb,id,j,i)
          rad_matrix(iu,ju) = -albw * fww(iurb,id,ju,iu-nz(iurb,id,j,i),j,i) *                   &
                              fr_wall(iurb,id,ju,j,i)
       ENDDO
!
!--    From east wall:
       DO  ju = 1+nz(iurb,id,j,i), 2*nz(iurb,id,j,i)
          rad_matrix(iu,ju) = 0.0_wp
       ENDDO
!
!-     Itsself:
       rad_matrix(iu,iu) = 1.0_wp
!
!--    From ground:
       rad_matrix(iu,2*nz(iurb,id,j,i)+1) = -albg * fgw(iurb,id,iu-nz(iurb,id,j,i),j,i)

       IF ( lrroofs )  THEN

          DO  ou = 1, nz(iurb,id,j,i)+1
!
!--          From west wall:
             DO  ju = 1, nz(iurb,id,j,i)
                rad_matrix(iu,ju) = rad_matrix(iu,ju) -                                            &
                                    albw * fwow(iurb,id,ju,ou,iu-nz(iurb,id,j,i),j,i) *            &
                                    fr_wall(iurb,id,ju,j,i) * fr_roof(iurb,id,ou,j,i)
             ENDDO

             rad_matrix(iu,2*nz(iurb,id,j,i)+1) = rad_matrix(iu,2*nz(iurb,id,j,i)+1) -             &
                                                  albg * fgow(iurb,id,ou,iu-nz(iurb,id,j,i),j,i) * &
                                                  fr_roof(iurb,id,ou,j,i)

          ENDDO
!
!--       From roof:
          DO  ju = 2*nz(iurb,id,j,i)+2, nz(iurb,id,j,i)+1+iu
             rad_matrix(iu,ju) =                                                                   &
                            -albr * frw(iurb,id,ju-(2*nz(iurb,id,j,i)+1),iu-nz(iurb,id,j,i),j,i) * &
                            fr_roof(iurb,id,ju-(2*nz(iurb,id,j,i)+1),j,i)
          ENDDO
!
!--       Rest:
          DO  ju = nz(iurb,id,j,i)+2+iu, ndim
             rad_matrix(iu,ju) = 0.0_wp
          ENDDO

       ENDIF

    ENDDO
!
!-- Ground.
!-- From west wall:
    DO  ju = 1, nz(iurb,id,j,i)
       rad_matrix(2*nz(iurb,id,j,i)+1,ju) = -albw * fwg(iurb,id,ju,j,i) * fr_wall(iurb,id,ju,j,i)
    ENDDO
!
!-- From east wall:
    DO  ju = nz(iurb,id,j,i)+1, 2*nz(iurb,id,j,i)
       rad_matrix(2*nz(iurb,id,j,i)+1,ju) = -albw * fwg(iurb,id,ju-nz(iurb,id,j,i),j,i) *          &
                                            fr_wall(iurb,id,ju-nz(iurb,id,j,i),j,i)
    ENDDO
!
!-- Itself:
    rad_matrix(2*nz(iurb,id,j,i)+1,2*nz(iurb,id,j,i)+1) = 1.0_wp

    IF ( lrroofs )  THEN
       DO  ju = 1, nz(iurb,id,j,i)
          DO  ou = 1, nz(iurb,id,j,i)+1
             rad_matrix(2*nz(iurb,id,j,i)+1,ju) = rad_matrix(2*nz(iurb,id,j,i)+1,ju) -             &
                                                  albw * fwog(iurb,id,ju,ou,j,i) *                 &
                                                  fr_wall(iurb,id,ju,j,i) * fr_roof(iurb,id,ou,j,i)
             rad_matrix(2*nz(iurb,id,j,i)+1,ju+nz(iurb,id,j,i)) =                                  &
                                                rad_matrix(2*nz(iurb,id,j,i)+1,ju+nz(iurb,id,j,i)) &
                                              - albw * fwog(iurb,id,ju,ou,j,i) *                   &
                                                fr_wall(iurb,id,ju,j,i) * fr_roof(iurb,id,ou,j,i)
          ENDDO
       ENDDO
!
!--    Nothing from roof:
       DO  ju = 2*nz(iurb,id,j,i)+2, ndim
          rad_matrix(2*nz(iurb,id,j,i)+1,ju) = 0.0_wp
       ENDDO
    ENDIF
!
!-- Roofs:
    IF ( lrroofs )  THEN
       DO  iu = 2*nz(iurb,id,j,i)+2, ndim
!
!--       Nothing from the low levels:
          DO  ju = 1, iu-(2*nz(iurb,id,j,i)+2)
             rad_matrix(iu,ju) = 0.0_wp
          ENDDO
!
!--       From west wall:
          DO  ju = iu-(2*nz(iurb,id,j,i)+1), nz(iurb,id,j,i)
             rad_matrix(iu,ju) = -albw * fwr(iurb,id,ju,iu-(2*nz(iurb,id,j,i)+1),j,i) *            &
                                 fr_wall(iurb,id,ju,j,i)
          ENDDO
!
!--       Nothing from the low levels:
          DO  ju = nz(iurb,id,j,i)+1, iu-nz(iurb,id,j,i)-2
             rad_matrix(iu,ju) = 0.0_wp
          ENDDO
!
!--       From east wall:
!--       @TOTO: explain following comment:
!--       nz(iurb,id,j,i)+iu-(2*nz(iurb,id,j,i)+1)
          DO  ju = iu-nz(iurb,id,j,i)-1, 2*nz(iurb,id,j,i)
             rad_matrix(iu,ju) = -albw *                                                           &
                                 fwr(iurb,id,ju-nz(iurb,id,j,i),iu-(2*nz(iurb,id,j,i)+1),j,i) *    &
                                 fr_wall(iurb,id,ju-nz(iurb,id,j,i),j,i)
          ENDDO
!
!--       From ground:
          rad_matrix(iu,2*nz(iurb,id,j,i)+1) = 0.0_wp
!
!--       Nothing from other roofs:
          DO  ju = 2*nz(iurb,id,j,i)+2, ndim
             rad_matrix(iu,ju) = 0.0_wp
          ENDDO
!
!--       Itself:
          rad_matrix(iu,iu) = 1.0_wp

       ENDDO

    ENDIF

 END FUNCTION rad_matrix


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function calculates the view factor from ground to wall.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION ground_otherwall_vf( i, j, id, iurb )

    INTEGER(iwp), INTENT(in) ::  i     !< index
    INTEGER(iwp), INTENT(in) ::  j     !< index
    INTEGER(iwp), INTENT(in) ::  id    !< index
    INTEGER(iwp), INTENT(in) ::  iurb  !< index

    INTEGER(iwp) ::  ri  !< roof index
    INTEGER(iwp) ::  wi  !< wall index

    REAL(wp) ::  ground_otherwall_vf(ke_uhl+1, ke_uhl)  !< function type
    REAL(wp) ::  temp                                   !<
    REAL(wp) ::  sfullvis                               !<
    REAL(wp) ::  snonvis                                !<


!
!-- Height of the building in middle.
    DO  ri = 1, nz(iurb,id,j,i)+1
       DO  wi = 1, nz(iurb,id,j,i)
!
!--       End of/from wall surface not visible fraction measured on the horizontal away from
!--       building with wall.
          IF ( z_uhl(ri) > 0.0_wp )  THEN
             IF ( z_uhl(ri) < z_uhl(wi+1) )  THEN
                snonvis = MIN( ( z_uhl(ri) / (z_uhl(wi+1) - z_uhl(ri) ) + 1.0_wp ) *               &
                                width_sglcan(iurb,id,j,i), width_dblcan(iurb,id,j,i) )
             ELSE
                snonvis = width_dblcan(iurb,id,j,i)
             ENDIF
          ELSE
             snonvis = width_sglcan(iurb,id,j,i)
          ENDIF
!
!--       Is there a visible ground fraction?
          IF ( snonvis < width_dblcan(iurb,id,j,i) )  THEN
!
!--          Begin of fully visible fraction
             IF ( z_uhl(ri) > 0.0_wp )  THEN
                IF ( z_uhl(ri) < z_uhl(wi) )  THEN
                   sfullvis = MIN( ( z_uhl(ri) / (z_uhl(wi) - z_uhl(ri) ) + 1.0_wp ) *             &
                                   width_sglcan(iurb,id,j,i), width_dblcan(iurb,id,j,i) )
                ELSE
                   sfullvis = width_dblcan(iurb,id,j,i)
                ENDIF
             ELSE
                sfullvis = width_sglcan(iurb,id,j,i)
             ENDIF

!
!--          Only the part that has varying visibility of receiving wall has to be integrated
!--          (sending area of that part is already included).
             ground_otherwall_vf(ri,wi) =                                                          &
                       integral( snonvis, sfullvis, ground_otherwall_vf_f, i, j, id, iurb, ri, wi )
!
!--          Plus part that is fully visible (sending area of that part is included).
             temp = fnrm14( sfullvis, width_dblcan(iurb,id,j,i), z_uhl(wi), z_uhl(wi + 1) )
             ground_otherwall_vf(ri,wi) = ground_otherwall_vf(ri,wi) + temp
!
!--          /sending * sending/receiving
             ground_otherwall_vf(ri,wi) = ground_otherwall_vf(ri,wi) / ( z_uhl(wi+1) - z_uhl(wi) )
          ELSE
             ground_otherwall_vf(ri,wi) = 0.0_wp
          ENDIF
       ENDDO
    ENDDO

 END FUNCTION ground_otherwall_vf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function calculates the view factor from ground to wall.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION wall_otherwall_vf( i, j, id, iurb )

    INTEGER(iwp), INTENT(in) ::  i     !<
    INTEGER(iwp), INTENT(in) ::  id    !<
    INTEGER(iwp), INTENT(in) ::  iurb  !<
    INTEGER(iwp), INTENT(in) ::  j     !<

    REAL(wp) ::  wall_otherwall_vf(ke_uhl,ke_uhl+1,ke_uhl)  !<  function type

!
!-- Local variables
    INTEGER(iwp) ::  ri   !<
    INTEGER(iwp) ::  wi   !<
    INTEGER(iwp) ::  wi2  !<

    REAL(wp) ::  sfullvis  !<
    REAL(wp) ::  snonvis   !<
    REAL(wp) ::  temp      !<


!
!-- Height of the building in middle.
    DO  ri = 1, nz(iurb,id,j,i)+1
       DO  wi = 1, nz(iurb,id,j,i)
!
!--       Height of non visible fraction and height of full visible.
          IF ( z_uhl(ri) > 0 )  THEN
             IF ( z_uhl(ri) < z_uhl(wi+1) )  THEN
                snonvis = MAX( ( width_dblcan(iurb,id,j,i) * z_uhl(ri) -                           &
                                 width_street(iurb,id,j,i) * z_uhl(wi+1)                           &
                                ) / width_sglcan(iurb,id,j,i), 0.0_wp )
             ELSE
                snonvis = MIN( width_dblcan(iurb,id,j,i) / width_street(iurb,id,j,i) *             &
                               ( z_uhl(ri) - z_uhl(wi+1) ) + z_uhl(wi+1), z_uhl(nz(iurb,id,j,i)+1) )
             ENDIF
             IF ( z_uhl(ri) < z_uhl(wi) )  THEN
                sfullvis = MAX( ( width_dblcan(iurb,id,j,i) * z_uhl(ri) -                          &
                                  width_street(iurb,id,j,i) * z_uhl(wi)                            &
                                ) / width_sglcan(iurb,id,j,i), 0.0_wp )
             ELSE
                sfullvis = MIN( width_dblcan(iurb,id,j,i) / width_street(iurb,id,j,i) *            &
                                (z_uhl(ri) - z_uhl(wi)) + z_uhl(wi), z_uhl(nz(iurb,id,j,i)+1) )
             ENDIF
          ELSE
             snonvis  = 0.0_wp
             sfullvis = 0.0_wp
          ENDIF

          DO  wi2 = wi, nz(iurb,id,j,i)
             IF ( snonvis < z_uhl(wi2+1) )  THEN
!
!--             Only the part that has varying visibility of receiving wall has to be integrated
!--             (sending area of that part is already included).
                IF ( sfullvis > z_uhl(wi2) )  THEN
                   wall_otherwall_vf(wi,ri,wi2) = integral( MAX( snonvis, z_uhl(wi2) ),            &
                                                            MIN( sfullvis, z_uhl(wi2+1) ),         &
                                                            wall_otherwall_vf_f, i, j, id, iurb,   &
                                                            ri, wi )
                ELSE
                   wall_otherwall_vf(wi,ri,wi2) = 0.0_wp
                ENDIF

                IF ( sfullvis < z_uhl(wi2+1) )  THEN
!
!--                Plus part that is fully visible (sending area of that part is included).
                   temp = fprl16( z_uhl(wi), z_uhl(wi+1), MAX( sfullvis, z_uhl(wi2) ),             &
                                  z_uhl(wi2+1),                                                    &
                                  2 * width_street(iurb,id,j,i) + width_build(iurb,id,j,i) )
                   wall_otherwall_vf(wi,ri,wi2) = wall_otherwall_vf(wi,ri,wi2) + temp
                ENDIF
!
!--             /sending * sending/receiving
                wall_otherwall_vf(wi,ri,wi2) = wall_otherwall_vf(wi,ri,wi2) /                      &
                                               ( z_uhl(wi2+1) - z_uhl(wi2) )
                wall_otherwall_vf(wi2,ri,wi) = wall_otherwall_vf(wi,ri,wi2) *                      &
                                               ( z_uhl(wi2+1) - z_uhl(wi2) ) /                     &
                                               ( z_uhl(wi+1)  - z_uhl(wi)  )
             ELSE
                wall_otherwall_vf(wi,ri,wi2) = 0.0_wp
                wall_otherwall_vf(wi2,ri,wi) = 0.0_wp
             ENDIF
          ENDDO
       ENDDO
    ENDDO

 END FUNCTION wall_otherwall_vf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function calculates the view factor from ground to sky.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION ground_othersky_vf( i, j, id, iurb )

    INTEGER(iwp), INTENT(in) ::  i     !<
    INTEGER(iwp), INTENT(in) ::  id    !<
    INTEGER(iwp), INTENT(in) ::  j     !<
    INTEGER(iwp), INTENT(in) ::  iurb  !<

    REAL(wp) ::  ground_othersky_vf(ke_uhl+1)  !< function type

!
!-- Local variables
    INTEGER(iwp) ::  ri  !< roof index

    REAL(wp) ::  spartvis  !<
    REAL(wp) ::  temp      !<


!
!-- Height of the building in middle.
    DO  ri = 1, nz(iurb,id,j,i)+1
!
!--    Partly visible length on the ground.
       IF ( z_uhl(ri) < z_uhl(nz(iurb,id,j,i)+1) )  THEN
          spartvis = MIN( z_uhl(ri) / ( z_uhl(nz(iurb,id,j,i)+1) - z_uhl(ri) ) *                   &
                          width_sglcan(iurb,id,j,i), width_street(iurb,id,j,i) )
       ELSE
          spartvis = width_street(iurb,id,j,i)
       ENDIF
!
!--    ONLY the part that has varying visibility of receiving wall has to be integrated (sending
!--    area of that part is already included).
       ground_othersky_vf(ri) = integral( width_street(iurb,id,j,i) - spartvis,                    &
                                          width_street(iurb,id,j,i), ground_othersky_vf_f, i, j,   &
                                          id, iurb, ri, -1 )
!
!--    Plus part that is fully visible (sending area of that part is included).
       temp = fprl134( width_street(iurb,id,j,i) - spartvis, width_dblcan(iurb,id,j,i),            &
                       z_uhl(nz(iurb,id,j,i)+1) )
       ground_othersky_vf(ri) = ground_othersky_vf(ri) + temp
!
!--    /sending * sending/receiving
       ground_othersky_vf(ri) = ground_othersky_vf(ri) / width_dblcan(iurb,id,j,i)
    ENDDO

 END FUNCTION ground_othersky_vf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function calculates the view factor from wall to sky.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION wall_othersky_vf( i, j, id, iurb )

    INTEGER(iwp), INTENT(in) ::  i     !<
    INTEGER(iwp), INTENT(in) ::  iurb  !<
    INTEGER(iwp), INTENT(in) ::  id    !<
    INTEGER(iwp), INTENT(in) ::  j     !<

    REAL(wp) ::  wall_othersky_vf(ke_uhl,ke_uhl+1)  !< function type

!
!-- Local variables
    INTEGER(iwp) ::  ri  !< roof index
    INTEGER(iwp) ::  wi  !< wall index

    REAL(wp) ::  sfullvis  !<
    REAL(wp) ::  temp      !<


!
!-- Height of the building in middle, highest case treated extra.
    DO  ri = 1, nz(iurb,id,j,i)
!
!--    Partly visible length.
       sfullvis = MAX( ( width_dblcan(iurb,id,j,i) * z_uhl(ri) -                                   &
                         width_street(iurb,id,j,i) * z_uhl(nz(iurb,id,j,i)+1) )                    &
                       / width_sglcan(iurb,id,j,i), 0.0_wp )
       DO  wi = 1, nz(iurb,id,j,i)
          IF ( sfullvis > z_uhl(wi) )  THEN
!
!--          Only the part that has varying visibility of receiving wall has to be integrated
!--          (sending area of that part is already included)
             wall_othersky_vf(wi,ri) = integral( z_uhl(wi), MIN( z_uhl(wi+1), sfullvis ),          &
                                                 wall_othersky_vf_f, i, j, id, iurb, ri, wi )
          ELSE
             wall_othersky_vf(wi,ri) = 0.0_wp
          ENDIF
          IF ( sfullvis < z_uhl(wi+1) )  THEN
!
!--          Plus part that is fully visible (sending area of that part is included).
             temp = fnrm13( z_uhl(nz(iurb,id,j,i)+1) - MAX( sfullvis, z_uhl(wi) ),                 &
                            z_uhl(nz(iurb,id,j,i)+1) - z_uhl(wi+1), width_dblcan(iurb,id,j,i) )
             wall_othersky_vf(wi,ri) = wall_othersky_vf(wi,ri) + temp
          ENDIF
!
!--       /sending * sending/receiving
          wall_othersky_vf(wi,ri) = wall_othersky_vf(wi,ri) / width_dblcan(iurb,id,j,i)
       ENDDO
    ENDDO
    DO wi = 1, nz(iurb,id,j,i)
       wall_othersky_vf(wi,nz(iurb,id,j,i)+1) = fnrm13( z_uhl(nz(iurb,id,j,i)+1) - z_uhl(wi),      &
                                                        z_uhl(nz(iurb,id,j,i)+1) - z_uhl(wi+1),    &
                                                        width_street(iurb,id,j,i) )
       wall_othersky_vf(wi,nz(iurb,id,j,i)+1) = wall_othersky_vf(wi,nz(iurb,id,j,i)+1) /           &
                                                width_dblcan(iurb,id,j,i)
    ENDDO

 END FUNCTION wall_othersky_vf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function calculates the view factor from wall to wall.
!--------------------------------------------------------------------------------------------------!
 PURE REAL(wp) FUNCTION wall_otherwall_vf_f( x, i, j, id, iurb, ri, wi )

    INTEGER(iwp), INTENT(in) ::  i     !<
    INTEGER(iwp), INTENT(in) ::  id    !<
    INTEGER(iwp), INTENT(in) ::  iurb  !<
    INTEGER(iwp), INTENT(in) ::  j     !<
    INTEGER(iwp), INTENT(in) ::  ri    !< roof index
    INTEGER(iwp), INTENT(in) ::  wi    !< wall index

    REAL(wp), INTENT(in) ::  x  !<

!
!-- Local variables
    REAL(wp) ::  hvis  !<
    REAL(wp) ::  t1    !<
    REAL(wp) ::  t2    !<


    IF ( z_uhl(ri) <= x )  THEN
       hvis = ( width_dblcan(iurb,id,j,i) * z_uhl(ri) - width_street(iurb,id,j,i) * x ) /        &
              width_sglcan(iurb,id,j,i)
       t1 = fprls_da( x - hvis,        width_dblcan(iurb,id,j,i) )
       t2 = fprls_da( x - z_uhl(wi+1), width_dblcan(iurb,id,j,i) )
    ELSE
       hvis = x + ( z_uhl(ri) - x ) * width_dblcan(iurb,id,j,i) / width_street(iurb,id,j,i)
       t1 = fprls_da( z_uhl(wi) - x, width_dblcan(iurb,id,j,i) )
       t2 = fprls_da( hvis - x,      width_dblcan(iurb,id,j,i) )
    ENDIF
    wall_otherwall_vf_f = t1 - t2

 END FUNCTION wall_otherwall_vf_f


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function calculates the view factor from wall to sky.
!--------------------------------------------------------------------------------------------------!
 PURE REAL(wp) FUNCTION wall_othersky_vf_f( x, i, j, id, iurb, ri, wi )

    INTEGER(iwp), INTENT(in) ::  i     !<
    INTEGER(iwp), INTENT(in) ::  id    !<
    INTEGER(iwp), INTENT(in) ::  iurb  !<
    INTEGER(iwp), INTENT(in) ::  j     !<
    INTEGER(iwp), INTENT(in) ::  ri    !< roof index
    INTEGER(iwp), INTENT(in) ::  wi    !< wi is handed over here only for consistency to other functions, i.e.
                                       !< ground_otherwall_vf_f, which are the argument f of the function integeral

    REAL(wp), INTENT(in) ::  x  !<
!
!-- Local variables
    REAL(wp) ::  l  !<

!
!-- To avoid compiler warning about unused variable.
    l = wi

    l = width_street(iurb,id,j,i) * ( z_uhl(nz(iurb,id,j,i)+1) - x ) / ( z_uhl(ri) - x )
    wall_othersky_vf_f = fnrms_da( l, z_uhl(nz(iurb,id,j,i)+1) - x )

 END FUNCTION wall_othersky_vf_f


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function calculates the view factor from ground to sky.
!--------------------------------------------------------------------------------------------------!
 PURE REAL(wp) FUNCTION ground_othersky_vf_f( x, i, j, id, iurb, ri, wi )

    INTEGER(iwp), INTENT(in) ::  i     !<
    INTEGER(iwp), INTENT(in) ::  id    !<
    INTEGER(iwp), INTENT(in) ::  iurb  !<
    INTEGER(iwp), INTENT(in) ::  j     !<
    INTEGER(iwp), INTENT(in) ::  ri    !<
    INTEGER(iwp), INTENT(in) ::  wi    !< wi is handed over here only for consistency to other functions, i.e.
                                       !< ground_otherwall_vf_f, which are the argument f of the function integeral
    REAL(wp), INTENT(in) ::  x  !<
!
!-- Local variables
    REAL(wp) ::  skyvis  !<
    REAL(wp) ::  temp1   !<
    REAL(wp) ::  temp2   !<


!
!-- To avoid compiler warning about unused variable.
    temp1 = wi

    IF ( z_uhl(ri) > 0.0_wp )  THEN
       skyvis = z_uhl(nz(iurb,id,j,i)+1) / z_uhl(ri) * ( width_street(iurb,id,j,i) - x )
    ELSE
       skyvis = width_dblcan(iurb,id,j,i) - x
    ENDIF

    temp1 = fprls_da( x,      z_uhl(nz(iurb,id,j,i)+1) )
    temp2 = fprls_da( skyvis, z_uhl(nz(iurb,id,j,i)+1) )

    ground_othersky_vf_f = temp1 + temp2

 END FUNCTION ground_othersky_vf_f


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function  calculates the view factor from ground to wall.
!--------------------------------------------------------------------------------------------------!
 PURE REAL(wp) FUNCTION ground_otherwall_vf_f( x, i, j, id, iurb, ri, wi )

    INTEGER(iwp), INTENT(in) ::  i     !<
    INTEGER(iwp), INTENT(in) ::  id    !<
    INTEGER(iwp), INTENT(in) ::  iurb  !<
    INTEGER(iwp), INTENT(in) ::  j     !<
    INTEGER(iwp), INTENT(in) ::  ri    !<
    INTEGER(iwp), INTENT(in) ::  wi    !<

    REAL(wp), INTENT(in) ::  x  !<

!
!-- Local variables
    REAL(wp) ::  rec                    !<
    REAL(wp) ::  temp1                  !<
    REAL(wp) ::  temp2                  !<


    rec = MAX( z_uhl(wi+1) - MAX( x * z_uhl(ri) / ( x - width_sglcan(iurb,id,j,i) ), z_uhl(wi) ),  &
               0.0_wp )

    temp1 = fnrms_da( z_uhl(wi+1),       x)
    temp2 = fnrms_da( z_uhl(wi+1) - rec, x)

    ground_otherwall_vf_f = temp1 - temp2

 END FUNCTION ground_otherwall_vf_f


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Function for numerical integration using Gauss-Legendre integration.
!--------------------------------------------------------------------------------------------------!
 PURE REAL(wp) FUNCTION integral( a, b, f, i, j, id, iurb, ri, wi )

    INTEGER(iwp), INTENT(in) ::  i      !< roof index, wall index
    INTEGER(iwp), INTENT(in) ::  id     !<
    INTEGER(iwp), INTENT(in) ::  j      !<
    INTEGER(iwp), INTENT(in) ::  iurb   !<
    INTEGER(iwp), INTENT(in) ::  ri     !<
    INTEGER(iwp), INTENT(in) ::  wi     !<

    REAL(wp), INTENT(in) ::  a  !<
    REAL(wp), INTENT(in) ::  b  !<

    INTERFACE
       PURE REAL(wp) FUNCTION f( x, i, j, id, iurb, ri, wi )
         USE kinds
         REAL(wp), INTENT(in) ::  x
         INTEGER(iwp), INTENT(in) ::  i     !<
         INTEGER(iwp), INTENT(in) ::  id    !<
         INTEGER(iwp), INTENT(in) ::  j     !<
         INTEGER(iwp), INTENT(in) ::  iurb  !<
         INTEGER(iwp), INTENT(in) ::  ri    !<
         INTEGER(iwp), INTENT(in) ::  wi    !<
       END FUNCTION f
    END INTERFACE
!
!-- Local variables
    INTEGER(iwp) ::  k  !<

    INTEGER(iwp), PARAMETER ::  maxk = 33554432  !<

    REAL(wp), PARAMETER ::  eps = 8.0E-16  !<

    REAL(wp) ::  errorfac     !<
    REAL(wp) ::  int_step_n   !<
    REAL(wp) ::  int_step_n1  !<


    IF ( a == b )  THEN
       integral = 0
       RETURN
    ENDIF

    int_step_n  = 0.0_wp
    int_step_n1 = 0.0_wp

    errorfac = 2.0_wp**(2*n_num_integ)

    k = 1
    int_step_n = int_step( a, b, f, k, i, j, id, iurb, ri, wi )

    integral = 0.0_wp

    DO WHILE ( k < maxk )
       k = k*2
       int_step_n1 = int_step( a, b, f, k, i, j, id, iurb, ri, wi )
!
!--    Error is N^(2*nx), so factor 2^(2*nx) between stepn and stepn1.
       integral = ( errorfac * int_step_n1 - int_step_n ) / ( errorfac - 1.0_wp )
       IF ( ABS( integral ) > 1.0E-15 )  THEN
          IF ( ABS( (integral - int_step_n1) / integral ) < eps )  RETURN
       ELSE
          IF ( ABS( integral - int_step_n1 ) < eps )  RETURN
       ENDIF
       int_step_n = int_step_n1
    ENDDO

 END FUNCTION integral


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> One step in the integration.
!--------------------------------------------------------------------------------------------------!
 PURE REAL(wp) FUNCTION int_step( a, b, f, n_intervalls, i, j, id, iurb, ri, wi )

    INTEGER(iwp), INTENT(in) :: n_intervalls  !<

    REAL(wp), INTENT(in) ::  a  !<
    REAL(wp), INTENT(in) ::  b  !<

    INTERFACE
       PURE REAL(wp) FUNCTION f( x, i, j, id, iurb, ri, wi )
         USE kinds
         INTEGER(iwp), INTENT(in) ::  i     !<
         INTEGER(iwp), INTENT(in) ::  id    !<
         INTEGER(iwp), INTENT(in) ::  j     !<
         INTEGER(iwp), INTENT(in) ::  iurb  !<
         INTEGER(iwp), INTENT(in) ::  ri    !<
         INTEGER(iwp), INTENT(in) ::  wi    !<
         REAL(wp), INTENT(in) ::  x
       END FUNCTION f
    END INTERFACE

    INTEGER(iwp), INTENT(in) ::  i     !< roof index, wall index
    INTEGER(iwp), INTENT(in) ::  id    !<
    INTEGER(iwp), INTENT(in) ::  iurb  !<
    INTEGER(iwp), INTENT(in) ::  j     !<
    INTEGER(iwp), INTENT(in) ::  ri    !<
    INTEGER(iwp), INTENT(in) ::  wi    !<
!
!-- Local variables
    INTEGER(iwp) ::  ix  !<
    INTEGER(iwp) ::  jx  !<

    REAL(wp) ::  h   !<
    REAL(wp) ::  h2  !<
    REAL(wp) ::  is  !<
    REAL(wp) ::  x   !<
    REAL(wp) ::  y   !<

    REAL(wp), DIMENSION(n_num_integ) :: modx  !<


    h = ( b - a ) / n_intervalls
    h2 = h / 2.0_wp

    DO  ix = 1, n_num_integ
       modx(ix) = legendrepoly_zero(ix) * h2
    ENDDO

    is = 0.0_wp
    x  = a + h2
    DO  ix = 1, n_intervalls
       DO  jx = 1, n_num_integ
          y = f( modx(jx) + x, i, j, id, iurb, ri, wi )
          is = is + weight(jx) * y
       ENDDO
       x = x + h
    ENDDO

    int_step = h2 * is

 END FUNCTION int_step

 END SUBROUTINE dcep_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate tendency due to urban area consideration using DCEP (Vector-optimized version).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_tendency( component )

    IMPLICIT NONE

    INTEGER(iwp) ::  component   !< prognostic variable
    INTEGER(iwp) ::  i           !< loop indices
    INTEGER(iwp) ::  ik          !< loop indices
    INTEGER(iwp) ::  j           !< loop indices
    INTEGER(iwp) ::  k           !< loop indices
    INTEGER(iwp) ::  k_topo      !< topography top index


    SELECT CASE ( component )
!
!--    u-component
       CASE ( 1 )
          IF ( lurbvel )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   k_topo = topo_top_ind(j,i,0)
                   DO  ik = 1, nz_meso(j,i)
                      k = k_topo + ik
                      tend(k,j,i) = tend(k,j,i) + fr_urb(j,i) * ut_urb(ik,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
!
!--    v-component
       CASE ( 2 )
          IF ( lurbvel )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   k_topo = topo_top_ind(j,i,0)
                   DO  ik = 1, nz_meso(j,i)
                      k = k_topo + ik
                      tend(k,j,i) = tend(k,j,i) + fr_urb(j,i) * vt_urb(ik,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
!
!--    w-component (not implemented yet)
       CASE ( 3 )
          CONTINUE
!
!--    pt-component
       CASE ( 4 )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                k_topo = topo_top_ind(j,i,0)
                DO  ik = 1, nz_meso(j,i)
                   k = k_topo + ik
                   tend(k,j,i) = tend(k,j,i) + fr_urb(j,i) * tt_urb(ik,j,i)
                ENDDO
             ENDDO
          ENDDO
!
!--    Humidity-component (not implemented yet)
       CASE ( 5 )
          CONTINUE
!
!--    Turbulence (not implemented yet)
       CASE ( 6 )
          CONTINUE
!
!--    Default case
       CASE DEFAULT
          WRITE( message_string, '(A,I7)' )                                                        &
                                'wrong component or not implemented yet for component: ', component
          CALL message( 'dcep_tendency', 'DCP0008', 1, 2, 0, 6, 0 )

    END SELECT

 END SUBROUTINE dcep_tendency


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate tendency due to urban area consideration using DCEP (Cache-optimized version).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_tendency_ji( i, j, component )

    IMPLICIT NONE

    INTEGER(iwp) ::  component    !< prognostic variable (???)
    INTEGER(iwp) ::  i            !< loop indices
    INTEGER(iwp) ::  ik           !< loop indices
    INTEGER(iwp) ::  j            !< loop indices
    INTEGER(iwp) ::  k            !< loop indices
    INTEGER(iwp) ::  k_topo       !< topography top index


    k_topo = topo_top_ind(j,i,0)

    SELECT CASE ( component )
!
!--    u-component
       CASE ( 1 )
          IF ( lurbvel )  THEN
             DO  ik = 1, nz_meso(j,i)
                k = k_topo + ik
                tend(k,j,i) = tend(k,j,i) + fr_urb(j,i) * ut_urb(ik,j,i)
             ENDDO
          ENDIF
!
!--    v-component
       CASE ( 2 )
          IF ( lurbvel )  THEN
             DO  ik = 1, nz_meso(j,i)
                k = k_topo + ik
                tend(k,j,i) = tend(k,j,i) + fr_urb(j,i) * vt_urb(ik,j,i)
             ENDDO
          ENDIF
!
!--    w-component (not implemented yet)
       CASE ( 3 )
          CONTINUE
!
!--    pt-component
       CASE ( 4 )
          DO  ik = 1, nz_meso(j,i)
             k = k_topo + ik
             tend(k,j,i) = tend(k,j,i) + fr_urb(j,i) * tt_urb(ik,j,i)
          ENDDO
!
!--    Humidity-component (not implemented yet)
       CASE ( 5 )
          CONTINUE
!
!--    Turbulence (not implemented yet)
       CASE ( 6 )
          CONTINUE
!
!--    Default case
       CASE DEFAULT
          WRITE( message_string, '(A,I7)' )                                                        &
                                'wrong component or not implemented yet for component: ', component
          CALL message( 'dcep_tendency_ji', 'DCP0008', 1, 2, 0, 6, 0 )

    END SELECT

  END SUBROUTINE dcep_tendency_ji


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Input dcep related inputs from a NetCDF file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_netcdf_input

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               ny

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  check_existence,                                                                    &
               close_input_file,                                                                   &
               get_dimension_length,                                                               &
               get_attribute,                                                                      &
               get_variable,                                                                       &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               open_read_file

    USE control_parameters,                                                                        &
        ONLY:  coupling_char

    IMPLICIT NONE

    CHARACTER(LEN=100) ::  input_file_dcep = 'PIDS_DCEP'         !< Name of file which comprises dcep input data
    CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< variable names

    LOGICAL ::  netcdf_dcep = .FALSE.  !< Flag: netcdf file exists

    INTEGER(iwp) ::  ke_uhl_f    !< urban leveles from PIDS_DCEP
    INTEGER(iwp) ::  id_dcep     !< NetCDF id of PIDS_DCEP
    INTEGER(iwp) ::  n_udir_f    !< no. street directions of PIDS_DCEP
    INTEGER(iwp) ::  num_vars    !< number of variables
    INTEGER(iwp) ::  nx_f        !< size of x direction of PIDS_DCEP
    INTEGER(iwp) ::  ny_f        !< size of y direction of PIDS_DCEP

    IF ( debug_output )  CALL debug_message( 'dcep_netcdf_input', 'start' )

    INQUIRE( FILE = TRIM( input_file_dcep ) //  TRIM( coupling_char ), EXIST = netcdf_dcep )

    IF ( netcdf_dcep )  THEN
!
!--    Open file in read-only mode.
       CALL open_read_file( TRIM( input_file_dcep ) // TRIM( coupling_char ), id_dcep )
!
!--    Inquire all variable names.
       CALL inquire_num_variables( id_dcep, num_vars )
!
!--    Check for correct horizontal and vertical dimension
       CALL get_dimension_length( id_dcep, nx_f, 'x'  )
       CALL get_dimension_length( id_dcep, ny_f, 'y'  )
       CALL get_dimension_length( id_dcep, ke_uhl_f, 'uheight1'  )
       CALL get_dimension_length( id_dcep, n_udir_f, 'streetdir'  )

       IF ( nx_f-1 /= nx  .OR.  ny_f-1 /= ny  .OR.                                                 &
            ke_uhl_f-1 /= ke_uhl  .OR.  n_udir_f /= n_udir_f )  THEN
          message_string = 'one or more dimension of in DCEP input file does not match the ' //    &
                           'inputs in parin file'
          CALL message( 'dcep_netcdf_input', 'DCP0009', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Allocate memory to store variable names.
       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_dcep, var_names )
!
!--    Read the street direction (ang_udir).
       IF ( check_existence( var_names, 'streetdir' ) )                                            &
            CALL get_variable( id_dcep, 'streetdir', ang_udir(1:n_udir) )
!
!--    Read the fraction of urban parts in a grid element (fr_urb).
       IF ( check_existence( var_names, 'FR_URBAN' ) )  THEN
          CALL get_variable( id_dcep, 'FR_URBAN', fr_urb, nxl, nxr, nys, nyn )
       ELSE
          WRITE( message_string, * ) 'missing urban fraction (FR_URBAN) in ' //                    &
                                     TRIM( input_file_dcep )
          CALL message( 'dcep_netcdf_input', 'DCP0010', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read the fraction of street direction (fr_udir).
       IF ( check_existence( var_names, 'FR_STREETD' ) )  THEN
          CALL get_variable( id_dcep, 'FR_STREETD', fr_udir, nxl, nxr, nys, nyn, 0, n_udir-1, 0,   &
                             n_uclass-1 )
       ELSE
          WRITE( message_string, * ) 'missing the fraction of street direction (FR_STREETD) in ' //&
                                     TRIM( input_file_dcep )
          CALL message( 'dcep_netcdf_input', 'DCP0010', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read the street width (width_street).
       IF ( check_existence( var_names, 'STREET_W' ) )  THEN
          CALL get_variable( id_dcep, 'STREET_W', width_street, nxl, nxr, nys, nyn, 0, n_udir-1, 0,&
                             n_uclass-1 )
       ELSE
          WRITE( message_string, * ) 'missing the street width (STREET_W) in ' //                  &
                                     TRIM( input_file_dcep )
          CALL message( 'dcep_netcdf_input', 'DCP0010', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read the building width (width_build).
       IF ( check_existence( var_names, 'BUILD_W' ) )  THEN
          CALL get_variable( id_dcep, 'BUILD_W', width_build, nxl, nxr, nys, nyn, 0, n_udir-1, 0,  &
                             n_uclass-1 )
       ELSE
          WRITE( message_string, * ) 'missing the building width (BUILD_W) in ' //                 &
                                     TRIM( input_file_dcep )
          CALL message( 'dcep_netcdf_input', 'DCP0010', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read the building probability (fr_roof).
       IF ( check_existence( var_names, 'BUILD_PROP' ) )  THEN
          CALL get_variable( id_dcep, 'BUILD_PROP', fr_roof, nxl, nxr, nys, nyn, 0, ke_uhl, 0,     &
                             n_udir-1, 0, n_uclass-1 )
       ELSE
          WRITE( message_string, * ) 'missing the building width (BUILD_PROP) in ' //              &
                                     TRIM( input_file_dcep )
          CALL message( 'dcep_netcdf_input', 'DCP0010', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Read the fraction of urban class in a grid element (fr_uclass).
       IF ( check_existence( var_names, 'FR_URBANCL' ) )  THEN
          CALL get_variable( id_dcep, 'FR_URBANCL', fr_uclass, nxl, nxr, nys, nyn, 0, n_uclass-1 )
       ELSE
          WRITE( message_string, * ) 'missing fraction of urban class (FR_URBANCL) in ' //         &
                                     TRIM( input_file_dcep )
          CALL message( 'dcep_netcdf_input', 'DCP0010', 1, 2, 0, 6, 0 )
       ENDIF

    ELSE

       message_string = 'input file "'// TRIM( input_file_dcep ) // TRIM( coupling_char ) //       &
                        '" for DCEP PIDS missing'
       CALL message( 'dcep_netcdf_input', 'DCP0011', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Close input file.
    CALL close_input_file( id_dcep )

    IF ( debug_output )  CALL debug_message( 'dcep_netcdf_input', 'end' )

 END SUBROUTINE dcep_netcdf_input


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to check data output for DCEP model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_check_data_output( variable, unit )

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  unit            !< unit of the variable
    CHARACTER(LEN=*) ::  variable        !< variable name
    CHARACTER(LEN=varnamelength) ::  var !< TRIM(variable)


    var = TRIM( variable )
!
!-- Check if variable is related to DCEP variables and DCEP is ON.
    IF ( var(1:5) == 'dcep_'  )  THEN
       IF ( .NOT. dcep )  THEN
          message_string = 'output of "' // TRIM( var ) // '" requires dcep = .TRUE.'
          CALL message( 'dcep_check_data_output', 'DCP0012', 1, 2, 0, 6, 0 )
       ENDIF
    ELSE
       unit = 'illegal'
       RETURN
    ENDIF
!
!-- Search for the variable.
    IF ( var  == 'dcep_nz_meso*'     ) unit = '1'
    IF ( var  == 'dcep_fr_wall'      ) unit = '1'
    IF ( var  == 'dcep_fr_roof'      ) unit = '1'
    IF ( var  == 'dcep_tt_urb'       ) unit = 'K s -1'
    IF ( var  == 'dcep_shfl_roof'    ) unit = 'W m -2'
    IF ( var  == 'dcep_sw_roof'      ) unit = 'W m -2'
    IF ( var  == 'dcep_t_g_urb*'     ) unit = 'K'
    IF ( var  == 'dcep_shfl_urb*'    ) unit = 'W m -2'
    IF ( var  == 'dcep_albedo_urb*'  ) unit = '1'
    IF ( var  == 'dcep_albedop_urb*' ) unit = '1'
    IF ( var  == 'dcep_emiss_urb*'   ) unit = '1'
    IF ( var  == 'dcep_t_grad_urb*'  ) unit = 'K'
    IF ( var  == 'dcep_rl_roof'      ) unit = 'W m -2'
    IF ( var  == 'dcep_rl_wallw'     ) unit = 'W m -2'
    IF ( var  == 'dcep_rl_walle'     ) unit = 'W m -2'
    IF ( var  == 'dcep_strfl_urb*'   ) unit = 'W m -2'
    IF ( var  == 'dcep_t_ground'     ) unit = 'K'
    IF ( var  == 'dcep_t_roof1'      ) unit = 'K'
    IF ( var  == 'dcep_t_walle'      ) unit = 'K'
    IF ( var  == 'dcep_t_wallw'      ) unit = 'K'

 END SUBROUTINE dcep_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 2D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )

    IMPLICIT NONE
!
!-- Input variables
    CHARACTER(LEN=*), INTENT(IN)  ::  variable !< Char identifier to select var for output

    INTEGER(iwp), INTENT(IN)      ::  av       !< Use averaged data? 0 = no, 1 = yes?
    INTEGER(iwp), INTENT(IN)      ::  nzb_do   !< Unused. 2D. nz bottom to nz top
    INTEGER(iwp), INTENT(IN)      ::  nzt_do   !< Unused.
!
!-- Output variables
    CHARACTER(LEN=*), INTENT(OUT)  ::  grid    !< Grid type (always "zu1" for biom)
    CHARACTER(LEN=*), INTENT(IN)   ::  mode    !<
    LOGICAL, INTENT(OUT)           ::  found   !< Output found?
    LOGICAL, INTENT(OUT)           ::  two_d   !< Flag parameter that indicates 2D variables,
!
!-- Horizontal cross sections, must be .TRUE. for
    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !< result grid to return
!
!-- Internal variables
    INTEGER(iwp) ::  i  !< Running index, x-dir
    INTEGER(iwp) ::  j  !< Running index, y-dir
    INTEGER(iwp) ::  k  !< Running index, z-dir


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

       CASE ( 'dcep_nz_meso*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = 1, 2
                      local_pf(i,j,nzb+1) = MAXVAL( width_build(:,:,j,i) )
                   ENDDO
                ENDDO
             ENDDO
             two_d = .TRUE.
             grid = 'zu1'
          ELSE
             found = .TRUE.
             grid  = 'none'
          ENDIF

       CASE ( 'dcep_t_g_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,nzb+1) = t_g_urb(j,i)
                ENDDO
             ENDDO
             two_d = .TRUE.
             IF ( mode == 'xy' )  grid = 'zw'
          ELSE
             found = .TRUE.
             grid  = 'none'
          ENDIF

       CASE ( 'dcep_shfl_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,nzb+1) = -shfl_urb(j,i)
                ENDDO
             ENDDO
             two_d = .TRUE.
             IF ( mode == 'xy' )  grid = 'zw'
          ELSE
             found = .TRUE.
             grid  = 'none'
          ENDIF

       CASE ( 'dcep_albedo_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,nzb+1) = albedo_urb(j,i)
                ENDDO
             ENDDO
             two_d = .TRUE.
             IF ( mode == 'xy' )  grid = 'zw'
          ELSE
             found = .TRUE.
             grid  = 'none'
          ENDIF

       CASE ( 'dcep_albedop_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,nzb+1) = albedop_urb(j,i)
                ENDDO
             ENDDO
             two_d = .TRUE.
             IF ( mode == 'xy' )  grid = 'zw'
          ELSE
             found = .TRUE.
             grid  = 'none'
          ENDIF

       CASE ( 'dcep_emiss_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,nzb+1) = emiss_urb(j,i)
                ENDDO
             ENDDO
             two_d = .TRUE.
             IF ( mode == 'xy' )  grid = 'zw'
          ELSE
             found = .TRUE.
             grid  = 'none'
          ENDIF

       CASE ( 'dcep_t_grad_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,nzb+1) = t_grad_urb(j,i)
                ENDDO
             ENDDO
             two_d = .TRUE.
             IF ( mode == 'xy' )  grid = 'zw'
          ELSE
             found = .TRUE.
             grid  = 'none'
          ENDIF

       CASE ( 'dcep_strfl_urb*_xy' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,nzb+1) = -strfl_urb(j,i)
                ENDDO
             ENDDO
             two_d = .TRUE.
             IF ( mode == 'xy' )  grid = 'zw'
          ELSE
             found = .TRUE.
             grid  = 'none'
          ENDIF

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT

 END SUBROUTINE dcep_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 3D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )

    USE indices

    USE arrays_3d,                                                                                 &
        ONLY: tpt_m

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)    ::  variable  !< variable name

    INTEGER(iwp),     INTENT(IN)    ::  av        !< flag for (non-)average output
    INTEGER(iwp),     INTENT(IN)    ::  nzb_do    !< vertical output index (bottom) (usually 0)
    INTEGER(iwp),     INTENT(IN)    ::  nzt_do    !< vertical output index (top) (usually nz_do3d)

    LOGICAL,          INTENT(INOUT) ::  found     !< flag if output variable is found

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf  !<
!
!-- Internal variables
    CHARACTER(len=varnamelength)  :: var  !<

    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  id    !<
    INTEGER(iwp) ::  iurb  !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  k     !<


    var = TRIM( variable )
!
!-- Return for non-DCEP variables.
    IF ( var(1:5) /= 'dcep_' )  THEN
       found = .FALSE.
       RETURN
    ENDIF
!
!-- Stop if average quantities are required (not implemented).
    IF ( av  == 1  )  THEN
       message_string = 'average of "'// TRIM(variable) // '" is not implemented'
       CALL message( 'dcep_data_output_3d', 'DCP0013', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Initialize.
    found = .TRUE.
!
!-- Select variable.
    SELECT CASE ( var )

       CASE ( 'dcep_fr_roof' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_uhl
                   local_pf(i,j,k) = fr_roof(1,1,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_fr_wall' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_uhl
                   local_pf(i,j,k) = fr_wall(1,1,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_tt_urb' )
          DO  i = nxl, nxr
             DO  j = nys, nyn
                IF ( fr_urb(j,i) > eps_urb )  THEN
                   DO  k = 1, nz_mesom
                      local_pf(i,j,k) = tt_urb(k,j,i)
                   ENDDO
                ELSE
                   DO  k = 1, nz_mesom
                      local_pf(i,j,k) = tpt_m(k,j,i)
                   ENDDO
                ENDIF
             ENDDO
          ENDDO

       CASE ( 'dcep_shfl_roof' )
          id = 1
          iurb = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_uhl
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,k) = shfl_roof(iurb,id,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_rl_roof' )
          id = 1
          iurb = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_uhl
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,k) = rl_roof(iurb,id,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_sw_roof' )
          id = 1
          iurb = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_uhl
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,k) = rs_roof(iurb,id,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_rl_wallw' )
          id = 1
          iurb = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_uhl
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,k) = rl_wall(iurb,2*id-1,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_rl_walle' )
          id = 1
          iurb = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_uhl
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,k) = rl_wall(iurb,2*id,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_t_ground' )
          id = 1
          iurb = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_ground
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,k) = t_ground(iurb,id,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_t_roof1' )
          id = 1
          iurb = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_roof
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,k) = t_roof(iurb,id,1,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_t_walle' )
          id = 1
          iurb = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_wall
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,k) = t_wall(iurb,2*id,1,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE ( 'dcep_t_wallw' )
          id = 1
          iurb = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                DO  k = 1, ke_wall
                   IF ( fr_urb(j,i) > eps_urb )  local_pf(i,j,k) = t_wall(iurb,2*id-1,1,k,j,i)
                ENDDO
             ENDDO
          ENDDO

       CASE DEFAULT
          found = .FALSE.

    END SELECT

 END SUBROUTINE dcep_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dcep_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN)  ::  variable    !<

    LOGICAL, INTENT(OUT)           ::  found       !<

    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x      !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y      !<
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z      !<
!
!-- Internal variables
    CHARACTER(len=varnamelength)  :: var  !<


    found  = .TRUE.
!
!-- Check for the grid.
    var = TRIM( variable )

    SELECT CASE ( TRIM( var ) )

       CASE ( 'dcep_nz_meso*_xy' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE ( 'dcep_fr_wall' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE ( 'dcep_fr_roof' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zuhl'

       CASE ( 'dcep_tt_urb' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zmeso'

       CASE ( 'dcep_shfl_roof' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zuhl'

       CASE ( 'dcep_rl_roof', 'dcep_rl_wallw', 'dcep_rl_walle' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zuhl'

       CASE ( 'dcep_t_ground' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zground'

       CASE ( 'dcep_t_roof1' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zroof'

       CASE ( 'dcep_t_walle', 'dcep_t_wallw' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zwall'

       CASE ( 'dcep_sw_roof' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zuhl'

       CASE ( 'dcep_t_g_urb*_xy', 'dcep_shfl_urb*_xy', 'dcep_albedo_urb*_xy',                      &
              'dcep_albedop_urb*_xy', 'dcep_emiss_urb*_xy', 'dcep_t_grad_urb*_xy',                 &
              'dcep_strfl_urb*_xy')
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE DEFAULT
           found  = .FALSE.
           grid_x = 'none'
           grid_y = 'none'
           grid_z = 'none'

    END SELECT

 END SUBROUTINE dcep_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> View_factors calculate the single canyon view factors and view factors with roofs.
!--------------------------------------------------------------------------------------------------!
 PURE SUBROUTINE view_factors( nz, z, ws, fww, fwg, fgw, fsg, fgs, fsw, fws, frw, fwr, fsr, frs,   &
                                    bs )

    IMPLICIT NONE
!
!-- Input variables
    INTEGER(iwp), INTENT(in) ::  nz          !< number of levels in the urban grid

    REAL(wp), INTENT(in)     ::  ws          !< street width
    REAL(wp), INTENT(in)     ::  z(ke_uhl+1) !< height of the urban grid levels
    REAL(wp), INTENT(in), OPTIONAL ::  bs    !< building width
!
!-- Output variables.
!-  fww,fwg,fgw,fsw,fsg are the view factors used to compute the long wave
!-  and the short wave radation. They are the part of radiation from a surface
!-  or from the sky to another surface.
    REAL(wp), INTENT(out) ::  fgs                 !< from ground to sky
    REAL(wp), INTENT(out) ::  fsg                 !< from sky to ground

    REAL(wp), INTENT(out) ::  fgw(ke_uhl)         !< from ground to wall
    REAL(wp), INTENT(out) ::  fsw(ke_uhl)         !< from sky to wall
    REAL(wp), INTENT(out) ::  fwg(ke_uhl)         !< from wall to ground
    REAL(wp), INTENT(out) ::  fws(ke_uhl)         !< from wall to sky

    REAL(wp), INTENT(out) ::  fww(ke_uhl,ke_uhl)  !< from wall to wall
!
!-  Optional view factor with roofs
    REAL(wp), INTENT(out), OPTIONAL ::  frs(ke_uhl+1)          !< roof to sky
    REAL(wp), INTENT(out), OPTIONAL ::  fsr(ke_uhl+1)          !< sky to roof

    REAL(wp), INTENT(out), OPTIONAL ::  frw(ke_uhl+1,ke_uhl)   !< roof to wall
    REAL(wp), INTENT(out), OPTIONAL ::  fwr(ke_uhl,ke_uhl+1)   !< wall to roof
!
!-- Internal variables
    INTEGER(iwp) ::  iz  !<
    INTEGER(iwp) ::  jz  !<

    REAL(wp) ::  a1     !<
    REAL(wp) ::  a12    !<
    REAL(wp) ::  a123   !<
    REAL(wp) ::  a2     !<
    REAL(wp) ::  a23    !<
    REAL(wp) ::  a3     !<
    REAL(wp) ::  a34    !<
    REAL(wp) ::  a4     !<
    REAL(wp) ::  a5     !<
    REAL(wp) ::  fprl   !<
    REAL(wp) ::  ftot   !<
    REAL(wp) ::  f1     !<
    REAL(wp) ::  f12    !<
    REAL(wp) ::  f123   !<
    REAL(wp) ::  f1234  !<
    REAL(wp) ::  f1245  !<
    REAL(wp) ::  f14    !<
    REAL(wp) ::  f2     !<
    REAL(wp) ::  f23    !<
    REAL(wp) ::  f234   !<
    REAL(wp) ::  hut    !<

!
!-- Highest height: this is where the sky is.
    hut = z(nz+1)
!
!-- Loop over receiving surface.
    DO  jz = 1, nz
!
!--    Radiation from wall to wall.
!--    Loop over sending surface.
       DO  iz = 1, nz
          f123 = fprls( ABS( z(jz+1) - z(iz)   ), ws )
          f23  = fprls( ABS( z(jz+1) - z(iz+1) ), ws )
          f12  = fprls( ABS( z(jz)   - z(iz)   ), ws )
          f2   = fprls( ABS( z(jz)   - z(iz+1) ), ws )
!
!--       Area segments on wall, street length cancels anyway.
          a123 = ABS( z(jz+1) - z(iz)   )
          a12  = ABS( z(jz)   - z(iz)   )
          a23  = ABS( z(jz+1) - z(iz+1) )
          a1   = ABS( z(iz+1) - z(iz)   )
          a2   = ABS( z(jz)   - z(iz+1) )
          a3   = ABS( z(jz+1) - z(jz)   )
!
!--       Sky view factor.
          ftot = 0.5_wp * ( a123 * f123 - a23  * f23 - a12 * f12 + a2 * f2 ) !/a1
!
!--       Effective svf taking sending and receiving area into account.
          fww(iz,jz) = ftot / a3 !*a1
       ENDDO

       DO  iz = nz+1, ke_uhl
          fww(iz,jz) = 0.0_wp
       ENDDO

       IF ( lrroofs )  THEN
!
!--       Radiation from wall to roof and roof to wall.
          DO  iz = 1, jz-1
             fwr(iz,jz) = 0.0_wp
             frw(jz,iz) = 0.0_wp
          ENDDO

          DO  iz = jz, nz
             a1  = z(iz+1) - z(iz)
             a12 = z(iz+1) - z(jz)
             a2  = z(iz)   - z(jz)
             a3  = ws
             a34 = ws + bs
             a4  = bs
             f1234 = fnrms( a34, a12 )
             f123  = fnrms(  a3, a12 )
             f234  = fnrms( a34,  a2 )
             f23   = fnrms(  a3,  a2 )
             ftot = ( a12 * ( f1234 - f123 ) + a2 * ( f23 - f234 ) )
             fwr(iz,jz) = ftot / a4
             frw(jz,iz) = ftot / a1
          ENDDO

          DO  iz = nz+1, ke_uhl
             fwr(iz,jz) = 0.0_wp
             fwr(jz,iz) = 0.0_wp
          ENDDO
!
!--       Radiation from sky to roof.
          a1   = ws
          a12  = ws + bs
          a5   = bs
          a123 = bs + 2.0_wp * ws
          f1245 = fprls( a12, hut-z(jz) )
          f14   = fprls( a1,  hut-z(jz) )
          ftot = a12 * f1245 - a1 * f14
          fsr(jz) = ftot / a5
          frs(jz) = ftot / a123
       ENDIF
!
!--    Radiation from ground to wall.
       f12 = fnrms( z(jz+1), ws )
       f1  = fnrms(   z(jz), ws )
       a1 = ws
       a4 = z(jz+1) - z(jz)
       ftot = f12 - f1
       fgw(jz) = ftot * a1 / a4
!
!--    Radiation from sky to wall.
       f12 = fnrms(   hut-z(jz), ws )
       f1  = fnrms( hut-z(jz+1), ws )
       a1 = ws
       a4 = z(jz+1) - z(jz)
       ftot = f12 - f1
       fsw(jz) = ftot * a1 / a4
       fws(jz) = ftot
    ENDDO
!
!-- Radiation from wall to ground.
    DO  iz = 1, nz
       f12 = fnrms( ws, z(iz+1) )
       f1  = fnrms( ws, z(iz)   )
       a2 = z(iz)
       a12= z(iz+1)
       a4 = ws
       ftot = (a12 * f12 - a2 * f1)   !/a1
       fwg(iz) = ftot / a4        !*a1
    ENDDO
!
!-- Radiation from sky to ground.
    fprl = fprls( ws, hut )
    fsg = fprl
    fgs = fprl

    DO  jz = nz+1, ke_uhl
       DO  iz = 1, ke_uhl
          fww(iz,jz) = 0.0_wp
       ENDDO
       fgw(jz) = 0.0_wp
       fsw(jz) = 0.0_wp
       fwg(jz) = 0.0_wp
    ENDDO

    IF ( lrroofs )  THEN
       fsr(nz+1) = 1.0_wp
       frs(nz+1) = bs / ( bs + 2.0_wp * ws )
       DO  jz = nz+2, ke_uhl+1
          fsr(jz) = 0.0_wp
          frs(jz) = 0.0_wp
       ENDDO
       DO  jz = nz+1, ke_uhl+1
          DO  iz = 1, ke_uhl
             fwr(iz,jz) = 0.0_wp
             frw(jz,iz) = 0.0_wp
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE view_factors


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> fprls calculates view factor of two parallel surfaces. One dimension uses
!> the limit infinity.
!>
!>      /|        /|
!>  inf/ |    inf/ |
!>    /  |      /  |
!>    |  |      |  |
!>   a|  /     a|  /
!>    | /       | /
!>    |/        |/
!>         c
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION fprls( a, c )

    IMPLICIT NONE

    REAL(wp) ::  fprls  !<

    REAL(wp), INTENT(in) ::  a  !<
    REAL(wp), INTENT(in) ::  c  !<

    REAL(wp) ::  y  !<

    IF ( c <= azero )  THEN
       fprls = 1.0_wp
    ELSE
       y = a / c
       IF ( y <= azero )  THEN
          fprls = 0.0_wp
       ELSE
          fprls = ( SQRT( 1.0_wp + y**2 ) - 1.0_wp ) / y
       ENDIF
    ENDIF

 END FUNCTION fprls


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> fnrms calculates view factor of two normal surfacse. One dimension uses
!> the limit infinity.
!>      /|
!>  inf/ |
!>    /  |
!>    |  |________
!>   a|  /       /
!>    | /       /
!>    |/_______/
!>         c
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION fnrms( a, c )

    IMPLICIT NONE

    REAL(wp) ::  fnrms  !<

    REAL(wp), INTENT(in) ::  a  !<
    REAL(wp), INTENT(in) ::  c  !<
!
!-- Internal variables.
    REAL(wp) ::  ratio  !<


    IF ( c <= azero )  THEN
       fnrms = 0.5_wp
    ELSE
       ratio = a / c
       fnrms = 0.5_wp * ( 1.0_wp + ratio - SQRT( 1 + ratio**2 ) )
    ENDIF

 END FUNCTION fnrms


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> fprls_da calculates view factor of two parallel surfaces, one with
!> infinitesimal height. One dimension uses the limit infinity.
!>              /|
!>          inf/ |
!>            /  |
!>     /      |  |
!> inf//     a|  /
!>   //       | /
!>   /        |/
!>         c
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION fprls_da( a, c )

    IMPLICIT NONE

    REAL(wp) ::  fprls_da  !<

    REAL(wp), INTENT(in) ::  a  !<
    REAL(wp), INTENT(in) ::  c  !<


    fprls_da = 0.5_wp * a / SQRT( a**2 + c**2 )

 END FUNCTION fprls_da


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> fnrms_da calculates view factor of two normal surfaces, one with
!> infinitesimal height. One dimension uses the limit infinity.
!>       /
!>   inf//
!>     //
!>     /  ________
!>   a   /       /
!>      /       /
!>     /_______/
!>         c
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION fnrms_da( a, c )

    IMPLICIT NONE

    REAL(wp) ::  fnrms_da  !<
!
!-- Internal variables
    REAL(wp), INTENT(in) ::  a  !<
    REAL(wp), INTENT(in) ::  c  !<

    fnrms_da = 0.5_wp * ( 1.0_wp - c / SQRT( a**2 + c**2 ) )

 END FUNCTION fnrms_da


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation function 
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION fnrm14( z1, z2, h1, h2 )

    IMPLICIT NONE

    REAL(wp) ::  fnrm14  !<

    REAL(wp), INTENT(in) ::  h1  !<
    REAL(wp), INTENT(in) ::  h2  !<
    REAL(wp), INTENT(in) ::  z1  !<
    REAL(wp), INTENT(in) ::  z2  !<
!
!-- Internal variables
    REAL(wp) ::  t1  !<
    REAL(wp) ::  t2  !<
    REAL(wp) ::  t3  !<
    REAL(wp) ::  t4  !<


    t1 = fnrms( z2, h2 )
    t2 = fnrms( z1, h2 )
    t3 = fnrms( z1, h1 )
    t4 = fnrms( z2, h1 )

    fnrm14 =  h2 * ( t1 - t2 ) + h1 * ( t3 - t4 )

 END FUNCTION fnrm14


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> help function
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION fprl16( s1, s2, r1, r2, d )

    IMPLICIT NONE

    REAL(wp) ::  fprl16  !<

    REAL(wp), INTENT(in) ::  d   !<
    REAL(wp), INTENT(in) ::  r1  !<
    REAL(wp), INTENT(in) ::  r2  !<
    REAL(wp), INTENT(in) ::  s1  !<
    REAL(wp), INTENT(in) ::  s2  !<
!
!-- Internal variables
    REAL(wp) ::  t1  !<
    REAL(wp) ::  t2  !<
    REAL(wp) ::  t3  !<
    REAL(wp) ::  t4  !<


    t1 = fprls( ABS( r2 - s1 ), d )
    t2 = fprls( ABS( r2 - s2 ), d )
    t3 = fprls( ABS( r1 - s1 ), d )
    t4 = fprls( ABS( r1 - s2 ), d )

    fprl16 =  0.5_wp * ( ABS( r2 - s1 ) * t1 - ABS( r2 - s2 ) * t2 - ABS( r1 - s1 ) * t3 +         &
                         ABS( r1 - s2 ) * t4 )

 END FUNCTION fprl16


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> help function
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION fprl134( g, s, h )

    IMPLICIT NONE

    REAL(wp), INTENT(in) ::  g  !<
    REAL(wp), INTENT(in) ::  h  !<
    REAL(wp), INTENT(in) ::  s  !<
!
!-- Internal variables
    REAL(wp) ::  fprl134  !<
    REAL(wp) ::  t1       !<
    REAL(wp) ::  t2       !<
    REAL(wp) ::  t3       !<


    t1 = fprls(   s, h )
    t2 = fprls(   g, h )
    t3 = fprls( s-g, h )

    fprl134 = 0.5_wp * ( s * t1 + g * t2 - (s - g) * t3 )

 END FUNCTION fprl134


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> help function
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION fnrm13( h1, h2, z )

    IMPLICIT NONE

    REAL(wp) ::  fnrm13  !<

    REAL(wp), INTENT(in) ::  h1  !<
    REAL(wp), INTENT(in) ::  h2  !<
    REAL(wp), INTENT(in) ::  z   !<
!
!-- Internal variables
    REAL(wp) ::  t1  !<
    REAL(wp) ::  t2  !<


    t1 = fnrms( z, h1 )
    t2 = fnrms( z, h2 )

    fnrm13 = h1 * t1 - h2 * t2

 END FUNCTION fnrm13


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate emitted radiation into the sky using energy conservation.
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION uprad( r, r_ground, r_wall, r_roof, albp_ground, albp_wall, albp_roof, heat )

    IMPLICIT NONE

    REAL(wp) ::  uprad(nys:nyn,nxl:nxr)  !< function type

    LOGICAL, INTENT(in)  ::  heat  !<

    REAL(wp), INTENT(in) ::  albp_ground(n_uclass)                             !<
    REAL(wp), INTENT(in) ::  albp_roof(n_uclass)                               !<
    REAL(wp), INTENT(in) ::  albp_wall(n_uclass)                               !<

    REAL(wp), INTENT(in) ::  r(nys:nyn,nxl:nxr)                                !<

    REAL(wp), INTENT(in) ::  r_ground(n_uclass,n_udir,nys:nyn,nxl:nxr)         !<

    REAL(wp), INTENT(in) ::  r_roof(n_uclass,n_udir,ke_uhl+1,nys:nyn,nxl:nxr)  !<
    REAL(wp), INTENT(in) ::  r_wall(n_uclass,2*n_udir,ke_uhl,nys:nyn,nxl:nxr)  !<
!
!-- Internal variables
    INTEGER ::  i     !<
    INTEGER ::  id    !<
    INTEGER ::  iu    !<
    INTEGER ::  iurb  !<
    INTEGER ::  j     !<

    REAL(wp) ::  r_emit_c  !<
    REAL(wp) ::  r_emit_d  !<
    REAL(wp) ::  r_inc_c   !<
    REAL(wp) ::  r_inc_d   !<
    REAL(wp) ::  r_up_c    !<
    REAL(wp) ::  r_up_d    !<


    DO  i = nxl, nxr
       DO  j = nys, nyn

          uprad(j,i) = 0.0_wp

          DO  iurb = 1, n_uclass
             IF ( fr_uclass(iurb,j,i) >= eps_urb )  THEN
!
!--             Upward longwave radiation for this street direction.
                r_up_c   = 0.0_wp
                r_inc_c  = 0.0_wp
                r_emit_c = 0.0_wp
!
!--             Number of street directions for urban class.
                DO  id = 1, n_udir
                   IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN
                      r_up_d = 0.0_wp
!
!--                   Ground:
                      r_emit_d = albp_ground(iurb) * r_ground(iurb,id,j,i) *                       &
                                 width_street(iurb,id,j,i)
                      r_inc_d  = r_ground(iurb,id,j,i) * width_street(iurb,id,j,i)
!
!--                   Roof:
                      DO  iu = 1, nz(iurb,id,j,i)+1
                         r_emit_d = r_emit_d + albp_roof(iurb) * r_roof(iurb,id,iu,j,i) *          &
                                               fr_roof(iurb,id,iu,j,i) * width_build(iurb,id,j,i)
                         r_inc_d = r_inc_d + r_roof(iurb,id,iu,j,i) * fr_roof(iurb,id,iu,j,i) *    &
                                             width_build(iurb,id,j,i)
                      ENDDO
!
!--                   Wall:
                      DO  iu = 1, nz(iurb,id,j,i)
                         r_emit_d = r_emit_d + albp_wall(iurb) *                                   &
                                       ( r_wall(iurb,2*id-1,iu,j,i) + r_wall(iurb,2*id,iu,j,i) ) * &
                                               dz_uhl(iu) * fr_wall(iurb,id,iu,j,i)
                         r_inc_d = r_inc_d +                                                       &
                                       ( r_wall(iurb,2*id-1,iu,j,i) + r_wall(iurb,2*id,iu,j,i) ) * &
                                       dz_uhl(iu) * fr_wall(iurb,id,iu,j,i)
                      ENDDO
                      r_up_c = r_up_c + r_up_d / width_sglcan(iurb,id,j,i) * fr_udir(iurb,id,j,i)
                      r_inc_c = r_inc_c + r_inc_d / width_sglcan(iurb,id,j,i) * fr_udir(iurb,id,j,i)
                      r_emit_c = r_emit_c +                                                        &
                                 r_emit_d / width_sglcan(iurb,id,j,i) * fr_udir(iurb,id,j,i)
                   ENDIF
                ENDDO
                r_up_c = r(j,i) - ( r_inc_c - r_emit_c )

                uprad(j,i) = uprad(j,i) + r_up_c * fr_uclass(iurb,j,i)

             ENDIF
          ENDDO

       ENDDO
    ENDDO

    IF ( heat )  THEN

       DO  i = nxl, nxr
          DO  j = nys, nyn

             DO  iurb = 1, n_uclass

                IF ( fr_uclass(iurb,j,i) >= eps_urb )  THEN

                   r_emit_c = 0.0_wp
!
!--                Number of street directions for urban class.
                   DO  id = 1, n_udir
                      IF ( fr_udir(iurb,id,j,i) > eps_urb )  THEN

                         r_emit_d = ( 1.0_wp - albp_ground(iurb) ) * sigma_sb *                    &
                                    t_ground(iurb,id,ke_ground,j,i)**4.0_wp *                      &
                                    width_street(iurb,id,j,i)
!
!--                      Roof:
                         DO  iu = 1, nz(iurb,id,j,i)+1
                            r_emit_d = r_emit_d + ( 1.0_wp - albp_roof(iurb) ) * sigma_sb *        &
                                                  t_roof(iurb,id,ke_roof,iu,j,i)**4.0_wp *         &
                                                  fr_roof(iurb,id,iu,j,i) * width_build(iurb,id,j,i)
                         ENDDO
!
!--                      Wall:
                         DO  iu = 1, nz(iurb,id,j,i)
                            r_emit_d = r_emit_d + ( 1.0_wp-albp_wall(iurb) ) * sigma_sb *          &
                                                  ( t_wall(iurb,2*id-1,ke_wall,iu,j,i)**4.0_wp +   &
                                                    t_wall(iurb,2*id,ke_wall,iu,j,i)**4.0_wp ) *   &
                                                  dz_uhl(iu) * fr_wall(iurb,id,iu,j,i)
                         ENDDO

                         r_emit_c = r_emit_c +                                                     &
                                    r_emit_d / width_sglcan(iurb,id,j,i) * fr_udir(iurb,id,j,i)
                      ENDIF
                   ENDDO
                   uprad(j,i) = uprad(j,i) + r_emit_c * fr_uclass(iurb,j,i)

                ENDIF

             ENDDO

          ENDDO
       ENDDO

    ENDIF

 END FUNCTION uprad


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine computes the effects of a shadow induced by a building of height hu, on a portion
!> of wall between z1 and z2. See equation A10, and correction described below formula A11, and
!> figure A1. Basically rd is the ratio between the horizontal surface illuminated and the portion
!> of wall. Referring to figure A1, multiplying radiation flux density on a horizontal surface
!> (rs1d) by x1-x2 we have the radiation energy per unit time. Dividing this by z2-z1, we obtain
!> the radiation flux density reaching the portion of the wall between z2 and z1 (everything is
!> assumed in 2D).
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION shade_wall( z1, z2, hu, tzr, saa, ws )

    IMPLICIT NONE

    REAL(wp) ::  shade_wall      !< Ratio between (x1-x2)/(z2-z1), see Fig. 1A. (Schubert 2013)
                                 !< Multiplying rd by rs1d (radiation flux density on
                                 !< a horizontal surface) gives the radiation flux
                                 !< density on the portion of wall between z1 and z2.

    REAL(wp), INTENT(in) ::  hu  !< Height of the building that generates the shadow
    REAL(wp), INTENT(in) ::  saa !< sin Angle between the sun direction and the face of the wall (A12)
    REAL(wp), INTENT(in) ::  tzr !< Solar zenith angle
    REAL(wp), INTENT(in) ::  ws  !< Width of the street
    REAL(wp), INTENT(in) ::  z1  !< Height of the level z_u(iz)
    REAL(wp), INTENT(in) ::  z2  !< Height of the level z_u(iz+1)
!
!-- Internal variables
    REAL(wp) ::  x1  !< see Fig. A1 in Schubert (2013)
    REAL(wp) ::  x2  !< see Fig. A1 in Schubert (2013)


    x1 = MIN( ( hu - z1 ) * tzr, MAX( 0.0_wp, ws / saa ) )

    x2 = MAX( ( hu - z2 ) * tzr, 0.0_wp )

    shade_wall = MAX( 0.0_wp, saa * MAX( 0.0_wp, x1 - x2 ) / ( z2 -z1 ) )

 END FUNCTION shade_wall


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> inhomogenities of linear system of equations for shortwave
!> radiation: includes diffuse radiation from the sky
!--------------------------------------------------------------------------------------------------!
 PURE FUNCTION short_rad_inhom( i, j, id, iurb, ndimfull, ndim, rsww, rswe, rsr, rsg, rsdiffww,    &
                                rsdiffwe, rsdiffr, rsdiffg, rc )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) ::  i         !<
    INTEGER(iwp), INTENT(in) ::  id        !<
    INTEGER(iwp), INTENT(in) ::  iurb      !<
    INTEGER(iwp), INTENT(in) ::  j         !<
    INTEGER(iwp), INTENT(in) ::  ndim      !<
    INTEGER(iwp), INTENT(in) ::  ndimfull  !<

    REAL(wp) ::  short_rad_inhom(ndimfull)  !<  function type

    REAL(wp), INTENT(in) ::  rc       !<
    REAL(wp), INTENT(in) ::  rsdiffg  !<
    REAL(wp), INTENT(in) ::  rsg      !<

    REAL(wp), DIMENSION(nz(iurb,id,j,i)+1), INTENT(in) ::  rsdiffr   !<
    REAL(wp), DIMENSION(nz(iurb,id,j,i)+1), INTENT(in) ::  rsr       !<
    REAL(wp), DIMENSION(nz(iurb,id,j,i))  , INTENT(in) ::  rsdiffwe  !<
    REAL(wp), DIMENSION(nz(iurb,id,j,i))  , INTENT(in) ::  rsdiffww  !<
    REAL(wp), DIMENSION(nz(iurb,id,j,i))  , INTENT(in) ::  rswe      !<
    REAL(wp), DIMENSION(nz(iurb,id,j,i))  , INTENT(in) ::  rsww      !<

!
!-- Internal variables
    INTEGER(iwp) ::  iu  !<
    INTEGER(iwp) ::  ju  !<
    INTEGER(iwp) ::  ou  !<


!
!-- West wall.
    DO  iu = 1, nz(iurb,id,j,i)
!
!--    Incoming from sky:
       short_rad_inhom(iu) = rsww(iu)

       IF ( lrroofs )  THEN
          DO  ou = 1, nz(iurb,id,j,i)+1
!
!--          Diffuse from sky:
             short_rad_inhom(iu) = short_rad_inhom(iu) +                                           &
                                   fsow(iurb,id,ou,iu,j,i) * rsdiffww(iu) * fr_roof(iurb,id,ou,j,i)
!
!--          Sky rad. from side other canyon:
             DO ju = 1,nz(iurb,id,j,i)
                short_rad_inhom(iu) = short_rad_inhom(iu) +                                        &
                                      fwow(iurb,id,ju,ou,iu,j,i) * rsdiffwe(ju) * rc *             &
                                      ( 1.0_wp - fr_wall(iurb,id,ju,j,i)) * fr_roof(iurb,id,ou,j,i )
             ENDDO
          ENDDO
       ELSE
!
!--       Diffuse from sky:
          short_rad_inhom(iu) = short_rad_inhom(iu) + fsw(iurb,id,iu,j,i) * rsdiffww(iu)
!
!--       Sky rad. from side:
          DO  ju = 1, nz(iurb,id,j,i)
             short_rad_inhom(iu) = short_rad_inhom(iu) + fww(iurb,id,ju,iu,j,i) * rsdiffwe(ju) *   &
                                                         rc * ( 1.0_wp - fr_wall(iurb,id,ju,j,i) )
          ENDDO
       ENDIF

    ENDDO
!
!-- East wall:
    DO  iu = 1+nz(iurb,id,j,i), 2*nz(iurb,id,j,i)
!
!--    From sky
       short_rad_inhom(iu) = rswe(iu-nz(iurb,id,j,i))

       IF ( lrroofs )  THEN
          DO  ou = 1, nz(iurb,id,j,i)+1
!
!--          Diffuse from sky:
             short_rad_inhom(iu) = short_rad_inhom(iu) + fsow(iurb,id,ou,iu-nz(iurb,id,j,i),j,i) * &
                                                         rsdiffwe(iu-nz(iurb,id,j,i)) *            &
                                                         fr_roof(iurb,id,ou,j,i)
!
!--          Sky rad. from side of other canyon:
             DO ju = 1, nz(iurb,id,j,i)
                short_rad_inhom(iu) = short_rad_inhom(iu) +                                        &
                                      fwow(iurb,id,ju,ou,iu-nz(iurb,id,j,i),j,i) * rsdiffwe(ju) *  &
                                      rc * ( 1.0_wp - fr_wall(iurb,id,ju,j,i) ) *                  &
                                      fr_roof(iurb,id,ou,j,i)
             ENDDO
          ENDDO
       ELSE
!
!--       Diffuse from sky:
          short_rad_inhom(iu) = short_rad_inhom(iu) +                                              &
                                fsw(iurb,id,iu-nz(iurb,id,j,i),j,i) * rsdiffwe(iu-nz(iurb,id,j,i))
!
!--       Sky rad. from side:
          DO  ju = 1, nz(iurb,id,j,i)
             short_rad_inhom(iu) = short_rad_inhom(iu) +                                           &
                                   fww(iurb,id,ju,iu-nz(iurb,id,j,i),j,i) * rsdiffwe(ju) * rc *    &
                                   ( 1.0_wp - fr_wall(iurb,id,ju,j,i) )
          ENDDO
       ENDIF

    ENDDO
!
!-- Ground
!-- From sky:
    short_rad_inhom(2*nz(iurb,id,j,i)+1) = rsg

    IF ( lrroofs )  THEN
!
!--    Diffuse from nearer side:
       DO  iu = 1, nz(iurb,id,j,i)
          short_rad_inhom(2*nz(iurb,id,j,i)+1) = short_rad_inhom(2*nz(iurb,id,j,i)+1) +            &
                                                 fwg(iurb,id,iu,j,i) *                             &
                                                 ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) * rsdiffg * rc
       ENDDO

       DO  ou = 1, nz(iurb,id,j,i)+1
!
!--       Diffuse from sky:
          short_rad_inhom(2*nz(iurb,id,j,i)+1) = short_rad_inhom(2*nz(iurb,id,j,i)+1) +            &
                                                 fsog(iurb,id,ou,j,i) * rsdiffg *                  &
                                                 fr_roof(iurb,id,ou,j,i)
!
!--       Diffuse from farer side:
          DO  iu = 1, nz(iurb,id,j,i)
             short_rad_inhom(2*nz(iurb,id,j,i)+1) = short_rad_inhom(2*nz(iurb,id,j,i)+1) +         &
                                                    fwog(iurb,id,iu,ou,j,i) * rc *                 &
                                                    ( 1.0_wp - fr_wall(iurb,id,iu,j,i) ) * rsdiffg &
                                                    * fr_roof(iurb,id,ou,j,i)
          ENDDO
       ENDDO
    ELSE
!
!--    Diffuse from sky:
       short_rad_inhom(2*nz(iurb,id,j,i)+1) = short_rad_inhom(2*nz(iurb,id,j,i)+1) +               &
                                              fsg(iurb,id,j,i) * rsdiffg
!
!--    Sky rad. from side:
       DO iu = 1, nz(iurb,id,j,i)
          short_rad_inhom(2*nz(iurb,id,j,i)+1) = short_rad_inhom(2*nz(iurb,id,j,i)+1) +            &
                                                 2.0_wp * fwg(iurb,id,iu,j,i) *                    &
                                                 ( 1.0_wp - fr_wall(iurb,id,iu,j,i)) * rsdiffg * rc
       ENDDO
    ENDIF
!
!-- Roofs
    IF ( lrroofs )  THEN
       DO  iu = 2*nz(iurb,id,j,i)+2, ndim
!
!--       Direct and diffuse from sky:
          short_rad_inhom(iu) = rsr(iu-(2*nz(iurb,id,j,i)+1)) +                                    &
                                fsr(iurb,id,iu-(2*nz(iurb,id,j,i)+1),j,i) *                        &
                                rsdiffr(iu-(2*nz(iurb,id,j,i)+1))
!
!--       Twice from side:
          DO ju = iu-2*nz(iurb,id,j,i)-1, nz(iurb,id,j,i)
             short_rad_inhom(iu) = short_rad_inhom(iu) +                                           &
                                   2.0_wp * fwr(iurb,id,ju,iu-(2*nz(iurb,id,j,i)+1),j,i) *         &
                                   ( 1.0_wp - fr_wall(iurb,id,ju,j,i) ) *                          &
                                   rsdiffr(iu-(2*nz(iurb,id,j,i)+1)) * rc
          ENDDO
       ENDDO
    ENDIF

  END FUNCTION short_rad_inhom


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> interpolate PALM array to the urban grid
!--------------------------------------------------------------------------------------------------!
  PURE FUNCTION interpol( nz, nz_u, z, z_u, c )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) ::  nz       !< number of mesoscale height levels
    INTEGER(iwp), INTENT(in) ::  nz_u     !< number of urban height levels

    REAL(wp) , INTENT(in) ::  c(1:nz)     !< field which has to be interpolated
    REAL(wp) , INTENT(in) ::  z(1:nz+1)   !< altitude of the mesoscale cell interface
    REAL(wp) , INTENT(in) ::  z_u(nz_u+2) !< altitude of the urban cell interface

    REAL(wp) ::  interpol(nz_u+1)  !< function type
!
!-- Internal variables
    INTEGER(iwp) ::  iz    !< help variables
    INTEGER(iwp) ::  iz_u  !< help variables

    REAL(wp) ::  c_tot  !< help variables
    REAL(wp) ::  dz     !< help variables


    DO  iz_u = 1, nz_u+1
       c_tot = 0.0_wp
       DO  iz = 1, nz
          dz = MAX( MIN( z(iz+1), z_u(iz_u+1) ) - MAX( z(iz), z_u(iz_u) ), 0.0_wp )
!
!--       Data starts is from top to down:
          c_tot = c_tot + c(nz-iz+1) * dz
       ENDDO
       interpol(iz_u) = c_tot / ( z_u(iz_u+1) - z_u(iz_u) )
    ENDDO

 END FUNCTION interpol

 END MODULE dcep_mod
