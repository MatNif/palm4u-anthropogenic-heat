! @file chem_emis_pollen_mod.f90
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
! Copyright 2018-2021 Leibniz Universitaet Hannover
! Copyright 2018-2021 Karlsruhe Institute of Technology
! Copyright 2018-2021 Freie Universitaet Berlin
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Parameterized pollen emission based on the EMPOL model (Zink et al, GMD 2013).
!--------------------------------------------------------------------------------------------------!
 MODULE chem_emis_pollen_mod

    USE chem_emis_generic_mod

    USE chem_emis_pollen_data_mod

    USE chem_emis_vsrc_mod

    USE chem_modules,                                                                              &
        ONLY:  nc_field_length

    USE control_parameters,                                                                        &
        ONLY:  message_string

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nys,                                                                                &
               nyn,                                                                                &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_top_ind

    USE kinds

    IMPLICIT NONE

    SAVE 
    PRIVATE
!
!-- Constant values used throughout the module.
    CHARACTER(LEN=*), PARAMETER ::  str_empty           = 'novalue'   !< fill value for empty string
    CHARACTER(LEN=*), PARAMETER ::  str_pollen_betula   = 'POL_BETU'  !< species string for betula (tree)
    CHARACTER(LEN=*), PARAMETER ::  str_pollen_alder    = 'POL_ALNU'  !< species string for adler (tree)
    CHARACTER(LEN=*), PARAMETER ::  str_pollen_ambrosia = 'POL_AMBR'  !< species string for ambrosia (vegetation)
    CHARACTER(LEN=*), PARAMETER ::  str_pollen_poaceae  = 'POL_POAC'  !< species string for poaceae (vegetation)

    INTEGER(iwp), PARAMETER ::  int_empty        = -127  !< fill value for integer variables

    REAL(wp), PARAMETER ::  count_scale      = 1.0_wp  !< 1.0E-6_wp   magnitude scale for pollen counts
    REAL(wp), PARAMETER ::  default_f_season = 0.9_wp  !< seasonal yield factor
    REAL(wp), PARAMETER ::  default_tuning   = 1.0_wp  !< tuning factor
!
!-- Derived type for mapping user to data module species index.
    TYPE ::  emis_pollen_species_attributes
       INTEGER(iwp) ::  crown      = 0               !< tree crown cell count
       INTEGER(iwp) ::  data_index = int_empty       !< index in species attributes and parameters
       INTEGER(iwp) ::  footprint  = 0               !< tree footprint cell count
       INTEGER(iwp) ::  patch      = 0               !< vegetation patch cell count

       REAL(wp)     ::  f_season   = 1.0_wp          !< yield factor due to seasonal conditions 
       REAL(wp)     ::  max_daily  = 0.0_wp          !< maximum daily available pollen
       REAL(wp)     ::  max_rate   = 0.0_wp          !< maximum pollen release rate
       REAL(wp)     ::  scale      = 1.0_wp          !< flux to source scaling factor
       REAL(wp)     ::  tuning     = default_tuning  !< tuning factor
    END TYPE emis_pollen_species_attributes
!
!-- Derived type for storing volume source status.
    TYPE :: emis_pollen_vsrc_attributes
       LOGICAL ::  is_alder    = .FALSE.  !< if volume contains ALNU species
       LOGICAL ::  is_ambrosia = .FALSE.  !< if volume contains AMBR species
       LOGICAL ::  is_betula   = .FALSE.  !< if volume contains BETU species
       LOGICAL ::  is_poaceae  = .FALSE.  !< if volume contains POAC species

       INTEGER(iwp) ::  species_index = int_empty  !< emission species index

       REAL(wp) ::  available_total = 0.0_wp  !< total pollen deposited on plant
       REAL(wp) ::  pollen_pool     = 0.0_wp  !< daily pollen pool
       REAL(wp) ::  release_rate    = 0.0_wp  !< nominal pollen release rate
       REAL(wp) ::  yield_release   = 1.0_wp  !< amalgamated release yield factor
       REAL(wp) ::  yield_emis      = 1.0_wp  !< amalgamated emission yield factor
    END TYPE emis_pollen_vsrc_attributes
!
!-- Interface to KPP chemical mechanism.
    INTEGER(iwp) ::  num_emis_species  !< num of emission species

    TYPE(chem_emis_species),              ALLOCATABLE, DIMENSION(:) ::  emis_species       !< user / KPP index mapping
    TYPE(emis_pollen_species_attributes), ALLOCATABLE, DIMENSION(:) ::  emis_species_info  !< species attributes
!
!-- Interface to chem_emis_vsrc_mod.
    INTEGER(iwp) ::  num_vsrc_pos  !< # volume sources
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::  vsrc_emis_value  !< volume source volues
    
    TYPE(emis_pollen_vsrc_attributes), ALLOCATABLE, DIMENSION(:) ::  vsrc_info  !< volume source attributes
    TYPE(chem_emis_vsrc_pos), ALLOCATABLE, DIMENSION(:)          ::  vsrc_pos   !< volume source positions
!
!-- Default base settings.
    INTEGER(iwp), PARAMETER ::  pollen_model_zink_default      = 0  !< default pollen model ID (zink)
    INTEGER(iwp), PARAMETER ::  pollen_pool_reset_hour_default = 0  !< default hour (UTC) 4 pollen pool reset
    INTEGER(iwp), PARAMETER ::  tke_scheme_default             = 0  !< default TKE scheme (default)
    INTEGER(iwp), PARAMETER ::  tke_scheme_adhoc               = 1  !< Holst TKE scheme
    INTEGER(iwp), PARAMETER ::  tke_scheme_dynamic             = 2  !< Seuhring TKE scheme

    REAL(wp), PARAMETER ::  model_update_interval_default = 300.0_wp !< default update interval
    REAL(wp), PARAMETER ::  tke_sgs_fraction_default      = 0.1_wp  !< default SGS TKE fraction
!
!-- Base settings.
    INTEGER(iwp) ::  pollen_model            !< pollen model
    INTEGER(iwp) ::  pollen_pool_reset_hour  !< UTC hour for pool pool reset
    INTEGER(iwp) ::  tke_scheme              !< TKE scheme

    LOGICAL ::  ignore_precipitation   !< if true all calculations are dry
    LOGICAL ::  ignore_solar_activity  !< if true pollen releases all day

    REAL(wp) ::  model_update_interval  !< default update interval
    REAL(wp) ::  tke_sgs_fraction       !< default SGS TKE fraction
!
!-- Triggers.
    LOGICAL ::  pollen_pool_reset_trigger_armed  !< whether trigger for pollen pool reset is armed

    REAL(wp) ::  time_since_last_update  !< time since last model update [s]
!
!-- Interface to mode-specific methods.
    INTERFACE chem_emis_pollen_init
       MODULE PROCEDURE chem_emis_pollen_init
    END INTERFACE chem_emis_pollen_init

    INTERFACE chem_emis_pollen_update
       MODULE PROCEDURE chem_emis_pollen_update
    END INTERFACE chem_emis_pollen_update

    INTERFACE chem_emis_pollen_update_settling
       MODULE PROCEDURE chem_emis_pollen_update_settling
    END INTERFACE chem_emis_pollen_update_settling
!
!-- Public methods specific to generic emission mode (i.e., this module).
    PUBLIC ::  chem_emis_pollen_init,                                                              &
               chem_emis_pollen_update,                                                            &
               chem_emis_pollen_update_settling

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Wrapper for intitializing emission mode.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_pollen_init( )

    IMPLICIT NONE


    CALL location_message( 'initializing pollen emission module', 'start' )
    CALL init_base_settings( )
    CALL init_species( ) 
    CALL init_volume_sources( )
    CALL location_message( 'initializing pollen emission module', 'finished'  )

 END SUBROUTINE chem_emis_pollen_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates emissions at every time step (called in time_integration).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_pollen_update( )

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_update_trigger

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_update_source

    USE radiation_model_mod,                                                                       &
        ONLY:  cos_zenith

    IMPLICIT NONE

    INTEGER(iwp) ::  k  !< generic index


!
!-- Trigger pollen reset at user-defined hour (UTC).
    IF ( pollen_pool_reset_trigger() )  THEN
       CALL reset_pollen_pool( )
    ENDIF
!
!-- Update pollen state and emissions only when it is time to release pollen.
    IF ( ignore_solar_activity  .OR.  ( cos_zenith >= 0.0_wp ) )  THEN
       IF ( chem_emis_update_trigger( time_since_last_update, model_update_interval ) )  THEN
          CALL update_pollen_source_attributes( )
       ENDIF
       CALL update_pollen_emissions( )
    ENDIF
!
!-- Transfer emission values to prognostic equation source terms.
    DO  k = 1, num_emis_species
       CALL chem_emis_vsrc_update_source( emis_species(k)%mech_index, vsrc_pos,                    &
                                          vsrc_emis_value(k,:) )
    ENDDO

 END SUBROUTINE chem_emis_pollen_update


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Corrects pollen concentration based due to vertical settling.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_pollen_update_settling( i, j )

    USE arrays_3d,                                                                                 &
        ONLY:  ddzu

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  g

    USE chem_modules,                                                                              &
        ONLY:  chem_species

    USE control_parameters,                                                                        &
        ONLY:  dt_3d

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  i  !< cell i index
    INTEGER(iwp), INTENT(IN) ::  j  !< cell j index

    INTEGER(iwp) ::  k_data  !< emission species index in pollen data module
    INTEGER(iwp) ::  k_emis  !< emission species index
    INTEGER(iwp) ::  k_mech  !< KPP species index

    REAL(wp) ::  mean_diameter  !< mean pollen diameter
    REAL(wp) ::  pollen_factor  !< scalar factor for settling velocity

    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  column  !< vertical column of concentration
    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  mu_air  !< dynamic viscosity of air
    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  ws      !< pollen settling velocity
    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  z_flux  !< downward momentum flux


    ALLOCATE( column(nzb:nzt+1) )
    ALLOCATE( mu_air(nzb:nzt+1) )
    ALLOCATE( ws(nzb:nzt+1) )
    ALLOCATE( z_flux(nzb:nzt+1) )
!
!-- Calculate viscosity of air.
    mu_air = get_dynamic_viscosity_sutherland( )

    DO  k_emis = 1, num_emis_species
       k_mech = emis_species(k_emis)%mech_index
       k_data = emis_species_info(k_emis)%data_index
!
!--    Make a working copy of veretical column.
       column = chem_species(k_mech)%conc(:,j,i)
!
!--    Calculate settling velocity.
!--    Note: the density term should be (rho_pollen - rho_air) instead of (rho_pollen)
!--          but since rho_pollen >> rho_air it can be dropped to save one
!--          array calculation without losing accuracy.
       mean_diameter = 1.0E-6_wp * pol_diameter_dry(k_data) / pol_diameter_ratio(k_data)
       pollen_factor = ( g / 18.0_wp ) * pol_density(k_data) * mean_diameter**2
       ws = pollen_factor / mu_air 
!
!--    Calculate vertical pollen flux over the column (cap maximum to amount in cell),
!--    displace it downwards (cap minimum to zero). Displacement on edge node (nzb) is ignored.
       z_flux(1:nzt+1) = MIN( 1.0_wp, dt_3d*ddzu(1:nzt+1)*ws(1:nzt+1) ) * column(1:nzt+1)
       column(nzb+1:nzt+1) = column(nzb+1:nzt+1) - z_flux(nzb+1:nzt+1)
       column(nzb+1:nzt) = MAX( 0.0_wp, column(nzb+1:nzt) + z_flux(nzb+2:nzt+1) )
!
!--    Reassign vertical column to prognostic equation.
       chem_species(k_mech)%conc(:,j,i) = column
    ENDDO

    DEALLOCATE( column )
    DEALLOCATE( mu_air )
    DEALLOCATE( ws )
    DEALLOCATE( z_flux )

 END SUBROUTINE chem_emis_pollen_update_settling


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> read namelist items and check validity
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_base_settings( )

    USE chem_modules,                                                                              &
        ONLY:  epol_ignore_precip,                                                                 &
               epol_ignore_solar,                                                                  &
               epol_model,                                                                         &
               epol_tke_scheme,                                                                    &
               epol_tke_sgs_fraction,                                                              &
               epol_pool_reset_hour,                                                               &
               epol_update_interval

    IMPLICIT NONE


!
!-- Module flags (ignore_precipitation is always foreced to true since there is currently no
!-- prognostic variable for it).
    ignore_precipitation = epol_ignore_precip
    IF ( .NOT. ignore_precipitation )  ignore_precipitation = .TRUE.

    ignore_solar_activity = epol_ignore_solar
!
!-- Pollen model (only Zink is available at the moment).
    SELECT CASE ( TRIM(epol_model) )

       CASE ( 'zink' )
          pollen_model = pollen_model_zink_default

       CASE DEFAULT
          message_string = 'Invalid pollen emission model selected. ' //                           &
                           'Reverting to the default Zink pollen model.'
          CALL message( 'init_base_setting', 'CHM0034', 0, 1, -1, 6, 0 )
          pollen_model = pollen_model_zink_default

    END SELECT
!
!-- TKE scheme.
    SELECT CASE ( TRIM(epol_tke_scheme) )

       CASE ( 'default' )
          tke_scheme = tke_scheme_default

       CASE ( 'adhoc'   )
          tke_scheme = tke_scheme_adhoc

       CASE ( 'dynamic' ) 
          tke_scheme = tke_scheme_dynamic 

       CASE DEFAULT
          message_string = 'unknown epol_tke_scheme = "' // TRIM( epol_tke_scheme ) // '"'
          CALL message( 'init_base_setting', 'CHM0035', 2, 2, 0, 6, 0 )

    END SELECT

    tke_sgs_fraction = epol_tke_sgs_fraction
    IF ( tke_sgs_fraction <= 0.0_wp )  tke_sgs_fraction = tke_sgs_fraction_default
!
!-- Update and reset timings.
    model_update_interval = epol_update_interval
    IF ( model_update_interval < 0.0_wp )  model_update_interval = model_update_interval_default

    pollen_pool_reset_hour = epol_pool_reset_hour
    IF ( pollen_pool_reset_hour < 0 )  pollen_pool_reset_hour = pollen_pool_reset_hour_default
!
!-- Initialize triggers.
    time_since_last_update          = 0.0_wp
    pollen_pool_reset_trigger_armed = .FALSE.
!
!-- Initialize array sizes.
    num_emis_species = 0
    num_vsrc_pos     = 0

 END SUBROUTINE init_base_settings


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initiate pollen species:
!>  - tag user species to KPP species
!>  - assign species attributes and parameters
!>  - pollen release attributes will be assigned in init_volume_sources
!>    when flux-to-source scaling are calculated
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_species( )

    USE chem_modules,                                                                              &
        ONLY:  epol_seasonal_factors,                                                              &
               epol_specs_names,                                                                   &
               epol_tuning_factors

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_species_match

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< generic indices
    INTEGER(iwp) ::  k  !< generic indices


!
!-- Tag user species to KPP mechanism.
    CALL chem_emis_species_match( num_emis_species, emis_species, epol_specs_names )

    IF ( num_emis_species == 0 )  THEN
       message_string = 'no valid pollen species defined'
       CALL message( 'init_species', 'CHM0036', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Assign species attributes and parameters after memory allocation,
!-- effective referencing to corresponding entries in chem_emis_pollen_data_mod.
    ALLOCATE( emis_species_info(num_emis_species) )


    DO  k = 1, num_emis_species
!
!--    Read and validate seasonal factors from namelist.
       emis_species_info(k)%f_season = epol_seasonal_factors(emis_species(k)%user_index)
       IF ( ( emis_species_info(k)%f_season < 0.0_wp )  .OR.                                       &
            ( emis_species_info(k)%f_season > 1.0_wp ) )                                           &
       THEN
          WRITE( message_string, * ) 'seasonal factor = ', emis_species_info(k)%f_season,          &
                                     ' reset to the default value 0.9'
          CALL message( 'init_species', 'CHM0037', 0, 1, 0, 6, 0 )

          emis_species_info(k)%f_season = default_f_season
       ENDIF
!
!--    Read and validate tuning factors from namelist.
       emis_species_info(k)%tuning = epol_tuning_factors(emis_species(k)%user_index)
       IF ( emis_species_info(k)%tuning < 0.0_wp )  THEN
           emis_species_info(k)%tuning = default_tuning
       ENDIF
!
!--    Create reference to data module (use default if species not found).
       emis_species_info(k)%data_index = pollen_default_offset

       DO  i = pollen_default_offset, pollen_num_total
          IF ( TRIM(emis_species(k)%name_str) /= TRIM(epol_species_names(i)) )  CYCLE
          emis_species_info(k)%data_index = i
          EXIT
       ENDDO

    ENDDO

 END SUBROUTINE init_species


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initiates volume sources:
!>  - reads static file (NC)
!>  - identifies volume source type (tree/vegetation)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_volume_sources( )

    USE arrays_3d,                                                                                 &
        ONLY:  dzu
    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_init_volume_source

    USE chem_modules,                                                                              &
        ONLY:  epol_tree_specs,                                                                    &
               epol_vegetation_specs 

    USE control_parameters,                                                                        &
        ONLY:  coupling_char

    USE indices,                                                                                   &
        ONLY:  nzb

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  close_input_file,                                                                   &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               input_file_static,                                                                  &
               open_read_file

    IMPLICIT NONE
!.
!-- NetCDF attribute parameters (dimension and variable names)
    CHARACTER(LEN=*), PARAMETER ::  nc_dim_nvegetation_pars = 'nvegetation_pars'
    CHARACTER(LEN=*), PARAMETER ::  nc_dim_nzlad            = 'zlad'
    CHARACTER(LEN=*), PARAMETER ::  nc_var_tree_type        = 'tree_type'
    CHARACTER(LEN=*), PARAMETER ::  nc_var_vegetation_type  = 'vegetation_type'
!
!-- NetCDF attributes.
    INTEGER(iwp) ::  ncid  !< netCDF handle for static file

    INTEGER(iwp) ::  nc_nvegetation_pars
    INTEGER(iwp) ::  nc_nzlad

    INTEGER(KIND=1), ALLOCATABLE, DIMENSION(:,:,:) ::  nc_tree_type        !< tree_type
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:)      ::  nc_vegetation_type  !< vegetation_type
!
!-- Local variables.
    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j
    INTEGER(iwp) ::  k          !< generic indices
    INTEGER(iwp) ::  k_data     !< emission species index in data module
    INTEGER(iwp) ::  k_emis     !< emission species index
    INTEGER(iwp) ::  k_species  !< tree species index
    INTEGER(iwp) ::  k_vsrc     !< volume source index
    INTEGER(iwp) ::  nz_tree    !< nc_nzlad - 1

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  tree_species        !< map 4 epol_tree_specs
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  vegetation_species  !< map 4 epol_vegetation_specs

!
!-- Get out if no species is tagged; there is nothing to do.
    IF ( num_emis_species == 0 )  RETURN
!
!-- Map user pollen species.
    ALLOCATE( tree_species(num_emis_species) )
    ALLOCATE( vegetation_species(num_emis_species) )

    DO  k = 1,num_emis_species
       tree_species(k)       = epol_tree_specs(emis_species(k)%user_index)
       vegetation_species(k) = epol_vegetation_specs(emis_species(k)%user_index)
    ENDDO
!
!-- Open static file and grab dimensions and variable data.
    CALL open_read_file( TRIM(input_file_static) // TRIM(coupling_char), ncid )

    CALL get_dimension_length( ncid, nc_nvegetation_pars, nc_dim_nvegetation_pars )
    CALL get_dimension_length( ncid, nc_nzlad, nc_dim_nzlad )

    nz_tree = nc_nzlad - 1
    ALLOCATE( nc_tree_type(0:nz_tree,nys:nyn,nxl:nxr) )
    ALLOCATE( nc_vegetation_type(nys:nyn,nxl:nxr) )

    CALL get_variable( ncid, nc_var_tree_type, nc_tree_type, nxl, nxr, nys, nyn, 0, nz_tree )
    CALL get_variable( ncid, nc_var_vegetation_type, nc_vegetation_type, nxl, nxr, nys, nyn )

    CALL close_input_file( ncid )
!
!-- Generate volume sources.
!-- First pass : tally up volume sources, tree crown and footprint size.
    k_vsrc = 0

    DO  j = nys, nyn
       DO  i = nxl, nxr
          IF ( ANY( nc_vegetation_type(j,i) == vegetation_species ) )  THEN
             k_vsrc = k_vsrc + 1
             k_species = locate_index( vegetation_species, nc_vegetation_type(j,i) )
             emis_species_info(k_species)%patch = emis_species_info(k_species)%patch + 1
          ENDIF

          IF ( ANY( nc_tree_type(0,j,i) == tree_species ) )  THEN
             k_species = locate_index( tree_species, INT( nc_tree_type(0,j,i), iwp ) )
             emis_species_info(k_species)%footprint = emis_species_info(k_species)%footprint + 1
          ENDIF
!
!--       Starting with 1 as z=0 is the tree footprint layer.
          DO  k = 1, nz_tree
             IF ( ANY( nc_tree_type(k,j,i) == tree_species ) )  THEN
                k_vsrc = k_vsrc + 1
                k_species = locate_index( tree_species, INT( nc_tree_type(k,j,i), iwp ) )
                emis_species_info(k_species)%crown = emis_species_info(k_species)%crown + 1
             ENDIF
          ENDDO 
       ENDDO
    ENDDO
!
!-- Calculate scale for unit conversion (aka "vertical distribution")
!-- pollen emissions produced by the model are expressed in units of m3 / m2
!-- which must be converted into units of                            m3 / m3
!-- (dx dy dz) / m2 x (cell horizontal area / cell volume)
!-- for tree       : cell horizontal area = footprint cell count x (dx dy)
!--                  cell volume          = crown cell count x (dx dy dz)
!-- for vegetation : cell horizontal area = (dx dy)
!--                  cell volume          = (dx dy dz)
!--                  (since crown = footprint )
!-- expression for scale
!-- for vegetation :  (dx dy)
!-- for tree       :  (dx dy) x ( footprint cell count / crown cell count ).
!-- Calculate emission species attributes:
    DO  k = 1,num_emis_species
       emis_species_info(k)%scale = 0.0_wp

       IF ( ( emis_species_info(k)%footprint > 0 )  .OR.                                           &
            ( emis_species_info(k)%patch > 0 ) )                                                   &
           emis_species_info(k)%scale = 1.0_wp / dzu(nzb+1)

       IF ( emis_species_info(k)%crown > 0 )                                                       &
           emis_species_info(k)%scale = emis_species_info(k)%scale *                               &
                                        REAL( emis_species_info(k)%footprint, wp ) /               &
                                        REAL( emis_species_info(k)%crown, wp ) 
!
!--    Calculate pollen properties (maximum daily availability and release rate).
       k_data = emis_species_info(k)%data_index
       emis_species_info(k)%max_daily = count_scale *                                              &
                                        emis_species_info(k)%scale  *                              &
                                        emis_species_info(k)%tuning *                              &
                                        emis_species_info(k)%f_season *                            &
                                        pol_qday(k_data)
       emis_species_info(k)%max_rate  = emis_species_info(k)%max_daily /                           &
                                        pol_qt_duration(k_data)
    ENDDO
!
!-- Update module volume source count.
    num_vsrc_pos = k_vsrc
!
!-- Second pass : assign volume source positions
!--               (sans topological correction)
    ALLOCATE( vsrc_emis_value(num_emis_species,num_vsrc_pos) )
    vsrc_emis_value = 0.0_wp

    ALLOCATE( vsrc_pos(num_vsrc_pos)  )
    ALLOCATE( vsrc_info(num_vsrc_pos) )

    k_vsrc = 0

    DO  j = nys, nyn
       DO  i = nxl, nxr

          IF ( ANY( nc_vegetation_type(j,i) == vegetation_species ) )  THEN
             k_vsrc = k_vsrc + 1
             vsrc_pos(k_vsrc)%i = i 
             vsrc_pos(k_vsrc)%j = j
             vsrc_pos(k_vsrc)%k = 0

             k_emis = MINLOC( ABS( vegetation_species - nc_vegetation_type(j,i) ), 1 )
             
             SELECT CASE ( TRIM(emis_species(k_emis)%name_str) )

                CASE  ( TRIM( str_pollen_ambrosia ) )
                   vsrc_info(k_vsrc)%is_ambrosia = .TRUE.

                CASE  ( TRIM( str_pollen_poaceae ) )
                   vsrc_info(k_vsrc)%is_poaceae  = .TRUE.

             END SELECT
          ENDIF
!
!--       Starting with 1 again (0 is footprint layer).
          DO  k = 1, nz_tree

             IF ( ANY( nc_tree_type(k,j,i) == tree_species ) )  THEN
                k_vsrc = k_vsrc + 1
                vsrc_pos(k_vsrc)%i = i
                vsrc_pos(k_vsrc)%j = j
                vsrc_pos(k_vsrc)%k = k

                k_emis = MINLOC( ABS( tree_species - nc_tree_type(k,j,i) ), 1 )

                SELECT CASE ( TRIM(emis_species(k_emis)%name_str) )

                   CASE ( TRIM(str_pollen_alder) )
                      vsrc_info(k_vsrc)%is_alder = .TRUE.

                   CASE ( TRIM(str_pollen_betula) )
                      vsrc_info(k_vsrc)%is_betula = .TRUE.

                END SELECT

             ENDIF
          ENDDO
       ENDDO
    ENDDO 
!
!-- Assign volume source attributes and perform topological correction.
    DO  k_vsrc = 1, num_vsrc_pos
       i = vsrc_pos(k_vsrc)%i
       j = vsrc_pos(k_vsrc)%j
       k = vsrc_pos(k_vsrc)%k
!
!--    Vegetation (tree is always at k>=1).
       IF ( k == 0 )  THEN
          vsrc_info(k_vsrc)%species_index =                                                        &
                   locate_index( vegetation_species, nc_vegetation_type(j,i) )
!
!--       First interior cell above terrain.
          vsrc_pos(k_vsrc)%k = 1
       ELSE
          vsrc_info(k_vsrc)%species_index =                                                        &
                   locate_index( tree_species, INT( nc_tree_type(k,j,i), iwp ) )
       ENDIF
!
!--    Terrain correction.
       vsrc_pos(k_vsrc)%k = vsrc_pos(k_vsrc)%k + topo_top_ind(j,i,0)
    ENDDO
!
!-- Register volume source to global volume source listing.
    CALL chem_emis_init_volume_source( vsrc_pos )
!
!-- Update volume source attributes.
       CALL reset_pollen_pool( )
       CALL update_pollen_source_attributes( )

    DEALLOCATE( tree_species )
    DEALLOCATE( vegetation_species )
    DEALLOCATE( nc_tree_type )
    DEALLOCATE( nc_vegetation_type )

 END SUBROUTINE init_volume_sources


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resets pollen pool.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE reset_pollen_pool( )

    IMPLICIT NONE

    INTEGER(iwp) :: k_data  !< species index in pollen data module
    INTEGER(iwp) :: k_emis  !< species index in emission species list
    INTEGER(iwp) :: k_vsrc  !< volume source index


    IF ( num_vsrc_pos == 0 )  RETURN

    DO  k_vsrc = 1, num_vsrc_pos
       CALL get_vsrc_indices( k_emis, k_data, k_vsrc )
       vsrc_info(k_vsrc)%pollen_pool = emis_species_info(k_emis)%max_daily
    ENDDO

 END SUBROUTINE reset_pollen_pool


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Refresh pollen attributes at volume source at update interval.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE update_pollen_source_attributes( )

    USE arrays_3d,                                                                                 &
        ONLY:  e,                                                                                  &
               exner,                                                                              &
               hyp,                                                                                &
               pt,                                                                                 &
               q,                                                                                  &
               rho_air

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  get_relative_humidity

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< cell position indices  x direction
    INTEGER(iwp) ::  j       !< cell position indices  y direction
    INTEGER(iwp) ::  k       !< cell position indices  z direction
    INTEGER(iwp) ::  k_data  !< species index in pollen data module
    INTEGER(iwp) ::  k_emis  !< species index in emission species list
    INTEGER(iwp) ::  k_vsrc  !< volume source index

    REAL(wp) ::  e_phi_rh           !< emission yield factor for relative humidity
    REAL(wp) ::  e_phi_tke          !< emission yield factor for tke
    REAL(wp) ::  r_phi_p            !< amalgamated release yield factor for precipitation
    REAL(wp) ::  ra_phi_meteo       !< amalgamated release yield factor for meteorology
    REAL(wp) ::  ra_phi_plant       !< amalgamated release yield factor for plant physiology
    REAL(wp) ::  r_phi_altitude     !< release yield factor for plant altitude
    REAL(wp) ::  r_phi_coverage     !< release yield factor for plant coverage
    REAL(wp) ::  r_phi_ta           !< release yield factor for temperature
    REAL(wp) ::  r_phi_rh           !< release yield factor for realtive humidity
    REAL(wp) ::  relative_humidity  !< relative humidity
    REAL(wp) ::  temperature        !< thermodynamic temperature
    REAL(wp) ::  tke_total          !< turbulent kinetic energy (resolved scale)


    IF ( num_vsrc_pos == 0 )  RETURN

    DO  k_vsrc = 1, num_vsrc_pos
       i = vsrc_pos(k_vsrc)%i
       j = vsrc_pos(k_vsrc)%j
       k = vsrc_pos(k_vsrc)%k
       CALL get_vsrc_indices( k_emis, k_data, k_vsrc ) 
!
!--    Get thermodynamic properties.
       temperature       = pt(k,j,i) * exner(k)
       relative_humidity = get_relative_humidity( q(k,j,i), temperature, hyp(k), rho_air(k) )
!
!--    Get resolved TKE.
       SELECT CASE ( tke_scheme )

          CASE ( tke_scheme_default )
             tke_total = get_tke_total_default( e(k,j,i), tke_sgs_fraction )

          CASE ( tke_scheme_adhoc )
             tke_total = get_tke_total_adhoc( k_vsrc, tke_sgs_fraction )

          CASE ( tke_scheme_dynamic )
             tke_total = get_tke_total_dynamic( e(k,j,i),k,j,i )

          CASE DEFAULT
             tke_total = get_tke_total_default( e(k,j,i), tke_sgs_fraction )

       END SELECT
!
!--    Update yield factors.
       e_phi_rh  = get_emission_yield_rh( relative_humidity )
       e_phi_tke = get_emission_yield_tke( tke_total )

       r_phi_p        = get_release_yield_precip( )
       r_phi_altitude = get_release_yield_altitude( )
       r_phi_coverage = get_release_yield_coverage( )

       IF ( vsrc_info(k_vsrc)%is_alder )  THEN
          r_phi_rh = get_release_yield_rh_alder( relative_humidity )
          r_phi_ta = get_release_yield_temperature_alder( temperature )

       ELSEIF ( vsrc_info(k_vsrc)%is_poaceae)  THEN
          r_phi_rh = get_release_yield_rh_poaceae( relative_humidity )
          r_phi_ta = get_release_yield_temperature_poaceae( temperature )


       ELSEIF ( vsrc_info(k_vsrc)%is_ambrosia)  THEN   
          r_phi_rh = get_release_yield_rh_ambrosia( relative_humidity )
          r_phi_ta = get_release_yield_temperature_ambrosia( temperature )

 
       ELSEIF ( vsrc_info(k_vsrc)%is_betula )  THEN
          r_phi_rh = get_release_yield_rh_betula( relative_humidity )
          r_phi_ta = get_release_yield_temperature_betula( temperature )

       ENDIF

       ra_phi_plant = r_phi_coverage * r_phi_altitude
       ra_phi_meteo = r_phi_ta * r_phi_rh

       vsrc_info(k_vsrc)%yield_emis    = r_phi_p * e_phi_tke * e_phi_rh
       vsrc_info(k_vsrc)%yield_release = ra_phi_plant
!
!--    If volume source is a tree.
       IF ( emis_species_info(k_emis)%footprint > 0 )  THEN
           vsrc_info(k_vsrc)%yield_release = vsrc_info(k_vsrc)%yield_release * ra_phi_meteo
       ENDIF
!
!--    Update pollen release rate.
       vsrc_info(k_vsrc)%release_rate = vsrc_info(k_vsrc)%yield_release *                          &
                                        emis_species_info(k_emis)%max_rate
    ENDDO

 END SUBROUTINE update_pollen_source_attributes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates pollen volume sources.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE update_pollen_emissions( )

    USE control_parameters,                                                                        &
        ONLY:  dt_3d

    IMPLICIT NONE

    INTEGER(iwp) ::  k_vsrc  !< volume source index
    INTEGER(iwp) ::  k_data  !< attribute index  
    INTEGER(iwp) ::  k_emis  !< species index in KPP mechanism

    REAL(wp) ::  d_pollen_total  !< total pollen not emitted
    REAL(wp) ::  e_pollen_total  !< total pollen emitted
    REAL(wp) ::  r_phi_random    !< release yield factor from random processes
    REAL(wp) ::  r_pollen_old    !< old pollen released
    REAL(wp) ::  r_pollen_new    !< new pollen released
    REAL(wp) ::  r_pollen_total  !< total pollen release


    IF ( num_vsrc_pos == 0 )  RETURN  
!
!-- Reset emission values for the time step.
    vsrc_emis_value = 0.0_wp

    DO  k_vsrc = 1, num_vsrc_pos

       CALL get_vsrc_indices( k_emis, k_data, k_vsrc )
!
!--    Get new pollen state.
       r_phi_random   = get_release_yield_random( pol_rand_coeff(k_data) )
       r_pollen_new   = MIN( vsrc_info(k_vsrc)%pollen_pool, dt_3d*vsrc_info(k_vsrc)%release_rate )
       r_pollen_old   = r_phi_random * vsrc_info(k_vsrc)%available_total
       r_pollen_total = r_pollen_new + r_pollen_old
       e_pollen_total = vsrc_info(k_vsrc)%yield_emis * r_pollen_total
       d_pollen_total = r_pollen_total - e_pollen_total
!
!--    Update volume source pollen attributes.
       vsrc_info(k_vsrc)%pollen_pool = MAX( 0.0_wp, vsrc_info(k_vsrc)%pollen_pool - r_pollen_new )
       vsrc_info(k_vsrc)%available_total = d_pollen_total
!
!--    Update volume source emission as volumetric rate.
       vsrc_emis_value(k_emis,k_vsrc) = e_pollen_total / dt_3d
    ENDDO

 END SUBROUTINE update_pollen_emissions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resets pollen pool.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE get_vsrc_indices( emis_index, data_index, vsrc_index )

    IMPLICIT NONE

    INTEGER(iwp)             ::  data_index  !< species index in pollen data module
    INTEGER(iwp)             ::  emis_index  !< species index in emission species list
    INTEGER(iwp), INTENT(IN) ::  vsrc_index  !< volume source index


    emis_index = vsrc_info(vsrc_index)%species_index
    data_index = emis_species_info(emis_index)%data_index

 END SUBROUTINE get_vsrc_indices


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Locate index in an array with matching integer.
!> integer version
!--------------------------------------------------------------------------------------------------!
 FUNCTION locate_index( array, target_value )  RESULT ( index )

    IMPLICIT NONE

    INTEGER(iwp)                           ::  index         !< array index; output
    INTEGER(iwp), INTENT(IN)               ::  target_value  !< target value

    INTEGER(iwp), INTENT(IN), DIMENSION(:) ::  array         !< input array


    index = MINLOC( ABS( array - target_value ), 1 )
    IF ( ( index < SIZE( array ) )  .AND.  ( array(index) /= target_value ) )  index = int_empty

 END FUNCTION locate_index


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Determine if trigger for pollen pool reset should be activated.
!--------------------------------------------------------------------------------------------------!
 FUNCTION pollen_pool_reset_trigger( )  RESULT ( reset_trigger )

    USE control_parameters,                                                                        &
        ONLY:  time_since_reference_point

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    IMPLICIT NONE

    INTEGER(iwp) ::  current_hour_utc  !< current model hour (UTC)

    LOGICAL      ::  reset_trigger  !< return value


    CALL get_date_time( time_since_reference_point, hour = current_hour_utc )
!
!-- Trigger reset only when it is armed and the hour is up. Rearm once the hour is past to prevent
!-- retrigger.
    reset_trigger = .FALSE.

    IF ( current_hour_utc == pollen_pool_reset_hour )  THEN
       IF ( pollen_pool_reset_trigger_armed )  THEN
          reset_trigger = .TRUE.
          pollen_pool_reset_trigger_armed = .FALSE.
       ENDIF
    ELSE
       pollen_pool_reset_trigger_armed = .TRUE.
    ENDIF

 END FUNCTION pollen_pool_reset_trigger


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Amalgamated release yield factor from precipitation. Currently a placeholder which returns unity.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_precip( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)                       ::  phi  !< yield factor; output
    REAL(wp), OPTIONAL, INTENT(IN) ::  xi   !< input; precipitation [mm]


    phi = 1.0_wp
    IF ( ( PRESENT( xi ) )  .OR.  ( .NOT. ignore_precipitation ) )  THEN
       phi = MAX( 0.0_wp, 1.0_wp - ( 2000.0_wp * xi ) )
    ENDIF

 END FUNCTION get_release_yield_precip


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Emission yield factor from turbulent kinetic energy.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_emission_yield_tke( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; turbulent kinetic energy [m2/s2]


    phi = logistic_function( xi, -2.1_wp, 4.0_wp )

 END FUNCTION get_emission_yield_tke


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Emission yield factor from relative humidity.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_emission_yield_rh( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; relative humidity [ ]


    phi = 1.0_wp
    IF ( xi > 0.90_wp )  THEN
       phi = 0.5_wp
    ELSEIF ( xi > 0.95_wp )  THEN
       phi = 0.0_wp
    ENDIF

 END FUNCTION get_emission_yield_rh


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Release yield factor from random mechanical processes.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_random( xi ) RESULT ( phi )

    USE control_parameters,                                                                        &
        ONLY:  dt_3d

    IMPLICIT NONE

    REAL(wp), PARAMETER  ::  log2 = 0.693147181   !< log(2) 
    REAL(wp), PARAMETER  ::  tiny = 1.0E-12       !< small number to prevent singularity

    REAL(wp)             ::  phi    !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi     !< input; species scaling factor [ ]


    phi = EXP( -( log2 * dt_3d ) / ( xi + tiny ) )

 END FUNCTION get_release_yield_random


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Release yield factor from plant coverage. Currently a placeholder which returns unity.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_coverage( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)                       ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN), OPTIONAL ::  xi   !< input; plant coverage


    IF ( ( .NOT. PRESENT( xi ) )  .OR.  ( PRESENT( xi ) ) )  phi = 1.0_wp

 END FUNCTION get_release_yield_coverage


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Release yield factor from plant altitude. Currently a placeholder which returns unity.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_altitude( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)                       ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN), OPTIONAL ::  xi   !< input; altitude [m]


    IF ( ( .NOT.  PRESENT(xi) )  .OR.  ( PRESENT(xi) ) )  phi = 1.0_wp

 END FUNCTION get_release_yield_altitude


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Release yield factor from temperature.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_temperature_betula( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; temperature [K]


    phi = 1.04_wp * logistic_function( xi, -0.277_wp,  76.0_wp ) *                                 &
                    logistic_function( xi,  0.45_wp, -137.0_wp )

 END FUNCTION get_release_yield_temperature_betula


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Release yield factor from temperature for alder trees.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_temperature_alder( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; temperature [K]


    phi = logistic_function( xi, -0.594_wp, 165.44_wp )  ! 2.2 * ( -0.27T + 75.2)

 END FUNCTION get_release_yield_temperature_alder


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Release yield factor from temperature for poaceae.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_temperature_poaceae( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; temperature [K]


    phi = logistic_function( xi,  -0.272_wp, 78.0_wp )

 END FUNCTION get_release_yield_temperature_poaceae


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Release yield factor from temperature for ambrosia. Empol 1.0 initial parameterization  
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_temperature_ambrosia( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; temperature [K]


    phi = 1.0_wp / ( 1.0_wp + EXP( -xi * 0.2_wp + 60.0_wp ) )

 END FUNCTION get_release_yield_temperature_ambrosia


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> release yield factor from relative humidity for betual and alder trees
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_rh_betula( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; relative humidity [ ]


    phi = logistic_function( xi, 21.0_wp, -15.0_wp )

 END FUNCTION get_release_yield_rh_betula


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> release yield factor from relative humidity for alder trees
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_rh_alder( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; relative humidity [ ]


    phi = logistic_function( xi, 21.0_wp, -15.0_wp )

 END FUNCTION get_release_yield_rh_alder


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> release yield factor from relative humidity for alder trees
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_rh_poaceae( xi ) RESULT ( phi )

    IMPLICIT NONE

    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; relative humidity [ ]

    REAL(wp) ::  w2 !< -(2w-1)^2

    w2 = -1.0_wp * ( 2.0_wp * xi - 1.0 ) * ( 2.0_wp * xi - 1.0_wp )

    phi = 2.0 * logistic_function( xi, -30.0_wp, 9.0_wp ) *                                        &
                EXP( w2 ) * ( 1.0 - 0.5 * EXP( w2 ) ) *                                            &
                logistic_function( xi, 45.0_wp, -40.0_wp )

 END FUNCTION get_release_yield_rh_poaceae


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> release yield factor from relative humidity for alder trees. Empol 1.0 initial parameterization  
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_release_yield_rh_ambrosia( xi ) RESULT ( phi )
 
    IMPLICIT NONE
 
    REAL(wp)             ::  phi  !< yield factor; output
    REAL(wp), INTENT(IN) ::  xi   !< input; relative humidity [K]


    phi = 1.0_wp / ( 1.0_wp + EXP( 20.0_wp * xi - 12.0_wp ) )

 END FUNCTION get_release_yield_rh_ambrosia


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates logistic function in the following form:
!> eta = 1 / [ 1 + exp( alpha1*xi + alpha0 ) ]
!--------------------------------------------------------------------------------------------------!
 FUNCTION logistic_function( xi, alpha1, alpha0 ) RESULT ( eta )

    IMPLICIT NONE

    REAL(wp) ::  alpha0  !< parameter for x^0
    REAL(wp) ::  alpha1  !< parameter for x^1
    REAL(wp) ::  eta     !< functional output
    REAL(wp) ::  xi      !< independent variable


    eta = 1.0 / ( 1.0 + EXP( alpha1 * xi + alpha0 ) )

 END FUNCTION logistic_function


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Estimates total TKE (resolved + SGS) from SGS TKE.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_tke_total_default( tke_sgs, phi_sgs ) RESULT ( tke_total )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) ::  tke_sgs    !< subgrid scale TKE; input
    REAL(wp)             ::  tke_total  !< resolved scale TKE; output
    REAL(wp), INTENT(IN) ::  phi_sgs    !< fraction of TKE in SGS


    tke_total = ( 1.0_wp / ( phi_sgs + 1.0E-12_wp ) ) * tke_sgs

 END FUNCTION get_tke_total_default


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Estimates total TKE (resolved + SGS) using ad hoc estimation of resolved scale TKE using
!> gradient transport.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_tke_total_adhoc( k_vsrc, phi_sgs ) RESULT ( tke_total )

    USE arrays_3d,                                                                                 &
        ONLY:  ddzu,                                                                               &
               e,                                                                                  &
               u,                                                                                  &
               v,                                                                                  &
               w 

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  k_vsrc  !< volume source inde ; input

    INTEGER(iwp) ::  i  !< cell indices in x direction
    INTEGER(iwp) ::  j  !< cell indices in y direction
    INTEGER(iwp) ::  k  !< cell indices in z direction

    REAL(wp) ::  dudx        !< u velocity gradient magnitude
    REAL(wp) ::  dvdy        !< v velocity gradient magnitude
    REAL(wp) ::  dwdz        !< w velocity gradient magnitude
    REAL(wp) ::  phi_uprime  !< random offset for u turbulent fluctuation
    REAL(wp) ::  phi_vprime  !< random offset for v turbulent fluctuation
    REAL(wp) ::  phi_wprime  !< random offset for w turbulent fluctuation
    REAL(wp) ::  tke_grad    !< estimate of total TKE using gradient transport
    REAL(wp) ::  tke_total   !< resolved scale TKE; output

    REAL(wp), INTENT(IN) ::  phi_sgs       !< fraction fo TKE in SGS


!
!-- Cell coordinates of volume source (limit z to first interior vertical layer).
    i = vsrc_pos(k_vsrc)%i
    j = vsrc_pos(k_vsrc)%j
    k = MAX( nzb+1,vsrc_pos(k_vsrc)%k )

    CALL RANDOM_NUMBER( phi_uprime )
    CALL RANDOM_NUMBER( phi_vprime )
    CALL RANDOM_NUMBER( phi_wprime )

    dudx = ddx     * ABS( u(k,j,i+1) - u(k,j,i) ) * ABS( 0.5_wp + phi_uprime )
    dvdy = ddy     * ABS( v(k,j+1,i) - v(k,j,i) ) * ABS( 0.5_wp + phi_vprime )
    dwdz = ddzu(k) * ABS( w(k+1,j,i) - w(k,j,i) ) * ABS( 0.5_wp + phi_wprime )
    
    tke_grad = 0.5_wp * ( dudx + dvdy + dwdz ) 
    tke_total = MIN( get_tke_total_default( e(k,j,i), phi_sgs ), tke_grad )

 END FUNCTION get_tke_total_adhoc


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Estimates total TKE (resolved + SGS) using ad hoc estimation of resolved scale TKE using
!> gradient transport.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_tke_total_dynamic( tke_sgs, k, j, i ) RESULT (tke_total)

    USE arrays_3d,                                                                                 &
        ONLY:  u,                                                                                  &
               v,                                                                                  &
               w

    USE indices,                                                                                   &
        ONLY:  ngp_2dh_outer,                                                                      &
               topo_flags

    USE statistics,                                                                                &
        ONLY:  flow_statistics_called,                                                             &
               hom,                                                                                &
               sums,                                                                               &
               sums_l

    IMPLICIT  NONE

    INTEGER(iwp), INTENT(IN) ::  i  !<
    INTEGER(iwp), INTENT(IN) ::  j  !<
    INTEGER(iwp), INTENT(IN) ::  k  !<

    REAL(wp) ::  flag1              !< flag to mask topography
    REAL(wp) ::  tke_rslvd_dom_avg  !<
    REAL(wp) ::  tke_sgs            !<
    REAL(wp) ::  tke_total          !<


    IF ( .NOT.  flow_statistics_called )  THEN
!
!--    First calculate horizontally averaged profiles of the horizontal velocities.
       sums_l(:,1,0) = 0.0_wp
       sums_l(:,2,0) = 0.0_wp
!
!--    Flag indicating vicinity of wall.
       flag1 = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 24 ) )
       sums_l(k,1,0) = sums_l(k,1,0) + u(k,j,i) * flag1
       sums_l(k,2,0) = sums_l(k,2,0) + v(k,j,i) * flag1
       sums(:,1) = sums_l(:,1,0)
       sums(:,2) = sums_l(:,2,0)
!
!--    Final values are obtained by division by the total number of grid points used for the
!--    summation.
       hom(:,1,1,0) = sums(:,1) / ngp_2dh_outer(:,0)   !< u velocity component
       hom(:,1,2,0) = sums(:,2) / ngp_2dh_outer(:,0)   !< v velocity component
!
!--    Now calculate the profiles of SGS TKE and the resolved-scale velocity variances.
       sums_l(:,8,0)  = 0.0_wp
       sums_l(:,30,0) = 0.0_wp
       sums_l(:,31,0) = 0.0_wp
       sums_l(:,32,0) = 0.0_wp
!
!--    Flag indicating vicinity of wall.
       flag1 = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 24 ) )

       sums_l(k,8,0)  = sums_l(k,8,0)  + tke_sgs                        * flag1
       sums_l(k,30,0) = sums_l(k,30,0) + ( u(k,j,i) - hom(k,1,1,0) )**2 * flag1
       sums_l(k,31,0) = sums_l(k,31,0) + ( v(k,j,i) - hom(k,1,2,0) )**2 * flag1
       sums_l(k,32,0) = sums_l(k,32,0) +   w(k,j,i)**2                  * flag1

       sums(:,8)  = sums_l(:,8,0)
       sums(:,30) = sums_l(:,30,0)
       sums(:,31) = sums_l(:,31,0)
       sums(:,32) = sums_l(:,32,0)
!
!--    Final values are obtained by division by the total number of grid points used for the
!--    summation.
       hom(:,1,8,0)  = sums(:,8)  / ngp_2dh_outer(:,0)   ! e 
       hom(:,1,30,0) = sums(:,30) / ngp_2dh_outer(:,0)   ! u*2
       hom(:,1,31,0) = sums(:,31) / ngp_2dh_outer(:,0)   ! v*2
       hom(:,1,32,0) = sums(:,32) / ngp_2dh_outer(:,0)   ! w*2
    ENDIF

    tke_rslvd_dom_avg = ( hom(k,1,30,0) + hom(k,1,31,0) + hom(k,1,32,0) ) * 0.5

    IF ( hom(k,1,8,0) <= 0 )  hom(k,1,8,0) = 1.0E-5
    
    tke_total = ( tke_rslvd_dom_avg / hom(k,1,8,0) * tke_sgs ) + tke_sgs
       
 END FUNCTION get_tke_total_dynamic


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates dynamic viscosity [Pa s] of a vertical column of air
!> as a function of temperature using Sutherland's Law.
!> NOTE 1: equation values are hard-coded for air
!> NOTE 2: mu is a very weak function of tempertaure so vertical thermal stratification
!>         has a much stronger effect than horizontal temperature variations
!>         using the ideal gas law on a single column is thus much more efficient than
!>         calculating pt x exner for every column in the domain with little impact on accuracy
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_dynamic_viscosity_sutherland( ) RESULT ( mu )

    USE arrays_3d,                                                                                 &
        ONLY:  hyp,                                                                                &
               rho_air

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  r_d

    IMPLICIT NONE

    REAL(wp), DIMENSION(nzb:nzt+1) ::  mu           !< settling velocity; output
    REAL(wp), DIMENSION(nzb:nzt+1) ::  temperature  !< thermodynamic temperature

!
!-- Calculate temperature.
    temperature = ( hyp / rho_air ) / r_d
!
!-- Sutherland's law ( T x sqrt(T) is faster than T**1.5 ).
    mu = 1.458E-6_wp * ( temperature * SQRT( temperature ) ) / ( 110.4_wp + temperature )

 END FUNCTION get_dynamic_viscosity_sutherland

 END MODULE chem_emis_pollen_mod
