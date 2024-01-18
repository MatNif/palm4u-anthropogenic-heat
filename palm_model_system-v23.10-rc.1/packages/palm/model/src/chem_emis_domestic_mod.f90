!> @file chem_emis_domestic_mod.f90
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
! Authors:
! --------
!> @author Edward C. Chan
!> @author Ilona Jäkel
!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Modules for generic emission mode for user-defined emission data.
!> Also contains interface functions and subroutines for other emission mode modules.
!--------------------------------------------------------------------------------------------------!
 MODULE chem_emis_domestic_mod

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  degc_to_k

    USE chem_emis_generic_mod  ! for derived type chem_emis_species

    USE chem_emis_vsrc_mod     ! for derived type chem_emis_vsrc_pos

    USE chem_modules,                                                                              &          
        ONLY:  emis_domestic_max_bld_types,                                                        &
               nc_field_length
    
    USE control_parameters,                                                                        &
        ONLY:  coupling_char
    
    USE kinds

    IMPLICIT NONE
    SAVE 
    PRIVATE

!
!-- General characteristics.
    CHARACTER(LEN=*), PARAMETER ::  input_nc_file = 'PIDS_EMIS_DOMESTIC'  !< emission mode input
    INTEGER(iwp)                ::  lod                                   !< nc lod

!
!-- Timestamps.
    CHARACTER(LEN=nc_field_length)                            ::  current_time  !< simulation time(time)
    CHARACTER(LEN=nc_field_length), ALLOCATABLE, DIMENSION(:) ::  timestamp     !< individual timestamps

    INTEGER(iwp) ::  num_timestamp  !< # timestamps
    INTEGER(iwp) ::  time_index     !< current timestamp

!
!-- Interface to KPP chemical mechanism.
    INTEGER(iwp)                                       ::  num_emis_species  !< # emission species
    TYPE(chem_emis_species), ALLOCATABLE, DIMENSION(:) ::  emis_species      !< emission species

!
!-- Interface to chem_emis_vsrc_mod.
    INTEGER(iwp)                                          ::  num_vsrc_pos     !< # volume sources
    INTEGER(iwp),             ALLOCATABLE, DIMENSION(:)   ::  vsrc_nc_index    !< location on nc file
    REAL(KIND=dp),            ALLOCATABLE, DIMENSION(:,:) ::  vsrc_emis_value  !< volume source volues
    TYPE(chem_emis_vsrc_pos), ALLOCATABLE, DIMENSION(:)   ::  vsrc_pos         !< volume source positions

!
!-- Update interval and trigger for LOD != 2.
    REAL(wp) ::  update_interval = 300.0_wp  !< time between updates [s]
    REAL(wp) ::  time_since_last_update      !< time since last update

!
!-- LOD 0 parameterization input variables (default values hard coded).
    INTEGER(iwp) ::  k_ambient = 10    !< number of vertical layers for determining domain ambient temperature

    REAL(wp) ::  annual_heating_degree = 2100.0_wp            !< heating degrees
    REAL(wp) ::  base_temperature      = 15.0_wp + degc_to_k  !< base temperature [K]

    REAL(wp), DIMENSION(emis_domestic_max_bld_types) ::  FC =   &  !< compactness factor for each building type
            (/ 0.23_wp, 0.28_wp, 0.28_wp, 0.26_wp, 0.29_wp, 0.29_wp /) 
    REAL(wp), DIMENSION(emis_domestic_max_bld_types) ::  kW =   &  !< energy demand for each building type
            (/ 130.0_wp, 100.0_wp, 100.0_wp, 110.0_wp, 89.0_wp, 89.0_wp /)
    
!
!-- LOD 0 parameterization derived variables .
    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  heat_degree_consumption   !< energy consumption per degree temperature [J/K]
    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  species_emission_factors  !< species emission factors [kg/J]


    INTERFACE chem_emis_domestic_init
       MODULE PROCEDURE chem_emis_domestic_init
    END INTERFACE chem_emis_domestic_init

    INTERFACE chem_emis_domestic_cleanup
       MODULE PROCEDURE chem_emis_domestic_cleanup
    END INTERFACE chem_emis_domestic_cleanup

    INTERFACE chem_emis_domestic_update
       MODULE PROCEDURE chem_emis_domestic_update
    END INTERFACE chem_emis_domestic_update

!
!-- Public methods specific to generic emission mode (i.e., this module).
    PUBLIC ::  chem_emis_domestic_cleanup, chem_emis_domestic_init,                                &
               chem_emis_domestic_update

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Wrapper for intitializing emission mode.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_domestic_init( )

    USE chem_modules,                                                                              &
        ONLY:  emis_domestic_lod

    IMPLICIT NONE


    lod = emis_domestic_lod
!
!-- Initialize emission mode based on LOD and update emission source for first iteration.
    SELECT CASE ( lod )
    CASE ( 0 )
       CALL emissions_domestic_init_lod0( )
       CALL emissions_vsrc_update_lod0( )
    CASE ( 1 )
       CALL emissions_domestic_init_lod1( )
       CALL emissions_vsrc_update_lod1( )
    CASE DEFAULT
       CALL emissions_domestic_init_lod2( )
       CALL emissions_vsrc_update_lod2( )
    END SELECT

 END SUBROUTINE chem_emis_domestic_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Extracts all volume source positions in subdomain from emissions netCDF file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_domestic_cleanup( )

    USE chem_emis_vsrc_mod

    IMPLICIT NONE


    IF ( ALLOCATED( timestamp                ) )  DEALLOCATE( timestamp                )
    IF ( ALLOCATED( emis_species             ) )  DEALLOCATE( emis_species             )
    IF ( ALLOCATED( vsrc_pos                 ) )  DEALLOCATE( vsrc_pos                 )
    IF ( ALLOCATED( vsrc_emis_value          ) )  DEALLOCATE( vsrc_emis_value          )
    IF ( ALLOCATED( vsrc_nc_index            ) )  DEALLOCATE( vsrc_nc_index            )
    IF ( ALLOCATED( species_emission_factors ) )  DEALLOCATE( species_emission_factors )
    IF ( ALLOCATED( heat_degree_consumption  ) )  DEALLOCATE( heat_degree_consumption  )

 END SUBROUTINE chem_emis_domestic_cleanup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates emissions at every time step (called in time_integratin).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_domestic_update( )

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_update_source

    IMPLICIT NONE

    INTEGER(iwp) ::  k  !< generic counter


    IF ( ( num_vsrc_pos == 0 )  .OR.  ( num_emis_species == 0 ) )  RETURN

    SELECT CASE ( lod )
       CASE ( 0 )
          CALL emissions_vsrc_update_lod0( )
       CASE ( 1 )
          CALL emissions_vsrc_update_lod1( )
       CASE DEFAULT
          CALL emissions_vsrc_update_lod2( )
    END SELECT

    DO  k = 1, num_emis_species
       CALL chem_emis_vsrc_update_source( emis_species(k)%mech_index, vsrc_pos,                    &
                                          vsrc_emis_value(k,:) )
    ENDDO

 END SUBROUTINE chem_emis_domestic_update


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize emission mode module under LOD 0.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emissions_domestic_init_lod0( )

    CALL location_message( 'reading LOD 0 domestic emissions data', 'start' )

    CALL read_namelist_lod0( )
    CALL read_nc_stack_data_lod0( )

!
!-- Allocate space for local volume source emissions.
    ALLOCATE( vsrc_emis_value( num_emis_species, num_vsrc_pos ) )
    vsrc_emis_value = 0.0

    time_since_last_update = update_interval
    
    CALL location_message( 'reading LOD 0 domestic emissions data', 'finished' )

 END SUBROUTINE emissions_domestic_init_lod0


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize emission mode module under LOD 1.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emissions_domestic_init_lod1( )

    ! right now it does nothing

 END SUBROUTINE emissions_domestic_init_lod1


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize emission mode module under LOD 2.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emissions_domestic_init_lod2( )

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_init_volume_source,                                                       &
               chem_emis_nc_get_species,                                                           &
               chem_emis_nc_get_time,                                                              &
               chem_emis_nc_get_vsrc_positions

    USE chem_modules

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               open_read_file

    IMPLICIT NONE

    INTEGER(iwp) ::  ncid  !< nc file handle


    CALL location_message( 'reading emissions data from ' // TRIM( input_nc_file ), 'start' )
!
!-- Open netCDF file and assign LOD (should be 2).
    CALL open_read_file( TRIM( input_nc_file ) //  TRIM( coupling_char ), ncid )
    CALL get_attribute( ncid, nc_att_lod, lod, .TRUE. )

!
!-- Extract timestamps and determine current time.
    CALL chem_emis_nc_get_time( num_timestamp, timestamp, ncid )
    time_index = chem_emis_get_current_time( current_time, timestamp )

!
!-- Extract emitting species in mechanism.
    CALL chem_emis_nc_get_species( num_emis_species, emis_species, ncid )

!
!-- Extract volume source locations.
    CALL chem_emis_nc_get_vsrc_positions( num_vsrc_pos, vsrc_pos, vsrc_nc_index, ncid )

!
!-- Initialize volume source data and global volume source data structure.
    CALL chem_emis_init_volume_source( vsrc_pos )

!
!-- Allocate space for local volume source emissions.
    ALLOCATE( vsrc_emis_value( num_emis_species, num_vsrc_pos ) )
    vsrc_emis_value = 0.0

    CALL close_input_file( ncid )  ! closes netCDF file
    CALL location_message( 'reading emissions data from ' // TRIM(input_nc_file), 'finished' )

 END SUBROUTINE emissions_domestic_init_lod2


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Update emission mode module under LOD 0.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emissions_vsrc_update_lod0( )

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_update_trigger

    IMPLICIT NONE

    INTEGER(iwp) ::  i      !< generic indices
    INTEGER(iwp) ::  k      !< generic indices

    REAL(wp), PARAMETER ::  day_to_seconds = 86400.0  !< seconds in a day
    REAL(wp)            ::  heating_rate              !< building heating rate   
    REAL(wp)            ::  temperature_deficit       !< indoor to ambient temperature deficit


!
!-- Check whether it is time to update source terms.
    IF ( .NOT. chem_emis_update_trigger( time_since_last_update, update_interval ) )  RETURN

!
!-- Assume no heating required unless temperature deficit existis.
    vsrc_emis_value = 0.0_wp
    temperature_deficit = base_temperature - get_domain_ambient_temperature( k_ambient )

    IF ( temperature_deficit > 0.0_wp )  THEN
       DO  k  = 1, num_vsrc_pos
          heating_rate = ( heat_degree_consumption(k) / temperature_deficit ) / day_to_seconds
          DO  i = 1, num_emis_species
              vsrc_emis_value(i,k) = species_emission_factors(i) * heating_rate
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE emissions_vsrc_update_lod0


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Update emission mode module under LOD 1.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emissions_vsrc_update_lod1( )

    ! right now it does nothing

 END SUBROUTINE emissions_vsrc_update_lod1


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates emission volume source under LOD 2 (only mode).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emissions_vsrc_update_lod2( )

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_nc_get_vsrc_species_var_name

    USE chem_modules

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  close_input_file,                                                                   &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               open_read_file

    IMPLICIT NONE

    CHARACTER(LEN=nc_field_length) ::  species_var_name  !< netCDF emission variable name

    INTEGER(iwp) ::  buf_nvsrc  !< # volume sources positions in nc file
    INTEGER(iwp) ::  i          !< counter
    INTEGER(iwp) ::  k          !< counter
    INTEGER(iwp) ::  ncid       !< nc file handle

    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) ::  buf_emissions  !< emission values

!
!-- Open netCDF file and get dimension.
    CALL open_read_file( TRIM( input_nc_file ) // TRIM( coupling_char ), ncid )
    CALL get_dimension_length( ncid, buf_nvsrc, nc_dim_nvsrc )
    ALLOCATE( buf_emissions( buf_nvsrc ) )

!
!-- For each emitting species appearing in the chemical mechanism determine the volume source
!-- variable name ( 'vsrc_' + species name ) and grab the emission values for the current time
!-- stamp.
    DO  k = 1, num_emis_species

       CALL chem_emis_nc_get_vsrc_species_var_name( species_var_name, emis_species(k)%name_str )
       CALL get_variable( ncid, species_var_name, buf_emissions, time_index, buf_nvsrc, .TRUE. )

       DO  i = 1, num_vsrc_pos
          vsrc_emis_value(k,i) = buf_emissions( vsrc_nc_index(i) )
       ENDDO

    ENDDO

!
!-- Cleaning up.
    DEALLOCATE( buf_emissions )
    CALL close_input_file( ncid )

 END SUBROUTINE emissions_vsrc_update_lod2


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads namelist items for LOD 0 parameters.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE read_namelist_lod0( )

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  degc_to_k

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_species_match

    USE chem_gasphase_mod,                                                                         &
        ONLY:  spc_names

    USE chem_modules,                                                                              &
        ONLY:  emis_domestic_base_temperature,                                                     &
               emis_domestic_compact_factors,                                                      &
               emis_domestic_energy_demands,                                                       &
               emis_domestic_heating_degree,                                                       &
               emis_domestic_max_bld_types,                                                        &
               emis_domestic_sampling_k,                                                           &
               emis_domestic_species_emission_factors,                                             &
               emis_domestic_species_names,                                                        &
               emis_domestic_update_interval

    USE indices,                                                                                   &
        ONLY:  nzt

    IMPLICIT NONE

!
!-- Constants.
    INTEGER(iwp), PARAMETER :: prefix_pm_len = 2  !< position of PM prefix

    CHARACTER(LEN=prefix_pm_len), PARAMETER ::  prefix_pm = 'PM'  !<

    REAL(wp), PARAMETER ::  molkg2umolmg = 1.0E+06_wp  !< mol/kg --> umol/mg
    REAL(wp), PARAMETER ::  j2tj         = 1.0E-12_wp  !< J --> TJ

!
!-- Local variables.
    INTEGER(iwp) ::  k       !< generic index
    INTEGER(iwp) ::  k_mech  !< KPP species index
    INTEGER(iwp) ::  k_user  !< user species index

    REAL(wp) ::  conv_factor  !< conversion factor


!
!-- Validate and transfer namelist parameters to module input variables.
    IF ( ( emis_domestic_sampling_k > 0 )  .AND.  ( emis_domestic_sampling_k < nzt ) )  THEN
       k_ambient = emis_domestic_sampling_k
    ENDIF

    IF ( emis_domestic_heating_degree > 0.0_wp  )                                                  &
       annual_heating_degree = emis_domestic_heating_degree

    IF ( emis_domestic_base_temperature > -degc_to_k )                                             &
       base_temperature = emis_domestic_base_temperature  + degc_to_k

    IF ( emis_domestic_update_interval > 0.0_wp )                                                  &
       update_interval = emis_domestic_update_interval

!
!-- Building type specific data (only uses namelist values if they are non-gevative).
    DO  k = 1, emis_domestic_max_bld_types
        IF ( emis_domestic_compact_factors(k) > 0.0_wp )  fc(k) = emis_domestic_compact_factors(k)
        IF ( emis_domestic_energy_demands(k)  > 0.0_wp )  kw(k) = emis_domestic_energy_demands(k)
    ENDDO

!
!-- Link volume source species to mechanism and user input.
    CALL chem_emis_species_match( num_emis_species, emis_species, emis_domestic_species_names )

    IF ( num_emis_species == 0 )  RETURN

!
!-- Assign emission factors with unit conversion to KPP molar units per TJ for gas phase species.
    ALLOCATE( species_emission_factors(num_emis_species) )
    DO  k = 1, num_emis_species
        k_mech = emis_species(k)%mech_index
        k_user = emis_species(k)%user_index
        conv_factor = j2tj
        IF ( .NOT. ( spc_names(k_mech)(1:prefix_pm_len) == prefix_pm ) )  THEN
           conv_factor = conv_factor * molkg2umolmg
        ENDIF
        species_emission_factors(k) = conv_factor * MAX( 0.0_wp,                                   &
                                                    emis_domestic_species_emission_factors(k_user) )
    ENDDO

 END SUBROUTINE read_namelist_lod0


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads static file items for LOD 0 parameters.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE read_nc_stack_data_lod0 ( )

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_init_volume_source

    USE control_parameters,                                                                        &
        ONLY:  coupling_char

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nyn,                                                                                &
               nxr,                                                                                &
               nys,                                                                                &
               topo_top_ind

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  close_input_file,                                                                   &
               get_variable,                                                                       &
               input_file_static,                                                                  &
               open_read_file

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !< generic indices
    INTEGER(iwp) ::  j     !< generic indices
    INTEGER(iwp) ::  k     !< generic indices
    INTEGER(iwp) ::  ncid  !< netcdf file handle
    INTEGER(iwp) ::  typ   !< stack building type

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  building_type  !< building type

    REAL(wp), ALLOCATABLE, DIMENSION(:,:)     ::  building_height  !< building height (building_2d)
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)     ::  stack_volume     !< building volume at stack location


!
!-- Allocate memory for static file variables.
    ALLOCATE( stack_volume(nys:nyn,nxl:nxr)    )
    ALLOCATE( building_height(nys:nyn,nxl:nxr) )
    ALLOCATE( building_type(nys:nyn,nxl:nxr)   )

!
!-- Transfer data from static file.
    CALL open_read_file( TRIM(input_file_static) // TRIM(coupling_char), ncid )
    CALL get_variable( ncid, 'stack_building_volume', stack_volume, nxl, nxr, nys, nyn )
    CALL get_variable( ncid, 'buildings_2d',  building_height, nxl, nxr, nys, nyn )
    CALL get_variable( ncid, 'building_type', building_type,   nxl, nxr, nys, nyn )
    CALL close_input_file( ncid )

!
!-- Tally up volume sources.
    num_vsrc_pos = 0

    DO  j = nys, nyn
       DO  i = nxl, nxr
          IF ( stack_volume(j,i) > 0.0_wp )  num_vsrc_pos = num_vsrc_pos + 1
       ENDDO
    ENDDO

!
!-- Quit (after cleaning up) when there are no stacks in domain.
    IF ( num_vsrc_pos == 0 )  THEN
       DEALLOCATE( stack_volume    )
       DEALLOCATE( building_height )
       DEALLOCATE( building_type   )
       RETURN
    ENDIF

!
!-- Allocate and populate volume source positions.
    ALLOCATE( vsrc_pos(num_vsrc_pos)                )
    ALLOCATE( heat_degree_consumption(num_vsrc_pos) )

    k = 0
    heat_degree_consumption = 0.0_wp

    DO  j = nys, nyn
       DO  i = nxl, nxr
           IF ( stack_volume(j,i) > 0.0_wp )  THEN
              k = k + 1
              vsrc_pos(k)%i = i
              vsrc_pos(k)%j = j
              vsrc_pos(k)%k = topo_top_ind(j,i,0) + 1
              typ           = building_type(j,i)
              heat_degree_consumption(k) = get_hdc( stack_volume(j,i), typ )
          ENDIF
       ENDDO
    ENDDO

!
!-- Initialize volume source data and global volume source data structure.
    CALL chem_emis_init_volume_source( vsrc_pos )

!
!-- Cleaning up.
    DEALLOCATE(building_type)
    DEALLOCATE(building_height)
    DEALLOCATE(stack_volume)

 END SUBROUTINE read_nc_stack_data_lod0


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Obtains heating degree consumption.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_hdc( volume, ktype )  RESULT ( hdc )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  ktype  !< building type

    REAL(wp), INTENT(IN) ::  volume  !< object volume [m3]

    REAL(wp) ::  hdc  !< output heat degree consumption [J/K]

    REAL(wp), PARAMETER ::  kwh2j  = 3.6E6_wp  !< from kWh to J

    hdc = kwh2j * KW(ktype) * FC(ktype) * volume / annual_heating_degree

 END FUNCTION get_hdc


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Obtains heating time from volume of object.
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_domain_ambient_temperature( ksample )  RESULT( domain_temperature )

    USE arrays_3d,                                                                                 &
        ONLY:  exner,                                                                              &
               pt

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nyn,                                                                                &
               nxr,                                                                                &
               nys,                                                                                &
               topo_top_ind

    INTEGER(iwp), INTENT(IN) ::  ksample  !< number vertical layers for sampling

    REAL(wp) ::  domain_temperature  !< domain mean temperature
    
    INTEGER(iwp) ::  i   !< generic indicies
    INTEGER(iwp) ::  j   !< generic indicies
    INTEGER(iwp) ::  k   !< generic indicies
    INTEGER(iwp) ::  kt  !< generic indicies

    REAL(wp) ::  ncells  !< number of cells
    
    ncells = 1.0E-10_wp
    domain_temperature = 0.0_wp

    DO  j = nys, nyn
       DO  i = nxl, nxr
          DO  k = 1, ksample
              kt = topo_top_ind(j,i,0) + k
              domain_temperature = domain_temperature + pt(kt,j,i) * exner(kt)
              ncells = ncells + 1.0_wp
          ENDDO
       ENDDO
    ENDDO

    domain_temperature = domain_temperature / ncells

 END FUNCTION get_domain_ambient_temperature

 END MODULE chem_emis_domestic_mod
