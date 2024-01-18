!> @file chem_emis_traffic_mod.f90
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
!> Also contains interface functions and subroutines for other emission mode modules
!--------------------------------------------------------------------------------------------------!
 MODULE chem_emis_traffic_mod

    USE chem_emis_generic_mod  ! for derived type chem_emis_species   
 
    USE chem_emis_vsrc_mod     ! for derived type chem_emis_vsrc_pos

    USE chem_modules,                                                                              &
        ONLY:  nc_field_length
    
    USE control_parameters,                                                                        &
        ONLY:  coupling_char

    USE kinds

    IMPLICIT NONE
    SAVE 
    PRIVATE

!
!-- general characteristics
    CHARACTER(LEN=*), PARAMETER ::  input_nc_file = 'PIDS_EMIS_TRAFFIC'         !< emission mode input
    INTEGER(iwp)                ::  lod                                         !< nc lod

!
!-- timestamps
    CHARACTER(LEN=nc_field_length)                            ::  current_time  !< simulation time(time)
    CHARACTER(LEN=nc_field_length), ALLOCATABLE, DIMENSION(:) ::  timestamp     !< individual timestamps

    INTEGER(iwp) :: num_timestamp     !< # timestamps
    INTEGER(iwp) :: time_index        !< current timestamp

!
!-- interface to KPP chemical mechanism
    INTEGER(iwp)                                       ::  num_emis_species     !< # emission species
    TYPE(chem_emis_species), ALLOCATABLE, DIMENSION(:) ::  emis_species         !< emission species

!
!-- interface to chem_emis_vsrc_mod
    INTEGER(iwp)                                          ::  num_vsrc_pos      !< # volume sources
    INTEGER(iwp),             ALLOCATABLE, DIMENSION(:)   ::  vsrc_nc_index     !< location on nc file
    REAL(KIND=dp),            ALLOCATABLE, DIMENSION(:,:) ::  vsrc_emis_value   !< volume source volues
    TYPE(chem_emis_vsrc_pos), ALLOCATABLE, DIMENSION(:)   ::  vsrc_pos          !< volume source positions

!
!-- module interface
    INTERFACE chem_emis_traffic_init
      MODULE PROCEDURE chem_emis_traffic_init
    END INTERFACE chem_emis_traffic_init

    INTERFACE chem_emis_traffic_cleanup
      MODULE PROCEDURE chem_emis_traffic_cleanup
    END INTERFACE chem_emis_traffic_cleanup

    INTERFACE chem_emis_traffic_update
      MODULE PROCEDURE chem_emis_traffic_update
    END INTERFACE chem_emis_traffic_update

!
!-- public methods specific to generic emission mode (i.e., this module)
    PUBLIC :: chem_emis_traffic_cleanup, chem_emis_traffic_init,                                   &
              chem_emis_traffic_update

 CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  PUBLIC METHODS SPECIFIC TO GENERIC EMISSION MODE
!!  WRAPPERS TO VARIOUS LOD FUNCTIONALITIES
!!  PREFIX chem_emis_domestic
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> wrapper for intitializing emission mode
!--------------------------------------------------------------------------------------------------!


 SUBROUTINE chem_emis_traffic_init( )

    USE chem_modules,                                                                              &
        ONLY:  emis_traffic_lod

    IMPLICIT NONE

    lod = emis_traffic_lod

    SELECT CASE ( lod )                            ! initialize emission mode based on LOD
    CASE ( 0 )
       CALL emissions_traffic_init_lod0 ( )
       CALL emissions_vsrc_update_lod0( )
    CASE ( 1 )
       CALL emissions_traffic_init_lod1( )
       CALL emissions_vsrc_update_lod1( )
    CASE DEFAULT
       CALL emissions_traffic_init_lod2( )
       CALL emissions_vsrc_update_lod2( )          ! update volume sources for first iteration
    END SELECT


 END SUBROUTINE chem_emis_traffic_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> extracts all volume source positions in subdomain from emissions netCDF file
!--------------------------------------------------------------------------------------------------!


 SUBROUTINE chem_emis_traffic_cleanup( )

    USE chem_emis_vsrc_mod

    IMPLICIT NONE

    IF ( ALLOCATED( timestamp)       )  DEALLOCATE( timestamp )
    IF ( ALLOCATED( emis_species)    )  DEALLOCATE( emis_species )
    IF ( ALLOCATED( vsrc_pos)        )  DEALLOCATE( vsrc_pos )
    IF ( ALLOCATED( vsrc_emis_value) )  DEALLOCATE( vsrc_emis_value )
    IF ( ALLOCATED( vsrc_nc_index)   )  DEALLOCATE( vsrc_nc_index )

 END SUBROUTINE chem_emis_traffic_cleanup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates emissions at every tiem step (called in time_integratin)
!--------------------------------------------------------------------------------------------------!


 SUBROUTINE chem_emis_traffic_update( )

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_update_source

    IMPLICIT NONE

    INTEGER(iwp) ::  k  !< generic counter

    SELECT CASE ( lod )
       CASE ( 0 )
          CALL emissions_vsrc_update_lod0( )
       CASE ( 1 )
          CALL emissions_vsrc_update_lod1( )
       CASE DEFAULT
          CALL emissions_vsrc_update_lod2( )
    END SELECT

    DO  k = 1, num_emis_species
       CALL chem_emis_vsrc_update_source ( emis_species(k)%mech_index, vsrc_pos, vsrc_emis_value(k,:) )
    END DO

 END SUBROUTINE chem_emis_traffic_update


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  PRIVATE METHODS FOR TRAFFIC EMISSION MODE
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> initialize emission mode module under LOD 0
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE emissions_traffic_init_lod0( )
    ! right now it does nothing
 END SUBROUTINE emissions_traffic_init_lod0

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> initialize emission mode module under LOD 1
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE emissions_traffic_init_lod1( )
    ! right now it does nothing
 END SUBROUTINE emissions_traffic_init_lod1

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> initialize emission mode module under LOD 2
!--------------------------------------------------------------------------------------------------!
 

 SUBROUTINE emissions_traffic_init_lod2( )

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


    CALL location_message( 'reading emissions data from ' // TRIM(input_nc_file), 'start' )
    CALL open_read_file( TRIM( input_nc_file ) // TRIM( coupling_char ), ncid )
    CALL get_attribute( ncid, nc_att_lod, lod, .TRUE. )   ! assign LOD (should be 2)

!
!-- extract timestamps and determine current time
    CALL chem_emis_nc_get_time( num_timestamp, timestamp, ncid )
    time_index = chem_emis_get_current_time( current_time, timestamp )

!
!-- extract emitting species in mechanism
    CALL chem_emis_nc_get_species( num_emis_species, emis_species, ncid )

!
!-- extract volume source locations
    CALL chem_emis_nc_get_vsrc_positions( num_vsrc_pos, vsrc_pos, vsrc_nc_index, ncid )

!
!-- initialize volume source data and global volume source data structure
    CALL chem_emis_init_volume_source( vsrc_pos )

!
!-- allocate space for local volume source emissions
    ALLOCATE( vsrc_emis_value( num_emis_species,num_vsrc_pos ) )
    vsrc_emis_value = 0.0

    CALL close_input_file( ncid )  ! closes netCDF file
    CALL location_message( 'reading emissions data from ' // TRIM(input_nc_file), 'finished' )

 END SUBROUTINE emissions_traffic_init_lod2


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> update emission mode module under LOD 0
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE emissions_vsrc_update_lod0( )
    ! right now it does nothing
 END SUBROUTINE emissions_vsrc_update_lod0

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> update emission mode module under LOD 1
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE emissions_vsrc_update_lod1( )
    ! right now it does nothing
 END SUBROUTINE emissions_vsrc_update_lod1

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates emission volume source under LOD 2 (only mode)
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

    INTEGER(iwp) ::  buf_nvsrc         !< # volume sources positions in nc file
    INTEGER(iwp) ::  i                 !< counter
    INTEGER(iwp) ::  k                 !< counter
    INTEGER(iwp) ::  ncid              !< nc file handle

    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) ::  buf_emissions     !< emission values

!
!-- setting everything up
    CALL open_read_file( TRIM( input_nc_file ) // TRIM( coupling_char ), ncid )
    CALL get_dimension_length( ncid, buf_nvsrc, nc_dim_nvsrc )  ! gets dimension
    ALLOCATE( buf_emissions( buf_nvsrc ) )

!
!-- for each emitting species appearing in the chemical mechanism
!-- determine the volume source variable name ( 'vsrc_' + species name )
!-- and grab the emission values for the current time stamp
    DO  k = 1, num_emis_species

       CALL chem_emis_nc_get_vsrc_species_var_name( species_var_name, emis_species(k)%name_str )
       CALL get_variable( ncid, species_var_name, buf_emissions, time_index, buf_nvsrc,  .TRUE.  )

       DO  i = 1, num_vsrc_pos
          vsrc_emis_value(k,i) = buf_emissions( vsrc_nc_index(i) )  ! transfer variable         
       END DO

    END DO

!
!-- cleaning up

    DEALLOCATE ( buf_emissions )
    CALL close_input_file( ncid )

 END SUBROUTINE emissions_vsrc_update_lod2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  PRIVATE METHODS FOR TRAFFIC EMISSION PHYSICS (LOD 0/1)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  END OF MODULE
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 END MODULE chem_emis_traffic_mod
