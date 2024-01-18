!> @file chem_emis_generic_mod.f90
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
! @author Edward C. Chan
! @author Ilona Jäkel
!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Modules for generic emission mode for user-defined emission data.
!> Also contains interface functions and subroutines for other emission mode modules.
!>
!> Public methods specific to generic emission mode, wrappers to various LOD functionalities have
!> prefix chem_emis_generic.
!>
!> Public methods applicable or all emission modes for reading netCDF files have prefix chem_emis_nc
!>
!> Public methods apllicable for all emission modes have prefix chem_emissions.
!>
!> NOTE: usually a switch structure will govern how module should be initialized based on LOD
!>       but for generic mode only LOD 2 is supported.
!--------------------------------------------------------------------------------------------------!
 MODULE chem_emis_generic_mod

    USE chem_emis_vsrc_mod  ! for derived type chem_emis_vsrc_pos

    USE chem_modules,                                                                              &     
        ONLY:  nc_field_length
    
    USE control_parameters,                                                                        &
        ONLY:  coupling_char

    USE kinds

    IMPLICIT NONE
    SAVE 
    PRIVATE

!
!-- data type for saving species and their linkage to KPP chemical mechanism
    TYPE, PUBLIC ::  chem_emis_species
       CHARACTER(LEN=nc_field_length) ::  name_str     !< emission species
       INTEGER(iwp)                   ::  mech_index   !< index per KPP mechanism
       INTEGER(iwp)                   ::  user_index   !< index per user input from namelist/nc-file
    END TYPE chem_emis_species
!
!-- general characteristics
    CHARACTER(LEN=*), PARAMETER ::  input_nc_file = 'PIDS_EMIS_GENERIC'  !< emission mode input
    CHARACTER(LEN=*), PARAMETER ::  str_empty     = 'novalue'            !< default empty string
    INTEGER(iwp)                ::  lod                                  !< nc lod (always 2)
!
!-- timestamps
    CHARACTER(LEN=nc_field_length)                            ::  current_time   !< simulation time(time)
    CHARACTER(LEN=nc_field_length), ALLOCATABLE, DIMENSION(:) ::  timestamp      !< individual timestamps

    INTEGER(iwp) ::  num_timestamp     !< # timestamps
    INTEGER(iwp) ::  time_index        !< current timestamp
!
!-- interface to KPP chemical mechanism
    INTEGER(iwp)                                       ::  num_emis_species      !< # emission species
    TYPE(chem_emis_species), ALLOCATABLE, DIMENSION(:) ::  emis_species          !< emission species
!
!-- interface to chem_emis_vsrc_mod
    INTEGER(iwp)                                          ::  num_vsrc_pos       !< # volume sources
    INTEGER(iwp),             ALLOCATABLE, DIMENSION(:)   ::  vsrc_nc_index      !< location on nc file

    REAL(KIND=dp),            ALLOCATABLE, DIMENSION(:,:) ::  vsrc_emis_value    !< volume source volues
    TYPE(chem_emis_vsrc_pos), ALLOCATABLE, DIMENSION(:)   ::  vsrc_pos           !< volume source positions

!
!-- interface to mode-specific methods
    INTERFACE chem_emis_generic_init
      MODULE PROCEDURE chem_emis_generic_init
    END INTERFACE chem_emis_generic_init

    INTERFACE chem_emis_generic_cleanup
      MODULE PROCEDURE chem_emis_generic_cleanup
    END INTERFACE chem_emis_generic_cleanup

    INTERFACE chem_emis_generic_update
      MODULE PROCEDURE chem_emis_generic_update
    END INTERFACE chem_emis_generic_update
!
!-- interface to netCDF methods
    INTERFACE chem_emis_nc_get_time
      MODULE PROCEDURE chem_emis_nc_get_time
    END INTERFACE chem_emis_nc_get_time

    INTERFACE chem_emis_nc_get_species
      MODULE PROCEDURE chem_emis_nc_get_species
    END INTERFACE chem_emis_nc_get_species

    INTERFACE chem_emis_nc_get_vsrc_positions
      MODULE PROCEDURE chem_emis_nc_get_vsrc_positions
    END INTERFACE chem_emis_nc_get_vsrc_positions

    INTERFACE chem_emis_nc_get_vsrc_species_var_name
      MODULE PROCEDURE chem_emis_nc_get_vsrc_specieS_var_name
    END INTERFACE chem_emis_nc_get_vsrc_species_var_name
!
!-- interface to auxilliary methods
    INTERFACE chem_emis_get_current_time
      MODULE PROCEDURE chem_emis_get_current_time
    END INTERFACE chem_emis_get_current_time

    INTERFACE chem_emis_init_volume_source
      MODULE PROCEDURE chem_emis_init_volume_source
    END INTERFACE chem_emis_init_volume_source

    INTERFACE chem_emis_species_match
      MODULE PROCEDURE chem_emis_species_match
    END INTERFACE chem_emis_species_match

    INTERFACE chem_emis_update_trigger
      MODULE PROCEDURE chem_emis_update_trigger
    END INTERFACE chem_emis_update_trigger

!
!-- public methods specific to generic emission mode (i.e., this module)
    PUBLIC ::  chem_emis_generic_cleanup, chem_emis_generic_init,                                  &
               chem_emis_generic_update
!
!-- public methods for all emisison modes for reading netCDF files
    PUBLIC ::  chem_emis_nc_get_species, chem_emis_nc_get_time,                                    &
               chem_emis_nc_get_vsrc_positions,                                                    &
               chem_emis_nc_get_vsrc_species_var_name
!
!-- public methods for all emission modes
    PUBLIC ::  chem_emis_get_current_time,                                                         &
               chem_emis_init_volume_source,                                                       &
               chem_emis_species_match,                                                            &
               chem_emis_update_trigger

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> wrapper for intitializing emission mode
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_generic_init( )

    IMPLICIT NONE

!
!-- This is where a select case based on lod would have to be called if
!-- the emission mode supports multiple LODs
    CALL emissions_generic_init_lod2( )  ! initialize emission module
    CALL emissions_vsrc_update_lod2( )   ! update volume sources for first iteration

 END SUBROUTINE chem_emis_generic_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> extracts all volume source positions in subdomain from emissions netCDF file
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_generic_cleanup( )

    USE chem_emis_vsrc_mod

    IMPLICIT NONE

    IF ( ALLOCATED( timestamp )       )  DEALLOCATE( timestamp       )
    IF ( ALLOCATED( emis_species )    )  DEALLOCATE( emis_species    )
    IF ( ALLOCATED( vsrc_pos )        )  DEALLOCATE( vsrc_pos        )
    IF ( ALLOCATED( vsrc_emis_value ) )  DEALLOCATE( vsrc_emis_value )
    IF ( ALLOCATED( vsrc_nc_index )   )  DEALLOCATE( vsrc_nc_index   )

 END SUBROUTINE chem_emis_generic_cleanup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates emissions at every tiem step (called in time_integratin)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_generic_update( )

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_update_source

    USE chem_modules,                                                                              &
        ONLY:  nc_field_length
 
    IMPLICIT NONE

    CHARACTER(LEN=nc_field_length) ::  current_timestamp    !< simulation time(time)

    INTEGER ::  current_time_index   !< time index for this 
    INTEGER ::  k                    !< counter

!
!-- checks if timestamp has advanced ...
    current_time_index = chem_emis_get_current_time( current_timestamp, timestamp )

!
!-- ... and update source terms if they are different
!-- (this is where a select case based on lod would have to be called if
!-- the emission mode supports multiple LODs)
    IF ( current_time_index > time_index )  THEN   ! we are not going back in time, n'est pas?   
       time_index = current_time_index             ! update time index
       CALL emissions_vsrc_update_lod2( )          ! update emission volume source
    ENDIF

!
!-- then update global volume source 
!-- Don't forget this gets reset at every time step so even IF there is no change
!-- the current values still need be updated / accumulated
    DO  k = 1, num_emis_species
       CALL chem_emis_vsrc_update_source( emis_species(k)%mech_index, vsrc_pos, vsrc_emis_value(k,:) )
    ENDDO

 END SUBROUTINE


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> reads emisison netCDF file to extract timestamps
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_nc_get_time( ntimestamp, timestamp, ncid )

    USE chem_modules,                                                                              &
        ONLY:  nc_dim_ntime,                                                                       &
               nc_field_length,                                                                    &
               nc_var_timestamp

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  get_dimension_length,                                                               &
               get_variable

    USE palm_date_time_mod,                                                                        &
        ONLY:  date_time_str_len

    IMPLICIT NONE
!
!-- arguments in order of appearance
    INTEGER(iwp)                                              ::  ntimestamp   !< number of timestamps
    CHARACTER(LEN=nc_field_length), ALLOCATABLE, DIMENSION(:) ::  timestamp    !< timestamps
    INTEGER(iwp),                   INTENT(IN)                ::  ncid         !< nc file handle

    CHARACTER(LEN=512), ALLOCATABLE, DIMENSION(:) ::  buf_timestamp   !< buffer variable
    INTEGER(iwp)                                  ::  k               !< counter

!
!-- grab timestamps from netCDF file into a buffer
!-- note:  buffer will be allocated in get_variable ( )
    CALL get_dimension_length( ncid, ntimestamp, nc_dim_ntime )  
    CALL get_variable( ncid, nc_var_timestamp, buf_timestamp, ntimestamp, date_time_str_len )

!
!-- reassign timestamps to clear NULL character in C-strings
    ALLOCATE( timestamp( ntimestamp ) ) 
    DO  k = 1, num_timestamp
       CALL copy_chars_trimmable( timestamp(k), buf_timestamp(k), nc_field_length )
    ENDDO

!
!-- cleaning up
    DEALLOCATE( buf_timestamp )  !< allocated in get_variable

 END SUBROUTINE chem_emis_nc_get_time


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> extracts all species names from emissions netCDF file
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_nc_get_species( nspecies, emis_species, ncid )

    USE chem_modules,                                                                              &
        ONLY:  nc_dim_nspecies,                                                                    &
               nc_field_length,                                                                    &
               nc_var_species

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  get_dimension_length,                                                               &
               get_variable
 
    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp)                                        ::  nspecies       !< number of emitting species
    TYPE(chem_emis_species),  ALLOCATABLE, DIMENSION(:) ::  emis_species   !< emitting species from netCDF
    INTEGER(iwp),             INTENT(IN)                ::  ncid           !< nc file handle

    CHARACTER(LEN=512),       ALLOCATABLE, DIMENSION(:) ::  buf_species    !< buffer (512 to match shape)
    INTEGER(iwp)                                        ::  buf_nspecies   !< buffer species count

!
!-- get species names from netCDF file and store them in a buffer
!-- note that buf_species will be allocated in get_variable
    CALL get_dimension_length( ncid, buf_nspecies, nc_dim_nspecies )
    CALL get_variable( ncid, nc_var_species, buf_species, buf_nspecies, nc_field_length )

!
!-- start species match
    CALL chem_emis_species_match( nspecies, emis_species, buf_species )

!
!-- cleaning up
    DEALLOCATE( buf_species )  ! allocated in get_variable  

 END SUBROUTINE chem_emis_nc_get_species


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> extracts all volume source positions in subdomain from emissions netCDF file
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_nc_get_vsrc_positions( nvsrc, vsrc_pos, nc_indices, ncid )

    USE chem_emis_vsrc_mod

    USE indices,                                                                                   &   
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               topo_top_ind

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  get_dimension_length,                                                               &
               get_variable
 
    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp)                                        ::  nvsrc       !< number of volumes source positions
    TYPE(chem_emis_vsrc_pos), ALLOCATABLE, DIMENSION(:) ::  vsrc_pos    !< volume source position vector
    INTEGER(iwp),             ALLOCATABLE, DIMENSION(:) ::  nc_indices  !< matching indices for netCDF file
    INTEGER(iwp),             INTENT(IN)                ::  ncid        !< nc file handle

    INTEGER(iwp) ::  buf_nvsrc             !< number of volume source positions
    INTEGER(iwp) ::  i                     !< counter
    INTEGER(iwp) ::  k                     !< counter

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  buf_indices           !< matching indices for netCDF file
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  buf_i                 !< volume source index from netCDF file
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  buf_j                 !< volume source index from netCDF file
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  buf_k                 !< volume source index from netCDF file

!
!-- get volume source vector data from netCDF file
    CALL get_dimension_length( ncid, buf_nvsrc, nc_dim_nvsrc )

    ALLOCATE( buf_i(buf_nvsrc) )
    ALLOCATE( buf_j(buf_nvsrc) )
    ALLOCATE( buf_k(buf_nvsrc) )

    CALL get_variable( ncid, nc_var_vsrc_i, buf_i )
    CALL get_variable( ncid, nc_var_vsrc_j, buf_j )
    CALL get_variable( ncid, nc_var_vsrc_k, buf_k )

!
!-- initialize and increment for every position located in subdomain
!-- assuming domain decomposition is in ij directions only
    ALLOCATE( buf_indices(buf_nvsrc) )
    buf_indices = 0

    nvsrc = 0  

    DO  k = 1, buf_nvsrc
       IF ( ( ( buf_i(k) >= nxl )  .AND.  ( buf_i(k) <= nxr) )  .AND.                                  &
            ( ( buf_j(k) >= nys )  .AND.  ( buf_j(k) <= nyn) ) )  THEN
          nvsrc = nvsrc + 1
          buf_indices(nvsrc) = k
       ENDIF
    ENDDO

!
!-- allocate space for return variables and assign
    ALLOCATE( vsrc_pos(nvsrc) )
    ALLOCATE( nc_indices(nvsrc) )

    DO  k = 1, nvsrc
       nc_indices(k) = buf_indices(k)
       i = nc_indices(k)
       vsrc_pos(k)%i = buf_i(i)
       vsrc_pos(k)%j = buf_j(i)
       vsrc_pos(k)%k = buf_k(i) + topo_top_ind(buf_j(i),buf_i(i),0)
    ENDDO

!
!-- clean up all buffers 
    DEALLOCATE( buf_indices )
    DEALLOCATE( buf_i )
    DEALLOCATE( buf_j )
    DEALLOCATE( buf_k )

 END SUBROUTINE chem_emis_nc_get_vsrc_positions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> given species name, return corresponding volume source variable name
!> (effectively prefix + species name)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_nc_get_vsrc_species_var_name( var_name, species_name )

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  nc_var_vsrc_prefix
    
    USE chem_modules,                                                                              &
        ONLY:  nc_field_length
 
    IMPLICIT NONE

!
!-- arguments in order of appearance
    CHARACTER(LEN=nc_field_length)             ::  var_name      !< volume source variable name
    CHARACTER(LEN=nc_field_length), INTENT(IN) ::  species_name  !< species name

    var_name = TRIM( nc_var_vsrc_prefix ) // TRIM( species_name )

 END SUBROUTINE chem_emis_nc_get_vsrc_species_var_name
 

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> gets current time stamp and closest non-exceeding index
!--------------------------------------------------------------------------------------------------!
 FUNCTION chem_emis_get_current_time( current_timestamp, timestamp )  RESULT ( k )

    USE chem_modules,                                                                              &    
        ONLY:  nc_field_length

    USE control_parameters,                                                                        &    
        ONLY:  time_since_reference_point

    USE palm_date_time_mod,                                                                        &    
        ONLY:  date_time_str_len,                                                                  &
               get_date_time

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp)                                             ::  k                   !< timestamp index
    CHARACTER(LEN=nc_field_length)                           ::  current_timestamp   !< current timestamp
    CHARACTER(LEN=nc_field_length), INTENT(IN), DIMENSION(:) ::  timestamp           !< listing of timestamps

    CALL get_date_time( time_since_reference_point, date_time_str=current_timestamp )

    DO  k = ( date_time_str_len+1 ), nc_field_length   ! get_date_time ( ) is returning dirty strings
       current_timestamp(k:k) = ' '                    ! so this di loop is necessary to clean it out
    ENDDO                                              ! so that timestamps can be compared

    k = locate_timestamp( timestamp, current_timestamp )

 END FUNCTION chem_emis_get_current_time


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> identifies emission species that appear in the KPP mechanism and user input from
!> namelist and netCDF files, for instance.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_species_match( nspecies, emis_species, ref_species )

    USE chem_emis_vsrc_mod,                                                                        &   
        ONLY:  default_index

    USE chem_gasphase_mod,                                                                         &   
        ONLY:  nvar,                                                                               &
               spc_names

    USE chem_modules,                                                                              &   
        ONLY:  nc_field_length

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp)                                       ::  nspecies      !< number of emitted species in mechanism
    TYPE(chem_emis_species), ALLOCATABLE, DIMENSION(:) ::  emis_species  !< emitted species in mechanism
    CHARACTER(LEN=*),        INTENT(IN),  DIMENSION(:) ::  ref_species   !< all emitted species listed from input
!
!-- local variables
    CHARACTER(LEN=nc_field_length)          ::  this_species      !< buffer for current species name 

    INTEGER(iwp)                            ::  i                 !< counter
    INTEGER(iwp)                            ::  k                 !< counter    
    INTEGER(iwp)                            ::  num_ref_species   !< total number of emitted species from input

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  mech_indices      !< indices of input species in melchanism
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  user_indices      !< indices of input species in user input

!
!-- initialize species indices and count
    nspecies        = 0
    num_ref_species = SIZE(ref_species)

    ALLOCATE( mech_indices(num_ref_species) )
    ALLOCATE( user_indices(num_ref_species) )
    mech_indices = default_index
    user_indices = default_index

!
!-- brute force search (no other way, spec_names is not sorted)
    DO  k = 1, num_ref_species

       IF ( TRIM( ref_species(k) ) == str_empty )  EXIT  ! assume rest of array is empty; exit loop
       CALL copy_chars_trimmable( this_species, ref_species(k), nc_field_length )

       DO  i = 1, nvar
          IF ( TRIM( this_species ) == TRIM( spc_names(i) ) )  THEN
             nspecies = nspecies + 1
             mech_indices(nspecies) = i
             user_indices(nspecies) = k
          ENDIF
       ENDDO

    ENDDO

!
!-- allocate space for emission species and assign values
    ALLOCATE( emis_species(nspecies) )

    DO  k = 1, nspecies
       emis_species(k)%mech_index = mech_indices(k)
       emis_species(k)%user_index = user_indices(k)
       CALL copy_chars_trimmable( emis_species(k)%name_str, ref_species(user_indices(k)),          &
                                  nc_field_length )
    ENDDO
       
    DEALLOCATE( mech_indices )
    DEALLOCATE( user_indices )
 
 END SUBROUTINE chem_emis_species_match


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> determine if update interval is reached
!
!> The time since last update will monotonically increase until it reaches the update interval
!> point, which it will decrease once before continuing to increase towards the next update
!> interval, thus signalling the update.  This approach guarantees that the update takes place
!> once and only once at each update interval irrespective of floating point precision.
!--------------------------------------------------------------------------------------------------!
 FUNCTION chem_emis_update_trigger( time_since_last_update, dt_update )  RESULT ( update_trigger )

    USE control_parameters,                                                                        &
        ONLY:  time_since_reference_point

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    IMPLICIT NONE

    LOGICAL ::  update_trigger  !< whether update is to be triggered

    REAL(WP)             ::  time_since_last_update  !< time since last update at previous time step
    REAL(wp), INTENT(IN) ::  dt_update               !< user-defined update interval

!
!-- Local variables.
    REAL(wp) ::  current_time_since_last_update  !< time since last update at current time step


!
!-- Get current time and time since last update.
    CALL get_date_time( time_since_reference_point )
    current_time_since_last_update = MOD( time_since_reference_point, dt_update )

!
!-- Always assume no update unless the update interval has been crossed.
    update_trigger = .FALSE.
    IF ( .NOT. ( current_time_since_last_update > time_since_last_update ) )  THEN
       update_trigger = .TRUE.
    ENDIF

!
!-- Update time from last update.
    time_since_last_update = current_time_since_last_update

 END FUNCTION chem_emis_update_trigger


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> initialize volume source hashkey and global volume source listing
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_init_volume_source( vsrc_pos )

    USE chem_emis_vsrc_mod

    IMPLICIT NONE

!
!-- arguments in order of appearance
    TYPE(chem_emis_vsrc_pos), ALLOCATABLE, DIMENSION(:) ::  vsrc_pos  !< volume source position vector
!
!-- local variables
    INTEGER(iwp) ::  k  !< counter

    IF ( SIZE( vsrc_pos ) == 0 )  RETURN

!
!-- calculate hash key value for each volume source location
    DO  k = 1, SIZE(vsrc_pos)
       vsrc_pos(k)%key = chem_emis_vsrc_key_from_indices( vsrc_pos(k)%i, vsrc_pos(k)%j,            &
                                                          vsrc_pos(k)%k )
    ENDDO

!
!-- append unique volume source to global volume source listing
    CALL chem_emis_vsrc_append_positions( vsrc_pos )

 END SUBROUTINE chem_emis_init_volume_source


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> initialize emission mode module under LOD 2 (only mode)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emissions_generic_init_lod2( )

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
    ALLOCATE( vsrc_emis_value( num_emis_species, num_vsrc_pos ) )
    vsrc_emis_value = 0.0

    CALL close_input_file( ncid )  ! closes netCDF file
    CALL location_message( 'reading emissions data from ' // TRIM(input_nc_file), 'finished' )

 END SUBROUTINE emissions_generic_init_lod2


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates emission volume source under LOD 2 (only mode)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emissions_vsrc_update_lod2( )

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
     
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) ::  buf_emissions     !< emission values


!
!-- setting everything up
    CALL open_read_file( TRIM( input_nc_file ) // TRIM( coupling_char ), ncid )
    CALL get_dimension_length( ncid, buf_nvsrc, nc_dim_nvsrc )  ! gets dimension
    ALLOCATE( buf_emissions(buf_nvsrc) )

!
!-- for each emitting species appearing in the chemical mechanism
!-- determine the volume source variable name ( 'vsrc_' + species name )
!-- and grab the emission values for the current time stamp
    DO  k = 1, num_emis_species

       CALL chem_emis_nc_get_vsrc_species_var_name( species_var_name, emis_species(k)%name_str )
       CALL get_variable( ncid, species_var_name, buf_emissions, time_index, buf_nvsrc,  .TRUE.  )

       DO  i = 1, num_vsrc_pos
          vsrc_emis_value(k,i) = buf_emissions( vsrc_nc_index(i) )  ! transfer variable
       ENDDO

    ENDDO

!
!-- cleaning up
    DEALLOCATE( buf_emissions )
    CALL close_input_file( ncid )

 END SUBROUTINE emissions_vsrc_update_lod2


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> given target timestamp locate most recent timestep in timestap array
!> timestamp is assumed to be sorted
!--------------------------------------------------------------------------------------------------!
 RECURSIVE FUNCTION locate_timestamp( x, xk, i0, i1 )  RESULT( k )

    USE chem_modules,                                                                              &
        ONLY:  nc_field_length

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp)                                             ::  k      !< index
    CHARACTER(LEN=nc_field_length), INTENT(IN), DIMENSION(:) ::  x      !< array of timestamps
    CHARACTER(LEN=nc_field_length), INTENT(IN)               ::  xk     !< target timestamp
    INTEGER(iwp),                   INTENT(IN), OPTIONAL     ::  i0     !< start index
    INTEGER(iwp),                   INTENT(IN), OPTIONAL     ::  i1     !< end index

    INTEGER(iwp) ::  k0  !< search bound
    INTEGER(iwp) ::  k1  !< search bound

!
!-- assign search bounds
    k0 = 1
    k1 = SIZE(x)
    IF ( PRESENT(i0) )  k0 = i0
    IF ( PRESENT(i1) )  k1 = i1

!
!-- Termination conditions.
    k = k0
!
!-- Return closest previous element.
    IF ( ( k1 - k0 ) <= 1 )  RETURN
    k = k1
!
!-- Special case to look for last element.
    IF ( TRIM( xk ) >= TRIM( x(k1) ) )  RETURN
    k = ( k0 + k1 ) / 2
!
!-- Or if timestamp is right on.
    IF ( TRIM( xk ) == TRIM( x(k)  ) )  RETURN

!
!-- Otherwise restrict search bounds and continue.
    IF ( TRIM( xk ) > TRIM( x(k) ) )  THEN
       k = locate_timestamp( x, xk, k, k1 )
    ELSE
       k = locate_timestamp( x, xk, k0, k  )
    ENDIF

 END FUNCTION locate_timestamp


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> copy string in a way that is trimmable (i.e., replace trailing null characters to white space)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE copy_chars_trimmable( des, src, len_des )

    IMPLICIT NONE

!
!-- arguments in order of appearance
    CHARACTER(LEN=*)             ::  des      !< destination string
    CHARACTER(LEN=*), INTENT(IN) ::  src      !< source string
    INTEGER(iwp),     INTENT(IN) ::  len_des  !< number of characters

    INTEGER(iwp)            ::  k             !< counter
    INTEGER(iwp), PARAMETER ::  NULL  =  0    !< null char in C string (i.e., from netCDF char variables)
    INTEGER(iwp), PARAMETER ::  SPACE = 32    !< "white space" - can be trimmed in Fortran


    DO  k = 1, len_des
      des(k:k) = CHAR( SPACE )  !< clear NULL character (in case of incoming C string)
      IF ( src(k:k) /= CHAR( NULL ) )  des(k:k) = src(k:k)
    ENDDO

 END SUBROUTINE copy_chars_trimmable

 END MODULE chem_emis_generic_mod
