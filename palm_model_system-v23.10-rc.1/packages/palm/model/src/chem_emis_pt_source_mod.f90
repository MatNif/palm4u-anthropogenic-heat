!> @file chem_emis_pt_source_mod.f90
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
! Copyright 2018-2022 Leibniz Universitaet Hannover
! Copyright 2018-2022 Karlsruhe Institute of Technology
! Copyright 2018-2022 Freie Universitaet Berlin
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
!> Module for emission mode for E-PRTR or GRETA point sources
!>
!> Wrappers to various LOD functionalities have prefix chem_emis_pt_source
!--------------------------------------------------------------------------------------------------!
 MODULE chem_emis_pt_source_mod

    USE chem_emis_generic_mod  ! for derived type chem_emis_species

    USE chem_emis_vsrc_mod     ! for derived type chem_emis_vsrc_pos

    USE chem_modules,                                                                              &          
        ONLY:  nc_field_length

    USE kinds

    IMPLICIT NONE
    SAVE 
    PRIVATE

!
!-- module-specific constants
    INTEGER(iwp), PARAMETER ::  default_k_spread = 2    !< maximum distribution
    INTEGER(iwp), PARAMETER ::  max_k_spread     = 5    !< maximum distribution
    INTEGER(iwp), PARAMETER ::  no_index         = -1   !< no assigned index
    INTEGER(iwp), PARAMETER ::  num_dim_ijk      = 3    !< size of 3D ijk dimensions
    INTEGER(iwp), PARAMETER ::  year_to_days     = 365  !< number of days in a year

    REAL(wp), PARAMETER ::  day_to_seconds = 86400.0    !< seconds in a day
!
!-- interface to KPP chemical mechanism
    INTEGER(iwp) ::  num_emis_species  !< emission species

    TYPE(chem_emis_species), ALLOCATABLE, DIMENSION(:) ::  emis_species  !< emission species
!
!-- interface to chem_emis_vsrc_mod
    INTEGER(iwp) ::  num_vsrc_pos  !< volume sources

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::  vsrc_emis_value  !< volume source volues

    TYPE(chem_emis_vsrc_pos), ALLOCATABLE, DIMENSION(:) ::  vsrc_pos  !< volume source positions
!
!-- consolidated namelist items
    INTEGER(iwp) ::  leap_year_day  !< leap year correction
    INTEGER(iwp) ::  num_k_weights  !< number of vertical distribution cells

    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  k_weights  !< normalized vertical distribution weights


    INTERFACE chem_emis_pt_source_init
       MODULE PROCEDURE chem_emis_pt_source_init
    END INTERFACE chem_emis_pt_source_init

    INTERFACE chem_emis_pt_source_cleanup
       MODULE PROCEDURE chem_emis_pt_source_cleanup
    END INTERFACE chem_emis_pt_source_cleanup

    INTERFACE chem_emis_pt_source_update
       MODULE PROCEDURE chem_emis_pt_source_update
    END INTERFACE chem_emis_pt_source_update

!
!-- public methods specific to generic emission mode (i.e., this module)
    PUBLIC ::  chem_emis_pt_source_cleanup, chem_emis_pt_source_init,                              &
               chem_emis_pt_source_update

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> wrapper for intitializing emission mode
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_pt_source_init( )

    IMPLICIT NONE

!
!-- wrapper for LOD 0 (anticipating new LODs)
    CALL emissions_pt_source_init_lod0( )

 END SUBROUTINE chem_emis_pt_source_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> extracts all volume source positions in subdomain from emissions netCDF file
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_pt_source_cleanup( )

    IMPLICIT NONE


    IF ( ALLOCATED( emis_species )    )  DEALLOCATE( emis_species )
    IF ( ALLOCATED( vsrc_pos )        )  DEALLOCATE( vsrc_pos )
    IF ( ALLOCATED( vsrc_emis_value ) )  DEALLOCATE( vsrc_emis_value )
    IF ( ALLOCATED( k_weights )       )  DEALLOCATE( k_weights )

 END SUBROUTINE chem_emis_pt_source_cleanup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates emissions at every time step (called in time_integratin)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_pt_source_update ( )

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_update_source

    IMPLICIT NONE

    INTEGER(iwp) ::  k  !< generic counter


    IF ( ( num_emis_species == 0 )  .OR.  ( num_vsrc_pos == 0 ) )  RETURN

    DO  k = 1, num_emis_species
       CALL chem_emis_vsrc_update_source( emis_species(k)%mech_index, vsrc_pos,                    &
                                          vsrc_emis_value(k,:) )
    ENDDO

 END SUBROUTINE chem_emis_pt_source_update


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> initialize emision mode module under LOD 0 (only mode)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emissions_pt_source_init_lod0( )

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_init_volume_source,                                                       &
               chem_emis_species_match  

    USE chem_modules,                                                                              &
        ONLY:  emis_pt_source_annual_values,                                                       &
               emis_pt_source_k_spread,                                                            &
               emis_pt_source_k_weights,                                                           &
               emis_pt_source_leap_year,                                                           &
               emis_pt_source_locations_ijk,                                                       &
               emis_pt_source_species_names 

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nyn,                                                                                &
               nxr,                                                                                &
               nys

    IMPLICIT NONE

    INTEGER(iwp) ::  i                !< generic index
    INTEGER(iwp) ::  j                !< generic index
    INTEGER(iwp) ::  k                !< generic index
    INTEGER(iwp) ::  k_user_pt        !< indices to map to volume source and namelist
    INTEGER(iwp) ::  k_user_species   !< indices to map to volume source and namelist
    INTEGER(iwp) ::  k_vsrc           !< indices to map to volume source and namelist

    INTEGER(iwp)                            ::  num_pt_sources    !< number of point sources defined in namelist
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  vsrc_user_indices !< user indices from namelist

    REAL(wp)                                ::  sum_k_weights     !< sum of spread weights (normalization factor)
    REAL(wp)                                ::  year_per_second   !< 1 / [seconds in a year]


    CALL location_message( 'reading LOD 0 point source emissions data', 'start' )

!
!-- get per second rate for converting annual output.
!-- assume 365 days per year unless leap year is indicated in namelist
    leap_year_day = 0
    IF ( emis_pt_source_leap_year )  leap_year_day = 1
    year_per_second = 1.0_wp / ( ( year_to_days + leap_year_day ) * day_to_seconds )
!
!-- get spread factors.  Note they are normalized so that
!-- sum of weights from 1 to k_spread = 1
    num_k_weights = emis_pt_source_k_spread
    IF ( ( num_k_weights < 1 )  .OR.  ( num_k_weights > max_k_spread ) )                           &
       num_k_weights = default_k_spread

    ALLOCATE( k_weights( num_k_weights ) )

    sum_k_weights = 0.0
    DO  k = 1, num_k_weights
       k_weights(k)  = emis_pt_source_k_weights(k)
       sum_k_weights = sum_k_weights + k_weights(k)      
    END DO

    IF ( .NOT. ( sum_k_weights > 0.0 ) )  sum_k_weights = 1.0_wp
    
    DO  k = 1, num_k_weights
        k_weights(k) = k_weights(k) / sum_k_weights
    END DO

!
!-- link volume source species to mechanism and user input
    CALL chem_emis_species_match( num_emis_species, emis_species, emis_pt_source_species_names )
!
!-- assign point source locations in domain
    ALLOCATE( vsrc_user_indices ( SIZE( emis_pt_source_locations_ijk,1 ) ) )
    vsrc_user_indices = no_index

    num_pt_sources = 0

    DO  k = 1, SIZE( emis_pt_source_locations_ijk,1 )
       IF ( ANY( emis_pt_source_locations_ijk(k,:) < 0 ) )  EXIT  ! assume no further source locations

       i = emis_pt_source_locations_ijk(k,1)
       j = emis_pt_source_locations_ijk(k,2)
!
!--    If cell is out of domain, proceed to next point.
       IF ( ( ( i < nxl )  .OR.  ( i > nxr ) )  .OR.  ( ( j < nys )  .OR.  ( j > nyn ) ) )  CYCLE

       num_pt_sources = num_pt_sources + 1
       vsrc_user_indices(num_pt_sources) = k
    ENDDO

!
!-- quit if there are no point sources
    num_vsrc_pos = num_k_weights * num_pt_sources

    IF ( num_vsrc_pos > 0 )  THEN
!
!--    iniatiate volume source and values
!--    it is assumed that k << top vertical layer though some kind of
!--    check would be nice as a prophylactic measure
       ALLOCATE( vsrc_pos(num_vsrc_pos) )
       ALLOCATE( vsrc_emis_value(num_emis_species,num_vsrc_pos) )

       vsrc_emis_value = 0.0_wp
       k_vsrc = 0

       DO  k = 1, num_pt_sources
          k_user_pt = vsrc_user_indices(k)                   ! point source index on namelist
          DO  j = 1, num_k_weights
             k_vsrc = k_vsrc + 1                             ! volume source index
             DO  i = 1, num_emis_species
                k_user_species     = emis_species(i)%user_index  ! species index on namelist
                vsrc_pos(k_vsrc)%i = emis_pt_source_locations_ijk(k_user_pt,1)
                vsrc_pos(k_vsrc)%j = emis_pt_source_locations_ijk(k_user_pt,2)
                vsrc_pos(k_vsrc)%k = emis_pt_source_locations_ijk(k_user_pt,3) + ( j - 1 )
                vsrc_emis_value(i,k_vsrc) = k_weights(j) * year_per_second *                       &
                               MAX( 0.0_wp, emis_pt_source_annual_values(k_user_pt,k_user_species) )
             ENDDO
          ENDDO
       ENDDO
!
!--    initialize volume source data and global volume source data structure
       CALL chem_emis_init_volume_source( vsrc_pos )

    ENDIF

    DEALLOCATE( vsrc_user_indices )

    CALL location_message( 'reading LOD 0 point source emissions data' , 'finished' )

 END SUBROUTINE emissions_pt_source_init_lod0

 END MODULE chem_emis_pt_source_mod
