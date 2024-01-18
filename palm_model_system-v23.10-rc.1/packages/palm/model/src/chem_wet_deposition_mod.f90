!> @file chem_wet_deposition_mod.f90
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
!> Module for the EMEP wet deposition (Simpson et al 2012, DOI 10.5194/acp-12-7825-2012).
!--------------------------------------------------------------------------------------------------!
 MODULE chem_wet_deposition_mod

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  precipitation

    USE chem_emis_generic_mod

    USE chem_modules,                                                                              &          
        ONLY: nc_field_length

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nys,                                                                                &
               nyn,                                                                                &
               nzb,                                                                                &
               nzt

    USE kinds

    IMPLICIT NONE

    SAVE 

    PRIVATE

    INTEGER(iwp), PARAMETER :: max_species = 9  !< maximum number of species

    INTEGER(iwp) ::  layer_cloud_lower  !< lower cloud layer
    INTEGER(iwp) ::  layer_cloud_upper  !< upper cloud layer
    INTEGER(iwp) ::  num_species        !< number of species

    LOGICAL :: model_override  !< switch to override bulk cloud model
!
!-- Species specific parameters (Supplement Table S20 of Simpson et al,
!-- DOI:10.5194/acp-12-7825-2012)
    CHARACTER(LEN=nc_field_length), DIMENSION(max_species), PARAMETER ::  species_names =          &  !< species names
        (/ 'SO2 ', 'HNO3', 'HONO', 'NH3 ', 'H2O2', 'HCHO', 'ROOH', 'PM25', 'PM10' /)                  !< xtra space 2 match lengths

    REAL(wp), DIMENSION(max_species), PARAMETER ::  collection_efficiency =                        &  !< particle collection efficiency
        (/  0.0,   0.0,    0.0,    0.0,   0.0,    0.0,    0.0,    0.02,   0.4 /)
    REAL(wp), DIMENSION(max_species), PARAMETER ::  scavenging_ratio_in   =                        &  !< in-could scavenging ratios
        (/  0.3,   1.4,    1.4,    1.4,   1.4,    0.1,    0.05,   1.0,    1.0 /)
    REAL(wp), DIMENSION(max_species), PARAMETER ::  scavenging_ratio_sub  =                        &  !< below-could scavenging ratios
        (/  0.15,  0.5,    0.5,    0.5,   0.5,    0.03,   0.015,  0.0,    0.0 /)  
!
!-- Single-valued constants (section 9 of Simpson et al, DOI:10.5194/acp-12-7825-2012).
    REAL(wp), PARAMETER ::  collection_coefficient  = 5.2_wp     !< Marshall-Palmer collection coefficient [m3/kg/s]
    REAL(wp), PARAMETER ::  scavenging_depth        = 1000.0_wp  !< scavenging depth [m]
    REAL(wp), PARAMETER ::  scavenging_ratio_factor = 1.0E+6_wp  !< scaling factor for scavenging ratios []
    REAL(wp), PARAMETER ::  raindrop_fall_speed     = 5.0_wp     !< raindrop fall speed [m/s]

    REAL(wp) ::  ppr_rain_override           !< precipitation
    REAL(wp) ::  update_interval = 300.0_wp  !< time between updates [s]
    REAL(wp) ::  time_since_last_update      !< time since last update

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: scavenging_rate    !< vertical column scavenging rates

    TYPE(chem_emis_species), ALLOCATABLE, DIMENSION(:) ::  species_info  !< link between data and mechanism

    INTERFACE chem_wet_deposition_cleanup
       MODULE PROCEDURE chem_wet_deposition_cleanup
    END INTERFACE chem_wet_deposition_cleanup

    INTERFACE chem_wet_deposition_init
       MODULE PROCEDURE chem_wet_deposition_init
    END INTERFACE chem_wet_deposition_init

    INTERFACE chem_wet_deposition_update_ij
       MODULE PROCEDURE chem_wet_deposition_update_ij
    END INTERFACE chem_wet_deposition_update_ij

    PUBLIC ::  chem_wet_deposition_cleanup, chem_wet_deposition_init, chem_wet_deposition_update_ij

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_wet_deposition_init( )

    USE chem_modules,                                                                              &          
        ONLY: chem_wet_deposition_cloud_level_lower,                                               &
              chem_wet_deposition_cloud_level_upper,                                               &
              chem_wet_deposition_model_override,                                                  &
              chem_wet_deposition_rain_rate,                                                       &
              chem_wet_deposition_update_interval

    USE chem_emis_generic_mod,                                                                     &          
        ONLY: chem_emis_species_match

    IMPLICIT NONE


    CALL location_message( 'initializing wet deposition module', 'start' )

    model_override = chem_wet_deposition_model_override

    IF ( .NOT. ( precipitation  .OR.  model_override )  )  THEN
        CALL location_message( 'wet deposition disabled - no precipitation model found',           &
                               'finished' )
        RETURN
    ENDIF
!
!-- Use update interval given in namelist if it is valid (i.e, non-negative).
    IF ( chem_wet_deposition_update_interval >= 0.0_wp )  THEN
       update_interval = chem_wet_deposition_update_interval
    ENDIF
!
!-- Match species in mechanism.
    CALL chem_emis_species_match( num_species, species_info, species_names )

    IF ( num_species == 0 )  THEN
        CALL location_message( 'wet deposition disabled - no active species found in mechanism',   &
                               'finished' )
        RETURN
    ENDIF

    IF ( model_override )  THEN
        layer_cloud_lower = chem_wet_deposition_cloud_level_lower
        layer_cloud_upper = chem_wet_deposition_cloud_level_upper
        ppr_rain_override = chem_wet_deposition_rain_rate
    ENDIF
!
!-- Init scavenging rate storage.
    ALLOCATE( scavenging_rate(num_species,nzb:nzt+1,nys:nyn,nxl:nxr) )
    scavenging_rate = 0.0_wp

    CALL location_message( 'initializing wet deposition module', 'finished' )

 END SUBROUTINE chem_wet_deposition_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Deallocate storage.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_wet_deposition_cleanup( )

    IMPLICIT NONE


    IF ( ALLOCATED( scavenging_rate ) )  DEALLOCATE( scavenging_rate )

 END SUBROUTINE chem_wet_deposition_cleanup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates emissions at every time step (called in time_integratin).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_wet_deposition_update_ij( i, j )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  i  !< horizontal grid index
    INTEGER(iwp), INTENT(IN) ::  j  !< horizontal grid index

    INTEGER(iwp) ::  k0  !< lower index bound for cloud layer
    INTEGER(iwp) ::  k1  !< upper index bound for cloud layer


!
!-- No precipitation.
    IF ( .NOT. ( precipitation  .OR.  model_override ) )  RETURN
!
!-- Mechanism contains no relevant species.
    IF ( num_species == 0 )  RETURN
!
!-- Update scavenging rates if necessary.
    IF ( chem_emis_update_trigger( time_since_last_update, update_interval ) )  THEN
        scavenging_rate = 0.0_wp
        CALL find_cloud_layer_ij( i, j, k0, k1 ) 
        CALL update_scavenging_rates_ij( i, j, k0, k1 )
    ENDIF

    CALL update_species_ij( i, j )
    
 END SUBROUTINE chem_wet_deposition_update_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Finds indices for the lowest cloud layer given a horizontal position.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE find_cloud_layer_ij( i, j, k0, k1 )

    USE arrays_3d,                                                                                 &
        ONLY: ql

    INTEGER(iwp), INTENT(IN)  ::  i   !< horizontal grid index
    INTEGER(iwp), INTENT(IN)  ::  j   !< horizontal grid index
    INTEGER(iwp), INTENT(OUT) ::  k0  !< lower index bound for cloud layer
    INTEGER(iwp), INTENT(OUT) ::  k1  !< upper index bound for cloud layer

    REAL(wp) ::  ql0  !< cloud threshold


    IF ( model_override )  THEN
        k0 = layer_cloud_lower
        k1 = layer_cloud_upper
    ELSE    
        ql0 = 1.0E-5_wp  ! as per van Zanten et. al., 2011
        k0  = MINLOC( ql(:,j,i),          DIM=1, MASK=( ql(:,j,i) > ql0 ) )
        k1  = MINLOC( ql(k0+1:nzt+1,j,i), DIM=1, MASK=( ql(:,j,i) < ql0 ) )
    ENDIF

 END SUBROUTINE find_cloud_layer_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates scavenging and collection coefficients over all species for each vertical column.
!> All equations referenced from Simpson et al DOI 10.5194/acp-12-7825-201
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE update_scavenging_rates_ij( i, j, k0, k1 )

    USE arrays_3d,                                                                                 &
        ONLY: rho_air, prr_rain

    USE basic_constants_and_equations_mod,                                                         &
        ONLY: rho_l

    USE indices,                                                                                   &
        ONLY: topo_top_ind

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  i   !< horizontal grid index
    INTEGER(iwp), INTENT(IN) ::  j   !< horizontal grid index
    INTEGER(iwp), INTENT(IN) ::  k0  !< lower index bound for cloud layer
    INTEGER(iwp), INTENT(IN) ::  k1  !< upper index bound for cloud layer

    INTEGER(iwp) ::  k          !< generic vertical index
    INTEGER(iwp) ::  k_data     !< species index in data arrays (species_names, scavenging_ratios, etc.)
    INTEGER(iwp) ::  k_species  !< species index in species_info
    INTEGER(iwp) ::  k_sub      !< vertical index directly below cloud layer
    INTEGER(iwp) ::  k_topo     !< vertical index directly above terrain topology

    REAL(wp) ::  c_fact  !< dimensionless collection prefactor
    REAL(wp) ::  e_bar   !< aerosol colletion effeciency
    REAL(wp) ::  s_fact  !< dimensionless scavenging prefactor
    REAL(wp) ::  w_in    !< in-cloud scavenging ratio
    REAL(wp) ::  w_sub   !< below-cloud scavenging ratio

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rain_rate  !< placeholder for ppr_rain


    k_sub  = k0 - 1
    k_topo = topo_top_ind(j,i,0)

    s_fact = scavenging_ratio_factor / ( scavenging_depth * rho_l )
    c_fact = collection_coefficient / raindrop_fall_speed

    ALLOCATE( rain_rate(nzb:nzt+1) )

    IF ( model_override )  THEN
        rain_rate = ppr_rain_override
    ELSE
        rain_rate = prr_rain(:,j,i)
    ENDIF

    DO  k_species = 1, num_species

       k_data = species_info(k_species)%user_index
       w_in   = scavenging_ratio_in(k_data)
       w_sub  = scavenging_ratio_sub(k_data)
       e_bar  = collection_efficiency(k_data)
!
!--    Eq. [71].
       DO  k = k0, k1
          scavenging_rate(k_species,k,j,i) = rho_air(k) * rain_rate(k) * s_fact * w_in
       ENDDO
!
!--    Eq. [72] + [73]
       DO  k = k_topo, k_sub
          scavenging_rate(k_species,k,j,i) =  rho_air(k) * rain_rate(k) * ( ( s_fact * w_sub )     &
                                                                          + ( c_fact * e_bar ) )
       ENDDO

    ENDDO

 END SUBROUTINE update_scavenging_rates_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Settle species from top of cloud layer downwards.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE update_species_ij( i, j )

    USE chem_modules,                                                                              &
        ONLY: chem_species

    USE control_parameters,                                                                        &
        ONLY: dt_3d

    INTEGER(iwp), INTENT(IN) ::  i  !< horizontal grid index
    INTEGER(iwp), INTENT(IN) ::  j  !< horizontal grid index

    INTEGER(iwp) ::  k_mech     !< species index in chemistry_gasphase_mod
    INTEGER(iwp) ::  k_species  !< species index in species_info

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  column  !< working copy of columnar concentration
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  flux    !< wet deposition factor with range (0,1)

    ALLOCATE( column(nzb:nzt+1) )
    ALLOCATE( flux(nzb:nzt+1) )

    DO  k_species = 1, num_species
       k_mech = species_info(k_species)%mech_index
       column = chem_species(k_mech)%conc(:,j,i)
       flux = 1.0_wp - MIN( dt_3d * scavenging_rate(k_species,:,j,i), 1.0_wp )
       chem_species(k_mech)%conc(:,j,i) = column * flux
    ENDDO

    DEALLOCATE( column )
    DEALLOCATE( flux   )

 END SUBROUTINE update_species_ij

 END MODULE chem_wet_deposition_mod
