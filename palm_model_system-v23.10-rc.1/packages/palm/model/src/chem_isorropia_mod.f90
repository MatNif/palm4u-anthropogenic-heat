!> @file chem_isorropia_mod.f90
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
!> Module for chemical production of secondary inorganic aerosols using the thermodynamical
!> equilibrium model ISORROPIA (Nenes, 1998) and ISORROPIA II (Fountoukis and Nenes, 2007).
!>
!> ISORPROPIA (II) subroutines will be imported from corresponding third-party shared libraries.
!> Also contains interface functions and subroutines for other emission mode modules.
!> Wrappers for these subroutines and default parameters for all variables used can be found in
!> the module chem_isorropia_interface_mod.
!--------------------------------------------------------------------------------------------------!
 MODULE chem_isorropia_mod

    USE chem_emis_generic_mod

    USE chem_isorropia_interface_mod 

    USE chem_modules,                                                                              &     
        ONLY:  nc_field_length

    USE kinds

    IMPLICIT NONE

    SAVE 

    PRIVATE

!
!-- Update control.
    REAL(wp)            ::  time_since_last_update              !< tracker for update interval

    REAL(wp), PARAMETER ::  default_update_interval = 300.0_wp  !< time interval between updates
!
!-- ISORROPIA parameters.
    INTEGER(iwp) ::  activity_coefficient_method  !< iacalci (0,1)
    INTEGER(iwp) ::  mass_conservation_mode       !< nadji (ISORROPIA II only)
    INTEGER(iwp) ::  max_activity_sweep           !< nsweepi
    INTEGER(iwp) ::  max_iteration                !< maxiti
    INTEGER(iwp) ::  mdr_weight_method            !< wftypi (0,1,2)
    INTEGER(iwp) ::  root_subdivisions            !< ndivi

    REAL(wp) ::  activity_tolerance  !< epsacti
    REAL(wp) ::  aerosol_state       !< cntrl[2] (0,1)
    REAL(wp) ::  problem_type        !< cntrl[1] (0,1)
    REAL(wp) ::  solver_tolerance    !< epsI
!
!-- Chemical species linkage to mechanism.
    INTEGER(iwp) ::  num_aqueous_aerosol_species  !< number of aqueous aerosols (aerliq) in mechanism
    INTEGER(iwp) ::  num_gas_species              !< number of gas species (gas) in mechanism
    INTEGER(iwp) ::  num_solid_aerosol_species    !< number of solid aerosols (aersld) in mechanism

    TYPE(chem_emis_species), ALLOCATABLE, DIMENSION(:) ::  aqueous_aerosol_species  !< aerliq
    TYPE(chem_emis_species), ALLOCATABLE, DIMENSION(:) ::  gas_species              !< gas
    TYPE(chem_emis_species), ALLOCATABLE, DIMENSION(:) ::  solid_aerosol_species    !< aersld

    INTERFACE chem_isorropia_init
      MODULE PROCEDURE chem_isorropia_init
    END INTERFACE chem_isorropia_init

    INTERFACE chem_isorropia_cleanup
      MODULE PROCEDURE chem_isorropia_cleanup
    END INTERFACE chem_isorropia_cleanup

    INTERFACE chem_isorropia_update
      MODULE PROCEDURE chem_isorropia_update
    END INTERFACE chem_isorropia_update

    PUBLIC ::  chem_isorropia_init, chem_isorropia_cleanup, chem_isorropia_update

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Wrapper for intitializing emission mode.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_isorropia_init( )

    USE chem_modules,                                                                              &
        ONLY:  chem_isorropia_activity_coefficient_method,                                         &
               chem_isorropia_activity_tolerance,                                                  &
               chem_isorropia_aerosol_state,                                                       &
               chem_isorropia_mass_conservation_mode,                                              &
               chem_isorropia_max_activity_sweep,                                                  &
               chem_isorropia_max_iteration,                                                       &
               chem_isorropia_mdr_weight_method,                                                   &
               chem_isorropia_problem_type,                                                        &
               chem_isorropia_root_subdivisions,                                                   &
               chem_isorropia_solver_tolerance,                                                    &
               chem_isorropia_update_interval

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_species_match

    IMPLICIT NONE


    CALL location_message( 'initializing ISORROPIA '//isr_versions_id, 'start' )

!
!-- Set update interval.
    IF ( chem_isorropia_update_interval < 0.0 )  THEN
       chem_isorropia_update_interval = default_update_interval
    ENDIF

!
!-- Problem setup (CNTRL).
    IF ( chem_isorropia_problem_type  < 0.0 )  chem_isorropia_problem_type  = isr_cntrl_problem_default
    problem_type  = chem_isorropia_problem_type
    IF ( chem_isorropia_aerosol_state < 0.0 )  chem_isorropia_aerosol_state = isr_cntrl_state_default
    aerosol_state = chem_isorropia_aerosol_state
    
!
!-- Solver parameters (SETPARM).
    mdr_weight_method           = chem_isorropia_mdr_weight_method
    activity_coefficient_method = chem_isorropia_activity_coefficient_method
    solver_tolerance            = chem_isorropia_solver_tolerance
    max_iteration               = chem_isorropia_max_iteration
    max_activity_sweep          = chem_isorropia_max_activity_sweep
    activity_tolerance          = chem_isorropia_activity_tolerance
    root_subdivisions           = chem_isorropia_root_subdivisions
    mass_conservation_mode      = chem_isorropia_mass_conservation_mode

!
!-- Pass parameters to ISORROPIA library.
    CALL isr_setparm( mdr_weight_method, activity_coefficient_method, solver_tolerance,            &
                      max_iteration, max_activity_sweep, activity_tolerance, root_subdivisions,    &
                      mass_conservation_mode )
!
!-- Where invalid entries will be replaced with default values.
    CALL isr_getparm( mdr_weight_method, activity_coefficient_method, solver_tolerance,            &
                      max_iteration, max_activity_sweep, activity_tolerance, root_subdivisions,    &
                      mass_conservation_mode )

!
!-- Repopulate entries back to namelist items (for chem_header).
    chem_isorropia_mdr_weight_method           = mdr_weight_method
    chem_isorropia_activity_coefficient_method = activity_coefficient_method
    chem_isorropia_solver_tolerance            = solver_tolerance
    chem_isorropia_max_iteration               = max_iteration
    chem_isorropia_max_activity_sweep          = max_activity_sweep
    chem_isorropia_activity_tolerance          = activity_tolerance
    chem_isorropia_root_subdivisions           = root_subdivisions
    chem_isorropia_mass_conservation_mode      = mass_conservation_mode
    

!
!-- Tag species in mechanism.
    CALL chem_emis_species_match( num_gas_species, gas_species, isr_gas_species    )
    CALL chem_emis_species_match( num_aqueous_aerosol_species, aqueous_aerosol_species,            &
                                  isr_aerliq_species )
    CALL chem_emis_species_match( num_solid_aerosol_species, solid_aerosol_species,                &
                                  isr_aersld_species )

    CALL location_message( 'initializing ISORROPIA '//isr_versions_id, 'finished' )

 END SUBROUTINE chem_isorropia_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Wrapper for deallocating all module-specific dynamic memory.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_isorropia_cleanup( )

    IMPLICIT NONE


    IF ( ALLOCATED( aqueous_aerosol_species ) )  DEALLOCATE( aqueous_aerosol_species )
    IF ( ALLOCATED( gas_species )             )  DEALLOCATE( gas_species )
    IF ( ALLOCATED( solid_aerosol_species )   )  DEALLOCATE( solid_aerosol_species )

 END SUBROUTINE chem_isorropia_cleanup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calls ISORROPIA off-line to obtain ion concentrations.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_isorropia_update( )

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nzb,                                                                                &
               nzt

    USE arrays_3d,                                                                                 &
        ONLY:  exner,                                                                              &
               hyp,                                                                                &
               pt,                                                                                 &
               q,                                                                                  &
               rho_air
              
    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  get_relative_humidity

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_update_trigger

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nspec

    USE chem_modules,                                                                              &
        ONLY:  chem_species,                                                                       &
               chem_isorropia_update_interval

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    IMPLICIT NONE

    INTEGER(iwp) ::  i          !< generic index
    INTEGER(iwp) ::  j          !< generic index
    INTEGER(iwp) ::  k          !< generic index
    INTEGER(iwp) ::  k_mech     !< mech index
    INTEGER(iwp) ::  k_species  !< species index counter
    INTEGER(iwp) ::  k_user     !< user (isorropia) index
!
!-- Derived thermodynamic state.
    REAL(KIND=dp) ::  relative_humidity  !< RH
    REAL(KIND=dp) ::  temperature        !< cell temperature
   
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) ::  gasphase_species  !< gasphase speceies concentrations
!
!-- ISORROPIA subroutine arguments.
    CHARACTER(LEN=isr_scasi_len) ::  case_identifier  !<

    REAL(KIND=dp), DIMENSION(isr_num_aerliq) ::  aqueous_aerosol_concentrations  !< aerliq
    REAL(KIND=dp), DIMENSION(isr_num_other)  ::  diagnostic_data                 !< other
    REAL(KIND=dp), DIMENSION(isr_num_gas)    ::  gas_concentrations              !< gas
    REAL(KIND=dp), DIMENSION(isr_num_w)      ::  precursors                      !< wi
    REAL(KIND=dp), DIMENSION(isr_num_aersld) ::  solid_aerosol_concentrations    !< aersld
    REAL(KIND=dp), DIMENSION(isr_num_w)      ::  total_concentrations            !< wt

    
!
!-- Update only when update interval is reached.
    IF ( .NOT. chem_emis_update_trigger( time_since_last_update, chem_isorropia_update_interval ) )&
    THEN
       RETURN
    ENDIF

    CALL cpu_log( log_point_s(56), 'chem.isorropia.update', 'start' )

    ALLOCATE( gasphase_species( nspec ) )

!
!-- Loop through all cells in domain.
    DO  j = nys, nyn
       DO  i = nxl, nxr
          DO  k = nzb, nzt
!
!--          Skip to next cell if no relevant trace or precursor species is found.
             IF ( .NOT. collect_species( gasphase_species, k, j, i ) )  CYCLE
             IF ( .NOT. collect_precursors( precursors, gasphase_species ) )  CYCLE

             temperature       = pt(k,j,i) * exner(k)
             relative_humidity = get_relative_humidity( q(k,j,i), temperature, hyp(k), rho_air(k) )

!
!--          Convert concentration unit from ppm to mol/m3 as required by ISORRPIA.
             DO  k_species = 1, isr_num_w
                precursors(k_species) = ppm_to_molm3( precursors(k_species), temperature, hyp(k) )
             ENDDO

             total_concentrations           = 0.0_dp
             gas_concentrations             = 0.0_dp
             aqueous_aerosol_concentrations = 0.0_dp
             solid_aerosol_concentrations   = 0.0_dp
             diagnostic_data                = 0.0_dp

!
!--          Call to ISORROPIA solver.
             CALL isr_isoropia( precursors, relative_humidity, temperature,                        &
                                (/ problem_type, aerosol_state /), total_concentrations,           &
                                gas_concentrations, aqueous_aerosol_concentrations,                &
                                solid_aerosol_concentrations, case_identifier, diagnostic_data )

!
!--          Revert concentration unit from mol/m3 to ppm for chemistry solver.
             DO  k_species = 1, num_aqueous_aerosol_species
                k_mech = aqueous_aerosol_species(k_species)%mech_index
                k_user = aqueous_aerosol_species(k_species)%user_index
                chem_species(k_mech)%conc(k,j,i) =                                                 &
                        molm3_to_ppm( aqueous_aerosol_concentrations(k_user), temperature, hyp(k) )
             ENDDO


             DO  k_species = 1, num_solid_aerosol_species
                k_mech = solid_aerosol_species(k_species)%mech_index
                k_user = solid_aerosol_species(k_species)%user_index
                chem_species(k_mech)%conc(k,j,i) =                                                 &
                        molm3_to_ppm( solid_aerosol_concentrations(k_user), temperature, hyp(k) )
             ENDDO

          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE( gasphase_species )

    CALL cpu_log( log_point_s(56), 'chem.isorropia.update', 'stop' )

 END SUBROUTINE chem_isorropia_update


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Assign gas phase species.
!> Returns FALSE if sum of all gas phase species concentration is zero.
!--------------------------------------------------------------------------------------------------!
 FUNCTION collect_species( concentrations, k, j, i )  RESULT ( has_species )

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nspec

    USE chem_modules,                                                                              &
        ONLY:  chem_species

    IMPLICIT NONE

!
!-- Arguments in order of appearance
    LOGICAL                     ::  has_species     !< whether cell contains gas phase species; output
    REAL(KIND=dp), DIMENSION(:) ::  concentrations  !< gas phase species concentrations
    INTEGER(iwp)                ::  k               !< cell index
    INTEGER(iwp)                ::  j               !< cell index
    INTEGER(iwp)                ::  i               !< cell index
!
!-- Local variables.
    INTEGER(iwp) ::  k_species  !< gas phase species index


!
!-- Reset gas psecies concentrations.
    concentrations = 0.0_dp

    DO  k_species = 1, nspec
        concentrations(k_species) = chem_species(k_species)%conc(k,j,i)
    ENDDO    

!
!-- Initally assume no gas phase species.
    has_species = .FALSE.
    IF ( SUM( concentrations ) > 0.0 )  has_species = .TRUE.

 END FUNCTION collect_species


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Assign species data to precursors.
!> Returns FALSE if sum of all precursor concentration is zero.
!--------------------------------------------------------------------------------------------------!
 FUNCTION collect_precursors( precursors, concentrations )  RESULT ( has_precursor )

    IMPLICIT NONE
!
!-- Arguents in order of appearance.
    LOGICAL                     ::  has_precursor   !< whether cell contains precursors; output
    REAL(KIND=dp), DIMENSION(:) ::  precursors      !< precursor concentrations
    REAL(KIND=dp), DIMENSION(:) ::  concentrations  !< gas phase species concentrations
!
!-- Local variables.
    INTEGER(iwp) ::  i        !< generic counter
    INTEGER(iwp) ::  k        !< generic counter
    INTEGER(iwp) ::  k_mech   !< species mechanism index (i.e. KPP)
    INTEGER(iwp) ::  k_mole   !< index for ISORROPIA molar concentrations


!
!-- Reset precursor concentrations.
    precursors = 0.0_dp
!
!-- Accumulate precursors.
    DO  k = 1, isr_num_w
!
!--    Gas phase species:
       DO  i = 1, num_gas_species
           k_mech = gas_species(i)%mech_index
           k_mole = isr_total_num_w * ( gas_species(i)%user_index - 1 ) + k
           precursors(k) = precursors(k) + ( concentrations(k_mech) *                              &
                                             isr_precursor_gas_molar_concentrations(k_mole) )
       ENDDO
!
!--    Aqueous aerosol species:
       DO  i = 1, num_aqueous_aerosol_species
           k_mech = aqueous_aerosol_species(i)%mech_index
           k_mole = isr_total_num_w * ( aqueous_aerosol_species(i)%user_index - 1 ) + k
           precursors(k) = precursors(k) + ( concentrations(k_mech) *                              &
                                             isr_precursor_aerliq_molar_concentrations(k_mole) )
       ENDDO
!
!--    Solid aersol species:
       DO  i = 1, num_solid_aerosol_species
           k_mech = solid_aerosol_species(i)%mech_index
           k_mole = isr_total_num_w * ( solid_aerosol_species(i)%user_index - 1 ) + k
           precursors(k) = precursors(k) + ( concentrations(k_mech) *                              &
                                             isr_precursor_aersld_molar_concentrations(k_mole) )
       ENDDO

    ENDDO
!
!-- Always assume no precursor.
    has_precursor = .FALSE.
    IF ( SUM( precursors ) > 0.0 )  has_precursor = .TRUE.

 END FUNCTION


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Converts concentration unit from ppm to mol/m3.
!--------------------------------------------------------------------------------------------------!
 FUNCTION ppm_to_molm3( ppm, temperature, pressure )  RESULT ( molm3 )

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  rgas_univ
    
    IMPLICIT NONE
!
!-- Arguments in order of appearance.
    REAL(KIND=dp) ::  molm3        !< concentration in mol/m3
    REAL(KIND=dp) ::  ppm          !< concentration in ppm
    REAL(KIND=dp) ::  temperature  !< thermodynamic state
    REAL(KIND=dp) ::  pressure  


    molm3 = 1.0E-6 * ppm * pressure / ( rgas_univ * temperature )

 END FUNCTION ppm_to_molm3


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Converts concentration unit from mol/m3 to ppm.
!--------------------------------------------------------------------------------------------------!
 FUNCTION molm3_to_ppm( molm3, temperature, pressure )  RESULT ( ppm )

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  rgas_univ

    IMPLICIT NONE
!
!-- Arguments in order of appearance.
    REAL(KIND=dp) ::  ppm          !< concentration in ppm
    REAL(KIND=dp) ::  molm3        !< concentration in mol/m3
    REAL(KIND=dp) ::  temperature  !< thermodynamic state
    REAL(KIND=dp) ::  pressure  


    ppm = 1.0E6 * molm3 * ( rgas_univ * temperature ) / pressure

 END FUNCTION molm3_to_ppm

 END MODULE chem_isorropia_mod
