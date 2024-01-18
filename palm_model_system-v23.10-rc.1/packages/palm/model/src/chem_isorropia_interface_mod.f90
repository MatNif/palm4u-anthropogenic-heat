!> @file chem_isorropia_interface_mod
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
!> Module for interfacing with third-party secondary inorganic aerosol library ISORROPIA
!> (Nenes, 1998) and ISORROPIA II (Fountoukis and Nenes, 2007).
!> Placeholder functions will be used in the event that no ISORROPIA library is supplied
!> during compilation.
!> All public symbols in the module will contain the prefix "isr_" to indicate namespace.
!--------------------------------------------------------------------------------------------------!
 MODULE chem_isorropia_interface_mod

    USE kinds

    USE chem_modules,                                                                              &
        ONLY:  nc_field_length

    IMPLICIT NONE

    SAVE
!-- @TODO: I don't understand the comment. The following statement will make every variable PUBLIC,
!--        so it does just the opposite of what the comment says
    PUBLIC  !-- only function propotypes and parameters are defined (no public variables)

!
!-- ISORROPIA version.
#if defined( __ISORROPIA1 )
    CHARACTER(LEN=*), PARAMETER ::  isr_versions_id  = 'I'
#elif defined( __ISORROPIA2 )
    CHARACTER(LEN=*), PARAMETER ::  isr_versions_id  = 'II'
#else
    CHARACTER(LEN=*), PARAMETER ::  isr_versions_id  = 'invalid'
#endif

!
!-- Absolute array sizes.
    INTEGER(iwp), PARAMETER ::  isr_errmsg_len        = 40  !< error message
    INTEGER(iwp), PARAMETER ::  isr_scasi_len         = 15  !< case identifier
    INTEGER(iwp), PARAMETER ::  isr_total_num_aerliq  = 15  !< liquid aerosol species
    INTEGER(iwp), PARAMETER ::  isr_total_num_aersld  = 19  !< solid aerosol species
    INTEGER(iwp), PARAMETER ::  isr_total_num_cntrl   =  2  !< program control
    INTEGER(iwp), PARAMETER ::  isr_total_num_errstki = 25  !< error stack size
    INTEGER(iwp), PARAMETER ::  isr_total_num_gas     =  3  !< gas phase species
    INTEGER(iwp), PARAMETER ::  isr_total_num_other   =  9  !< misc. quantities
    INTEGER(iwp), PARAMETER ::  isr_total_num_w       =  8  !< precursor (WI) and total (WT)
    INTEGER(iwp), PARAMETER ::  isr_version_len       = 14  !< version string
!
!-- Version-specific array sizes (default to ISORROPIA 1).
#if defined( __ISORROPIA2 )
    INTEGER(iwp), PARAMETER ::  isr_num_aerliq  = isr_total_num_aerliq
    INTEGER(iwp), PARAMETER ::  isr_num_aersld  = isr_total_num_aersld
    INTEGER(iwp), PARAMETER ::  isr_num_cntrl   = isr_total_num_cntrl 
    INTEGER(iwp), PARAMETER ::  isr_num_errstki = isr_total_num_errstki
    INTEGER(iwp), PARAMETER ::  isr_num_gas     = isr_total_num_gas
    INTEGER(iwp), PARAMETER ::  isr_num_other   = isr_total_num_other
    INTEGER(iwp), PARAMETER ::  isr_num_w       = isr_total_num_w
#else
    INTEGER(iwp), PARAMETER ::  isr_num_aerliq  = 12
    INTEGER(iwp), PARAMETER ::  isr_num_aersld  =  9
    INTEGER(iwp), PARAMETER ::  isr_num_cntrl   = isr_total_num_cntrl
    INTEGER(iwp), PARAMETER ::  isr_num_errstki = isr_total_num_errstki
    INTEGER(iwp), PARAMETER ::  isr_num_gas     = isr_total_num_gas
    INTEGER(iwp), PARAMETER ::  isr_num_other   =  7
    INTEGER(iwp), PARAMETER ::  isr_num_w       =  5
#endif

!
!-- Indices for array cntrl.
!-- ISORROPIA solver control.
    INTEGER(iwp), PARAMETER ::  isr_cntrl_problem = 1   !< problem type
    INTEGER(iwp), PARAMETER ::  isr_cntrl_state   = 2   !< aerosol state
!
!-- Possible values for cntrl (isoropia).
!-- ISORROPIA solver control.
    REAL(KIND=dp), PARAMETER ::  isr_cntrl_problem_forward    = 0.0_dp  !< forward problem (aerosol+gas)
    REAL(KIND=dp), PARAMETER ::  isr_cntrl_problem_reverse    = 1.0_dp  !< reverse problem (aerosol only)
    REAL(KIND=dp), PARAMETER ::  isr_cntrl_state_aqueous      = 1.0_dp  !< aqueous only aerosol state
    REAL(KIND=dp), PARAMETER ::  isr_cntrl_state_deliquescent = 0.0_dp  !< deliquescent aerosol state
!
!-- Possible values for wftypi (setparm / getparm).
!-- MDR weighting method.
    INTEGER(iwp), PARAMETER ::  isr_wftypi_dry  = 0  !< mutual deliquescence regions are assumed dry
    INTEGER(iwp), PARAMETER ::  isr_wftypi_half = 1  !< solution assumed half dry / half wet
    INTEGER(iwp), PARAMETER ::  isr_wftypi_mean = 2  !< RH weighted mean of dry / wet solutions
!
!-- Default parameters for isoropia (see subroutine definitions for details).
    REAL(KIND=dp), PARAMETER ::  isr_cntrl_problem_default = isr_cntrl_problem_forward
    REAL(KIND=dp), PARAMETER ::  isr_cntrl_state_default   = isr_cntrl_state_deliquescent
!
!-- Possible values for iacalci (setparm / getparm).
!-- Activity coefficient calculation method.
    INTEGER(iwp), PARAMETER ::  isr_iacalci_lookup  = 1  !< use precalculated tables
    INTEGER(iwp), PARAMETER ::  isr_iacalci_runtime = 0  !< calculate activity coefficients during runtime
!
!-- Possible values for nadji (setparm / getparm) for ISORROPIA II.
!-- Mass conservation mode.
    INTEGER(iwp), PARAMETER ::  isr_nadji_adjust = 1  !< automatic adjustment for mass conservation
    INTEGER(iwp), PARAMETER ::  isr_nadji_normal = 0  !< normal calculation
!
!-- Default parameters for setparm (see subroutine definitions for details).
    INTEGER(iwp),  PARAMETER ::  isr_iacalci_default = isr_iacalci_lookup
    INTEGER(iwp),  PARAMETER ::  isr_maxiti_default  = 100
    INTEGER(iwp),  PARAMETER ::  isr_nadji_default   = isr_nadji_normal    ! ISORROPIA II only
    INTEGER(iwp),  PARAMETER ::  isr_ndivi_default   = 5
    INTEGER(iwp),  PARAMETER ::  isr_nsweepi_default = 4
    INTEGER(iwp),  PARAMETER ::  isr_wftypi_default  = isr_wftypi_mean
    
    REAL(KIND=dp), PARAMETER ::  isr_epsacti_default = 0.05_dp
    REAL(KIND=dp), PARAMETER ::  isr_epsi_default    = 1.0E-6

!
!-- Descriptiosn of chemical species (w, gas, aerliq, aersld), must be of length nc_field_length
!-- array element alignment.
    CHARACTER(LEN=nc_field_length), PARAMETER, DIMENSION(isr_total_num_w) ::                       &
       isr_w_species      = (/ 'NA                                                              ', &
                               'SO4                                                             ', &
                               'NH4                                                             ', &
                               'NO3                                                             ', &
                               'CL                                                              ', &
                               'CA                                                              ', &
                               'K                                                               ', &
                               'MG                                                              ' /)
                                !---5----1----5----2----5----3----5----4----5----5----5----6---!

    CHARACTER(LEN=nc_field_length), PARAMETER, DIMENSION(isr_total_num_gas) ::                     &
       isr_gas_species    = (/ 'NH3                                                             ', &
                               'HNO3                                                            ', &
                               'HCL                                                             ' /)
                                !---5----1----5----2----5----3----5----4----5----5----5----6---!

    CHARACTER(LEN=nc_field_length), PARAMETER, DIMENSION(isr_total_num_aerliq) ::                  &
       isr_aerliq_species = (/ 'H_aq                                                            ', &
                               'NA_aq                                                           ', &
                               'NH4_aq                                                          ', &
                               'CL_aq                                                           ', &
                               'SO4_aq                                                          ', &
                               'HSO4_aq                                                         ', &
                               'NO3_aq                                                          ', &
                               'H2O_aq                                                          ', &
                               'NH3_aq                                                          ', &
                               'HCL_aq                                                          ', &
                               'HNO3_aq                                                         ', &
                               'OH_aq                                                           ', &
                               'CA_aq                                                           ', &
                               'K_aq                                                            ', &
                               'MG_aq                                                           ' /)
                                !---5----1----5----2----5----3----5----4----5----5----5----6---!

    CHARACTER(LEN=nc_field_length), PARAMETER, DIMENSION(isr_total_num_aersld) ::                  &
       isr_aersld_species = (/ 'NANO3_s                                                         ', &
                               'NH4NO3_s                                                        ', &
                               'NACL_s                                                          ', &
                               'NH4CL_s                                                         ', &
                               'NA2SO4_s                                                        ', &
                               'NH42SO4_s                                                       ', &
                               'NAHSO4_s                                                        ', &
                               'NH4HSO4_s                                                       ', &
                               'NH44HSO42_s                                                     ', &
                               'CASO4_s                                                         ', &
                               'CANO32_s                                                        ', &
                               'CACL2_s                                                         ', &
                               'K2SO4_s                                                         ', &
                               'KHSO4_s                                                         ', &
                               'KNO3_s                                                          ', &
                               'KCL_s                                                           ', &
                               'MGSO4_s                                                         ', &
                               'MGNO32_s                                                        ', &
                               'MGCL2_s                                                         ' /)
                                !---5----1----5----2----5----3----5----4----5----5----5----6---!

!
!-- Descriptions of diagnostic variables (cntrl, other).
    CHARACTER(LEN=nc_field_length), PARAMETER, DIMENSION(isr_total_num_cntrl) ::                   &
       isr_cntrl_names    = (/ 'problem type                                                    ', &
                               'aerosol state                                                   ' /)
                                !---5----1----5----2----5----3----5----4----5----5----5----6---!
                                                               
    CHARACTER(LEN=nc_field_length), PARAMETER, DIMENSION(isr_total_num_other) ::                   &
       isr_other_names    = (/ 'aqueous phase indicator                                         ', &
                               'total sulfate molar ratio                                       ', &
                               'aerosol sulfate molar ratio                                     ', &
                               'total sodium molar ratio                                        ', &
                               'ionic strength of aqueous aerosol                               ', &
                               'total calls to activity coefficient calculation subroutine      ', &
                               'total sulphate molar ratio with crustal species                 ', &
                               'crustals to sodium molar ratio                                  ', &
                               'crustal species molar ratio                                     ' /)
                                !---5----1----5----2----5----3----5----4----5----5----5----6---!

!
!-- Molar precursor concentrations (viz. Ch 6 ISORROPIA v 2.1 reference manual).
    REAL(KIND=dp), PARAMETER, DIMENSION(isr_total_num_gas*isr_total_num_w) ::                      &
       isr_precursor_gas_molar_concentrations =                                                    &
                            !  Na   SO4  NH4  NO3  Cl   Ca   K    Mg
                            (/ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! NH3
                               0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,           &  ! HNO3
                               0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 /)            ! HCl
                            !  Na   SO4  NH4  NO3  Cl   Ca   K    Mg

    REAL(KIND=dp), PARAMETER, DIMENSION(isr_total_num_aerliq*isr_total_num_w) ::                   &
       isr_precursor_aerliq_molar_concentrations =                                                 &
                            !  Na   SO4  NH4  NO3  Cl   Ca   K    Mg
                            (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! H
                               1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! Na
                               0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! NH4
                               0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,           &  ! Cl
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! SO4
                               0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! HSO4
                               0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,           &  ! NO3
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! H2O
                               0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! NH3
                               0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,           &  ! HCl
                               0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,           &  ! HNO3
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! OH
                               0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,           &  ! Ca
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,           &  ! K
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 /)            ! Mg
                            !  Na   SO4  NH4  NO3  Cl   Ca   K    Mg

    REAL(KIND=dp), PARAMETER, DIMENSION(isr_total_num_aersld*isr_total_num_w) ::                   &
       isr_precursor_aersld_molar_concentrations =                                                 &
                            !  Na   SO4  NH4  NO3  Cl   Ca   K    Mg
                            (/ 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,           &  ! NaNO3
                               0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,           &  ! NH4NO3
                               1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,           &  ! NaCl
                               0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,           &  ! NH4Cl
                               2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! Na2SO4
                               0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! (NH4)2SO4
                               1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! NaHSO4
                               0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! (NH4)HSO4
                               0.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0,           &  ! (NH4)4H(SO4)2
                               0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,           &  ! CaSO4
                               0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0,           &  ! Ca(NO3)2
                               0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0,           &  ! CaCl2
                               0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0,           &  ! K2SO4
                               0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,           &  ! KHSO4
                               0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,           &  ! KNO3
                               0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,           &  ! KCl
                               0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,           &  ! MgSO4
                               0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0,           &  ! Mg(NO3)2
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0 /)            ! MgCl2
                            !  Na   SO4  NH4  NO3  Cl   Ca   K    Mg

    INTERFACE isr_isoropia             ! ISORROPIA driver (note spelling)
       MODULE PROCEDURE isr_isoropia
    END INTERFACE isr_isoropia

    INTERFACE isr_setparm              ! set execution parameters
       MODULE PROCEDURE isr_setparm    ! can be specified in namelist
    END INTERFACE isr_setparm

    INTERFACE isr_getparm              ! obtain execution parameters
       MODULE PROCEDURE isr_getparm
    END INTERFACE isr_getparm

    INTERFACE isr_iserrinf             ! diagnostic functions (print error stack)
        MODULE PROCEDURE isr_iserrinf
    END INTERFACE isr_iserrinf

    INTERFACE isr_errstat              ! print error message
       MODULE PROCEDURE isr_errstat
    END INTERFACE isr_errstat

    INTERFACE isr_isorinf              ! obtain general information on ISORROPIA library
       MODULE PROCEDURE isr_isorinf
    END INTERFACE isr_isorinf

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> ISORROPIA driver
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE isr_isoropia( wi, rhi, tempi, cntrl, wt, gas, aerliq, aersld, scasi, other )

    IMPLICIT NONE

!
!-- Arguments in order of appearance.
    REAL(KIND=dp), DIMENSION(isr_num_w)      ::  wi     !< precursor conentrations
    REAL(KIND=dp)                            ::  rhi    !< relative humidity (0-1)
    REAL(KIND=dp)                            ::  tempi  !< ambient temperature [K]
    REAL(KIND=dp), DIMENSION(isr_num_cntrl)  ::  cntrl  !< problem control
    REAL(KIND=dp), DIMENSION(isr_num_w)      ::  wt     !< total species concentrations
    REAL(KIND=dp), DIMENSION(isr_num_gas)    ::  gas    !< gaeous species concentrations
    REAL(KIND=dp), DIMENSION(isr_num_aerliq) ::  aerliq !< liquid aerosol concentrations
    REAL(KIND=dp), DIMENSION(isr_num_aersld) ::  aersld !< solid aerosol concentrations
    CHARACTER(LEN=isr_scasi_len)             ::  scasi  !< case identifier
    REAL(KIND=dp), DIMENSION(isr_num_other)  ::  other  !< other output quantities
 

#if defined( __ISORROPIA1 ) || defined( __ISORROPIA2 )
!
!-- Check run control.
    IF ( ( cntrl(isr_cntrl_problem) < isr_cntrl_problem_forward )  .OR.                            &
         ( cntrl(isr_cntrl_problem) > isr_cntrl_problem_reverse ) )                                &
    THEN
       cntrl(isr_cntrl_problem) = isr_cntrl_problem_default
    ENDIF

    IF ( ( cntrl(isr_cntrl_state) < isr_cntrl_state_deliquescent )  .OR.                           &
         ( cntrl(isr_cntrl_state) > isr_cntrl_state_aqueous ) )                                    &
    THEN
       cntrl(isr_cntrl_state) = isr_cntrl_state_default
    ENDIF

    CALL ISOROPIA( wi, rhi, tempi, cntrl, wt, gas, aerliq, aersld, scasi, other )

#else
!
!-- Dummy call (does nothing).
    wi     = wi
    rhi    = rhi
    tempi  = tempi
    cntrl  = cntrl
    wt     = wt
    gas    = gas
    aerliq = aerliq
    aersld = aersld
    scasi  = scasi
    other  = other

#endif    

 END SUBROUTINE isr_isoropia


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set ISORROPIA solver PARAMETERs.
!> Note difference in arguments between ISORROPIA 1 and 2.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE isr_setparm( wftypi, iacalci, epsi, maxiti, nsweepi, epsacti, ndivi, nadji )

    IMPLICIT NONE

!
!-- Arguments in order of appearance.
    INTEGER(iwp)  ::  wftypi  !< mutual deliquescence region weighting
    INTEGER(iwp)  ::  iacalci !< activity coefficeint calculation method
    REAL(KIND=dp) ::  epsi    !< general convergence criterion
    INTEGER(iwp)  ::  maxiti  !< general maximum iterations
    INTEGER(iwp)  ::  nsweepi !< maximum iterations for activity coefficient calculation
    REAL(KIND=dp) ::  epsacti !< convergence criterion for activity coefficient calculation
    INTEGER(iwp)  ::  ndivi   !< subdivisions for root tracking
    INTEGER(iwp)  ::  nadji   !< enforcement of mass conservation (ISORROPIA II only)


#if defined( __ISORROPIA1 ) || defined( __ISORROPIA2 )
!
!-- Input check.
    IF ( ( wftypi < isr_wftypi_dry )  .OR.  ( wftypi > isr_wftypi_mean ) )  THEN
       wftypi = isr_wftypi_default
    ENDIF

    IF ( ( iacalci < isr_iacalci_runtime )  .OR.  ( iacalci > isr_iacalci_lookup ) )  THEN
       iacalci = isr_iacalci_default
    ENDIF

    IF ( epsi    < 0.0 )  epsi    = isr_epsi_default
    IF ( maxiti  <=  0 )  maxiti  = isr_maxiti_default
    IF ( nsweepi <=  0 )  nsweepi = isr_nsweepi_default
    IF ( epsacti < 0.0 )  epsacti = isr_epsacti_default
    IF ( ndivi   <=  0 )  ndivi   = isr_ndivi_default

    IF ( ( nadji < isr_nadji_normal )  .OR.  ( nadji > isr_nadji_adjust ) )  THEN
       nadji = isr_nadji_default
    ENDIF

!
!-- Select subroutine based on ISORROPIA version
#if defined( __ISORROPIA2 )
    CALL SETPARM( wftypi, iacalci, epsi, maxiti, nsweepi, epsacti, ndivi, nadji )
#else
    CALL SETPARM( wftypi, iacalci, epsi, maxiti, nsweepi, epsacti, ndivi )
#endif

#else
!
!-- Dummy call (does nothing).
    wftypi  = wftypi
    iacalci = iacalci
    epsi    = epsi
    maxiti  = maxiti
    nsweepi = nsweepi
    epsacti = epsacti
    ndivi   = ndivi
    nadji   = nadji
#endif

 END SUBROUTINE isr_setparm


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set ISORROPIA solver PARAMETERs.
!> Note difference in arguments between ISORROPIA 1 and 2.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE isr_getparm( wftypi, iacalci, epsi, maxiti, nsweepi, epsacti, ndivi, nadji )

    IMPLICIT NONE

!
!-- Arguments in order of appearance.
    INTEGER(iwp)  ::  wftypi   !< mutual deliquescence region weighting
    INTEGER(iwp)  ::  iacalci  !< activity coefficeint calculation method
    REAL(KIND=dp) ::  epsi     !< general convergence criterion
    INTEGER(iwp)  ::  maxiti   !< general maximum iterations
    INTEGER(iwp)  ::  nsweepi  !< maximum iterations for activity coefficient calculation
    REAL(KIND=dp) ::  epsacti  !< convergence criterion for activity coefficient calculation
    INTEGER(iwp)  ::  ndivi    !< subdivisions for root tracking
    INTEGER(iwp)  ::  nadji    !< enforcement of mass conservation (ISORROPIA II only)


#if defined( __ISORROPIA1 ) || defined( __ISORROPIA2 )
!
!-- Select subroutine based on ISORROPIA version.
#if defined( __ISORROPIA2 )
    CALL GETPARM( wftypi, iacalci, epsi, maxiti, nsweepi, epsacti, ndivi, nadji )
#else
    CALL GETPARM( wftypi, iacalci, epsi, maxiti, nsweepi, epsacti, ndivi )
    nadji = nadji
#endif

#else
!
!-- Dummy call (does nothing).
    wftypi  = wftypi
    iacalci = iacalci
    epsi    = epsi
    maxiti  = maxiti
    nsweepi = nsweepi
    epsacti = epsacti
    ndivi   = ndivi
    nadji   = nadji
#endif

 END SUBROUTINE isr_getparm


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Retrieves a copy of the error stack.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE isr_iserrinf( errstki, errmsgi, noferi, stkofl )

    IMPLICIT NONE

!
!-- Arguments in order of appearance.
    INTEGER(iwp),                  DIMENSION(isr_num_errstki) ::  errstki  !< error stack
    CHARACTER(LEN=isr_errmsg_len), DIMENSION(isr_num_errstki) ::  errmsgi  !< error message
    INTEGER(iwp)                                              ::  noferi   !< number of errors
    LOGICAL                                                   ::  stkofl   !< error stack overflow


#if defined( __ISORROPIA1 ) || defined( __ISORROPIA2 )
    CALL ISERRINF( errstki, errmsgi, noferi, stkofl )
#else
!
!-- Dummy call (does nothing).
    errstki = errstki
    errmsgi = errmsgi
    noferi  = noferi
    stkofl  = stkofl
#endif

 END SUBROUTINE isr_iserrinf


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Prints error message.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE isr_errstat( io, ierr, errinf )

    IMPLICIT NONE

!
!-- Arguments in order of appearance.
    INTEGER(iwp)     ::  io      !< unit of already opened file
    INTEGER(iwp)     ::  ierr    !< error code
    
    CHARACTER(LEN=*) ::  errinf  !< error message

    
#if defined( __ISORROPIA1 ) || defined( __ISORROPIA2 )
    CALL ERRSTAT( io, ierr, errinf )
#else
!
!-- Dummy call (does nothing).
    io     = io
    ierr   = ierr
    errinf = errinf
#endif

 END SUBROUTINE isr_errstat


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Returns basic information of ISORROPIA.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE isr_isorinf( versi, ncmp, nion, naqgas, nsol, nerr, tin, grt )

    IMPLICIT NONE

!
!-- Arguments in order of appearance.
    CHARACTER(LEN=isr_version_len) ::  versi   !< version string
    
    INTEGER(iwp)                   ::  ncmp    !< number of components in array WI / WT
    INTEGER(iwp)                   ::  nion    !< max number of aqueous ionic species
    INTEGER(iwp)                   ::  naqgas  !< max number of undissociated aqueuous species
    INTEGER(iwp)                   ::  nsol    !< max number of solid aerosol species
    INTEGER(iwp)                   ::  nerr    !< max size of error stack
    
    REAL(KIND=dp)                  ::  tin     !< value of numerical zero (tiny)
    REAL(KIND=dp)                  ::  grt     !< value of numerical infinity (great)


#if defined( __ISORROPIA1 ) || defined( __ISORROPIA2 )
    CALL ISORINF( versi, ncmp, nion, naqgas, nsol, nerr, tin, grt )
#else
!
!-- Dummy call (does nothing).
    versi  = versi
    ncmp   = ncmp
    nion   = nion
    naqgas = naqgas
    nsol   = nsol
    nerr   = nerr
    tin    = tin
    grt    = grt
#endif
  
 END SUBROUTINE isr_isorinf

 END MODULE chem_isorropia_interface_mod
