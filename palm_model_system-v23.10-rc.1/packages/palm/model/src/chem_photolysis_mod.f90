!> @file chem_photolysis_mod.f90
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
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Renate Forkel
!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> photolysis models and interfaces (Adapted from photolysis_model_mod.f90)
!> @todo more complex scheme, add shading
!--------------------------------------------------------------------------------------------------!
 MODULE chem_photolysis_mod

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nphot,                                                                              &
               phot,                                                                               &
               phot_names

    USE chem_modules,                                                                              &
        ONLY:  phot_frequen,                                                                       &
               photolysis_scheme,                                                                  &
               photolysis_shading

    USE control_parameters,                                                                        &
        ONLY:  time_since_reference_point

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nzb,                                                                                &
               nzt

    USE kinds

    IMPLICIT NONE


!   LOGICAL ::  unscheduled_photolysis_calls = .TRUE., & !< flag parameter indicating whether additional calls of the photolysis
!                                                        !< code are allowed
!               constant_albedo = .FALSE.,             & !< flag parameter indicating whether the albedo may change depending on
!                                                        !< zenith
!               force_photolysis_call = .FALSE.,       & !< flag parameter for unscheduled photolysis calls
!               photolysis = .FALSE.,                  & !< flag parameter indicating whether the photolysis model is used
!               sun_up    = .TRUE.,                    & !< flag parameter indicating whether the sun is up or down
!               photolysis = .TRUE.,                   & !< flag parameter indicing whether photolysis shall be calculated
!               sun_direction = .FALSE.                  !< flag parameter indicing whether solar direction shall be calculated

!
!-- Parameters for constant photolysis frequencies
    INTEGER,PARAMETER :: nconst = 15               !< available predefined photolysis prequencies for constant
!
!-- Names for predefined fixed photolysis frequencies at zenith angle 0
    CHARACTER(LEN=10), PARAMETER, DIMENSION(nconst) :: names_c =  (/                               &
                     'J_O31D    ','J_O33P    ','J_NO2     ','J_HNO3    ','J_RCHO    ',             &
                     'J         ','J         ','J         ','J         ','J         ',             &
                     'J         ','J         ','J         ','J         ','J         ' /)
!
!-- Photolysis frequency at zenith angle 0 degrees in 1/s
    REAL(wp), PARAMETER, DIMENSION(nconst) :: phot0 =  (/                                          &
                      2.489E-05_wp, 3.556E-04_wp, 8.89E-03_wp,5.334E-07_wp, 3.734E-05_wp,          &
                      0.0000E00_wp, 0.0000E00_wp, 0.0000E00_wp,0.0000E00_wp, 0.0000E00_wp,         &
                      0.0000E00_wp, 0.0000E00_wp, 0.0000E00_wp,0.0000E00_wp, 0.0000E00_wp /)
!
!-- Parameters for simple photolysis frequencies from MCM (http://mcm.leeds.ac.uk/MCM)
!-- Saunders et al., 2003, Atmos. Chem. Phys., 3, 161-180
    INTEGER,PARAMETER :: nsimple = 15               !< available predefined photolysis prequencies for simple parameterisation
!
!-- Names for simple photolysis frequencies parameterisation (
    CHARACTER(LEN=10), PARAMETER, DIMENSION(nsimple) :: names_s =  (/                              &
                     'J_O31D    ','J_O33P    ','J_H2O2    ','J_NO2     ','J_NO3_A   ',             &
                     'J_NO3_B   ','J_HONO    ','J_HNO3    ','J_HCHO_A  ','J_HCHO_B  ',             &
                     'J_CH3CHO  ','J         ','J         ','J         ','J_RCHO    ' /)
!
!-- Species dependent parameters for simple photolysis frequencies from MCM
!-- (http://mcm.leeds.ac.uk/MCM)
!-- J = l*COSx@m*EXP(-n*SECx)  with l,m,n named par_l etc., x is the zenith angle
    REAL(wp), PARAMETER, DIMENSION(nconst) :: par_l =  (/                                          &
                       6.073E-05_wp, 4.775E-04_wp, 1.041E-05_wp, 1.165E-02_wp, 2.485E-02_wp,       &
                       1.747E-01_wp, 2.644E-03_wp, 9.312E-07_wp, 4.642E-05_wp, 6.853E-05_wp,       &
                       7.344E-06_wp, 0.0000E00_wp, 0.0000E00_wp, 0.000E00_wp,  6.853E-05_wp /)

    REAL(wp), PARAMETER, DIMENSION(nconst) :: par_m =  (/                                          &
                           1.743_wp,    0.298_wp,    0.723_wp,    0.244_wp,    0.168_wp,           &
                           0.155_wp,    0.261_wp,    1.230_wp,    0.762_wp,    0.477_wp,           &
                           1.202_wp,    0.000_wp,    0.000_wp,    0.000_wp,    0.477_wp /)

    REAL(wp), PARAMETER, DIMENSION(nconst) :: par_n =  (/                                          &
                           0.474_wp,    0.080_wp,    0.279_wp,    0.267_wp,    0.108_wp,           &
                           0.125_wp,    0.288_wp,    0.307_wp,    0.353_wp,    0.323_wp,           &
                           0.417_wp,    0.000_wp,    0.000_wp,    0.000_wp,    0.323_wp /)


    REAL(wp)     :: cosz = 0.7_wp                   !< cosine of fixed zenith angle (45 deg, if not
                                                    !< specified otherwise)


    INTERFACE photolysis_constant
       MODULE PROCEDURE photolysis_constant
    END INTERFACE photolysis_constant

    INTERFACE photolysis_simple
       MODULE PROCEDURE photolysis_simple
    END INTERFACE photolysis_simple
!
!   INTERFACE photolysis_fastj
!      MODULE PROCEDURE photolysis_fastj
!   END INTERFACE photolysis_fastj
!
    INTERFACE photolysis_control
       MODULE PROCEDURE photolysis_control
    END INTERFACE photolysis_control

    SAVE

    PRIVATE

    PUBLIC  photolysis_control

    PUBLIC  photolysis_scheme

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine controls the calls of the photolysis schemes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE photolysis_control

    IMPLICIT NONE

    SELECT CASE ( TRIM( photolysis_scheme ) )

       CASE ( 'constant' )
          CALL photolysis_constant

       CASE ( 'simple' )
          CALL photolysis_simple

!      CASE ( 'fastj' )
!         CALL photolysis_fastj

       CASE DEFAULT

    END SELECT


 END SUBROUTINE photolysis_control


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This scheme keeps the prescribed net radiation constant during the run.
!> Default zenith angle is 45 deg
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE photolysis_constant

    IMPLICIT NONE

    INTEGER(iwp) :: iphot,iav !< loop index for photolysis reaction

    DO  iphot = 1, nphot
       DO  iav = 1, nconst
          IF ( TRIM( names_c(iav) ) == TRIM( phot_names(iphot) ) )  THEN
!--             Prescribe fixed photolysis frequencies  [1/s]
                phot_frequen(iphot)%freq(nzb+1:nzt,:,:) =  phot0(iav) * cosz
          ENDIF
       ENDDO
    ENDDO


 END SUBROUTINE photolysis_constant


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This scheme applies a simple parameterisation for clear sky photolysis frequencies from the
!> Master Chemical Mechanism, MCM v3.2 (http://mcm.leeds.ac.uk/MCM).
!> Reference: Saunders et al., Atmos. Chem. Phys., 3, 161, 2003
!> J = l*COSx@m*EXP(-n*SECx)  with l,m,n named par_l etc., x is the zenith angle
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE photolysis_simple

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    USE radiation_model_mod,                                                                       &
        ONLY:  calc_zenith,                                                                        &
               cos_zenith,                                                                         &
               rad_sw_in_diff,                                                                     &
               rad_sw_in_dir,                                                                      &
               rad_shade_h,                                                                        &
               solar_constant

    IMPLICIT NONE

    INTEGER(iwp) ::  day_of_year  !< day of the year
    INTEGER(iwp) ::  i            !< loop index x-direction
    INTEGER(iwp) ::  iav          !< loop indix for photolysis reaction
    INTEGER(iwp) ::  iphot        !< loop index for photolysis reaction
    INTEGER(iwp) ::  j            !< loop index y-direction
    INTEGER(iwp) ::  nhmax        !< level of maximum height of shade

    REAL(wp)     ::  cos_exp          !< pre-calculatated cosine power of
    REAL(wp)     ::  coszi            !< 1./cosine of zenith angle
    REAL(wp)     ::  cloud_factor     !< ratio between clear-sky and current value of solar radiation
    REAL(wp)     ::  exp_cos          !< pre-calculatated exponential function
    REAL(wp)     ::  second_of_day    !< second of the day
    REAL(wp)     ::  rad_sw_in_clear  !< incoming clear-sky solar radiation  (simple)
    REAL(wp)     ::  sky_trans        !< 1./cosine of zenith angle

    DO  iphot = 1, nphot
       phot_frequen(iphot)%freq = 0.0_wp
    ENDDO

    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year,                     &
                        second_of_day=second_of_day )
    CALL calc_zenith( day_of_year, second_of_day )

    IF ( cos_zenith > 0.0_wp )  THEN
       coszi = 1.0_wp / cos_zenith
       sky_trans = 0.6_wp + 0.2_wp * cos_zenith
       rad_sw_in_clear = solar_constant * sky_trans * cos_zenith

       DO  iphot = 1, nphot
          DO  iav = 1, nsimple

             IF ( TRIM( names_s(iav) ) == TRIM( phot_names(iphot) ) )  THEN
                exp_cos = EXP( - par_n(iav) * coszi )
                cos_exp = par_l(iav) * cos_zenith**par_m(iav)

                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      cloud_factor = ( rad_sw_in_dir(j,i) + rad_sw_in_diff(j,i) ) / rad_sw_in_clear
                      phot_frequen(iphot)%freq(nzb+1:nzt,j,i) = cos_exp * exp_cos * cloud_factor

                   ENDDO
                ENDDO
!
!--             Reduce photolysis frequency inside shadows (simple approach based on maximum shadow
!--             height only: reduce down to 20 %).
                IF ( photolysis_shading )  THEN
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         nhmax = rad_shade_h(j,i)
                         phot_frequen(iphot)%freq(nzb+1:nhmax,j,i) =                               &
                                                phot_frequen(iphot)%freq(nzb+1:nhmax,j,i) * 0.2_wp
                      ENDDO
                   ENDDO
                ENDIF

             ENDIF

          ENDDO
       ENDDO

    ENDIF

 END SUBROUTINE photolysis_simple

 END MODULE chem_photolysis_mod
