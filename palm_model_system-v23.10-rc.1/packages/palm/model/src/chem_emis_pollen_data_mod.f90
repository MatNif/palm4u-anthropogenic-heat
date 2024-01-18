!> @file chem_emis_pollen_data_mod.f90
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
! Copyright 1997-2021 Leibniz Universitaet Hannover
! Copyright 2017-2021 Karlsruhe Institute of Technology
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Module for defnition of physical and parameteric attributes of all possible pollen species
!> Prefix epol
!--------------------------------------------------------------------------------------------------!
 MODULE chem_emis_pollen_data_mod

    USE chem_modules,                                                                              &
        ONLY:  nc_field_length

    USE kinds

    IMPLICIT NONE

    SAVE

    PUBLIC
!
!-- Number of pollens defined in this module.
    INTEGER(iwp), PARAMETER ::  pollen_default_offset = 1  !< default data index
    INTEGER(iwp), PARAMETER ::  pollen_num_species    = 4  !< num species (sans default)
    INTEGER(iwp), PARAMETER ::  pollen_num_total      =                                            &
                                pollen_num_species + pollen_default_offset  !< # species (w/ default)
!
!-- Instantiate pollen attributes (first element used as default).
!-- epol_species_names  : name of species per KPP mechanism
!-- pol_diameter_dry   : dry grain diameter                           [ 1e-6 m ]
!-- pol_diameter_ratio : dry diameter / mean diameter                 [  ]
!-- pol_density        : pollen density                               [ kg / m3 ]
!-- pol_qday           : flux for maximum available pollen            [ N / (m2 day) ] 
!-- pol_qt_duration    : duration of max pollen release               [ s ]
!-- pol_rand_coeff     : coefficients for random mechanical processes [ s ]
    CHARACTER(LEN=nc_field_length), PARAMETER,                                                     &
       DIMENSION(pollen_num_total) ::  epol_species_names =                                        &
       (/ '_DEFAULT',     'POL_BETU',    'POL_AMBR',    'POL_POAC',    'POL_ALNU'              /)

    REAL(wp), PARAMETER,                                                                           &
       DIMENSION(pollen_num_total) ::  pol_diameter_dry =                                          &
       (/  20.0,          23.0,          20.0,          36.0,          23.0                    /)

    REAL(wp), PARAMETER,                                                                           &
       DIMENSION(pollen_num_total) ::  pol_diameter_ratio =                                        &
       (/  1.0,           1.0,           1.0,           1.0,           0.987                   /)

    REAL(wp), PARAMETER,                                                                           &
       DIMENSION(pollen_num_total) ::  pol_density =                                               &
       (/  810.0,         810.0,         830.0,         980.0,         860.0                   /)

    REAL(wp), PARAMETER,                                                                           &
       DIMENSION(pollen_num_total) ::  pol_qday =                                                  &
       (/  1.0E5,         1.0166E7,       7.535E2,       1.09E5,        9.8E5                   /)

    REAL(wp), PARAMETER,                                                                           &
       DIMENSION(pollen_num_total) ::  pol_qt_duration =                                           &
       (/  57600.0,       57600.0,       32400.0,       57600.0,       57600.0                 /)

    REAL(wp), PARAMETER,                                                                           &
       DIMENSION(pollen_num_total) ::  pol_rand_coeff =                                            &
       (/  811.0,         43200.0,       811.0,         811.0,         43200.0                 /)

 END MODULE chem_emis_pollen_data_mod
