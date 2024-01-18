!> @file src/inifor_defs.f90
!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2017-2021 Leibniz Universitaet Hannover
! Copyright 2017-2021 Deutscher Wetterdienst Offenbach
!------------------------------------------------------------------------------!
!
! Authors:
! --------
!> @author Eckhard Kadasch (Deutscher Wetterdienst, Offenbach)
!
! Description:
! ------------
!> The defs module provides global constants used in INIFOR.
!------------------------------------------------------------------------------!
 MODULE inifor_defs
 
 !USE kinds,                                                                    &
 !    ONLY :  wp, iwp

 IMPLICIT NONE

! 
!-- Parameters for type definitions
 INTEGER, PARAMETER  ::  iwp = 8
 INTEGER, PARAMETER  ::  wp = 8
 INTEGER, PARAMETER  ::  PATH  = 140 !< length of file path strings
 INTEGER, PARAMETER  ::  LNAME = 150 !< length of long name strings
 INTEGER, PARAMETER  ::  SNAME = 40  !< length of short name strings
 INTEGER, PARAMETER  ::  DATE  = 10  !< length of date strings

! 
!-- Trigonomentry
 REAL(wp), PARAMETER ::  PI = 3.14159265358979323846264338_wp !< Ratio of a circle's circumference to its diamter [-]
 REAL(wp), PARAMETER ::  TO_RADIANS = PI / 180.0_wp           !< Conversion factor from degrees to radiant [-]
 REAL(wp), PARAMETER ::  TO_DEGREES = 180.0_wp / PI           !< Conversion factor from radians to degrees [-]

! 
!-- COSMO parameters
 INTEGER(iwp), PARAMETER  ::  WATER_ID = 9                !< Integer corresponding to the water soil type in COSMO-DE [-]
 REAL(wp), PARAMETER ::  EARTH_RADIUS = 6371229.0_wp !< Earth radius used in COSMO-DE [m]
 REAL(wp), PARAMETER ::  P_SL = 1e5_wp               !< Reference pressure for computation of COSMO-DE's basic state pressure [Pa]
 REAL(wp), PARAMETER ::  T_SL = 288.15_wp            !< Reference temperature for computation of COSMO-DE's basic state pressure [K]
 REAL(wp), PARAMETER ::  BETA = 42.0_wp              !< logarithmic lapse rate, dT / d ln(p), for computation of COSMO-DE's basic
                                                     !< state pressure [K]
 REAL(wp), PARAMETER ::  RD   = 287.05_wp            !< specific gas constant of dry air, used in computation of COSMO-DE's basic
                                                     !< state [J/kg/K]
 REAL(wp), PARAMETER ::  RV   = 461.51_wp            !< specific gas constant of water vapor [J/kg/K]
 REAL(wp), PARAMETER ::  G    = 9.80665_wp           !< acceleration of Earth's gravity, used in computation of COSMO-DE's basic
                                                     !< state [m/s/s]
 REAL(wp), PARAMETER ::  RHO_L = 1e3_wp              !< density of liquid water, used to convert W_SO from [kg/m^2] to [m^3/m^3],
                                                     !< in [kg/m^3]
 REAL(wp), PARAMETER ::  HECTO = 100_wp              !< unit conversion factor from hPa to Pa

!
!-- PALM-4U parameters
 REAL(wp), PARAMETER ::  OMEGA   = 7.29e-5_wp !< angular velocity of Earth's rotation [s^-1]
 REAL(wp), PARAMETER ::  P_REF   = 1e5_wp     !< Reference pressure for potential temperature [Pa]
 REAL(wp), PARAMETER ::  RD_PALM = 287.0_wp   !< specific gas constant of dry air, used in computation of PALM-4U's potential temperature [J/kg/K]
 REAL(wp), PARAMETER ::  CP_PALM = 1005.0_wp  !< heat capacity of dry air at constant pressure, used in computation of PALM-4U's potential temperature [J/kg/K]

!
!-- PALM static driver attribute names (PIDS 1.9)
 CHARACTER(SNAME), PARAMETER ::  PIDS_ORIGIN_LON = 'origin_lon'
 CHARACTER(SNAME), PARAMETER ::  PIDS_ORIGIN_LAT = 'origin_lat'
 CHARACTER(SNAME), PARAMETER ::  PIDS_ORIGIN_Z   = 'origin_z'

! 
!-- COSMO netCDF parameters
 INTEGER, PARAMETER          ::  NC_DEPTH_DIM_IDX = 3
 CHARACTER(SNAME), PARAMETER ::  NC_DEPTH_NAME = 'depth_2'
 CHARACTER(SNAME), PARAMETER ::  NC_HHL_NAME = 'HHL'
 CHARACTER(SNAME), PARAMETER ::  NC_RLAT_NAME = 'rlat'
 CHARACTER(SNAME), PARAMETER ::  NC_RLON_NAME = 'rlon'
 CHARACTER(SNAME), PARAMETER ::  NC_ROTATED_POLE_NAME = 'rotated_pole'
 CHARACTER(SNAME), PARAMETER ::  NC_POLE_LATITUDE_NAME = 'grid_north_pole_latitude'
 CHARACTER(SNAME), PARAMETER ::  NC_POLE_LONGITUDE_NAME = 'grid_north_pole_longitude'

!
!-- INIFOR parameters
 CHARACTER(LEN=*), PARAMETER ::  CFG_INIT_PROFILE = 'profile'
 CHARACTER(LEN=*), PARAMETER ::  CFG_INIT_VOLUME = 'volume'
 CHARACTER(LEN=*), PARAMETER ::  CFG_INIT_SOIL_PROFILE = 'profile'
 CHARACTER(LEN=*), PARAMETER ::  CFG_INIT_SOIL_VOLUME = 'volume'
 CHARACTER(LEN=*), PARAMETER ::  CFG_FORCING_HETERO = 'hetero'
 CHARACTER(LEN=*), PARAMETER ::  CFG_FORCING_HOMO = 'homo'
 CHARACTER(LEN=*), PARAMETER ::  CFG_FORCING_NUDGING = 'nudging'
 INTEGER(iwp), PARAMETER     ::  FILL_ITERATIONS = 5          !< Number of iterations for extrapolating soil data into COSMO-DE
                                                              !< water cells [-]
 INTEGER(iwp), PARAMETER     ::  FORCING_STEP = 1             !< Number of hours between forcing time steps [h]
 REAL(wp), PARAMETER         ::  NUDGING_TAU = 21600.0_wp     !< Nudging relaxation time scale [s]
 CHARACTER(LEN=*), PARAMETER ::  COPYRIGHT = 'Copyright 2017-2021 Leibniz Universitaet Hannover' // &
    ACHAR( 10 ) // ' Copyright 2017-2021 Deutscher Wetterdienst Offenbach' !< Copyright notice
 CHARACTER(LEN=*), PARAMETER ::  LOG_FILE_NAME = 'inifor.log' !< Name of INIFOR's log file
 CHARACTER(LEN=*), PARAMETER ::  VERSION = '2.2.1'            !< INIFOR version number
 
 END MODULE inifor_defs
