!> @file local_tremain_ini.f90
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
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Initialization of CPU-time measurements for different operating systems
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE local_tremain_ini


    USE cpulog,                                                                                    &
        ONLY:  initial_wallclock_time

    USE kinds

    IMPLICIT NONE

    INTEGER(idp)     ::  count      !<
    INTEGER(idp)     ::  count_rate !<

!
!-- Get initial wall clock time
    CALL SYSTEM_CLOCK( count, count_rate )
    initial_wallclock_time = REAL( count, KIND=wp ) / REAL( count_rate, KIND=wp )

 END SUBROUTINE local_tremain_ini
