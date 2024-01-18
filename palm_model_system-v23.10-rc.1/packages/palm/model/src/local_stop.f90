!> @file local_stop.f90
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
!> Stop program execution
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE local_stop

#if defined( __parallel )

    USE control_parameters,                                                                        &
        ONLY:  abort_mode,                                                                         &
               nested_run

    USE MPI

    USE pegrid

    USE pmc_interface,                                                                             &
        ONLY:  atmosphere_ocean_coupled_run


    IF ( nested_run  .OR.  atmosphere_ocean_coupled_run )  THEN
!
!--    Workaround: If any of the nested/coupled models crashes, it aborts the whole run with
!--                MPI_ABORT, regardless of the reason given by abort_mode.
       CALL MPI_ABORT( MPI_COMM_WORLD, 9999, ierr )
    ELSE
       IF ( abort_mode == 1 )  THEN
          CALL MPI_FINALIZE( ierr )
          STOP 1
       ELSEIF ( abort_mode == 2 )  THEN
          CALL MPI_ABORT( comm2d, 9999, ierr )
       ELSEIF ( abort_mode == 3 )  THEN
          CALL MPI_ABORT( MPI_COMM_WORLD, 9999, ierr )
       ENDIF
    ENDIF

#else

    STOP 1

#endif

 END SUBROUTINE local_stop
