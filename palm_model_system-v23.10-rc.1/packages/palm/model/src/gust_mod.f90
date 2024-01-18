!> @file gust_mod.f90
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
!> Gust model.
!>
!> @todo This is just a dummy module. The actual module ist not released yet.
!--------------------------------------------------------------------------------------------------!
 MODULE gust_mod

    USE control_parameters,                                                                        &
        ONLY:  restart_data_format_output

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nzb,                                                                                &
               nzt

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  idum  !< dummy variable used to avoid compiler warnings about unused variables

    LOGICAL ::  dummy_logical = .FALSE.        !< switch to avoid compiler warnings about unused variables
    LOGICAL ::  gust_module_enabled = .FALSE.  !< switch, if the entire module is used at all

    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC                                                                                         &
       gust_parin,                                                                                 &
       gust_check_parameters,                                                                      &
       gust_check_data_output_pr,                                                                  &
       gust_check_data_output,                                                                     &
       gust_init_arrays,                                                                           &
       gust_init,                                                                                  &
       gust_define_netcdf_grid,                                                                    &
       gust_header,                                                                                &
       gust_actions,                                                                               &
       gust_prognostic_equations,                                                                  &
       gust_swap_timelevel,                                                                        &
       gust_3d_data_averaging,                                                                     &
       gust_data_output_2d,                                                                        &
       gust_data_output_3d,                                                                        &
       gust_statistics,                                                                            &
       gust_rrd_global,                                                                            &
       gust_wrd_global,                                                                            &
       gust_rrd_local,                                                                             &
       gust_wrd_local
!
!-- Public parameters, constants and initial values
    PUBLIC                                                                                         &
       gust_module_enabled


    INTERFACE gust_parin
       MODULE PROCEDURE gust_parin
    END INTERFACE gust_parin

    INTERFACE gust_check_parameters
       MODULE PROCEDURE gust_check_parameters
    END INTERFACE gust_check_parameters

    INTERFACE gust_check_data_output_pr
       MODULE PROCEDURE gust_check_data_output_pr
    END INTERFACE gust_check_data_output_pr

    INTERFACE gust_check_data_output
       MODULE PROCEDURE gust_check_data_output
    END INTERFACE gust_check_data_output

    INTERFACE gust_init_arrays
       MODULE PROCEDURE gust_init_arrays
    END INTERFACE gust_init_arrays

    INTERFACE gust_init
       MODULE PROCEDURE gust_init
    END INTERFACE gust_init

    INTERFACE gust_define_netcdf_grid
       MODULE PROCEDURE gust_define_netcdf_grid
    END INTERFACE gust_define_netcdf_grid

    INTERFACE gust_header
       MODULE PROCEDURE gust_header
    END INTERFACE gust_header

    INTERFACE gust_actions
       MODULE PROCEDURE gust_actions
       MODULE PROCEDURE gust_actions_ij
    END INTERFACE gust_actions

    INTERFACE gust_prognostic_equations
       MODULE PROCEDURE gust_prognostic_equations
       MODULE PROCEDURE gust_prognostic_equations_ij
    END INTERFACE gust_prognostic_equations

    INTERFACE gust_swap_timelevel
       MODULE PROCEDURE gust_swap_timelevel
    END INTERFACE gust_swap_timelevel

    INTERFACE gust_3d_data_averaging
       MODULE PROCEDURE gust_3d_data_averaging
    END INTERFACE gust_3d_data_averaging

    INTERFACE gust_data_output_2d
       MODULE PROCEDURE gust_data_output_2d
    END INTERFACE gust_data_output_2d

    INTERFACE gust_data_output_3d
       MODULE PROCEDURE gust_data_output_3d
    END INTERFACE gust_data_output_3d

    INTERFACE gust_statistics
       MODULE PROCEDURE gust_statistics
    END INTERFACE gust_statistics

    INTERFACE gust_rrd_global
       MODULE PROCEDURE gust_rrd_global_ftn
       MODULE PROCEDURE gust_rrd_global_mpi
    END INTERFACE gust_rrd_global

    INTERFACE gust_wrd_global
       MODULE PROCEDURE gust_wrd_global
    END INTERFACE gust_wrd_global

    INTERFACE gust_rrd_local
       MODULE PROCEDURE gust_rrd_local_ftn
       MODULE PROCEDURE gust_rrd_local_mpi
    END INTERFACE gust_rrd_local

    INTERFACE gust_wrd_local
       MODULE PROCEDURE gust_wrd_local
    END INTERFACE gust_wrd_local

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &gust_parameters for gust module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_parin

    IMPLICIT NONE

    CHARACTER(LEN=100)  ::  line  !< dummy string that contains the current line of the parameter
                                  !< file

    INTEGER(iwp)  ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file


    NAMELIST /gust_parameters/  switch_off_module

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, gust_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    gust_parameters namelist was found and read correctly. Set flag that indicates that the gust
!--    module is switched on.
       IF ( .NOT. switch_off_module )  gust_module_enabled = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    gust_parameters namelist was found, but contained errors. Print an error message including
!--    the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'gust_parameters', line )

    ENDIF

 END SUBROUTINE gust_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for gust module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_check_parameters

    IMPLICIT NONE

 END SUBROUTINE gust_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output of profiles for gust module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_check_data_output_pr( variable, var_count, unit, dopr_unit )

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  dopr_unit !< local value of dopr_unit
    CHARACTER(LEN=*) ::  unit      !<
    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  var_count      !<

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = LEN( unit ) + LEN( variable ) + LEN( dopr_unit ) + var_count

 END SUBROUTINE gust_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for gust module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_check_data_output( var, unit )

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  unit  !<
    CHARACTER(LEN=*) ::  var   !<

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = LEN( var ) + LEN( unit )

 END SUBROUTINE gust_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate gust module arrays and define pointers
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_init_arrays

    IMPLICIT NONE

 END SUBROUTINE gust_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the gust module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_init

    IMPLICIT NONE

 END SUBROUTINE gust_init


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  grid_x      !<
    CHARACTER(LEN=*), INTENT(IN) ::  grid_y      !<
    CHARACTER(LEN=*), INTENT(IN) ::  grid_z      !<
    CHARACTER(LEN=*), INTENT(IN) ::  var         !<

    LOGICAL, INTENT(IN)           ::  found       !<


!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( found )  idum = LEN( var ) + LEN( grid_x ) + LEN( grid_y ) + LEN( grid_z )

 END SUBROUTINE gust_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for gust module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_header ( io )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  io  !< Unit of the output file


!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = io

 END SUBROUTINE gust_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_actions( location )

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  location !<

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = LEN( location )

 END SUBROUTINE gust_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_actions_ij( i, j, location )

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  location

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = i + j + LEN( location )

 END SUBROUTINE gust_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_prognostic_equations()

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = 1

 END SUBROUTINE gust_prognostic_equations


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_prognostic_equations_ij( i, j, i_omp_start, tn )

    INTEGER(iwp), INTENT(IN) ::  i            !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  i_omp_start  !< first loop index of i-loop in prognostic_equations
    INTEGER(iwp), INTENT(IN) ::  j            !< grid index in y-direction
    INTEGER(iwp), INTENT(IN) ::  tn           !< task number of openmp task

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = i + j + i_omp_start + tn

 END SUBROUTINE gust_prognostic_equations_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of timelevels
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_swap_timelevel ( mod_count )

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  mod_count


!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = mod_count

 END SUBROUTINE gust_swap_timelevel


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine for averaging 3D data
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_3d_data_averaging( mode, variable )

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  mode    !<
    CHARACTER(LEN=*) ::  variable !<

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = LEN( mode ) + LEN( variable )

 END SUBROUTINE gust_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 2D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(INOUT) ::  grid       !< name of vertical grid
    CHARACTER(LEN=*), INTENT(IN)    ::  mode       !< either 'xy', 'xz' or 'yz'
    CHARACTER(LEN=*), INTENT(IN)    ::  variable   !< name of variable

    INTEGER(iwp), INTENT(IN) ::  av        !< flag for (non-)average output
    INTEGER(iwp), INTENT(IN) ::  nzb_do    !< vertical output index (bottom)
    INTEGER(iwp), INTENT(IN) ::  nzt_do    !< vertical output index (top)

    LOGICAL, INTENT(INOUT) ::  found   !< flag if output variable is found
    LOGICAL, INTENT(INOUT) ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf !< local
                                                                                   !< array to which output data is resorted to

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( found .AND. two_d )  THEN
       idum = av + LEN( variable ) + LEN( grid // mode ) + local_pf(nxl,nys,nzb_do)
    ENDIF

 END SUBROUTINE gust_data_output_2d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine defining 3D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  variable   !< name of variable

    INTEGER(iwp), INTENT(IN) ::  av        !< flag for (non-)average output
    INTEGER(iwp), INTENT(IN) ::  nzb_do    !< lower limit of the data output (usually 0)
    INTEGER(iwp), INTENT(IN) ::  nzt_do    !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL, INTENT(INOUT) ::  found     !< flag if output variable is found

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf   !< local
                                                                                     !< array to which output data is resorted to

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( found )  idum = av + LEN( variable ) + local_pf(nxl,nys,nzb_do)

 END SUBROUTINE gust_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine computes profile and timeseries data for the gust module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_statistics( mode, sr, tn, dots_max )

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  mode  !<

    INTEGER(iwp) ::  dots_max   !<
    INTEGER(iwp) ::  sr         !<
    INTEGER(iwp) ::  tn         !<

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( dummy_logical )  idum = dots_max + sr + tn + LEN( mode )

 END SUBROUTINE gust_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_rrd_global_ftn( found )

    USE control_parameters,                                                                     &
        ONLY:  length,                                                                          &
               restart_string

    IMPLICIT NONE

    LOGICAL, INTENT(OUT)  ::  found


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'global_paramter' )
!          READ ( 13 )  global_parameter

       CASE DEFAULT
          found = .FALSE.

    END SELECT

 END SUBROUTINE gust_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_rrd_global_mpi

!    CALL rrd_mpi_io( 'global_parameter', global_parameter )
!          READ ( 13 )  global_parameter

 END SUBROUTINE gust_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync,   &
                                nyn_on_file, nysf, nysc, nys_on_file, tmp_2d, tmp_3d, found )

    USE control_parameters

    USE indices

    USE kinds

    USE pegrid


    IMPLICIT NONE

    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  nxlc            !<
    INTEGER(iwp) ::  nxlf            !<
    INTEGER(iwp) ::  nxl_on_file     !<
    INTEGER(iwp) ::  nxrc            !<
    INTEGER(iwp) ::  nxrf            !<
    INTEGER(iwp) ::  nxr_on_file     !<
    INTEGER(iwp) ::  nync            !<
    INTEGER(iwp) ::  nynf            !<
    INTEGER(iwp) ::  nyn_on_file     !<
    INTEGER(iwp) ::  nysc            !<
    INTEGER(iwp) ::  nysf            !<
    INTEGER(iwp) ::  nys_on_file     !<

    LOGICAL, INTENT(OUT)  ::  found

    REAL(wp), DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_2d   !<
    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d   !<


!
!-- Next lines are just to avoid compiler warnings about unused variables in case of empty user
!-- interface routine. You may remove them.
    IF ( dummy_logical )  THEN
       idum = k + nxlc + nxlf + nxrc + nxrf + nync + nynf + nysc + nysf +                          &
              tmp_2d(nys_on_file,nxl_on_file) + tmp_3d(nzb,nys_on_file,nxl_on_file)
    ENDIF

!
!-- Here the reading of user-defined restart data follows:
!-- Sample for user-defined output
    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( '.......' )
!          IF ( .NOT. ALLOCATED( u2_av ) ) THEN
!               ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!          ENDIF
!          IF ( k == 1 )  READ ( 13 )  tmp_3d
!             u2_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =         &
!                tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
       CASE DEFAULT
          found = .FALSE.

       END SELECT

 END SUBROUTINE gust_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_rrd_local_mpi

!    CALL rrd_mpi_io( 'local_array', local_array )

 END SUBROUTINE gust_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the gust module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_wrd_global

    IMPLICIT NONE


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

!       CALL wrd_write_string( 'global_parameter' )
!       WRITE ( 14 )  global_parameter

!       IF ( ALLOCATED( inflow_damping_factor ) )  THEN
!          CALL wrd_write_string( 'inflow_damping_factor' )
!          WRITE ( 14 )  inflow_damping_factor
!       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

!       CALL wrd_mpi_io( 'global_parameter', global_parameter )
!       IF ( ALLOCATED( inflow_damping_factor ) )  THEN
!          CALL wrd_mpi_io_global_array( 'inflow_damping_factor', inflow_damping_factor )
!       ENDIF

    ENDIF

 END SUBROUTINE gust_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data for the gust module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gust_wrd_local

    IMPLICIT NONE


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

!       IF ( ALLOCATED( u2_av ) )  THEN
!          CALL wrd_write_string( 'u2_av' )
!          WRITE ( 14 )  u2_av
!       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

!       IF ( ALLOCATED( u2_av ) )  CALL wrd_mpi_io( 'u2_av', u2_av )

    ENDIF

 END SUBROUTINE gust_wrd_local

 END MODULE gust_mod
