!> @file src/inifor_control.f90
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
!> The control module provides routines for timing INIFOR and writing runtime
!> feedback to the terminal and a log file.
!------------------------------------------------------------------------------!
 MODULE inifor_control

    USE inifor_defs,                                                           &
        ONLY:  COPYRIGHT, LNAME, LOG_FILE_NAME, PATH, VERSION, iwp, wp
    USE inifor_util,                                                           &
        ONLY:  real_to_str, real_to_str_f, str

    IMPLICIT NONE

    INTEGER(iwp), SAVE         ::  u                     !< Fortran file unit for the log file
    INTEGER(iwp), PARAMETER    ::  n_max_wrngs = 512     !< Fortran file unit for the log file
    INTEGER(iwp), SAVE         ::  n_wrngs = 0           !< Fortran file unit for the log file
    CHARACTER (LEN=5000)       ::  message = ''          !< log message buffer
    CHARACTER (LEN=5000)       ::  tip     = ''          !< optional log message buffer for tips on how to rectify encountered errors
    CHARACTER (LEN=5000), SAVE ::  warnings(n_max_wrngs) !< log of warnings

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> 
!> report() is INIFOR's general logging routine. It prints the given 'message'
!> to the terminal and writes it to the INIFOR log file.
!> 
!> You can use this routine to log events across INIFOR's code to log. For
!> warnings and abort messages, you may use the dedicated routines warn() and
!> inifor_abort() in this module. Both use report() and add specific behaviour
!> to it.
!------------------------------------------------------------------------------!
 SUBROUTINE report( routine, message, debug )

    CHARACTER(LEN=*), INTENT(IN)  ::  routine !< name of calling subroutine of function
    CHARACTER(LEN=*), INTENT(IN)  ::  message !< log message
    LOGICAL, OPTIONAL, INTENT(IN) ::  debug   !< flag the current message as debugging message

    LOGICAL, SAVE                 ::  is_first_run = .TRUE. !< control flag for file opening mode
    LOGICAL                       ::  suppress_message      !< control falg for additional debugging log

    IF ( is_first_run )  THEN
       OPEN( NEWUNIT=u, FILE=LOG_FILE_NAME, STATUS='replace' )
       is_first_run = .FALSE.
    ENDIF
       

    suppress_message = .FALSE.
    IF ( PRESENT( debug ) )  THEN
       IF ( .NOT. debug )  suppress_message = .TRUE.
    ENDIF

    IF ( .NOT. suppress_message )  THEN
       CALL write_to_sdtout_and_logfile(                                       &
          TRIM( message ) // "  [ " // TRIM( routine ) // " ]"                 &
       )
    ENDIF

 END SUBROUTINE report


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine is used to report whether a Boolean setting is enabled or
!> disabled as indicated by the 'toggle' attribute.
!------------------------------------------------------------------------------!
 SUBROUTINE report_toggle( routine, message, toggle )

    CHARACTER(LEN=*), INTENT(IN)  ::  routine !< name of calling subroutine of function
    CHARACTER(LEN=*), INTENT(IN)  ::  message !< log message
    LOGICAL, INTENT(IN)           ::  toggle  !< flag indicating whether to repot a setting as enabled or disabled

    IF ( toggle )  THEN
       CALL report( routine, message // " enabled" )
    ELSE
       CALL report( routine, message // " disabled" )
    ENDIF

 END SUBROUTINE report_toggle


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the given message to SDTOUT as well as to the INIFOR log
!> file.
!------------------------------------------------------------------------------!
 SUBROUTINE write_to_sdtout_and_logfile( message )

    CHARACTER(LEN=*), INTENT(IN)  ::  message

    WRITE(*, '(A)') "inifor: " // TRIM( message )
    WRITE(u, '(A)') TRIM( message )

 END SUBROUTINE write_to_sdtout_and_logfile


!------------------------------------------------------------------------------!
! Description:
! ------------
!> 
!> warn() prepends "WARNING:" the given 'message' and prints the result to the
!> terminal and writes it to the INIFOR logfile.
!> 
!> You can use this routine for messaging issues, that still allow INIFOR to
!> continue.
!------------------------------------------------------------------------------!
 SUBROUTINE warn( routine, message )

    CHARACTER(LEN=*), INTENT(IN) ::  routine !< name of calling subroutine or function
    CHARACTER(LEN=*), INTENT(IN) ::  message !< log message

    CALL cache_warning( routine, message )
    CALL report( routine, "WARNING: " // TRIM( message ) )

 END SUBROUTINE warn


 SUBROUTINE cache_warning( routine, message )

    CHARACTER(LEN=*), INTENT(IN) ::  routine !< name of calling subroutine or function
    CHARACTER(LEN=*), INTENT(IN) ::  message !< log message

    n_wrngs = n_wrngs + 1
    warnings(n_wrngs) = "  WARNING: " // TRIM( message ) //                      &
                        "  [ " // TRIM( routine ) // " ]"

 END SUBROUTINE cache_warning


!------------------------------------------------------------------------------!
! Description:
! ------------
!> 
!> This routine writes all warnings cached with cache_warning() to STDOUT
!> and the INIFOR log file.
!------------------------------------------------------------------------------!
 SUBROUTINE report_warnings()

    INTEGER(iwp)        ::  warning_idx
    CHARACTER (LEN=500) ::  warning = ''

    IF (n_wrngs > 0)  THEN
       warning = 'Encountered the following '// TRIM( str( n_wrngs ) ) // " warning(s) during this run:"
       CALL report( 'report_warnings', warning)

       DO warning_idx = 1, n_wrngs
          CALL write_to_sdtout_and_logfile( warnings(warning_idx) )
       ENDDO
    ENDIF

 END SUBROUTINE report_warnings

!------------------------------------------------------------------------------!
! Description:
! ------------
!> 
!> Report successful run. To be called at the end of the main loop.
!------------------------------------------------------------------------------!
 SUBROUTINE report_success( output_file_name )

    CHARACTER(LEN=PATH), INTENT(IN) ::  output_file_name

    message = "Finished writing dynamic driver '" // TRIM( output_file_name )
    message = TRIM( message ) // "' successfully."
    IF (n_wrngs > 0)  THEN
       message = TRIM( message ) // " Some warnings were encountered, see above."
    ENDIF
    CALL report( 'main loop', message )

 END SUBROUTINE report_success
   

!------------------------------------------------------------------------------!
! Description:
! ------------
!> 
!> Report runtime statistics
!------------------------------------------------------------------------------!
 SUBROUTINE report_runtime()

    CALL log_runtime( 'report', 'void' )

 END SUBROUTINE report_runtime


!------------------------------------------------------------------------------!
! Description:
! ------------
!> 
!> inifor_abort() prepends "ERROR:" the given 'message' and prints the result to
!> stdout, writes it to the INIFOR logfile, and exits INIFOR.
!> 
!> You can use this routine for messaging issues, that are critical and prevent
!> INIFOR from continueing.
!------------------------------------------------------------------------------!
 SUBROUTINE inifor_abort( routine , message )

    CHARACTER(LEN=*), INTENT(IN) ::  routine !< name of calling subroutine or function
    CHARACTER(LEN=*), INTENT(IN) ::  message !< log message

    CALL report_warnings
    CALL report( routine, "ERROR: " // TRIM( message ) // " Stopping." )
    CALL close_log
    CALL EXIT(1)

 END SUBROUTINE inifor_abort


 SUBROUTINE close_log()

    CLOSE( u )

 END SUBROUTINE close_log


!------------------------------------------------------------------------------!
! Description:
! ------------
!> 
!> print_version() prints the INIFOR version number and copyright notice.
!------------------------------------------------------------------------------!
 SUBROUTINE print_version()
    PRINT *, "INIFOR " // VERSION
    PRINT *, COPYRIGHT
 END SUBROUTINE print_version


!------------------------------------------------------------------------------!
! Description:
! ------------
!> 
!> log_runtime() measures the run times of various parts of INIFOR and
!> accumulates them in timing budgets.
!------------------------------------------------------------------------------!
 SUBROUTINE log_runtime( mode, budget )

    CHARACTER(LEN=*), INTENT(IN) ::  mode   !< name of the calling mode
    CHARACTER(LEN=*), INTENT(IN) ::  budget !< name of the timing budget

    REAL(wp), SAVE ::  t0               !< begin of timing interval
    REAL(wp), SAVE ::  t1               !< end of timing interval
    REAL(wp), SAVE ::  t_comp  = 0.0_wp !< computation timing budget
    REAL(wp), SAVE ::  t_alloc = 0.0_wp !< allocation timing budget
    REAL(wp), SAVE ::  t_init  = 0.0_wp !< initialization timing budget
    REAL(wp), SAVE ::  t_read  = 0.0_wp !< reading timing budget
    REAL(wp), SAVE ::  t_total = 0.0_wp !< total time
    REAL(wp), SAVE ::  t_write = 0.0_wp !< writing timing budget

    CHARACTER(LEN=*), PARAMETER  ::  fmt='(F6.2)' !< floating-point output format


    SELECT CASE( TRIM( mode ) )

    CASE( 'init' )
       CALL CPU_TIME(t0)

    CASE( 'time' )

       CALL CPU_TIME(t1)

       SELECT CASE( TRIM( budget ) )

          CASE( 'alloc' )
             t_alloc = t_alloc + t1 - t0

          CASE( 'init' )
             t_init = t_init + t1 - t0

          CASE( 'read' )
             t_read = t_read + t1 - t0

          CASE( 'write' )
             t_write = t_write + t1 - t0

          CASE( 'comp' )
             t_comp = t_comp + t1 - t0

          CASE DEFAULT
             CALL inifor_abort(                                                &
                'log_runtime',                                                 &
                "Time Budget '" // TRIM( mode ) // "' is not supported."       &
             )

       END SELECT

       t0 = t1

    CASE( 'report' )
        t_total = t_init + t_read + t_write + t_comp

        CALL report( 'log_runtime', "*** CPU time ***" )

        CALL report( 'log_runtime', "Initialization:  " // TRIM( real_to_str( t_init ) ) // &
                     " s  (" // TRIM( real_to_str( 100 * t_init / t_total, fmt ) ) // " %)" )

        CALL report( 'log_runtime', "(De-)Allocation: " // TRIM( real_to_str( t_alloc ) ) // &
                     " s  (" // TRIM( real_to_str( 100 * t_alloc / t_total, fmt ) ) // " %)" )

        CALL report( 'log_runtime', "Reading data:    " // TRIM( real_to_str( t_read ) )  // &
                     " s  (" // TRIM( real_to_str( 100 * t_read / t_total, fmt ) ) // " %)" )

        CALL report( 'log_runtime', "Writing data:    " // TRIM( real_to_str( t_write ) ) // &
                     " s  (" // TRIM( real_to_str( 100 * t_write / t_total, fmt ) ) // " %)" )

        CALL report( 'log_runtime', "Computation:     " // TRIM( real_to_str( t_comp ) )  // &
                     " s  (" // TRIM( real_to_str( 100 * t_comp / t_total, fmt) ) // " %)" )

        CALL report( 'log_runtime', "Total:           " // TRIM( real_to_str( t_total ) ) // &
                     " s  (" // TRIM( real_to_str( 100 * t_total / t_total, fmt ) ) // " %)" )

    CASE DEFAULT
       CALL inifor_abort( 'log_runtime', "Mode '" // TRIM(mode) // "' is not supported." )

    END SELECT

 END SUBROUTINE log_runtime


 END MODULE inifor_control
