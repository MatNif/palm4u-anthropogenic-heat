!> @file cpulog_mod.f90
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
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> CPU-time measurements for any program part whatever. Results of the measurements are output at
!> the end of the run in local file CPU_MEASURES.
!>
!> To measure the CPU-time (better to say the wallclock time) of a specific code segment, two calls
!> of cpu_log have to be used as brackets in front and at the end of the segment:
!>
!>     CALL cpu_log( log_point(n), 'any identifier', 'start' )
!>       ... code segment ...
!>     CALL cpu_log( log_point(n), 'any identifier', 'stop' )
!>
!> Parts of the code segment can be excluded from the measurement by additional call of cpu_log:
!>
!>       ... first segment to be measured
!>     CALL cpu_log( log_point(n), 'any identifier', 'pause' )
!>       ... oart of segment to be excluded from measurement
!>     CALL cpu_log( log_point(n), 'any identifier', 'continue' )
!>       ... second segment to be mesasured
!>
!> n is an INTEGER within the interval [1,100] defining the id of the specific code segment,
!> 'any identifier' is a string describing the code segment to be measured. It can be freely chosen
!> and results will appear under this name in file CPU_MEASURES. ids can only be used once. If you
!> like to do a measurement of a new part of the code, please look for an id which is unused so far.
!>
!> runtime_parameters-parameter cpu_log_barrierwait can be used to set an MPI barrier at the
!> beginning of the measurement (modus 'start' or 'continue'), to avoid that idle times (due to MPI
!> calls in the code segment, which are waiting for other processes to be finished) affect the
!> measurements.
!> If barriers shall not be used at all, a fourth, optional parameter has to be
!> given:
!>
!>     CALL cpu_log( ..., ..., 'start', cpu_log_nowait )
!>
!> Variable log_point should be used for non-overlapping code segments, and they should sum up to
!> the total cpu-time required by the complete run.
!> Variable log_point_s can be used for any other special (s) measurements.
!--------------------------------------------------------------------------------------------------!
 MODULE cpulog

#if defined( __parallel )
    USE MPI
#endif

    USE control_parameters,                                                                        &
        ONLY:  message_string,                                                                     &
               nr_timesteps_this_run,                                                              &
               output_3d_file_size,                                                                &
               restart_data_format_output,                                                         &
               restart_file_size,                                                                  &
               run_description_header,                                                             &
               synchronous_exchange,                                                               &
               write_binary

    USE indices,                                                                                   &
        ONLY:  ngp_3d,                                                                             &
               nx,                                                                                 &
               ny,                                                                                 &
               nz

    USE kinds

    USE pegrid

    IMPLICIT NONE

    PRIVATE
    PUBLIC   cpu_log, cpu_log_barrierwait, cpu_log_nowait, cpu_statistics, initial_wallclock_time, &
             log_point, log_point_s

    INTERFACE cpu_log
       MODULE PROCEDURE cpu_log
    END INTERFACE cpu_log

    INTERFACE cpu_statistics
       MODULE PROCEDURE cpu_statistics
    END INTERFACE cpu_statistics

    INTEGER(iwp), PARAMETER ::  cpu_log_continue = 0  !<
    INTEGER(iwp), PARAMETER ::  cpu_log_pause = 1     !<
    INTEGER(iwp), PARAMETER ::  cpu_log_start = 2     !<
    INTEGER(iwp), PARAMETER ::  cpu_log_stop = 3      !<

    LOGICAL, PARAMETER ::  cpu_log_nowait = .FALSE.       !<

    LOGICAL            ::  cpu_log_barrierwait = .FALSE.  !<

    REAL(dp) ::  initial_wallclock_time  !<

    TYPE logpoint
       REAL(dp)           ::  isum       !<
       REAL(dp)           ::  ivect      !<
       REAL(dp)           ::  mean       !<
       REAL(dp)           ::  mtime      !<
       REAL(dp)           ::  mtimevec   !<
       REAL(dp)           ::  sum        !<
       REAL(dp)           ::  vector     !<
       INTEGER(iwp)       ::  counts     !<
       CHARACTER (LEN=30) ::  place      !<
    END TYPE logpoint

    TYPE(logpoint), DIMENSION(200) ::  log_point = logpoint( 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,       &
                                                             0.0_dp, 0.0_dp, 0.0_dp, 0, ' ' )
    TYPE(logpoint), DIMENSION(200) ::  log_point_s = logpoint( 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,     &
                                                               0.0_dp, 0.0_dp, 0.0_dp, 0, ' ' )

    SAVE

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE cpu_log( log_event, place, modus, barrierwait )

       IMPLICIT NONE

       CHARACTER (LEN=*) ::  modus              !<
       CHARACTER (LEN=*) ::  place              !<

       LOGICAL           ::  wait_allowed       !<
       LOGICAL, OPTIONAL ::  barrierwait        !<
       LOGICAL, SAVE     ::  first = .TRUE.     !<

       REAL(dp)          ::  mtime = 0.0_dp     !<
       REAL(dp)          ::  mtimevec = 0.0_dp  !<
       TYPE(logpoint)    ::  log_event          !<

       INTEGER(idp)     ::  count        !<
       INTEGER(idp)     ::  count_rate   !<


!
!--    Initialize and check, respectively, point of measurement
       IF ( log_event%place == ' ' )  THEN
          log_event%place = place
       ELSEIF ( log_event%place /= place )  THEN
          WRITE( message_string, * ) 'wrong argument: "', TRIM( place ), '", expected: "',         &
                                     TRIM(log_event%place), '"'
          CALL message( 'cpu_log', 'PAC0179', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Determine, if barriers are allowed to set
       IF ( PRESENT( barrierwait ) )  THEN
          wait_allowed = barrierwait
       ELSE
          wait_allowed = .TRUE.
       ENDIF

!
!--    MPI barrier, if requested, in order to avoid measuring wait times
!--    caused by MPI routines waiting for other MPI routines of other
!--    PEs that have not yet finished
#if defined( __parallel )
       IF ( cpu_log_barrierwait  .AND.  wait_allowed  .AND.                                        &
            ( modus == 'start'  .OR.  modus == 'continue' ) )  THEN
          CALL MPI_BARRIER( comm2d, ierr )
       ENDIF
#endif

!
!--    Take current time
       CALL SYSTEM_CLOCK( count, count_rate )
       mtime = REAL( count, KIND=dp ) / REAL( count_rate, KIND=dp )

!
!--    Start, stop or pause measurement
       IF ( modus == 'start'  .OR.  modus == 'continue' )  THEN
          log_event%mtime    = mtime
          log_event%mtimevec = mtimevec
       ELSEIF ( modus == 'pause' )  THEN
          IF ( ( mtime - log_event%mtime ) < 0.0  .AND.  first )  THEN
             WRITE( message_string, * ) 'negative time interval occured',                          &
                                        '&PE',myid,' L=PAUSE "',TRIM(log_event%place),             &
                                        '" new=', mtime,' last=',log_event%mtime
             CALL message( 'cpu_log', 'PAC0180', 0, 1, -1, 6, 0 )
             first = .FALSE.
          ENDIF
          log_event%isum     = log_event%isum + mtime - log_event%mtime
          log_event%ivect    = log_event%ivect  + mtimevec - log_event%mtimevec
       ELSEIF ( modus == 'stop' )  THEN
          IF ( ( mtime - log_event%mtime + log_event%isum ) < 0.0  .AND. first )  THEN
             WRITE( message_string, * ) 'negative time interval occured',                          &
                                        '&PE',myid,' L=STOP "',TRIM(log_event%place),'" new=',     &
                                        mtime,' last=',log_event%mtime,' isum=',log_event%isum
             CALL message( 'cpu_log', 'PAC0181', 0, 1, -1, 6, 0 )
             first = .FALSE.
          ENDIF
          log_event%mtime    = mtime    - log_event%mtime    + log_event%isum
          log_event%mtimevec = mtimevec - log_event%mtimevec + log_event%ivect
          log_event%sum      = log_event%sum  + log_event%mtime
          IF ( log_event%sum < 0.0  .AND.  first )  THEN
             WRITE( message_string, * ) 'negative time interval occured',                          &
                                        '&PE',myid,' L=STOP "',TRIM(log_event%place),'" sum=',     &
                                        log_event%sum,' mtime=',log_event%mtime
             CALL message( 'cpu_log', 'PAC0182', 0, 1, -1, 6, 0 )
             first = .FALSE.
          ENDIF
          log_event%vector   = log_event%vector + log_event%mtimevec
          log_event%counts   = log_event%counts + 1
          log_event%isum     = 0.0_dp
          log_event%ivect    = 0.0_dp
       ELSE
          message_string = 'unknown modus of time measurement: "' // TRIM( modus ) // '"'
          CALL message( 'cpu_log', 'PAC0183', 0, 1, -1, 6, 0 )
       ENDIF

    END SUBROUTINE cpu_log


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Analysis and output of the cpu-times measured. All PE results are collected
!> on PE0 in order to calculate the mean cpu-time over all PEs and other
!> statistics. The output is sorted according to the amount of cpu-time consumed
!> and output on PE0.
!------------------------------------------------------------------------------!

    SUBROUTINE cpu_statistics

       IMPLICIT NONE

       INTEGER(iwp)    ::  i               !<
       INTEGER(iwp)    ::  ii(1)           !<
#if defined( __parallel )
       INTEGER(iwp)    ::  iii             !<
       INTEGER(iwp)    ::  sender          !<
#endif
       REAL(dp)       ::  average_cputime  !<
       REAL(dp), SAVE ::  norm = 1.0_dp    !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  pe_max          !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  pe_max_s        !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  pe_min          !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  pe_min_s        !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  pe_rms          !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  pe_rms_s        !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  pe_tmp          !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  pe_tmp_s        !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  sum             !<
       REAL(dp), DIMENSION(:),   ALLOCATABLE ::  sum_s           !<
       REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  pe_log_points   !<
       REAL(dp), DIMENSION(:,:), ALLOCATABLE ::  pe_log_points_s !<


       CALL location_message( 'calculating cpu statistics', 'start' )

!
!--    Compute CPU-times in seconds for the global non-overlapping measurements.
       log_point%mtime  = log_point%mtime  / norm
       log_point%sum    = log_point%sum    / norm
       log_point%vector = log_point%vector / norm
       WHERE ( log_point%counts /= 0 )
          log_point%mean = log_point%sum / log_point%counts
       END WHERE


!
!--    Collect CPU-times from all PEs and calculate statistics
       IF ( myid == 0 )  THEN
!
!--       Allocate and initialize temporary arrays needed for statistics
          ALLOCATE( pe_max( SIZE( log_point ) ), pe_min( SIZE( log_point ) ),                      &
                    pe_rms( SIZE( log_point ) ), pe_tmp( SIZE( log_point ) ),                      &
                    pe_log_points( SIZE( log_point ), 0:numprocs-1 ) )
          pe_min = log_point%sum
          pe_max = log_point%sum    ! need to be set in case of 1 PE
          pe_rms = 0.0_dp
          pe_tmp = 0.0_dp

#if defined( __parallel )
!
!--       Receive data from all PEs
          DO  i = 1, numprocs-1
             CALL MPI_RECV( pe_tmp(1), SIZE( log_point ), MPI_DOUBLE_PRECISION, i, i, comm2d,      &
                            status, ierr )
             sender = status(MPI_SOURCE)
             pe_log_points(:,sender) = pe_tmp
          ENDDO
          pe_log_points(:,0) = log_point%sum   ! Results from PE0
!
!--       Calculate mean of all PEs, store it on log_point%sum
!--       and find minimum and maximum
          DO  iii = 1, SIZE( log_point )
             DO  i = 1, numprocs-1
                log_point(iii)%sum = log_point(iii)%sum + pe_log_points(iii,i)
                pe_min(iii) = MIN( pe_min(iii), pe_log_points(iii,i) )
                pe_max(iii) = MAX( pe_max(iii), pe_log_points(iii,i) )
             ENDDO
             log_point(iii)%sum = log_point(iii)%sum / numprocs
!
!--          Calculate rms
             DO  i = 0, numprocs-1
                pe_rms(iii) = pe_rms(iii) + ( pe_log_points(iii,i) - log_point(iii)%sum )**2
             ENDDO
             pe_rms(iii) = SQRT( pe_rms(iii) / numprocs )
          ENDDO
       ELSE
!
!--       Send data to PE0 (pe_max is used as temporary storage to send
!--       the data in order to avoid sending the data type log)
          ALLOCATE( pe_max( SIZE( log_point ) ) )
          pe_max = log_point%sum
          CALL MPI_SEND( pe_max(1), SIZE( log_point ), MPI_DOUBLE_PRECISION, 0, myid, comm2d, ierr )
#endif

       ENDIF

!
!--    The same procedure again for the special measurements.
!
!--    Compute cpu-times in seconds
       log_point_s%mtime  = log_point_s%mtime  / norm
       log_point_s%sum    = log_point_s%sum    / norm
       log_point_s%vector = log_point_s%vector / norm
       WHERE ( log_point_s%counts /= 0 )
          log_point_s%mean = log_point_s%sum / log_point_s%counts
       END WHERE

!
!--    Collect CPU-times from all PEs and calculate statistics
#if defined( __parallel )
!
!--    Set barrier in order to avoid that PE0 receives log_point_s-data
!--    while still busy with receiving log_point-data (see above)
       CALL MPI_BARRIER( comm2d, ierr )
#endif
       IF ( myid == 0 )  THEN
!
!--       Allocate and initialize temporary arrays needed for statistics
          ALLOCATE( pe_max_s( SIZE( log_point_s ) ), pe_min_s( SIZE( log_point_s ) ),              &
                    pe_rms_s( SIZE( log_point_s ) ), pe_tmp_s( SIZE( log_point_s ) ),              &
                    pe_log_points_s( SIZE( log_point_s ), 0:numprocs-1 ) )
          pe_min_s = log_point_s%sum
          pe_max_s = log_point_s%sum    ! need to be set in case of 1 PE
          pe_rms_s = 0.0_dp

#if defined( __parallel )
!
!--       Receive data from all PEs
          DO  i = 1, numprocs-1
             CALL MPI_RECV( pe_tmp_s(1), SIZE( log_point_s ), MPI_DOUBLE_PRECISION,                &
                            MPI_ANY_SOURCE, MPI_ANY_TAG, comm2d, status, ierr )
             sender = status(MPI_SOURCE)
             pe_log_points_s(:,sender) = pe_tmp_s
          ENDDO
          pe_log_points_s(:,0) = log_point_s%sum   ! Results from PE0
!
!--       Calculate mean of all PEs, store it on log_point_s%sum
!--       and find minimum and maximum
          DO  iii = 1, SIZE( log_point_s )
             DO  i = 1, numprocs-1
                log_point_s(iii)%sum = log_point_s(iii)%sum + pe_log_points_s(iii,i)
                pe_min_s(iii) = MIN( pe_min_s(iii), pe_log_points_s(iii,i) )
                pe_max_s(iii) = MAX( pe_max_s(iii), pe_log_points_s(iii,i) )
             ENDDO
             log_point_s(iii)%sum = log_point_s(iii)%sum / numprocs
!
!--          Calculate rms
             DO  i = 0, numprocs-1
                pe_rms_s(iii) = pe_rms_s(iii) + ( pe_log_points_s(iii,i) - log_point_s(iii)%sum )**2
             ENDDO
             pe_rms_s(iii) = SQRT( pe_rms_s(iii) / numprocs )
          ENDDO
       ELSE
!
!--       Send data to PE0 (pe_max is used as temporary storage to send
!--       the data in order to avoid sending the data type log)
          ALLOCATE( pe_max_s( SIZE( log_point_s ) ) )
          pe_max_s = log_point_s%sum
          CALL MPI_SEND( pe_max_s(1), SIZE( log_point_s ), MPI_DOUBLE_PRECISION, 0, 0, comm2d,     &
                         ierr )
#endif

       ENDIF
!
!--    Write CPU-times for the non-overlapping measurements.
       IF ( myid == 0 )  THEN
!
!--       Re-store sums
          ALLOCATE( sum( SIZE( log_point ) ) )
          WHERE ( log_point%counts /= 0 )
             sum = log_point%sum
          ELSEWHERE
             sum = -1.0_dp
          ENDWHERE

!
!--       Get total time in order to calculate CPU-time per gridpoint and timestep.
!--       Time for initialization (2) and spinup (41) are not included here!
          IF ( nr_timesteps_this_run /= 0 )  THEN
             average_cputime = log_point_s(10)%sum / REAL( ngp_3d(0), KIND=dp ) /                  &
                               REAL( nr_timesteps_this_run, KIND=dp ) * 1E6_dp     ! in micro-sec
          ELSE
             average_cputime = -1.0_dp
          ENDIF

!
!--       Write cpu-times sorted by size
          CALL check_open( 18 )
#if defined( __parallel )
          WRITE ( 18, 100 )  TRIM( run_description_header ), numprocs * threads_per_task, npex,    &
                             npey, threads_per_task, nx+1, ny+1, nz, nr_timesteps_this_run,        &
                             average_cputime

          WRITE ( 18, 110 )
#else
          WRITE ( 18, 100 )  TRIM( run_description_header ), numprocs * threads_per_task, 1, 1,    &
                             threads_per_task, nx+1, ny+1, nz, nr_timesteps_this_run,              &
                             average_cputime

          WRITE ( 18, 110 )
#endif
          DO
             ii = MAXLOC( sum )
             i = ii(1)
             IF ( sum(i) /= -1.0_dp )  THEN
                WRITE ( 18, 102 )  log_point(i)%place, log_point(i)%sum,                           &
                                   log_point(i)%sum / log_point(1)%sum * 100.0_dp,                 &
                                   log_point(i)%counts, pe_min(i), pe_max(i), pe_rms(i)
                sum(i) = -1.0_dp
             ELSE
                EXIT
             ENDIF
          ENDDO
       ENDIF

!
!--    Write CPU-times for special measurements.
       IF ( myid == 0 )  THEN
!
!--       Re-store sums
          ALLOCATE( sum_s( SIZE( log_point_s ) ) )
          WHERE ( log_point_s%counts /= 0 )
             sum_s = log_point_s%sum
          ELSEWHERE
             sum_s = -1.0_dp
          ENDWHERE

!
!--       Write cpu-times sorted by size
          WRITE ( 18, 101 )
          DO
             ii = MAXLOC( sum_s )
             i = ii(1)
             IF ( sum_s(i) /= -1.0_dp )  THEN
                WRITE ( 18, 102 )  log_point_s(i)%place, log_point_s(i)%sum,                       &
                                   log_point_s(i)%sum / log_point(1)%sum * 100.0_dp,               &
                                   log_point_s(i)%counts, pe_min_s(i), pe_max_s(i), pe_rms_s(i)
                sum_s(i) = -1.0_dp
             ELSE
                EXIT
             ENDIF
          ENDDO

!
!--       Output of handling of MPI operations
          IF ( collective_wait )  THEN
             WRITE ( 18, 103 )
          ELSE
             WRITE ( 18, 104 )
          ENDIF
          IF ( cpu_log_barrierwait )  WRITE ( 18, 111 )
          IF ( synchronous_exchange )  THEN
             WRITE ( 18, 105 )
          ELSE
             WRITE ( 18, 106 )
          ENDIF
!
!--       Output of restart data transfer speed
          IF ( write_binary  .AND.  restart_data_format_output(1:3) == 'mpi' )  THEN
             WRITE ( 18, 107 )  restart_file_size, restart_file_size / log_point(22)%sum
          ENDIF
!
!--       Output of NetCDF data transfer speed
          IF ( output_3d_file_size /= 0.0_wp  .AND.                                                &
               ( log_point_s(58)%sum /= 0.0_wp  .OR.  log_point_s(59)%sum /= 0.0_wp ) )            &
          THEN
             WRITE ( 18, 108 )  output_3d_file_size,                                               &
                                output_3d_file_size / ( log_point_s(58)%sum + log_point_s(59)%sum )
          ENDIF
!
!--       Empty lines in order to create a gap to the results from the next restart run
          WRITE ( 18, 109 )

!
!--       Unit 18 is not needed anymore
          CALL close_file( 18 )

       ENDIF

       CALL location_message( 'calculating cpu statistics', 'finished' )


   100 FORMAT (A/11('-')//'CPU measures for ',I5,' PEs (',I5,'(x) * ',I5,'(y',                     &
               &') tasks *',I5,' threads):'//                                                      &
               'gridpoints (x/y/z): ',20X,I5,' * ',I5,' * ',I5/                                    &
               'nr of timesteps: ',22X,I6/                                                         &
               'cpu time per grid point and timestep: ',5X,F8.5,' * 10**-6 s')

   101 FORMAT (/'special measures:'/ &
               &'-----------------------------------------------------------',                     &
               &'---------------------')

   102 FORMAT (A25,1X,F10.3,1X,F6.2,1X,I7,2(1X,F9.2),1X,F8.2)
   103 FORMAT (/'Barriers are set in front of collective operations')
   104 FORMAT (/'No barriers are set in front of collective operations')
   105 FORMAT (/'Exchange of ghostpoints via MPI_SENDRCV')
   106 FORMAT (/'Exchange of ghostpoints via MPI_ISEND/MPI_IRECV')
   107 FORMAT (/'Restart file size:   ',F12.1,' MByte'/                                            &
               &'I/O transfer speed:  ',F12.1,' MByte / sec')
   108 FORMAT (/'3d-NetCDF file size: ',F12.1,' MByte'/                                            &
               &'I/O transfer speed:  ',F12.1,' MByte / sec')
   109 FORMAT (//)
   110 FORMAT ('------------------------------------------------------------',                     &
               &'----------'//                                                                     &
               &'place:                             mean       counts     ',                       &
               &' min       max      rms'/                                                         &
               &'                                sec.    %                ',                       &
               &'sec.      sec.     sec.'/                                                         &
               &'-----------------------------------------------------------',                     &
               &'---------------------')
   111 FORMAT (/'Barriers are set at beginning (start/continue) of measurements')

    END SUBROUTINE cpu_statistics

 END MODULE cpulog
