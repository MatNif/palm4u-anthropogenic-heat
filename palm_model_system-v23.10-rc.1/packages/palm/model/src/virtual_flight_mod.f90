!> @file virtual_flights_mod.f90
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
! @author Matthias Suehring
!
!
! Description:
! ------------
!> Module for virtual flight measurements.
!> @todo Err msg VFL0005: flight can be inside topography -> extra check?
!--------------------------------------------------------------------------------------------------!
 MODULE flight_mod

#if defined( __parallel )
    USE MPI
#endif

    USE control_parameters,                                                                        &
        ONLY:  debug_output,                                                                       &
               fl_max, num_leg,                                                                    &
               num_var_fl,                                                                         &
               num_var_fl_user,                                                                    &
               restart_data_format_output,                                                         &
               virtual_flight

    USE kinds

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array,                                                              &
               rrd_mpi_io_global_array,                                                            &
               wrd_mpi_io_global_array

    USE user_init_flight_mod,                                                                      &
        ONLY:  user_init_flight


    CHARACTER(LEN=6), DIMENSION(fl_max) ::  leg_mode = 'cyclic'  !< flight mode through the model domain, either 'cyclic' or
                                                                 !<'return'

    INTEGER(iwp) ::  l           !< index for flight leg
    INTEGER(iwp) ::  var_index   !< index for measured variable

    LOGICAL, DIMENSION(:), ALLOCATABLE  ::  cyclic_leg  !< flag to identify fly mode

    REAL(wp) ::  flight_end = 9999999.9_wp  !< end time of virtual flight
    REAL(wp) ::  flight_begin = 0.0_wp      !< end time of virtual flight

    REAL(wp), DIMENSION(fl_max) ::  flight_angle = 45.0_wp    !< angle determining the horizontal flight direction
    REAL(wp), DIMENSION(fl_max) ::  flight_level = 100.0_wp   !< flight level
    REAL(wp), DIMENSION(fl_max) ::  max_elev_change = 0.0_wp  !< maximum elevation change for the respective flight leg
    REAL(wp), DIMENSION(fl_max) ::  rate_of_climb = 0.0_wp    !< rate of climb or descent
    REAL(wp), DIMENSION(fl_max) ::  speed_agl = 25.0_wp       !< absolute horizontal flight speed above ground level (agl)
    REAL(wp), DIMENSION(fl_max) ::  x_start = 999999999.0_wp  !< start x position
    REAL(wp), DIMENSION(fl_max) ::  x_end   = 999999999.0_wp  !< end x position
    REAL(wp), DIMENSION(fl_max) ::  y_start = 999999999.0_wp  !< start y position
    REAL(wp), DIMENSION(fl_max) ::  y_end   = 999999999.0_wp  !< end y position

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_agl  !< u-component of flight speed
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_agl  !< v-component of flight speed
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  w_agl  !< w-component of flight speed
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  x_pos  !< aircraft x-position
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  y_pos  !< aircraft y-position
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  z_pos  !< aircraft z-position

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sensor_l  !< measured data on local PE
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sensor    !< measured data

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  var_u  !< dummy array for possibly user-defined quantities

    SAVE

    PRIVATE

    INTERFACE flight_header
       MODULE PROCEDURE flight_header
    END INTERFACE flight_header

    INTERFACE flight_init
       MODULE PROCEDURE flight_init
    END INTERFACE flight_init

    INTERFACE flight_init_output
       MODULE PROCEDURE flight_init_output
    END INTERFACE flight_init_output

    INTERFACE flight_check_parameters
       MODULE PROCEDURE flight_check_parameters
    END INTERFACE flight_check_parameters

    INTERFACE flight_parin
       MODULE PROCEDURE flight_parin
    END INTERFACE flight_parin

    INTERFACE interpolate_xyz
       MODULE PROCEDURE interpolate_xyz
    END INTERFACE interpolate_xyz

    INTERFACE flight_measurement
       MODULE PROCEDURE flight_measurement
    END INTERFACE flight_measurement

    INTERFACE flight_rrd_global
       MODULE PROCEDURE flight_rrd_global_ftn
       MODULE PROCEDURE flight_rrd_global_mpi
    END INTERFACE flight_rrd_global

    INTERFACE flight_wrd_global
       MODULE PROCEDURE flight_wrd_global
    END INTERFACE flight_wrd_global

!
!-- Private interfaces
    PRIVATE flight_check_parameters,                                                               &
            flight_init_output,                                                                    &
            interpolate_xyz
!
!-- Public interfaces
    PUBLIC flight_init,                                                                            &
           flight_header,                                                                          &
           flight_parin,                                                                           &
           flight_measurement,                                                                     &
           flight_wrd_global,                                                                      &
           flight_rrd_global
!
!-- Public variables
    PUBLIC fl_max,                                                                                 &
           sensor,                                                                                 &
           x_pos,                                                                                  &
           y_pos,                                                                                  &
           z_pos

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for flight module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE flight_header ( io )


    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  io  !< Unit of the output file

    WRITE ( io, 1  )
    WRITE ( io, 2  )
    WRITE ( io, 3  ) num_leg
    WRITE ( io, 4  ) flight_begin
    WRITE ( io, 5  ) flight_end

    DO  l=1, num_leg
       WRITE ( io, 6   ) l
       WRITE ( io, 7   ) speed_agl(l)
       WRITE ( io, 8   ) flight_level(l)
       WRITE ( io, 9   ) max_elev_change(l)
       WRITE ( io, 10  ) rate_of_climb(l)
       WRITE ( io, 11  ) leg_mode(l)
    ENDDO


 1   FORMAT (' Virtual flights: ----------------' )
 2   FORMAT ('       Output every timestep' )
 3   FORMAT ('       Number of flight legs:',    I3 )
 4   FORMAT ('       Begin of measurements:',    F10.1    , ' s' )
 5   FORMAT ('       End of measurements:',      F10.1    , ' s' )
 6   FORMAT ('       Leg', I3/, '       ------' )
 7   FORMAT ('          Flight speed            : ', F5.1, ' m/s' )
 8   FORMAT ('          Flight level            : ', F5.1, ' m' )
 9   FORMAT ('          Maximum elevation change: ', F5.1, ' m/s' )
 10  FORMAT ('          Rate of climb / descent : ', F5.1, ' m/s' )
 11  FORMAT ('          Leg mode                : ', A/ )

 END SUBROUTINE flight_header

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads the namelist flight_par.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE flight_parin

    IMPLICIT NONE

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /virtual_flight_parameters/  flight_angle,                                            &
                                          flight_begin,                                            &
                                          flight_end,                                              &
                                          flight_level,                                            &
                                          leg_mode,                                                &
                                          max_elev_change,                                         &
                                          rate_of_climb,                                           &
                                          speed_agl,                                               &
                                          switch_off_module,                                       &
                                          x_end,                                                   &
                                          x_start,                                                 &
                                          y_end,                                                   &
                                          y_start


!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, virtual_flight_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    virtual_flight_parameters namelist was found and read correctly. Set switch that virtual
!--    flights are carried out.
       IF ( .NOT. switch_off_module )  virtual_flight = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    virtual_flight_parameters namelist was found, but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'virtual_flight_parameters', line )

    ENDIF

 END SUBROUTINE flight_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Inititalization of required arrays, number of legs and flags. Moreover, initialize flight speed
!> and -direction, as well as start positions.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE flight_init

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  pi

    USE control_parameters,                                                                        &
        ONLY:  initializing_actions

    USE indices,                                                                                   &
        ONLY:  nxlg,                                                                               &
               nxrg,                                                                               &
               nysg,                                                                               &
               nyng,                                                                               &
               nzb,                                                                                &
               nzt

    IMPLICIT NONE

    REAL(wp) ::  distance  !< distance between start and end position of a flight leg


    IF ( debug_output )  CALL debug_message( 'flight_init', 'start' )
!
!-- Determine the number of flight legs
    l = 1
    DO  WHILE ( x_start(l) /= 999999999.0_wp  .AND.  l <= SIZE(x_start) )
       l       = l + 1
    ENDDO
    num_leg = l-1
!
!-- Check for proper parameter settings
    CALL flight_check_parameters
!
!-- Allocate and initialize logical array for flight pattern
    ALLOCATE( cyclic_leg(1:num_leg) )
!
!-- Initialize flags for cyclic/return legs
    DO  l = 1, num_leg
       cyclic_leg(l) = MERGE( .TRUE., .FALSE., TRIM( leg_mode(l) ) == 'cyclic' )
    ENDDO
!
!-- Allocate and initialize arraxs for flight position and speed. In case of restart runs these data
!-- are read by the routine read_flight_restart_data instead.
    IF (  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

       ALLOCATE( x_pos(1:num_leg), y_pos(1:num_leg ), z_pos(1:num_leg) )
!
!--    Inititalize x-, y-, and z-positions with initial start position
       x_pos(1:num_leg) = x_start(1:num_leg)
       y_pos(1:num_leg) = y_start(1:num_leg)
       z_pos(1:num_leg) = flight_level(1:num_leg)
!
!--    Allocate arrays for flight-speed components
       ALLOCATE( u_agl(1:num_leg),                                                                 &
                 v_agl(1:num_leg),                                                                 &
                 w_agl(1:num_leg) )
!
!--    Inititalize u-, v- and w-component.
       DO  l = 1, num_leg
!
!--       In case of return-legs, the flight direction, i.e. the horizontal flight-speed components,
!--       are derived from the given start/end positions.
          IF (  .NOT.  cyclic_leg(l) )  THEN
             distance = SQRT( ( x_end(l) - x_start(l) )**2 + ( y_end(l) - y_start(l) )**2 )
             u_agl(l) = speed_agl(l) * ( x_end(l) - x_start(l) ) / distance
             v_agl(l) = speed_agl(l) * ( y_end(l) - y_start(l) ) / distance
             w_agl(l) = rate_of_climb(l)
!
!--       In case of cyclic-legs, flight direction is directly derived from the given flight angle.
          ELSE
             u_agl(l) = speed_agl(l) * COS( flight_angle(l) * pi / 180.0_wp )
             v_agl(l) = speed_agl(l) * SIN( flight_angle(l) * pi / 180.0_wp )
             w_agl(l) = rate_of_climb(l)
          ENDIF

       ENDDO

    ENDIF
!
!-- Initialized data output
    CALL flight_init_output
!
!-- Allocate array required for user-defined quantities if necessary.
    IF ( num_var_fl_user  > 0 )  ALLOCATE( var_u(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!-- Allocate and initialize arrays containing the measured data
    ALLOCATE( sensor_l(1:num_var_fl,1:num_leg) )
    ALLOCATE( sensor(1:num_var_fl,1:num_leg)   )
    sensor_l = 0.0_wp
    sensor   = 0.0_wp

    IF ( debug_output )  CALL debug_message( 'flight_init', 'end' )

 END SUBROUTINE flight_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of output-variable names and units.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE flight_init_output

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model

    USE control_parameters,                                                                        &
        ONLY:  cloud_droplets,                                                                     &
               humidity,                                                                           &
               neutral,                                                                            &
               passive_scalar

    USE netcdf_interface

    IMPLICIT NONE

    CHARACTER(LEN=10) ::  label_leg  !< dummy argument to convert integer to string

    INTEGER(iwp) ::  i         !< loop variable
    INTEGER(iwp) ::  id_pt     !< identifyer for labeling
    INTEGER(iwp) ::  id_q      !< identifyer for labeling
    INTEGER(iwp) ::  id_ql     !< identifyer for labeling
    INTEGER(iwp) ::  id_s      !< identifyer for labeling
    INTEGER(iwp) ::  id_u = 1  !< identifyer for labeling
    INTEGER(iwp) ::  id_v = 2  !< identifyer for labeling
    INTEGER(iwp) ::  id_w = 3  !< identifyer for labeling
    INTEGER(iwp) ::  k         !< index variable

    LOGICAL      ::  init = .TRUE.  !< flag to distiquish calls of user_init_flight
!
!-- Define output quanities, at least three variables are measured (u,v,w)
    num_var_fl = 3
    IF ( .NOT. neutral                     )  THEN
       num_var_fl = num_var_fl + 1
       id_pt      = num_var_fl
    ENDIF
    IF ( humidity                          )  THEN
       num_var_fl = num_var_fl + 1
       id_q       = num_var_fl
    ENDIF
    IF ( bulk_cloud_model .OR. cloud_droplets )  THEN
       num_var_fl = num_var_fl + 1
       id_ql      = num_var_fl
    ENDIF
    IF ( passive_scalar                    )  THEN
       num_var_fl = num_var_fl + 1
       id_s       = num_var_fl
    ENDIF
!
!-- Write output strings for the spatial variables x, y, z
    DO  l=1, num_leg

       IF ( l < 10                    )  WRITE( label_leg, '(I1)' )  l
       IF ( l >= 10   .AND.  l < 100  )  WRITE( label_leg, '(I2)' )  l
       IF ( l >= 100  .AND.  l < 1000 )  WRITE( label_leg, '(I3)' )  l

       dofl_label_x(l)  = 'x_' // TRIM( label_leg )
       dofl_label_y(l)  = 'y_' // TRIM( label_leg )
       dofl_label_z(l)  = 'z_' // TRIM( label_leg )

    ENDDO

!
!-- Call user routine to initialize further variables
    CALL user_init_flight( init )
!
!-- Write output labels and units for the quanities
    k = 1
    DO  l=1, num_leg

       IF ( l < 10                    )  WRITE( label_leg, '(I1)' )  l
       IF ( l >= 10   .AND.  l < 100  )  WRITE( label_leg, '(I2)' )  l
       IF ( l >= 100  .AND.  l < 1000 )  WRITE( label_leg, '(I3)' )  l

       label_leg = 'leg_' // TRIM(label_leg)
       DO  i=1, num_var_fl

          IF ( i == id_u      )  THEN
             dofl_label(k) = TRIM( label_leg ) // '_u'
             dofl_unit(k)  = 'm/s'
             k             = k + 1
          ELSEIF ( i == id_v  )  THEN
             dofl_label(k) = TRIM( label_leg ) // '_v'
             dofl_unit(k)  = 'm/s'
             k             = k + 1
          ELSEIF ( i == id_w  )  THEN
             dofl_label(k) = TRIM( label_leg ) // '_w'
             dofl_unit(k)  = 'm/s'
             k             = k + 1
          ELSEIF ( i == id_pt )  THEN
             dofl_label(k) = TRIM( label_leg ) // '_theta'
             dofl_unit(k)  = 'K'
             k             = k + 1
          ELSEIF ( i == id_q  )  THEN
             dofl_label(k) = TRIM( label_leg ) // '_q'
             dofl_unit(k)  = 'kg/kg'
             k             = k + 1
          ELSEIF ( i == id_ql )  THEN
             dofl_label(k) = TRIM( label_leg ) // '_ql'
             dofl_unit(k)  = 'kg/kg'
             k             = k + 1
          ELSEIF ( i == id_s  )  THEN
             dofl_label(k) = TRIM( label_leg ) // '_s'
             dofl_unit(k)  = 'kg/kg'
             k             = k + 1
          ENDIF
       ENDDO

       DO i=1, num_var_fl_user
          CALL user_init_flight( init, k, i, label_leg )
       ENDDO

    ENDDO
!
!-- Finally, set the total number of flight-output quantities.
    num_var_fl = num_var_fl + num_var_fl_user

 END SUBROUTINE flight_init_output

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine calculates the current flight positions and calls the respective interpolation
!> routine to measure the data at the current flight position.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE flight_measurement

    USE arrays_3d,                                                                                 &
        ONLY:  ddzu,                                                                               &
               ddzw,                                                                               &
               pt,                                                                                 &
               q,                                                                                  &
               ql,                                                                                 &
               s,                                                                                  &
               u,                                                                                  &
               v,                                                                                  &
               w,                                                                                  &
               zu,                                                                                 &
               zw

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model

    USE control_parameters,                                                                        &
        ONLY:  cloud_droplets,                                                                     &
               dt_3d,                                                                              &
               humidity,                                                                           &
               neutral,                                                                            &
               passive_scalar,                                                                     &
               time_since_reference_point

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               nxl,                                                                                &
               nxr,                                                                                &
               ny,                                                                                 &
               nys,                                                                                &
               nyn

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< index of current grid box along x
    INTEGER(iwp) ::  j  !< index of current grid box along y
    INTEGER(iwp) ::  n  !< loop variable for number of user-defined output quantities

    LOGICAL  ::  on_pe  !< flag to check if current flight position is on current PE

    CALL cpu_log( log_point(65), 'flight_measurement', 'start' )
!
!-- Perform flight measurement if start time is reached.
    IF ( time_since_reference_point >= flight_begin  .AND.  time_since_reference_point <= flight_end )  THEN

       sensor_l = 0.0_wp
       sensor   = 0.0_wp
!
!--    Loop over all flight legs
       DO  l = 1, num_leg
!
!--       Update location for each leg
          x_pos(l) = x_pos(l) + u_agl(l) * dt_3d
          y_pos(l) = y_pos(l) + v_agl(l) * dt_3d
          z_pos(l) = z_pos(l) + w_agl(l) * dt_3d
!
!--       Check if location must be modified for return legs.
!--       Carry out horizontal reflection if required.
          IF ( .NOT. cyclic_leg(l) )  THEN

             IF ( x_start(l) <= x_end(l) )  THEN
!
!--             Outward flight, i.e. from start to end
                IF ( u_agl(l) >= 0.0_wp  .AND.  x_pos(l) > x_end(l)      )  THEN
                   x_pos(l) = 2.0_wp * x_end(l)   - x_pos(l)
                   u_agl(l) = - u_agl(l)
!
!--             Return flight, i.e. from end to start
                ELSEIF ( u_agl(l) < 0.0_wp  .AND.  x_pos(l) < x_start(l) )  THEN
                   x_pos(l) = 2.0_wp * x_start(l) - x_pos(l)
                   u_agl(l) = - u_agl(l)
                ENDIF
             ELSE
!
!--             Outward flight, i.e. from start to end
                IF ( u_agl(l) < 0.0_wp  .AND.  x_pos(l) < x_end(l)      )  THEN
                   x_pos(l) = 2.0_wp * x_end(l)   - x_pos(l)
                   u_agl(l) = - u_agl(l)
!
!--             Return flight, i.e. from end to start
                ELSEIF ( u_agl(l) >= 0.0_wp  .AND.  x_pos(l) > x_start(l) )  THEN
                   x_pos(l) = 2.0_wp * x_start(l) - x_pos(l)
                   u_agl(l) = - u_agl(l)
                ENDIF
             ENDIF

             IF ( y_start(l) <= y_end(l) )  THEN
!
!--             Outward flight, i.e. from start to end
                IF ( v_agl(l) >= 0.0_wp  .AND.  y_pos(l) > y_end(l)      )  THEN
                   y_pos(l) = 2.0_wp * y_end(l)   - y_pos(l)
                   v_agl(l) = - v_agl(l)
!
!--             Return flight, i.e. from end to start
                ELSEIF ( v_agl(l) < 0.0_wp  .AND.  y_pos(l) < y_start(l) )  THEN
                   y_pos(l) = 2.0_wp * y_start(l) - y_pos(l)
                   v_agl(l) = - v_agl(l)
                ENDIF
             ELSE
!
!--             Outward flight, i.e. from start to end
                IF ( v_agl(l) < 0.0_wp  .AND.  y_pos(l) < y_end(l)      )  THEN
                   y_pos(l) = 2.0_wp * y_end(l)   - y_pos(l)
                   v_agl(l) = - v_agl(l)
!
!--             Return flight, i.e. from end to start
                ELSEIF ( v_agl(l) >= 0.0_wp  .AND.  y_pos(l) > y_start(l) )  THEN
                   y_pos(l) = 2.0_wp * y_start(l) - y_pos(l)
                   v_agl(l) = - v_agl(l)
                ENDIF
             ENDIF
!
!--       Check if flight position is outside the model domain and apply cyclic conditions if required
          ELSEIF ( cyclic_leg(l) )  THEN
!
!--          Check if aircraft leaves the model domain at the right boundary
             IF ( ( flight_angle(l) >= 0.0_wp     .AND.                                            &
                    flight_angle(l) <= 90.0_wp )  .OR.                                             &
                  ( flight_angle(l) >= 270.0_wp   .AND.                                            &
                    flight_angle(l) <= 360.0_wp ) )  THEN
                IF ( x_pos(l) >= ( nx + 0.5_wp ) * dx )                                            &
                     x_pos(l) = x_pos(l) - ( nx + 1 ) * dx
!
!--          Check if aircraft leaves the model domain at the left boundary
             ELSEIF ( flight_angle(l) > 90.0_wp  .AND.  flight_angle(l) < 270.0_wp )  THEN
                IF ( x_pos(l) < -0.5_wp * dx )                                                     &
                     x_pos(l) = ( nx + 1 ) * dx + x_pos(l)
             ENDIF
!
!--          Check if aircraft leaves the model domain at the north boundary
             IF ( flight_angle(l) >= 0.0_wp  .AND.  flight_angle(l) <= 180.0_wp )  THEN
                IF ( y_pos(l) >= ( ny + 0.5_wp ) * dy )                                            &
                     y_pos(l) = y_pos(l) - ( ny + 1 ) * dy
!
!--          Check if aircraft leaves the model domain at the south boundary
             ELSEIF ( flight_angle(l) > 180.0_wp  .AND.  flight_angle(l) < 360.0_wp )  THEN
                IF ( y_pos(l) < -0.5_wp * dy )                                                     &
                     y_pos(l) = ( ny + 1 ) * dy + y_pos(l)
             ENDIF

          ENDIF
!
!--       Check if maximum elevation change is already reached. If required reflect vertically.
!--       One have to distinguish between a flight that starts with a descent or a climb.
          IF ( rate_of_climb(l) > 0.0_wp )  THEN
!
!--          First check if aircraft is too high
             IF (  w_agl(l) > 0.0_wp  .AND.  z_pos(l) - flight_level(l) > max_elev_change(l) )  THEN
                z_pos(l) = 2.0_wp * ( flight_level(l) + max_elev_change(l) ) - z_pos(l)
                w_agl(l) = - w_agl(l)
!
!--          Check if aircraft is too low
             ELSEIF (  w_agl(l) < 0.0_wp  .AND.  z_pos(l) < flight_level(l) )  THEN
                z_pos(l) = 2.0_wp * flight_level(l) - z_pos(l)
                w_agl(l) = - w_agl(l)
             ENDIF

          ELSEIF ( rate_of_climb(l) < 0.0_wp )  THEN
!
!--          First check if aircraft is too high
             IF (  w_agl(l) > 0.0_wp  .AND.  z_pos(l) > flight_level(l) )  THEN
                z_pos(l) = 2.0_wp * flight_level(l) - z_pos(l)
                w_agl(l) = - w_agl(l)
!
!--          Check if aircraft is too low
             ELSEIF (  w_agl(l) < 0.0_wp  .AND.  flight_level(l) - z_pos(l) > max_elev_change(l) )  THEN
                z_pos(l) = 2.0_wp * ( flight_level(l) - max_elev_change(l) ) - z_pos(l)
                w_agl(l) = - w_agl(l)

             ENDIF
          ENDIF
!
!--       Determine grid indices for flight position along x- and y-direction. Please note, there is
!--       a special treatment for the index along z-direction, which is due to vertical grid stretching.
          i = x_pos(l) * ddx
          j = y_pos(l) * ddy
!
!--       Check if indices are on current PE
          on_pe = ( i >= nxl  .AND.  i <= nxr  .AND.  j >= nys  .AND.  j <= nyn )

          IF ( on_pe )  THEN

             var_index = 1
!
!--          Recalculate indices for u on the xu-y grid. Indicies are calculated so that no
!--          case differentiation need to be made but variables only need to be interpolate
!--          between indices i and i+1, j and j+1, as well as k and k+1.
             i =   x_pos(l) * ddx
             j = ( y_pos(l) - 0.5_wp * dy ) * ddy
!
!--          Interpolate u-component onto current flight position.
             CALL interpolate_xyz( u, zu, ddzu, var_index, j, i )
             var_index = var_index + 1
!
!--          Recalculate indices for v on the x-yv grid.
             i = ( x_pos(l) - 0.5_wp * dx ) * ddx
             j =   y_pos(l) * ddy
!
!--          Interpolate v-component onto current flight position.
             CALL interpolate_xyz( v, zu, ddzu, var_index, j, i )
             var_index = var_index + 1
!
!--          Interpolate w and scalar quantities. Recalculate indices on the x-y grid.
             i  = ( x_pos(l) - 0.5_wp * dx ) * ddx
             j  = ( y_pos(l) - 0.5_wp * dy ) * ddy
!
!--          Interpolate w-velocity component.
             CALL interpolate_xyz( w, zw, ddzw, var_index, j, i )
             var_index = var_index + 1
!
!--          Potential temerature
             IF ( .NOT. neutral )  THEN
                CALL interpolate_xyz( pt, zu, ddzu, var_index, j, i )
                var_index = var_index + 1
             ENDIF
!
!--          Humidity
             IF ( humidity )  THEN
                CALL interpolate_xyz( q, zu, ddzu, var_index, j, i )
                var_index = var_index + 1
             ENDIF
!
!--          Liquid water content
             IF ( bulk_cloud_model  .OR.  cloud_droplets )  THEN
                CALL interpolate_xyz( ql, zu, ddzu, var_index, j, i )
                var_index = var_index + 1
             ENDIF
!
!--          Passive scalar
             IF ( passive_scalar )  THEN
                CALL interpolate_xyz( s, zu, ddzu, var_index, j, i )
                var_index = var_index + 1
             ENDIF
!
!--          Treat user-defined variables if required
             DO  n = 1, num_var_fl_user
                CALL user_flight( var_u, n )
                CALL interpolate_xyz( var_u, zu, ddzu, var_index, j, i )
                var_index = var_index + 1
             ENDDO
          ENDIF

       ENDDO
!
!--    Write local data on global array.
#if defined( __parallel )
       CALL MPI_ALLREDUCE( sensor_l(1,1), sensor(1,1), num_var_fl * num_leg, MPI_REAL, MPI_SUM,    &
                           comm2d, ierr )
#else
       sensor     = sensor_l
#endif
    ENDIF

    CALL cpu_log( log_point(65), 'flight_measurement', 'stop' )

 END SUBROUTINE flight_measurement

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine bi-linearly interpolates the respective data onto the current flight position.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE interpolate_xyz( var, z_uw, ddz_uw, var_ind, j, i )

    USE grid_variables,                                                                            &
       ONLY:  ddx,                                                                                 &
              ddy,                                                                                 &
              dx,                                                                                  &
              dy

    USE indices,                                                                                   &
        ONLY:  nzb,                                                                                &
               nzt,                                                                                &
               nxlg,                                                                               &
               nxrg,                                                                               &
               nysg,                                                                               &
               nyng

    IMPLICIT NONE

    INTEGER(iwp) ::  i        !< index along x
    INTEGER(iwp) ::  j        !< index along y
    INTEGER(iwp) ::  k        !< index along z
    INTEGER(iwp) ::  k1       !< dummy variable
    INTEGER(iwp) ::  var_ind  !< index variable for output quantity

    REAL(wp) ::  var_int      !< interpolated variable at current probe position
    REAL(wp) ::  var_int_l    !< horizontally interpolated variable at k-level
    REAL(wp) ::  var_int_ly1  !< horizontally interpolated variable at k-level at index j
    REAL(wp) ::  var_int_ly2  !< horizontally interpolated variable at k-level at index j+1
    REAL(wp) ::  var_int_u    !< horizontally interpolated variable at (k+1)-level
    REAL(wp) ::  var_int_uy1  !< horizontally interpolated variable at (k+1)-level at index j
    REAL(wp) ::  var_int_uy2  !< horizontally interpolated variable at (k+1)-level at index j+1

    REAL(wp), DIMENSION(1:nzt+1)   ::  ddz_uw  !< inverse vertical grid spacing
    REAL(wp), DIMENSION(nzb:nzt+1) ::  z_uw    !< height level

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  var  !< treated quantity
!
!-- Obtain vertical index by searching. This is required due to the vertical grid stretching
!-- as well as due to the asymmetric vertical grid at the lowest grid level between k = 0 and k = 1.
    DO  k1 = nzb, nzt
       IF ( z_pos(l) >= z_uw(k1) .AND. z_pos(l) < z_uw(k1+1) )  THEN
          k = k1
          EXIT
       ENDIF
    ENDDO
!
!-- Bi-linearly interpolate the required variable onto its x-y sensor position at discrete levels
!-- k and (k+1). Therefore, first interpolate the variable along x at its discrete y-locations
!-- j and (j+1). In a second step interpolate onto the x-y sensor position along y.
    var_int_ly1 = var(k,j,i)   + ( var(k,j,i+1)   - var(k,j,i)   ) * ddx * ( x_pos(l) - i * dx )
    var_int_ly2 = var(k,j+1,i) + ( var(k,j+1,i+1) - var(k,j+1,i) ) * ddx * ( x_pos(l) - i * dx )
    var_int_l   = var_int_ly1 + ( var_int_ly2 - var_int_ly1 ) * ddy * ( y_pos(l) - j * dy )

    var_int_uy1 = var(k+1,j,i)   + ( var(k+1,j,i+1)   - var(k+1,j,i)   ) * ddx * ( x_pos(l) - i * dx )
    var_int_uy2 = var(k+1,j+1,i) + ( var(k+1,j+1,i+1) - var(k+1,j+1,i) ) * ddx * ( x_pos(l) - i * dx )
    var_int_u   = var_int_uy1 + ( var_int_uy2 - var_int_uy1 ) * ddy * ( y_pos(l) - j * dy )
!
!-- Now interpolate linearly onto the exact z-position.
    var_int = var_int_l + (var_int_u - var_int_l ) * ddz_uw(k+1) * ( z_pos(l) - z_uw(k) )
!
!-- Store on local data array
    sensor_l(var_ind,l) = var_int

 END SUBROUTINE interpolate_xyz

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform parameter checks.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE flight_check_parameters

    USE arrays_3d,                                                                                 &
        ONLY:  zu

    USE control_parameters,                                                                        &
        ONLY:  bc_lr_cyc,                                                                          &
               bc_ns_cyc,                                                                          &
               message_string

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               ny,                                                                                 &
               nz

    IMPLICIT NONE

!
!-- Check if start positions are properly set.
    DO  l = 1, num_leg
       IF ( x_start(l) < 0.0_wp  .OR.  x_start(l) > ( nx + 1 ) * dx )  THEN
          message_string = 'start x position is outside the model domain'
          CALL message( 'flight_check_parameters', 'VFL0001', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( y_start(l) < 0.0_wp  .OR.  y_start(l) > ( ny + 1 ) * dy )  THEN
          message_string = 'start y position is outside the model domain'
          CALL message( 'flight_check_parameters', 'VFL0001', 1, 2, 0, 6, 0 )
       ENDIF

    ENDDO
!
!-- Check for leg mode
    DO  l = 1, num_leg
!
!--    Check if leg mode matches the overall lateral model boundary conditions.
       IF ( TRIM( leg_mode(l) ) == 'cyclic' )  THEN
          IF ( .NOT. bc_lr_cyc  .OR.  .NOT. bc_ns_cyc )  THEN
             message_string = 'cyclic flight leg does not match lateral boundary condition'
             CALL message( 'flight_check_parameters', 'VFL0002', 1, 2, 0, 6, 0 )
          ENDIF
!
!--    Check if end-positions are inside the model domain in case of return-legs.
       ELSEIF ( TRIM( leg_mode(l) ) == 'return' )  THEN
          IF ( x_end(l) > ( nx + 1 ) * dx  .OR.  y_end(l) > ( ny + 1 ) * dx )  THEN
             message_string = 'flight leg or parts of it are outside the model domain'
             CALL message( 'flight_check_parameters', 'VFL0003', 1, 2, 0, 6, 0 )
          ENDIF
       ELSE
          message_string = 'unknown flight mode "' // TRIM( leg_mode(l) ) // '"'
          CALL message( 'flight_check_parameters', 'VFL0004', 1, 2, 0, 6, 0 )
       ENDIF

    ENDDO
!
!-- Check if given flight object remains inside model domain if a rate of climb / descent is
!-- prescribed.
    DO  l = 1, num_leg
       IF ( flight_level(l) + max_elev_change(l) > zu(nz) .AND. rate_of_climb(l) > 0.0_wp .OR.     &
            flight_level(l) - max_elev_change(l) <= 0.0_wp .AND. rate_of_climb(l) < 0.0_wp )  THEN
          message_string = 'flight level is outside the model domain '
          CALL message( 'flight_check_parameters', 'VFL0005', 1, 2, 0, 6, 0 )
       ENDIF
    ENDDO


 END SUBROUTINE flight_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE flight_rrd_global_ftn( found )


    USE control_parameters,                                                                        &
        ONLY: length,                                                                              &
              restart_string


    IMPLICIT NONE

    LOGICAL, INTENT(OUT)  ::  found  !< flag indicating if a variable string is found


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'u_agl' )
          IF ( .NOT. ALLOCATED( u_agl ) )  ALLOCATE( u_agl(1:num_leg) )
          READ ( 13 )  u_agl
       CASE ( 'v_agl' )
          IF ( .NOT. ALLOCATED( v_agl ) )  ALLOCATE( v_agl(1:num_leg) )
          READ ( 13 )  v_agl
       CASE ( 'w_agl' )
          IF ( .NOT. ALLOCATED( w_agl ) )  ALLOCATE( w_agl(1:num_leg) )
          READ ( 13 )  w_agl
       CASE ( 'x_pos' )
          IF ( .NOT. ALLOCATED( x_pos ) )  ALLOCATE( x_pos(1:num_leg) )
          READ ( 13 )  x_pos
       CASE ( 'y_pos' )
          IF ( .NOT. ALLOCATED( y_pos ) )  ALLOCATE( y_pos(1:num_leg) )
          READ ( 13 )  y_pos
       CASE ( 'z_pos' )
          IF ( .NOT. ALLOCATED( z_pos ) )  ALLOCATE( z_pos(1:num_leg) )
          READ ( 13 )  z_pos

       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE flight_rrd_global_ftn



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!------------------------------------------------------------------------------!
 SUBROUTINE flight_rrd_global_mpi


    IMPLICIT NONE

    LOGICAL  ::  array_found  !< flag indicating if respective array is found in restart file


    CALL rd_mpi_io_check_array( 'u_agl', found = array_found )
    IF ( array_found)  THEN
       IF ( .NOT. ALLOCATED( u_agl ) )  ALLOCATE( u_agl(1:num_leg) )
       CALL rrd_mpi_io_global_array( 'u_agl', u_agl )
    ENDIF
    CALL rd_mpi_io_check_array( 'v_agl', found = array_found )
    IF ( array_found)  THEN
       IF ( .NOT. ALLOCATED( v_agl ) )  ALLOCATE( v_agl(1:num_leg) )
       CALL rrd_mpi_io_global_array( 'v_agl', v_agl )
    ENDIF
    CALL rd_mpi_io_check_array( 'w_agl', found = array_found )
    IF ( array_found)  THEN
       IF ( .NOT. ALLOCATED( w_agl ) )  ALLOCATE( w_agl(1:num_leg) )
       CALL rrd_mpi_io_global_array( 'w_agl', w_agl )
    ENDIF
    CALL rd_mpi_io_check_array( 'x_pos', found = array_found )
    IF ( array_found)  THEN
       IF ( .NOT. ALLOCATED( x_pos ) )  ALLOCATE( x_pos(1:num_leg) )
       CALL rrd_mpi_io_global_array( 'x_pos', x_pos )
    ENDIF
    CALL rd_mpi_io_check_array( 'y_pos', found = array_found )
    IF ( array_found)  THEN
       IF ( .NOT. ALLOCATED( y_pos ) )  ALLOCATE( y_pos(1:num_leg) )
       CALL rrd_mpi_io_global_array( 'y_pos', y_pos )
    ENDIF
    CALL rd_mpi_io_check_array( 'z_pos', found = array_found )
    IF ( array_found)  THEN
       IF ( .NOT. ALLOCATED( z_pos ) )  ALLOCATE( z_pos(1:num_leg) )
       CALL rrd_mpi_io_global_array( 'z_pos', z_pos )
    ENDIF

 END SUBROUTINE flight_rrd_global_mpi


    !--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE flight_wrd_global


    IMPLICIT NONE


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL wrd_write_string( 'u_agl' )
       WRITE ( 14 )  u_agl

       CALL wrd_write_string( 'v_agl' )
       WRITE ( 14 )  v_agl

       CALL wrd_write_string( 'w_agl' )
       WRITE ( 14 )  w_agl

       CALL wrd_write_string( 'x_pos' )
       WRITE ( 14 )  x_pos

       CALL wrd_write_string( 'y_pos' )
       WRITE ( 14 )  y_pos

       CALL wrd_write_string( 'z_pos' )
       WRITE ( 14 )  z_pos

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       CALL wrd_mpi_io_global_array( 'u_agl', u_agl )
       CALL wrd_mpi_io_global_array( 'v_agl', v_agl )
       CALL wrd_mpi_io_global_array( 'w_agl', w_agl )
       CALL wrd_mpi_io_global_array( 'x_pos', x_pos )
       CALL wrd_mpi_io_global_array( 'y_pos', y_pos )
       CALL wrd_mpi_io_global_array( 'z_pos', z_pos )

    ENDIF

 END SUBROUTINE flight_wrd_global


 END MODULE flight_mod
