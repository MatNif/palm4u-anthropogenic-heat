!> @file fastv8_updata.f90
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
! Copyright 2017-2023 Carl von Ossietzky Universität Oldenburg
! Copyright 2022-2023 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Additional routines used in C Codebase of FASTv8 Coupler
!--------------------------------------------------------------------------------------------------!
 MODULE fastv8_updata

    USE fastv8_coupler_mod
    USE kinds
   
    IMPLICIT NONE

    PRIVATE

!
!- Public functions
    PUBLIC                                                                                         &
       palm_bld_data,                                                                              &
       palm_sim_env,                                                                               &
       palm_turb_par,                                                                              &
       palm_force_table_bld,                                                                       &
       palm_position_table,                                                                        &
       palm_vel_value


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates PALM simulation/environment variables based on FAST information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_sim_env(turbine_id, dt_fast_in, shaft_height_fast_in)                             &
            BIND(C,NAME="palm_sim_env")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID

    REAL(kind=C_DOUBLE), INTENT(IN) ::  dt_fast_in            !< time steps size FAST
    REAL(kind=C_DOUBLE), INTENT(IN) ::  shaft_height_fast_in  !< shaft height from FAST

    dt_fast = dt_fast_in
    shaft_height_fast(turbine_id) = shaft_height_fast_in

 END SUBROUTINE palm_sim_env


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates PALM turbine parameters based on FAST information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_turb_par(turbine_id, simtimef, vel_rot, xs_x, xs_y, xs_z)                         &
            BIND(C,NAME="palm_turb_par")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID

    REAL(kind=C_DOUBLE), INTENT(IN) ::  simtimef  !< current FAST simulation time
    REAL(kind=C_DOUBLE), INTENT(IN) ::  vel_rot   !< rotational velocity from FAST
    REAL(kind=C_DOUBLE), INTENT(IN) ::  xs_x      !< shaft coordinate in X direction
    REAL(kind=C_DOUBLE), INTENT(IN) ::  xs_y      !< shaft coordinate in Y direction
    REAL(kind=C_DOUBLE), INTENT(IN) ::  xs_z      !< shaft coordinate in Z direction
   
    rotspeed(turbine_id) = vel_rot
    current_time_fast(turbine_id) = simtimef

    shaft_coordinates(turbine_id,1) = xs_x
    shaft_coordinates(turbine_id,2) = xs_y
    shaft_coordinates(turbine_id,3) = xs_z

 END SUBROUTINE palm_turb_par


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write the FAST blade element and hub positions into the PALM lookup table.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE palm_position_table(turbine_id, max_blade_elem_id, blade_elem_id,            &
                                          pos_x, pos_y, pos_z, error_status)                       &
                      BIND(C,NAME="palm_position_table")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id      !< turbine ID
    INTEGER(kind=C_INT), INTENT(IN) ::  max_blade_elem_id  !<
    INTEGER(kind=C_INT), INTENT(IN) ::  blade_elem_id     !<

    REAL(kind=C_DOUBLE), INTENT(IN) ::  pos_x  !< input information: positions X
    REAL(kind=C_DOUBLE), INTENT(IN) ::  pos_y  !< input information: positions Y
    REAL(kind=C_DOUBLE), INTENT(IN) ::  pos_z  !< input information: positions Z

    INTEGER(kind=C_INT), INTENT(OUT) :: error_status  !< determines if an error has occured

    INTEGER ::  cflag
    INTEGER ::  fast_blade_elem_idx
 
    cflag = 0
    error_status = 0
!
!-- check which blade the current blade element is
    IF ( blade_elem_id <= (max_blade_elem_id-1) / 3 )  THEN
       cflag = 1
    ELSEIF ( ( blade_elem_id > (max_blade_elem_id-1) / 3 )  .AND.                                  &
             ( blade_elem_id <= (max_blade_elem_id-1) / 3 * 2 ) )  THEN
       cflag = 2
    ELSEIF ( ( blade_elem_id > (max_blade_elem_id-1) / 3 * 2 )  .AND.                              &
             ( blade_elem_id <= (max_blade_elem_id-1) ) )  THEN
       cflag = 3
!
!-- the last element is the hub position
    ELSEIF ( blade_elem_id == max_blade_elem_id )  THEN
       cflag = 4
    ELSE
       error_status = 1
    END IF
!
!-- update blade element positions
    IF ( (cflag >= 1)  .AND.  (cflag <= 3) )  THEN

       fast_blade_elem_idx = blade_elem_id - (cflag-1)*((max_blade_elem_id-1)/3)

       fast_position_blade(fast_blade_elem_idx, cflag, turbine_id, 1) = pos_x
       fast_position_blade(fast_blade_elem_idx, cflag, turbine_id, 2) = pos_y
       fast_position_blade(fast_blade_elem_idx, cflag, turbine_id, 3) = pos_z
!
!-- update hub position
    ELSEIF ( cflag == 4 )  THEN
       fast_hub_center_pos(turbine_id, 1) = pos_x
       fast_hub_center_pos(turbine_id, 2) = pos_y
       fast_hub_center_pos(turbine_id, 3) = pos_z

    ELSE
       error_status = 1
    END IF

 END SUBROUTINE palm_position_table


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write the FAST blade element forces into the PALM lookup table.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE palm_force_table_bld(turbine_id, max_blade_elem_id, blade_elem_id,           &
                                           force_x, force_y, force_z, error_status)                &
                      BIND(C,NAME="palm_force_table_bld")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id      !< turbine ID
    INTEGER(kind=C_INT), INTENT(IN) ::  max_blade_elem_id  !< 
    INTEGER(kind=C_INT), INTENT(IN) ::  blade_elem_id     !< blade element ID

    REAL(kind=C_DOUBLE), INTENT(IN) ::  force_x  !< input information: force in X direction
    REAL(kind=C_DOUBLE), INTENT(IN) ::  force_y  !< input information: force in Y direction
    REAL(kind=C_DOUBLE), INTENT(IN) ::  force_z  !< input information: force in Z direction

    INTEGER(kind=C_INT), INTENT(OUT)  :: error_status  !< determines if an error has occured

    INTEGER ::  cflag
    INTEGER ::  fast_blade_elem_idx
 
    cflag = 0
    error_status = 0
!
!-- check which blade the current blade element is
    IF ( blade_elem_id <= (max_blade_elem_id-1) / 3 )  THEN
       cflag = 1
    ELSEIF ( ( blade_elem_id >  (max_blade_elem_id-1) / 3 )  .AND.                                 &
             ( blade_elem_id <= (max_blade_elem_id-1) / 3 * 2 ) )  THEN
       cflag = 2
    ELSEIF ( ( blade_elem_id >  (max_blade_elem_id-1) / 3 * 2 )  .AND.                             &
             ( blade_elem_id <= (max_blade_elem_id-1) ) )  THEN
       cflag = 3
    ELSE
       error_status = 1
    END IF
!
!-- update blade element forces
    IF ( (cflag >= 1)  .AND.  (cflag <= 3) )  THEN

       fast_blade_elem_idx = blade_elem_id - (cflag-1)*((max_blade_elem_id-1)/3)

       fast_force_blade(fast_blade_elem_idx, cflag, turbine_id, 1) = force_x
       fast_force_blade(fast_blade_elem_idx, cflag, turbine_id, 2) = force_y
       fast_force_blade(fast_blade_elem_idx, cflag, turbine_id, 3) = force_z

    ELSE
       error_status = 1
    END IF

 END SUBROUTINE palm_force_table_bld


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write the FAST blade element data into the PALM lookup table.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE palm_bld_data(turbine_id, numbld, numbldelem, rotrad)                                  &
            BIND(C,NAME="palm_bld_data")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN) ::  turbine_id  !< turbine ID
    INTEGER(kind=C_INT), INTENT(IN) ::  numbld      !< number of blades
    INTEGER(kind=C_INT), INTENT(IN) ::  numbldelem  !< number of blade elements

    REAL(kind=C_DOUBLE), INTENT(IN) ::  rotrad      !< rotor radius (length of blades)

    fast_n_blades(turbine_id) =     numbld
    fast_n_blade_elem(turbine_id) = numbldelem
    fast_n_radius(turbine_id) =     rotrad

 END SUBROUTINE palm_bld_data


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads the velocities for fast out of the PALM data structures.
!--------------------------------------------------------------------------------------------------!
 RECURSIVE SUBROUTINE palm_vel_value(turbine_id, bladeid, elemid,                                  &
                                     compu, compv, compw, error_status)                            &
                      BIND(C,NAME="palm_vel_value")

    USE ISO_C_BINDING,                                                                             &
        ONLY:  C_INT,                                                                              &
               C_DOUBLE

    INTEGER(kind=C_INT), INTENT(IN)  ::  turbine_id  !< turbine ID
    INTEGER(kind=C_INT), INTENT(IN)  ::  bladeid     !< blade ID (index)
    INTEGER(kind=C_INT), INTENT(IN)  ::  elemid      !< blade element ID (index)

    REAL(kind=C_DOUBLE), INTENT(OUT) ::  compu  !< wind velocity component in X direction
    REAL(kind=C_DOUBLE), INTENT(OUT) ::  compv  !< wind velocity component in Y direction
    REAL(kind=C_DOUBLE), INTENT(OUT) ::  compw  !< wind velocity component in Z direction

    INTEGER(kind=C_INT), INTENT(OUT) ::  error_status  !< determines if an error has been encountered

    error_status = 0

    IF ( (turbine_id > 0)  .AND.  (turbine_id <= nturbines) )  THEN
      
       IF ( (bladeid > 0)  .AND.  (bladeid <= fast_n_blades(turbine_id)) )  THEN

          IF ( (elemid > 0)  .AND.  (elemid <= fast_n_blade_elem(turbine_id)) )  THEN 
             compu = palm_vel_blade(elemid, bladeid, turbine_id, 1)
             compv = palm_vel_blade(elemid, bladeid, turbine_id, 2)
             compw = palm_vel_blade(elemid, bladeid, turbine_id, 3)
          ELSE
             error_status = 1
          ENDIF
          
       ELSEIF ( (bladeid == 0)  .AND.  (elemid == 0) )  THEN
          compu = palm_hub_center_vel(turbine_id, 1)
          compv = palm_hub_center_vel(turbine_id, 2)
          compw = palm_hub_center_vel(turbine_id, 3)
       ELSE
          error_status = 1
       ENDIF

    ELSE
       error_status = 1
    ENDIF

 END SUBROUTINE palm_vel_value

 END MODULE fastv8_updata
