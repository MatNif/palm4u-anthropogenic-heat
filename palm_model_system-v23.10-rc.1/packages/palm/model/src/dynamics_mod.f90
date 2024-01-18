!> @file dynamics_mod.f90
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
!> This module contains the dynamics of PALM.
!--------------------------------------------------------------------------------------------------!
 MODULE dynamics_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  diss,                                                                               &
               diss_p,                                                                             &
               dzu,                                                                                &
               e,                                                                                  &
               e_p,                                                                                &
               exner,                                                                              &
               hyp,                                                                                &
               pt,                                                                                 &
               pt_1,                                                                               &
               pt_2,                                                                               &
               pt_init,                                                                            &
               pt_p,                                                                               &
               q, q_1,                                                                             &
               q_2,                                                                                &
               q_p,                                                                                &
               rho_air,                                                                            &
               s,                                                                                  &
               s_1,                                                                                &
               s_2,                                                                                &
               s_p,                                                                                &
               u,                                                                                  &
               u_1,                                                                                &
               u_2,                                                                                &
               u_init,                                                                             &
               u_p,                                                                                &
               v,                                                                                  &
               v_1,                                                                                &
               v_2,                                                                                &
               v_p,                                                                                &
               v_init,                                                                             &
               w,                                                                                  &
               w_1,                                                                                &
               w_2,                                                                                &
               w_p,                                                                                &
               zu

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  magnus,                                                                             &
               rd_d_rv

    USE control_parameters,                                                                        &
        ONLY:  bc_dirichlet_l,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               bc_pt_b,                                                                            &
               bc_pt_t_val,                                                                        &
               bc_q_t_val,                                                                         &
               bc_s_t_val,                                                                         &
               check_realistic_q,                                                                  &
               child_domain,                                                                       &
               constant_diffusion,                                                                 &
               cyclic_fill_initialization,                                                         &
               dt_3d,                                                                              &
               homogenize_surface_temperature,                                                     &
               humidity,                                                                           &
               ibc_pt_b,                                                                           &
               ibc_pt_t,                                                                           &
               ibc_q_b,                                                                            &
               ibc_q_t,                                                                            &
               ibc_s_b,                                                                            &
               ibc_s_t,                                                                            &
               ibc_uv_b,                                                                           &
               ibc_uv_t,                                                                           &
               initializing_actions,                                                               &
               intermediate_timestep_count,                                                        &
               land_surface,                                                                       &
               length,                                                                             &
               message_string,                                                                     &
               nested_run,                                                                         &
               nesting_offline,                                                                    &
               neutral,                                                                            &
               nudging,                                                                            &
               passive_scalar,                                                                     &
               pt_surface_heating_rate,                                                            &
               pt_surface_initial_change,                                                          &
               q_surface_initial_change,                                                           &
               restart_string,                                                                     &
               rans_mode,                                                                          &
               rans_tke_e,                                                                         &
               s_surface_initial_change,                                                           &
               time_since_reference_point,                                                         &
               tsc,                                                                                &
               urban_surface

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz


    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nys,                                                                                &
               nysg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nzb,                                                                                &
               nzt

    USE kinds

    USE pegrid

    USE pmc_interface,                                                                             &
        ONLY:  root_model

    USE surface_mod,                                                                               &
        ONLY:  bc_hv,                                                                              &
               surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_usm

    USE statistics,                                                                                &
        ONLY:  weight_substep


    IMPLICIT NONE

    LOGICAL ::  dynamics_module_enabled = .FALSE.   !<

    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC                                                                                         &
       dynamics_parin,                                                                             &
       dynamics_check_parameters,                                                                  &
       dynamics_check_data_output_ts,                                                              &
       dynamics_check_data_output_pr,                                                              &
       dynamics_check_data_output,                                                                 &
       dynamics_init_masks,                                                                        &
       dynamics_define_netcdf_grid,                                                                &
       dynamics_init_arrays,                                                                       &
       dynamics_init,                                                                              &
       dynamics_init_checks,                                                                       &
       dynamics_header,                                                                            &
       dynamics_actions,                                                                           &
       dynamics_non_advective_processes,                                                           &
       dynamics_exchange_horiz,                                                                    &
       dynamics_prognostic_equations,                                                              &
       dynamics_boundary_conditions,                                                               &
       dynamics_swap_timelevel,                                                                    &
       dynamics_3d_data_averaging,                                                                 &
       dynamics_data_output_2d,                                                                    &
       dynamics_data_output_3d,                                                                    &
       dynamics_statistics,                                                                        &
       dynamics_rrd_global,                                                                        &
       dynamics_rrd_local,                                                                         &
       dynamics_wrd_global,                                                                        &
       dynamics_wrd_local,                                                                         &
       dynamics_last_actions

!
!-- Public parameters, constants and initial values
    PUBLIC                                                                                         &
       dynamics_module_enabled

    INTERFACE dynamics_parin
       MODULE PROCEDURE dynamics_parin
    END INTERFACE dynamics_parin

    INTERFACE dynamics_check_parameters
       MODULE PROCEDURE dynamics_check_parameters
    END INTERFACE dynamics_check_parameters

    INTERFACE dynamics_check_data_output_ts
       MODULE PROCEDURE dynamics_check_data_output_ts
    END INTERFACE dynamics_check_data_output_ts

    INTERFACE dynamics_check_data_output_pr
       MODULE PROCEDURE dynamics_check_data_output_pr
    END INTERFACE dynamics_check_data_output_pr

    INTERFACE dynamics_check_data_output
       MODULE PROCEDURE dynamics_check_data_output
    END INTERFACE dynamics_check_data_output

    INTERFACE dynamics_init_masks
       MODULE PROCEDURE dynamics_init_masks
    END INTERFACE dynamics_init_masks

    INTERFACE dynamics_define_netcdf_grid
       MODULE PROCEDURE dynamics_define_netcdf_grid
    END INTERFACE dynamics_define_netcdf_grid

    INTERFACE dynamics_init_arrays
       MODULE PROCEDURE dynamics_init_arrays
    END INTERFACE dynamics_init_arrays

    INTERFACE dynamics_init
       MODULE PROCEDURE dynamics_init
    END INTERFACE dynamics_init

    INTERFACE dynamics_init_checks
       MODULE PROCEDURE dynamics_init_checks
    END INTERFACE dynamics_init_checks

    INTERFACE dynamics_header
       MODULE PROCEDURE dynamics_header
    END INTERFACE dynamics_header

    INTERFACE dynamics_actions
       MODULE PROCEDURE dynamics_actions
       MODULE PROCEDURE dynamics_actions_ij
    END INTERFACE dynamics_actions

    INTERFACE dynamics_non_advective_processes
       MODULE PROCEDURE dynamics_non_advective_processes
       MODULE PROCEDURE dynamics_non_advective_processes_ij
    END INTERFACE dynamics_non_advective_processes

    INTERFACE dynamics_exchange_horiz
       MODULE PROCEDURE dynamics_exchange_horiz
    END INTERFACE dynamics_exchange_horiz

    INTERFACE dynamics_prognostic_equations
       MODULE PROCEDURE dynamics_prognostic_equations
       MODULE PROCEDURE dynamics_prognostic_equations_ij
    END INTERFACE dynamics_prognostic_equations

    INTERFACE dynamics_boundary_conditions
       MODULE PROCEDURE dynamics_boundary_conditions
    END INTERFACE dynamics_boundary_conditions

    INTERFACE dynamics_swap_timelevel
       MODULE PROCEDURE dynamics_swap_timelevel
    END INTERFACE dynamics_swap_timelevel

    INTERFACE dynamics_3d_data_averaging
       MODULE PROCEDURE dynamics_3d_data_averaging
    END INTERFACE dynamics_3d_data_averaging

    INTERFACE dynamics_data_output_2d
       MODULE PROCEDURE dynamics_data_output_2d
    END INTERFACE dynamics_data_output_2d

    INTERFACE dynamics_data_output_3d
       MODULE PROCEDURE dynamics_data_output_3d
    END INTERFACE dynamics_data_output_3d

    INTERFACE dynamics_statistics
       MODULE PROCEDURE dynamics_statistics
    END INTERFACE dynamics_statistics

    INTERFACE dynamics_rrd_global
       MODULE PROCEDURE dynamics_rrd_global_ftn
       MODULE PROCEDURE dynamics_rrd_global_mpi
    END INTERFACE dynamics_rrd_global

    INTERFACE dynamics_rrd_local
       MODULE PROCEDURE dynamics_rrd_local_ftn
       MODULE PROCEDURE dynamics_rrd_local_mpi
    END INTERFACE dynamics_rrd_local

    INTERFACE dynamics_wrd_global
       MODULE PROCEDURE dynamics_wrd_global
    END INTERFACE dynamics_wrd_global

    INTERFACE dynamics_wrd_local
       MODULE PROCEDURE dynamics_wrd_local
    END INTERFACE dynamics_wrd_local

    INTERFACE dynamics_last_actions
       MODULE PROCEDURE dynamics_last_actions
    END INTERFACE dynamics_last_actions


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific namelist
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_parin

    CHARACTER(LEN=100)  ::  line  !< dummy string that contains the current line of the parameter
                                  !< file
    INTEGER(iwp)  ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /dynamics_parameters/  switch_off_module


!
!-- For the time beeing (unless the dynamics module is further developed), set default module
!-- switch to true.
    dynamics_module_enabled = .TRUE.

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, dynamics_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    dynamics_parameters namelist was found and read correctly.
       IF ( .NOT. switch_off_module )  dynamics_module_enabled = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    dynamics_parameters namelist was found, but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'dynamics_parameters', line )

    ENDIF

 END SUBROUTINE dynamics_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check control parameters and deduce further quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_check_parameters

!
!-- Checks concering homogenization of surface temperature.
    IF ( homogenize_surface_temperature )  THEN
!
!--    Homogenization of the surface temperature is only allowed if the temperature equation is
!--    enabled and if a Dirichlet boundary conditions for the pot. temperature is set.
       IF ( neutral  .OR.  bc_pt_b /= 'dirichlet' )  THEN
          message_string = 'homogenize_surface_temperature = .T. is only allowed in&' //           &
                           'combination with neutral = .F. and bc_pt_b = "dirichlet".'
          CALL message( 'dynamics_check_parameters', 'PAC0197', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Homogenization of surface temperature is not allowed in combination with the land- and urban-
!--    surface model, where the soil- and surface temperature is determined by a prognostic equation.
!--    Note, this covers also the situation where the soil/surface temperature can be input via
!--    the dynamic input file.
       IF ( land_surface  .OR.  urban_surface )  THEN
          message_string = 'homogenize_surface_temperature = .T. is not allowed in combination ' //&
                           'with the land- or urban-surface model.'
          CALL message( 'dynamics_check_parameters', 'PAC0198', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

 END SUBROUTINE dynamics_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set module-specific timeseries units and labels
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )

    INTEGER(iwp),      INTENT(IN)     ::  dots_max

    CHARACTER (LEN=*), DIMENSION(dots_max), INTENT(INOUT)  :: dots_label
    CHARACTER (LEN=*), DIMENSION(dots_max), INTENT(INOUT)  :: dots_unit

    INTEGER(iwp),      INTENT(INOUT)  ::  dots_num

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( dots_num == 0  .OR.  dots_label(1)(1:1) == ' '  .OR.  dots_unit(1)(1:1) == ' ' )  CONTINUE


 END SUBROUTINE dynamics_check_data_output_ts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of module-specific profile output quantities. For those variables not recognized,
!> the parameter unit is set to "illegal", which tells the calling routine that the output variable
!> is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_check_data_output_pr( variable, var_count, unit, dopr_unit )


    CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit
    CHARACTER (LEN=*) ::  unit     !<
    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  var_count     !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( unit(1:1) == ' '  .OR.  dopr_unit(1:1) == ' '  .OR.  var_count == 0 )  CONTINUE

    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'var_name' )

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE dynamics_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of module-specific output quantities. For those variables not recognized,
!> the parameter unit is set to "illegal", which tells the calling routine that the output variable
!< is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_check_data_output( variable, unit )


    CHARACTER (LEN=*) ::  unit     !<
    CHARACTER (LEN=*) ::  variable !<

    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2' )

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE dynamics_check_data_output


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Initialize module-specific masked output
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_init_masks( variable, unit )


    CHARACTER (LEN=*) ::  unit     !<
    CHARACTER (LEN=*) ::  variable !<


    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2' )

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE dynamics_init_masks


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize module-specific arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_init_arrays


 END SUBROUTINE dynamics_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of module-specific initializing actions
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_init

    INTEGER(iwp) ::  i !< grid index in x direction
    INTEGER(iwp) ::  j !< grid index in y direction
    INTEGER(iwp) ::  k !< grid index in z direction
    INTEGER(iwp) ::  m !< running index over boundary grid points

    REAL(wp) ::  pt_surf_mean !< mean surface temperature

    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.                                &
         .NOT. cyclic_fill_initialization )  THEN
!
!--    If required, change the surface temperature at the start of the 3D run. Further, directly
!--    write the surface temperature onto the corresponding surface attribute %pt_surface, for
!--    sake of data output before the very first timestep will be done.
       IF ( pt_surface_initial_change /= 0.0_wp )  THEN
          DO  m = 1, bc_hv%ns_bgp
             i = bc_hv%i_bgp(m)
             j = bc_hv%j_bgp(m)
             k = bc_hv%k_bgp(m)
             pt(k,j,i) = pt(k,j,i) + pt_surface_initial_change
          ENDDO

          DO  m = 1, surf_def%ns
             i = surf_def%i(m) + surf_def%ioff(m)
             j = surf_def%j(m) + surf_def%joff(m)
             k = surf_def%k(m) + surf_def%koff(m)
             surf_def%pt_surface(m) = pt(k,j,i)
          ENDDO
          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m) + surf_lsm%ioff(m)
             j = surf_lsm%j(m) + surf_lsm%joff(m)
             k = surf_lsm%k(m) + surf_lsm%koff(m)
             surf_lsm%pt_surface(m) = pt(k,j,i)
          ENDDO
          DO  m = 1, surf_usm%ns
             i = surf_usm%i(m) + surf_usm%ioff(m)
             j = surf_usm%j(m) + surf_usm%joff(m)
             k = surf_usm%k(m) + surf_usm%koff(m)
             surf_usm%pt_surface(m) = pt(k,j,i)
          ENDDO
       ENDIF
!
!--    If required, change the surface humidity/scalar at the start of the 3D run.
!--    Note, to store boundary values of passive scalar on the surface types is not require at
!--    at the moment.
       IF ( humidity  .AND.  q_surface_initial_change /= 0.0_wp )  THEN
          DO  m = 1, bc_hv%ns_bgp
             i = bc_hv%i_bgp(m)
             j = bc_hv%j_bgp(m)
             k = bc_hv%k_bgp(m)
             q(k,j,i) = q(k,j,i) + q_surface_initial_change
          ENDDO

          DO  m = 1, surf_def%ns
             i = surf_def%i(m) + surf_def%ioff(m)
             j = surf_def%j(m) + surf_def%joff(m)
             k = surf_def%k(m) + surf_def%koff(m)
             surf_def%q_surface(m) = q(k,j,i)
          ENDDO
          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m) + surf_lsm%ioff(m)
             j = surf_lsm%j(m) + surf_lsm%joff(m)
             k = surf_lsm%k(m) + surf_lsm%koff(m)
             surf_lsm%q_surface(m) = q(k,j,i)
          ENDDO
          DO  m = 1, surf_usm%ns
             i = surf_usm%i(m) + surf_usm%ioff(m)
             j = surf_usm%j(m) + surf_usm%joff(m)
             k = surf_usm%k(m) + surf_usm%koff(m)
             surf_usm%q_surface(m) = q(k,j,i)
          ENDDO
       ENDIF

       IF ( passive_scalar  .AND.  s_surface_initial_change /= 0.0_wp )  THEN
          DO  m = 1, bc_hv%ns_bgp
             i = bc_hv%i_bgp(m)
             j = bc_hv%j_bgp(m)
             k = bc_hv%k_bgp(m)
             s(k,j,i) = s(k,j,i) + s_surface_initial_change
          ENDDO
       ENDIF
    ENDIF

!
!-- If required, homogenize the surface temperature. This is useful to remove horizontal
!-- differences from the surface temperature in case the surface boundary condition for the
!-- temperature changes between two runs. More precisely, if the previous run used Neumann
!-- conditions but the current run uses Dirichlet conditions (e.g. to consider a rapid
!-- stabilization of the atmosphere in an idealized way). Therefore, compute the mean surface
!-- temperature and map this onto all boundary grid points. In case of a nested run, the surface
!-- temperature from the root domain is taken and send to all child domains, in order to
!-- guarantee the same surface temperature in all domains.
    IF ( homogenize_surface_temperature )  THEN
!
!--    Calculate mean surface potential temperature in the outermost domain. Loop over all
!--    boundary grid points.
       IF ( root_model )  THEN

          pt_surf_mean = 0.0_wp
          DO  m = 1, bc_hv%ns_bgp
             i = bc_hv%i_bgp(m)
             j = bc_hv%j_bgp(m)
             k = bc_hv%k_bgp(m)

             pt_surf_mean = pt_surf_mean + pt(k,j,i)
          ENDDO

#if defined( __parallel )
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, pt_surf_mean, 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#endif
          pt_surf_mean = pt_surf_mean / bc_hv%ns_bgp_tot
!
!--       Apply initial change of surface temperature in the restart run. Note,
!--       pt_surface_initial_change is not part of the restart IO, so that the change in the
!--       surface temperature needs to be explicitly set in the runtime_parameters namelist.
          IF ( TRIM( initializing_actions ) == 'read_restart_data'  .AND.                          &
               pt_surface_initial_change /= 0.0_wp )  THEN
             pt_surf_mean = pt_surf_mean + pt_surface_initial_change
          ENDIF
       ENDIF

       IF ( nested_run )  THEN
#if defined( __parallel )
          CALL MPI_BCAST( pt_surf_mean, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr )
#endif
       ENDIF
!
!--    Map surface temperature onto all boundary grid points. Note, surface temperature is also
!--    written onto the default surface array (land- and urban-surfaces are not allowed). This
!--    is only required to obtain the correct initial output of the surface temperature, which
!--    is stored in the %pt_surface attribute.
!--    Moreover, please note that the content of pt is copied onto pt_p at the end of init_3d_model,
!--    so that this is not necessary here.
       DO  m = 1, bc_hv%ns_bgp
          i = bc_hv%i_bgp(m)
          j = bc_hv%j_bgp(m)
          k = bc_hv%k_bgp(m)
          pt(k,j,i) = pt_surf_mean
       ENDDO

       DO  m = 1, surf_def%ns
          surf_def%pt_surface(m) = pt_surf_mean
       ENDDO
       DO  m = 1, surf_lsm%ns
          surf_lsm%pt_surface(m) = pt_surf_mean
       ENDDO
       DO  m = 1, surf_usm%ns
          surf_usm%pt_surface(m) = pt_surf_mean
       ENDDO
    ENDIF

 END SUBROUTINE dynamics_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific post-initialization checks
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_init_checks

    INTEGER(iwp) ::  i !< loop index in x-direction
    INTEGER(iwp) ::  j !< loop index in y-direction
    INTEGER(iwp) ::  k !< loop index in z-direction

    LOGICAL      ::  realistic_q = .TRUE. !< flag indicating realistic mixing ratios

    REAL(wp)     ::  e_s !< saturation water vapor pressure
    REAL(wp)     ::  q_s !< saturation mixing ratio
    REAL(wp)     ::  t_l !< actual temperature
    REAL(wp)     ::  rh_check = 9999999.9_wp !< relative humidity
    REAL(wp)     ::  rh_min = 9999999.9_wp !< max relative humidity
    REAL(wp)     ::  height = 9999999.9_wp !< height of supersaturated regions
    REAL(wp)     ::  min_height = 9999999.9_wp !< height of supersaturated regions

!
!-- Check for realistic initial mixing ratio. This must be in a realistic phyiscial range and must
!-- not exceed the saturation mixing ratio by more than 2 percent. Please note, the check is
!-- performed for each grid point (not just for a vertical profile), in order to cover also
!-- three-dimensional initialization. Note, this check gives an error only for the initial run not
!-- for a restart run. In case there are no cloud physics considered, the mixing ratio can exceed
!-- the saturation moisture. This case a warning is given.
    IF ( humidity  .AND.  .NOT. neutral  .AND.  check_realistic_q )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             Calculate actual temperature, water vapor saturation pressure, and based on this
!--             the saturation mixing ratio.
                t_l = exner(k) * pt(k,j,i)
                e_s = magnus( t_l )
                q_s = rd_d_rv * e_s / ( hyp(k) - e_s )

                IF ( q(k,j,i) > 1.02_wp * q_s )  THEN
                   realistic_q = .FALSE.
                   rh_check = q(k,j,i) / q_s * 100.0_wp
                   height = zu(k)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Since the check is performed locally, merge the logical flag from all mpi ranks,
!--    in order to do not print the error message multiple times.
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, realistic_q, 1, MPI_LOGICAL, MPI_LAND, comm2d, ierr)
       CALL MPI_ALLREDUCE( rh_check, rh_min, 1, MPI_REAL, MPI_MIN, comm2d, ierr )
       CALL MPI_ALLREDUCE( height, min_height, 1, MPI_REAL, MPI_MIN, comm2d, ierr )
#endif

       IF ( .NOT. realistic_q  .AND.  TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          WRITE( message_string, * ) 'initial mixing ratio exceeds the saturation mixing ratio' // &
                                     ', with rh = ', rh_min, '% at a height of ', min_height,      &
                                     'm for the first time (initial run)'
          CALL message( 'dynamic_init_checks', 'PAC0199', 2, 2, 0, 6, 0 )
       ELSEIF ( .NOT. realistic_q  .AND.  TRIM( initializing_actions ) == 'read_restart_data' ) THEN
          WRITE( message_string, * ) 'initial mixing ratio exceeds the saturation mixing ratio' // &
                                     ', with rh = ', rh_min, '% at a height of ', min_height,      &
                                     'm for the first time (restart run)'
          CALL message( 'dynamic_init_checks', 'PAC0200', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF

 END SUBROUTINE dynamics_init_checks


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the grids on which module-specific output quantities are defined. Allowed values for
!> grid_x are "x" and "xu", for grid_y "y" and "yv", and for grid_z "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )


    CHARACTER (LEN=*) ::  grid_x     !<
    CHARACTER (LEN=*) ::  grid_y     !<
    CHARACTER (LEN=*) ::  grid_z     !<
    CHARACTER (LEN=*) ::  variable   !<

    LOGICAL ::  found   !<


    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2' )

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT


 END SUBROUTINE dynamics_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print a header with module-specific information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_header( io )


    INTEGER(iwp) ::  io   !< output-file id

!
!-- Write dynamics module header
!-- NOTE: Deactivated because no relevant information to write so far
    IF ( .FALSE. )  WRITE ( io, 100 )

!
!-- Format-descriptors
100 FORMAT (//' Dynamics module information:'/                                              &
              ' -----------------------------------'//)

 END SUBROUTINE dynamics_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute module-specific actions for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_actions( location )


    CHARACTER (LEN=*) ::  location !<

    INTEGER(iwp) ::  i !< running grid index in x-direction
    INTEGER(iwp) ::  j !< running grid index in y-direction
    INTEGER(iwp) ::  k !< running grid index in z-direction

    LOGICAL, SAVE ::  rh_valid = .TRUE.

    REAL(wp) ::  e_s !< saturation water vapor pressure
    REAL(wp) ::  q_s !< saturation mixing ratio


!
!-- Here actions of the dynamics module follow.
!-- No calls for single grid points are allowed at locations before and after the timestep, since
!-- these calls are not within an i,j-loop
    SELECT CASE ( location )

       CASE ( 'before_timestep' )


       CASE ( 'before_prognostic_equations' )


       CASE ( 'after_integration' )
!
!--       Check values of the mixing ratio. They must lie in a realistic physical
!--       range and must not exceed the saturation mixing ratio by more than 2 percent, else a
!--       warning message is triggered. Note, this warning is only triggered one time to
!--       make users aware of potential issues. For example, in case of LSM/USM runs, the energy
!--       balance solvers might yield unrealistic fluxes causing strong numerical oscillations
!--       in the surface temperature resulting in numerical instabilities.
!--       Too high relative humidities may occur if no cloud water condensation is considered.
          IF ( humidity  .AND.  .NOT. neutral  .AND.  check_realistic_q  .AND.  rh_valid )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb+1, nzt
!
!--                   Calculate water vapor saturation pressure, and based on this, the saturation
!--                   mixing ratio.
                      e_s = magnus( exner(k) * pt(k,j,i) )
                      q_s = rd_d_rv * e_s / ( hyp(k) - e_s )

                      IF ( q(k,j,i) > 1.02_wp * q_s )  rh_valid = .FALSE.

                   ENDDO
                ENDDO
             ENDDO
!
!--          Since the check is performed locally, merge the logical flag from all mpi ranks.
#if defined( __parallel )
             CALL MPI_ALLREDUCE( MPI_IN_PLACE, rh_valid, 1, MPI_LOGICAL, MPI_LAND, comm2d, ierr)
#endif
             IF ( .NOT. rh_valid )  THEN
                WRITE( message_string, *) 'Mixing ratio exceeds the saturation mixing ratio for ', &
                                          'the first time after ', time_since_reference_point,     &
                                          ' s of simulated time.&No further warning will be ',     &
                                          'given. Please check your results carefully.'
                CALL message( 'dynamics_actions', 'PAC0201', 0, 1, 0, 6, 0 )
             ENDIF
          ENDIF

       CASE ( 'after_timestep' )


       CASE ( 'u-tendency' )


       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE dynamics_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute module-specific actions for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_actions_ij( i, j, location )


    CHARACTER (LEN=*) ::  location

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j

!
!-- Here the user-defined actions follow
    SELECT CASE ( location )

       CASE ( 'u-tendency' )

!
!--       Next line is to avoid compiler warning about unused variables. Please remove.
          IF ( i +  j < 0 )  CONTINUE

       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE dynamics_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific non-advective processes for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_non_advective_processes



 END SUBROUTINE dynamics_non_advective_processes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific non-advective processes for grid points i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_non_advective_processes_ij( i, j )


    INTEGER(iwp) ::  i                 !<
    INTEGER(iwp) ::  j                 !<

!
!--    Next line is just to avoid compiler warnings about unused variables. You may remove it.
       IF ( i + j < 0 )  CONTINUE


 END SUBROUTINE dynamics_non_advective_processes_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific horizontal boundary exchange
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_exchange_horiz( location )

       CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

       SELECT CASE ( location )

          CASE ( 'before_prognostic_equation' )

          CASE ( 'after_prognostic_equation' )

             CALL exchange_horiz( u_p, nbgp )
             CALL exchange_horiz( v_p, nbgp )
             CALL exchange_horiz( w_p, nbgp )
             CALL exchange_horiz( pt_p, nbgp )
             IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e_p, nbgp )
             IF ( rans_tke_e  )               CALL exchange_horiz( diss_p, nbgp )
             IF ( humidity )                  CALL exchange_horiz( q_p, nbgp )
             IF ( passive_scalar )            CALL exchange_horiz( s_p, nbgp )

          CASE ( 'after_anterpolation' )

             CALL exchange_horiz( u, nbgp )
             CALL exchange_horiz( v, nbgp )
             CALL exchange_horiz( w, nbgp )
             IF ( .NOT. neutral )             CALL exchange_horiz( pt, nbgp )
             IF ( humidity )                  CALL exchange_horiz( q, nbgp )
             IF ( passive_scalar )            CALL exchange_horiz( s, nbgp )
             IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e, nbgp )
             IF ( .NOT. constant_diffusion  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
                CALL exchange_horiz( diss, nbgp )
             ENDIF

       END SELECT

 END SUBROUTINE dynamics_exchange_horiz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific prognostic equations for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_prognostic_equations



 END SUBROUTINE dynamics_prognostic_equations


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific prognostic equations for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_prognostic_equations_ij( i, j, i_omp_start, tn )


    INTEGER(iwp), INTENT(IN) ::  i            !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  i_omp_start  !< first loop index of i-loop in prognostic_equations
    INTEGER(iwp), INTENT(IN) ::  j            !< grid index in y-direction
    INTEGER(iwp), INTENT(IN) ::  tn           !< task number of openmp task

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( i + j + i_omp_start + tn < 0 )  CONTINUE

 END SUBROUTINE dynamics_prognostic_equations_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute boundary conditions of dynamics model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_boundary_conditions

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index x direction
    INTEGER(iwp) ::  j  !< grid index y direction
    INTEGER(iwp) ::  k  !< grid index z direction
    INTEGER(iwp) ::  m  !< running index surface elements

!
!-- Bottom boundary
    IF ( ibc_uv_b == 1 )  THEN
       u_p(nzb,:,:) = u_p(nzb+1,:,:)
       v_p(nzb,:,:) = v_p(nzb+1,:,:)
    ENDIF
!
!-- Set zero vertical velocity at topography adjacent topography grid points.
!-- Is this really necessary?
    !$OMP PARALLEL DO PRIVATE( i, j, k )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
    !$ACC PRESENT(bc_hv, w_p)
    DO  m = 1, bc_hv%ns
       i = bc_hv%i(m)
       j = bc_hv%j(m)
       k = bc_hv%k(m)
       w_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = 0.0_wp
    ENDDO
!
!-- Top boundary. A nested domain ( ibc_uv_t = 3 ) does not require settings.
    IF ( ibc_uv_t == 0 )  THEN
        !$ACC KERNELS PRESENT(u_p, v_p, u_init, v_init)
        u_p(nzt+1,:,:) = u_init(nzt+1)
        v_p(nzt+1,:,:) = v_init(nzt+1)
        !$ACC END KERNELS
    ELSEIF ( ibc_uv_t == 1 )  THEN
        u_p(nzt+1,:,:) = u_p(nzt,:,:)
        v_p(nzt+1,:,:) = v_p(nzt,:,:)
    ENDIF
!
!-- Vertical velocity is zero at the model top, but maybe non-zero in child domains and
!-- in case of offline nesting.
    IF ( .NOT. child_domain  .AND.  .NOT. nesting_offline )  THEN
       !$ACC KERNELS PRESENT(w_p)
       w_p(nzt:nzt+1,:,:) = 0.0_wp  !< nzt is not a prognostic level (but cf. pres)
       !$ACC END KERNELS
    ENDIF

!
!-- Temperature at bottom and top boundary.
!-- In case of coupled runs (ibc_pt_b = 2) the temperature is given by the sea surface temperature
!-- of the coupled ocean model.
    IF ( .NOT. neutral )  THEN
!
!--    Dirichlet
       IF ( ibc_pt_b == 0 )  THEN
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m) + bc_hv%ioff(m)
             j = bc_hv%j(m) + bc_hv%joff(m)
             k = bc_hv%k(m) + bc_hv%koff(m)
             pt_p(k,j,i) = pt(k,j,i)
          ENDDO
!
!--       Increase temperature pt(0) according to pt_surface_heating_rate (convert from K/h to K/s).
!--       Do this each Runge-Kutta substep.
          IF ( pt_surface_heating_rate /= 0.0_wp ) THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k )
             DO  m = 1, bc_hv%ns_bgp
                i = bc_hv%i_bgp(m)
                j = bc_hv%j_bgp(m)
                k = bc_hv%k_bgp(m)
                pt_p(k,j,i) = pt_p(k,j,i) + dt_3d * pt_surface_heating_rate / 3600.0_wp            &
                                          * weight_substep(intermediate_timestep_count)
             ENDDO
          ENDIF
!
!--    Neumann, zero-gradient
       ELSEIF ( ibc_pt_b == 1 )  THEN
          !$OMP PARALLEL DO PRIVATE( i, j, k )
          !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
          !$ACC PRESENT(bc_hv, pt_p)
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m)
             j = bc_hv%j(m)
             k = bc_hv%k(m)
             pt_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = pt_p(k,j,i)
          ENDDO
       ENDIF

!
!--    Temperature at top boundary
       IF ( ibc_pt_t == 0 )  THEN
           pt_p(nzt+1,:,:) = pt(nzt+1,:,:)
!
!--        In case of nudging adjust top boundary to pt which is
!--        read in from NUDGING-DATA
           IF ( nudging )  THEN
              pt_p(nzt+1,:,:) = pt_init(nzt+1)
           ENDIF
       ELSEIF ( ibc_pt_t == 1 )  THEN
           pt_p(nzt+1,:,:) = pt_p(nzt,:,:)
       ELSEIF ( ibc_pt_t == 2 )  THEN
           !$ACC KERNELS PRESENT(pt_p, dzu)
           pt_p(nzt+1,:,:) = pt_p(nzt,:,:) + bc_pt_t_val * dzu(nzt+1)
           !$ACC END KERNELS
       ENDIF
    ENDIF
!
!-- Boundary conditions for total water content, bottom and top boundary (see also temperature)
    IF ( humidity )  THEN
!
!--    Surface conditions for constant_humidity_flux
!--    Run loop over all non-natural and natural walls. Note, in wall-datatype the k coordinate
!--    belongs to the atmospheric grid point, therefore, set q_p at k-1
       IF ( ibc_q_b == 0 ) THEN

          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m)
             j = bc_hv%j(m)
             k = bc_hv%k(m)
             q_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) =                                &
                                             q(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m))
          ENDDO

       ELSE

          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m)
             j = bc_hv%j(m)
             k = bc_hv%k(m)
             q_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = q_p(k,j,i)
          ENDDO
       ENDIF
!
!--    Top boundary
       IF ( ibc_q_t == 0 ) THEN
          q_p(nzt+1,:,:) = q(nzt+1,:,:)
       ELSEIF ( ibc_q_t == 1 ) THEN
          q_p(nzt+1,:,:) = q_p(nzt,:,:) + bc_q_t_val * dzu(nzt+1)
       ENDIF
    ENDIF
!
!-- Boundary conditions for scalar, bottom and top boundary (see also temperature)
    IF ( passive_scalar )  THEN
!
!--    Surface conditions for constant_humidity_flux
!--    Run loop over all non-natural and natural walls. Note, in wall-datatype the k coordinate
!--    belongs to the atmospheric grid point, therefore, set s_p at k-1
       IF ( ibc_s_b == 0 ) THEN

          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m)
             j = bc_hv%j(m)
             k = bc_hv%k(m)
             s_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) =                                &
                                             s(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m))
          ENDDO

       ELSE

          !$OMP PARALLEL DO PRIVATE( i, j, k )
          DO  m = 1, bc_hv%ns
             i = bc_hv%i(m)
             j = bc_hv%j(m)
             k = bc_hv%k(m)
             s_p(k+bc_hv%koff(m),j+bc_hv%joff(m),i+bc_hv%ioff(m)) = s_p(k,j,i)
          ENDDO
       ENDIF
!
!--    Top boundary condition
       IF ( ibc_s_t == 0 )  THEN
          s_p(nzt+1,:,:) = s(nzt+1,:,:)
       ELSEIF ( ibc_s_t == 1 )  THEN
          s_p(nzt+1,:,:) = s_p(nzt,:,:)
       ELSEIF ( ibc_s_t == 2 )  THEN
          s_p(nzt+1,:,:) = s_p(nzt,:,:) + bc_s_t_val * dzu(nzt+1)
       ENDIF

    ENDIF
!
!-- In case of inflow or nest boundary at the south boundary the boundary for v is at nys and in
!-- case of inflow or nest boundary at the left boundary the boundary for u is at nxl. Since in
!-- prognostic_equations (cache optimized version) these levels are handled as a prognostic level,
!-- boundary values have to be restored here.
    IF ( bc_dirichlet_s )  THEN
       v_p(:,nys,:) = v_p(:,nys-1,:)
    ELSEIF ( bc_dirichlet_l ) THEN
       u_p(:,:,nxl) = u_p(:,:,nxl-1)
    ENDIF

!
!-- Lateral boundary conditions for scalar quantities at the outflow.
!-- Lateral oundary conditions for TKE and dissipation are set in tcm_boundary_conds.
    IF ( bc_radiation_s )  THEN
       pt_p(:,nys-1,:)     = pt_p(:,nys,:)
       IF ( humidity )  THEN
          q_p(:,nys-1,:) = q_p(:,nys,:)
       ENDIF
       IF ( passive_scalar )  s_p(:,nys-1,:) = s_p(:,nys,:)
    ELSEIF ( bc_radiation_n )  THEN
       pt_p(:,nyn+1,:)     = pt_p(:,nyn,:)
       IF ( humidity )  THEN
          q_p(:,nyn+1,:) = q_p(:,nyn,:)
       ENDIF
       IF ( passive_scalar )  s_p(:,nyn+1,:) = s_p(:,nyn,:)
    ELSEIF ( bc_radiation_l )  THEN
       pt_p(:,:,nxl-1)     = pt_p(:,:,nxl)
       IF ( humidity )  THEN
          q_p(:,:,nxl-1) = q_p(:,:,nxl)
       ENDIF
       IF ( passive_scalar )  s_p(:,:,nxl-1) = s_p(:,:,nxl)
    ELSEIF ( bc_radiation_r )  THEN
       pt_p(:,:,nxr+1)     = pt_p(:,:,nxr)
       IF ( humidity )  THEN
          q_p(:,:,nxr+1) = q_p(:,:,nxr)
       ENDIF
       IF ( passive_scalar )  s_p(:,:,nxr+1) = s_p(:,:,nxr)
    ENDIF

!
!-- Radiation boundary conditions for the velocities at the respective outflow.
!-- The phase velocity is set to the maximum phase velocity that ensures numerical
!-- stability (CFL-condition), i.e. a Courant number of one is assumed.
    IF ( bc_radiation_s )  THEN
       u_p(:,-1,:) = u(:,0,:)
       v_p(:,0,:)  = v(:,1,:)
       w_p(:,-1,:) = w(:,0,:)
    ENDIF

    IF ( bc_radiation_n )  THEN
       u_p(:,ny+1,:) = u(:,ny,:)
       v_p(:,ny+1,:) = v(:,ny,:)
       w_p(:,ny+1,:) = w(:,ny,:)
    ENDIF

    IF ( bc_radiation_l )  THEN
       u_p(:,:,0)  = u(:,:,1)
       v_p(:,:,-1) = v(:,:,0)
       w_p(:,:,-1) = w(:,:,0)
    ENDIF

    IF ( bc_radiation_r )  THEN
       u_p(:,:,nx+1) = u(:,:,nx)
       v_p(:,:,nx+1) = v(:,:,nx)
       w_p(:,:,nx+1) = w(:,:,nx)
    ENDIF

 END SUBROUTINE dynamics_boundary_conditions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Swap timelevels of module-specific array pointers
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_swap_timelevel ( mod_count )


    INTEGER, INTENT(IN) :: mod_count


    SELECT CASE ( mod_count )

       CASE ( 0 )

          u  => u_1;   u_p  => u_2
          v  => v_1;   v_p  => v_2
          w  => w_1;   w_p  => w_2
          IF ( .NOT. neutral )  THEN
             pt => pt_1;  pt_p => pt_2
          ENDIF
          IF ( humidity )  THEN
             q => q_1;    q_p => q_2
          ENDIF
          IF ( passive_scalar )  THEN
             s => s_1;    s_p => s_2
          ENDIF

       CASE ( 1 )

          u  => u_2;   u_p  => u_1
          v  => v_2;   v_p  => v_1
          w  => w_2;   w_p  => w_1
          IF ( .NOT. neutral )  THEN
             pt => pt_2;  pt_p => pt_1
          ENDIF
          IF ( humidity )  THEN
             q => q_2;    q_p => q_1
          ENDIF
          IF ( passive_scalar )  THEN
             s => s_2;    s_p => s_1
          ENDIF

    END SELECT

 END SUBROUTINE dynamics_swap_timelevel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average module-specific output quantities as well as allocate the array necessary
!> for storing the average.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_3d_data_averaging( mode, variable )


    CHARACTER (LEN=*) ::  mode    !<
    CHARACTER (LEN=*) :: variable !<


    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

!          CASE ( 'u2' )

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

!          CASE ( 'u2' )

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

!          CASE ( 'u2' )

       END SELECT

    ENDIF

 END SUBROUTINE dynamics_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the module-specific output quantity with indices (k,j,i) to a temporary array with
!> indices (i,j,k) and sets the grid on which it is defined.
!> Allowed values for grid are "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do,     &
                                     nzt_do )


    CHARACTER (LEN=*)             ::  grid       !<
    CHARACTER (LEN=*), INTENT(IN) ::  mode       !< either 'xy', 'xz' or 'yz'
    CHARACTER (LEN=*)             ::  variable   !<

    INTEGER(iwp) ::  av     !< flag to control data output of instantaneous or time-averaged data
!    INTEGER(iwp) ::  i      !< grid index along x-direction
!    INTEGER(iwp) ::  j      !< grid index along y-direction
!    INTEGER(iwp) ::  k      !< grid index along z-direction
!    INTEGER(iwp) ::  m      !< running index surface elements
    INTEGER(iwp) ::  nzb_do !< lower limit of the domain (usually nzb)
    INTEGER(iwp) ::  nzt_do !< upper limit of the domain (usually nzt+1)

    LOGICAL      ::  found !<
    LOGICAL      ::  two_d !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !<

!
!-- Next line is just to avoid compiler warnings about unused variables. You may remove it.
    IF ( two_d .AND. av + LEN( mode ) + local_pf(nxl,nys,nzb_do) < 0.0 )  CONTINUE

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2_xy', 'u2_xz', 'u2_yz' )

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT


 END SUBROUTINE dynamics_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the module-specific output quantity with indices (k,j,i) to a temporary array with
!> indices (i,j,k).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av    !<
!    INTEGER(iwp) ::  i     !<
!    INTEGER(iwp) ::  j     !<
!    INTEGER(iwp) ::  k     !<
    INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found !<

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av + local_pf(nxl,nys,nzb_do) < 0.0 )  CONTINUE


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!       CASE ( 'u2' )

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE dynamics_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of module-specific statistics, i.e. horizontally averaged profiles and time series.
!> This is called for every statistic region sr, but at least for the region "total domain" (sr=0).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_statistics( mode, sr, tn )


    CHARACTER (LEN=*) ::  mode   !<
!    INTEGER(iwp) ::  i    !<
!    INTEGER(iwp) ::  j    !<
!    INTEGER(iwp) ::  k    !<
    INTEGER(iwp) ::  sr   !<
    INTEGER(iwp) ::  tn   !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( sr == 0  .OR.  tn == 0 )  CONTINUE

    IF ( mode == 'profiles' )  THEN

    ELSEIF ( mode == 'time_series' )  THEN

    ENDIF

 END SUBROUTINE dynamics_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_rrd_global_ftn( found )

    LOGICAL, INTENT(OUT)  ::  found


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'global_paramter' )
!          READ ( 13 )  global_parameter

       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE dynamics_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_rrd_global_mpi

!    CALL rrd_mpi_io( 'global_parameter', global_parameter )
    CONTINUE

 END SUBROUTINE dynamics_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!> Subdomain index limits on file are given by nxl_on_file, etc.
!> Indices nxlc, etc. indicate the range of gridpoints to be mapped from the subdomain on file (f)
!> to the subdomain of the current PE (c). They have been calculated in routine rrd_local.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf,     &
                                    nync, nyn_on_file, nysf, nysc, nys_on_file, tmp_2d, tmp_3d,    &
                                    found )


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
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( k + nxlc + nxlf + nxrc + nxrf + nync + nynf + nysc + nysf +                               &
         tmp_2d(nys_on_file,nxl_on_file) +                                                         &
         tmp_3d(nzb,nys_on_file,nxl_on_file) < 0.0 )  CONTINUE
!
!-- Here the reading of user-defined restart data follows:
!-- Sample for user-defined output

    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

!       CASE ( 'u2_av' )

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE dynamics_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_rrd_local_mpi

    IMPLICIT NONE

!    LOGICAL ::  array_found  !<

!
!--  Restart input of time-averaged quantities is skipped in case of cyclic-fill initialization.
!--  This case, input of time-averaged data is useless and can lead to faulty averaging.
!    IF ( .NOT. cyclic_fill_initialization )  THEN

!       CALL rd_mpi_io_check_array( 'u2_av' , found = array_found )
!       IF ( array_found )  THEN
!          IF ( .NOT. ALLOCATED( u2_av ) )  ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!          CALL rrd_mpi_io( 'u2_av', u2_av )
!       ENDIF

!    ENDIF

    CONTINUE

 END SUBROUTINE dynamics_rrd_local_mpi



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes global module-specific restart data into binary file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_wrd_global


 END SUBROUTINE dynamics_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes processor specific and module-specific restart data into binary file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_wrd_local


 END SUBROUTINE dynamics_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execute module-specific actions at the very end of the program.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dynamics_last_actions


 END SUBROUTINE dynamics_last_actions

 END MODULE dynamics_mod
