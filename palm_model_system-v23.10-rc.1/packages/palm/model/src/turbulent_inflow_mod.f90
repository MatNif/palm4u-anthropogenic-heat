!> @file turbulent_inflow_mod.f90
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
! Copyright 1997-2021  Leibniz Universitaet Hannover
! Copyright 2022-2022  pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Turbulent inflow module. This module includes two approaches to guarantee a turbulent inflow,
!> i.e. the turbulence recycling method according to Kataoka and Mizuno (2002), and the so-called
!> read-from-file method where the turbulent inflow data has been sampled in a precursor simulation.
!--------------------------------------------------------------------------------------------------!
 MODULE turbulent_inflow_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  drho_air,                                                                           &
               dzw,                                                                                &
               e,                                                                                  &
               e_p,                                                                                &
               mean_inflow_profiles,                                                               &
               pt,                                                                                 &
               pt_p,                                                                               &
               q,                                                                                  &
               q_p,                                                                                &
               rho_air,                                                                            &
               s,                                                                                  &
               s_p,                                                                                &
               u,                                                                                  &
               u_init,                                                                             &
               u_p,                                                                                &
               v,                                                                                  &
               v_init,                                                                             &
               v_p,                                                                                &
               w,                                                                                  &
               w_p,                                                                                &
               zu,                                                                                 &
               zw

    USE control_parameters,                                                                        &
        ONLY:  bc_dirichlet_l,                                                                     &
               bc_radiation_r,                                                                     &
               bc_lr,                                                                              &
               child_domain,                                                                       &
               constant_diffusion,                                                                 &
               coupling_char,                                                                      &
               cyclic_fill_initialization,                                                         &
               dt_3d,                                                                              &
               end_time,                                                                           &
               humidity,                                                                           &
               initializing_actions,                                                               &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               length,                                                                             &
               message_string,                                                                     &
               neutral,                                                                            &
               num_mean_inflow_profiles,                                                           &
               nesting_offline,                                                                    &
               passive_scalar,                                                                     &
               restart_data_format_output,                                                         &
               restart_string,                                                                     &
               terrain_following_mapping,                                                          &
               spinup_time,                                                                        &
               syn_turb_gen,                                                                       &
               time_since_reference_point,                                                         &
               turbulent_inflow,                                                                   &
               turbulent_outflow,                                                                  &
               use_prescribed_profile_data,                                                        &
               y_shift

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE general_utilities,                                                                         &
        ONLY:  interpolate_linear

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxlu,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nys,                                                                                &
               nysg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nzb,                                                                                &
               nz,                                                                                 &
               nzt,                                                                                &
               topo_flags,                                                                         &
               topo_top_ind

    USE kinds

    USE pegrid

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  char_fill,                                                                          &
               check_existence,                                                                    &
               close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               input_pids_dynamic,                                                                 &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               input_file_dynamic,                                                                 &
               num_var_pids,                                                                       &
               open_read_file,                                                                     &
               pids_id

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rrd_mpi_io,                                                                         &
               rd_mpi_io_check_array,                                                              &
               rrd_mpi_io_global_array,                                                            &
               wrd_mpi_io,                                                                         &
               wrd_mpi_io_global_array

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               hom_sum,                                                                            &
               pr_palm

    IMPLICIT NONE

    CHARACTER(LEN=80) ::  turbulent_inflow_method = 'recycle_turbulent_fluctuation'   !< namelist parameter

    INTEGER(iwp) ::  id_inflow = 0         !< myidx of procs at inflow (turbulent inflow method)
    INTEGER(iwp) ::  id_recycling = 0      !< myidx of procs containing the recycling plane (turbulence recycling method)
    INTEGER(iwp) ::  input_block_size = 30 !< namelist parameter: number of input timesteps per IO call
    INTEGER(iwp) ::  recycling_plane       !< position of recycling plane along x (in grid points) in case of turbulence recycling

    LOGICAL ::  turbulent_inflow_rec = .FALSE. !< internal switch for turbulence recycling
    LOGICAL ::  turbulent_inflow_rff = .FALSE. !< internal switch for the read-from-file inflow method

    REAL(wp) ::  fac_dt                                !< factor in time interpolation
    REAL(wp) ::  inflow_damping_height = 9999999.9_wp  !< namelist parameter
    REAL(wp) ::  inflow_damping_width = 9999999.9_wp   !< namelist parameter
    REAL(wp) ::  recycling_width = HUGE( 1.0_wp )      !< namelist parameter

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  inflow_damping_factor !< used for turbulent inflow (non-cyclic boundary conditions)

    TYPE inflow
!
!--    Define input variable strings.
       CHARACTER(LEN=14) ::  char_e    = 'inflow_plane_e'  !< variable name for SGS-TKE input
       CHARACTER(LEN=15) ::  char_pt   = 'inflow_plane_pt' !< variable name for potential temperature input
       CHARACTER(LEN=14) ::  char_q    = 'inflow_plane_q'  !< variable name for mixing ratio input
       CHARACTER(LEN=11) ::  char_time = 'time_inflow'     !< variable name for time input
       CHARACTER(LEN=1)  ::  char_y    = 'y'               !< variable name for y input
       CHARACTER(LEN=14) ::  char_u    = 'inflow_plane_u'  !< variable name for u input
       CHARACTER(LEN=14) ::  char_v    = 'inflow_plane_v'  !< variable name for v input
       CHARACTER(LEN=14) ::  char_w    = 'inflow_plane_w'  !< variable name for w input
       CHARACTER(LEN=2)  ::  char_zu   = 'z'               !< variable name for zu input  !!!!!!! should be 'z' !!!!
       CHARACTER(LEN=2)  ::  char_zw   = 'zw'              !< variable name for zw input
!
!--    Define different kinds of indices.
       INTEGER(iwp) ::  js   !< start index in y-dimension
       INTEGER(iwp) ::  je   !< end index in y-dimension
       INTEGER(iwp) ::  nt   !< dimension length of time_inflow
       INTEGER(iwp) ::  nts  !< number of timesteps to be read per IO call
       INTEGER(iwp) ::  ny   !< dimension length of y
       INTEGER(iwp) ::  nzu  !< dimension length of zu
       INTEGER(iwp) ::  nzw  !< dimension length of zw
       INTEGER(iwp) ::  t    !< time index at current time
       INTEGER(iwp) ::  tp   !< time index at current time + dt_3d
       INTEGER(iwp) ::  ts   !< start index of currently loaded time array

       REAL(wp) ::  fill_time  !< _FillValue for time variable
       REAL(wp) ::  fill_var   !< _FillValue for physical input variables
!
!--    Define dimensions of inflow data.
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  time  !< time levels in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  y     !< y levels in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu    !< vertical levels at scalar grid in dynamic input file
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw    !< vertical levels at w grid in dynamic input file
!
!--    Define arrays for linearly interpolated boundary values.
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  e_in   !< interpolated SGS-TKE at left boundary
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pt_in  !< interpolated potentital temperature at left boundary
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  q_in   !< interpolated mixing ratio at left boundary
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  u_in   !< interpolated u-component at left boundary
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  v_in   !< interpolated v-component at left boundary
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  w_in   !< interpolated w-component at left boundary
!
!--    Define arrays for inflow quantities.
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  e_l   !< input SGS-TKE at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  pt_l  !< input potentital temperature at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  q_l   !< input mixing ratio at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u_l   !< input u-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  v_l   !< input v-component at left boundary
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  w_l   !< input w-component at left boundary

    END TYPE inflow

    TYPE(inflow) ::  inflow_data  !< derived type which contains data for read-from-file inflow

    PRIVATE
!
!-- Public subroutines.
    PUBLIC turbulent_inflow_check_parameters,                                                      &
           turbulent_inflow_header,                                                                &
           turbulent_inflow_impose,                                                                &
           turbulent_inflow_init,                                                                  &
           turbulent_inflow_init_arrays,                                                           &
           turbulent_inflow_parin,                                                                 &
           turbulent_inflow_rrd_global,                                                            &
           turbulent_inflow_wrd_global
!
!-- Public variables.
    PUBLIC turbulent_inflow_method

    INTERFACE turbulent_inflow_check_parameters
       MODULE PROCEDURE turbulent_inflow_check_parameters
    END INTERFACE turbulent_inflow_check_parameters

    INTERFACE turbulent_inflow_header
       MODULE PROCEDURE turbulent_inflow_header
    END INTERFACE turbulent_inflow_header

    INTERFACE turbulent_inflow_impose
       MODULE PROCEDURE turbulent_inflow_impose
    END INTERFACE turbulent_inflow_impose

    INTERFACE turbulent_inflow_init
       MODULE PROCEDURE turbulent_inflow_init
    END INTERFACE turbulent_inflow_init

    INTERFACE turbulent_inflow_init_arrays
       MODULE PROCEDURE turbulent_inflow_init_arrays
    END INTERFACE turbulent_inflow_init_arrays

    INTERFACE turbulent_inflow_parin
       MODULE PROCEDURE turbulent_inflow_parin
    END INTERFACE turbulent_inflow_parin

    INTERFACE turbulent_inflow_rrd_global
       MODULE PROCEDURE turbulent_inflow_rrd_global_ftn
       MODULE PROCEDURE turbulent_inflow_rrd_global_mpi
    END INTERFACE turbulent_inflow_rrd_global

    INTERFACE turbulent_inflow_wrd_global
       MODULE PROCEDURE turbulent_inflow_wrd_global
    END INTERFACE turbulent_inflow_wrd_global

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Performs consistency checks.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_check_parameters

!
!-- Check for correct input of the inflow method.
    IF ( TRIM( turbulent_inflow_method ) /= 'read_from_file'  .AND.                                &
         TRIM( turbulent_inflow_method ) /= 'recycle_turbulent_fluctuation'  .AND.                 &
         TRIM( turbulent_inflow_method ) /= 'recycle_absolute_value'         .AND.                 &
         TRIM( turbulent_inflow_method ) /= 'recycle_absolute_value_thermodynamic' )               &
    THEN
       WRITE( message_string, * )  'unknown turbulent inflow method: ',                            &
                                   TRIM( turbulent_inflow_method )
       CALL message( 'turbulent_inflow_check_parameters', 'TUI0001', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set control switches for the recycling and the read-from-file method.
    turbulent_inflow_rff = ( INDEX( turbulent_inflow_method, 'read_from_file' ) /= 0 )
    turbulent_inflow_rec = ( INDEX( turbulent_inflow_method, 'recycle' ) /= 0 )

!
!-- A turbulent inflow requires Dirichlet conditions at the respective inflow boundary (so far, a
!-- turbulent inflow is realized from the left side only).
    IF ( turbulent_inflow_rec  .OR.  turbulent_inflow_rff )  THEN
!
!--    Turbulent inflow is not allowed in conjunction with mesoscale offline nesting.
       IF ( nesting_offline )  THEN
          message_string = 'turbulent inflow is not allowed in combination with mesoscale nesting'
          CALL message( 'turbulent_inflow_check_parameters', 'TUI0002', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Turbulent inflow is not allowed in nested child domains.
       IF ( child_domain )  THEN
          message_string = 'turbulent_inflow_method = "' // TRIM( turbulent_inflow_method ) //     &
                           '"& is not allowed in nested child domains'
          CALL message( 'turbulent_inflow_check_parameters', 'TUI0003', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( bc_lr /= 'dirichlet/radiation' )  THEN
          message_string = 'turbulent_inflow_method = "' // TRIM( turbulent_inflow_method ) //     &
                           '"& requires a Dirichlet condition at the inflow boundary'
          CALL message( 'turbulent_inflow_check_parameters', 'TUI0004', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    read-from-file inflow is not allowed in conjunction with a synthetic turbulence generator.
       IF ( syn_turb_gen )  THEN
          message_string = 'turbulent_inflow_method = "' // TRIM( turbulent_inflow_method ) //     &
                           '"& is not allowed in combination with the synthetic turbulence ' //    &
                           'generator'
          CALL message( 'turbulent_inflow_check_parameters', 'TUI0005', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF
!
!-- Specific checks for the turbulence recycling method.
    IF ( turbulent_inflow_rec )  THEN
!
!--    Turbulence recycling requires that 3d arrays have been cyclically filled with data from
!--    prerun in the first main run.
       IF ( .NOT. cyclic_fill_initialization  .AND.  initializing_actions /= 'read_restart_data' ) &
       THEN
          message_string = 'turbulence recycling requires ' //                                     &
                           'initializing_actions = ''cyclic_fill'' or ' //                         &
                           'initializing_actions = ''read_restart_data'' '
          CALL message( 'turbulent_inflow_check_parameters', 'TUI0006', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check setting of the recycling plane.
       IF ( recycling_width <= dx  .OR.  recycling_width >= nx * dx )  THEN
          WRITE( message_string, * )  'illegal value for recycling_width: ', recycling_width
          CALL message( 'turbulent_inflow_check_parameters', 'TUI0007', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF
!
!-- Specific checks for the read-from-file method.
    IF ( turbulent_inflow_rff )  THEN
!
!--    read-from-file inflow method requires a dynamic input file.
       IF ( .NOT. input_pids_dynamic )  THEN
          message_string = 'turbulent_inflow_method = "' // TRIM( turbulent_inflow_method ) //     &
                           '" requires a dynamic input file'
          CALL message( 'turbulent_inflow_check_parameters', 'TUI0008', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    read-from-file inflow cannot be combined with the outflow turbulence method as the mass
!--    conservation cannot be guaranteed.
       IF ( turbulent_outflow )  THEN
          message_string = 'turbulent_inflow_method = "' // TRIM( turbulent_inflow_method ) //     &
                           '"& is not allowed in combination with turbulent_outflow = .T.'
          CALL message( 'turbulent_inflow_check_parameters', 'TUI0009', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

 END SUBROUTINE turbulent_inflow_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_header( io )

    INTEGER(iwp), INTENT(IN) ::  io !< Unit of the output file


    WRITE ( io, 1 )

    IF ( turbulent_inflow_rec )  THEN
       IF ( y_shift == 0 )  THEN
          WRITE ( io, 2 )  recycling_width, recycling_plane, inflow_damping_height,                &
                           inflow_damping_width
       ELSE
          WRITE ( io, 3 )  y_shift, recycling_width, recycling_plane, inflow_damping_height,       &
                           inflow_damping_width

       ENDIF
    ENDIF

    IF ( turbulent_inflow_rff )  WRITE ( io, 4 )  inflow_data%time(inflow_data%nt-1)

1 FORMAT (//' Turbulent inflow at left lateral boundary is enabled.')

2 FORMAT (  '       Turbulence recycling at inflow switched on'/                                   &
            '       width of recycling domain: ',F7.1,' m   grid index: ',I4/                      &
            '       inflow damping height: ',F6.1,' m   width: ',F6.1,' m')

3 FORMAT (  '       turbulence recycling at inflow switched on'/                                   &
            '       y-shift of the recycled inflow turbulence is',I3,' PE'/                        &
            '       width of recycling domain: ',F7.1,' m   grid index: ',I4/                      &
            '       inflow damping height: ',F6.1,' m   width: ',F6.1,' m'/)

4 FORMAT (  '       Turbulent inflow data is read from dynamic input file.'/                       &
            '       Inflow data covers a time range of: ', F7.1, ' s.' )

 END SUBROUTINE turbulent_inflow_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Wrapper routine to call desired initialization routine.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_init

!
!-- Initialize read-from-file turbulent inflow method.
    IF ( turbulent_inflow_rff )  CALL turbulent_inflow_rff_init
!
!-- Initialize turbulence recycling method.
    IF ( turbulent_inflow_rec )  CALL turbulent_inflow_recycling_init

 END SUBROUTINE turbulent_inflow_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_init_arrays

    CONTINUE

 END SUBROUTINE turbulent_inflow_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Wrapper routine to call desired routine to impose either a recycled signal or boundary data read
!> from the dynamic input file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_impose

!
!-- Turbulent inflow with external data read from file. Map boundary data onto domain boundaries.
    IF ( turbulent_inflow_rff )  THEN
       CALL turbulent_inflow_rff_bc
       CALL turbulent_inflow_mass_flux_conservation
    ENDIF
!
!-- Turbulence recycling method.
    IF ( turbulent_inflow_rec )  CALL turbulent_inflow_recycling

 END SUBROUTINE turbulent_inflow_impose


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Mass-flux conservation used for the turbulent inflow method with external data read from file.
!> A mass-flux conservation is required to maintain stability, i.e. to guarantee a divergence-free
!> velocity field at each timestep. Since the read-from-file inflow method allows instantionary
!> flows, the mass flux can also change between two timesteps.
!> The mass flux is corrected via the right (eastern) outflow boundary.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_mass_flux_conservation

    INTEGER(iwp) ::  i !< grid index in x-direction
    INTEGER(iwp) ::  j !< grid index in y-direction
    INTEGER(iwp) ::  k !< grid index in z-direction

    REAL(wp), DIMENSION(nzb+1:nzt) ::  u_corr     !< correction velocity required to balance mass flux
    REAL(wp), DIMENSION(nzb+1:nzt) ::  vol_flux_l !< volume flux at west domain boundary
    REAL(wp), DIMENSION(nzb+1:nzt) ::  vol_flux_r !< volume flux at west domain boundary

    vol_flux_l = 0.0_wp
    vol_flux_r = 0.0_wp
!
!-- Compute mass flux at eastern outflow boundary.
    IF ( bc_radiation_r )  THEN
       i = nxr + 1
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             vol_flux_r(k) = vol_flux_r(k) - u(k,j,i) * dy * dzw(k) * rho_air(k)                   &
                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
       ENDDO
    ENDIF
!
!-- Compute mass flux at western inflow boundary.
    IF ( bc_dirichlet_l )  THEN
       i = nxlu-1
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             vol_flux_l(k) = vol_flux_l(k) + u(k,j,i) * dy * dzw(k) * rho_air(k)                   &
                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
       ENDDO
    ENDIF
!
!-- Sum-up over all subdomains.
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, vol_flux_l, nzt, MPI_REAL, MPI_SUM, comm2d, ierr )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, vol_flux_r, nzt, MPI_REAL, MPI_SUM, comm2d, ierr )
#endif
!
!-- Determine correction velocity required to balance the total mass flux.
    DO  k = nzb+1, nzt
       u_corr(k) = ( vol_flux_l(k) + vol_flux_r(k) ) * drho_air(k) / ( ( ny + 1 ) * dy * dzw(k) )
    ENDDO
!
!-- Finally, add the correction velocitiy at the eastern outflow boundary to balance the total
!-- mass flux.
    IF ( bc_radiation_r )  THEN
       i = nxr + 1
       DO  j = nysg, nyng
          DO  k = nzb+1, nzt
             u(k,j,i) = u(k,j,i) + u_corr(k)                                                       &
                                 * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
          u(nzt+1,j,i) = u(nzt,j,i)
       ENDDO
    ENDIF

 END SUBROUTINE turbulent_inflow_mass_flux_conservation


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads the parameter list turbulent_inflow_parameters.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /turbulent_inflow_parameters/  inflow_damping_height,                                 &
                                            inflow_damping_width,                                  &
                                            input_block_size,                                      &
                                            recycling_width,                                       &
                                            switch_off_module,                                     &
                                            turbulent_inflow_method

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, turbulent_inflow_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    turbulent_inflow_parameters namelist was found and read correctly. Enable turbulent inflow
!--    module.
       IF ( .NOT. switch_off_module )  turbulent_inflow = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    turbulent_inflow_parameters namelist was found but contained errors. Print an error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'turbulent_inflow_parameters', line )

    ENDIF

 END SUBROUTINE turbulent_inflow_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Imposing turbulence at the respective inflow using the turbulence recycling method of
!> Kataoka and Mizuno (2002).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_recycling

    INTEGER(iwp) ::  i        !< loop index
    INTEGER(iwp) ::  j        !< loop index
    INTEGER(iwp) ::  k        !< loop index
    INTEGER(iwp) ::  l        !< loop index
    INTEGER(iwp) ::  ngp_ifd  !< number of grid points stored in avpr
    INTEGER(iwp) ::  ngp_pr   !< number of grid points stored in inflow_dist

    REAL(wp), DIMENSION(nzb:nzt+1,num_mean_inflow_profiles,nbgp) ::  avpr   !< stores averaged profiles at recycling plane
    REAL(wp), DIMENSION(nzb:nzt+1,num_mean_inflow_profiles,nbgp) ::  avpr_l !< auxiliary variable to calculate avpr

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,num_mean_inflow_profiles,nbgp) ::  inflow_dist !< turbulence signal,
                                                                                           !< added at inflow boundary

    CALL cpu_log( log_point(40), 'turbulent_inflow_recycl', 'start' )

!
!-- Carry out spanwise averaging in the recycling plane
    avpr_l = 0.0_wp
    ngp_pr = ( nzt - nzb + 2 ) * num_mean_inflow_profiles * nbgp
    ngp_ifd = ngp_pr * ( nyn - nys + 1 + 2 * nbgp )

!
!-- First, local averaging within the recycling domain (not required if absolute values are
!-- recycled).
    i = recycling_plane

    IF ( myidx == id_recycling )  THEN

       DO  l = 1, nbgp
          DO  j = nys, nyn
             DO  k = nzb, nzt + 1

                avpr_l(k,1,l) = avpr_l(k,1,l) + u(k,j,i)
                avpr_l(k,2,l) = avpr_l(k,2,l) + v(k,j,i)
                avpr_l(k,3,l) = avpr_l(k,3,l) + w(k,j,i)
                IF ( .NOT. neutral )             avpr_l(k,4,l) = avpr_l(k,4,l) + pt(k,j,i)
                IF ( .NOT. constant_diffusion )  avpr_l(k,5,l) = avpr_l(k,5,l) + e(k,j,i)
                IF ( humidity )                  avpr_l(k,6,l) = avpr_l(k,6,l) + q(k,j,i)
                IF ( passive_scalar )            avpr_l(k,7,l) = avpr_l(k,7,l) + s(k,j,i)

             ENDDO
          ENDDO
          i = i + 1
       ENDDO

    ENDIF

#if defined( __parallel )
!
!-- Now, averaging over all PEs
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( avpr_l(nzb,1,1), avpr(nzb,1,1), ngp_pr, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    avpr = avpr_l
#endif

    avpr = avpr / ( ny + 1 )

!
!-- Calculate the disturbances at the recycling plane. In case of recycling of absolute quantities,
!-- the disturbance is defined as the absolute value (and not as the deviation from the mean
!-- profile).
    IF ( myidx == id_recycling )  THEN
!
!--    Recycling of turbulent fluctuations.
       IF ( TRIM( turbulent_inflow_method ) == 'recycle_turbulent_fluctuation' )  THEN

          i = recycling_plane
          DO  l = 1, nbgp
             DO  j = nysg, nyng
                DO  k = nzb, nzt + 1

                   inflow_dist(k,j,1,l) = u(k,j,i+1) - avpr(k,1,l)
                   inflow_dist(k,j,2,l) = v(k,j,i)   - avpr(k,2,l)
                   inflow_dist(k,j,3,l) = w(k,j,i)   - avpr(k,3,l)
                   inflow_dist(k,j,4,l) = pt(k,j,i)  - avpr(k,4,l)
                   inflow_dist(k,j,5,l) = e(k,j,i)   - avpr(k,5,l)
                   IF ( humidity )        inflow_dist(k,j,6,l) = q(k,j,i) - avpr(k,6,l)
                   IF ( passive_scalar )  inflow_dist(k,j,7,l) = s(k,j,i) - avpr(k,7,l)

                ENDDO
             ENDDO
             i = i + 1
          ENDDO
!
!--    Recycling of absolute values. Only the inflow profile for w demands zero mean.
       ELSEIF ( TRIM( turbulent_inflow_method ) == 'recycle_absolute_value' )  THEN

          i = recycling_plane
          DO  l = 1, nbgp
             DO  j = nysg, nyng
                DO  k = nzb, nzt + 1

                   inflow_dist(k,j,1,l) = u(k,j,i+1)
                   inflow_dist(k,j,2,l) = v(k,j,i)
                   inflow_dist(k,j,3,l) = w(k,j,i) - avpr(k,3,l)
                   inflow_dist(k,j,4,l) = pt(k,j,i)
                   inflow_dist(k,j,5,l) = e(k,j,i)
                   IF ( humidity )        inflow_dist(k,j,6,l) = q(k,j,i)
                   IF ( passive_scalar )  inflow_dist(k,j,7,l) = s(k,j,i)

                ENDDO
             ENDDO
             i = i + 1
          ENDDO
!
!--    Recycling of absolute values for thermodynamic quantities + passive scalar.
!--    For all other variables, only turbulent fluctuations are recycled.
       ELSEIF ( TRIM( turbulent_inflow_method ) == 'recycle_absolute_value_thermodynamic' )  THEN

          i = recycling_plane
          DO  l = 1, nbgp
             DO  j = nysg, nyng
                DO  k = nzb, nzt + 1

                   inflow_dist(k,j,1,l) = u(k,j,i+1) - avpr(k,1,l)
                   inflow_dist(k,j,2,l) = v(k,j,i)   - avpr(k,2,l)
                   inflow_dist(k,j,3,l) = w(k,j,i)   - avpr(k,3,l)
                   inflow_dist(k,j,4,l) = pt(k,j,i)
                   inflow_dist(k,j,5,l) = e(k,j,i)   - avpr(k,5,l)
                   IF ( humidity )        inflow_dist(k,j,6,l) = q(k,j,i)
                   IF ( passive_scalar )  inflow_dist(k,j,7,l) = s(k,j,i)

                ENDDO
             ENDDO
             i = i + 1
          ENDDO
       ENDIF

    ENDIF

!
!-- For parallel runs, send the disturbances to the respective inflow PE
#if defined( __parallel )
    IF ( myidx == id_recycling  .AND.  myidx /= id_inflow )  THEN

       CALL MPI_SEND( inflow_dist(nzb,nysg,1,1), ngp_ifd, MPI_REAL, id_inflow, 1, comm1dx, ierr )

    ELSEIF ( myidx /= id_recycling  .AND.  myidx == id_inflow )  THEN

       inflow_dist = 0.0_wp
       CALL MPI_RECV( inflow_dist(nzb,nysg,1,1), ngp_ifd, MPI_REAL, id_recycling, 1, comm1dx,      &
                      status, ierr )

    ENDIF
#endif
!
!-- If y-shift is set, exchange the inflow_dist array along the y-dimension.
    CALL  turbulent_inflow_y_shift( inflow_dist = inflow_dist )

!
!-- Add the disturbance at the inflow boundary.
    IF ( nxl == 0 )  THEN
!
!--    Impose turbulent fluctuations at the inflow plane.
       IF ( TRIM( turbulent_inflow_method ) == 'recycle_turbulent_fluctuation' )  THEN
          DO  j = nysg, nyng
             DO  k = nzb, nzt + 1
                u(k,j,-nbgp+1:0) = mean_inflow_profiles(k,1) +                                     &
                                   inflow_dist(k,j,1,1:nbgp) * inflow_damping_factor(k)
                v(k,j,-nbgp:-1)  = mean_inflow_profiles(k,2) +                                     &
                                   inflow_dist(k,j,2,1:nbgp) * inflow_damping_factor(k)
                w(k,j,-nbgp:-1)  = inflow_dist(k,j,3,1:nbgp) * inflow_damping_factor(k)
                pt(k,j,-nbgp:-1) = mean_inflow_profiles(k,4) +                                     &
                                   inflow_dist(k,j,4,1:nbgp) * inflow_damping_factor(k)
                e(k,j,-nbgp:-1)  = mean_inflow_profiles(k,5) +                                     &
                                   inflow_dist(k,j,5,1:nbgp) * inflow_damping_factor(k)
                IF ( humidity )  THEN
                   q(k,j,-nbgp:-1) = mean_inflow_profiles(k,6) +                                   &
                                     inflow_dist(k,j,6,1:nbgp) * inflow_damping_factor(k)
                ENDIF
                IF ( passive_scalar )  THEN
                   s(k,j,-nbgp:-1)  = mean_inflow_profiles(k,7) +                                  &
                                      inflow_dist(k,j,7,1:nbgp) * inflow_damping_factor(k)
                ENDIF
             ENDDO
          ENDDO
!
!--    Impose absolute values at at the inflow plane, except for w, where turbulent
!--    fluctuations are imposed.
       ELSEIF ( TRIM( turbulent_inflow_method ) == 'recycle_absolute_value' )  THEN
          DO  j = nysg, nyng
             DO  k = nzb, nzt + 1
                u(k,j,-nbgp+1:0) = inflow_dist(k,j,1,1:nbgp)
                v(k,j,-nbgp:-1)  = inflow_dist(k,j,2,1:nbgp)
                w(k,j,-nbgp:-1)  = inflow_dist(k,j,3,1:nbgp)
                pt(k,j,-nbgp:-1) = inflow_dist(k,j,4,1:nbgp)
                e(k,j,-nbgp:-1)  = inflow_dist(k,j,5,1:nbgp)
                IF ( humidity )        q(k,j,-nbgp:-1) = inflow_dist(k,j,6,1:nbgp)
                IF ( passive_scalar )  s(k,j,-nbgp:-1) = inflow_dist(k,j,7,1:nbgp)
             ENDDO
          ENDDO
!
!--    Impose turbulent perturbations at recycling plane for all dynamic quantities, except for
!--    the thermodynamic quantities (pt, q) and the passive scalar where absolute values are
!--    imposed.
       ELSEIF ( TRIM( turbulent_inflow_method ) == 'recycle_absolute_value_thermodynamic' )  THEN
          DO  j = nysg, nyng
             DO  k = nzb, nzt + 1
                u(k,j,-nbgp+1:0) = mean_inflow_profiles(k,1) +                                     &
                                   inflow_dist(k,j,1,1:nbgp) * inflow_damping_factor(k)
                v(k,j,-nbgp:-1)  = mean_inflow_profiles(k,2) +                                     &
                                   inflow_dist(k,j,2,1:nbgp) * inflow_damping_factor(k)
                w(k,j,-nbgp:-1)  = inflow_dist(k,j,3,1:nbgp) * inflow_damping_factor(k)
                pt(k,j,-nbgp:-1) = inflow_dist(k,j,4,1:nbgp)
                e(k,j,-nbgp:-1)  = mean_inflow_profiles(k,5) +                                     &
                                   inflow_dist(k,j,5,1:nbgp) * inflow_damping_factor(k)
                IF ( humidity )        q(k,j,-nbgp:-1) = inflow_dist(k,j,6,1:nbgp)
                IF ( passive_scalar )  s(k,j,-nbgp:-1) = inflow_dist(k,j,7,1:nbgp)
             ENDDO
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point(40), 'turbulent_inflow_recycl', 'stop' )

 END SUBROUTINE turbulent_inflow_recycling


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize the turbulence recycling method.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_recycling_init

    INTEGER(iwp) ::  i               !< running index in x-direction
    INTEGER(iwp) ::  j               !< running index in y-direction
    INTEGER(iwp) ::  k               !< running index in z-direction
    INTEGER(iwp) ::  nz_s_shift      !< topography-top index on scalar-grid, used to vertically shift initial profiles
    INTEGER(iwp) ::  nz_u_shift      !< topography-top index on u-grid, used to vertically shift initial profiles
    INTEGER(iwp) ::  nz_v_shift      !< topography-top index on v-grid, used to vertically shift initial profiles
    INTEGER(iwp) ::  nz_w_shift      !< topography-top index on w-grid, used to vertically shift initial profiles


!
!-- Calculate index of the recycling layer.
    recycling_plane = recycling_width / dx
!
!-- Broadcast the id of the inflow PE.
    IF ( bc_dirichlet_l )  THEN
       id_inflow = myidx
    ELSE
       id_inflow = 0
    ENDIF
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, id_inflow, 1, MPI_INTEGER, MPI_SUM, comm1dx, ierr )
#endif
!
!-- Broadcast the id of the recycling plane.
!-- WARNING: needs to be adjusted in case of inflows other than from left side!
    IF ( NINT( recycling_width / dx, KIND=idp ) >= nxl  .AND.                                      &
         NINT( recycling_width / dx, KIND=idp ) <= nxr )  THEN
       id_recycling = myidx
    ELSE
       id_recycling = 0
    ENDIF
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, id_recycling, 1, MPI_INTEGER, MPI_SUM, comm1dx, ierr )
#endif

!
!-- Initialization of the turbulence recycling method.
    IF ( cyclic_fill_initialization )  THEN
!
!--    First store the profiles to be used at the inflow.
!--    These profiles are the (temporally) and horizontally averaged vertical profiles from the
!--    prerun. Alternatively, prescribed profiles for u,v-components can be used.
       ALLOCATE( mean_inflow_profiles(nzb:nzt+1,1:num_mean_inflow_profiles) )

       IF ( use_prescribed_profile_data )  THEN
          mean_inflow_profiles(:,1) = u_init            ! u
          mean_inflow_profiles(:,2) = v_init            ! v
       ELSE
          mean_inflow_profiles(:,1) = hom_sum(:,1,0)    ! u
          mean_inflow_profiles(:,2) = hom_sum(:,2,0)    ! v
       ENDIF
       mean_inflow_profiles(:,4) = hom_sum(:,4,0)       ! pt
       IF ( .NOT. constant_diffusion )  mean_inflow_profiles(:,5) = hom_sum(:,8,0)   ! e
       IF ( humidity )  mean_inflow_profiles(:,6) = hom_sum(:,41,0)                  ! q
       IF ( passive_scalar )  mean_inflow_profiles(:,7) = hom_sum(:,115,0)           ! s

!
!--    In case of complex terrain, determine vertical displacement at inflow boundary and adjust
!--    mean inflow profiles.
       IF ( terrain_following_mapping )  THEN

          IF ( myidx == id_inflow )  THEN
             nz_u_shift = MINVAL( topo_top_ind(nys:nyn,nxl,1)   )
             nz_v_shift = MINVAL( topo_top_ind(nys:nyn,nxl-1,2) )
             nz_w_shift = MINVAL( topo_top_ind(nys:nyn,nxl-1,3) )
             nz_s_shift = MINVAL( topo_top_ind(nys:nyn,nxl-1,0) )
          ELSE
             nz_u_shift = 0
             nz_v_shift = 0
             nz_w_shift = 0
             nz_s_shift = 0
          ENDIF

#if defined( __parallel )
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, nz_u_shift, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, nz_v_shift, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, nz_w_shift, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, nz_s_shift, 1, MPI_INTEGER, MPI_MIN, comm2d, ierr )
#endif

          mean_inflow_profiles(:,1) = 0.0_wp
          mean_inflow_profiles(nz_u_shift:nzt+1,1) = hom_sum(0:nzt+1-nz_u_shift,1,0)  ! u

          mean_inflow_profiles(:,2) = 0.0_wp
          mean_inflow_profiles(nz_v_shift:nzt+1,2) = hom_sum(0:nzt+1-nz_v_shift,2,0)  ! v

          mean_inflow_profiles(nz_s_shift:nzt+1,4) = hom_sum(0:nzt+1-nz_s_shift,4,0)  ! pt

          IF ( .NOT. constant_diffusion )  THEN
             mean_inflow_profiles(nz_s_shift:nzt+1,5) = hom_sum(0:nzt+1-nz_s_shift,8,0)  ! e
          ENDIF

       ENDIF

!
!--    If necessary, adjust the horizontal flow field to the prescribed profiles.
       IF ( use_prescribed_profile_data )  THEN
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = nzb, nzt+1
                   u(k,j,i) = u(k,j,i) - hom_sum(k,1,0) + u_init(k)
                   v(k,j,i) = v(k,j,i) - hom_sum(k,2,0) + v_init(k)
                ENDDO
             ENDDO
          ENDDO
       ENDIF

!
!--    Use these mean profiles at the inflow (provided that Dirichlet conditions are used).
       IF ( bc_dirichlet_l )  THEN
          DO  j = nysg, nyng
             DO  k = nzb, nzt+1
                u(k,j,nxlg:-1)  = mean_inflow_profiles(k,1)
                v(k,j,nxlg:-1)  = mean_inflow_profiles(k,2)
                w(k,j,nxlg:-1)  = 0.0_wp
                pt(k,j,nxlg:-1) = mean_inflow_profiles(k,4)
                IF ( .NOT. constant_diffusion )  e(k,j,nxlg:-1)  = mean_inflow_profiles(k,5)
                IF ( humidity )  q(k,j,nxlg:-1)  = mean_inflow_profiles(k,6)
                IF ( passive_scalar )  s(k,j,nxlg:-1)  = mean_inflow_profiles(k,7)
             ENDDO
          ENDDO
       ENDIF

       u  = MERGE( u, 0.0_wp, BTEST( topo_flags, 1 ) )
       v  = MERGE( v, 0.0_wp, BTEST( topo_flags, 2 ) )
       w  = MERGE( w, 0.0_wp, BTEST( topo_flags, 3 ) )

!
!--    Calculate the damping factors to be used at the inflow. For a turbulent inflow the
!--    turbulent fluctuations have to be limited vertically because otherwise the turbulent
!--    inflow layer will grow in time.
       IF ( inflow_damping_height == 9999999.9_wp )  THEN
!
!--       Default: use the inversion height calculated by the prerun; if this is zero,
!--       inflow_damping_height must be explicitly specified.
          IF ( hom_sum(nzb+6,pr_palm,0) /= 0.0_wp )  THEN
             inflow_damping_height = hom_sum(nzb+6,pr_palm,0)
          ELSE
             WRITE( message_string, * ) 'inflow_damping_height must be ',                       &
                                        'explicitly specified because&the inversion height ',   &
                                        'calculated by the prerun is zero.'
             CALL message( 'init_3d_model', 'TUI0010', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       IF ( inflow_damping_width == 9999999.9_wp )  THEN
!
!--       Default for the transition range: one tenth of the undamped layer.
          inflow_damping_width = 0.1_wp * inflow_damping_height
       ENDIF

       ALLOCATE( inflow_damping_factor(nzb:nzt+1) )

       DO  k = nzb, nzt+1

          IF ( zu(k) <= inflow_damping_height )  THEN
             inflow_damping_factor(k) = 1.0_wp
          ELSEIF ( zu(k) <= ( inflow_damping_height + inflow_damping_width ) )  THEN
             inflow_damping_factor(k) = 1.0_wp -                                                &
                                        ( zu(k) - inflow_damping_height ) / inflow_damping_width
          ELSE
             inflow_damping_factor(k) = 0.0_wp
          ENDIF

       ENDDO

    ENDIF
!
!-- Initialize new time levels again (only done in order to set boundary values including ghost
!-- points).
    pt_p = pt; u_p = u; v_p = v; w_p = w
    IF ( .NOT. constant_diffusion )  e_p = e
    IF ( humidity )                  q_p = q
    IF ( passive_scalar )            s_p = s

 END SUBROUTINE turbulent_inflow_recycling_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Impose boundary values in case of the read-from-file inflow method.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_rff_bc

    INTEGER(iwp) ::  i       !< grid index in x-direction
    INTEGER(iwp) ::  iu      !< grid index in x-direction on u-grid
    INTEGER(iwp) ::  j       !< grid index in y-direction
    INTEGER(iwp) ::  k       !< grid index in z-direction
    INTEGER(iwp) ::  te      !< lower index of time dimension
    INTEGER(iwp) ::  ts      !< upper index of time dimension
    INTEGER(iwp) ::  tmp_nts !< dummy for the block size


    CALL cpu_log( log_point_s(66), 'turbulent_inflow_rff', 'start' )

!
!-- Check if input is required. This is only required at last RK3-substep.
    IF ( intermediate_timestep_count == intermediate_timestep_count_max )  THEN
       IF ( time_since_reference_point + dt_3d > inflow_data%time(inflow_data%ts+inflow_data%nts) )&
       THEN
!
!--       Determine start index. Decrement start index by one if simulation time is smaller than
!--       start time.
          inflow_data%ts = MINLOC( ABS( inflow_data%time - time_since_reference_point ), DIM = 1 ) &
                           - 1

          IF ( inflow_data%time(inflow_data%ts) > time_since_reference_point )  THEN
             inflow_data%ts = inflow_data%ts - 1
          ENDIF
!
!--       Determine the block size, i.e. the number of time layers to be read. This must not be
!--       larger than the total number of time layers as well as the number of remaining time layers.
!--       However, it must be sufficiently large to cover the time interval t+dt (this
!--       scenario might happen if the time resolution of the precursor is much higher than the
!--       main run, e.g. during the model spinup phase). The following do loop increases the block
!--       size as long as this condition is fulfilled. Note, this assumes that sufficient data is
!--       provided by the dynamic input file, which is checked during the initialization.
          tmp_nts = inflow_data%nts
          tmp_nts = MIN( inflow_data%nts, inflow_data%nt - 1 )
          tmp_nts = MIN( tmp_nts,         inflow_data%nt - 1 - inflow_data%ts )

          DO WHILE ( inflow_data%time(inflow_data%ts+tmp_nts) < time_since_reference_point + dt_3d )
             IF ( inflow_data%ts + tmp_nts  >= UBOUND( inflow_data%time, 1 ) )  EXIT
             tmp_nts = tmp_nts + 1
          ENDDO
!
!--       Deallocate arrays and allocate them again with updated lower/upper bounds for the
!--       time dimension and possibly different dimension size if only a few timesteps are left
!--       on file.
          inflow_data%nts = tmp_nts

          ts = inflow_data%ts
          te = inflow_data%ts + inflow_data%nts

          IF ( ALLOCATED( inflow_data%e_l  ) )  DEALLOCATE( inflow_data%e_l  )
          IF ( ALLOCATED( inflow_data%u_l  ) )  DEALLOCATE( inflow_data%u_l  )
          IF ( ALLOCATED( inflow_data%v_l  ) )  DEALLOCATE( inflow_data%v_l  )
          IF ( ALLOCATED( inflow_data%w_l  ) )  DEALLOCATE( inflow_data%w_l  )
          IF ( ALLOCATED( inflow_data%pt_l ) )  DEALLOCATE( inflow_data%pt_l )
          IF ( ALLOCATED( inflow_data%q_l  ) )  DEALLOCATE( inflow_data%q_l  )

          ALLOCATE( inflow_data%e_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je) )
          ALLOCATE( inflow_data%u_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je) )
          ALLOCATE( inflow_data%v_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je) )
          ALLOCATE( inflow_data%w_l(ts:te,nzb+1:nzt-1,inflow_data%js:inflow_data%je) )
          IF ( .NOT. neutral )  THEN
             ALLOCATE( inflow_data%pt_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je) )
          ENDIF
          IF ( humidity )  THEN
             ALLOCATE( inflow_data%q_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je)  )
          ENDIF
!
!--       Read the new boundary data.
          CALL turbulent_inflow_rff_input
       ENDIF
!
!--    Determine current time indices that refer to the time before and after the actual (relative)
!--    model time.
       inflow_data%tp = MINLOC( ABS( inflow_data%time - ( time_since_reference_point + dt_3d ) ),  &
                                DIM = 1 ) - 1
!
!--    Future time index must refer to a time larger/equal the actual model time + model timestep.
!--    Loop is required here because the time resolution of the precursor data might be
!--    temporarily higher than the actual model timestep.
       DO WHILE ( inflow_data%time(inflow_data%tp) < time_since_reference_point + dt_3d )
          IF ( inflow_data%tp >= UBOUND( inflow_data%time, 1 ) )  EXIT
          inflow_data%tp = inflow_data%tp + 1
       ENDDO
!
!--    Previous time index must refer to a time smaller/equal than the interpolated time.
       inflow_data%t  = inflow_data%tp
       IF ( inflow_data%time(inflow_data%t) > time_since_reference_point + dt_3d )  THEN
          inflow_data%t = inflow_data%t - 1
       ENDIF
!
!--    Make sure that the determined timestep indices also cover the range of allocated data.
!--    In the special case that e.g. only one data layer is loaded at the end, t and tp can
!--    easily being out of the allocated array sizes.
       inflow_data%tp = MIN( inflow_data%tp, UBOUND( inflow_data%u_l, 1 ) )
       inflow_data%t  = MAX( inflow_data%t,  LBOUND( inflow_data%u_l, 1 ) )
!
!--    Determine interpolation factor and limit it to 1.
!--    This is because t+dt can slightly exceed time(t+1) before boundary data is updated again.
       IF ( inflow_data%time(inflow_data%tp) - inflow_data%time(inflow_data%t) > 0.0_wp )  THEN
          fac_dt = ( time_since_reference_point + dt_3d - inflow_data%time(inflow_data%t) ) /      &
                   ( inflow_data%time(inflow_data%tp) - inflow_data%time(inflow_data%t) )
       ELSE
          fac_dt = 0.0_wp
       ENDIF

       fac_dt = MIN( 1.0_wp, fac_dt )
!
!--    Precalculate boundary data, i.e. linearly interpolate in time.
       IF ( bc_dirichlet_l )  THEN
!
!--       e, u, v
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                inflow_data%e_in(k,j) = interpolate_linear( inflow_data%e_l(inflow_data%t,k,j),    &
                                                            inflow_data%e_l(inflow_data%tp,k,j),   &
                                                            fac_dt )
                inflow_data%u_in(k,j) = interpolate_linear( inflow_data%u_l(inflow_data%t,k,j),    &
                                                            inflow_data%u_l(inflow_data%tp,k,j),   &
                                                            fac_dt )
                inflow_data%v_in(k,j) = interpolate_linear( inflow_data%v_l(inflow_data%t,k,j),    &
                                                            inflow_data%v_l(inflow_data%tp,k,j),   &
                                                            fac_dt )
             ENDDO
          ENDDO
!
!--       w
          DO  j = nys, nyn
             DO  k = nzb+1, nzt-1
                inflow_data%w_in(k,j) = interpolate_linear( inflow_data%w_l(inflow_data%t,k,j),    &
                                                            inflow_data%w_l(inflow_data%tp,k,j),   &
                                                            fac_dt )
             ENDDO
          ENDDO
!
!--       theta
          IF ( .NOT. neutral )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   inflow_data%pt_in(k,j) =                                                        &
                                        interpolate_linear( inflow_data%pt_l(inflow_data%t,k,j),   &
                                                            inflow_data%pt_l(inflow_data%tp,k,j),  &
                                                            fac_dt )
                ENDDO
             ENDDO
          ENDIF
!
!--       q
          IF ( humidity )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   inflow_data%q_in(k,j) =                                                         &
                                        interpolate_linear( inflow_data%q_l(inflow_data%t,k,j),    &
                                                            inflow_data%q_l(inflow_data%tp,k,j),   &
                                                            fac_dt )
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDIF
!
!-- Impose boundary data onto the prognostic arrays. This is done every RK3-substep.
    IF ( bc_dirichlet_l )  THEN
       i = nxl - 1
       iu = nxlu - 1
!
!--    e, u, v
       DO  j = nys, nyn
          DO  k = nzb+1, nzt
             e(k,j,i)    = inflow_data%e_in(k,j) *                                                 &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             u(k,j,iu)   = inflow_data%u_in(k,j) *                                                 &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,iu), 1 ) )
             u(k,j,iu-1) = u(k,j,iu)
             v(k,j,i)    = inflow_data%v_in(k,j) *                                                 &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO
!
!--       For sake of consistency, set also top boundary condition.
          e(nzt+1,j,i)  = e(nzt,j,i)
          u(nzt+1,j,iu) = u(nzt,j,iu)
          v(nzt+1,j,i)  = v(nzt,j,i)
       ENDDO
!
!--    w
       DO  j = nys, nyn
          DO  k = nzb+1, nzt-1
             w(k,j,i) = inflow_data%w_in(k,j) *                                                    &
                        MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
          ENDDO
!
!--       For sake of consistency, set also top boundary condition.
          w(nzt:nzt+1,j,i) = w(nzt-1,j,i)
       ENDDO
!
!--    theta
       IF ( .NOT. neutral )  THEN
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                pt(k,j,i) = inflow_data%pt_in(k,j)
             ENDDO
!
!--          For sake of consistency, set also top boundary condition.
             pt(nzt+1,j,i) = pt(nzt,j,i)
          ENDDO
       ENDIF
!
!--    q
       IF ( humidity )  THEN
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                q(k,j,i) = inflow_data%q_in(k,j)
             ENDDO
!
!--          For sake of consistency, set also top boundary condition.
             q(nzt+1,j,i) = q(nzt,j,i)
          ENDDO
       ENDIF
    ENDIF

    CALL cpu_log( log_point_s(66), 'turbulent_inflow_rff', 'stop' )

 END SUBROUTINE turbulent_inflow_rff_bc


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize turbulent inflow method where data is read from the dynamic input file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_rff_init

    CHARACTER(LEN=500) ::  variable_not_found = ''              !< string to gather information of missing
                                                                !< dimensions and variables
    CHARACTER(LEN=500), DIMENSION(:), ALLOCATABLE ::  var_names !< list of variable names in dynamic input file

    INTEGER(iwp) ::  i  !< grid index in x-direction
    INTEGER(iwp) ::  iu !< grid index in x-direction on u-grid
    INTEGER(iwp) ::  j  !< grid index in y-direction
    INTEGER(iwp) ::  k  !< grid index in z-direction
    INTEGER(iwp) ::  te !< lower index of time dimension
    INTEGER(iwp) ::  ts !< upper index of time dimension

    LOGICAL ::  trigger_error_message !< control flag to check whether all required input data is available


    trigger_error_message = .FALSE.
    variable_not_found = ''

#if defined( __netcdf )
!
!-- Open the file and read metadata.
    CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ), pids_id )

    CALL inquire_num_variables( pids_id, num_var_pids )
!
!-- Allocate memory to store variable names and read them.
    ALLOCATE( var_names(1:num_var_pids) )
    CALL inquire_variable_names( pids_id, var_names )
!
!-- Read time dimension, allocate memory and finally read time array. Further, trigger
!-- error message in case the dimension is not given on file.
    IF ( check_existence( var_names, inflow_data%char_time ) )  THEN
       CALL get_attribute( pids_id, char_fill, inflow_data%fill_time, .FALSE.,                     &
                           inflow_data%char_time, ignore_error = .TRUE. )
       CALL get_dimension_length( pids_id, inflow_data%nt, inflow_data%char_time )
       ALLOCATE( inflow_data%time(0:inflow_data%nt-1) )
       CALL get_variable( pids_id, inflow_data%char_time, inflow_data%time )
    ELSE
       trigger_error_message = .TRUE.
       variable_not_found = TRIM( inflow_data%char_time ) // ', '
    ENDIF
!
!-- Read horizontal dimension as well as vertical dimension of scalar und w grid.
!-- Further, trigger error message in case the dimension is not given on file.
    IF ( check_existence( var_names, inflow_data%char_y ) )  THEN
       CALL get_dimension_length( pids_id, inflow_data%ny, inflow_data%char_y )
       ALLOCATE( inflow_data%y(1:inflow_data%ny) )
       CALL get_variable( pids_id, inflow_data%char_y, inflow_data%y )
    ELSE
       trigger_error_message = .TRUE.
       variable_not_found = TRIM( variable_not_found ) // inflow_data%char_y // ', '
    ENDIF

    IF ( check_existence( var_names, inflow_data%char_zu ) )  THEN
       CALL get_dimension_length( pids_id, inflow_data%nzu, inflow_data%char_zu )
       ALLOCATE( inflow_data%zu(1:inflow_data%nzu) )
       CALL get_variable( pids_id, inflow_data%char_zu, inflow_data%zu )
    ELSE
       trigger_error_message = .TRUE.
       variable_not_found = TRIM( variable_not_found ) // inflow_data%char_zu // ', '
    ENDIF

    IF ( check_existence( var_names, inflow_data%char_zw ) )  THEN
       CALL get_dimension_length( pids_id, inflow_data%nzw, inflow_data%char_zw )
       ALLOCATE( inflow_data%zw(1:inflow_data%nzw) )
       CALL get_variable( pids_id, inflow_data%char_zw, inflow_data%zw )
    ELSE
       trigger_error_message = .TRUE.
       variable_not_found = TRIM( variable_not_found ) // inflow_data%char_zw // ', '
    ENDIF
!
!-- Close input file.
    CALL close_input_file( pids_id )
#endif

    IF ( trigger_error_message )  THEN
       message_string = 'turbulent_inflow_method = "' // TRIM( turbulent_inflow_method ) //       &
                        '" - dimension(s): &"' //                                                  &
                        TRIM( variable_not_found ) // '" not found in dynamic driver'
       CALL message( 'turbulent_inflow_init', 'TUI0011', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Now, check if all required input variables are available.
    trigger_error_message = .FALSE.
    variable_not_found = ''

    IF ( .NOT. check_existence( var_names, inflow_data%char_e ) )  THEN
       trigger_error_message = .TRUE.
       variable_not_found = inflow_data%char_e // ', '
    ENDIF

    IF ( .NOT. check_existence( var_names, inflow_data%char_u ) )  THEN
       trigger_error_message = .TRUE.
       variable_not_found = TRIM( variable_not_found ) // inflow_data%char_u // ', '
    ENDIF

    IF ( .NOT. check_existence( var_names, inflow_data%char_v ) )  THEN
       trigger_error_message = .TRUE.
       variable_not_found = TRIM( variable_not_found ) // inflow_data%char_v // ', '
    ENDIF

    IF ( .NOT. check_existence( var_names, inflow_data%char_w ) )  THEN
       trigger_error_message = .TRUE.
       variable_not_found = TRIM( variable_not_found ) // inflow_data%char_w // ', '
    ENDIF

    IF ( .NOT. check_existence( var_names, inflow_data%char_pt )  .AND.  .NOT. neutral )  THEN
       trigger_error_message = .TRUE.
       variable_not_found = TRIM( variable_not_found ) // inflow_data%char_pt // ', '
    ENDIF

    IF ( .NOT. check_existence( var_names, inflow_data%char_q )  .AND.  humidity )  THEN
       trigger_error_message = .TRUE.
       variable_not_found = TRIM( variable_not_found ) // inflow_data%char_q // ', '
    ENDIF

    IF ( trigger_error_message )  THEN
       message_string = 'turbulent_inflow_method = "' // TRIM( turbulent_inflow_method ) //        &
                        '" - variable(s): &"' //                                                   &
                        TRIM( variable_not_found ) // '" not found in dynamic input file.'
       CALL message( 'turbulent_inflow_init', 'TUI0012', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- After variable checks, further check if dimensions are defined properly.
    IF ( inflow_data%ny-1 /= ny )  THEN
       message_string = 'number of y-grid points in dynamic driver does not match ' //             &
                        'the number of numeric grid points'
       CALL message( 'turbulent_inflow_init', 'TUI0013', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( inflow_data%nzu /= nz )  THEN
       message_string = 'number of zu-grid points in dynamic driver does not match ' //            &
                        'the number of numeric grid points'
       CALL message( 'turbulent_inflow_init', 'TUI0014', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( inflow_data%nzw /= nz-1 )  THEN
       message_string = 'number of zw-grid points in dynamic driver does not match ' //            &
                        'the number of numeric grid points'
       CALL message( 'turbulent_inflow_init', 'TUI0015', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check if sufficient timesteps have been provided in the dynamic input file.
!-- Note, if surface spinup is employed, the spinup time must be subtracted from end_time!
    IF ( inflow_data%time(inflow_data%nt-1) < end_time - spinup_time )  THEN
       WRITE( message_string, *) 'dynamic driver provides too few time levels for ',               &
                                 'turbulent inflow: &', inflow_data%time(inflow_data%nt-1),        &
                                 ' is less than ', end_time - spinup_time
       CALL message( 'turbulent_inflow_init', 'TUI0016', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check if time starts at zero (plus/minus a small offset).
    IF ( ABS( inflow_data%time(0) ) > 10E-5_wp )  THEN
       message_string = 'dimension "time_inflow" in dynamic driver must start at 0.0s'
       CALL message( 'turbulent_inflow_init', 'TUI0017', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check if time dimension contains any _FillValues.
    IF ( ANY( inflow_data%time == inflow_data%fill_time ) )  THEN
       message_string = 'dimension "time_inflow" in dynamic driver contains _FillValues'
       CALL message( 'turbulent_inflow_init', 'TUI0018', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Obtain time index for current point in time. Note, the time coordinate in the input file is
!-- always relative to zero. Further, since time_since_reference_point is negativ when spinup
!-- is applied, use MAX function to obtain correct time index.
    inflow_data%ts = MINLOC( ABS( inflow_data%time - MAX( time_since_reference_point, 0.0_wp ) ),  &
                             DIM = 1 ) - 1
!
!-- Note, in case of restart runs, the time index for the boundary data may indicate a time in
!-- the future. This needs to be checked and corrected.
    IF ( TRIM( initializing_actions ) == 'read_restart_data'  .AND.                                &
         inflow_data%time(inflow_data%ts) > time_since_reference_point )  THEN
       inflow_data%ts = inflow_data%ts - 1
    ENDIF
!
!-- Check input block size for too large or negative values. For example, input block size should
!-- not exceed 2x the subdomain size in x-direction, so that per input call the maximum data size
!-- is limited to two additional 3D arrays per variable.
    IF ( input_block_size < 2  .OR.  input_block_size > 2 * ( nxr - nxl + 1 ) )  THEN
       WRITE( message_string, *) 'input_block_size must be smaller or equal than 2 x the ',        &
                                 'subdomain size in x-direction, which is: ', nxr - nxl + 1,       &
                                 '&and larger or equal than 2'
       CALL message( 'turbulent_inflow_init', 'TUI0019', 1, 2, 0, 6, 0 )
    ENDIF
    inflow_data%nts = input_block_size

!
!-- Determine the block size, i.e. the number of time layers to be read. This must not be
!-- larger than the total number of time layers as well as the number of remaining time layers.
!-- The storage of multiple time layers avoids that NetCDF IO is required every LES timestep.
    inflow_data%nts = MIN( inflow_data%nts, inflow_data%nt - 1 )
    inflow_data%nts = MIN( inflow_data%nts, inflow_data%nt - 1 - inflow_data%ts )
!
!-- Allocate memory for the variables. Since also inner subdomains participate in the
!-- NetCDF collective read, allocate dummy arrays with size 0. Furthermore, boundary arrays
!-- are defined over ghost points, in order to avoid additional ghost point exchange. An exception
!-- is only made for the left-south and left-north boundary core, where the lowest and largest
!-- index are limited to 0 and ny, respectively.
    inflow_data%js = MAX( MERGE( nysg, 1, bc_dirichlet_l ), 0  )
    inflow_data%je = MIN( MERGE( nyng, 1, bc_dirichlet_l ), ny )
!
!-- Determine lower and upper allocation bounds for time dimension.
    ts = inflow_data%ts
    te = inflow_data%ts + inflow_data%nts

    ALLOCATE( inflow_data%e_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je) )
    ALLOCATE( inflow_data%u_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je) )
    ALLOCATE( inflow_data%v_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je) )
    ALLOCATE( inflow_data%w_l(ts:te,nzb+1:nzt-1,inflow_data%js:inflow_data%je) )
    IF ( .NOT. neutral ) ALLOCATE( inflow_data%pt_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je) )
    IF ( humidity )      ALLOCATE( inflow_data%q_l(ts:te,nzb+1:nzt,inflow_data%js:inflow_data%je)  )
!
!-- Read first block of time layers.
    CALL turbulent_inflow_rff_input
!
!-- Allocate arrays where linearly interpolated data is stored on.
    ALLOCATE( inflow_data%e_in(nzb+1:nzt,inflow_data%js:inflow_data%je) )
    ALLOCATE( inflow_data%u_in(nzb+1:nzt,inflow_data%js:inflow_data%je) )
    ALLOCATE( inflow_data%v_in(nzb+1:nzt,inflow_data%js:inflow_data%je) )
    ALLOCATE( inflow_data%w_in(nzb+1:nzt-1,inflow_data%js:inflow_data%je) )
    IF ( .NOT. neutral )  ALLOCATE( inflow_data%pt_in(nzb+1:nzt,inflow_data%js:inflow_data%je) )
    IF ( humidity      )  ALLOCATE( inflow_data%q_in(nzb+1:nzt,inflow_data%js:inflow_data%je)  )
!
!-- Finally, set the boundary values. Also store them on internal arrays. Later-on, these are used
!-- to set the boundary values during the time-integration. This detour is necessary to have
!-- full control over the boundary values. Else, for some variables, Neumann conditions would be
!-- set.
    IF ( bc_dirichlet_l )  THEN
       i = nxl - 1
       iu = nxlu - 1
!
!--    Loop over j corresponds to j = MAX( 0, nysg ), MIN( nyn, ny )
       DO  j = inflow_data%js, inflow_data%je
          DO  k = nzb+1, nzt
             e(k,j,i) =    inflow_data%e_l(inflow_data%ts,k,j) *                                   &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
             u(k,j,iu) =   inflow_data%u_l(inflow_data%ts,k,j) *                                   &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,iu), 1 ) )
             u(k,j,iu-1) = u(k,j,iu)
             v(k,j,i) =    inflow_data%v_l(inflow_data%ts,k,j) *                                   &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )

             inflow_data%e_in(k,j) = e(k,j,i)
             inflow_data%u_in(k,j) = u(k,j,iu)
             inflow_data%v_in(k,j) = v(k,j,i)
          ENDDO
!
!--       For sake of consistency, set also top boundary condition.
          e(nzt+1,j,i)  = e(nzt,j,i)
          u(nzt+1,j,iu) = u(nzt,j,iu)
          v(nzt+1,j,i)  = v(nzt,j,i)
       ENDDO

       DO  j = inflow_data%js, inflow_data%je
          DO  k = nzb+1, nzt-1
             w(k,j,i) = inflow_data%w_l(inflow_data%ts,k,j) *                                      &
                        MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
             inflow_data%w_in(k,j) = w(k,j,i)
          ENDDO
!
!--       For sake of consistency, set also top boundary condition.
          w(nzt:nzt+1,j,i) = w(nzt-1,j,i)
       ENDDO

       IF ( .NOT. neutral )  THEN
          DO  j = inflow_data%js, inflow_data%je
             DO  k = nzb+1, nzt
                pt(k,j,i) = inflow_data%pt_l(inflow_data%ts,k,j) *                                 &
                            MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                inflow_data%pt_in(k,j) = pt(k,j,i)
             ENDDO
!
!--          For sake of consistency, set also top boundary condition.
             pt(nzt+1,j,i) = pt(nzt,j,i)
          ENDDO
       ENDIF

       IF ( humidity )  THEN
          DO  j = inflow_data%js, inflow_data%je
             DO  k = nzb+1, nzt
                q(k,j,i) = inflow_data%q_l(inflow_data%ts,k,j) *                                   &
                           MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
                inflow_data%q_in(k,j) = q(k,j,i)
             ENDDO
!
!--          For sake of consistency, set also top boundary condition.
             q(nzt+1,j,i) = q(nzt,j,i)
          ENDDO
       ENDIF
    ENDIF
!
!-- Take care of mass-flux conservation.
    CALL turbulent_inflow_mass_flux_conservation
!
!-- Initialize prognostic levels as well.
    pt_p = pt
    u_p = u
    v_p = v
    w_p = w
    e_p = e
    IF ( humidity )  q_p = q
    IF ( passive_scalar )  s_p  = s

 END SUBROUTINE turbulent_inflow_rff_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> NetCDF input.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_rff_input

#if defined( __netcdf )
!
!-- Open the file and read metadata.
    CALL open_read_file( TRIM( input_file_dynamic ) // TRIM( coupling_char ), pids_id )
!
!-- Note, even though only lateral data is requied, input data is accessed by parallel IO because
!-- collective parallel access shows better performance than just a conditional access. This is the
!-- reason why different arguments are passed depending on the boundary control flags.
!-- Cores that do not belong to the respective boundary, only perform a dummy read with
!-- count = 0, just in order to participate the collective operation.
    CALL get_variable( pids_id, inflow_data%char_e,                                                &
                       inflow_data%e_l,                                                            &
                       MERGE( inflow_data%js+1, 1, bc_dirichlet_l ),                               &
                       MERGE( nzb+1, 1, bc_dirichlet_l ),                                          &
                       MERGE( inflow_data%ts+1, 1, bc_dirichlet_l ),                               &
                       MERGE( inflow_data%je-inflow_data%js+1, 0, bc_dirichlet_l ),                &
                       MERGE( inflow_data%nzu, 0, bc_dirichlet_l ),                                &
                       MERGE( inflow_data%nts+1, 0, bc_dirichlet_l ),                              &
                       .TRUE. )

    CALL get_variable( pids_id, inflow_data%char_u,                                                &
                       inflow_data%u_l,                                                            &
                       MERGE( inflow_data%js+1, 1, bc_dirichlet_l ),                               &
                       MERGE( nzb+1, 1, bc_dirichlet_l ),                                          &
                       MERGE( inflow_data%ts+1, 1, bc_dirichlet_l ),                               &
                       MERGE( inflow_data%je-inflow_data%js+1, 0, bc_dirichlet_l ),                &
                       MERGE( inflow_data%nzu, 0, bc_dirichlet_l ),                                &
                       MERGE( inflow_data%nts+1, 0, bc_dirichlet_l ),                              &
                       .TRUE. )

    CALL get_variable( pids_id, inflow_data%char_v,                                                &
                       inflow_data%v_l,                                                            &
                       MERGE( inflow_data%js+1, 1, bc_dirichlet_l ),                               &
                       MERGE( nzb+1, 1, bc_dirichlet_l ),                                          &
                       MERGE( inflow_data%ts+1, 1, bc_dirichlet_l ),                               &
                       MERGE( inflow_data%je-inflow_data%js+1, 0, bc_dirichlet_l ),                &
                       MERGE( inflow_data%nzu, 0, bc_dirichlet_l ),                                &
                       MERGE( inflow_data%nts+1, 0, bc_dirichlet_l ),                              &
                       .TRUE. )

    CALL get_variable( pids_id, inflow_data%char_w,                                                &
                       inflow_data%w_l,                                                            &
                       MERGE( inflow_data%js+1, 1, bc_dirichlet_l ),                               &
                       MERGE( nzb+1, 1, bc_dirichlet_l ),                                          &
                       MERGE( inflow_data%ts+1, 1, bc_dirichlet_l ),                               &
                       MERGE( inflow_data%je-inflow_data%js+1, 0, bc_dirichlet_l ),                &
                       MERGE( inflow_data%nzw, 0, bc_dirichlet_l ),                                &
                       MERGE( inflow_data%nts+1, 0, bc_dirichlet_l ),                              &
                       .TRUE. )

    IF ( .NOT. neutral )  THEN
       CALL get_variable( pids_id, inflow_data%char_pt,                                            &
                          inflow_data%pt_l,                                                        &
                          MERGE( inflow_data%js+1, 1, bc_dirichlet_l ),                            &
                          MERGE( nzb+1, 1, bc_dirichlet_l ),                                       &
                          MERGE( inflow_data%ts+1, 1, bc_dirichlet_l ),                            &
                          MERGE( inflow_data%je-inflow_data%js+1, 0, bc_dirichlet_l ),             &
                          MERGE( inflow_data%nzu, 0, bc_dirichlet_l ),                             &
                          MERGE( inflow_data%nts+1, 0, bc_dirichlet_l ),                           &
                          .TRUE. )
    ENDIF

    IF ( humidity )  THEN
       CALL get_variable( pids_id, inflow_data%char_q,                                             &
                          inflow_data%q_l,                                                         &
                          MERGE( inflow_data%js+1, 1, bc_dirichlet_l ),                            &
                          MERGE( nzb+1, 1, bc_dirichlet_l ),                                       &
                          MERGE( inflow_data%ts+1, 1, bc_dirichlet_l ),                            &
                          MERGE( inflow_data%je-inflow_data%js+1, 0, bc_dirichlet_l ),             &
                          MERGE( inflow_data%nzu, 0, bc_dirichlet_l ),                             &
                          MERGE( inflow_data%nts+1, 0, bc_dirichlet_l ),                           &
                          .TRUE. )
    ENDIF
!
!-- Close input file
    CALL close_input_file( pids_id )
#endif

 END SUBROUTINE turbulent_inflow_rff_input


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_rrd_global_ftn( found )

    LOGICAL, INTENT(OUT) ::  found  !< control flag to indicate whether a variable has been found in the module or not


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'inflow_damping_factor' )
          IF ( .NOT. ALLOCATED( inflow_damping_factor ) )  THEN
             ALLOCATE( inflow_damping_factor(0:nz+1) )
          ENDIF
          READ( 13 )  inflow_damping_factor

       CASE ( 'inflow_damping_height' )
          READ( 13 )  inflow_damping_height

       CASE ( 'inflow_damping_width' )
          READ( 13 )  inflow_damping_width

       CASE ( 'recycling_width' )
          READ( 13 )  recycling_width

       CASE ( 'turbulent_inflow_method' )
          READ( 13 )  turbulent_inflow_method

       CASE DEFAULT
          found = .FALSE.

    END SELECT

 END SUBROUTINE turbulent_inflow_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_rrd_global_mpi

    LOGICAL ::  array_found  !< flag to indicate that an array which is searched for exists


    CALL rd_mpi_io_check_array( 'inflow_damping_factor', found = array_found )

    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( inflow_damping_factor ) )  THEN
          ALLOCATE( inflow_damping_factor(0:nz+1) )
       ENDIF
       CALL rrd_mpi_io_global_array( 'inflow_damping_factor', inflow_damping_factor )
    ENDIF
    CALL rrd_mpi_io( 'inflow_damping_height',   inflow_damping_height      )
    CALL rrd_mpi_io( 'inflow_damping_width',    inflow_damping_width       )
    CALL rrd_mpi_io( 'recycling_width',         recycling_width            )
    CALL rrd_mpi_io( 'turbulent_inflow_method', turbulent_inflow_method    )

 END SUBROUTINE turbulent_inflow_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes global restart data
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_wrd_global


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       IF ( ALLOCATED( inflow_damping_factor ) )  THEN
          CALL wrd_write_string( 'inflow_damping_factor' )
          WRITE ( 14 )  inflow_damping_factor
       ENDIF

       CALL wrd_write_string( 'inflow_damping_height' )
       WRITE ( 14 )  inflow_damping_height

       CALL wrd_write_string( 'inflow_damping_width' )
       WRITE ( 14 )  inflow_damping_width

       CALL wrd_write_string( 'recycling_width' )
       WRITE ( 14 )  recycling_width

       CALL wrd_write_string( 'turbulent_inflow_method' )
       WRITE ( 14 )  turbulent_inflow_method

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       IF ( ALLOCATED( inflow_damping_factor ) )  THEN
          CALL wrd_mpi_io_global_array( 'inflow_damping_factor', inflow_damping_factor )
       ENDIF
       CALL wrd_mpi_io( 'inflow_damping_height',   inflow_damping_height   )
       CALL wrd_mpi_io( 'inflow_damping_width',    inflow_damping_width    )
       CALL wrd_mpi_io( 'recycling_width',         recycling_width         )
       CALL wrd_mpi_io( 'turbulent_inflow_method', turbulent_inflow_method )

    ENDIF

 END SUBROUTINE turbulent_inflow_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Performs y-shift of the input data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE turbulent_inflow_y_shift( inflow_dist, inflow_var )

#if defined( __parallel )
    INTEGER(iwp) ::  next     !< ID of receiving PE for y-shift
    INTEGER(iwp) ::  npg      !< number of grid data send during shifting
    INTEGER(iwp) ::  prev     !< ID of sending PE for y-shift
#endif

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng), INTENT(INOUT), OPTIONAL ::  inflow_var  !< inflow variable

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,num_mean_inflow_profiles,nbgp), INTENT(INOUT),         &
              OPTIONAL ::  inflow_dist !< turbulence signal
#if defined( __parallel )
    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  local_var  !< auxiliary variable, used for y-shift
#endif


#if defined( __parallel )
    IF ( ( y_shift /= 0 )  .AND.  myidx == id_inflow )  THEN
!
!--    Calculate the ID of the PE which sends data to this PE (prev) and of the PE which receives
!--    data from this PE (next).
       prev = MODULO( myidy - y_shift, npey )
       next = MODULO( myidy + y_shift, npey )
!
!--    For turbulence recycling method, shift inflow_dist in positive y direction by a number of PEs
!--    equal to y_shift.
       IF ( turbulent_inflow_rec )  THEN
          ALLOCATE( local_var(nzb:nzt+1,nysg:nyng,num_mean_inflow_profiles,nbgp) )
          local_var = 0.0_wp
!
!--       Determine number of grid data send during shifting.
          npg = ( nzt - nzb + 2 ) * num_mean_inflow_profiles * nbgp * ( nyn - nys + 1 + 2 * nbgp )

          CALL MPI_SENDRECV( inflow_dist(nzb,nysg,1,1), npg, MPI_REAL, next, 1,                    &
                             local_var(nzb,nysg,1,1),   npg, MPI_REAL, prev, 1, comm1dy,           &
                             status, ierr )

          inflow_dist = local_var

          DEALLOCATE( local_var )
!
!--    For the read-from file inflow method, shift only a single variable along y-direction.
       ELSEIF ( turbulent_inflow_rff )  THEN

          ALLOCATE( local_var(nzb:nzt+1,nysg:nyng,1:1,1:1) )
          local_var = 0.0_wp
!
!--       Determine number of grid data send during shifting.
          npg = ( nzt - nzb + 2 ) * ( nyn - nys + 1 + 2 * nbgp )

          CALL MPI_SENDRECV( inflow_var(nzb,nysg),    npg, MPI_REAL, next, 1,                      &
                             local_var(nzb,nysg,1,1), npg, MPI_REAL, prev, 1, comm1dy,             &
                             status, ierr )

          inflow_var(nzb:nzt+1,nysg:nyng) = local_var(nzb:nzt+1,nysg:nyng,1,1)

          DEALLOCATE( local_var )

       ENDIF

    ENDIF

#endif
!
!-- Dummy check to avoid compiler errors concerning unused variables in serial mode.
    IF ( PRESENT( inflow_dist )  .OR.  .NOT.  PRESENT( inflow_dist ) )  CONTINUE
    IF ( PRESENT( inflow_var )   .OR.  .NOT.  PRESENT( inflow_var  ) )  CONTINUE

 END SUBROUTINE turbulent_inflow_y_shift

 END MODULE turbulent_inflow_mod
