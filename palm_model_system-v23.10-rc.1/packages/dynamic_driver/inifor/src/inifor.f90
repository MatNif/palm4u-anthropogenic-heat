!> @file src/inifor.f90
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
!> @author Eckhard Kadasch, (Deutscher Wetterdienst, Offenbach)
!
! Description:
! ------------
!> INIFOR is an interpolation tool for generating meteorological initialization
!> and forcing data for the urban climate model PALM-4U. The required 
!> meteorological fields are interpolated from output data of the mesoscale
!> model COSMO. This is the main program file.
!------------------------------------------------------------------------------!
 PROGRAM inifor

#if defined ( __netcdf )

    USE inifor_control
    USE inifor_defs
    USE inifor_grid,                                                           &
        ONLY:  averaging_width_ns,                                             &
               averaging_width_ew,                                             &
               cfg,                                                            &    
               cosmo_grid,                                                     &
               f3,                                                             &
               fini_grids,                                                     &
               fini_io_groups,                                                 &
               fini_variables,                                                 &
               fini_file_lists,                                                &
               io_group_list,                                                  &
               lam_centre,                                                     &
               lambda_n,                                                       &
               ls_forcing_variables_required,                                  &
               nx, ny, nz,                                                     &
               origin_lat,                                                     &
               origin_lon,                                                     &
               output_file,                                                    &
               p0,                                                             &
               phi_centre,                                                     &
               phi_n,                                                          &
               preprocess,                                                     &
               palm_grid,                                                      &
               setup_grids,                                                    &
               setup_parameters,                                               &
               setup_variable_tables,                                          &
               setup_io_groups
    USE inifor_io
    USE inifor_transform,                                                      &
        ONLY:  average_pressure_perturbation,                                  &
               average_profile,                                                & 
               geostrophic_winds,                                              &
               get_surface_pressure,                                           &
               interp_average_profile,                                         &
               interpolate_1d,                                                 &
               interpolate_1d_arr,                                             & 
               interpolate_2d,                                                 &
               interpolate_3d
    USE inifor_types
    
    IMPLICIT NONE
    
    INTEGER(iwp) ::  igroup !< loop index for IO groups
    INTEGER(iwp) ::  ivar   !< loop index for output variables
    INTEGER(iwp) ::  iter   !< loop index for time steps

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)     ::  output_arr !< array buffer for interpolated quantities
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  rho_centre_cosmo !< density profile of the centre averaging domain
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  ug_cosmo   !< profile of the geostrophic wind in x direction on COSMO levels
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  vg_cosmo   !< profile of the geostrophic wind in y direction on COSMO levels
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  ug_palm    !< profile of the geostrophic wind in x direction interpolated onto PALM levels
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  vg_palm    !< profile of the geostrophic wind in y direction interpolated onto PALM levels
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  rho_north_cosmo  !< density profile of the northern averaging domain
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  rho_south_cosmo  !< density profile of the southern averaging domain
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  rho_east_cosmo   !< density profile of the eastern averaging domain
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  rho_west_cosmo   !< density profile of the western averaging domain
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  p_north_cosmo    !< pressure profile of the northern averaging domain
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  p_south_cosmo    !< pressure profile of the southern averaging domain
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  p_east_cosmo     !< pressure profile of the eastern averaging domain
    REAL(wp), ALLOCATABLE, DIMENSION(:), TARGET ::  p_west_cosmo     !< pressure profile of the western averaging domain

    REAL(wp), POINTER, DIMENSION(:) ::  internal_cosmo_arr !< pointer to the currently processed internal array (density, pressure)
    REAL(wp), POINTER, DIMENSION(:) ::  ug_vg_palm   !< pointer to the currently processed geostrophic wind component

    TYPE(nc_var), POINTER        ::  output_var      !< pointer to the currently processed output variable
    TYPE(io_group), POINTER      ::  group           !< pointer to the currently processed IO group
    TYPE(container), ALLOCATABLE ::  input_buffer(:) !< buffer of the current IO group's input arrays

    LOGICAL, SAVE ::  ug_vg_have_been_computed = .FALSE. !< flag for managing geostrophic wind allocation and computation
    !LOGICAL, SAVE ::  debugging_output = .TRUE.          !< flag controllging output of internal variables
    
!> \mainpage About INIFOR
!>  ...
!
!------------------------------------------------------------------------------
!- Section 1: Initialization
!------------------------------------------------------------------------------
    CALL log_runtime( 'init', 'void' )

!
!-- Initialize INIFOR's parameters from command-line interface and namelists
    CALL setup_parameters

!
!-- Initialize all grids, including interpolation neighbours and weights
    CALL setup_grids
    CALL log_runtime( 'time', 'init' )

!
!-- Initialize the netCDF output file and define dimensions
    CALL setup_netcdf_dimensions( output_file, palm_grid, cfg%start_date,      &
                                  origin_lon, origin_lat )
    CALL log_runtime( 'time', 'write' )

!
!-- Set up the tables containing the input and output variables and set
!-- the corresponding netCDF dimensions for each output variable
    CALL setup_variable_tables
    CALL setup_io_groups
    CALL log_runtime( 'time', 'init' )
!
!-- Add the output variables to the netCDF output file
    CALL setup_netcdf_variables( output_file%name, io_group_list)
    CALL log_runtime( 'time', 'write' )

!------------------------------------------------------------------------------
!-- Section 2: Main loop
!------------------------------------------------------------------------------
!
!-- Input and output variables are organized into IO groups. For instance, the
!-- 'thermodynamics' IO group bundles the input variaebls T, P, QV and the
!-- output variables p, theta, rho, and qv.
!-- An IO group bunldes variables that are physically dependent on each other.
!-- In case of the 'thermodynamics' group, theta = f(P,T), rho = f(P,T,QV).
    DO  igroup = 1, SIZE( io_group_list )

       group => io_group_list(igroup)
       IF ( group%to_be_processed )  THEN
          
!--       Loop over all output time steps for the current group.
          DO  iter = 1, group%nt 

!------------------------------------------------------------------------------
!-- Section 2.1: Read and preprocess input data
!------------------------------------------------------------------------------
             CALL read_input_variables( group, iter, input_buffer )
             CALL log_runtime( 'time', 'read' )

!--          Carry out all required physical conversion of the input variables
!--          of the current IO group on the input (COSMO) grid. For instance,
!--          horizontal velocities are rotated to the PALM coordinate system and
!--          potential temperature is computed from the absolute temperature and
!--          pressure.
             CALL preprocess( group, input_buffer, cosmo_grid )
             CALL log_runtime( 'time', 'comp' )

             !TODO: move this assertion into 'preprocess'.
             IF ( .NOT. ALL(input_buffer(:)%is_preprocessed .AND. .TRUE.) )  THEN
                message = "Input buffers for group '" // TRIM( group%kind ) // &
                          "' could not be preprocessed sucessfully."
                CALL inifor_abort( 'main loop', message )
             ENDIF

!------------------------------------------------------------------------------
!-- Section 2.2: Interpolate each output variable of the group
!------------------------------------------------------------------------------
             DO  ivar = 1, group%nv

                output_var => group%out_vars(ivar)

                IF ( output_var%to_be_processed .AND.                          &
                     iter .LE. output_var%nt )  THEN

                   message = "Processing '" // TRIM( output_var%name ) //      &
                             "' (" // TRIM( output_var%kind ) //               &
                             "), iteration " // TRIM( str( iter ) ) //" of " //&
                             TRIM( str( output_var%nt ) )
                   CALL report( 'main loop', message )

                   SELECT CASE( TRIM( output_var%task ) )

!--                   2D horizontal interpolation
                      CASE( 'interpolate_2d' ) 
                      
                         SELECT CASE( TRIM( output_var%kind ) )
                          
                         CASE( 'init soil' )
   
                            ALLOCATE( output_arr(0:output_var%grid%nx,         &
                                                 0:output_var%grid%ny,         &
                                                 SIZE( output_var%grid%depths )) )
   
                         CASE ( 'surface forcing' )
   
                            ALLOCATE( output_arr(0:output_var%grid%nx,         &
                                                 0:output_var%grid%ny, 1) )
   
                         CASE DEFAULT
   
                             message = "'" // TRIM( output_var%kind ) // "' is not a soil variable"
                             CALL inifor_abort( "main loop", message )
   
                         END SELECT
                         CALL log_runtime( 'time', 'alloc' )
   
                         CALL interpolate_2d( input_buffer(output_var%input_id)%array(:,:,:), &
                                 output_arr(:,:,:), output_var%intermediate_grid, output_var )
                         CALL log_runtime( 'time', 'comp' )
   
   
!--                   Interpolation in 3D, used for atmospheric initial and
!--                   boundary conditions.
                      CASE ( 'interpolate_3d' )
   
                         ALLOCATE( output_arr(0:output_var%grid % nx,          &
                                              0:output_var%grid % ny,          &
                                              1:output_var%grid % nz) )
   
                         CALL log_runtime( 'time', 'alloc' )
                         CALL interpolate_3d(                                  &
                            input_buffer(output_var%input_id)%array(:,:,:),    &
                            output_arr(:,:,:),                                 &
                            output_var%intermediate_grid,                      &
                            output_var%grid)
                         CALL log_runtime( 'time', 'comp' )
   
!--                   Compute initial avaerage profiles (if --init-mode profile
!--                   is used) 
                      CASE ( 'average profile' )
   
                         ALLOCATE( output_arr(0:0, 0:0, 1:output_var%averaging_grid%nz) )
                         CALL log_runtime( 'time', 'alloc' )
                         CALL interp_average_profile(                          &
                            input_buffer(output_var%input_id)%array(:,:,:),    &
                            output_arr(0,0,:),                                 &
                            output_var%averaging_grid )
   
                         IF ( TRIM( output_var%name ) ==                       &
                              'surface_forcing_surface_pressure' )  THEN
   
                            IF ( cfg%p0_is_set )  THEN
                               output_arr(0,0,1) = p0
                            ELSE
                               CALL get_surface_pressure(                      &
                                  output_arr(0,0,:),                           &
                                  rho_centre_cosmo,                            &
                                  output_var%averaging_grid )
                            ENDIF
   
                         ENDIF
                         CALL log_runtime( 'time', 'comp' )

                      CASE ( 'average levels' )

                         ALLOCATE( output_arr(0:0, 0:0, 1:output_var%averaging_grid%nz) )
                         CALL log_runtime( 'time', 'alloc' )
                         CALL average_profile(                                 &
                            input_buffer(output_var%input_id)%array(:,:,:),    &
                            output_arr(0,0,:),                                 &
                            output_var%averaging_grid                          &
                         )
                         CALL log_runtime( 'time', 'comp' )
   
!--                   Compute internal profiles, required for differentiation of
!--                   geostrophic wind
                      CASE ( 'internal profile' )
   
                         message = "Averaging of internal profile for variable '" //&
                            TRIM( output_var%name ) // "' is not supported."
   
                         SELECT CASE ( TRIM( output_var%name ) )
   
                         CASE( 'internal_density_centre' )
                            ALLOCATE( rho_centre_cosmo(1:cosmo_grid%nz) )
                            internal_cosmo_arr => rho_centre_cosmo
   
                         CASE( 'internal_density_north' )
                            ALLOCATE( rho_north_cosmo(1:cosmo_grid%nz) )
                            internal_cosmo_arr => rho_north_cosmo
   
                         CASE( 'internal_density_south' )
                            ALLOCATE( rho_south_cosmo(1:cosmo_grid%nz) )
                            internal_cosmo_arr => rho_south_cosmo
   
                         CASE( 'internal_density_east' )
                            ALLOCATE( rho_east_cosmo(1:cosmo_grid%nz) )
                            internal_cosmo_arr => rho_east_cosmo
   
                         CASE( 'internal_density_west' )
                            ALLOCATE( rho_west_cosmo(1:cosmo_grid%nz) )
                            internal_cosmo_arr => rho_west_cosmo
   
                         CASE( 'internal_pressure_north' )
                            ALLOCATE( p_north_cosmo(1:cosmo_grid%nz) )
                            internal_cosmo_arr => p_north_cosmo
   
                         CASE( 'internal_pressure_south' )
                            ALLOCATE( p_south_cosmo(1:cosmo_grid%nz) )
                            internal_cosmo_arr => p_south_cosmo
   
                         CASE( 'internal_pressure_east' )
                            ALLOCATE( p_east_cosmo(1:cosmo_grid%nz) )
                            internal_cosmo_arr => p_east_cosmo
   
                         CASE( 'internal_pressure_west' )
                            ALLOCATE( p_west_cosmo(1:cosmo_grid%nz) )
                            internal_cosmo_arr => p_west_cosmo
   
                         CASE DEFAULT
                            CALL inifor_abort( 'main loop', message )
   
                         END SELECT
                         CALL log_runtime( 'time', 'alloc' )
   
   
                         SELECT CASE( TRIM( output_var%name ) )
   
                         CASE( 'internal_pressure_north',                      &
                               'internal_pressure_south',                      &
                               'internal_pressure_east',                       &
                               'internal_pressure_west' )
   
                            CALL average_pressure_perturbation(                &
                               input_buffer(output_var%input_id)%array(:,:,:), &
                               internal_cosmo_arr(:),                          &
                               cosmo_grid, output_var%averaging_grid           &
                            )
   
                         CASE DEFAULT
   
                            CALL average_profile(                              &
                               input_buffer(output_var%input_id)%array(:,:,:), &
                               internal_cosmo_arr(:),                          &
                               output_var%averaging_grid                       &
                            )

                      END SELECT

                      CALL log_runtime( 'time', 'comp' )

!
!--                This case gets called twice, the first time for ug, the
!--                second time for vg. We compute ug and vg at the first call
!--                and keep both of them around for the second call.
                   CASE ( 'geostrophic winds' )

                      IF (.NOT. ug_vg_have_been_computed )  THEN
                         ALLOCATE( ug_palm(output_var%grid%nz) )
                         ALLOCATE( vg_palm(output_var%grid%nz) )
                         ALLOCATE( ug_cosmo(cosmo_grid%nz) )
                         ALLOCATE( vg_cosmo(cosmo_grid%nz) )

                         IF ( cfg%ug_defined_by_user )  THEN
                            ug_palm = cfg%ug
                            vg_palm = cfg%vg
                         ELSE
                            CALL geostrophic_winds(                            &
                               p_north_cosmo, p_south_cosmo, p_east_cosmo,     &
                               p_west_cosmo, rho_centre_cosmo, f3,             &
                               averaging_width_ew,                             &
                               averaging_width_ns,                             &
                               phi_n, lambda_n,                                &
                               phi_centre, lam_centre,                         &
                               ug_cosmo, vg_cosmo                              &
                            )

                            CALL interpolate_1d( ug_cosmo, ug_palm,            &
                                                 output_var%grid )

                            CALL interpolate_1d( vg_cosmo, vg_palm,            &
                                                 output_var%grid )
                         ENDIF

                         ug_vg_have_been_computed = .TRUE.

                      ENDIF

!
!--                   Select output array of current geostrophic wind component
                      SELECT CASE( TRIM( output_var%name ) )
                      CASE ( 'ls_forcing_ug' )
                         ug_vg_palm => ug_palm
                      CASE ( 'ls_forcing_vg' )
                         ug_vg_palm => vg_palm
                      END SELECT

                      ALLOCATE( output_arr(1,1,output_var%grid%nz) )
                      output_arr(1,1,:) = ug_vg_palm(:)

!--                User defined constant profiles
                   CASE ( 'set profile' )
                      
                      ALLOCATE( output_arr(1,1,1:nz) )
                      CALL log_runtime( 'time', 'alloc' )

                      SELECT CASE ( TRIM( output_var%name ) )

                      CASE ( 'nudging_tau' )
                          output_arr(1, 1, :) = NUDGING_TAU

                      CASE DEFAULT
                          message = "'" // TRIM( output_var%name ) //          &
                             "' is not a valid '" // TRIM( output_var%kind ) //&
                             "' variable kind."
                          CALL inifor_abort( 'main loop', message )
                      END SELECT
                      CALL log_runtime( 'time', 'comp' )

                   CASE DEFAULT
                      message = "Processing task '" // TRIM( output_var%task ) //&
                               "' not recognized."
                      CALL inifor_abort( '', message )

                   END SELECT
                   CALL log_runtime( 'time', 'comp' )

!------------------------------------------------------------------------------
!- Section 2.3: Write current time step of current variable
!------------------------------------------------------------------------------
!
!--                Output of geostrophic pressure profiles (with --debug
!--                option) is currently deactivated, since they are now
!--                defined on averaged COSMO levels instead of PALM levels
!--                (requires definiton of COSMO levels in netCDF output.)
                   !IF (.NOT. output_var%is_internal .OR. debugging_output)  THEN

                   IF ( .NOT. output_var%is_internal )  THEN
                      message = "Writing variable '" // TRIM( output_var%name ) // "'."
                      CALL report( 'main loop', message )
                      CALL update_output( output_var, output_arr, iter,        &
                                          output_file, cfg )
                      CALL log_runtime( 'time', 'write' )
                   ENDIF

                   IF ( ALLOCATED( output_arr ) )  DEALLOCATE( output_arr )
                   CALL log_runtime( 'time', 'alloc' )

                ENDIF

!
!--          output variable loop
             ENDDO

             ug_vg_have_been_computed = .FALSE.
             IF ( group%kind == 'thermodynamics' )  THEN
                DEALLOCATE( rho_centre_cosmo )
                IF ( ls_forcing_variables_required                             &
                     .OR. cfg%ug_defined_by_user )  THEN
                   DEALLOCATE( ug_palm )
                   DEALLOCATE( vg_palm )
                   DEALLOCATE( ug_cosmo )
                   DEALLOCATE( vg_cosmo )
                ENDIF
                IF ( .NOT. cfg%ug_defined_by_user )  THEN
                   DEALLOCATE( rho_north_cosmo )
                   DEALLOCATE( rho_south_cosmo )
                   DEALLOCATE( rho_east_cosmo )
                   DEALLOCATE( rho_west_cosmo )
                   DEALLOCATE( p_north_cosmo )
                   DEALLOCATE( p_south_cosmo )
                   DEALLOCATE( p_east_cosmo )
                   DEALLOCATE( p_west_cosmo )
                ENDIF
             ENDIF

!
!--          Keep input buffer around for averaged (radiation) and 
!--          accumulated COSMO quantities (precipitation).
             IF ( group%kind == 'running average' .OR. &
                  group%kind == 'accumulated' )  THEN
             ELSE
                CALL report( 'main loop', 'Deallocating input buffer', cfg%debug )
                DEALLOCATE( input_buffer )
             ENDIF
             CALL log_runtime( 'time', 'alloc' )

!
!--       time steps / input files loop
          ENDDO

          IF ( ALLOCATED( input_buffer ) )  THEN
             CALL report( 'main loop', 'Deallocating input buffer', cfg%debug )
             DEALLOCATE( input_buffer )
          ENDIF
          CALL log_runtime( 'time', 'alloc' )

       ELSE 

          message = "Skipping IO group " // TRIM( str( igroup ) ) // " '" // TRIM( group%kind ) // "'"
          IF ( ALLOCATED( group%in_var_list ) )  THEN
              message = TRIM( message ) // " with input variable '" //         &
              TRIM( group%in_var_list(1)%name ) // "'."
          ENDIF

          CALL report( 'main loop', message, cfg%debug )

!
!--    IO group%to_be_processed conditional
       ENDIF

!
!-- IO groups loop
    ENDDO

!------------------------------------------------------------------------------
!- Section 3: Clean up.
!------------------------------------------------------------------------------
    CALL fini_file_lists
    CALL fini_io_groups
    CALL fini_variables
    !CALL fini_grids
    CALL log_runtime( 'time', 'alloc' )

    CALL report_warnings
    CALL report_success( output_file%name )
    CALL report_runtime
    CALL close_log

#else

    USE inifor_control
    IMPLICIT NONE
    
    message = "INIFOR was compiled without netCDF support, which is required for it to run. "  //     &
              "To use INIFOR, recompile PALM with netCDF support by adding the -D__netcdf " //        &
              "precompiler flag to your .palm.config file."
    CALL inifor_abort( 'main loop', message )
 
#endif

 END PROGRAM inifor
