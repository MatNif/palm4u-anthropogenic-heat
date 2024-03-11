!> @file wind_turbine_model_mod.f90
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
! Copyright 2009-2021 Carl von Ossietzky Universitaet Oldenburg
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> This module introduces anthropogenic heat emissions from external sources into PALM. It considers
!> heat from buildings, traffic, and other point sources. 
!> AH-profiles for BUILDINGS are typically obtained from a building energy model (BEM) 
!> that replicates PALM's buildings in terms of geometry and materials, but applies a different 
!> heating or cooling system to them.
!> POINT SORUCES are a more generic data fromat that can be used to represent any other source of 
!> anthropogenic heat: e.g. incineration plants, power plants, or industrial complexes.
!>
!> *potential future development*: The module could be extended to include heat from traffic, if:
!>  a) a traffic model is available that provides heat profiles, and
!>  b) PALM4U started working with unique street IDs, identifying each street in the model domain.
!
!--------------------------------------------------------------------------------------------------!
 MODULE anthropogenic_heat_mod

#if defined( __parallel )
    USE MPI
#endif
        
    USE basic_constants_and_equations_mod,                                                             &
       ONLY:  pi
    
    USE control_parameters,                                                                            &
       ONLY:  coupling_char,                                                                           &
              external_anthropogenic_heat,                                                             &
              initializing_actions,                                                                    &
              message_string,                                                                          &
              origin_date_time,                                                                        &
              time_since_reference_point                                                           

    USE cpulog,                                                                                        &
       ONLY:  cpu_log,                                                                                 &
              log_point_s
   
    USE grid_variables,                                                                                &
       ONLY:  ddx,                                                                                     &
              ddy
    
    USE indices,                                                                                       &
       ONLY:  nx,                                                                                      &
              ny

    USE kinds
    
    USE netcdf_data_input_mod,                                                                         &
       ONLY:  char_fill,                                                                               &
              check_existence,                                                                         &
              close_input_file,                                                                        &
              get_variable,                                                                            &
              get_attribute,                                                                           &
              get_dimension_length,                                                                    &
              init_model,                                                                              &
              input_pids_ah,                                                                           &
              inquire_num_variables,                                                                   &
              inquire_variable_names,                                                                  &
              input_file_ah,                                                                           &
              open_read_file
    
    USE pegrid

    USE surface_mod,                                                                                   &
       ONLY:  surf_usm,                                                                                &
              surf_lsm
    
    !
    !-- Define data types
    TYPE int_dimension
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: val      !< 1d integer-value array
       LOGICAL ::  from_file = .FALSE.                     !< flag indicating whether an input variable is available and read from file 
                                                           !< or default value is used
    END TYPE int_dimension
    
    TYPE real_2d_matrix
       REAL(wp) ::  fill                                  !< fill value
       REAL(wp), DIMENSION(:,:), ALLOCATABLE :: val       !< 2d real-value matrix
       LOGICAL ::  from_file = .FALSE.                    !< flag indicating whether an input variable is available and read from file 
                                                          !< or default value is used
    END TYPE real_2d_matrix

    TYPE point_coordinates
       REAL(wp) ::  x_fill                                !< fill value for x-coordinates
       REAL(wp) ::  y_fill                                !< fill value for y-coordinates
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  x          !< x-coordinates of point sources
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  y          !< y-coordinates of point sources
       LOGICAL ::  from_file = .FALSE.                    !< flag indicating whether an input variable is available and read from file 
                                                          !< or default value is used
    END TYPE point_coordinates

       
    IMPLICIT NONE
    
    PRIVATE

    TYPE(int_dimension) :: building_ids             !< ids of buildings with anthropogenic heat profiles
    TYPE(int_dimension) :: street_ids               !< ids of streets with anthropogenic heat profiles
    TYPE(int_dimension) :: point_ids                !< ids of point sources with anthropogenic heat profiles
    REAL(wp), DIMENSION(:), ALLOCATABLE :: ah_time  !< time steps


    TYPE(real_2d_matrix) :: building_ah       !< anthropogenic heat profiles for buildings
    TYPE(real_2d_matrix) :: street_ah         !< anthropogenic heat profiles for streets
    TYPE(real_2d_matrix) :: point_ah          !< anthropogenic heat profiles for point sources
    

    TYPE(point_coordinates) :: point_coords      !< exact coordinates of point sources
    
    
    SAVE
    
    INTERFACE ah_parin
       MODULE PROCEDURE ah_parin
    END INTERFACE ah_parin
    
    INTERFACE ah_profiles_netcdf_data_input
       MODULE PROCEDURE ah_profiles_netcdf_data_input
    END INTERFACE ah_profiles_netcdf_data_input
    
    INTERFACE ah_apply_to_surfaces
       MODULE PROCEDURE ah_apply_to_surfaces
    END INTERFACE ah_apply_to_surfaces
    
    
    PUBLIC                                                                                         &
       building_ah,                                                                                &
       street_ah,                                                                                  &
       point_ah
    
    PUBLIC                                                                                         &
       ah_parin,                                                                                   &
       ah_profiles_netcdf_data_input,                                                              &
       ah_apply_to_surfaces
    
    
    CONTAINS

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Parin for &anthropogenic_heat_par for external anthropogenic heat sources model
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_parin
    
       IMPLICIT NONE
  
       CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file
   
       INTEGER(iwp) ::  io_status   !< status after reading the namelist file
  
       LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                                !< although the respective module namelist appears in
                                                !< the namelist file
  
       NAMELIST /anthropogenic_heat_parameters/  switch_off_module                                  
  
       !
       !-- Move to the beginning of the namelist file and try to find and read the namelist.
       REWIND( 11 )
       READ( 11, anthropogenic_heat_parameters, IOSTAT=io_status )
  
       !
       !-- Action depending on the READ status
       IF ( io_status == 0 )  THEN
          !
          !--    anthropogenic_heat_parameters namelist was found and read correctly. Enable the
          !--    external anthropogenic heat module.
          IF ( .NOT. switch_off_module ) external_anthropogenic_heat = .TRUE.
  
       ELSEIF ( io_status > 0 )  THEN
          !
          !--    anthropogenic_heat_parameters namelist was found but contained errors. Print an error message
          !--    including the line that caused the problem.
          BACKSPACE( 11 )
          READ( 11 , '(A)' ) line
          CALL parin_fail_message( 'anthropogenic_heat_turbines', line )
  
       ENDIF
  
    END SUBROUTINE ah_parin


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Reads anthropogenic heat profiles emitted from building and ground surfaces from a NetCDF file.
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_profiles_netcdf_data_input
        
       IMPLICIT NONE

       INTEGER(iwp) ::  id_netcdf                                    !< NetCDF id of input file
       INTEGER(iwp) ::  num_vars                                     !< number of variables in input file
       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names   !< variable names in static input file   
       
       
       INTEGER(iwp) ::  n_buildings = 0           !< number of buildings (for array allocation)
       INTEGER(iwp) ::  n_streets = 0             !< number of streets (for array allocation)
       INTEGER(iwp) ::  n_points = 0              !< number of point sources (for array allocation)
       INTEGER(iwp) ::  n_timesteps = 0           !< number of time steps (for array allocation)
        
    !-- If no anthropogenic heat input file is available, skip this routine
       IF ( .NOT. input_pids_ah )  RETURN
    !
    !-- Measure CPU time
       CALL cpu_log( log_point_s(82), 'NetCDF input', 'start' )
    !
    !-- Skip the following if no expternal anthropogenic heat profiles are to be considered in the calculations.
       IF ( .NOT. external_anthropogenic_heat )  RETURN
        
    #if defined ( __netcdf )
    !
    !-- Open file in read-only mode
       CALL open_read_file( TRIM( input_file_ah ) // TRIM( coupling_char ) , id_netcdf )
    !
    !-- Inquire all variable names.
    !-- This will be used to check whether an optional input variable exists or not.
       CALL inquire_num_variables( id_netcdf, num_vars )
        
       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_netcdf, var_names )
    !
    !-- Read dimensions from file
       CALL get_dimension_length( id_netcdf, n_buildings, 'building_id' )
       CALL get_dimension_length( id_netcdf, n_streets, 'street_id' )
       CALL get_dimension_length( id_netcdf, n_points, 'point_id' )
       CALL get_dimension_length( id_netcdf, n_timesteps, 'time' )
    !
    !-- Read building ids from file 
       IF ( n_buildings > 0 .AND. check_existence( var_names, 'building_id' ) )  THEN
          building_ids%from_file = .TRUE.
        
          ALLOCATE( building_ids%val(1:n_buildings) ) 
        
          CALL get_variable( id_netcdf, 'building_id', building_ids%val )
       ELSE
          building_ids%from_file = .FALSE.
       ENDIF
    !
    !-- Read streed ids from file 
       IF ( n_streets > 0 .AND. check_existence( var_names, 'street_id' ) )  THEN
          street_ids%from_file = .TRUE.
      
          ALLOCATE( street_ids%val(1:n_streets) ) 
      
          CALL get_variable( id_netcdf, 'street_id', street_ids%val )
       ELSE
          street_ids%from_file = .FALSE.
       ENDIF
    !
    !-- Read point ids from file 
       IF ( n_points > 0 .AND. check_existence( var_names, 'point_id' ) )  THEN
         point_ids%from_file = .TRUE.
     
         ALLOCATE( point_ids%val(1:n_points) ) 
     
         CALL get_variable( id_netcdf, 'point_id', point_ids%val )
      ELSE
         point_ids%from_file = .FALSE.
      ENDIF
    !
    !-- Read timesteps from file
       IF ( n_timesteps > 0 .AND. check_existence( var_names, 'time' ) .AND. 
           ( building_ids%from_file .OR. street_ids%from_file .OR. point_ids%from_file ) )  THEN   
          ALLOCATE( ah_time(1:n_timesteps) ) 
         
          CALL get_variable( id_netcdf, 'time', ah_time )
       ENDIF

    !
    !-- Read anthrpogenic heat profiles from buildings from file
       IF ( building_ids%from_file .AND. check_existence( var_names, 'building_ah') ) THEN
          building_ah%from_file = .TRUE.
          CALL get_attribute( id_netcdf, char_fill, building_ah%fill, .FALSE., 'building_ah', .FALSE. )
        
          ALLOCATE( building_ah%val(1:n_buildings,1:n_timesteps) )
        
          CALL get_variable( id_netcdf, 'building_ah', building_ah%val )
          !-- TODO:(DONE) Output warning if the anthropogenic heat profile is not fully defined (use message.f90 module, use message_string from control_parameters under modules.f90)
          CALL ah_check_input_profiles('buildings', building_ah)
       ELSE
          building_ah%from_file = .FALSE.
       ENDIF
    !
    !-- Read anthrpogenic heat profiles from streets from file
       IF ( street_ids%from_file .AND. check_existence( var_names, 'street_ah') ) THEN
          street_ah%from_file = .TRUE.
          CALL get_attribute( id_netcdf, char_fill, street_ah%fill, .FALSE., 'street_ah', .FALSE. )
       
          ALLOCATE( street_ah%val(1:n_buildings,1:n_timesteps) )
       
          CALL get_variable( id_netcdf, 'street_ah', street_ah%val )
          !-- TODO:(DONE) Output warning if the anthropogenic heat profile is not fully defined (use message.f90 module, use message_string from control_parameters under modules.f90)
          CALL ah_check_input_profiles('streets', street_ah)
       ELSE
          street_ah%from_file = .FALSE.
       ENDIF
    !
    !-- Read anthropogenic heat profiles from points from file
       IF ( point_ids%from_file .AND. check_existence( var_names, 'point_ah') ) THEN
          point_ah%from_file = .TRUE.
          CALL get_attribute( id_netcdf, char_fill, point_ah%fill, .FALSE., 'point_ah', .FALSE. )
       
          ALLOCATE( point_ah%val(1:n_buildings,1:n_timesteps) )
       
          CALL get_variable( id_netcdf, 'point_ah', point_ah%val )
          !-- TODO:(DONE) Output warning if the anthropogenic heat profile is not fully defined (use message.f90 module)
          CALL ah_check_input_profiles('points', point_ah)
       ELSE
          point_ah%from_file = .FALSE.
       ENDIF
   
    !
    !-- Read coordinates of point sources from file
       IF ( point_ah%from_file .AND. check_existence( var_names, 'point_x') .AND. check_existence( var_names, 'point_y') ) THEN
         point_coords%from_file = .TRUE.
         CALL get_attribute( id_netcdf, char_fill, point_coords%x_fill, .FALSE., 'point_x', .FALSE. )
         CALL get_attribute( id_netcdf, char_fill, point_coords%y_fill, .FALSE., 'point_y', .FALSE. )
      
         ALLOCATE( point_coords%x(1:n_points) )
         ALLOCATE( point_coords%y(1:n_points) )
      
         CALL get_variable( id_netcdf, 'point_x', point_coords%x )
         CALL get_variable( id_netcdf, 'point_y', point_coords%y )

         !-- TODO:(DONE) Implement error handling if point coorinates are not fully defined (use message.f90 module)
         CALL ah_check_point_source_locations(point_coords)
      ELSE
         point_coords%from_file = .FALSE.
      ENDIF
    
    !
    !-- Finally, close input file
       CALL close_input_file( id_netcdf )
    #endif
    !
    !-- End of CPU measurement
       CALL cpu_log( log_point_s(82), 'NetCDF input', 'stop' )
    
    END SUBROUTINE ah_profiles_netcdf_data_input
    

    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Check anthropogenic heat profiles from input for buildings and point sources
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_check_input_profiles(source, ah_profile)

       IMPLICIT NONE

       CHARACTER(len=20), INTENT(IN) :: source            !< source of the anthropogenic heat profile
       TYPE(real_2d_matrix), INTENT(IN) :: ah_profile     !< anthropogenic heat profile

       !-- Check if the anthropogenic heat profile contains only non-negative values
       IF ( ANY( ( ah_profile%val < 0.0_wp ) .AND. ( ah_profile%val > ah_profile%fill )) ) THEN

         message_string = 'Anthropogenic heat profiles from ' // source // ' must be non-negative.'
         CALL message( 'netcdf_data_input_anthro_heat_profiles', 'AH0001', 1, 2, 0, 6, 0 )
       
       !-- Check if all values of the input anthropogenic heat profile have been defined
       ELSE IF ( ANY( ah_profile%val == ah_profile%fill ) ) THEN

         message_string = 'Some timesteps in the anthropogenic heat profile from ' // source // ' are not fully defined.'
         CALL message( 'netcdf_data_input_anthro_heat_profiles', 'AH0002', 0, 1, 0, 6, 0 ) 
         
       ENDIF 

    END SUBROUTINE ah_check_input_profiles


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Check anthropogenic heat profiles from input for point sources
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_check_point_source_locations(point_coordinates)

       IMPLICIT NONE

       TYPE(point_coords_t), INTENT(IN) :: point_coordinates     !< point source coordinates

       INTEGER(iwp) :: p                                         !< auxiliary index
       INTEGER(iwp) :: i_grid, j_grid                            !< domain grid indices of the point source

       !-- Check if any of the point source coordinates are not fully defined
       IF ( ANY( point_coordinates%x == point_coordinates%x_fill ) .OR. &
             ANY( point_coordinates%y == point_coordinates%y_fill ) ) THEN
   
          message_string = 'Some point sources lack fully defined coordinates.' // ACHAR(10) // & 
                            'Please make sure every point source is associated to a valid x and y coordinate.'
          CALL message( 'netcdf_data_input_anthro_heat_profiles', 'AH0003', 1, 2, 0, 6, 0 )

       ENDIF
      
       ! -- TODO: (Done) Check if the point source coordinates are within the model domain
       CALL metric_coords_to_grid_indices(point_coordinates%x, point_coordinates%y, i_grid, j_grid)
       IF ( i_grid < 0 .OR. i_grid > nx .OR. j_grid < 0 .OR. j_grid > ny ) THEN
          message_string = 'Some point sources are located outside the model domain.'
          CALL message( 'netcdf_data_input_anthro_heat_profiles', 'AH0004', 1, 2, 0, 6, 0 )
       ENDIF

    END SUBROUTINE ah_check_point_source_locations

    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This routine applies the imported anthropogenic heat profiles to the corresponding urban surfaces.
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_apply_to_surfaces

       IMPLICIT NONE

       INTEGER(iwp) :: t_step, b, m, p                            !< auxiliary indices
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE :: b_surf_indexes  !< array of building surface indexes
       INTEGER(iwp) :: b_surf_index                               !< building surface index
       INTEGER(iwp), ALLOCATABLE :: p_surf_index                  !< point source surface index

       ! -- Identify the time step to which the anthropogenic heat profiles will be applied
       t_exact = time_since_reference_point
       t_step = MINLOC( ABS( ah_time - MAX( t_exact, 0.0_wp ) ), DIM = 1 ) - 1

       ! TODO: Go over this routine again and make sure it is implemented correctly 
       !       (Can the time step really only be 1 ahead of the current time upon restart?)
       ! -- 
       ! -- Note, in case of restart runs, the time step for the boundary data may indicate a time in
       ! -- the future. This needs to be checked and corrected.
       IF ( TRIM( initializing_actions ) == 'read_restart_data'  .AND.                                &
          ah_time(t_step) > t_exact )  THEN
          t_step = t_step - 1
       ENDIF

       !
       ! Identify the relevant surfaces and apply the anthropogenic heat profiles to them
       ! -- for buildings
       DO  b = LBOUND(building_ids%val), UBOUND(building_ids%val)
          CALL ah_building_id_to_surfaces(building_ids%val(b), b_surf_indexes)
          ! -- distribute the anthropogenic heat profile of the building to the corresponding surfaces
          DO  m = LBOUND(b_surf_indexes), UBOUND(b_surf_indexes)
             b_surf_index = b_surf_indexes(m)
             surf_usm%waste_heat(surf_ind) = ( building_ah%val(b, t_step + 1) * ( time_since_reference_point - ah_time(t_step) ) +   &
                                               building_ah%val(b, t_step) * ( ah_time(t_step + 1) - time_since_reference_point ) )   &
                                             / ( ah_time(t_int + 1) - ah_time(t_step) ) &
                                             / SIZE(surf_indexes)
          ENDDO
       ENDDO

       ! -- for point sources
       DO p = LBOUND(point_ids%val), UBOUND(point_ids%val)
          CALL ah_point_id_to_surfaces(point_ids%val(p), p_surf_index)
          surf_usm%waste_heat(p_surf_index) = ( point_ah%val(p, t_step + 1) * ( time_since_reference_point - ah_time(t_step) ) +   &
                                                point_ah%val(p, t_step) * ( ah_time(t_step + 1) - time_since_reference_point ) )   &
                                              / ( ah_time(t_int + 1) - ah_time(t_step) )
       ENDDO


    END SUBROUTINE ah_apply_to_surfaces


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This routine fetches the roof surface tiles based on building ids
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_building_id_to_surfaces(building_id, surf_indexes)

       USE indoor_model_mod, ONLY: buildings

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) :: building_id                       !< building id
       INTEGER(iwp), DIMENSION(:), INTENT(OUT) :: surf_indexes       !< surface index

       INTEGER(iwp) :: b                                             !< auxiliary building index
       

       DO b = lbound(buildings), ubound(buildings)
          IF ( buildings(b)%id == building_id ) THEN
             ALLOCATE( surf_indexes(lbound(buildings(b)%m):ubound(buildings(b)%m)) )
             surf_indexes = buildings(b)%m
             EXIT
          ENDIF
       END DO

    END SUBROUTINE ah_building_id_to_surfaces


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> This routine fetches the ground surface tiles based on point-source ids
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_point_id_to_surface(point_id, surf_index)

       IMPLICIT NONE

       INTEGER(iwp), INTENT(IN) :: point_id     !< point id
       INTEGER(iwp), INTENT(OUT) :: surf_index  !< surface index

       INTEGER(iwp) :: p                        !< auxiliary point source index
       REAL(wp) :: x_coord_abs, y_coord_abs     !< absolute metric coordinates of the point source
       INTEGER(iwp) :: i, j, m, is, js          !< auxiliary surface indices
       LOGICAL :: found = .FALSE.               !< flag to indicate if a match was found for the point source

       IF ( point_coords%from_file ) THEN
          !-- Retrieve metric coordinates of point sources from file
          x_coord_abs = point_coords%x(point_id)
          y_coord_abs = point_coords%y(point_id)

          !-- Transform metric coordinates to grid indices
          CALL metric_coords_to_grid_indices(x_coord_abs, y_coord_abs, is, js)

          ! TODO:(DONE) Finish writing comments
          ! -- Find the surface index of the grid cell in which the point source is located.
          ! -- First, search the land surfaces (lsm) for a match with the coordinates of the point source.
          DO m = 1, surf_lsm%ns
             i = surf_lsm%i(m)
             j = surf_lsm%j(m)
             IF ( i == is .AND. j == js ) THEN
                surf_index = m
                found = .TRUE.
                EXIT
             ENDIF
          END DO

          ! -- If no match was found in the land surfaces, search the urban surfaces (usm).
          IF (.NOT. found) THEN
             DO m = 1, surf_usm%ns   ! TODO: Limit this looping to horizontal surfaces
                i = surf_usm%i(m)
                j = surf_usm%j(m)
                IF ( i == is .AND. j == js ) THEN
                   surf_index = m
                   found = .TRUE.
                   EXIT
                ENDIF
             END DO
          ENDIF
         
       ENDIF

    END SUBROUTINE ah_point_id_to_surfaces


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Transform metric point coordinates (x,y) to their equivalent grid point indices in the 
    !> model domain.
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE metric_coords_to_grid_indices(x_coord_metric, y_coord_metric, i, j)

       IMPLICIT NONE

       REAL(wp), INTENT(IN) :: x_coord_metric  !< x-coordinate of the point source
       REAL(wp), INTENT(IN) :: y_coord_metric  !< y-coordinate of the point source
       INTEGER(iwp), INTENT(OUT) :: i          !< x-index of the grid cell
       INTEGER(iwp), INTENT(OUT) :: j          !< y-index of the grid cell

       REAL(wp) :: x_coord_rel, y_coord_rel    !< metric point coordinates relative to the model origin
       REAL(wp) :: x_coord, y_coord            !< metric point coordinates after translation and rotation 
                                               !< to match the model domain grid

       !
       ! -- offset metric point coordinates to reflect their relative position to the model origin
       x_coord_rel = x_coord_metric - init_model%origin_x
       y_coord_rel = y_coord_metric - init_model%origin_y
       
       !
       ! -- Rotate the metric point coordinates to align with the model grid
       x_coord = COS( init_model%rotation_angle * pi / 180.0_wp ) * x_coord_rel(t)            &
               - SIN( init_model%rotation_angle * pi / 180.0_wp ) * y_coord_rel(t)
       y_coord = SIN( init_model%rotation_angle * pi / 180.0_wp ) * x_coord_rel(t)            &
               + COS( init_model%rotation_angle * pi / 180.0_wp ) * y_coord_rel(t)

       !
       ! -- Then, compute the indices of the grid cell in which the point source is located.
       i = INT( x_coord * ddx, KIND = iwp )
       j = INT( y_coord * ddy, KIND = iwp )

    END SUBROUTINE metric_coords_to_grid_indices
   
   
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Check namelist parameter TODO: Check if this method is needed in AH mod
    !--------------------------------------------------------------------------------------------------!
    SUBROUTINE ah_check_parameters
   
    END SUBROUTINE ah_check_parameters
    
    
    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Initialization of the anthropogenic heat model TODO: Check if this method is needed in AH mod
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE ah_init
    
     END SUBROUTINE ah_init


    !--------------------------------------------------------------------------------------------------!
    ! Description:
    ! ------------
    !> Print header for the anthropogenic heat model TODO: Check if this method is needed in AH mod
    !--------------------------------------------------------------------------------------------------!
     SUBROUTINE ah_header
    
     END SUBROUTINE ah_header

    
 END MODULE anthropogenic_heat_mod
    