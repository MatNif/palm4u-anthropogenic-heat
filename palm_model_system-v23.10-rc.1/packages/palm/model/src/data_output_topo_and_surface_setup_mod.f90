!> @file data_output_topo_and_surface_setup_mod.f90
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
! Authors:
! --------
! @author Matthias Suehring
!
!
! Description:
! ------------
!> Output of filtered and modified topography, final surface classification and surface types.
!--------------------------------------------------------------------------------------------------!
 MODULE data_output_topo_and_surface_setup_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  x,                                                                                  &
               y,                                                                                  &
               zu,                                                                                 &
               zw

    USE control_parameters,                                                                        &
        ONLY:  coupling_char,                                                                      &
               message_string,                                                                     &
               topography

    USE data_output_module

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               nxl,                                                                                &
               nxr,                                                                                &
               ny,                                                                                 &
               nyn,                                                                                &
               nys,                                                                                &
               nzb,                                                                                &
               nzb_max,                                                                            &
               nzt,                                                                                &
               topo_flags

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  buildings_f,                                                                        &
               building_type_f,                                                                    &
               building_id_f,                                                                      &
               init_model,                                                                         &
               input_pids_static,                                                                  &
               pavement_type_f,                                                                    &
               terrain_height_f,                                                                   &
               vegetation_type_f,                                                                  &
               water_type_f

    USE pegrid

    USE surface_mod,                                                                               &
        ONLY:  surf_lsm,                                                                           &
               surf_usm

    IMPLICIT NONE

    CHARACTER(LEN=14)  ::  filename = 'DATA_TOPO_SURF' !< file name
    CHARACTER(LEN=10)  ::  char_fill = '_FillValue'    !< attribute name for fill value
    CHARACTER(LEN=30)  ::  nc_filename                 !< output file name (file name + coupling_char)
    CHARACTER(LEN=100) ::  variable_name               !< name of output variable
    CHARACTER(LEN=800) ::  dom_error_message           !< error message returned by the data-output module

    INTEGER(iwp) ::  k_rel_max_buildings !< max. building height in index range
    INTEGER(iwp) ::  k_ttop_ij           !< terrain-top index at location i,j
    INTEGER(iwp) ::  return_value        !< returned status value of called data-output function

    SAVE

    INTERFACE do_topo_and_surface_actions
       MODULE PROCEDURE do_topo_and_surface_actions
    END INTERFACE do_topo_and_surface_actions

    INTERFACE do_topo_and_surface_init_output
       MODULE PROCEDURE do_topo_and_surface_init_output
    END INTERFACE do_topo_and_surface_init_output
!
!-- Public subroutines
    PUBLIC do_topo_and_surface_actions,                                                            &
           do_topo_and_surface_init_output

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write topography and surface setup to netCDF file.
!> @NOTE So far, topography output is only necessary at the beginning of a simulation. Later on, if
!> topography may include moving objects, this may change and topography may become time-dependent.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE do_topo_and_surface_actions( location )

       CHARACTER(LEN=*) ::  location !< call location string

       INTEGER(iwp) ::  i  !< running index in x-direction
       INTEGER(iwp) ::  j  !< running index in y-direction
       INTEGER(iwp) ::  k  !< running index in z-direction
       INTEGER(iwp) ::  m  !< running index surface elemets

       INTEGER(ibp), DIMENSION(:,:,:), POINTER              ::  output_values_3d_int8_pointer   !< pointer to output array
       INTEGER(ibp), DIMENSION(:,:,:), ALLOCATABLE, TARGET  ::  output_values_3d_int8_target    !< output array

       INTEGER(iwp), DIMENSION(:,:),   POINTER              ::  output_values_2d_int32_pointer  !< pointer to output array
       INTEGER(iwp), DIMENSION(:,:),   ALLOCATABLE, TARGET  ::  output_values_2d_int32_target   !< output array
       INTEGER(iwp), DIMENSION(:,:,:), POINTER              ::  output_values_3d_int32_pointer  !< pointer to output array
       INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE, TARGET  ::  output_values_3d_int32_target   !< output array

       LOGICAL ::  output_required  !< flag indicating if netCDF output of topograhy and surface classification is required

       REAL(sp), DIMENSION(:,:), POINTER             ::  output_values_2d_real32_pointer  !< pointer to output array
       REAL(sp), DIMENSION(:,:), ALLOCATABLE, TARGET ::  output_values_2d_real32_target   !< output array

       SELECT CASE ( location )
!
!--       Initial output of topography and surface-classification information.
          CASE ( 'before_integration' )
!
!--          Check if output file has been created. If it has not been created (no topography
!--          or surface classification), skip the following code.
             INQUIRE( FILE = TRIM( filename // TRIM( coupling_char ) ), EXIST = output_required )

             IF ( .NOT. output_required )  RETURN
!
!--          Output of terrain height information.
             IF ( terrain_height_f%from_file )  THEN
                variable_name = 'zt'
                ALLOCATE( output_values_2d_real32_target(nys:nyn,nxl:nxr) )
                output_values_2d_real32_target = terrain_height_f%fill
!
!--             Obtain filtered and modified terrain height from flag array.
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         output_values_2d_real32_target(j,i) = MERGE(                              &
                                                             REAL( zw(k), KIND = sp ),             &
                                                             output_values_2d_real32_target(j,i),  &
                                                             BTEST( topo_flags(k,j,i), 5 ) )
                      ENDDO
                   ENDDO
                ENDDO

                output_values_2d_real32_pointer => output_values_2d_real32_target
                return_value = dom_write_var( nc_filename, variable_name,                          &
                                              values_real32_2d = output_values_2d_real32_pointer,  &
                                              bounds_start = (/nxl, nys/),                         &
                                              bounds_end   = (/nxr, nyn/) )
                DEALLOCATE( output_values_2d_real32_target )

             ENDIF
!
!--          Output of building information.
             IF ( buildings_f%from_file )  THEN
!
!--             2D buildings
                IF ( buildings_f%lod == 1 )  THEN
                   variable_name = 'buildings_2d'
                   ALLOCATE( output_values_2d_real32_target(nys:nyn,nxl:nxr) )
                   output_values_2d_real32_target = buildings_f%fill1
!
!--                Obtain filtered and modified 2D building height from flag array. The building
!--                height is output relative to the terrain height here. Therefore, first determine
!--                top index of uppermost terrain grid point.
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         k_ttop_ij = 0
!
!--                      Determine terrain-top index. Note, in case of ASCII topography skip this
!--                      part as the entire topography is assumed building (input_pids_static = .F.,
!--                      i.e. no static input file).
                         IF ( input_pids_static )  THEN
                            DO  k = nzb, nzt+1
                               k_ttop_ij = MERGE( k, k_ttop_ij, BTEST( topo_flags(k,j,i), 5 ) )
                            ENDDO
                         ENDIF
                         DO  k = k_ttop_ij+1, nzt+1
                            output_values_2d_real32_target(j,i) = MERGE(                           &
                                                          REAL( zw(k) - zw(k_ttop_ij), KIND = sp ),&
                                                          output_values_2d_real32_target(j,i),     &
                                                          BTEST( topo_flags(k,j,i), 6 ) )
                         ENDDO
                      ENDDO
                   ENDDO
                   output_values_2d_real32_pointer => output_values_2d_real32_target
                   return_value = dom_write_var( nc_filename, variable_name,                       &
                                                 values_real32_2d=output_values_2d_real32_pointer, &
                                                 bounds_start = (/nxl, nys/),                      &
                                                 bounds_end   = (/nxr, nyn/) )
                   DEALLOCATE( output_values_2d_real32_target )
                ENDIF
!
!--             3D buildings
                IF ( buildings_f%lod == 2 )  THEN
                   variable_name = 'buildings_3d'
                   ALLOCATE( output_values_3d_int8_target(nzb:k_rel_max_buildings,nys:nyn,nxl:nxr) )
                   output_values_3d_int8_target = buildings_f%fill2
!
!--                Obtain filtered and modified 2D building height from flag array. The building
!--                height is output relative to the terrain height here. Therefore, first determine
!--                top index of uppermost terrain grid point.
                   DO  i = nxl, nxr
                      DO  j = nys, nyn
                         k_ttop_ij = 0
                         DO  k = nzb, nzt+1
                            k_ttop_ij = MERGE( k, k_ttop_ij, BTEST( topo_flags(k,j,i), 5 ) )
                         ENDDO
                         DO  k = k_ttop_ij+1, k_ttop_ij + k_rel_max_buildings
                            output_values_3d_int8_target(k-k_ttop_ij,j,i) = INT(                   &
                                                    MERGE( 1, 0,  BTEST( topo_flags(k,j,i), 6 ) ), &
                                                                                 KIND = ibp )
                         ENDDO
                         output_values_3d_int8_target(nzb,j,i) =                                   &
                                                        output_values_3d_int8_target(nzb+1,j,i)
                      ENDDO
                   ENDDO
                   output_values_3d_int8_pointer => output_values_3d_int8_target
                   return_value = dom_write_var( nc_filename, variable_name,                       &
                                                 values_int8_3d = output_values_3d_int8_pointer,   &
                                                 bounds_start = (/nxl, nys, nzb/),                 &
                                                 bounds_end   = (/nxr, nyn, k_rel_max_buildings/) )
                   DEALLOCATE( output_values_3d_int8_target )
                ENDIF
             ENDIF
!
!--          Output of merged topography information (terrain + buildings).
             IF ( topography /= 'flat' )  THEN
                variable_name = 'topo_all'
                ALLOCATE( output_values_3d_int32_target(nzb:nzb_max,nys:nyn,nxl:nxr) )
!
!--             Initialize with a fill value. Take the one from building-id.
                output_values_3d_int32_target = building_id_f%fill
!
!--             Set total topography array. 0 - no obstacle, 1 - building, 2 - terrain. In contrast
!--             to output of buildings, output of entire topography is relative to absolute
!--             coordinates.
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzb_max
                         IF ( BTEST( topo_flags(k,j,i), 0 ) )  THEN
                            output_values_3d_int32_target(k,j,i) = 0
                         ELSEIF ( BTEST( topo_flags(k,j,i), 6 ) )  THEN
                            output_values_3d_int32_target(k,j,i) = 1
                         ELSEIF ( BTEST( topo_flags(k,j,i), 5 ) )  THEN
                            output_values_3d_int32_target(k,j,i) = 2
                         ELSE
                            output_values_3d_int32_target(k,j,i) = 3
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
                output_values_3d_int32_pointer => output_values_3d_int32_target
                return_value = dom_write_var( nc_filename, variable_name,                          &
                                                 values_int32_3d = output_values_3d_int32_pointer, &
                                                 bounds_start = (/nxl, nys, nzb/),                 &
                                                 bounds_end   = (/nxr, nyn, nzb_max/) )
                DEALLOCATE( output_values_3d_int32_target )
             ENDIF
!
!--          Output of building ids.
             IF ( building_id_f%from_file )  THEN
                variable_name = 'building_id'
                ALLOCATE( output_values_2d_int32_target(nys:nyn,nxl:nxr) )
                output_values_2d_int32_target = building_id_f%fill

                DO  m = 1, surf_usm%ns
                   IF ( surf_usm%upward(m) )  THEN
                      i = surf_usm%i(m)
                      j = surf_usm%j(m)
                      output_values_2d_int32_target(j,i) = building_id_f%var(j,i)
                   ENDIF
                ENDDO

                output_values_2d_int32_pointer => output_values_2d_int32_target
                return_value = dom_write_var( nc_filename, variable_name,                          &
                                              values_int32_2d = output_values_2d_int32_pointer,    &
                                              bounds_start = (/nxl, nys/),                         &
                                              bounds_end   = (/nxr, nyn/) )
                DEALLOCATE( output_values_2d_int32_target )
             ENDIF
!
!--          Output of building types.
             IF ( building_type_f%from_file )  THEN
                variable_name = 'building_type'
                ALLOCATE( output_values_2d_int32_target(nys:nyn,nxl:nxr) )
                output_values_2d_int32_target = building_type_f%fill

                DO  m = 1, surf_usm%ns
                   IF ( surf_usm%upward(m) )  THEN
                      i = surf_usm%i(m)
                      j = surf_usm%j(m)
                      output_values_2d_int32_target(j,i) = surf_usm%building_type(m)
                   ENDIF
                ENDDO

                output_values_2d_int32_pointer => output_values_2d_int32_target
                return_value = dom_write_var( nc_filename, variable_name,                          &
                                              values_int32_2d = output_values_2d_int32_pointer,    &
                                              bounds_start = (/nxl, nys/),                         &
                                              bounds_end   = (/nxr, nyn/) )
                DEALLOCATE( output_values_2d_int32_target )
             ENDIF
!
!--          Output of vegetation types.
             IF ( vegetation_type_f%from_file )  THEN
                variable_name = 'vegetation_type'
                ALLOCATE( output_values_2d_int32_target(nys:nyn,nxl:nxr) )
                output_values_2d_int32_target = vegetation_type_f%fill

                DO  m = 1, surf_lsm%ns
                   IF ( surf_lsm%upward(m)  .AND.  surf_lsm%vegetation_surface(m) )  THEN
                      i = surf_lsm%i(m)
                      j = surf_lsm%j(m)
                      output_values_2d_int32_target(j,i) = surf_lsm%vegetation_type(m)
                   ENDIF
                ENDDO

                output_values_2d_int32_pointer => output_values_2d_int32_target
                return_value = dom_write_var( nc_filename, variable_name,                          &
                                              values_int32_2d = output_values_2d_int32_pointer,    &
                                              bounds_start = (/nxl, nys/),                         &
                                              bounds_end   = (/nxr, nyn/) )
                DEALLOCATE( output_values_2d_int32_target )
             ENDIF
!
!--          Output of pavement types.
             IF ( pavement_type_f%from_file )  THEN
                variable_name = 'pavement_type'
                ALLOCATE( output_values_2d_int32_target(nys:nyn,nxl:nxr) )
                output_values_2d_int32_target = pavement_type_f%fill

                DO  m = 1, surf_lsm%ns
                   IF ( surf_lsm%upward(m)  .AND.  surf_lsm%pavement_surface(m) )  THEN
                      i = surf_lsm%i(m)
                      j = surf_lsm%j(m)
                      output_values_2d_int32_target(j,i) = surf_lsm%pavement_type(m)
                   ENDIF
                ENDDO

                output_values_2d_int32_pointer => output_values_2d_int32_target
                return_value = dom_write_var( nc_filename, variable_name,                          &
                                              values_int32_2d = output_values_2d_int32_pointer,    &
                                              bounds_start = (/nxl, nys/),                         &
                                              bounds_end   = (/nxr, nyn/) )
                DEALLOCATE( output_values_2d_int32_target )
             ENDIF
!
!--          Output of water types.
             IF ( water_type_f%from_file )  THEN
                variable_name = 'water_type'
                ALLOCATE( output_values_2d_int32_target(nys:nyn,nxl:nxr) )
                output_values_2d_int32_target = water_type_f%fill

                DO  m = 1, surf_lsm%ns
                   IF ( surf_lsm%upward(m)  .AND.  surf_lsm%water_surface(m) )  THEN
                      i = surf_lsm%i(m)
                      j = surf_lsm%j(m)
                      output_values_2d_int32_target(j,i) = surf_lsm%water_type(m)
                   ENDIF
                ENDDO

                output_values_2d_int32_pointer => output_values_2d_int32_target
                return_value = dom_write_var( nc_filename, variable_name,                          &
                                              values_int32_2d = output_values_2d_int32_pointer,    &
                                              bounds_start = (/nxl, nys/),                         &
                                              bounds_end   = (/nxr, nyn/) )
                DEALLOCATE( output_values_2d_int32_target )
             ENDIF


          CASE DEFAULT
             CONTINUE

       END SELECT

!
!--    Check if DOM reported any error.
!--    Attention: In case of DATA_TOPO_SURF DOM occasionally reports the error "Parallel operation
!--    on file opened for non-parallel access", although output seems to be fine. For this reason
!--    the status is declared as WARNING instead of ERROR.
       dom_error_message = dom_get_error_message()
       IF ( TRIM( dom_error_message ) /= '' )  THEN
          message_string = 'issue while writing output: "' // dom_error_message // '"'
          CALL message( 'do_topo_and_surface_actions', 'PAC0191', 0, 1, 0, 6, 0 )
       ENDIF

    END SUBROUTINE do_topo_and_surface_actions

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize netCDF output file, set global attributes, define dimensions and variables.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE do_topo_and_surface_init_output

       INTEGER(iwp) ::  i  !< running index in x-direction
       INTEGER(iwp) ::  j  !< running index in y-direction
       INTEGER(iwp) ::  k  !< running index in z-direction
!
!--    Check if output is required. Output is only required when topography is defined or
!--    any surfaces are classified. The latter is only the case if a static input file is
!--    available. If there is no need to write this static information, skip the initialization
!--    of the topography and surface output file.
       IF ( topography == 'flat'  .AND.  .NOT. input_pids_static )  RETURN
!
!--    Topography and surface output is only realized for parallel netCDF. So, the output file is
!--    only defined in this case. If parallel netCDF is not available, no file is defined and the
!--    further output definition must be skipped.
#if defined( __netcdf4_parallel )
!
!--    Define output file. Note, DOM automatically adds the coupling_char to the filename!
       nc_filename = TRIM( filename )
       return_value = dom_def_file( nc_filename, 'netcdf4-parallel' )
#else
       message_string = 'Topography and surface-setup output requires parallel netCDF. ' //        &
                        'No output file will be created.'
       CALL message( 'do_topo_and_surface_setup', 'PAC0192', 0, 1, 0, 6, 0 )

       RETURN
#endif
!
!--    Define global attributes.
       return_value = dom_def_att( nc_filename, attribute_name = 'title',                          &
                                   value = 'Topography and surface classification')
       return_value = dom_def_att( nc_filename, attribute_name = 'source', value = 'PALM-4U')
       return_value = dom_def_att( nc_filename, attribute_name = 'origin_x',                       &
                                   value = init_model%origin_x )
       return_value = dom_def_att( nc_filename, attribute_name = 'origin_y',                       &
                                   value = init_model%origin_y )
       return_value = dom_def_att( nc_filename, attribute_name = 'origin_lon',                     &
                                   value = init_model%longitude )
       return_value = dom_def_att( nc_filename, attribute_name = 'origin_lat',                     &
                                   value = init_model%latitude )
       return_value = dom_def_att( nc_filename, attribute_name = 'origin_z',                       &
                                   value = init_model%origin_z )
       return_value = dom_def_att( nc_filename, attribute_name = 'rotation_angle',                 &
                                   value = init_model%rotation_angle )
!
!--    Output of z-shift to bring the global minimum in terrain height to zero. Note, this
!--    will be only enabled when topography initialization is modularized.
!        return_value = dom_def_att( nc_filename, attribute_name = 'z_shift', value = oro_min )

!
!--    Define spatial dimensions. Note, z-dimension is only defined up to topography top
!--    (indicated by nzb_max).
!--    zb is dimensioned according to the maximum building height.
!--    Before dimensions will be defined, determine maximum building height in order to properly
!--    set dimension bounds for building output.
       k_rel_max_buildings = 0
       k_ttop_ij = 0
       DO  i = nxl, nxr
          DO  j = nys, nyn
             k_ttop_ij = 0
             DO  k = nzb, nzt+1
                k_ttop_ij = MERGE( k, k_ttop_ij, BTEST( topo_flags(k,j,i), 5 ) )
             ENDDO
             DO  k = k_ttop_ij + 1, nzt+1
                k_rel_max_buildings = MERGE( k-k_ttop_ij, k_rel_max_buildings,                     &
                                             BTEST( topo_flags(k,j,i), 6 ) )
             ENDDO
          ENDDO
       ENDDO
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, k_rel_max_buildings, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr)
#endif
!
!--    Define grid coordinates and its attributes.
!--    x-coordinate
       return_value = dom_def_dim( nc_filename, dimension_name = 'x', output_type = 'real32',      &
                                   bounds = (/0, nx/), values_realwp = x(0:nx) )
       return_value = dom_def_att( nc_filename, variable_name = 'x', attribute_name = 'long_name', &
                                   value = 'distance to origin in x-direction' )
       return_value = dom_def_att( nc_filename, variable_name = 'x', attribute_name = 'axis',      &
                                   value = 'X' )
       return_value = dom_def_att( nc_filename, variable_name = 'x', attribute_name = 'units',     &
                                   value = 'm' )
!
!--    y-coordinate
       return_value = dom_def_dim( nc_filename, dimension_name = 'y', output_type = 'real32',      &
                                   bounds = (/0, ny/), values_realwp = y(0:ny) )
       return_value = dom_def_att( nc_filename, variable_name = 'y', attribute_name = 'long_name', &
                                   value = 'distance to origin in y-direction' )
       return_value = dom_def_att( nc_filename, variable_name = 'y', attribute_name = 'axis',      &
                                   value = 'Y' )
       return_value = dom_def_att( nc_filename, variable_name = 'y', attribute_name = 'units',     &
                                   value = 'm' )
!
!--    z-coordinate up to topography top
       return_value = dom_def_dim( nc_filename, dimension_name = 'z', output_type = 'real32',      &
                                   bounds = (/0, nzb_max/), values_realwp = zu(0:nzb_max) )
       return_value = dom_def_att( nc_filename, variable_name = 'z', attribute_name = 'long_name', &
                                   value = 'height above origin' )
       return_value = dom_def_att( nc_filename, variable_name = 'z', attribute_name = 'axis',      &
                                   value = 'Z' )
       return_value = dom_def_att( nc_filename, variable_name = 'z', attribute_name = 'positive',  &
                                   value = 'up' )
       return_value = dom_def_att( nc_filename, variable_name = 'z', attribute_name = 'units',     &
                                   value = 'm' )
!
!--    z-coordinate for buildings only
       return_value = dom_def_dim( nc_filename, dimension_name = 'zb', output_type = 'real32',     &
                                   bounds = (/0, k_rel_max_buildings/),                            &
                                   values_realwp = zu(0:k_rel_max_buildings) )
       return_value = dom_def_att( nc_filename, variable_name = 'zb', attribute_name = 'long_name',&
                                   value = 'z-coordinate for buildings' )
       return_value = dom_def_att( nc_filename, variable_name = 'zb', attribute_name = 'units',    &
                                   value = 'm' )
!
!--    Define output for terrain height.
       IF ( terrain_height_f%from_file )  THEN
          variable_name = 'zt'
          return_value = dom_def_var( nc_filename, variable_name = variable_name,                  &
                                      dimension_names = (/ 'x', 'y'/),                             &
                                      output_type = 'real32' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = char_fill,                                  &
                                      value = REAL( terrain_height_f%fill, KIND = sp ) )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'long_name',                                &
                                      value = 'ground level altitude' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'standard_name',                            &
                                      value = 'ground_level_altitude' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'units',                                    &
                                      value = 'm' )
       ENDIF
!
!--    Define output for buildings. Depending on the given LOD, either 2D building heights or
!--    3D building information will be output. Assure that also deprecated ASCII topography
!--    as well as generic topography works.
       IF ( buildings_f%from_file )  THEN
          IF ( buildings_f%lod == 1 )  THEN
             variable_name = 'buildings_2d'
             return_value = dom_def_var( nc_filename, variable_name = variable_name,               &
                                         dimension_names = (/ 'x', 'y'/),                          &
                                         output_type = 'real32' )
             return_value = dom_def_att( nc_filename, variable_name = variable_name,               &
                                         attribute_name = char_fill,                               &
                                         value = REAL( buildings_f%fill1, KIND = sp ) )
             return_value = dom_def_att( nc_filename, variable_name = variable_name,               &
                                         attribute_name = 'long_name',                             &
                                         value = 'building height' )
             return_value = dom_def_att( nc_filename, variable_name = variable_name,               &
                                         attribute_name = 'units',                                 &
                                         value = 'm' )
          ENDIF
!
!--       Note, output of 3D building information is terrain-following.
          IF ( buildings_f%lod == 2 )  THEN
             variable_name = 'buildings_3d'
             return_value = dom_def_var( nc_filename, variable_name = variable_name,               &
                                         dimension_names = (/ 'x ', 'y ', 'zb'/),                  &
                                         output_type = 'int8' )
             return_value = dom_def_att( nc_filename, variable_name = variable_name,               &
                                         attribute_name = char_fill,                               &
                                         value = buildings_f%fill2 )
             return_value = dom_def_att( nc_filename, variable_name = variable_name,               &
                                         attribute_name = 'long_name',                             &
                                         value = 'building flag' )
             return_value = dom_def_att( nc_filename, variable_name = variable_name,               &
                                         attribute_name = 'units',                                 &
                                         value = '1' )
             return_value = dom_def_att( nc_filename, variable_name = variable_name,               &
                                         attribute_name = 'values',                                &
                                         value = '0: non topography, 1: building' )
          ENDIF
       ENDIF
!
!--    Define output for total topography (terrain + buildings).
       IF ( topography /= 'flat' )  THEN
          variable_name = 'topo_all'
          return_value = dom_def_var( nc_filename, variable_name = variable_name,                  &
                                      dimension_names = (/ 'x', 'y', 'z'/),                        &
                                      output_type = 'int32' )
!
!--       Set some _FillValue. Simply take it from building_id_f structure.
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = char_fill,                                  &
                                      value = building_id_f%fill )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'long_name',                                &
                                      value = 'all topography' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'units',                                    &
                                      value = '1' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'values',                                   &
                                      value = '0: non topography, 1: building, 2: terrain, ' //    &
                                              '3: non-classified topography' )
       ENDIF
!
!--    Define output for building id.
       IF ( building_id_f%from_file )  THEN
          variable_name = 'building_id'
          return_value = dom_def_var( nc_filename, variable_name = variable_name,                  &
                                      dimension_names = (/ 'x', 'y'/),                             &
                                      output_type = 'int32' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = char_fill,                                  &
                                      value = building_id_f%fill )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'long_name',                                &
                                      value = 'building id number' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'units',                                    &
                                      value = '1' )
       ENDIF
!
!--    Define output for building type.
       IF ( building_type_f%from_file )  THEN
          variable_name = 'building_type'
          return_value = dom_def_var( nc_filename, variable_name = variable_name,                  &
                                      dimension_names = (/ 'x', 'y'/),                             &
                                      output_type = 'int32' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = char_fill,                                  &
                                      value = INT( building_type_f%fill, KIND = iwp ) )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'long_name',                                &
                                      value = 'building type classification' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'units',                                    &
                                      value = '1' )
       ENDIF
!
!--    Define output for vegetation type.
       IF ( vegetation_type_f%from_file )  THEN
          variable_name = 'vegetation_type'
          return_value = dom_def_var( nc_filename, variable_name = variable_name,                  &
                                      dimension_names = (/ 'x', 'y'/),                             &
                                      output_type = 'int32' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = char_fill,                                  &
                                      value = INT( vegetation_type_f%fill, KIND = iwp ) )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'long_name',                                &
                                      value = 'vegetation type classification' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'units',                                    &
                                      value = '1' )
       ENDIF
!
!--    Define output for pavement type.
       IF ( pavement_type_f%from_file )  THEN
          variable_name = 'pavement_type'
          return_value = dom_def_var( nc_filename, variable_name = variable_name,                  &
                                      dimension_names = (/ 'x', 'y'/),                             &
                                      output_type = 'int32' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = char_fill,                                  &
                                      value = INT( pavement_type_f%fill, KIND = iwp ) )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'long_name',                                &
                                      value = 'pavement type classification' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'units',                                    &
                                      value = '1' )
       ENDIF
!
!--    Define output for water type.
       IF ( water_type_f%from_file )  THEN
          variable_name = 'water_type'
          return_value = dom_def_var( nc_filename, variable_name = variable_name,                  &
                                      dimension_names = (/ 'x', 'y'/),                             &
                                      output_type = 'int32' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = char_fill,                                  &
                                      value = INT( water_type_f%fill, KIND = iwp ) )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'long_name',                                &
                                      value = 'water type classification' )
          return_value = dom_def_att( nc_filename, variable_name = variable_name,                  &
                                      attribute_name = 'units',                                    &
                                      value = '1' )
       ENDIF

!
!--    Check if DOM reported any error
!--    Attention: In case of DATA_TOPO_SURF DOM occasionally reports the error "Parallel operation
!--    on file opened for non-parallel access", although output seems to be fine. For this reason
!--    the status is declared as WARNING instead of ERROR.
       dom_error_message = dom_get_error_message()
       IF ( TRIM( dom_error_message ) /= '' )  THEN
          message_string = 'issue while defining output: "' // dom_error_message // '"'
          CALL message( 'do_topo_and_surface_init_output', 'PAC0193', 0, 1, 0, 6, 0 )
       ENDIF

    END SUBROUTINE do_topo_and_surface_init_output

 END MODULE data_output_topo_and_surface_setup_mod
