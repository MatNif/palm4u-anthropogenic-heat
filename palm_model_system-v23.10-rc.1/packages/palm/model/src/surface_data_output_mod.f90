!> @file surface_data_output_mod.f90
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
!
! Authors:
! --------
! @author Klaus Ketelsen, Matthias Suehring, Tobias Gronemeier
!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Generate output for surface data.
!>
!> @todo Create namelist file for post-processing tool.
!--------------------------------------------------------------------------------------------------!
MODULE surface_data_output_mod

#if defined( __parallel )
    USE MPI
#endif

   USE kinds

   USE arrays_3d,                                                                                  &
       ONLY:  heatflux_output_conversion,                                                          &
              momentumflux_output_conversion,                                                      &
              scalarflux_output_conversion,                                                        &
              waterflux_output_conversion,                                                         &
              zu,                                                                                  &
              zw

   USE control_parameters,                                                                         &
       ONLY:  coupling_char,                                                                       &
              cyclic_fill_initialization,                                                          &
              data_output_during_spinup,                                                           &
              end_time,                                                                            &
              flux_output_mode,                                                                    &
              message_string,                                                                      &
              origin_date_time,                                                                    &
              restart_data_format_output,                                                          &
              run_description_header,                                                              &
              simulated_time_at_begin,                                                             &
              spinup_time,                                                                         &
              surface_output

   USE grid_variables,                                                                             &
       ONLY: dx,                                                                                   &
             dy

   USE indices,                                                                                    &
       ONLY: nxl,                                                                                  &
             nxr,                                                                                  &
             nys,                                                                                  &
             nyn,                                                                                  &
             nzb,                                                                                  &
             nzt

   USE module_interface,                                                                           &
       ONLY:  module_interface_check_data_output_surf,                                             &
              module_interface_data_output_surf,                                                   &
              module_interface_surface_data_averaging

#if defined( __netcdf )
   USE NETCDF
#endif

   USE netcdf_interface,                                                                           &
       ONLY:  nc_stat,                                                                             &
              netcdf_create_att,                                                                   &
              netcdf_create_dim,                                                                   &
              netcdf_create_global_atts,                                                           &
              netcdf_create_var,                                                                   &
              netcdf_data_format,                                                                  &
              netcdf_handle_error

   USE pegrid

   USE restart_data_mpi_io_mod,                                                                    &
       ONLY:  rrd_mpi_io,                                                                          &
              rd_mpi_io_check_array,                                                               &
              rrd_mpi_io_surface,                                                                  &
              rd_mpi_io_surface_filetypes,                                                         &
              wrd_mpi_io,                                                                          &
              wrd_mpi_io_surface

   USE surface_mod,                                                                                &
       ONLY:  ind_pav_green,                                                                       &
              ind_veg_wall,                                                                        &
              ind_wat_win,                                                                         &
              surf_out,                                                                            &
              surf_def,                                                                            &
              surf_lsm,                                                                            &
              surf_usm

   IMPLICIT NONE

   CHARACTER(LEN=100), DIMENSION(300)     ::  data_output_surf = ' '  !< namelist variable which describes the output variables
   CHARACTER(LEN=100), DIMENSION(0:1,300) ::  dosurf = ' '            !< internal variable which describes the output variables
                                                                      !< and separates averaged from non-averaged output
   CHARACTER(LEN=100), DIMENSION(0:1,300) ::  dosurf_unit = ' '       !< internal variable which holds the unit of the given output
                                                                      !< variable

   INTEGER(iwp) ::  average_count_surf = 0  !< number of ensemble members used for averaging
   INTEGER(iwp) ::  dosurf_no(0:1)     = 0  !< number of surface output quantities
#if defined( __netcdf4_parallel )
   INTEGER(iwp) ::  oldmode                 !< save old set-fill-mode of netcdf file (not needed, but required for routine call)

   INTEGER(iwp), DIMENSION(0:1) ::  dosurf_time_count = 0  !< count of output time steps
   INTEGER(iwp), DIMENSION(0:1) ::  id_dim_s_surf          !< netcdf ID for dimension s
   INTEGER(iwp), DIMENSION(0:1) ::  id_dim_time_surf       !< netcdf ID for dimension time
   INTEGER(iwp), DIMENSION(0:1) ::  id_set_surf            !< netcdf ID for file
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_azimuth_surf    !< netcdf ID for variable azimuth
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_time_surf       !< netcdf ID for variable time
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_s_surf          !< netcdf ID for variable s
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_xs_surf         !< netcdf ID for variable xs
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_ys_surf         !< netcdf ID for variable ys
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_zenith_surf     !< netcdf ID for variable zenith
   INTEGER(iwp), DIMENSION(0:1) ::  id_var_zs_surf         !< netcdf ID for variable zs
   INTEGER(iwp), DIMENSION(0:1) ::  ntdim_surf             !< number of output time steps

   INTEGER(iwp), DIMENSION(0:1,300) ::  id_var_dosurf      !< netcdf ID for output variables
#endif

   LOGICAL ::  first_output(0:1) = .FALSE.  !< true if first output was already called
   LOGICAL ::  to_netcdf         = .FALSE.  !< flag indicating parallel NetCDF output
   LOGICAL ::  to_vtk            = .FALSE.  !< flag indicating binary surface-data output that can be further
                                            !< processed to VTK format

   REAL(wp) ::  averaging_interval_surf  = 9999999.9_wp  !< averaging interval
   REAL(wp) ::  dt_dosurf                = 9999999.9_wp  !< time interval for instantaneous data output
   REAL(wp) ::  dt_dosurf_av             = 9999999.9_wp  !< time interval for averaged data output
   REAL(wp) ::  skip_time_dosurf         = 0.0_wp        !< skip time for instantaneous data output
   REAL(wp) ::  skip_time_dosurf_av      = 0.0_wp        !< skip time for averaged data output
   REAL(wp) ::  time_dosurf              = 0.0_wp        !< internal counter variable to check for instantaneous data output
   REAL(wp) ::  time_dosurf_av           = 0.0_wp        !< internal counter variable to check for averaged data output

   SAVE

   PRIVATE

   INTERFACE  surface_data_output
      MODULE PROCEDURE surface_data_output
   END INTERFACE  surface_data_output

   INTERFACE  surface_data_output_averaging
      MODULE PROCEDURE surface_data_output_averaging
   END INTERFACE  surface_data_output_averaging

   INTERFACE  surface_data_output_check_parameters
      MODULE PROCEDURE surface_data_output_check_parameters
   END INTERFACE  surface_data_output_check_parameters

   INTERFACE  surface_data_output_init
      MODULE PROCEDURE surface_data_output_init
   END INTERFACE  surface_data_output_init

   INTERFACE  surface_data_output_init_arrays
      MODULE PROCEDURE surface_data_output_init_arrays
   END INTERFACE  surface_data_output_init_arrays

   INTERFACE  surface_data_output_last_action
      MODULE PROCEDURE surface_data_output_last_action
   END INTERFACE  surface_data_output_last_action

   INTERFACE  surface_data_output_parin
      MODULE PROCEDURE surface_data_output_parin
   END INTERFACE  surface_data_output_parin

   INTERFACE  surface_data_output_rrd_global
      MODULE PROCEDURE surface_data_output_rrd_global_ftn
      MODULE PROCEDURE surface_data_output_rrd_global_mpi
   END INTERFACE  surface_data_output_rrd_global

   INTERFACE  surface_data_output_rrd_local
      MODULE PROCEDURE surface_data_output_rrd_local_ftn
      MODULE PROCEDURE surface_data_output_rrd_local_mpi
   END INTERFACE  surface_data_output_rrd_local

   INTERFACE  surface_data_output_wrd_global
      MODULE PROCEDURE surface_data_output_wrd_global
   END INTERFACE  surface_data_output_wrd_global

   INTERFACE  surface_data_output_wrd_local
      MODULE PROCEDURE surface_data_output_wrd_local
   END INTERFACE  surface_data_output_wrd_local

   INTERFACE  surface_data_output_sum_up
      MODULE PROCEDURE surface_data_output_sum_up_1d
      MODULE PROCEDURE surface_data_output_sum_up_2d
   END INTERFACE  surface_data_output_sum_up

   INTERFACE  surface_data_output_collect
      MODULE PROCEDURE surface_data_output_collect_1d
      MODULE PROCEDURE surface_data_output_collect_2d
   END INTERFACE  surface_data_output_collect

!
!--Public subroutines
   PUBLIC surface_data_output,                                                                     &
          surface_data_output_averaging,                                                           &
          surface_data_output_check_parameters,                                                    &
          surface_data_output_init,                                                                &
          surface_data_output_init_arrays,                                                         &
          surface_data_output_last_action,                                                         &
          surface_data_output_parin,                                                               &
          surface_data_output_rrd_global,                                                          &
          surface_data_output_rrd_local,                                                           &
          surface_data_output_wrd_local,                                                           &
          surface_data_output_wrd_global
!
!--Public variables
   PUBLIC average_count_surf,                                                                      &
          averaging_interval_surf,                                                                 &
          dt_dosurf,                                                                               &
          dt_dosurf_av,                                                                            &
          skip_time_dosurf,                                                                        &
          skip_time_dosurf_av,                                                                     &
          time_dosurf,                                                                             &
          time_dosurf_av

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine counts the number of surfaces on each core and allocates arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_init_arrays

    IMPLICIT NONE

!
!-- Determine the number of surface elements on subdomain
    surf_out%ns = surf_def%ns + surf_lsm%ns + surf_usm%ns
!
!--  Determine the total number of surfaces in the model domain
#if defined( __parallel )
     CALL MPI_ALLREDUCE( surf_out%ns, surf_out%ns_total, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
     surf_out%ns_total = surf_out%ns
#endif
!
!-- Allocate output variable and set to _FillValue attribute
    ALLOCATE ( surf_out%var_out(1:surf_out%ns) )
    surf_out%var_out = surf_out%fillvalue
!
!-- If there is an output of time average output variables, allocate the required array.
    IF ( dosurf_no(1) > 0 )  THEN
       ALLOCATE ( surf_out%var_av(1:surf_out%ns,1:dosurf_no(1)) )
       surf_out%var_av = 0.0_wp
    ENDIF

 END SUBROUTINE surface_data_output_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization surface-data output data structure: calculation of vertices and polygon data for
!> the surface elements, allocation of required arrays.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_init

    IMPLICIT NONE

#if defined( __netcdf4_parallel )
    CHARACTER (LEN=100)  ::  filename           !< name of output file
    CHARACTER (LEN=80)   ::  time_average_text  !< string written to file attribute time_avg
    CHARACTER (LEN=4000) ::  var_list           !< list of variables written to NetCDF file

    INTEGER(iwp) ::  av                 !< flag for averaged (=1) and non-averaged (=0) data
#endif
    INTEGER(iwp) ::  i                  !< grid index in x-direction, also running variable for counting non-average data output
    INTEGER(iwp) ::  j                  !< grid index in y-direction, also running variable for counting average data output
    INTEGER(iwp) ::  k                  !< grid index in z-direction
    INTEGER(iwp) ::  kw                 !< surface grid index on w-grid
    INTEGER(iwp) ::  m                  !< running index for surface elements
    INTEGER(iwp) ::  mm                 !< local counting variable for surface elements
    INTEGER(iwp) ::  npg                !< counter variable for all surface elements ( or polygons )
    INTEGER(iwp) ::  point_index_count  !< local counter variable for point index
    INTEGER(iwp) ::  start_count        !< local start counter for the surface index

    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_points_on_pe    !< array which contains the number of points on all mpi ranks
    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_surfaces_on_pe  !< array which contains the number of surfaces on all mpi ranks
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:,:) ::  point_index   !< dummy array used to check where the reference points for
                                                                  !< surface polygons are located

    REAL(wp) ::  az     !< azimuth angle, indicated the vertical orientation of a surface element
    REAL(wp) ::  off_x  !< grid offset in x-direction between the stored grid index and the actual wall
    REAL(wp) ::  off_y  !< grid offset in y-direction between the stored grid index and the actual wall
    REAL(wp) ::  ze     !< zenith angle, indicated the horizontal orientation of a surface element
    REAL(wp) ::  zpos   !< z-position of a surface element

#if defined( __netcdf4_parallel )
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  netcdf_data_1d  !< dummy array to output 1D data into netcdf file
#endif

!
!-- If output to VTK format is enabled, initialize point and polygon data.
!-- In a first step, count the number of points which are defining the surfaces and the polygons.
    IF ( to_vtk )  THEN
       ALLOCATE( point_index(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1) )
       point_index = -1
!
!--    Default surfaces
       surf_out%npoints = 0
       DO  m = 1, surf_def%ns
!
!--       Determine the indices of the respective grid cell inside the topography. Please note,
!--       j-index for north-facing surfaces is identical to the reference j-index outside the
!--       grid box (surface is defined at ((i-0.5)*dx,(j-0.5)*dy). Equivalent for east-facing
!--       surfaces and i-index. This means that ioff and joff are only added if their values is 1.
          i = surf_def%i(m) + MERGE( surf_def%ioff(m), 0, surf_def%westward(m)  )
          j = surf_def%j(m) + MERGE( surf_def%joff(m), 0, surf_def%southward(m) )
          k = surf_def%k(m) + surf_def%koff(m)
!
!--       Check if the vertices that define the surface element are already defined, if not,
!--       increment the counter. Check vertices that define a horizontal surface.
          IF ( surf_def%upward(m)  .OR.  surf_def%downward(m) )  THEN
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = surf_out%npoints - 1
             ENDIF
             IF ( point_index(k,j,i+1) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j,i+1) = surf_out%npoints - 1
             ENDIF
             IF ( point_index(k,j+1,i+1) < 0 )  THEN
                surf_out%npoints       = surf_out%npoints + 1
                point_index(k,j+1,i+1) = surf_out%npoints - 1
             ENDIF
             IF ( point_index(k,j+1,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j+1,i) = surf_out%npoints - 1
             ENDIF
!
!--       Now check for vertices that define a vertical surface.
          ELSE
!
!--          Lower left /front vertex
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = surf_out%npoints - 1
             ENDIF
!
!--          Upper left / front vertex
             IF ( point_index(k+1,j,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k+1,j,i) = surf_out%npoints - 1
             ENDIF
!
!--          Upper / lower right index for north- and south-facing surfaces
             IF ( surf_def%northward(m)  .OR.  surf_def%southward(m) )  THEN
                IF ( point_index(k,j,i+1) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j,i+1) = surf_out%npoints - 1
                ENDIF
                IF ( point_index(k+1,j,i+1) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j,i+1) = surf_out%npoints - 1
                ENDIF
!
!--          Upper / lower front index for east- and west-facing surfaces
             ELSEIF ( surf_def%eastward(m)  .OR.  surf_def%westward(m) )  THEN
                IF ( point_index(k,j+1,i) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j+1,i) = surf_out%npoints - 1
                ENDIF
                IF ( point_index(k+1,j+1,i) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j+1,i) = surf_out%npoints - 1
                ENDIF
             ENDIF

          ENDIF
       ENDDO
       DO  m = 1, surf_lsm%ns
!
!--       Determine the indices of the respective grid cell inside the topography. Please note,
!--       j-index for north-facing surfaces is identical to the reference j-index outside the
!--       grid box (surface is defined at ((i-0.5)*dx,(j-0.5)*dy). Equivalent for east-facing
!--       surfaces and i-index. This means that ioff and joff are only added if their values is 1.
          i = surf_lsm%i(m) + MERGE( surf_lsm%ioff(m), 0, surf_lsm%westward(m)  )
          j = surf_lsm%j(m) + MERGE( surf_lsm%joff(m), 0, surf_lsm%southward(m) )
          k = surf_lsm%k(m) + surf_lsm%koff(m)
!
!--       Check if the vertices that define the surface element are already defined, if not,
!--       increment the counter. Check vertices that define a horizontal surface.
          IF ( surf_lsm%upward(m)  .OR.  surf_lsm%downward(m) )  THEN
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = surf_out%npoints - 1
             ENDIF
             IF ( point_index(k,j,i+1) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j,i+1) = surf_out%npoints - 1
             ENDIF
             IF ( point_index(k,j+1,i+1) < 0 )  THEN
                surf_out%npoints       = surf_out%npoints + 1
                point_index(k,j+1,i+1) = surf_out%npoints - 1
             ENDIF
             IF ( point_index(k,j+1,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j+1,i) = surf_out%npoints - 1
             ENDIF
!
!--       Now check for vertices that define a vertical surface.
          ELSE
!
!--          Lower left /front vertex
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = surf_out%npoints - 1
             ENDIF
!
!--          Upper left / front vertex
             IF ( point_index(k+1,j,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k+1,j,i) = surf_out%npoints - 1
             ENDIF
!
!--          Upper / lower right index for north- and south-facing surfaces
             IF ( surf_lsm%northward(m)  .OR.  surf_lsm%southward(m) )  THEN
                IF ( point_index(k,j,i+1) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j,i+1) = surf_out%npoints - 1
                ENDIF
                IF ( point_index(k+1,j,i+1) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j,i+1) = surf_out%npoints - 1
                ENDIF
!
!--          Upper / lower front index for east- and west-facing surfaces
             ELSEIF ( surf_lsm%eastward(m)  .OR.  surf_lsm%westward(m) )  THEN
                IF ( point_index(k,j+1,i) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j+1,i) = surf_out%npoints - 1
                ENDIF
                IF ( point_index(k+1,j+1,i) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j+1,i) = surf_out%npoints - 1
                ENDIF
             ENDIF

          ENDIF
       ENDDO
       DO  m = 1, surf_usm%ns
!
!--       Determine the indices of the respective grid cell inside the topography. Please note,
!--       j-index for north-facing surfaces is identical to the reference j-index outside the
!--       grid box (surface is defined at ((i-0.5)*dx,(j-0.5)*dy). Equivalent for east-facing
!--       surfaces and i-index. This means that ioff and joff are only added if their values is 1.
          i = surf_usm%i(m) + MERGE( surf_usm%ioff(m), 0, surf_usm%westward(m)  )
          j = surf_usm%j(m) + MERGE( surf_usm%joff(m), 0, surf_usm%southward(m) )
          k = surf_usm%k(m) + surf_usm%koff(m)
!
!--       Check if the vertices that define the surface element are already defined, if not,
!--       increment the counter. Check vertices that define a horizontal surface.
          IF ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) )  THEN
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = surf_out%npoints - 1
             ENDIF
             IF ( point_index(k,j,i+1) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j,i+1) = surf_out%npoints - 1
             ENDIF
             IF ( point_index(k,j+1,i+1) < 0 )  THEN
                surf_out%npoints       = surf_out%npoints + 1
                point_index(k,j+1,i+1) = surf_out%npoints - 1
             ENDIF
             IF ( point_index(k,j+1,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j+1,i) = surf_out%npoints - 1
             ENDIF
!
!--       Now check for vertices that define a vertical surface.
          ELSE
!
!--          Lower left /front vertex
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = surf_out%npoints - 1
             ENDIF
!
!--          Upper left / front vertex
             IF ( point_index(k+1,j,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k+1,j,i) = surf_out%npoints - 1
             ENDIF
!
!--          Upper / lower right index for north- and south-facing surfaces
             IF ( surf_usm%northward(m)  .OR.  surf_usm%southward(m) )  THEN
                IF ( point_index(k,j,i+1) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j,i+1) = surf_out%npoints - 1
                ENDIF
                IF ( point_index(k+1,j,i+1) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j,i+1) = surf_out%npoints - 1
                ENDIF
!
!--          Upper / lower front index for east- and west-facing surfaces
             ELSEIF ( surf_usm%eastward(m)  .OR.  surf_usm%westward(m) )  THEN
                IF ( point_index(k,j+1,i) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j+1,i) = surf_out%npoints - 1
                ENDIF
                IF ( point_index(k+1,j+1,i) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j+1,i) = surf_out%npoints - 1
                ENDIF
             ENDIF

          ENDIF
       ENDDO
!
!--    Allocate the number of points and polygons. Note, the number of polygons is identical to the
!--    number of surfaces elements, whereas the number of points (vertices), which define the
!--    polygons, can be larger.
       ALLOCATE( surf_out%points(3,1:surf_out%npoints) )
       ALLOCATE( surf_out%polygons(5,1:surf_out%ns)    )
!
!--    Note, PARAVIEW expects consecutively ordered points, in order to unambiguously identify
!--    surfaces. Hence, all PEs should know where they start counting, depending on the number of
!--    points on the other PE's with lower MPI rank.
#if defined( __parallel )
       CALL MPI_ALLGATHER( surf_out%npoints, 1, MPI_INTEGER, num_points_on_pe, 1, MPI_INTEGER,     &
                           comm2d, ierr  )
#else
       num_points_on_pe = surf_out%npoints
#endif

!
!--    After the number of vertices is counted, repeat the loops and define the vertices. Start with
!--    the horizontal default surfaces. First, however, determine the offset where couting of points
!--    should be started, which is the sum of points of all PE's with lower MPI rank.
       i                 = 0
       point_index_count = 0
       DO WHILE  ( i < myid  .AND.  i <= SIZE( num_points_on_pe ) )
          point_index_count = point_index_count + num_points_on_pe(i)
          i                 = i + 1
       ENDDO

       surf_out%npoints = 0
       point_index      = -1
       npg              = 0

       DO  m = 1, surf_def%ns
!
!--       Determine the indices of the respective grid cell inside the topography. Please note,
!--       j-index for north-facing surfaces is identical to the reference j-index outside the
!--       grid box (surface is defined at ((i-0.5)*dx,(j-0.5)*dy). Equivalent for east-facing
!--       surfaces and i-index. This means that ioff and joff are only added if their values is 1.
          i  = surf_def%i(m) + MERGE( surf_def%ioff(m), 0, surf_def%westward(m) )
          j  = surf_def%j(m) + MERGE( surf_def%joff(m), 0, surf_def%southward(m) )
          k  = surf_def%k(m) + surf_def%koff(m)
          kw = MERGE( k, surf_def%k(m), surf_def%upward(m) )
!
!--       Check if the vertices that define the surface element are already defined, if not,
!--       increment the counter. Check vertices that define a horizontal surface.
          IF ( surf_def%upward(m)  .OR.  surf_def%downward(m) )  THEN
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = point_index_count
                point_index_count  = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             IF ( point_index(k,j,i+1) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j,i+1) = point_index_count
                point_index_count    = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j     - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             IF ( point_index(k,j+1,i+1) < 0 )  THEN
                surf_out%npoints       = surf_out%npoints + 1
                point_index(k,j+1,i+1) = point_index_count
                point_index_count      = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             IF ( point_index(k,j+1,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j+1,i) = point_index_count
                point_index_count    = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i     - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             npg                      = npg + 1
             surf_out%polygons(1,npg) = 4
             surf_out%polygons(2,npg) = point_index(k,j,i)
             surf_out%polygons(3,npg) = point_index(k,j,i+1)
             surf_out%polygons(4,npg) = point_index(k,j+1,i+1)
             surf_out%polygons(5,npg) = point_index(k,j+1,i)
!
!--       Now check for vertices that define a vertical surface.
          ELSE
!
!--          Lower left /front vertex
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = point_index_count
                point_index_count  = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(k-1)
             ENDIF
!
!--          Upper left / front vertex
             IF ( point_index(k+1,j,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k+1,j,i) = point_index_count
                point_index_count    = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(k)
             ENDIF
!
!--          Upper / lower right index for north- and south-facing surfaces
             IF ( surf_def%northward(m)  .OR.  surf_def%southward(m) )  THEN
                IF ( point_index(k,j,i+1) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j,i+1) = point_index_count
                   point_index_count    = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j     - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k-1)
                ENDIF
                IF ( point_index(k+1,j,i+1) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j,i+1) = point_index_count
                   point_index_count      = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j     - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k)
                ENDIF

                npg                      = npg + 1
                surf_out%polygons(1,npg) = 4
                surf_out%polygons(2,npg) = point_index(k,j,i)
                surf_out%polygons(3,npg) = point_index(k,j,i+1)
                surf_out%polygons(4,npg) = point_index(k+1,j,i+1)
                surf_out%polygons(5,npg) = point_index(k+1,j,i)
!
!--          Upper / lower front index for east- and west-facing surfaces
             ELSEIF ( surf_def%northward(m)  .OR.  surf_def%southward(m) )  THEN
                IF ( point_index(k,j+1,i) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j+1,i) = point_index_count
                   point_index_count    = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i     - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k-1)
                ENDIF
                IF ( point_index(k+1,j+1,i) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j+1,i) = point_index_count
                   point_index_count      = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i     - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k)
                ENDIF

                npg                      = npg + 1
                surf_out%polygons(1,npg) = 4
                surf_out%polygons(2,npg) = point_index(k,j,i)
                surf_out%polygons(3,npg) = point_index(k,j+1,i)
                surf_out%polygons(4,npg) = point_index(k+1,j+1,i)
                surf_out%polygons(5,npg) = point_index(k+1,j,i)
             ENDIF
          ENDIF
       ENDDO
       DO  m = 1, surf_lsm%ns
!
!--       Determine the indices of the respective grid cell inside the topography. Please note,
!--       j-index for north-facing surfaces is identical to the reference j-index outside the
!--       grid box (surface is defined at ((i-0.5)*dx,(j-0.5)*dy). Equivalent for east-facing
!--       surfaces and i-index. This means that ioff and joff are only added if their values is 1.
          i  = surf_lsm%i(m) + MERGE( surf_lsm%ioff(m), 0, surf_lsm%westward(m) )
          j  = surf_lsm%j(m) + MERGE( surf_lsm%joff(m), 0, surf_lsm%southward(m) )
          k  = surf_lsm%k(m) + surf_lsm%koff(m)
          kw = MERGE( k, surf_lsm%k(m), surf_lsm%upward(m) )
!
!--       Check if the vertices that define the surface element are already defined, if not,
!--       increment the counter. Check vertices that define a horizontal surface.
          IF ( surf_lsm%upward(m)  .OR.  surf_lsm%downward(m) )  THEN
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = point_index_count
                point_index_count  = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             IF ( point_index(k,j,i+1) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j,i+1) = point_index_count
                point_index_count    = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j     - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             IF ( point_index(k,j+1,i+1) < 0 )  THEN
                surf_out%npoints       = surf_out%npoints + 1
                point_index(k,j+1,i+1) = point_index_count
                point_index_count      = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             IF ( point_index(k,j+1,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j+1,i) = point_index_count
                point_index_count    = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i     - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF

             npg                      = npg + 1
             surf_out%polygons(1,npg) = 4
             surf_out%polygons(2,npg) = point_index(k,j,i)
             surf_out%polygons(3,npg) = point_index(k,j,i+1)
             surf_out%polygons(4,npg) = point_index(k,j+1,i+1)
             surf_out%polygons(5,npg) = point_index(k,j+1,i)
!
!--       Now check for vertices that define a vertical surface.
          ELSE
!
!--          Lower left /front vertex
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = point_index_count
                point_index_count  = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(k-1)
             ENDIF
!
!--          Upper left / front vertex
             IF ( point_index(k+1,j,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k+1,j,i) = point_index_count
                point_index_count    = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(k)
             ENDIF
!
!--          Upper / lower right index for north- and south-facing surfaces
             IF ( surf_lsm%northward(m)  .OR.  surf_lsm%southward(m) )  THEN
                IF ( point_index(k,j,i+1) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j,i+1) = point_index_count
                   point_index_count    = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j     - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k-1)
                ENDIF
                IF ( point_index(k+1,j,i+1) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j,i+1) = point_index_count
                   point_index_count      = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j     - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k)
                ENDIF

                npg                      = npg + 1
                surf_out%polygons(1,npg) = 4
                surf_out%polygons(2,npg) = point_index(k,j,i)
                surf_out%polygons(3,npg) = point_index(k,j,i+1)
                surf_out%polygons(4,npg) = point_index(k+1,j,i+1)
                surf_out%polygons(5,npg) = point_index(k+1,j,i)
!
!--          Upper / lower front index for east- and west-facing surfaces
             ELSEIF ( surf_lsm%northward(m)  .OR.  surf_lsm%southward(m) )  THEN
                IF ( point_index(k,j+1,i) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j+1,i) = point_index_count
                   point_index_count    = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i     - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k-1)
                ENDIF
                IF ( point_index(k+1,j+1,i) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j+1,i) = point_index_count
                   point_index_count      = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i     - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k)
                ENDIF

                npg                      = npg + 1
                surf_out%polygons(1,npg) = 4
                surf_out%polygons(2,npg) = point_index(k,j,i)
                surf_out%polygons(3,npg) = point_index(k,j+1,i)
                surf_out%polygons(4,npg) = point_index(k+1,j+1,i)
                surf_out%polygons(5,npg) = point_index(k+1,j,i)
             ENDIF
          ENDIF
       ENDDO

       DO  m = 1, surf_usm%ns
!
!--       Determine the indices of the respective grid cell inside the topography. Please note,
!--       j-index for north-facing surfaces is identical to the reference j-index outside the
!--       grid box (surface is defined at ((i-0.5)*dx,(j-0.5)*dy). Equivalent for east-facing
!--       surfaces and i-index. This means that ioff and joff are only added if their values is 1.
          i  = surf_usm%i(m) + MERGE( surf_usm%ioff(m), 0, surf_usm%westward(m) )
          j  = surf_usm%j(m) + MERGE( surf_usm%joff(m), 0, surf_usm%southward(m) )
          k  = surf_usm%k(m) + surf_usm%koff(m)
          kw = MERGE( k, surf_usm%k(m), surf_usm%upward(m) )
!
!--       Check if the vertices that define the surface element are already defined, if not,
!--       increment the counter. Check vertices that define a horizontal surface.
          IF ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) )  THEN
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = point_index_count
                point_index_count  = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             IF ( point_index(k,j,i+1) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j,i+1) = point_index_count
                point_index_count    = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j     - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             IF ( point_index(k,j+1,i+1) < 0 )  THEN
                surf_out%npoints       = surf_out%npoints + 1
                point_index(k,j+1,i+1) = point_index_count
                point_index_count      = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             IF ( point_index(k,j+1,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k,j+1,i) = point_index_count
                point_index_count    = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i     - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(kw)
             ENDIF
             npg                      = npg + 1
             surf_out%polygons(1,npg) = 4
             surf_out%polygons(2,npg) = point_index(k,j,i)
             surf_out%polygons(3,npg) = point_index(k,j,i+1)
             surf_out%polygons(4,npg) = point_index(k,j+1,i+1)
             surf_out%polygons(5,npg) = point_index(k,j+1,i)
!
!--       Now check for vertices that define a vertical surface.
          ELSE
!
!--          Lower left /front vertex
             IF ( point_index(k,j,i) < 0 )  THEN
                surf_out%npoints   = surf_out%npoints + 1
                point_index(k,j,i) = point_index_count
                point_index_count  = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(k-1)
             ENDIF
!
!--          Upper left / front vertex
             IF ( point_index(k+1,j,i) < 0 )  THEN
                surf_out%npoints     = surf_out%npoints + 1
                point_index(k+1,j,i) = point_index_count
                point_index_count    = point_index_count + 1
                surf_out%points(1,surf_out%npoints) = ( i - 0.5_wp ) * dx
                surf_out%points(2,surf_out%npoints) = ( j - 0.5_wp ) * dy
                surf_out%points(3,surf_out%npoints) = zw(k)
             ENDIF
!
!--          Upper / lower right index for north- and south-facing surfaces
             IF ( surf_usm%northward(m)  .OR.  surf_usm%southward(m) )  THEN
                IF ( point_index(k,j,i+1) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j,i+1) = point_index_count
                   point_index_count    = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j     - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k-1)
                ENDIF
                IF ( point_index(k+1,j,i+1) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j,i+1) = point_index_count
                   point_index_count      = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i + 1 - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j     - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k)
                ENDIF

                npg                      = npg + 1
                surf_out%polygons(1,npg) = 4
                surf_out%polygons(2,npg) = point_index(k,j,i)
                surf_out%polygons(3,npg) = point_index(k,j,i+1)
                surf_out%polygons(4,npg) = point_index(k+1,j,i+1)
                surf_out%polygons(5,npg) = point_index(k+1,j,i)
!
!--          Upper / lower front index for east- and west-facing surfaces
             ELSEIF ( surf_usm%northward(m)  .OR.  surf_usm%southward(m) )  THEN
                IF ( point_index(k,j+1,i) < 0 )  THEN
                   surf_out%npoints     = surf_out%npoints + 1
                   point_index(k,j+1,i) = point_index_count
                   point_index_count    = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i     - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k-1)
                ENDIF
                IF ( point_index(k+1,j+1,i) < 0 )  THEN
                   surf_out%npoints       = surf_out%npoints + 1
                   point_index(k+1,j+1,i) = point_index_count
                   point_index_count      = point_index_count + 1
                   surf_out%points(1,surf_out%npoints) = ( i     - 0.5_wp ) * dx
                   surf_out%points(2,surf_out%npoints) = ( j + 1 - 0.5_wp ) * dy
                   surf_out%points(3,surf_out%npoints) = zw(k)
                ENDIF

                npg                      = npg + 1
                surf_out%polygons(1,npg) = 4
                surf_out%polygons(2,npg) = point_index(k,j,i)
                surf_out%polygons(3,npg) = point_index(k,j+1,i)
                surf_out%polygons(4,npg) = point_index(k+1,j+1,i)
                surf_out%polygons(5,npg) = point_index(k+1,j,i)
             ENDIF
          ENDIF
       ENDDO
!
!--    Deallocate temporary dummy variable
       DEALLOCATE ( point_index )
!
!--    Sum-up total number of vertices on domain. This will be needed for post-processing.
       surf_out%npoints_total = 0
#if defined( __parallel )
        CALL MPI_ALLREDUCE( surf_out%npoints, surf_out%npoints_total, 1, MPI_INTEGER, MPI_SUM,     &
                            comm2d, ierr )
#else
        surf_out%npoints_total = surf_out%npoints
#endif
     ENDIF
!
!--  If output to netcdf is enabled, set-up the coordinate arrays that unambiguously describe the
!--  position and orientation of each surface element.
     IF ( to_netcdf )  THEN
!
!--     Allocate local coordinate arrays
        ALLOCATE( surf_out%s(1:surf_out%ns) )
        ALLOCATE( surf_out%xs(1:surf_out%ns) )
        ALLOCATE( surf_out%ys(1:surf_out%ns) )
        ALLOCATE( surf_out%zs(1:surf_out%ns) )
        ALLOCATE( surf_out%azimuth(1:surf_out%ns) )
        ALLOCATE( surf_out%zenith(1:surf_out%ns) )
!
!--     Gather the number of surface on each processor, in order to number the surface elements in
!--     ascending order with respect to the total number of surfaces in the domain.
#if defined( __parallel )
        CALL MPI_ALLGATHER( surf_out%ns, 1, MPI_INTEGER, num_surfaces_on_pe, 1, MPI_INTEGER,       &
                            comm2d, ierr  )
#else
        num_surfaces_on_pe = surf_out%ns
#endif
!
!--     First, however, determine the offset where couting of the surfaces should start (the sum of
!--     surfaces on all PE's with lower MPI rank).
        i           = 0
        start_count = 1
        DO WHILE  ( i < myid  .AND.  i <= SIZE( num_surfaces_on_pe ) )
           start_count = start_count + num_surfaces_on_pe(i)
           i           = i + 1
        ENDDO
!
!--     Set coordinate arrays. For horizontal surfaces, azimuth angles are not defined (fill value).
!--     Zenith angle is 0 (180) for upward (downward)-facing surfaces, respectively.
!--     Azimuth angle: northward (0), eastward (90), southward (180), westward (270).
!--     For vertical surfaces, zenith angles are 90.
        i  = start_count
        mm = 1
        DO  m = 1, surf_def%ns
!
!--        Pre-define azimuth, zenith and offset values.
!--        Uward-facing
           IF ( surf_def%upward(m) )  THEN
              az    = surf_out%fillvalue
              ze    = 0.0_wp
              off_x = 0.5_wp
              off_y = 0.5_wp
              zpos  = zw(surf_def%k(m)+surf_def%koff(m))
!
!--        Downward-facing
           ELSEIF ( surf_def%downward(m) )  THEN
              az    = surf_out%fillvalue
              ze    = 180.0_wp
              off_x = 0.5_wp
              off_y = 0.5_wp
!
!--           In contrast to upward-facing surfaces, no offset is added.
              zpos  = zw(surf_def%k(m))
!
!--        Northward-facing
           ELSEIF ( surf_def%northward(m) )  THEN
              az    = 0.0_wp
              ze    = 90.0_wp
              off_x = 0.5_wp
              off_y = 0.0_wp
              zpos  = zu(surf_def%k(m))
!
!--        Southward-facing
           ELSEIF ( surf_def%southward(m) )  THEN
              az    = 180.0_wp
              ze    = 90.0_wp
              off_x = 0.5_wp
              off_y = 1.0_wp
              zpos  = zu(surf_def%k(m))
!
!--        Eastward-facing
           ELSEIF ( surf_def%eastward(m) )  THEN
              az    = 90.0_wp
              ze    = 90.0_wp
              off_x = 0.0_wp
              off_y = 0.5_wp
              zpos  = zu(surf_def%k(m))
!
!--        Westward-facing
           ELSEIF ( surf_def%westward(m) )  THEN
              az    = 270.0_wp
              ze    = 90.0_wp
              off_x = 1.0_wp
              off_y = 0.5_wp
              zpos  = zu(surf_def%k(m))
           ENDIF

           surf_out%s(mm)       = i
           surf_out%xs(mm)      = ( surf_def%i(m) + off_x ) * dx
           surf_out%ys(mm)      = ( surf_def%j(m) + off_y ) * dy
           surf_out%zs(mm)      = zpos
           surf_out%azimuth(mm) = az
           surf_out%zenith(mm)  = ze
           i                    = i + 1
           mm                   = mm + 1
        ENDDO
        DO  m = 1, surf_lsm%ns
!
!--        Pre-define azimuth, zenith and offset values.
!--        Uward-facing
           IF ( surf_lsm%upward(m) )  THEN
              az    = surf_out%fillvalue
              ze    = 0.0_wp
              off_x = 0.5_wp
              off_y = 0.5_wp
              zpos  = zw(surf_lsm%k(m)+surf_lsm%koff(m))
!
!--        Downward-facing
           ELSEIF ( surf_lsm%downward(m) )  THEN
              az    = surf_out%fillvalue
              ze    = 180.0_wp
              off_x = 0.5_wp
              off_y = 0.5_wp
!
!--           In contrast to upward-facing surfaces, no offset is added.
              zpos  = zw(surf_lsm%k(m))
!
!--        Northward-facing
           ELSEIF ( surf_lsm%northward(m) )  THEN
              az    = 0.0_wp
              ze    = 90.0_wp
              off_x = 0.5_wp
              off_y = 0.0_wp
              zpos  = zu(surf_lsm%k(m))
!
!--        Southward-facing
           ELSEIF ( surf_lsm%southward(m) )  THEN
              az    = 180.0_wp
              ze    = 90.0_wp
              off_x = 0.5_wp
              off_y = 1.0_wp
              zpos  = zu(surf_lsm%k(m))
!
!--        Eastward-facing
           ELSEIF ( surf_lsm%eastward(m) )  THEN
              az    = 90.0_wp
              ze    = 90.0_wp
              off_x = 0.0_wp
              off_y = 0.5_wp
              zpos  = zu(surf_lsm%k(m))
!
!--        Westward-facing
           ELSEIF ( surf_lsm%westward(m) )  THEN
              az    = 270.0_wp
              ze    = 90.0_wp
              off_x = 1.0_wp
              off_y = 0.5_wp
              zpos  = zu(surf_lsm%k(m))
           ENDIF

           surf_out%s(mm)       = i
           surf_out%xs(mm)      = ( surf_lsm%i(m) + off_x ) * dx
           surf_out%ys(mm)      = ( surf_lsm%j(m) + off_y ) * dy
           surf_out%zs(mm)      = zpos
           surf_out%azimuth(mm) = az
           surf_out%zenith(mm)  = ze
           i                    = i + 1
           mm                   = mm + 1
        ENDDO
        DO  m = 1, surf_usm%ns
!
!--        Pre-define azimuth, zenith and offset values.
!--        Uward-facing
           IF ( surf_usm%upward(m) )  THEN
              az    = surf_out%fillvalue
              ze    = 0.0_wp
              off_x = 0.5_wp
              off_y = 0.5_wp
              zpos  = zw(surf_usm%k(m)+surf_usm%koff(m))
!
!--        Downward-facing
           ELSEIF ( surf_usm%downward(m) )  THEN
              az    = surf_out%fillvalue
              ze    = 180.0_wp
              off_x = 0.5_wp
              off_y = 0.5_wp
!
!--           In contrast to upward-facing surfaces, no offset is added.
              zpos  = zw(surf_usm%k(m))
!
!--        Northward-facing
           ELSEIF ( surf_usm%northward(m) )  THEN
              az    = 0.0_wp
              ze    = 90.0_wp
              off_x = 0.5_wp
              off_y = 0.0_wp
              zpos  = zu(surf_usm%k(m))
!
!--        Southward-facing
           ELSEIF ( surf_usm%southward(m) )  THEN
              az    = 180.0_wp
              ze    = 90.0_wp
              off_x = 0.5_wp
              off_y = 1.0_wp
              zpos  = zu(surf_usm%k(m))
!
!--        Eastward-facing
           ELSEIF ( surf_usm%eastward(m) )  THEN
              az    = 90.0_wp
              ze    = 90.0_wp
              off_x = 0.0_wp
              off_y = 0.5_wp
              zpos  = zu(surf_usm%k(m))
!
!--        Westward-facing
           ELSEIF ( surf_usm%westward(m) )  THEN
              az    = 270.0_wp
              ze    = 90.0_wp
              off_x = 1.0_wp
              off_y = 0.5_wp
              zpos  = zu(surf_usm%k(m))
           ENDIF

           surf_out%s(mm)       = i
           surf_out%xs(mm)      = ( surf_usm%i(m) + off_x ) * dx
           surf_out%ys(mm)      = ( surf_usm%j(m) + off_y ) * dy
           surf_out%zs(mm)      = zpos
           surf_out%azimuth(mm) = az
           surf_out%zenith(mm)  = ze
           i                    = i + 1
           mm                   = mm + 1
        ENDDO
!
!--     Initialize NetCDF data output. Please note, local start position for the surface elements in
!--     the NetCDF file is surfaces%s(1), while the number of surfaces on the subdomain is given by
!--     surfaces%ns.
#if defined( __netcdf4_parallel )

!
!--     Calculate number of time steps to be output
        ntdim_surf(0) = dosurf_time_count(0) + CEILING( ( end_time - MAX(                          &
                          MERGE( skip_time_dosurf, skip_time_dosurf + spinup_time,                 &
                                 data_output_during_spinup ), simulated_time_at_begin )            &
                                                        ) / dt_dosurf )

        ntdim_surf(1) = dosurf_time_count(1) + CEILING( ( end_time - MAX(                          &
                          MERGE( skip_time_dosurf_av, skip_time_dosurf_av + spinup_time,           &
                                 data_output_during_spinup ), simulated_time_at_begin )            &
                                                        ) / dt_dosurf_av )

!
!--     Create NetCDF4 files for parallel writing
        DO  av = 0, 1
!
!--        If there is no instantaneous data (av=0) or averaged data (av=1) requested, do not create
!--        the corresponding NetCDF file
           IF ( dosurf_no(av) == 0 ) CYCLE

           IF ( av == 0 )  THEN
              filename = 'SURFACE_DATA_NETCDF' // TRIM( coupling_char )
           ELSE
              filename = 'SURFACE_DATA_AV_NETCDF' // TRIM( coupling_char )
           ENDIF
!
!--        Open file using netCDF4/HDF5 format, parallel
           nc_stat = NF90_CREATE( TRIM(filename),                                                  &
                                  IOR( NF90_NOCLOBBER, IOR( NF90_NETCDF4, NF90_MPIIO ) ),          &
                                  id_set_surf(av), COMM = comm2d, INFO = MPI_INFO_NULL )
           CALL netcdf_handle_error( 'surface_data_output_mod', 5550 )

           !- Write some global attributes
           IF ( av == 0 )  THEN
              CALL netcdf_create_global_atts( id_set_surf(av), 'surface-data',                     &
                                              TRIM( run_description_header ), 5551 )
              time_average_text = ' '
           ELSE
              CALL netcdf_create_global_atts( id_set_surf(av), 'surface-data_av',                  &
                                              TRIM( run_description_header ), 5552 )
              WRITE ( time_average_text,'(F7.1,'' s avg'')' )  averaging_interval_surf
              nc_stat = NF90_PUT_ATT( id_set_surf(av), NF90_GLOBAL, 'time_avg',                    &
                                      TRIM( time_average_text ) )
              CALL netcdf_handle_error( 'surface_data_output_mod', 5553 )
           ENDIF


!
!--        Define time coordinate for surface data.
!--        For parallel output the time dimension has to be limited (ntdim_surf), otherwise the
!--        performance drops significantly.
           CALL netcdf_create_dim( id_set_surf(av), 'time', ntdim_surf(av), id_dim_time_surf(av),  &
                                   5554 )

           CALL netcdf_create_var( id_set_surf(av), (/ id_dim_time_surf(av) /),                    &
                                   'time', NF90_DOUBLE, id_var_time_surf(av),                      &
                                   'seconds since '// TRIM( origin_date_time ),                    &
                                   'time', 5555, 5555, 5555 )

           CALL netcdf_create_att( id_set_surf(av), id_var_time_surf(av), 'standard_name', 'time', &
                                   5556)

           CALL netcdf_create_att( id_set_surf(av), id_var_time_surf(av), 'axis', 'T', 5557)
!
!--        Define spatial dimensions and coordinates:
!--        Define index of surface element
           CALL netcdf_create_dim( id_set_surf(av), 's', surf_out%ns_total, id_dim_s_surf(av),     &
                                   5558 )
           CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), 's', NF90_DOUBLE,     &
                                   id_var_s_surf(av), '1', 'number of surface element', 5559,      &
                                   5559, 5559 )
!
!--        Define x coordinate
           CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), 'xs', NF90_DOUBLE,    &
                                   id_var_xs_surf(av), 'meters',                                   &
                                   'distance to origin in x-direction', 5561, 5561, 5561 )
!
!--         Define y coordinate
           CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), 'ys', NF90_DOUBLE,    &
                                   id_var_ys_surf(av), 'meters',                                   &
                                   'distance to origin in y-direction', 5562, 5562, 5562 )
!
!--        Define z coordinate
           CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), 'zs', NF90_DOUBLE,    &
                                   id_var_zs_surf(av), 'meters', 'height', 5560, 5560, 5560 )
           CALL netcdf_create_att( id_set_surf(av), id_var_zs_surf(av), 'standard_name', 'height', &
                                   5583 )

!
!--        Define angles
           CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), 'azimuth',            &
                                   NF90_DOUBLE, id_var_azimuth_surf(av), 'degree',                 &
                                   'azimuth angle', 5577, 5578, 5579, fill = .TRUE. )
           CALL netcdf_create_att( id_set_surf(av), id_var_azimuth_surf(av), 'standard_name',      &
                                   'surface_azimuth_angle', 5584 )

           CALL netcdf_create_var( id_set_surf(av), (/ id_dim_s_surf(av) /), 'zenith',             &
                                   NF90_DOUBLE, id_var_zenith_surf(av), 'degree', '', 5580, 5581,  &
                                   5582, fill = .TRUE. )
!
!--        Define the variables
           var_list = ';'
           i = 1

           DO WHILE ( dosurf(av,i)(1:1) /= ' ' )

              CALL netcdf_create_var( id_set_surf(av),                                             &
                                      (/ id_dim_s_surf(av), id_dim_time_surf(av) /), dosurf(av,i), &
                                      NF90_REAL4, id_var_dosurf(av,i), dosurf_unit(av,i),          &
                                      dosurf(av,i), 5565, 5565, 5565, .TRUE. )
!
!--           Set no fill for every variable to increase performance.
              nc_stat = NF90_DEF_VAR_FILL( id_set_surf(av), id_var_dosurf(av,i), NF90_NOFILL, 0 )
              CALL netcdf_handle_error( 'surface_data_output_init', 5566 )
!
!--           Set collective io operations for parallel io
              nc_stat = NF90_VAR_PAR_ACCESS( id_set_surf(av), id_var_dosurf(av,i), NF90_COLLECTIVE )
              CALL netcdf_handle_error( 'surface_data_output_init', 5567 )
              var_list = TRIM( var_list ) // TRIM( dosurf(av,i) ) // ';'

              i = i + 1

           ENDDO
!
!--        Write the list of variables as global attribute (this is used by restart runs and by
!--        combine_plot_fields)
           nc_stat = NF90_PUT_ATT( id_set_surf(av), NF90_GLOBAL, 'VAR_LIST', var_list )
           CALL netcdf_handle_error( 'surface_data_output_init', 5568 )

!
!--        Set general no fill, otherwise the performance drops significantly for parallel output.
           nc_stat = NF90_SET_FILL( id_set_surf(av), NF90_NOFILL, oldmode )
           CALL netcdf_handle_error( 'surface_data_output_init', 5569 )

!
!--        Leave netCDF define mode
           nc_stat = NF90_ENDDEF( id_set_surf(av) )
           CALL netcdf_handle_error( 'surface_data_output_init', 5570 )

!
!--        These data are only written by PE0 for parallel output to increase the performance.
           IF ( myid == 0 )  THEN
!
!--           Write data for surface indices
              ALLOCATE( netcdf_data_1d(1:surf_out%ns_total) )

              DO  i = 1, surf_out%ns_total
                 netcdf_data_1d(i) = i
              ENDDO

              nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_s_surf(av), netcdf_data_1d,          &
                                      start = (/ 1 /), count = (/ surf_out%ns_total /) )
              CALL netcdf_handle_error( 'surface_data_output_init', 5571 )

              DEALLOCATE( netcdf_data_1d )

           ENDIF

!
!--        Write surface positions to file
           nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_xs_surf(av),                            &
                                   surf_out%xs, start = (/ surf_out%s(1) /),                       &
                                   count = (/ surf_out%ns /) )
           CALL netcdf_handle_error( 'surface_data_output_init', 5572 )

           nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_ys_surf(av), surf_out%ys,               &
                                   start = (/ surf_out%s(1) /), count = (/ surf_out%ns /) )
           CALL netcdf_handle_error( 'surface_data_output_init', 5573 )

           nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_zs_surf(av), surf_out%zs,               &
                                   start = (/ surf_out%s(1) /), count = (/ surf_out%ns /) )
           CALL netcdf_handle_error( 'surface_data_output_init', 5574 )

           nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_azimuth_surf(av), surf_out%azimuth,     &
                                   start = (/ surf_out%s(1) /), count = (/ surf_out%ns /) )
           CALL netcdf_handle_error( 'surface_data_output_init', 5585 )

           nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_zenith_surf(av), surf_out%zenith,       &
                                   start = (/ surf_out%s(1) /), count = (/ surf_out%ns /) )
           CALL netcdf_handle_error( 'surface_data_output_init', 5586 )

          ENDDO
#endif

     ENDIF

 END SUBROUTINE surface_data_output_init

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for controlling the data output. Surface data is collected from different types of
!> surfaces (default, natural, urban) and different orientation and written to one 1D-output array.
!> Further, NetCDF routines are called to write the surface data in the respective NetCDF files.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output( av )

    USE control_parameters,                                                                        &
        ONLY:  io_blocks,                                                                          &
               io_group,                                                                           &
               time_since_reference_point

#if defined( __parallel )
    USE pegrid,                                                                                    &
        ONLY:  comm2d,                                                                             &
               ierr
#endif


    IMPLICIT NONE

    CHARACTER(LEN=100) ::  trimvar = ' '  !< dummy for single output variable

    INTEGER(iwp) ::  av     !< id indicating average or non-average data output
    INTEGER(iwp) ::  i      !< loop index
    INTEGER(iwp) ::  m      !< running index for surface elements
    INTEGER(iwp) ::  n_out  !< counter variables for surface output

    LOGICAL      ::  found  !< flag if output variable is found

!
!-- Return, if nothing to output
    IF ( dosurf_no(av) == 0 )  RETURN
!
!-- In case of VTK output, check if binary files are open and write coordinates.
    IF ( to_vtk )  THEN

       CALL check_open( 25 + av )

       IF ( .NOT. first_output(av) )  THEN
          DO  i = 0, io_blocks - 1
             IF ( i == io_group )  THEN
                WRITE ( 25 + av )  surf_out%npoints
                WRITE ( 25 + av )  surf_out%npoints_total
                WRITE ( 25 + av )  surf_out%ns
                WRITE ( 25 + av )  surf_out%ns_total
                WRITE ( 25 + av )  surf_out%points
                WRITE ( 25 + av )  surf_out%polygons
             ENDIF
#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             first_output(av) = .TRUE.
          ENDDO
       ENDIF
    ENDIF
!
!-- In case of NetCDF output, check if enough time steps are available in file and update time axis.
    IF ( to_netcdf )  THEN
#if defined( __netcdf4_parallel )
       IF ( dosurf_time_count(av) + 1 > ntdim_surf(av) )  THEN
          WRITE ( message_string, * ) 'output of surface data is not given at t=',                 &
                                      time_since_reference_point, 's because the maximum ',        &
                                      'number of output time levels is exceeded'
          CALL message( 'surface_data_output', 'SDO0001', 0, 1, 0, 6, 0 )

          RETURN

       ENDIF
!
!--    Update the netCDF time axis
!--    In case of parallel output, this is only done by PE0 to increase the performance.
       dosurf_time_count(av) = dosurf_time_count(av) + 1
       IF ( myid == 0 )  THEN
          nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_time_surf(av),                           &
                                  (/ time_since_reference_point /),                                &
                                  start = (/ dosurf_time_count(av) /), count = (/ 1 /) )
          CALL netcdf_handle_error( 'surface_data_output', 6666 )
       ENDIF
#endif
    ENDIF

!
!-- Cycle through output quantities and write them to file.
    n_out = 0
    DO  WHILE ( dosurf(av,n_out+1)(1:1) /= ' ' )

       n_out   = n_out + 1
       trimvar = TRIM( dosurf(av,n_out) )
!
!--    Set the output array to the _FillValue in case it is not defined for each type of surface.
       surf_out%var_out = surf_out%fillvalue
       SELECT CASE ( trimvar )

          CASE ( 'us' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%us, surf_lsm%us, surf_usm%us )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'ts' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%ts, surf_lsm%ts, surf_usm%ts )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'qs' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%qs, surf_lsm%qs, surf_usm%qs )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'ss' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%ss, surf_lsm%ss, surf_usm%ss )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'qcs' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%qcs, surf_lsm%qcs, surf_usm%qcs )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'ncs' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%ncs, surf_lsm%ncs, surf_usm%ncs )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'qis' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%qis, surf_lsm%qis, surf_usm%qis )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'nis' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%nis, surf_lsm%nis, surf_usm%nis )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'qrs' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%qrs, surf_lsm%qrs, surf_usm%qrs )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'nrs' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%nrs, surf_lsm%nrs, surf_usm%nrs )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'ol' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%ol, surf_lsm%ol, surf_usm%ol )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'z0' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%z0, surf_lsm%z0, surf_usm%z0 )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'z0h' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%z0h, surf_lsm%z0h, surf_usm%z0h )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'z0q' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%z0q, surf_lsm%z0q, surf_usm%z0q )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'theta1' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%pt1, surf_lsm%pt1, surf_usm%pt1 )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'qv1' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%qv1, surf_lsm%qv1, surf_usm%qv1 )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'thetav1' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%vpt1, surf_lsm%vpt1, surf_usm%vpt1 )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'usws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%usws, surf_lsm%usws, surf_usm%usws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'vsws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%vsws, surf_lsm%vsws, surf_usm%vsws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'shf' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%shf, surf_lsm%shf, surf_usm%shf )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp
             ENDIF

          CASE ( 'qsws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%qsws, surf_lsm%qsws, surf_usm%qsws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'ssws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%ssws, surf_lsm%ssws, surf_usm%ssws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'qcsws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%qcsws, surf_lsm%qcsws, surf_usm%qcsws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'ncsws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%ncsws, surf_lsm%ncsws, surf_usm%ncsws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF


          CASE ( 'qisws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%qisws, surf_lsm%qisws, surf_usm%qisws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'nisws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%nisws, surf_lsm%nisws, surf_usm%nisws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'qrsws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%qrsws, surf_lsm%qrsws, surf_usm%qrsws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'nrsws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%nrsws, surf_lsm%nrsws, surf_usm%nrsws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'sasws' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%sasws, surf_lsm%sasws, surf_usm%sasws )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'q_surface' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%q_surface, surf_lsm%q_surface,          &
                                                  surf_usm%q_surface )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'theta_surface' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%pt_surface, surf_lsm%pt_surface,        &
                                                  surf_usm%pt_surface )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'thetav_surface' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%vpt_surface, surf_lsm%vpt_surface,      &
                                                  surf_usm%vpt_surface )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'theta_10cm' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%pt_10cm, surf_lsm%pt_10cm, surf_usm%pt_10cm )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_net' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_net, surf_lsm%rad_net,              &
                                                  surf_usm%rad_net )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_lw_in' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_lw_in, surf_lsm%rad_lw_in,          &
                                                  surf_usm%rad_lw_in )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_lw_out' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_lw_out, surf_lsm%rad_lw_out,        &
                                                  surf_usm%rad_lw_out )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_sw_in' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_sw_in, surf_lsm%rad_sw_in,          &
                                                  surf_usm%rad_sw_in )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_sw_out' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_sw_out, surf_lsm%rad_sw_out,        &
                                                  surf_usm%rad_sw_out )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'ghf' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
!
!--             Sum up ground / wall heat flux. Note, for urban surfaces the wall heat flux is
!--             aggregated from the different green, window and wall tiles.
                DO  m = 1, surf_usm%ns
                   surf_usm%ghf(m) = surf_usm%frac(m,ind_veg_wall)  * surf_usm%wghf_eb(m) +        &
                                     surf_usm%frac(m,ind_pav_green) * surf_usm%wghf_eb_green(m) +  &
                                     surf_usm%frac(m,ind_wat_win)   * surf_usm%wghf_eb_window(m)
                ENDDO

                CALL surface_data_output_collect( surf_def%ghf, surf_lsm%ghf, surf_usm%ghf )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'r_a' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%r_a, surf_lsm%r_a, surf_usm%r_a )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'r_soil' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%r_soil, surf_lsm%r_soil, surf_usm%r_soil )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'r_canopy' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%r_canopy, surf_lsm%r_canopy, surf_usm%r_canopy )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'r_s' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%r_s, surf_lsm%r_s, surf_usm%r_s )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_sw_dir' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_sw_dir, surf_lsm%rad_sw_dir,        &
                                                  surf_usm%rad_sw_dir )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_sw_dif' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_sw_dif, surf_lsm%rad_sw_dif,        &
                                                  surf_usm%rad_sw_dif )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_sw_ref' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_sw_ref, surf_lsm%rad_sw_ref,        &
                                                  surf_usm%rad_sw_ref )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_sw_res' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_sw_res, surf_lsm%rad_sw_res,        &
                                                  surf_usm%rad_sw_res )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_lw_dif' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_lw_dif, surf_lsm%rad_lw_dif,        &
                                                  surf_usm%rad_lw_dif )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_lw_ref' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_lw_ref, surf_lsm%rad_lw_ref,        &
                                                  surf_usm%rad_lw_ref )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'rad_lw_res' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%rad_lw_res, surf_lsm%rad_lw_res,        &
                                                  surf_usm%rad_lw_res )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF

          CASE ( 'uvw1' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%uvw_abs, surf_lsm%uvw_abs,              &
                                                  surf_usm%uvw_abs )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF
!
!--       Waste heat from indoor model
          CASE ( 'waste_heat' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%waste_heat, surf_lsm%waste_heat,        &
                                                  surf_usm%waste_heat )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF
!
!--       Innermost building wall flux from indoor model
          CASE ( 'im_hf' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%iwghf_eb, surf_lsm%iwghf_eb,            &
                                                  surf_usm%iwghf_eb )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF
!
!--       Surface albedo (tile approach)
          CASE ( 'albedo' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%albedo, surf_lsm%albedo, surf_usm%albedo )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF
!
!--       Surface emissivity (tile approach)
          CASE ( 'emissivity' )
!
!--          Output of instantaneous data
             IF ( av == 0 )  THEN
                CALL surface_data_output_collect( surf_def%emissivity, surf_lsm%emissivity,        &
                                                  surf_usm%emissivity )
             ELSE
!
!--             Output of averaged data
                surf_out%var_out(:) = surf_out%var_av(:,n_out) / REAL( average_count_surf, KIND=wp )
                surf_out%var_av(:,n_out) = 0.0_wp

             ENDIF
!
!--          Add further variables:
!--          'css', 'cssws', 'qsws_liq', 'qsws_soil', 'qsws_veg'

          CASE DEFAULT
!
!--          Try other modules
             found = .FALSE.

             CALL module_interface_data_output_surf( av, trimvar, found )

       END SELECT
!
!--    Write to binary file:
!--    - surfaces%points ( 3, 1-npoints )
!--    - surfaces%polygons ( 5, 1-ns )
!--    - surfaces%var_out ( 1-ns, time )
!--    - Dimension: 1-nsurfaces, 1-npoints - can be ordered consecutively
!--    - Distinguish between average and non-average data
       IF ( to_vtk )  THEN
          DO  i = 0, io_blocks - 1
             IF ( i == io_group )  THEN
                WRITE ( 25 + av )  LEN_TRIM( 'time' )
                WRITE ( 25 + av )  'time'
                WRITE ( 25 + av )  time_since_reference_point
                WRITE ( 25 + av )  LEN_TRIM( trimvar )
                WRITE ( 25 + av )  TRIM( trimvar )
                WRITE ( 25 + av )  surf_out%var_out
             ENDIF
#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
          ENDDO
       ENDIF

       IF ( to_netcdf )  THEN
#if defined( __netcdf4_parallel )
!
!--       Write output array to file
          nc_stat = NF90_PUT_VAR( id_set_surf(av), id_var_dosurf(av,n_out), surf_out%var_out,      &
                                  start = (/ surf_out%s(1), dosurf_time_count(av) /),              &
                                  count = (/ surf_out%ns, 1 /) )
          CALL netcdf_handle_error( 'surface_data_output', 6667 )
#endif
       ENDIF

    ENDDO

!
!-- If averaged output was written to NetCDF file, set the counter to zero
    IF ( av == 1 )  average_count_surf = 0

 END SUBROUTINE surface_data_output

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for controlling the data averaging.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_averaging

    IMPLICIT NONE

    CHARACTER(LEN=100) ::  trimvar  !< dummy variable for current output variable

    INTEGER(iwp) ::  m      !< running index for surface elements
    INTEGER(iwp) ::  n_out  !< counter variables for surface output

    n_out = 0
    DO  WHILE ( dosurf(1,n_out+1)(1:1) /= ' ' )

       n_out   = n_out + 1
       trimvar = TRIM( dosurf(1,n_out) )

       SELECT CASE ( trimvar )

          CASE ( 'us' )
             CALL surface_data_output_sum_up( surf_def%us, surf_lsm%us, surf_usm%us, n_out )

          CASE ( 'ts' )
             CALL surface_data_output_sum_up( surf_def%ts, surf_lsm%ts, surf_usm%ts, n_out )

          CASE ( 'qs' )
             CALL surface_data_output_sum_up( surf_def%qs, surf_lsm%qs, surf_usm%qs, n_out )

          CASE ( 'ss' )
             CALL surface_data_output_sum_up( surf_def%ss, surf_lsm%ss, surf_usm%ss, n_out )

          CASE ( 'qcs' )
             CALL surface_data_output_sum_up( surf_def%qcs, surf_lsm%qcs, surf_usm%qcs, n_out )

          CASE ( 'ncs' )
             CALL surface_data_output_sum_up( surf_def%ncs, surf_lsm%ncs, surf_usm%ncs, n_out )

          CASE ( 'qis' )
             CALL surface_data_output_sum_up( surf_def%qis, surf_lsm%qis, surf_usm%qis, n_out )

          CASE ( 'nis' )
             CALL surface_data_output_sum_up( surf_def%nis, surf_lsm%nis, surf_usm%nis, n_out )

          CASE ( 'qrs' )
             CALL surface_data_output_sum_up( surf_def%qrs, surf_lsm%qrs, surf_usm%qrs, n_out )

          CASE ( 'nrs' )
             CALL surface_data_output_sum_up( surf_def%nrs, surf_lsm%nrs, surf_usm%nrs, n_out )

          CASE ( 'ol' )
             CALL surface_data_output_sum_up( surf_def%ol, surf_lsm%ol, surf_usm%ol, n_out )

          CASE ( 'z0' )
             CALL surface_data_output_sum_up( surf_def%z0, surf_lsm%z0, surf_usm%z0, n_out )

          CASE ( 'z0h' )
             CALL surface_data_output_sum_up( surf_def%z0h, surf_lsm%z0h, surf_usm%z0h, n_out )

          CASE ( 'z0q' )
             CALL surface_data_output_sum_up( surf_def%z0q, surf_lsm%z0q, surf_usm%z0q, n_out )

          CASE ( 'theta1' )
             CALL surface_data_output_sum_up( surf_def%pt1, surf_lsm%pt1, surf_usm%pt1, n_out )

          CASE ( 'qv1' )
             CALL surface_data_output_sum_up( surf_def%qv1, surf_lsm%qv1, surf_usm%qv1, n_out )

          CASE ( 'thetav1' )
             CALL surface_data_output_sum_up( surf_def%vpt1, surf_lsm%vpt1, surf_usm%vpt1, n_out )

          CASE ( 'usws' )
             CALL surface_data_output_sum_up( surf_def%usws, surf_lsm%usws, surf_usm%usws, n_out,  &
                                              momentumflux_output_conversion )

          CASE ( 'vsws' )
             CALL surface_data_output_sum_up( surf_def%vsws, surf_lsm%vsws, surf_usm%vsws, n_out,  &
                                              momentumflux_output_conversion )

          CASE ( 'shf' )
             CALL surface_data_output_sum_up( surf_def%shf, surf_lsm%shf, surf_usm%shf, n_out,     &
                                              heatflux_output_conversion )

          CASE ( 'qsws' )
             CALL surface_data_output_sum_up( surf_def%qsws, surf_lsm%qsws, surf_usm%qsws, n_out,  &
                                              waterflux_output_conversion )

          CASE ( 'ssws' )
             CALL surface_data_output_sum_up( surf_def%ssws, surf_lsm%ssws, surf_usm%ssws, n_out,  &
                                              scalarflux_output_conversion )

          CASE ( 'qcsws' )
             CALL surface_data_output_sum_up( surf_def%qcsws, surf_lsm%qcsws, surf_usm%qcsws, n_out )

          CASE ( 'ncsws' )
             CALL surface_data_output_sum_up( surf_def%ncsws, surf_lsm%ncsws, surf_usm%ncsws, n_out )

          CASE ( 'qisws' )
             CALL surface_data_output_sum_up( surf_def%qisws, surf_lsm%qisws, surf_usm%qisws, n_out )

          CASE ( 'nisws' )
             CALL surface_data_output_sum_up( surf_def%nisws, surf_lsm%nisws, surf_usm%nisws, n_out )

          CASE ( 'qrsws' )
             CALL surface_data_output_sum_up( surf_def%qrsws, surf_lsm%qrsws, surf_usm%qrsws, n_out )

          CASE ( 'nrsws' )
             CALL surface_data_output_sum_up( surf_def%nrsws, surf_lsm%nrsws, surf_usm%nrsws, n_out )

          CASE ( 'sasws' )
             CALL surface_data_output_sum_up( surf_def%sasws, surf_lsm%sasws, surf_usm%sasws, n_out )

          CASE ( 'q_surface' )
             CALL surface_data_output_sum_up( surf_def%q_surface, surf_lsm%q_surface,              &
                                              surf_usm%q_surface, n_out )

          CASE ( 'theta_surface' )
             CALL surface_data_output_sum_up( surf_def%pt_surface, surf_lsm%pt_surface,            &
                                              surf_usm%pt_surface, n_out )

          CASE ( 'thetav_surface' )
             CALL surface_data_output_sum_up( surf_def%vpt_surface, surf_lsm%vpt_surface,          &
                                              surf_usm%vpt_surface, n_out )

          CASE ( 'theta_10cm' )
             CALL surface_data_output_sum_up( surf_def%pt_10cm, surf_lsm%pt_10cm,                  &
                                              surf_usm%pt_10cm , n_out )

          CASE ( 'rad_net' )
             CALL surface_data_output_sum_up( surf_def%rad_net, surf_lsm%rad_net,                  &
                                              surf_usm%rad_net, n_out )

          CASE ( 'rad_lw_in' )
             CALL surface_data_output_sum_up( surf_def%rad_lw_in, surf_lsm%rad_lw_in,              &
                                              surf_usm%rad_lw_in, n_out )

          CASE ( 'rad_lw_out' )
             CALL surface_data_output_sum_up( surf_def%rad_lw_out, surf_lsm%rad_lw_out,            &
                                              surf_usm%rad_lw_out, n_out )

          CASE ( 'rad_sw_in' )
             CALL surface_data_output_sum_up( surf_def%rad_sw_in, surf_lsm%rad_sw_in,              &
                                              surf_usm%rad_sw_in, n_out )

          CASE ( 'rad_sw_out' )
             CALL surface_data_output_sum_up( surf_def%rad_sw_out, surf_lsm%rad_sw_out,            &
                                              surf_usm%rad_sw_out, n_out )

          CASE ( 'ghf' )
!
!--          Sum up ground / wall heat flux. Note, for urban surfaces the wall heat flux is
!--          aggregated from the different green, window and wall tiles.
             DO  m = 1, surf_usm%ns
                surf_usm%ghf(m) = surf_usm%frac(m,ind_veg_wall)  * surf_usm%wghf_eb(m) +           &
                                  surf_usm%frac(m,ind_pav_green) * surf_usm%wghf_eb_green(m) +     &
                                  surf_usm%frac(m,ind_wat_win)   * surf_usm%wghf_eb_window(m)
             ENDDO

             CALL surface_data_output_sum_up( surf_def%ghf, surf_lsm%ghf, surf_usm%ghf, n_out )


          CASE ( 'r_a' )
             CALL surface_data_output_sum_up( surf_def%r_a, surf_lsm%r_a, surf_usm%r_a, n_out )

          CASE ( 'r_soil' )
             CALL surface_data_output_sum_up( surf_def%r_soil, surf_lsm%r_soil, surf_usm%r_soil,   &
                                              n_out )

          CASE ( 'r_canopy' )
             CALL surface_data_output_sum_up( surf_def%r_canopy, surf_lsm%r_canopy,                &
                                              surf_usm%r_canopy, n_out )

          CASE ( 'r_s' )
             CALL surface_data_output_sum_up( surf_def%r_s, surf_lsm%r_s, surf_usm%r_s, n_out )


          CASE ( 'rad_sw_dir' )
             CALL surface_data_output_sum_up( surf_def%rad_sw_dir, surf_lsm%rad_sw_dir,            &
                                              surf_usm%rad_sw_dir, n_out )
          CASE ( 'rad_sw_dif' )
             CALL surface_data_output_sum_up( surf_def%rad_sw_dif, surf_lsm%rad_sw_dif,            &
                                              surf_usm%rad_sw_dif, n_out )

          CASE ( 'rad_sw_ref' )
             CALL surface_data_output_sum_up( surf_def%rad_sw_ref, surf_lsm%rad_sw_ref,            &
                                              surf_usm%rad_sw_ref, n_out )

          CASE ( 'rad_sw_res' )
             CALL surface_data_output_sum_up( surf_def%rad_sw_res, surf_lsm%rad_sw_res,            &
                                              surf_usm%rad_sw_res, n_out )

          CASE ( 'rad_lw_dif' )
             CALL surface_data_output_sum_up( surf_def%rad_lw_dif, surf_lsm%rad_lw_dif,            &
                                              surf_usm%rad_lw_dif, n_out )

          CASE ( 'rad_lw_ref' )
             CALL surface_data_output_sum_up( surf_def%rad_lw_ref, surf_lsm%rad_lw_ref,            &
                                              surf_usm%rad_lw_ref, n_out )

          CASE ( 'rad_lw_res' )
             CALL surface_data_output_sum_up( surf_def%rad_lw_res, surf_lsm%rad_lw_res,            &
                                              surf_usm%rad_lw_res, n_out )

          CASE ( 'uvw1' )
             CALL surface_data_output_sum_up( surf_def%uvw_abs, surf_lsm%uvw_abs,                  &
                                              surf_usm%uvw_abs, n_out )

          CASE ( 'waste_heat' )
             CALL surface_data_output_sum_up( surf_def%waste_heat, surf_lsm%waste_heat,            &
                                              surf_usm%waste_heat, n_out )

          CASE ( 'im_hf' )
             CALL surface_data_output_sum_up( surf_def%iwghf_eb, surf_lsm%iwghf_eb,                &
                                              surf_usm%iwghf_eb, n_out )

          CASE ( 'albedo' )
             CALL surface_data_output_sum_up( surf_def%albedo, surf_lsm%albedo,                    &
                                              surf_usm%albedo, n_out )


          CASE ( 'emissivity' )
             CALL surface_data_output_sum_up( surf_def%emissivity, surf_lsm%emissivity,            &
                                              surf_usm%emissivity, n_out )

          CASE DEFAULT
!
!--          Process variables from other modules
             CALL module_interface_surface_data_averaging( trimvar, n_out )

       END SELECT
    ENDDO


 END SUBROUTINE surface_data_output_averaging

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum-up the surface data for average output variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_sum_up_1d( var_def, var_lsm, var_usm, n_out, fac )

    IMPLICIT NONE

    INTEGER(iwp) ::  k       !< height index of surface element
    INTEGER(iwp) ::  m       !< running index for surface elements
    INTEGER(iwp) ::  n_out   !< index for output variable
    INTEGER(iwp) ::  n_surf  !< running index for surface elements

    REAL(wp), DIMENSION(:), OPTIONAL                ::  fac                !< passed output conversion factor for heatflux output
    REAL(wp), DIMENSION(nzb:nzt+1)                  ::  conversion_factor  !< effective array for output conversion factor

    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def  !< output variable for default-type surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm  !< output variable for natural-type surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm  !< output variable for urban-type surfaces

!
!-- Set conversion factor to one if not present
    IF ( .NOT. PRESENT( fac ) )  THEN
       conversion_factor = 1.0_wp
    ELSE
       conversion_factor = fac
    ENDIF
!
!-- Set counter variable to zero before the variable is written to the output array.
    n_surf = 0
!
!-- Write the horizontal surfaces.
!-- Before each variable is written to the output data structure, first check if the variable
!-- for the respective surface type is defined. If a variable is not defined, skip the block and
!-- increment the counter variable by the number of surface elements of this type. Usually this is
!-- zero, however, there might be the situation that e.g. urban surfaces are defined but the
!-- respective variable is not allocated for this surface type. To write the data on the exact
!-- position, increment the counter.
    IF ( ALLOCATED( var_def ) )  THEN
       DO  m = 1, surf_def%ns
          n_surf                        = n_surf + 1
          k                             = surf_def%k(m)
          surf_out%var_av(n_surf,n_out) = surf_out%var_av(n_surf,n_out) + var_def(m) *             &
                                          conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_def%ns
    ENDIF
    IF ( ALLOCATED( var_lsm ) )  THEN
       DO  m = 1, surf_lsm%ns
          n_surf                        = n_surf + 1
          k                             = surf_lsm%k(m)
          surf_out%var_av(n_surf,n_out) = surf_out%var_av(n_surf,n_out) + var_lsm(m) *             &
                                          conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_lsm%ns
    ENDIF
    IF ( ALLOCATED( var_usm ) )  THEN
       DO  m = 1, surf_usm%ns
          n_surf                        = n_surf + 1
          k                             = surf_usm%k(m)
          surf_out%var_av(n_surf,n_out) = surf_out%var_av(n_surf,n_out) + var_usm(m) *             &
                                          conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_usm%ns
    ENDIF

 END SUBROUTINE surface_data_output_sum_up_1d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum-up the surface data for average output variables for properties which are defined using tile
!> approach.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_sum_up_2d( var_def, var_lsm, var_usm, n_out, fac )

    IMPLICIT NONE

    INTEGER(iwp) ::  k       !< height index of surface element
    INTEGER(iwp) ::  m       !< running index for surface elements
    INTEGER(iwp) ::  n_out   !< index for output variable
    INTEGER(iwp) ::  n_surf  !< running index for surface elements

    REAL(wp), DIMENSION(:), OPTIONAL                  ::  fac                !< passed output conversion factor for heatflux output
    REAL(wp), DIMENSION(nzb:nzt+1)                    ::  conversion_factor  !< effective array for output conversion factor

    REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) ::  var_def  !< output variable for default-type surfaces
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) ::  var_lsm  !< output variable for natural-type surfaces
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) ::  var_usm  !< output variable for urban-type surfaces
!
!-- Set conversion factor to one if not present
    IF ( .NOT. PRESENT( fac ) )  THEN
       conversion_factor = 1.0_wp
    ELSE
       conversion_factor = fac
    ENDIF
!
!-- Set counter variable to zero before the variable is written to the output array.
    n_surf = 0

!
!-- Write the horizontal surfaces.
!-- Before each variable is written to the output data structure, first check if the variable
!-- for the respective surface type is defined. If a variable is not defined, skip the block and
!-- increment the counter variable by the number of surface elements of this type. Usually this is
!-- zero, however, there might be the situation that e.g. urban surfaces are defined but the
!-- respective variable is not allocated for this surface type. To write the data on the exact
!-- position, increment the counter.
    IF ( ALLOCATED( var_def ) )  THEN
       DO  m = 1, surf_def%ns
          n_surf = n_surf + 1
          k      = surf_def%k(m)
          surf_out%var_av(n_surf,n_out) = surf_out%var_av(n_surf,n_out)                            &
                                          + SUM ( surf_def%frac(m,:) * var_def(m,:) )              &
                                          * conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_def%ns
    ENDIF
    IF ( ALLOCATED( var_lsm ) )  THEN
       DO  m = 1, surf_lsm%ns
          n_surf = n_surf + 1
          k      = surf_lsm%k(m)
          surf_out%var_av(n_surf,n_out) = surf_out%var_av(n_surf,n_out)                            &
                                          + SUM ( surf_lsm%frac(m,:) * var_lsm(m,:) )              &
                                          * conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_lsm%ns
    ENDIF
    IF ( ALLOCATED( var_usm ) )  THEN
       DO  m = 1, surf_usm%ns
          n_surf = n_surf + 1
          k      = surf_usm%k(m)
          surf_out%var_av(n_surf,n_out) = surf_out%var_av(n_surf,n_out)                            &
                                          + SUM ( surf_usm%frac(m,:) * var_usm(m,:) )              &
                                          * conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_usm%ns
    ENDIF

 END SUBROUTINE surface_data_output_sum_up_2d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collect the surface data from different types and different orientation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_collect_1d( var_def, var_lsm, var_usm, fac )

    IMPLICIT NONE

    INTEGER(iwp) ::  k       !< height index of surface element
    INTEGER(iwp) ::  m       !< running index for surface elements
    INTEGER(iwp) ::  n_surf  !< running index for surface elements

    REAL(wp), DIMENSION(:), OPTIONAL                ::  fac                !< passed output conversion factor for heatflux output
    REAL(wp), DIMENSION(nzb:nzt+1)                  ::  conversion_factor  !< effective array for output conversion factor

    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_def  !< output variable for default-type surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_lsm  !< output variable for natural-type surfaces
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(IN) ::  var_usm  !< output variable for urban-type surfaces

!
!-- Set conversion factor to one if not present
    IF ( .NOT. PRESENT( fac ) )  THEN
       conversion_factor = 1.0_wp
    ELSE
       conversion_factor = fac
    ENDIF
!
!-- Set counter variable to zero before the variable is written to the output array.
    n_surf = 0
!
!-- Write the horizontal surfaces.
!-- Before each variable is written to the output data structure, first check if the variable
!-- for the respective surface type is defined. If a variable is not defined, skip the block and
!-- increment the counter variable by the number of surface elements of this type. Usually this is
!-- zero, however, there might be the situation that e.g. urban surfaces are defined but the
!-- respective variable is not allocated for this surface type. To write the data on the exact
!-- position, increment the counter.
    IF ( ALLOCATED( var_def ) )  THEN
       DO  m = 1, surf_def%ns
          n_surf = n_surf + 1
          k      = surf_def%k(m)
          surf_out%var_out(n_surf) = var_def(m) * conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_def%ns
    ENDIF
    IF ( ALLOCATED( var_lsm ) )  THEN
       DO  m = 1, surf_lsm%ns
          n_surf = n_surf + 1
          k      = surf_lsm%k(m)
          surf_out%var_out(n_surf) = var_lsm(m) * conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_lsm%ns
    ENDIF
    IF ( ALLOCATED( var_usm ) )  THEN
       DO  m = 1, surf_usm%ns
          n_surf = n_surf + 1
          k      = surf_usm%k(m)
          surf_out%var_out(n_surf) = var_usm(m) * conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_usm%ns
    ENDIF

 END SUBROUTINE surface_data_output_collect_1d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Collect the surface data from different types and different orientation for properties which are
!> defined using tile approach.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_collect_2d( var_def, var_lsm, var_usm, fac )

    IMPLICIT NONE

    INTEGER(iwp) ::  k       !< height index of surface element
    INTEGER(iwp) ::  m       !< running index for surface elements
    INTEGER(iwp) ::  n_surf  !< running index for surface elements

    REAL(wp), DIMENSION(:), OPTIONAL                  ::  fac                !< passed output conversion factor for heatflux output
    REAL(wp), DIMENSION(nzb:nzt+1)                    ::  conversion_factor  !< effective array for output conversion factor

    REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) ::  var_def  !< output variable for default-type surfaces
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) ::  var_lsm  !< output variable for natural-type surfaces
    REAL(wp), DIMENSION(:,:), ALLOCATABLE, INTENT(IN) ::  var_usm  !< output variable for urban-type surfaces
!
!-- Set conversion factor to one if not present
    IF ( .NOT. PRESENT( fac ) )  THEN
       conversion_factor = 1.0_wp
    ELSE
       conversion_factor = fac
    ENDIF
!
!-- Set counter variable to zero before the variable is written to the output array.
    n_surf = 0
!
!-- Write the horizontal surfaces.
!-- Before each variable is written to the output data structure, first check if the variable
!-- for the respective surface type is defined. If a variable is not defined, skip the block and
!-- increment the counter variable by the number of surface elements of this type. Usually this is
!-- zero, however, there might be the situation that e.g. urban surfaces are defined but the
!-- respective variable is not allocated for this surface type. To write the data on the exact
!-- position, increment the counter.
    IF ( ALLOCATED( var_def ) )  THEN
       DO  m = 1, surf_def%ns
          n_surf = n_surf + 1
          k      = surf_def%k(m)
          surf_out%var_out(n_surf) = SUM( surf_def%frac(m,:) * var_def(m,:) ) * conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_def%ns
    ENDIF

    IF ( ALLOCATED( var_lsm ) )  THEN
       DO  m = 1, surf_lsm%ns
          n_surf = n_surf + 1
          k      = surf_lsm%k(m)
          surf_out%var_out(n_surf) = SUM( surf_lsm%frac(m,:) * var_lsm(m,:) ) * conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_lsm%ns
    ENDIF

    IF ( ALLOCATED( var_usm ) )  THEN
       DO  m = 1, surf_usm%ns
          n_surf = n_surf + 1
          k      = surf_usm%k(m)
          surf_out%var_out(n_surf) = SUM( surf_usm%frac(m,:) * var_usm(m,:) ) * conversion_factor(k)
       ENDDO
    ELSE
       n_surf = n_surf + surf_usm%ns
    ENDIF

 END SUBROUTINE surface_data_output_collect_2d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for output of surface parameters
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_parin

    IMPLICIT NONE

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /surface_data_output_parameters/  averaging_interval_surf,                            &
                                               data_output_surf,                                   &
                                               dt_dosurf,                                          &
                                               dt_dosurf_av,                                       &
                                               skip_time_dosurf,                                   &
                                               skip_time_dosurf_av,                                &
                                               switch_off_module,                                  &
                                               to_netcdf,                                          &
                                               to_vtk

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, surface_data_output_parameters, IOSTAT=io_status )
!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    surface_data_output_parameters namelist was found and read correctly. Set flag that indicates
!--    that surface data output is switched on.
       IF ( .NOT. switch_off_module )  surface_output = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    surface_data_output_parameters namelist was found but contained errors. Print an error
!--    message including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'surface_data_output_parameters', line )

    ENDIF

 END SUBROUTINE surface_data_output_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check the input parameters for consistency. Further pre-process the given output variables, i.e.
!> separate them into average and non-average output variables and map them onto internal output
!> array.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_check_parameters

    USE control_parameters,                                                                        &
        ONLY:  averaging_interval,                                                                 &
               dt_data_output,                                                                     &
               dt_data_output_av,                                                                  &
               indoor_model,                                                                       &
               message_string

    IMPLICIT NONE

    CHARACTER(LEN=100) ::  trimvar  !< dummy for single output variable
    CHARACTER(LEN=100) ::  unit     !< dummy for unit of output variable

    INTEGER(iwp) ::  av     !< id indicating average or non-average data output
    INTEGER(iwp) ::  ilen   !< string length
    INTEGER(iwp) ::  ivar   !< loop index
    INTEGER(iwp) ::  n_out  !< running index for number of output variables

    LOGICAL ::  is_duplicate  !< true if string has duplicates in list
!
!-- Check if any output file type is selected
    IF ( .NOT. to_vtk  .AND.  .NOT. to_netcdf )  THEN
       WRITE( message_string, * ) 'No output file type selected for surface-data output.&' //      &
                                  'Set at least either "to_vtk" or "to_netcdf" to .TRUE.'
       CALL message( 'surface_data_output_check_parameters', 'SDO0002', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check the average interval
    IF ( averaging_interval_surf == 9999999.9_wp )  THEN
       averaging_interval_surf = averaging_interval
    ENDIF
!
!-- Set the default data-output intervals to dt_data_output or dt_data_output_av if necessary
    IF ( dt_dosurf    == 9999999.9_wp )  dt_dosurf    = dt_data_output
    IF ( dt_dosurf_av == 9999999.9_wp )  dt_dosurf_av = dt_data_output_av

    IF ( averaging_interval_surf > dt_dosurf_av )  THEN
       WRITE( message_string, * )  'averaging_interval_surf = ', averaging_interval_surf,          &
                                   ' must be <= dt_dosurf_av = ', dt_dosurf_av
       CALL message( 'surface_data_output_check_parameters', 'SDO0003', 1, 2, 0, 6, 0 )
    ENDIF

#if ! defined( __netcdf4_parallel )
!
!-- Surface output via NetCDF requires parallel NetCDF
    IF ( to_netcdf )  THEN
       message_string = 'to_netcdf = .True. requires parallel netCDF'
       CALL message( 'surface_data_output_check_parameters', 'SDO0004', 1, 2, 0, 6, 0 )
    ENDIF
#endif
!
!-- In case of parallel NetCDF output the output timestep must not be zero. This is because the
!-- number of requiered output timesteps is pre-calculated, which is not possible with zero output
!-- timestep.
    IF ( netcdf_data_format > 4 )  THEN
       IF ( dt_dosurf == 0.0_wp )  THEN
          message_string = 'dt_dosurf = 0.0 while using a variable timestep and parallel ' //      &
                           'netCDF4 is not allowed.'
          CALL message( 'surface_data_output_check_parameters', 'SDO0005', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( dt_dosurf_av == 0.0_wp )  THEN
          message_string = 'dt_dosurf_av = 0.0 while using a variable timestep and parallel ' //   &
                           'netCDF4 is not allowed.'
          CALL message( 'surface_data_output_check_parameters', 'SDO0005', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- Count number of output variables and separate output strings for average and non-average output
!-- variables.
    n_out = 0
    DO WHILE ( data_output_surf(n_out+1)(1:1) /= ' ' )

       n_out = n_out + 1
       ilen = LEN_TRIM( data_output_surf(n_out) )
       trimvar = TRIM( data_output_surf(n_out) )

!
!--    Check for duplications
       is_duplicate = .FALSE.
       DO ivar = 1, n_out-1
          IF ( trimvar == TRIM( data_output_surf(ivar) ) ) THEN
             message_string = 'The variable ' // TRIM( trimvar ) //                                &
                              ' is defined more than once. Duplications are removed.'
             CALL message( 'surface_data_output_check_parameters', 'SDO0006', 0, 1, 0, 6, 0 )
             is_duplicate = .TRUE.
             CYCLE
          ENDIF
       ENDDO
       IF ( is_duplicate ) CYCLE
!
!--    Check for data averaging
       av = 0
       IF ( ilen > 3 )  THEN
          IF ( data_output_surf(n_out)(ilen-2:ilen) == '_av' )  THEN
             trimvar = data_output_surf(n_out)(1:ilen-3)
             av      = 1
          ENDIF
       ENDIF

       dosurf_no(av) = dosurf_no(av) + 1
       dosurf(av,dosurf_no(av)) = TRIM( trimvar )

!
!--    Check if all output variables are known and assign a unit.
!--    Note, in case of shf and qsws the unit depends on the setting of flux_output_mode.
!--    Further, please note, for the passive scalar as well as for the momentum fluxes
!--    the output is independent on the flux_output_mode, because the density cancels out.
       unit = 'not set'
       SELECT CASE ( TRIM( trimvar ) )

          CASE ( 'css', 'cssws', 'qsws_liq', 'qsws_soil', 'qsws_veg' )
             message_string = '"' // TRIM( trimvar ) // '" is not yet implemented in the ' //      &
                              'surface output'
             CALL message( 'surface_data_output_check_parameters', 'SDO0007', 1, 2, 0, 6, 0 )

          CASE ( 'us', 'uvw1' )
             unit = 'm/s'

          CASE ( 'ss', 'qcs', 'ncs', 'qis', 'nis', 'qrs', 'nrs' )
             unit = '1'

          CASE ( 'z0', 'z0h', 'z0q', 'ol' )
             unit = 'm'

          CASE ( 'ts', 'theta1', 'thetav1', 'theta_surface', 'thetav_surface' )
             unit = 'K'

          CASE ( 'usws', 'vsws' )
             unit = 'm2/s2'

          CASE ( 'qcsws', 'ncsws', 'qisws', 'nisws', 'qrsws', 'nrsws', 'sasws' )

          CASE ( 'shf' )
             IF ( TRIM( flux_output_mode ) == 'kinematic' )  THEN
                unit = 'K m/s'
             ELSE
                unit = 'W/m2'
             ENDIF

          CASE ( 'qsws' )
             IF ( TRIM( flux_output_mode ) == 'kinematic' )  THEN
                unit = 'kg/kg m/s'
             ELSE
                unit = 'W/m2'
             ENDIF

          CASE ( 'ssws' )
             unit = 'kg/m2/s'

          CASE ( 'qs', 'q_surface', 'qv1' )
             unit = 'kg/kg'

          CASE ( 'rad_net' )
             unit = 'W/m2'

          CASE ( 'rad_lw_in', 'rad_lw_out', 'rad_lw_dif', 'rad_lw_ref', 'rad_lw_res' )
             unit = 'W/m2'

          CASE ( 'rad_sw_in', 'rad_sw_out', 'rad_sw_dif', 'rad_sw_ref', 'rad_sw_res', 'rad_sw_dir' )
             unit = 'W/m2'

          CASE ( 'ghf' )
             unit = 'W/m2'

          CASE ( 'r_a', 'r_canopy', 'r_soil', 'r_s' )
             unit = 's/m'

          CASE ( 'waste_heat', 'im_hf' )
             IF ( .NOT. indoor_model )  THEN
                message_string = '"' // TRIM( trimvar ) // '" requires the indoor model'
                CALL message( 'surface_data_output_check_parameters', 'SDO0009', 1, 2, 0, 6, 0 )
             ENDIF

             unit = 'W/m2'

           CASE ( 'theta_10cm' )
             IF ( .NOT. indoor_model )  THEN
                message_string = '"' // TRIM( trimvar ) // '" requires the indoor model'
                CALL message( 'surface_data_output_check_parameters', 'SDO0009', 1, 2, 0, 6, 0 )
             ENDIF

             unit = 'K'

          CASE ( 'albedo', 'emissivity' )
             unit = '1'

          CASE DEFAULT
!
!--          Check for other modules
             CALL module_interface_check_data_output_surf( trimvar, unit, av )

             IF ( unit == 'illegal' )  THEN
                message_string = '"' // TRIM( trimvar ) // '" is not part of the surface output'
                CALL message( 'surface_data_output_check_parameters', 'SDO0010', 1, 2, 0, 6, 0 )
             ENDIF

       END SELECT

       dosurf_unit(av,dosurf_no(av)) = unit

    ENDDO

 END SUBROUTINE surface_data_output_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Last action.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_last_action( av )

   USE control_parameters,                                                                         &
       ONLY:  io_blocks,                                                                           &
              io_group

#if defined( __parallel )
   USE pegrid,                                                                                     &
       ONLY:  comm2d,                                                                              &
              ierr
#endif

   IMPLICIT NONE

   INTEGER(iwp) ::  av  !< id indicating average or non-average data output
   INTEGER(iwp) ::  i   !< loop index

!
!--Return, if nothing to output
   IF ( dosurf_no(av) == 0 )  RETURN
!
!--If output to VTK files is enabled, check if files are open and write an end-of-file statement.
   IF ( to_vtk )  THEN
      CALL check_open( 25 + av )
!
!--   Write time coordinate
      DO  i = 0, io_blocks - 1
         IF ( i == io_group )  THEN
            WRITE ( 25 + av )  LEN_TRIM( 'END' )
            WRITE ( 25 + av )  'END'
         ENDIF
#if defined( __parallel )
         CALL MPI_BARRIER( comm2d, ierr )
#endif
      ENDDO
   ENDIF

 END SUBROUTINE surface_data_output_last_action


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!---------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_rrd_global_ftn( found )


    USE control_parameters,                                                                        &
        ONLY:  length,                                                                             &
               restart_string

    IMPLICIT NONE

    LOGICAL, INTENT(OUT) ::  found  !< flag indicating if variable was found

    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'average_count_surf' )
          READ ( 13 )  average_count_surf
       CASE ( 'time_dosurf' )
          READ ( 13 )  time_dosurf
       CASE ( 'time_dosurf_av' )
          READ ( 13 )  time_dosurf_av

       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE surface_data_output_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_rrd_global_mpi

    CALL rrd_mpi_io( 'average_count_surf', average_count_surf )
    CALL rrd_mpi_io( 'time_dosurf', time_dosurf )
    CALL rrd_mpi_io( 'time_dosurf_av', time_dosurf_av )

 END SUBROUTINE surface_data_output_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_rrd_local_ftn( found )


    USE control_parameters,                                                                        &
        ONLY:  length,                                                                             &
               restart_string

    IMPLICIT NONE

    LOGICAL, INTENT(OUT) ::  found


    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'surf_out%var_av' )
          READ ( 13 )  surf_out%var_av

       CASE DEFAULT

          found = .FALSE.

       END SELECT


 END SUBROUTINE surface_data_output_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_rrd_local_mpi

    IMPLICIT NONE

    CHARACTER(LEN=3) ::  dum !< dummy string to create output-variable name

    INTEGER(iwp) ::  i  !< grid index in x-direction
    INTEGER(iwp) ::  j  !< grid index in y-direction
    INTEGER(iwp) ::  m  !< running index surface elements
    INTEGER(iwp) ::  n  !< counting variable
    INTEGER(iwp) ::  nv !< running index over number of variables

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  end_index           !< end index of surface data at (j,i)
    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE ::  global_end_index    !< end index array for surface data (MPI-IO)
    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE ::  global_start_index  !< start index array for surface data (MPI-IO)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  num_surf            !< number of surface data at (j,i)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  start_index         !< start index of surface data at (j,i)

    LOGICAL ::  array_found  !< flag indicating whether variable is in restart data or not
    LOGICAL ::  ldum         !< dummy variable

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surf_in !< input array in expected restart format

!
!-- Skip restart input in case of cyclic-fill initialization. This case, input of time-
!-- averaged data is useless and can lead to faulty averaging.
    IF ( cyclic_fill_initialization )  RETURN
!
!-- Note, surface data which is written to file is organized in a different way than
!-- the output surface data. The output surface data is a concatenated array of the
!-- different surface types and orientations, while the mpi-io expects surface data that
!-- is consecutive in terms of start- and end-index, i.e. organized along the (j,i)
!-- grid index. Hence, data need to be tranformed back to the output surface data.
    ALLOCATE( end_index(nys:nyn,nxl:nxr)          )
    ALLOCATE( global_end_index(nys:nyn,nxl:nxr)   )
    ALLOCATE( global_start_index(nys:nyn,nxl:nxr) )
    ALLOCATE( num_surf(nys:nyn,nxl:nxr)           )
    ALLOCATE( start_index(nys:nyn,nxl:nxr)        )
    ALLOCATE( surf_in(1:surf_out%ns) )

    CALL rd_mpi_io_check_array( 'surf_out%global_start_index', found = array_found )
    IF ( array_found )  CALL rrd_mpi_io( 'surf_out%global_start_index', global_start_index )

    CALL rd_mpi_io_check_array( 'surf_out%global_end_index', found = array_found )
    IF ( array_found )  CALL rrd_mpi_io( 'surf_out%global_end_index', global_end_index )
!
!-- Check if data input for surface-output variables is required. Note, only invoke routine if
!-- restart data is on file. To check this, use the array_found control flag.
    IF ( array_found )  THEN
       CALL rd_mpi_io_surface_filetypes( start_index, end_index, ldum, global_start_index,         &
                                         global_end_index )
    ENDIF

    DO  nv = 1, dosurf_no(1)
       WRITE( dum, '(I3.3)' )  nv

       CALL rd_mpi_io_check_array( 'surf_out%var_av' // TRIM( dum ), found = array_found )

       IF ( array_found )  THEN

          CALL rrd_mpi_io_surface( 'surf_out%var_av' // TRIM(dum), surf_in )
!
!--       Write temporary input variable back to surface-output data array.
          n = 0
          num_surf = 0
          DO  m = 1, surf_def%ns
             i = surf_def%i(m)
             j = surf_def%j(m)
             n = n + 1
             surf_out%var_av(n,nv) = surf_in(start_index(j,i)+num_surf(j,i))
             num_surf(j,i)         = num_surf(j,i) + 1
          ENDDO
          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m)
             j = surf_lsm%j(m)
             n = n + 1
             surf_out%var_av(n,nv) = surf_in(start_index(j,i)+num_surf(j,i))
             num_surf(j,i)         = num_surf(j,i) + 1
          ENDDO
          DO  m = 1, surf_usm%ns
             i = surf_usm%i(m)
             j = surf_usm%j(m)
             n = n + 1
             surf_out%var_av(n,nv) = surf_in(start_index(j,i)+num_surf(j,i))
             num_surf(j,i)         = num_surf(j,i) + 1
          ENDDO
       ENDIF
    ENDDO

    DEALLOCATE( end_index )
    DEALLOCATE( global_end_index )
    DEALLOCATE( global_start_index )
    DEALLOCATE( num_surf )
    DEALLOCATE( start_index )

 END SUBROUTINE surface_data_output_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_wrd_global

    IMPLICIT NONE

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL wrd_write_string( 'average_count_surf' )
       WRITE ( 14 )  average_count_surf

       CALL wrd_write_string( 'time_dosurf' )
       WRITE ( 14 )  time_dosurf

       CALL wrd_write_string( 'time_dosurf_av' )
       WRITE ( 14 )  time_dosurf_av

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

      CALL wrd_mpi_io( 'average_count_surf', average_count_surf )
      CALL wrd_mpi_io( 'time_dosurf', time_dosurf )
      CALL wrd_mpi_io( 'time_dosurf_av', time_dosurf_av )

    ENDIF

 END SUBROUTINE surface_data_output_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes restart data which individual on each PE
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_data_output_wrd_local

    IMPLICIT NONE

    CHARACTER(LEN=3) ::  dum !< dummy string to create output-variable name

    INTEGER(iwp) ::  i                      !< grid index in x-direction
    INTEGER(iwp) ::  j                      !< grid index in y-direction
    INTEGER(iwp) ::  m                      !< running index surface elements
    INTEGER(iwp) ::  n                      !< counting variable
    INTEGER(iwp) ::  nv                     !< running index over number of variables
    INTEGER(iwp) ::  start_index_aggregated !< sum of start-index at (j,i) over all surface types

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  end_index           !< end index of surface data at (j,i)
    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE ::  global_end_index    !< end index array for surface data (MPI-IO)
    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE ::  global_start_index  !< start index array for surface data (MPI-IO)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  num_surf            !< number of surface data at (j,i)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  start_index         !< start index of surface data at (j,i)

    LOGICAL ::  surface_data_to_write  !< switch for MPI-I/O if PE has surface data to write

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  surf_data  !< surface data in expected restart format


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       IF ( ALLOCATED( surf_out%var_av ) )  THEN
          CALL wrd_write_string( 'surf_out%var_av' )
          WRITE ( 14 )  surf_out%var_av
       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
!
!--    Note, surface data which is written to file is organized in a different way than
!--    the output surface data. The output surface data is a concatenated array of the
!--    different surface types and orientations, while the mpi-io expects surface data that
!--    is consecutive in terms of start- and end-index, i.e. organized along the (j,i)
!--    grid index. Hence, data need to be tranformed before it can be written to file.
       IF ( ALLOCATED( surf_out%var_av ) )  THEN

          ALLOCATE( end_index(nys:nyn,nxl:nxr)          )
          ALLOCATE( global_end_index(nys:nyn,nxl:nxr)   )
          ALLOCATE( global_start_index(nys:nyn,nxl:nxr) )
          ALLOCATE( num_surf(nys:nyn,nxl:nxr)           )
          ALLOCATE( start_index(nys:nyn,nxl:nxr)        )
          ALLOCATE( surf_data(1:surf_out%ns)            )
!
!--       Determine the start and end index at each (j,i)-pair and resort the surface data
          start_index            = 1
          end_index              = 0
          start_index_aggregated = 1
          num_surf = 0
          DO  m = 1, surf_def%ns
             i = surf_def%i(m)
             j = surf_def%j(m)
             num_surf(j,i) = num_surf(j,i) + 1
          ENDDO
          DO  m = 1, surf_lsm%ns
             i = surf_lsm%i(m)
             j = surf_lsm%j(m)
             num_surf(j,i) = num_surf(j,i) + 1
          ENDDO
          DO  m = 1, surf_usm%ns
             i = surf_usm%i(m)
             j = surf_usm%j(m)
             num_surf(j,i) = num_surf(j,i) + 1
          ENDDO

          start_index = 0
          end_index   = 0
          start_index_aggregated = 1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                start_index(j,i)       = start_index_aggregated
                end_index(j,i)         = start_index(j,i) + num_surf(j,i) - 1
                start_index_aggregated = start_index_aggregated + num_surf(j,i)
             ENDDO
          ENDDO

          CALL rd_mpi_io_surface_filetypes( start_index, end_index, surface_data_to_write,         &
                                            global_start_index, global_end_index )

          CALL wrd_mpi_io( 'surf_out%global_start_index', global_start_index )
          CALL wrd_mpi_io( 'surf_out%global_end_index', global_end_index )

          DO  nv = 1, dosurf_no(1)
             n = 0
             num_surf = 0
             DO  m = 1, surf_def%ns
                i = surf_def%i(m)
                j = surf_def%j(m)
                n = n + 1
                surf_data(start_index(j,i)+num_surf(j,i)) = surf_out%var_av(n,nv)
                num_surf(j,i)                             = num_surf(j,i) + 1
             ENDDO
             DO  m = 1, surf_lsm%ns
                i = surf_lsm%i(m)
                j = surf_lsm%j(m)
                n = n + 1
                surf_data(start_index(j,i)+num_surf(j,i)) = surf_out%var_av(n,nv)
                num_surf(j,i)                             = num_surf(j,i) + 1
             ENDDO
             DO  m = 1, surf_usm%ns
                i = surf_usm%i(m)
                j = surf_usm%j(m)
                n = n + 1
                surf_data(start_index(j,i)+num_surf(j,i)) = surf_out%var_av(n,nv)
                num_surf(j,i)                             = num_surf(j,i) + 1
             ENDDO

             WRITE( dum, '(I3.3)' )  nv

             CALL wrd_mpi_io_surface( 'surf_out%var_av' // TRIM( dum ), surf_data )
          ENDDO

          DEALLOCATE( end_index )
          DEALLOCATE( global_end_index )
          DEALLOCATE( global_start_index )
          DEALLOCATE( num_surf )
          DEALLOCATE( start_index )
          DEALLOCATE( surf_data )

       ENDIF

    ENDIF

 END SUBROUTINE surface_data_output_wrd_local


END MODULE surface_data_output_mod
