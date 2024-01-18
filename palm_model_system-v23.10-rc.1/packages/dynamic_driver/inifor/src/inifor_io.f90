!> @file src/inifor_io.f90
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
!> The io module contains the functions needed to read and write netCDF data in
!> INIFOR.
!------------------------------------------------------------------------------!
#if defined ( __netcdf )
 MODULE inifor_io

    USE inifor_control
    USE inifor_defs,                                                           &
        ONLY:  CFG_INIT_PROFILE, CFG_INIT_VOLUME,                              &
               CFG_INIT_SOIL_PROFILE, CFG_INIT_SOIL_VOLUME,                    &
               CFG_FORCING_HETERO, CFG_FORCING_HOMO, CFG_FORCING_NUDGING,      &
               DATE, SNAME, PATH, PI, TO_RADIANS, TO_DEGREES, VERSION,         &
               NC_DEPTH_DIM_IDX, NC_DEPTH_NAME, NC_HHL_NAME, NC_RLAT_NAME,     &
               NC_RLON_NAME, NC_ROTATED_POLE_NAME, NC_POLE_LATITUDE_NAME,      &
               NC_POLE_LONGITUDE_NAME, RHO_L, iwp, wp,                         &
               PIDS_ORIGIN_LON, PIDS_ORIGIN_LAT, PIDS_ORIGIN_Z
    USE inifor_types
    USE inifor_util,                                                           &
        ONLY:  add_hours_to, nearly_equal, reverse, str, real_to_str
    USE netcdf

    IMPLICIT NONE

!------------------------------------------------------------------------------!
! Description:
! ------------
!> get_netcdf_variable() reads the netCDF data and metadate for the netCDF
!> variable 'in_var%name' from the file 'in_file'. The netCDF data array is
!> stored in the 'buffer' array and metadata added to the respective members of
!> the given 'in_var'.
!------------------------------------------------------------------------------!
    INTERFACE get_netcdf_variable
       MODULE PROCEDURE get_netcdf_variable_int
       MODULE PROCEDURE get_netcdf_variable_real
    END INTERFACE get_netcdf_variable

    PRIVATE ::  get_netcdf_variable_int, get_netcdf_variable_real

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> get_netcdf_variable_int() implements the integer variant for the
!> get_netcdf_variable interface.
!------------------------------------------------------------------------------!
 SUBROUTINE get_netcdf_variable_int(in_file, in_var, buffer)

    CHARACTER(LEN=PATH), INTENT(IN)          ::  in_file
    TYPE(nc_var), INTENT(INOUT)              ::  in_var
    INTEGER(iwp), ALLOCATABLE, INTENT(INOUT) ::  buffer(:,:,:)

    INTEGER               ::  ncid
    INTEGER, DIMENSION(3) ::  start, count

    IF ( nf90_open( TRIM(in_file), NF90_NOWRITE, ncid ) .EQ. NF90_NOERR .AND. &
         nf90_inq_varid( ncid, in_var%name, in_var%varid ) .EQ. NF90_NOERR )  THEN

       CALL get_input_dimensions(in_var, ncid)

       CALL get_netcdf_start_and_count(in_var, start, count)
       CALL log_runtime('time', 'read')

       ALLOCATE( buffer( count(1), count(2), count(3) ) )
       CALL log_runtime('time', 'alloc')

       CALL check(nf90_get_var( ncid, in_var%varid, buffer,                  &
                                start = start,                                 &
                                count = count ))

    ELSE

       message = "Failed to read '" // TRIM(in_var%name) // &
          "' from file '" // TRIM(in_file) // "'."
       CALL inifor_abort('get_netcdf_variable', message)

    ENDIF

    CALL check(nf90_close(ncid))
    CALL log_runtime('time', 'read')

 END SUBROUTINE get_netcdf_variable_int


!------------------------------------------------------------------------------!
! Description:
! ------------
!> get_netcdf_variable_real() implements the real variant for the
!> get_netcdf_variable interface.
!------------------------------------------------------------------------------!
 SUBROUTINE get_netcdf_variable_real(in_file, in_var, buffer)

    CHARACTER(LEN=PATH), INTENT(IN)      ::  in_file
    TYPE(nc_var), INTENT(INOUT)          ::  in_var
    REAL(wp), ALLOCATABLE, INTENT(INOUT) ::  buffer(:,:,:)

    INTEGER               ::  ncid
    INTEGER, DIMENSION(3) ::  start, count

    IF ( nf90_open( TRIM(in_file), NF90_NOWRITE, ncid ) .EQ. NF90_NOERR .AND. &
         nf90_inq_varid( ncid, in_var%name, in_var%varid ) .EQ. NF90_NOERR )  THEN

       CALL get_input_dimensions(in_var, ncid)

       CALL get_netcdf_start_and_count(in_var, start, count)
       CALL log_runtime('time', 'read')

       ALLOCATE( buffer( count(1), count(2), count(3) ) )
       CALL log_runtime('time', 'alloc')

       CALL check(nf90_get_var( ncid, in_var%varid, buffer,                  &
                                start = start,                                 &
                                count = count ))

    ELSE

       message = "Failed to read '" // TRIM(in_var%name) // &
          "' from file '" // TRIM(in_file) // "'."
       CALL inifor_abort('get_netcdf_variable', message)

    ENDIF

    CALL check(nf90_close(ncid))
    CALL log_runtime('time', 'read')

 END SUBROUTINE get_netcdf_variable_real


!------------------------------------------------------------------------------!
! Description:
! ------------
!> get_netcdf_dim_vector() reads the coordinate array 'coordname' from the
!> netCDF file 'filename'.
!------------------------------------------------------------------------------!
 SUBROUTINE get_netcdf_dim_vector(filename, coordname, coords)

    CHARACTER(LEN=*), INTENT(IN)         ::  filename
    CHARACTER(LEN=*), INTENT(IN)         ::  coordname
    REAL(wp), ALLOCATABLE, INTENT(INOUT) ::  coords(:)

    INTEGER ::  ncid, varid, dimlen
    INTEGER ::  dimids(NF90_MAX_VAR_DIMS)

    IF ( nf90_open( TRIM(filename), NF90_NOWRITE, ncid ) .EQ. NF90_NOERR .AND. &
         nf90_inq_varid( ncid, coordname, varid ) .EQ. NF90_NOERR )  THEN

       CALL check(nf90_inquire_variable( ncid, varid, dimids = dimids ))
       CALL check(nf90_inquire_dimension( ncid, dimids(1), len = dimlen ))

       ALLOCATE(coords(dimlen))
       CALL check(nf90_get_var( ncid, varid, coords))

    ELSE

       message = "Failed to read '" // TRIM(coordname) // &
          "' from file '" // TRIM(filename) // "'."
       CALL inifor_abort('get_netcdf_dim_vector', message)

    ENDIF

 END SUBROUTINE get_netcdf_dim_vector


!------------------------------------------------------------------------------!
! Description:
! ------------
!> get_input_dimensions() reads dimensions metadata of the netCDF variable given
!> by 'in_var%name' into 'in_var' data structure.
!------------------------------------------------------------------------------!
 SUBROUTINE get_input_dimensions(in_var, ncid)

    TYPE(nc_var), INTENT(INOUT) ::  in_var
    INTEGER, INTENT(IN)         ::  ncid

    INTEGER ::  i

    CALL check(nf90_get_att( ncid, in_var%varid, "long_name",             &
                             in_var%long_name))

    CALL check(nf90_get_att( ncid, in_var%varid, "units",                 &
                             in_var%units ))

    CALL check(nf90_inquire_variable( ncid, in_var%varid,                 &
                                      ndims  = in_var%ndim,               &
                                      dimids = in_var%dimids ))

    DO  i = 1, in_var%ndim
       CALL check(nf90_inquire_dimension( ncid, in_var%dimids(i),         &
                                          name = in_var%dimname(i),       &
                                          len  = in_var%dimlen(i) ))
    ENDDO

 END SUBROUTINE get_input_dimensions


!------------------------------------------------------------------------------!
! Description:
! ------------
!> get_netcdf_start_and_count() gets the start position and element counts for
!> the given netCDF file. This information is used in get_netcdf_variable_int()
!> and _real() for reading input variables..
!------------------------------------------------------------------------------!
 SUBROUTINE get_netcdf_start_and_count(in_var, start, count)

    TYPE(nc_var), INTENT(INOUT)        ::  in_var
    INTEGER, DIMENSION(3), INTENT(OUT) ::  start, count

    INTEGER ::  ndim

    IF ( in_var%ndim .LT. 2  .OR.  in_var%ndim .GT. 4 )  THEN

       message = "Failed reading NetCDF variable " //                       &
          TRIM(in_var%name) // " with " //                                  &
          TRIM(str(INT(in_var%ndim, kind=iwp))) //                          &
          " dimensions because only two- and and three-dimensional" //      &
          " variables are supported."
       CALL inifor_abort('get_netcdf_start_and_count', message)

    ENDIF

    start = (/ 1, 1, 1 /)
    IF ( TRIM(in_var%name) .EQ. 'T_SO' .AND.                                &
         in_var%has_redundant_first_level )  THEN
       
!
!--    Skip depth = 0.0 for T_SO and reduce number of depths from 9 to 8
       in_var%dimlen(3) = in_var%dimlen(3) - 1

!
!--    Start reading from second level, e.g. depth = 0.005 instead of 0.0
       start(3) = 2
    ENDIF

    IF (in_var%ndim .EQ. 2)  THEN
       in_var%dimlen(3) = 1
    ENDIF

    ndim = MIN(in_var%ndim, 3)
    count = (/ 1, 1, 1 /)
    count(1:ndim) = in_var%dimlen(1:ndim)

 END SUBROUTINE get_netcdf_start_and_count


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Routine for defining netCDF variables in the dynamic driver, INIFOR's netCDF
!> output.
!------------------------------------------------------------------------------!
 SUBROUTINE netcdf_define_variable(var, ncid)

    TYPE(nc_var), INTENT(INOUT) ::  var
    INTEGER, INTENT(IN)         ::  ncid

    CALL check(nf90_def_var(ncid, var%name, NF90_FLOAT,       var%dimids(1:var%ndim), var%varid))
    CALL check(nf90_put_att(ncid, var%varid, "long_name",     var%long_name))
    CALL check(nf90_put_att(ncid, var%varid, "units",         var%units))
    IF ( var%lod .GE. 0 )  THEN
       CALL check(nf90_put_att(ncid, var%varid, "lod",           var%lod))
    ENDIF 
    CALL check(nf90_put_att(ncid, var%varid, "source",        var%source))
    CALL check(nf90_put_att(ncid, var%varid, "_FillValue",    NF90_FILL_REAL))

 END SUBROUTINE netcdf_define_variable
    

!------------------------------------------------------------------------------!
! Description:
! ------------
!> netcdf_get_dimensions() reads in all dimensions and their lengths and stores
!> them in the given the 'var' data structure. This information is used later 
!> for writing output variables in update_output().
!------------------------------------------------------------------------------!
 SUBROUTINE netcdf_get_dimensions(var, ncid)

    TYPE(nc_var), INTENT(INOUT) ::  var
    INTEGER, INTENT(IN)         ::  ncid
    INTEGER                     ::  i
    CHARACTER(SNAME)            ::  null

    DO  i = 1, var%ndim
       CALL check(nf90_inquire_dimension(ncid, var%dimids(i), &
                                         name = null, &
                                         len  = var%dimlen(i)  ) )
    ENDDO

 END SUBROUTINE netcdf_get_dimensions


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine parses and interpretes the command-line options and stores the
!> resulting settings in the 'cfg' data structure.
!------------------------------------------------------------------------------!
 SUBROUTINE parse_command_line_arguments( cfg )

    TYPE(inifor_config), INTENT(INOUT) ::  cfg

    CHARACTER(LEN=PATH)                ::  option, arg
    INTEGER                            ::  arg_count, i

    cfg%flow_prefix_is_set = .FALSE.
    cfg%input_prefix_is_set = .FALSE.
    cfg%map_terrain = .FALSE.
    cfg%p0_is_set = .FALSE.
    cfg%precipitation_prefix_is_set = .FALSE.
    cfg%process_precipitation = .FALSE.
    cfg%radiation_prefix_is_set = .FALSE.
    cfg%soil_prefix_is_set = .FALSE.
    cfg%static_driver_is_set = .FALSE.
    cfg%ug_defined_by_user = .FALSE.
    cfg%vg_defined_by_user = .FALSE.
    cfg%z0_is_set = .FALSE.

    arg_count = COMMAND_ARGUMENT_COUNT()
    IF (arg_count .GT. 0)  THEN

       i = 1
       DO  WHILE (i .LE. arg_count)

          CALL GET_COMMAND_ARGUMENT( i, option )

          SELECT CASE( TRIM(option) )

             CASE( '--averaging-mode' )
                CALL get_option_argument( i, arg )
                cfg%averaging_mode = TRIM(arg)

             CASE( '-date', '-d', '--date' )
                CALL get_option_argument( i, arg )
                cfg%start_date = TRIM(arg)

             CASE( '--debug' )
                cfg%debug = .TRUE.

             CASE( '-z0', '-z', '--elevation' )
                cfg%z0_is_set = .TRUE.
                CALL get_option_argument( i, arg )
                READ(arg, *) cfg%z0

             CASE( '-p0', '-r', '--surface-pressure' )
                cfg%p0_is_set = .TRUE.
                CALL get_option_argument( i, arg )
                READ(arg, *) cfg%p0

             CASE( '-ug', '-u', '--geostrophic-u' )
                cfg%ug_defined_by_user = .TRUE.
                CALL get_option_argument( i, arg )
                READ(arg, *) cfg%ug

             CASE( '-vg', '-v', '--geostrophic-v' )
                cfg%vg_defined_by_user = .TRUE.
                CALL get_option_argument( i, arg )
                READ(arg, *) cfg%vg

             CASE( '-path', '-p', '--path' )
                CALL get_option_argument( i, arg )
                 cfg%input_path = TRIM(arg)

             CASE( '-hhl', '-l', '--hhl-file' )
                CALL get_option_argument( i, arg )
                cfg%hhl_file = TRIM(arg)

             CASE( '--input-prefix')
                CALL get_option_argument( i, arg )
                cfg%input_prefix = TRIM(arg)
                cfg%input_prefix_is_set = .TRUE.
   
             CASE( '-a', '--averaging-angle' )
                CALL get_option_argument( i, arg )
                READ(arg, *) cfg%averaging_angle

             CASE( '-static', '-t', '--static-driver' )
                cfg%static_driver_is_set = .TRUE.
                CALL get_option_argument( i, arg )
                cfg%static_driver_file = TRIM(arg)

             CASE( '-soil', '-s', '--soil-file')
                CALL get_option_argument( i, arg )
                cfg%soiltyp_file = TRIM(arg)

             CASE( '--flow-prefix')
                CALL get_option_argument( i, arg )
                cfg%flow_prefix = TRIM(arg)
                cfg%flow_prefix_is_set = .TRUE.
   
             CASE( '--radiation-prefix')
                CALL get_option_argument( i, arg )
                cfg%radiation_prefix = TRIM(arg)
                cfg%radiation_prefix_is_set = .TRUE.
   
             CASE( '--soil-prefix')
                CALL get_option_argument( i, arg )
                cfg%soil_prefix = TRIM(arg)
                cfg%soil_prefix_is_set = .TRUE.
   
             CASE( '--precipitation-prefix')
                CALL get_option_argument( i, arg )
                cfg%precipitation_prefix = TRIM(arg)
                cfg%precipitation_prefix_is_set = .TRUE.

             CASE( '--precipitation')
                cfg%process_precipitation = .TRUE.

             CASE( '-m', '--map-terrain' )
                cfg%map_terrain = .TRUE.

             CASE( '-n', '--namelist' )
                CALL get_option_argument( i, arg )
                cfg%namelist_file = TRIM(arg)

             CASE( '-o', '--output' )
                CALL get_option_argument( i, arg )
                cfg%output_file = TRIM(arg)

             CASE( '-i', '--init-mode' )
                CALL get_option_argument( i, arg )
                cfg%ic_mode = TRIM(arg)

             CASE( '--soil-init-mode' )
                CALL get_option_argument( i, arg )
                cfg%isc_mode = TRIM(arg)

             CASE( '-f', '--forcing-mode' )
                CALL get_option_argument( i, arg )
                cfg%bc_mode = TRIM(arg)

             CASE( '--version' )
                CALL print_version
                STOP

             CASE( '--help' )
                CALL print_version
                PRINT *, ""
                PRINT *, &
                   "For documentation and a list of available command-line options " // NEW_LINE( " " ) // &
                   " please visit https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/inifor."
                STOP

             CASE DEFAULT
                message = "unknown option '" // TRIM(option) // "'."
                CALL inifor_abort('parse_command_line_arguments', message)

          END SELECT

          i = i + 1

       ENDDO

    ELSE
         
       CALL print_version
       CALL report( 'parse_command_line_arguments', 'No arguments present, exiting.' )
       STOP

    ENDIF

 END SUBROUTINE parse_command_line_arguments



 SUBROUTINE get_datetime_file_list( start_date_string, start_hour, end_hour, &
                                    step_hour, input_path, prefix, suffix,   &
                                    file_list )

    CHARACTER (LEN=DATE), INTENT(IN) ::  start_date_string
    CHARACTER (LEN=*),    INTENT(IN) ::  prefix, suffix, input_path
    INTEGER(iwp),         INTENT(IN) ::  start_hour, end_hour, step_hour
    CHARACTER(LEN=*), ALLOCATABLE, INTENT(INOUT) ::  file_list(:)

    INTEGER(iwp)        ::  number_of_intervals, hour, i
    CHARACTER(LEN=DATE) ::  date_string

    number_of_intervals = CEILING( REAL(end_hour - start_hour) / step_hour )
    ALLOCATE( file_list(number_of_intervals + 1) )

    DO  i = 0, number_of_intervals

       hour = start_hour + i * step_hour
       date_string = add_hours_to(start_date_string, hour)

       file_list(i+1) = TRIM(input_path) // TRIM(prefix) //                  &
                        TRIM(date_string) // TRIM(suffix) // '.nc'

    ENDDO

 END SUBROUTINE get_datetime_file_list

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Establish a list of files based on the given start and end times and file
!> prefixes and suffixes.
!------------------------------------------------------------------------------!
 SUBROUTINE get_input_file_list( start_date_string, start_hour, end_hour,    &
                                 step_hour, input_path, prefix, suffix,      &
                                 file_list )

    CHARACTER (LEN=DATE), INTENT(IN) ::  start_date_string
    CHARACTER (LEN=*),    INTENT(IN) ::  prefix, suffix, input_path
    INTEGER(iwp),         INTENT(IN) ::  start_hour, end_hour, step_hour
    CHARACTER(LEN=*), ALLOCATABLE, INTENT(INOUT) ::  file_list(:)

    CALL get_datetime_file_list( start_date_string, start_hour, end_hour,    &
                                 step_hour, input_path, prefix, suffix,      &
                                 file_list )

 END SUBROUTINE get_input_file_list


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Abort INIFOR if the given file is not present.
!------------------------------------------------------------------------------!
LOGICAL FUNCTION file_is_present( filename, file_kind, message )

    CHARACTER(LEN=*), INTENT(IN)            ::  filename, file_kind
    CHARACTER(LEN=*), INTENT(OUT)           ::  message

    file_is_present = file_present(filename)

    IF (.NOT. file_is_present)  THEN

       IF (LEN( TRIM( filename ) ) == 0)  THEN
          message = "No name was given for the " // TRIM( file_kind ) // " file."
       ELSE
          message = "The " // TRIM( file_kind ) // " file '" //                  &
                    TRIM( filename ) // "' was not found."
       ENDIF

    ENDIF

END FUNCTION file_is_present


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Abort INIFOR if the given file is not present.
!------------------------------------------------------------------------------!
 SUBROUTINE verify_file( file_name, file_kind )

    CHARACTER(LEN=*), INTENT(IN)           ::  file_name, file_kind

    IF (.NOT. file_is_present( file_name, file_kind, message ))  THEN

       CALL inifor_abort( 'verify_file', message )

    ENDIF

    message = "Set up input file name '" // TRIM( file_name ) // "'"
    CALL report( 'verify_file', message )

 END SUBROUTINE verify_file
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Get the argument of the i'th command line option, which is at the location
!> i+1 of the argument list.
!------------------------------------------------------------------------------!
 SUBROUTINE get_option_argument( i, arg )
    CHARACTER(LEN=PATH), INTENT(INOUT) ::  arg
    INTEGER, INTENT(INOUT)             ::  i

    i = i + 1
    CALL GET_COMMAND_ARGUMENT( i, arg )

 END SUBROUTINE


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks the INIFOR configuration 'cfg' for plausibility.
!------------------------------------------------------------------------------!
 SUBROUTINE validate_config( cfg )
    TYPE(inifor_config), INTENT(IN) ::  cfg

    CALL verify_file( cfg%hhl_file, 'HHL' )
    CALL verify_file( cfg%namelist_file, 'NAMELIST' )
    CALL verify_file( cfg%soiltyp_file, 'SOILTYP' )

!
!-- Only check optional static driver file name, if it has been given.
    IF (TRIM( cfg%static_driver_file ) .NE. '')  THEN
       CALL verify_file( cfg%static_driver_file, 'static driver' )
    ENDIF

    SELECT CASE( TRIM( cfg%ic_mode ) )
       CASE( CFG_INIT_PROFILE, CFG_INIT_VOLUME )
       CASE DEFAULT
          message = "Initialization mode '" // TRIM( cfg%ic_mode ) //&
                    "' is not supported. " //&
                    "Please select either '" // CFG_INIT_PROFILE //"' or '" // &
                    CFG_INIT_VOLUME //"', " //&
                    "or omit the -i/--init-mode option entirely, which corresponds "//&
                    "to the latter."
          CALL inifor_abort( 'validate_config', message )
    END SELECT

    SELECT CASE( TRIM( cfg%isc_mode ) )
       CASE( CFG_INIT_SOIL_PROFILE, CFG_INIT_SOIL_VOLUME )
       CASE DEFAULT
          message = "Soil initialization mode '" // TRIM( cfg%isc_mode ) //&
                    "' is not supported. " //&
                    "Please select either '" // CFG_INIT_SOIL_PROFILE //"' or '" // &
                    CFG_INIT_SOIL_VOLUME //"', " //&
                    "or omit the --soil-init-mode option entirely, which corresponds "//&
                    "to the latter."
          CALL inifor_abort( 'validate_config', message )
    END SELECT

    SELECT CASE( TRIM(cfg%bc_mode) )
       CASE( CFG_FORCING_HOMO, CFG_FORCING_HETERO, CFG_FORCING_NUDGING )
       CASE DEFAULT
          message = "Forcing mode '" // TRIM( cfg%bc_mode ) //& 
                    "' is not supported. " //&
                    "Please select either '" // CFG_FORCING_NUDGING //          &
                    "', '" // CFG_FORCING_HOMO // "', or '" //                 &
                    CFG_FORCING_HETERO // "' " //&
                    "or omit the -f/--forcing-mode option entirely, which corresponds "//&
                    "to the latter."
          CALL inifor_abort( 'validate_config', message )
    END SELECT

    SELECT CASE( TRIM( cfg%averaging_mode ) )
       CASE( 'level' )
       CASE( 'height' )
          message = "Averaging mode '" // TRIM( cfg%averaging_mode ) //&
                    "' is currently not supported. " //&
                    "Please use level-based averaging by selecting 'level', " //&
                    "or by omitting the --averaging-mode option entirely."
          CALL inifor_abort( 'validate_config', message )
       CASE DEFAULT
          message = "Averaging mode '" // TRIM( cfg%averaging_mode ) //&
                    "' is not supported. " //&
          !          "Please select either 'height' or 'level', " //&
          !          "or omit the --averaging-mode option entirely, which corresponds "//&
          !          "to the latter."
                    "Please use level-based averaging by selecting 'level', " //&
                    "or by omitting the --averaging-mode option entirely."
          CALL inifor_abort( 'validate_config', message )
    END SELECT

    IF ( cfg%ug_defined_by_user .NEQV. cfg%vg_defined_by_user )  THEN
       message = "You specified only one component of the geostrophic " // &
                 "wind. Please specify either both or none."
       CALL inifor_abort( 'validate_config', message )
    ENDIF

    IF ( .NOT. cfg%static_driver_is_set .AND. .NOT. cfg%z0_is_set )  THEN
       message =                                                               &
          "The vertical origin of the PALM grid has not been defined. " // NEW_LINE( " " ) // &
          "Please specify the right value for your setup by either " // NEW_LINE( " " ) // &
          "  - using the command-line option --elevation <height above sea level>, or by" // NEW_LINE( " " ) // &
          "  - specifying a static driver file using --static <filename> in order to use " // NEW_LINE( " " ) // &
          "    use the value of origin_z (and origin_lon and origin_lat) specifed therein."
       CALL inifor_abort( 'validate_config', message )
    ENDIF

    IF ( cfg%map_terrain .AND. .NOT. cfg%static_driver_is_set )  THEN
       message =                                                               &
          "You switched on mesoscale-microscale terrain mapping " //           &
          "(--map-terrain). Please also provide the microscale terrain by " // &
          "specifying the static driver file (--static <filename>)."
       CALL inifor_abort( 'validate_config', message )
    ENDIF

 END SUBROUTINE validate_config


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks wheather the COSMO grid matches the shape of the meteorological input
!> data by comparing the number of netCDF dimensions and their lengths in the
!> hhl.nc and the first of the *-flow files.
!------------------------------------------------------------------------------!
 SUBROUTINE validate_dataset(flow_files, hhl_file)
    CHARACTER(LEN=PATH), INTENT(IN) ::  flow_files(:) !< paths to files containing atmospheric variables
    CHARACTER(LEN=PATH), INTENT(IN) ::  hhl_file      !< path to file containing the HHL variable (height of half layers)

    CHARACTER(SNAME), PARAMETER ::  NC_W_NAME = 'W'
    TYPE(nc_var)                ::  hhl_var, flow_var
    INTEGER                     ::  dim_idx, ncid, ndims_flow, ndims_hhl, varid
    REAL(wp), ALLOCATABLE       ::  hhl_dim_vector(:), flow_dim_vector(:)
    LOGICAL                     ::  dims_have_same_length

    IF ( nf90_open( TRIM( flow_files(1) ), NF90_NOWRITE, ncid ) .EQ. NF90_NOERR )  THEN

       CALL check( nf90_inq_varid( ncid, NC_W_NAME, varid=varid ) )
       CALL check( nf90_inquire_variable( ncid, varid, ndims=ndims_flow ) )
       CALL netcdf_get_dimensions( flow_var, ncid )
       CALL check( nf90_close( ncid ) )

    ELSE

       message = "Failed to read netCDF dimensions'" //                        &
                 "' from file '" // TRIM( flow_files(1) ) // "'."
       CALL inifor_abort( 'validate_dataset', message )

    ENDIF

    IF ( nf90_open( TRIM( hhl_file ), NF90_NOWRITE, ncid ) .EQ. NF90_NOERR )  THEN

       CALL check( nf90_inq_varid( ncid, NC_HHL_NAME, varid=varid ) )
       CALL check( nf90_inquire_variable( ncid, varid, ndims=ndims_hhl ) )
       CALL netcdf_get_dimensions( hhl_var, ncid )
       CALL check( nf90_close( ncid ) )

    ELSE

       message = "Failed to read netCDF dimensions'" //                        &
                 "' from file '" // TRIM(hhl_file) // "'."
       CALL inifor_abort( 'validate_dataset', message )

    ENDIF

!
!-- Compare number dimensions of 'HHL' in hhl file and 'W' in first flow file
    IF  ( .NOT. ndims_flow .EQ. ndims_hhl )  THEN
       message = "Mesoscale data inconsistent. Number of dimensions in the " //&
                 "hhl file does not match with the meteorologial fields " //   &
                 "in the *-flow files (" //                                    &
                 "HHL: ndims = " // TRIM( str( ndims_hhl ) )  // ", " //       &
                 "W: ndims = "   // TRIM( str( ndims_flow ) ) // ")."
       CALL inifor_abort( 'validate_dataset', message )
    ENDIF


!
!-- Compare lengths of each dimension, ignoring time (dim_idx = 1)
    DO dim_idx = 2, ndims_hhl

       CALL get_dimension_vector_of_variable(                                  &
          NC_HHL_NAME,                                                         &
          dim_idx = dim_idx,                                                   &
          filename = hhl_file,                                                 &
          dim_vector = hhl_dim_vector                                          &
       )

       CALL get_dimension_vector_of_variable(                                  &
          NC_W_NAME,                                                           &
          dim_idx = dim_idx,                                                   &
          filename = flow_files(1),                                            &
          dim_vector = flow_dim_vector                                         &
       )

       dims_have_same_length = SIZE( flow_dim_vector ) .EQ. SIZE( hhl_dim_vector )
       IF  ( .NOT. dims_have_same_length )  THEN
          message = &
             "Mesoscale data inconsistent. Number of grid points " //          &
             "in dimension #" // TRIM( str( dim_idx ) ) //                     &
             " do not match in the hhl and *-flow files (" //                  &
             "HHL: len = " // TRIM( str( SIZE( hhl_dim_vector ) ) ) // ", " // &
             "W: len = "   // TRIM( str( SIZE( flow_dim_vector ) ) )// ")."
          CALL inifor_abort( 'validate_dataset', message )
       ENDIF

    ENDDO

 END SUBROUTINE validate_dataset

 SUBROUTINE get_cosmo_grid( hhl_file, soil_file, rlon, rlat, hhl, hfl, depths, &
                            d_depth, d_depth_rho_inv, phi_n, lambda_n,         &
                            phi_equat,                                         &
                            lonmin_cosmo, lonmax_cosmo,                        &
                            latmin_cosmo, latmax_cosmo,                        &
                            nlon, nlat, nlev, ndepths )

    CHARACTER(LEN=PATH), INTENT(IN)                      ::  hhl_file  !< path to file containing the HHL variable (height of half layers)
    CHARACTER(LEN=PATH), INTENT(IN)                      ::  soil_file !< path to one of the soil input files for reading soil layer depths
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)     ::  rlon      !< longitudes of COSMO-DE's rotated-pole grid
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)     ::  rlat      !< latitudes of COSMO-DE's rotated-pole grid
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) ::  hhl       !< heights of half layers (cell faces) above sea level in COSMO-DE, read in from external file
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) ::  hfl       !< heights of full layers (cell centres) above sea level in COSMO-DE, computed as arithmetic average of hhl
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)     ::  depths    !< COSMO-DE's TERRA-ML soil layer depths
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)     ::  d_depth
    REAL(wp), DIMENSION(:), ALLOCATABLE, INTENT(OUT)     ::  d_depth_rho_inv
    REAL(wp), INTENT(OUT)                                ::  phi_n
    REAL(wp), INTENT(OUT)                                ::  phi_equat
    REAL(wp), INTENT(OUT)                                ::  lambda_n
    REAL(wp), INTENT(OUT)                                ::  lonmin_cosmo !< Minimunm longitude of COSMO-DE's rotated-pole grid [COSMO rotated-pole rad]
    REAL(wp), INTENT(OUT)                                ::  lonmax_cosmo !< Maximum longitude of COSMO-DE's rotated-pole grid [COSMO rotated-pole rad]
    REAL(wp), INTENT(OUT)                                ::  latmin_cosmo !< Minimunm latitude of COSMO-DE's rotated-pole grid [COSMO rotated-pole rad]
    REAL(wp), INTENT(OUT)                                ::  latmax_cosmo !< Maximum latitude of COSMO-DE's rotated-pole grid [COSMO rotated-pole rad]
    INTEGER(iwp), INTENT(OUT)                            ::  nlon, nlat, nlev, ndepths

    TYPE(nc_var) ::  cosmo_var !< COSMO dummy variable, used for reading HHL, rlon, rlat
    INTEGER(iwp) ::  k

!
!-- Read in COSMO's heights of half layers (vertical cell faces)
    cosmo_var%name = NC_HHL_NAME
    CALL get_netcdf_variable( hhl_file, cosmo_var, hhl )
    CALL get_netcdf_dim_vector( hhl_file, NC_RLON_NAME, rlon )
    CALL get_netcdf_dim_vector( hhl_file, NC_RLAT_NAME, rlat )
    CALL get_netcdf_dim_vector( soil_file, NC_DEPTH_NAME, depths)
    CALL log_runtime( 'time', 'read' )

    CALL reverse( hhl )
    nlon = SIZE( hhl, 1 )
    nlat = SIZE( hhl, 2 )
    nlev = SIZE( hhl, 3 )
    ndepths = SIZE( depths )

    CALL log_runtime( 'time', 'comp' )

    ALLOCATE( hfl( nlon, nlat, nlev-1 ) )
    ALLOCATE( d_depth( ndepths ), d_depth_rho_inv( ndepths ) )
    CALL log_runtime('time', 'alloc')

    CALL get_soil_layer_thickness( depths, d_depth )

    d_depth_rho_inv(:) = 0.0_wp
    DO  k = 1, ndepths
       IF ( d_depth(k) /= 0.0_wp )                                           &
          d_depth_rho_inv(k) = 1.0_wp / ( d_depth(k) * RHO_L )
    ENDDO

!
!-- Compute COSMO's heights of full layers (cell centres)
    DO  k = 1, nlev-1
       hfl(:,:,k) = 0.5_wp * ( hhl(:,:,k) +                                  &
                               hhl(:,:,k+1) )
    ENDDO
!
!-- COSMO rotated pole coordinates
    phi_n = TO_RADIANS                                                       &
          * get_netcdf_variable_attribute( hhl_file,                         &
                                           NC_ROTATED_POLE_NAME,             &
                                           NC_POLE_LATITUDE_NAME )

    lambda_n = TO_RADIANS                                                    &
             * get_netcdf_variable_attribute( hhl_file,                      &
                                              NC_ROTATED_POLE_NAME,          &
                                              NC_POLE_LONGITUDE_NAME )

    phi_equat = 90.0_wp * TO_RADIANS - phi_n

    lonmin_cosmo = MINVAL( rlon ) * TO_RADIANS
    lonmax_cosmo = MAXVAL( rlon ) * TO_RADIANS
    latmin_cosmo = MINVAL( rlat ) * TO_RADIANS
    latmax_cosmo = MAXVAL( rlat ) * TO_RADIANS
    CALL log_runtime('time', 'comp')

 END SUBROUTINE get_cosmo_grid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Fills the thickness array of the COSMO soil layers. Since COSMO's (i.e. 
!> TERRA_ML's [1]) soil layer boundaries follow the rule
!>
!>    depth(0) = 0.0, and
!>    depth(k) = 0.01 * 3**(k-1), k in [1,2,3,...,7]
!>
!> and full levels are defined as the midpoints between two layer boundaries, 
!> all except the first layer thicknesses equal the depth of the midpoint.
!> 
!> [1] A Description of the Nonhydrostatic Regional COSMO Model Part II :
!>     Physical Parameterization*, Sect. 11 TERRA_ML.
!>     http://www.cosmo-model.org/content/model/documentation/core/cosmoPhysParamtr.pdf)
!> 
!> Input parameters:
!> -----------------
!>
!> depths: array of full soil layer depths (cell midpoints)
!>
!>
!> Output parameters:
!> ------------------
!>
!> d_depth: array of soil layer thicknesses
!>
!------------------------------------------------------------------------------!
 SUBROUTINE get_soil_layer_thickness( depths, d_depth )

    REAL(wp), INTENT(IN)  ::  depths(:)
    REAL(wp), INTENT(OUT) ::  d_depth(:)

    d_depth(:) = depths(:)
    d_depth(1) = 2.0_wp * depths(1)

 END SUBROUTINE get_soil_layer_thickness
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check whether the given file is present on the filesystem.
!------------------------------------------------------------------------------!
 LOGICAL FUNCTION file_present( filename )
    CHARACTER(LEN=PATH), INTENT(IN) ::  filename

    INQUIRE( FILE=filename, EXIST=file_present )

 END FUNCTION file_present


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine initializes the dynamic driver file, i.e. INIFOR's netCDF output
!> file.
!> 
!> Besides writing metadata, such as global attributes, coordinates, variables,
!> in the netCDF file, various netCDF IDs are saved for later, when INIFOR
!> writes the actual data.
!------------------------------------------------------------------------------!
 SUBROUTINE setup_netcdf_dimensions( output_file, palm_grid,                  &
                                     start_date_string, origin_lon, origin_lat )

    TYPE(nc_file), INTENT(INOUT)      ::  output_file
    TYPE(grid_definition), INTENT(IN) ::  palm_grid
    CHARACTER (LEN=DATE), INTENT(IN)  ::  start_date_string
    REAL(wp), INTENT(IN)              ::  origin_lon, origin_lat

    CHARACTER (LEN=8)     ::  date_string
    CHARACTER (LEN=10)    ::  time_string
    CHARACTER (LEN=5)     ::  zone_string
    CHARACTER (LEN=SNAME) ::  history_string
    INTEGER               ::  nx, ny, nz, nt
    INTEGER               ::  ncid, dimids(3), dimvarids(3)
    REAL(wp)              ::  z0

    message = "Initializing PALM-4U dynamic driver file '" //               &
              TRIM(output_file%name) // "' and setting up dimensions."
    CALL report('setup_netcdf_dimensions', message)

!
!-- Create the netCDF file as in netCDF-4/HDF5 format if __netcdf4 preprocessor flag is given
#if defined( __netcdf4 )
    CALL check(nf90_create(TRIM(output_file%name), OR(NF90_CLOBBER, NF90_HDF5), ncid))
#else
    CALL check(nf90_create(TRIM(output_file%name), NF90_CLOBBER, ncid))
#endif

!------------------------------------------------------------------------------
!- Section 1: Define NetCDF dimensions and coordinates
!------------------------------------------------------------------------------
    nt = SIZE(output_file%time)
    nx = palm_grid%nx
    ny = palm_grid%ny
    nz = palm_grid%nz
    z0 = palm_grid%z0


!
!------------------------------------------------------------------------------
!- Section 2: Write global NetCDF attributes
!------------------------------------------------------------------------------
    CALL date_and_time(DATE=date_string, TIME=time_string, ZONE=zone_string)
    history_string =                                                        &
        'Created on '// date_string      //                                 &
        ' at '       // time_string(1:2) // ':' // time_string(3:4) //      &
        ' (UTC'      // zone_string // ')'

    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'title',          'PALM input file for scenario ...'))
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'institution',    'Deutscher Wetterdienst, Offenbach'))
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'author',         'Eckhard Kadasch, eckhard.kadasch@dwd.de'))
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'history',        TRIM(history_string)))
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'references',     '--'))
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'comment',        '--'))
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'origin_lat',     TRIM(real_to_str(origin_lat*TO_DEGREES, '(F18.13)'))))
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'origin_lon',     TRIM(real_to_str(origin_lon*TO_DEGREES, '(F18.13)'))))
!
!-- FIXME: This is the elevation relative to COSMO-DE/D2 sea level and does
!-- FIXME: not necessarily comply with DHHN2016 (c.f. PALM Input Data
!-- FIXME: Standard v1.9., origin_z)
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'origin_z',       TRIM(real_to_str(z0, '(F18.13)'))))
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'inifor_version', TRIM(VERSION)))
    CALL check(nf90_put_att(ncid, NF90_GLOBAL, 'palm_version',   '--'))

!
!
!------------------------------------------------------------------------------
!- Section 2a: Define dimensions for cell centers (scalars in soil and atmosph.)
!------------------------------------------------------------------------------
!
!-- reset dimids first
    dimids = (/0, 0, 0/)
    CALL check( nf90_def_dim(ncid, "x", nx+1, dimids(1)) )
    CALL check( nf90_def_dim(ncid, "y", ny+1, dimids(2)) )
    CALL check( nf90_def_dim(ncid, "z", nz, dimids(3)) )
!
!-- save dimids for later
    output_file%dimids_scl = dimids 

!
!-- reset dimvarids first
    dimvarids = (/0, 0, 0/)
    CALL check(nf90_def_var(ncid, "x", NF90_FLOAT, dimids(1), dimvarids(1)))
    CALL check(nf90_put_att(ncid, dimvarids(1), "standard_name", "x coordinate of cell centers"))
    CALL check(nf90_put_att(ncid, dimvarids(1), "units", "m"))

    CALL check(nf90_def_var(ncid, "y", NF90_FLOAT, dimids(2), dimvarids(2)))
    CALL check(nf90_put_att(ncid, dimvarids(2), "standard_name", "y coordinate of cell centers"))
    CALL check(nf90_put_att(ncid, dimvarids(2), "units", "m"))

    CALL check(nf90_def_var(ncid, "z", NF90_FLOAT, dimids(3), dimvarids(3)))
    CALL check(nf90_put_att(ncid, dimvarids(3), "standard_name", "z coordinate of cell centers"))
    CALL check(nf90_put_att(ncid, dimvarids(3), "units", "m"))
!
!-- save dimvarids for later
    output_file%dimvarids_scl = dimvarids

!
!-- overwrite third dimid with the one of depth
    CALL check(nf90_def_dim(ncid, "zsoil", SIZE(palm_grid%depths), dimids(3)) )
!
!-- save dimids for later
    output_file%dimids_soil = dimids

!
!-- overwrite third dimvarid with the one of depth
    CALL check(nf90_def_var(ncid, "zsoil", NF90_FLOAT, output_file%dimids_soil(3), dimvarids(3)))
    CALL check(nf90_put_att(ncid, dimvarids(3), "standard_name", "depth_below_land"))
    CALL check(nf90_put_att(ncid, dimvarids(3), "positive", "down"))
    CALL check(nf90_put_att(ncid, dimvarids(3), "units", "m"))
!
!-- save dimvarids for later
    output_file%dimvarids_soil = dimvarids
!
!------------------------------------------------------------------------------
!- Section 2b: Define dimensions for cell faces/velocities
!------------------------------------------------------------------------------
!
!-- reset dimids first
    dimids = (/0, 0, 0/)
    CALL check(nf90_def_dim(ncid, "xu", nx, dimids(1)) )
    CALL check(nf90_def_dim(ncid, "yv", ny, dimids(2)) )
    CALL check(nf90_def_dim(ncid, "zw", nz-1, dimids(3)) )
!
!-- save dimids for later
    output_file%dimids_vel = dimids

!
!-- reset dimvarids first
    dimvarids = (/0, 0, 0/)
    CALL check(nf90_def_var(ncid, "xu", NF90_FLOAT, dimids(1), dimvarids(1)))
    CALL check(nf90_put_att(ncid, dimvarids(1), "standard_name", "x coordinate of cell faces"))
    CALL check(nf90_put_att(ncid, dimvarids(1), "units", "m"))

    CALL check(nf90_def_var(ncid, "yv", NF90_FLOAT, dimids(2), dimvarids(2)))
    CALL check(nf90_put_att(ncid, dimvarids(2), "standard_name", "y coordinate of cell faces"))
    CALL check(nf90_put_att(ncid, dimvarids(2), "units", "m"))

    CALL check(nf90_def_var(ncid, "zw", NF90_FLOAT, dimids(3), dimvarids(3)))
    CALL check(nf90_put_att(ncid, dimvarids(3), "standard_name", "z coordinate of cell faces"))
    CALL check(nf90_put_att(ncid, dimvarids(3), "units", "m"))
!
!-- save dimvarids for later
    output_file%dimvarids_vel = dimvarids

!
!------------------------------------------------------------------------------
!- Section 2c: Define time dimension
!------------------------------------------------------------------------------
    CALL check(nf90_def_dim(ncid, "time", nt, output_file%dimid_time) )
    CALL check(nf90_def_var(ncid, "time", NF90_FLOAT, &
                                          output_file%dimid_time, &
                                          output_file%dimvarid_time))
    CALL check(nf90_put_att(ncid, output_file%dimvarid_time, "standard_name", "time"))
    CALL check(nf90_put_att(ncid, output_file%dimvarid_time, "long_name", "time"))
    CALL check(nf90_put_att(ncid, output_file%dimvarid_time, "units",     &
                            "seconds since " // start_date_string // " UTC"))

    CALL check(nf90_enddef(ncid))

!
!------------------------------------------------------------------------------
!- Section 3: Write grid coordinates
!------------------------------------------------------------------------------
    CALL check(nf90_put_var(ncid, output_file%dimvarids_scl(1), palm_grid%x))
    CALL check(nf90_put_var(ncid, output_file%dimvarids_scl(2), palm_grid%y))
    CALL check(nf90_put_var(ncid, output_file%dimvarids_scl(3), palm_grid%z))

    CALL check(nf90_put_var(ncid, output_file%dimvarids_vel(1), palm_grid%xu))
    CALL check(nf90_put_var(ncid, output_file%dimvarids_vel(2), palm_grid%yv))
    CALL check(nf90_put_var(ncid, output_file%dimvarids_vel(3), palm_grid%zw))

!
!-- TODO Read in soil depths from input file before this.
    CALL check(nf90_put_var(ncid, output_file%dimvarids_soil(3), palm_grid%depths))

!
!-- Write time vector
    CALL check(nf90_put_var(ncid, output_file%dimvarid_time, output_file%time))

!
!-- Close the file
    CALL check(nf90_close(ncid))

 END SUBROUTINE setup_netcdf_dimensions


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Defines the netCDF variables to be written to the dynamic driver file
!------------------------------------------------------------------------------!
 SUBROUTINE setup_netcdf_variables(filename, io_group_list)

    CHARACTER (LEN=*), INTENT(IN)       ::  filename
    TYPE(io_group), INTENT(IN), TARGET  ::  io_group_list(:)

    TYPE(io_group), POINTER             ::  group
    TYPE(nc_var), POINTER               ::  var
    INTEGER(iwp)                        ::  group_idx, var_idx, n_var
    INTEGER                             ::  ncid
    LOGICAL                             ::  to_be_written

    message = "Defining variables in dynamic driver '" // TRIM(filename) // "'."
    CALL report('setup_netcdf_variables', message)

    CALL check(nf90_open(TRIM(filename), NF90_WRITE, ncid))
    CALL check(nf90_redef(ncid))

    n_var = 0
    DO  group_idx = 1, SIZE( io_group_list )

       group => io_group_list(group_idx)
       DO var_idx = 1, SIZE( group%out_vars )

          to_be_written = .FALSE.

          IF (ALLOCATED( group%out_vars ))  THEN
             var => group%out_vars(var_idx)
             n_var = n_var + 1

             to_be_written = (                                                 &
                group%to_be_processed  .AND.  var%to_be_processed              &
                .AND. .NOT. var%is_internal                                    &
             )
          ENDIF

          IF ( to_be_written )  THEN
             message = "  variable #" // TRIM(str(n_var)) // " '" // TRIM(var%name) // "'."
             CALL report('setup_netcdf_variables', message)

             CALL netcdf_define_variable(var, ncid)
             CALL netcdf_get_dimensions(var, ncid)
          ENDIF

       ENDDO
        
    ENDDO

    CALL check(nf90_enddef(ncid))
    CALL check(nf90_close(ncid))

    message = "Dynamic driver '" // TRIM(filename) // "' initialized successfully."
    CALL report('setup_netcdf_variables', message)

 END SUBROUTINE setup_netcdf_variables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine reads and returns all input variables of the given IO group
!> It accomodates the data by allocating a container variable that represents a
!> list of arrays of the same length as the groups variable list. (This list
!> will typically contain one or two items.) After the container, its members
!> are allocated one by one with the appropriate, possibly different,
!> dimensions.
!>
!> The 'group' is an INTENT(INOUT) variable so that 'get_netcdf_variable()' can
!> record netCDF IDs in the 'in_var_list()' member variable.
!------------------------------------------------------------------------------!
 SUBROUTINE read_input_variables(group, iter, buffer)
    TYPE(io_group), INTENT(INOUT), TARGET       ::  group
    INTEGER(iwp), INTENT(IN)                    ::  iter
    TYPE(container), ALLOCATABLE, INTENT(INOUT) ::  buffer(:)
    INTEGER(iwp)                                ::  hour, buf_id
    TYPE(nc_var), POINTER                       ::  input_var
    CHARACTER(LEN=PATH), POINTER                ::  input_file
    INTEGER(iwp)                                ::  ivar, nbuffers

    message = "Reading data for I/O group '" // TRIM(group%in_var_list(1)%name) // "'."
    CALL report('read_input_variables', message)

    input_file => group%in_files(iter)

!
!------------------------------------------------------------------------------
!- Section 1: Load input buffers for accumulated variables
!------------------------------------------------------------------------------
!
!-- radiation budgets, precipitation
    IF ( group%kind == 'running average' .OR.                                  &
         group%kind == 'accumulated' )  THEN

       IF ( SIZE( group%in_var_list ) .GT. 1 ) THEN
          message = "I/O groups may not contain more than one " // & 
                    "accumulated variable. Group '" // TRIM(group%kind) //&
                    "' contains " //                                           &
                    TRIM( str( SIZE( group%in_var_list, kind=iwp ) ) ) // "."
          CALL inifor_abort( 'read_input_variables | accumulation', message )
       ENDIF

!
!--    use two buffer arrays
       nbuffers = 2
       IF ( .NOT. ALLOCATED( buffer ) )  ALLOCATE( buffer(nbuffers) )

!
!--    hour of the day
       hour = iter - 1
!
!--    chose correct buffer array
       buf_id = select_buffer(hour)

       CALL log_runtime('time', 'read')
       IF ( ALLOCATED(buffer(buf_id)%array) )  DEALLOCATE(buffer(buf_id)%array)
       CALL log_runtime('time', 'alloc')

       input_var => group%in_var_list(1)
       CALL get_netcdf_variable(input_file, input_var, buffer(buf_id)%array)
       CALL report('read_input_variables', "Read accumulated " // TRIM(group%in_var_list(1)%name)) 

       IF ( input_var%is_upside_down )  CALL reverse(buffer(buf_id)%array)
       CALL log_runtime('time', 'comp')
          
!------------------------------------------------------------------------------
!- Section 2: Load input buffers for normal I/O groups
!------------------------------------------------------------------------------
    ELSE

!
!--    Allocate one input buffer per output quantity. If more quantities
!--    have to be computed than input variables exist in this group,
!--    allocate more buffers. For instance, in the thermodynamics group,
!--    there are three input variabels (PP, T, Qv) and four quantities
!--    necessart (P, Theta, Rho, qv) for the corresponding output fields
!--    (p0, Theta, qv, ug, and vg)
       ALLOCATE( buffer(group%n_output_quantities) )
       CALL log_runtime('time', 'alloc')
       
!
!--    Read in all input variables, leave extra buffers-if any-untouched.
       DO  ivar = 1, group%n_inputs

          input_var => group%in_var_list(ivar)

          IF ( input_var%to_be_processed )  THEN
!            Check wheather P or PP is present in input file
             IF (input_var%name == 'P')  THEN
                input_var%name = TRIM( get_pressure_varname(input_file) )
             CALL log_runtime('time', 'read')
             ENDIF

             CALL get_netcdf_variable(input_file, input_var, buffer(ivar)%array)

             IF ( input_var%is_upside_down )  CALL reverse(buffer(ivar)%array)
             CALL log_runtime('time', 'comp')
          ENDIF

       ENDDO
    ENDIF

 END SUBROUTINE read_input_variables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Select the appropriate buffer ID for accumulated COSMO input variables
!> depending on the current hour.
!------------------------------------------------------------------------------!
 INTEGER(iwp) FUNCTION select_buffer(hour)
    INTEGER(iwp), INTENT(IN) ::  hour
    INTEGER(iwp)             ::  step

    select_buffer = 0_iwp
    step = MODULO(hour, 3_iwp) + 1_iwp

    SELECT CASE(step)
       CASE(1, 3)
           select_buffer = 1_iwp
       CASE(2)
           select_buffer = 2_iwp
       CASE DEFAULT
           message = "Invalid step '" // TRIM(str(step))
           CALL inifor_abort('select_buffer', message)
    END SELECT
 END FUNCTION select_buffer


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if the input_file contains the absolute pressure, 'P', or the pressure
!> perturbation, 'PP', and returns the appropriate string.
!------------------------------------------------------------------------------!
 CHARACTER(LEN=2) FUNCTION get_pressure_varname(input_file) RESULT(var)
    CHARACTER(LEN=*) ::  input_file
    INTEGER          ::  ncid, varid

    CALL check(nf90_open( TRIM(input_file), NF90_NOWRITE, ncid ))
    IF ( nf90_inq_varid( ncid, 'P', varid ) .EQ. NF90_NOERR )  THEN

       var = 'P'

    ELSE IF ( nf90_inq_varid( ncid, 'PP', varid ) .EQ. NF90_NOERR )  THEN

       var = 'PP'
       CALL report('get_pressure_var', 'Using PP instead of P')

    ELSE

       message = "Failed to read '" // TRIM(var) // &
                 "' from file '" // TRIM(input_file) // "'."
       CALL inifor_abort('get_pressure_var', message)

    ENDIF

    CALL check(nf90_close(ncid))

 END FUNCTION get_pressure_varname


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read the given global attribute form the given netCDF file.
!------------------------------------------------------------------------------!
 FUNCTION get_netcdf_attribute(filename, attribute) RESULT(attribute_value)

    CHARACTER(LEN=*), INTENT(IN) ::  filename, attribute
    REAL(wp)                     ::  attribute_value

    INTEGER                      ::  ncid

    IF ( nf90_open( TRIM(filename), NF90_NOWRITE, ncid ) == NF90_NOERR )  THEN

       CALL check(nf90_get_att(ncid, NF90_GLOBAL, TRIM(attribute), attribute_value))
       CALL check(nf90_close(ncid))

    ELSE

       message = "Failed to read '" // TRIM(attribute) // &
                 "' from file '" // TRIM(filename) // "'."
       CALL inifor_abort('get_netcdf_attribute', message)

    ENDIF

 END FUNCTION get_netcdf_attribute


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read the attribute of the given variable form the given netCDF file.
!------------------------------------------------------------------------------!
 FUNCTION get_netcdf_variable_attribute(filename, varname, attribute)       &
    RESULT(attribute_value)

    CHARACTER(LEN=*), INTENT(IN) ::  filename, varname, attribute
    REAL(wp)                     ::  attribute_value

    INTEGER                      ::  ncid, varid

    IF ( nf90_open( TRIM(filename), NF90_NOWRITE, ncid ) == NF90_NOERR )  THEN

       CALL check( nf90_inq_varid( ncid, TRIM( varname ), varid ) )
       CALL check( nf90_get_att( ncid, varid, TRIM( attribute ),            &
                   attribute_value ) )
       CALL check( nf90_close( ncid ) )

    ELSE

       message = "Failed to read '" // TRIM( varname ) // ":" //            &
                 TRIM( attribute ) // "' from file '" // TRIM(filename) // "'."
       CALL inifor_abort('get_netcdf_variable_attribute', message)

    ENDIF

 END FUNCTION get_netcdf_variable_attribute

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates the dynamic driver with the interpolated field of the current
!> variable at the current time step.
!------------------------------------------------------------------------------!
 SUBROUTINE update_output(var, array, iter, output_file, cfg)
    TYPE(nc_var), INTENT(IN)  ::  var
    REAL(wp), INTENT(IN)      ::  array(:,:,:)
    INTEGER(iwp), INTENT(IN)  ::  iter
    TYPE(nc_file), INTENT(IN) ::  output_file
    TYPE(inifor_config)       ::  cfg

    INTEGER      ::  ncid, ndim, start(4), count(4)
    LOGICAL      ::  var_is_time_dependent

    var_is_time_dependent = (                                                  &
       var%dimids( var%ndim ) == output_file%dimid_time                        &
    )

!
!-- Skip time dimension for output
    ndim = var%ndim
    IF ( var_is_time_dependent )  ndim = var%ndim - 1

    start(:)      = (/1,1,1,1/)
    start(ndim+1) = iter
    count(1:ndim) = var%dimlen(1:ndim)

    CALL check(nf90_open(output_file%name, NF90_WRITE, ncid))

!
!-- Reduce dimension of output array according to variable kind
    SELECT CASE (TRIM(var%kind))
       
       CASE( 'init scalar profile', 'init u profile', 'init v profile',        &
             'init w profile', 'init soil profile' )

          CALL check(nf90_put_var( ncid, var%varid, array(1,1,:) ) )

       CASE( 'init soil', 'init scalar', 'init u', 'init v', 'init w' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,:,:) ) )

       CASE( 'left scalar', 'right scalar', 'left w', 'right w' ) 

          CALL check(nf90_put_var( ncid, var%varid, array(1,:,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )
          

          IF (.NOT. SIZE(array, 2) .EQ. var%dimlen(1))  THEN
             PRINT *, "inifor: update_output: Dimension ", 1, " of variable ", &
                 TRIM(var%name), " (", var%dimlen(1),                          &
                 ") does not match the dimension of the output array (",       &
                 SIZE(array, 2), ")."
             STOP
          ENDIF
          

       CASE( 'north scalar', 'south scalar', 'north w', 'south w' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,1,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )
          

       CASE( 'surface forcing', 'top scalar', 'top w' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,:,1),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )
          
       CASE( 'left u', 'right u', 'left v', 'right v' )

          CALL check(nf90_put_var( ncid, var%varid, array(1,:,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )

       CASE( 'north u', 'south u', 'north v', 'south v' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,1,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )

       CASE( 'top u', 'top v' )

          CALL check(nf90_put_var( ncid, var%varid, array(:,:,1),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )

       CASE( 'time series' )

          CALL check(nf90_put_var( ncid, var%varid, array(1,1,1),              &
                                   start=start(1:ndim+1) ) )

       CASE( 'constant scalar profile', 'geostrophic',                         &
             'left scalar profile', 'right scalar profile',                    &
             'north scalar profile', 'south scalar profile',                   &
             'left u profile', 'right u profile',                              &
             'north u profile', 'south u profile',                             &
             'left v profile', 'right v profile',                              &
             'north v profile', 'south v profile',                             &
             'top scalar profile', 'top u profile', 'top v profile',           &
             'top w profile',                                                  &
             'left w profile', 'right w profile',                              &
             'north w profile', 'south w profile' )

          CALL check(nf90_put_var( ncid, var%varid, array(1,1,:),              &
                                   start=start(1:ndim+1),                      &
                                   count=count(1:ndim) ) )

       CASE ( 'internal profile' )

          IF ( cfg%debug )  THEN
             CALL check(nf90_put_var( ncid, var%varid, array(1,1,:),           &
                                      start=start(1:ndim+1),                   &
                                      count=count(1:ndim) ) )
          ENDIF

       CASE ( 'large-scale scalar forcing', 'large-scale w forcing' )

           message = "Doing nothing in terms of writing large-scale forings."
           CALL report('update_output', message)

       CASE DEFAULT

           message = "Variable kind '" // TRIM(var%kind) //                  &
                    "' not recognized."
           CALL inifor_abort('update_output', message)

    END SELECT

    CALL check(nf90_close(ncid))

 END SUBROUTINE update_output


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks the status of a netCDF API call and aborts if an error occured
!------------------------------------------------------------------------------!
 SUBROUTINE check(status)

    INTEGER, INTENT(IN) ::  status

    IF (status /= nf90_noerr)  THEN
       message = "NetCDF API call failed with error: " //                     &
                 TRIM( nf90_strerror(status) )
       CALL inifor_abort('io.check', message)
    ENDIF

 END SUBROUTINE check


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Setup the origin of the PALM coordinate system based on what is given in the
!> INIFOR namelist file and via an optional static driver.
!------------------------------------------------------------------------------!
 SUBROUTINE set_palm_origin( cfg, namelist_longitude, namelist_latitude,       &
                             origin_lon, origin_lat, z0 )

    TYPE(inifor_config), INTENT(IN) ::  cfg
    REAL(wp), INTENT(IN)            ::  namelist_longitude, namelist_latitude
    REAL(wp), INTENT(OUT)           ::  origin_lon, origin_lat, z0

    message = 'Reading PALM-4U origin from'
    IF ( TRIM( cfg%static_driver_file ) .NE. '' )  THEN

       origin_lon = get_netcdf_attribute( cfg%static_driver_file, TRIM( PIDS_ORIGIN_LON ) )
       origin_lat = get_netcdf_attribute( cfg%static_driver_file, TRIM( PIDS_ORIGIN_LAT ) )
       z0         = get_netcdf_attribute( cfg%static_driver_file, TRIM( PIDS_ORIGIN_Z ) )

       message = TRIM(message) // " static driver file '"                      &
                               // TRIM( cfg%static_driver_file ) // "'"


    ELSE

       origin_lon = namelist_longitude
       origin_lat = namelist_latitude

       message = TRIM( message ) // " namlist file '"                          &
                                 // TRIM( cfg%namelist_file ) // "'"

    ENDIF
    origin_lon = origin_lon * TO_RADIANS
    origin_lat = origin_lat * TO_RADIANS

    CALL report('set_palm_origin', message)

    IF ( cfg%z0_is_set )  THEN
       z0 = cfg%z0
       IF ( TRIM( cfg%static_driver_file ) .NE. '' )  THEN
          message = 'You specified both --static-driver/-t and --elevation/-z0. ' // &
                    'Using the command line value (' // TRIM( real_to_str_f( cfg%z0 ) ) // &
                    ') instead of static driver value (' // TRIM( real_to_str_f( z0 ) ) // ').'
          CALL warn( 'set_palm_origin', message )
       ENDIF
    ENDIF

 END SUBROUTINE set_palm_origin


!------------------------------------------------------------------------------!
! Description:
! ------------
! This function is meant to check weather a COSMO soil variable has an
! additional and redunant surface value at depth = 0.0. For instance operational
! DWD COSMO output contains the surface temperature in T_SO as a copy of the
! values in the first soil layer.
!------------------------------------------------------------------------------!
 LOGICAL FUNCTION has_surface_value( soil_var, filename ) 

    TYPE(nc_var), INTENT(IN)     ::  soil_var
    CHARACTER(LEN=*), INTENT(IN) ::  filename
    
    REAL(wp), ALLOCATABLE        ::  depths(:)

    CALL get_dimension_vector_of_variable(                                     &
       soil_var%name,                                                          &
       dim_idx = NC_DEPTH_DIM_IDX,                                             &
       filename = filename,                                                    &
       dim_vector = depths                                                     &
    )

    has_surface_value = nearly_equal( depths(1), 0.0_wp, 10 * EPSILON(1.0_wp) )

 END FUNCTION has_surface_value 


!------------------------------------------------------------------------------!
! Description:
! ------------
! This routine reads the dim_idx-th dimension vector of the variable varname
! from netCDF file filename. It is used for finding the depth coordinate vector
! of COSMO soil variables without knowing its name.
!------------------------------------------------------------------------------!
 SUBROUTINE get_dimension_vector_of_variable( varname, dim_idx, filename, dim_vector )
    CHARACTER(LEN=*), INTENT(IN)              ::  varname, filename
    INTEGER, INTENT(IN)                       ::  dim_idx

    REAL(wp), INTENT(OUT), ALLOCATABLE ::  dim_vector(:)

    INTEGER                            ::  dimids(NF90_MAX_VAR_DIMS)
    INTEGER                            ::  varid
    CHARACTER(LEN=NF90_MAX_NAME)       ::  dimname

    INTEGER ::  ncid

    IF ( nf90_open( TRIM( filename ), NF90_NOWRITE, ncid ) .EQ. NF90_NOERR )  THEN

!
!--    get id of variable varname
       CALL check( nf90_inq_varid( ncid, TRIM( varname ), varid ) )
       
!
!--    get dimension ids of variable with varid
       CALL check( nf90_inquire_variable( ncid, varid, dimids = dimids ) )

!
!--    get name of dim_idx-th dimension variable
       CALL check( nf90_inquire_dimension( ncid, dimids(dim_idx), name = dimname ) )
       CALL check( nf90_close( ncid ) )

    ELSE

       message = "Failed to open file '" // TRIM(filename) // "'."
       CALL inifor_abort('get_netcdf_variable', message)

    ENDIF

    ! get dimension vector with dimname
    CALL get_netcdf_dim_vector( filename, dimname, dim_vector )

 END SUBROUTINE get_dimension_vector_of_variable


 LOGICAL FUNCTION netcdf_variable_present_in_file( varname, filename )
    CHARACTER(LEN=SNAME), INTENT(IN) ::  varname
    CHARACTER(LEN=PATH), INTENT(IN)  ::  filename

    INTEGER ::  ncid, varid
    
    netcdf_variable_present_in_file = (                                        &
       nf90_open( TRIM( filename ), NF90_NOWRITE, ncid ) .EQ. NF90_NOERR .AND. &
       nf90_inq_varid( ncid, varname, varid ) .EQ. NF90_NOERR .AND.            &
       nf90_close( ncid ) .EQ. NF90_NOERR                                      &
    )
 
 END FUNCTION netcdf_variable_present_in_file


 END MODULE inifor_io
#endif
