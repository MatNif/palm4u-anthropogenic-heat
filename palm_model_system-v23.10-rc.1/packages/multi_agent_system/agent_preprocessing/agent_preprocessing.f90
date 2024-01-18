!> @agent_preprocessing.f90
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
! Copyright 1997-2018 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Reads topography and building data and converts this information into a
!> navigation mesh (NavMesh, a visibility graph). This mesh is necessary for
!> the use of the Multi Agent System in PALM for the agents to navigate in
!> complex (urban) terrain.
!------------------------------------------------------------------------------!

 MODULE kinds

    IMPLICIT NONE

!
!-- Floating point kinds
    INTEGER, PARAMETER ::  sp = 4           !< single precision (32 bit)
    INTEGER, PARAMETER ::  dp = 8           !< double precision (64 bit)

!
!-- Integer kinds
    INTEGER, PARAMETER ::  isp = SELECTED_INT_KIND(  9 )   !< single precision (32 bit)
    INTEGER, PARAMETER ::  idp = SELECTED_INT_KIND( 14 )   !< double precision (64 bit)

!
!-- Set kinds to be used as defaults
    INTEGER, PARAMETER ::   wp =  dp          !< default real kind
    INTEGER, PARAMETER ::  iwp = isp          !< default integer kind

    SAVE

 END MODULE kinds

 MODULE variables

    USE kinds

    CHARACTER(LEN=3)  ::  char_lod  = 'lod'         !< name of level-of-detail attribute in NetCDF file
    CHARACTER(LEN=10) ::  char_fill = '_FillValue'  !< name of fill value attribute in NetCDF file
    CHARACTER(LEN=128) ::  runname                  !< Run name

    LOGICAL ::  internal_buildings = .FALSE.  !< Flag that indicates whether buildings within closed courtyards should be deleted
    LOGICAL ::  flag_2d            = .FALSE.  !< Flag that indicates that 2d buildings will be used in all cases

    INTEGER(iwp) ::  i                          !< Index along x
    INTEGER(iwp) ::  j                          !< Index along y
    INTEGER(iwp) ::  nx = 99999                 !< Number of grid points in x-direction
    INTEGER(iwp) ::  ny = 99999                 !< Number of grid points in x-direction
    INTEGER(iwp) ::  nov                        !< Number of vertices
    INTEGER(iwp) ::  polygon_counter            !< Iterator for the number of building polygons
    INTEGER(iwp) ::  number_of_connections = 0  !< Counter for number of connections in mesh
    INTEGER(iwp) ::  i_cn                       !< Min number of corners left in polygons after Douglas Poiker algorithm
    INTEGER(iwp) ::  i_sc                       !< Cycle number for Douglas-Peucker algorithm
    INTEGER(iwp) ::  nc_stat                    !< return value of nf90 function call
    INTEGER(iwp) ::  vertex_counter             !< Counter: total number of vertices

    INTEGER, DIMENSION(:,:), ALLOCATABLE ::  wall_flags_0  !< Bit-array containing surface information
    INTEGER, DIMENSION(:,:), ALLOCATABLE ::  polygon_id    !< Identifies each grid point as part of exactly one building

    REAL(wp)    ::  ddx              !< inverse of dx
    REAL(wp)    ::  ddy              !< inverse of dy
    REAL(wp)    ::  dx = 99999.9_wp  !< grid spacing in x-direction
    REAL(wp)    ::  dy = 99999.9_wp  !< grid spacing in x-direction
    REAL(wp)    ::  dz = 99999.9_wp  !< grid spacing in x-direction
    REAL(wp)    ::  finish           !< variable for CPU time measurement
    REAL(wp)    ::  start            !< variable for CPU time measurement

    REAL(wp), DIMENSION(0:2) ::  tolerance_dp = 999999.0_wp  !< tolerance in Douglas-Peucker algorithm

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  obstacle_height  !< height of obstacles

!
!-- Define data structure where the dimension and type of the input depends 
!-- on the given level of detail. 
!-- For buildings, the input is either 2D float, or 3d byte. 
       TYPE build_in
          INTEGER(iwp)    ::  lod = 1                                  !< level of detail                  
          INTEGER(KIND=1) ::  fill2 = -127                             !< fill value for lod = 2
          INTEGER(iwp)    ::  nz                                       !< number of vertical layers in file
          INTEGER(KIND=1), DIMENSION(:,:,:), ALLOCATABLE ::  var_3d    !< 3d variable (lod = 2)
          REAL(wp), DIMENSION(:), ALLOCATABLE ::  z                    !< vertical coordinate for 3D building, used for consistency check
          LOGICAL ::  from_file = .FALSE.                              !< flag indicating whether an input variable is available and read from file or default values are used  
          REAL(wp)                              ::  fill1 = -9999.9_wp !< fill values for lod = 1
          REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var_2d             !< 2d variable (lod = 1)
       END TYPE build_in

!
!-- Topography grid point
    TYPE  grid_point
        LOGICAL      ::  checked     !< Flag to indicate whether this grid point has been evaluated already
        INTEGER(iwp) ::  i           !< x-index
        INTEGER(iwp) ::  j           !< y-index
        INTEGER(iwp) ::  polygon_id  !< ID of the polygon this grid point belongs to
    END TYPE grid_point

!
!-- Node in the visibility graph navigation mesh
    TYPE  mesh_point
        INTEGER(iwp)                            ::  polygon_id          !< Polygon the point belongs to
        INTEGER(iwp)                            ::  vertex_id           !< Vertex in the polygon
        INTEGER(iwp)                            ::  noc                 !< number of connections
        INTEGER(iwp)                            ::  origin_id           !< ID of previous mesh point on path (A*)
        REAL(wp)                                ::  cost_so_far         !< Cost to reach this mesh point (A*)
        REAL(wp)                                ::  x                   !< x-coordinate
        REAL(wp)                                ::  y                   !< y-coordinate
        REAL(wp)                                ::  x_s                 !< corner shifted outward from building by 1m (x)
        REAL(wp)                                ::  y_s                 !< corner shifted outward from building by 1m (y)
        INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  connected_vertices  !< Index of connected vertices
        REAL(wp), DIMENSION(:), ALLOCATABLE     ::  distance_to_vertex  !< Distance to each vertex
    END TYPE mesh_point

!
!-- Vertex of a polygon
    TYPE  vertex_type
        LOGICAL               ::  delete  !< Flag to mark vertex for deletion
        REAL(wp)              ::  x       !< x-coordinate
        REAL(wp)              ::  y       !< y-coordinate
    END TYPE vertex_type

!
!-- Polygon containing a number of vertices
    TYPE  polygon_type
        INTEGER(iwp)                                 ::  nov       !< Number of vertices in this polygon
        TYPE(vertex_type), DIMENSION(:), ALLOCATABLE ::  vertices  !< Array of vertices
    END TYPE polygon_type

!
!-- Define data type to read 2D real variables
    TYPE real_2d
       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used  
       REAL(wp) ::  fill = -9999.9_wp                !< fill value
       REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  var !< respective variable
    END TYPE real_2d

!
!-- Define data type to read 2D real variables
    TYPE real_3d
       LOGICAL ::  from_file = .FALSE.  !< flag indicating whether an input variable is available and read from file or default values are used  
       INTEGER(iwp) ::  nz   !< number of grid points along vertical dimension
       REAL(wp) ::  fill = -9999.9_wp                  !< fill value
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  var !< respective variable
    END TYPE real_3d

    TYPE(grid_point), DIMENSION(:,:), ALLOCATABLE ::  grid  !< 2d Topography grid

    TYPE(mesh_point), DIMENSION(:), ALLOCATABLE, TARGET ::  mesh  !< Navigation mesh

    TYPE(real_2d) ::  terrain_height_f       !< input variable for terrain height

    TYPE(vertex_type) ::  dummy_vertex  !< placeholder vertex used for data copying
    TYPE(vertex_type) ::  null_vertex   !< placeholder vertex used for initialisation

    TYPE(vertex_type), DIMENSION(:), ALLOCATABLE ::  dummy_v_list  !< Dummy for reallocation of polygon array

    TYPE(polygon_type), POINTER ::  polygon  !< Current polygon

    TYPE(polygon_type), DIMENSION(:), ALLOCATABLE, TARGET ::  polygons  !< Building polygons

    TYPE(build_in) ::  buildings_f  !< input variable for buildings


 END MODULE variables
 
 MODULE data_input

    USE kinds
    USE variables

#if defined ( __netcdf )
    USE NETCDF
#endif

    INTERFACE get_variable
       MODULE PROCEDURE get_variable_1d_int
       MODULE PROCEDURE get_variable_1d_real
       MODULE PROCEDURE get_variable_2d_int8
       MODULE PROCEDURE get_variable_2d_int32
       MODULE PROCEDURE get_variable_2d_real
       MODULE PROCEDURE get_variable_3d_int8
       MODULE PROCEDURE get_variable_3d_real
       MODULE PROCEDURE get_variable_4d_real
    END INTERFACE get_variable

    INTERFACE get_attribute
       MODULE PROCEDURE get_attribute_real
       MODULE PROCEDURE get_attribute_int8
       MODULE PROCEDURE get_attribute_int32
       MODULE PROCEDURE get_attribute_string
    END INTERFACE get_attribute
    
    CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads orography and building information.
!------------------------------------------------------------------------------!
    SUBROUTINE netcdf_data_input_topo ( input_trunk )

       IMPLICIT NONE

       CHARACTER(LEN=*) ::  input_trunk       !< run path
       CHARACTER(LEN=200) ::  input_filename  !< filename

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< variable names in static input file


       INTEGER(iwp) ::  i            !< running index along x-direction
       INTEGER(iwp) ::  k_head       !< minimum k index for agents to walk underneath overhanging buildings
       INTEGER(iwp) ::  id_topo      !< NetCDF id of topograhy input file
       INTEGER(iwp) ::  j            !< running index along y-direction
       INTEGER(iwp) ::  num_vars     !< number of variables in netcdf input file

       LOGICAL ::  netcdf_flag = .FALSE.     !< indicates whether netcdf file is used for input
       LOGICAL ::  lod_flag = .FALSE.        !< true if 3d building data is used
       LOGICAL ::  topo_file_flag = .FALSE.  !< true if 3d building data is used

       WRITE(*,'((1X,A,/))') 'Looking for topography/building information'
       INQUIRE( FILE = TRIM( input_trunk )//'_static', EXIST = netcdf_flag )

       IF ( netcdf_flag ) THEN
          input_filename = TRIM( input_trunk )//'_static'
          WRITE(*,'(2(3X,A,/))') 'Topography/building data will be used from', & 
          TRIM( input_trunk )//'_static'
       ELSE
          WRITE(*,'(2(3X,A,/))') 'No static driver was found.',                &
                      'Trying to read building data from _topo file.'
          input_filename = TRIM( input_trunk )//'_topo'
          INQUIRE( FILE = TRIM( input_filename ), EXIST = topo_file_flag )
          IF ( .NOT. topo_file_flag ) THEN
             WRITE(*,'(6(3X,A,/))')                                            &
                 'No ASCII topography file was found in INPUT directory.',     &
                 'Make sure you provided building data in the form of either', &
                 '   (A) a static driver (<runname>_static) or',               &
                 '   (B) an ASCII topography file (<runname>_topo).',          &
                 NEW_LINE('A')//'Aborting nav_mesh program...'
             STOP
          ENDIF
          WRITE(*,'(2(3X,A,/))') 'Topography/building data will be used from', & 
          TRIM( input_trunk )//'_topo'
       ENDIF

!
!--    Input via palm-input data standard 
       IF ( netcdf_flag )  THEN
#if defined ( __netcdf )
!
!--       Open file in read-only mode
          CALL open_read_file( TRIM(input_filename) , id_topo )

!
!--       At first, inquire all variable names. 
!--       This will be used to check whether an input variable exist or not. 
          nc_stat = NF90_INQUIRE( id_topo, NVARIABLES = num_vars )
          CALL handle_error( 'inquire_num_variables', 534 )
!
!--       Allocate memory to store variable names and inquire them.
          ALLOCATE( var_names(1:num_vars) )
          CALL inquire_variable_names( id_topo, var_names )
!
!--       Terrain height. First, get variable-related _FillValue attribute
          IF ( check_existence( var_names, 'zt' ) )  THEN
             terrain_height_f%from_file = .TRUE. 
             CALL get_attribute( id_topo, char_fill,                           &
                                 terrain_height_f%fill,                        &
                                 .FALSE., 'zt' ) 
!
!--          PE-wise reading of 2D terrain height.
             ALLOCATE ( terrain_height_f%var(0:ny,0:nx)  )
             DO  i = 0, nx
                CALL get_variable( id_topo, 'zt',                    &
                                   i, terrain_height_f%var(:,i) )  
             ENDDO
          ELSE
             terrain_height_f%from_file = .FALSE. 
          ENDIF

!
!--       Read building height. First, read its _FillValue attribute, 
!--       as well as lod attribute
          buildings_f%from_file = .FALSE.
          IF ( check_existence( var_names, 'buildings_2d' ) )  THEN
             buildings_f%from_file = .TRUE. 
             CALL get_attribute( id_topo, char_lod, buildings_f%lod,           &
                                 .FALSE., 'buildings_2d' )

             CALL get_attribute( id_topo, char_fill,                           &
                                 buildings_f%fill1,                            &
                                 .FALSE., 'buildings_2d' )

!
!--             Read 2D topography
             IF ( buildings_f%lod == 1 )  THEN
                ALLOCATE ( buildings_f%var_2d(0:ny,0:nx) )
                DO  i = 0, nx
                   CALL get_variable( id_topo, 'buildings_2d',                 &
                                      i, buildings_f%var_2d(:,i) ) 
                ENDDO
             ELSE 
                WRITE(*,'(A)') 'NetCDF attribute lod ' //                      &
                                 '(level of detail) is not set properly.'
             ENDIF
          ENDIF
!
!--       If available, also read 3D building information. If both are
!--       available, use 3D information. Do this only if the flag that indicates
!--       that 2d buildings shall be used no matter what is false.
          IF ( check_existence( var_names, 'buildings_3d' )                    &
               .AND. .NOT. flag_2d )                                           &
          THEN
             lod_flag = .TRUE.
             buildings_f%from_file = .TRUE. 
             CALL get_attribute( id_topo, char_lod, buildings_f%lod,           &
                                 .FALSE., 'buildings_3d' )     

             CALL get_attribute( id_topo, char_fill,                           &
                                 buildings_f%fill2,                            &
                                 .FALSE., 'buildings_3d' )

             CALL get_dimension_length( id_topo, buildings_f%nz, 'z' )

             IF ( buildings_f%lod == 2 )  THEN
                ALLOCATE( buildings_f%z(0:buildings_f%nz-1) )
                CALL get_variable( id_topo, 'z', buildings_f%z )

                ALLOCATE( buildings_f%var_3d(0:buildings_f%nz-1,               &
                                             0:ny,0:nx) )
                buildings_f%var_3d = 0
!
!--                Read data PE-wise. Read yz-slices.
                DO  i = 0, nx
                   DO  j = 0, ny
                      CALL get_variable( id_topo, 'buildings_3d',              &
                                         i, j,                                 &
                                         buildings_f%var_3d(:,j,i) )
                   ENDDO
                ENDDO
             ELSE 
                WRITE(*,'(A)') 'NetCDF attribute lod ' //                      &
                                 '(level of detail) is not set properly.'
             ENDIF
          ENDIF

!
!--       Close topography input file
          CALL close_input_file( id_topo )
#endif
!
!--    ASCII input
       ELSE

          OPEN( 90, FILE= input_filename,                                      &
                STATUS='OLD', FORM='FORMATTED' )
!
!--             Read data from nyn to nys and nxl to nxr. Therefore, skip 
!--             column until nxl-1 is reached
          ALLOCATE ( buildings_f%var_2d(0:ny,0:nx) )
          DO  j = ny, 0, -1
             READ( 90, *, ERR=11, END=11 )                                     &
                             ( buildings_f%var_2d(j,i), i = 0, nx )
          ENDDO

          GOTO 12

11             WRITE(*,'(2A)') 'errors in file ',input_filename

12             CLOSE( 90 )
          buildings_f%from_file = .TRUE.

       ENDIF
!
!--    In case no terrain height is provided by static input file, allocate 
!--    array nevertheless and set terrain height to 0, which simplifies 
!--    topography initialization.
       IF ( .NOT. terrain_height_f%from_file )  THEN
          ALLOCATE ( terrain_height_f%var(0:ny,0:nx) )
          terrain_height_f%var = 0.0_wp
       ENDIF

!
!--    Transfer read data to uniform format: For agents the only relevant 
!--    information is whether they can walk or not at ground level.
       k_head = CEILING(2./dz)
       IF ( buildings_f%from_file ) THEN
          IF ( lod_flag ) THEN
             obstacle_height(0:nx,0:ny) = 1.
             DO j = 0, ny
                DO i = 0, nx
!
!--                For this purpose, an overhanging structure that an angent
!--                can walk beneath (e.g. a doorway) is not considered an
!--                obstacle.
                   IF ( ALL( buildings_f%var_3d(0:k_head,j,i) == 0 ) ) THEN
                      obstacle_height(i,j) = 0.
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             DO j = 0, ny
                DO i = 0, nx
                   obstacle_height(i,j) = buildings_f%var_2d(j,i)
                ENDDO
             ENDDO
          ENDIF
       ELSE
          WRITE(*,*) 'No building data was read from file. There will be no' //&
                     'navigation data available to agents.'
       ENDIF
       
    END SUBROUTINE netcdf_data_input_topo

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks if a given variables is on file
!------------------------------------------------------------------------------!
    FUNCTION check_existence( vars_in_file, var_name )

       IMPLICIT NONE

       CHARACTER(LEN=*) ::  var_name                   !< variable to be checked
       CHARACTER(LEN=*), DIMENSION(:) ::  vars_in_file !< list of variables in file

       INTEGER(iwp) ::  i                              !< loop variable

       LOGICAL ::  check_existence                     !< flag indicating whether a variable exist or not - actual return value

       i = 1
       check_existence = .FALSE.
       DO  WHILE ( i <= SIZE( vars_in_file ) )
          check_existence = TRIM( vars_in_file(i) ) == TRIM( var_name )  .OR.  &
                            check_existence
          i = i + 1
       ENDDO

       RETURN

    END FUNCTION check_existence


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Closes an existing netCDF file.
!------------------------------------------------------------------------------!
    SUBROUTINE close_input_file( id )
#if defined( __netcdf )

       IMPLICIT NONE

       INTEGER(iwp), INTENT(INOUT)        ::  id        !< file id

       nc_stat = NF90_CLOSE( id )
       CALL handle_error( 'close', 537 )
#endif
    END SUBROUTINE close_input_file

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Opens an existing netCDF file for reading only and returns its id.
!------------------------------------------------------------------------------!
    SUBROUTINE open_read_file( filename, id )
#if defined( __netcdf )

       IMPLICIT NONE

       CHARACTER (LEN=*), INTENT(IN) ::  filename  !< filename
       INTEGER(iwp), INTENT(INOUT)   ::  id        !< file id

       nc_stat = NF90_OPEN( filename, NF90_NOWRITE, id )

       CALL handle_error( 'open_read_file', 536 )

#endif
    END SUBROUTINE open_read_file



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Get dimension array for a given dimension
!------------------------------------------------------------------------------!
     SUBROUTINE get_dimension_length( id, dim_len, variable_name )
#if defined( __netcdf )

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name    !< dimension name 
       CHARACTER(LEN=100)          ::  dum              !< dummy variable to receive return character

       INTEGER(iwp)                ::  dim_len          !< dimension size
       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_dim           !< dimension id

!
!--    First, inquire dimension ID
       nc_stat = NF90_INQ_DIMID( id, TRIM( variable_name ), id_dim )
       CALL handle_error( 'get_dimension_length', 526 )
!
!--    Inquire dimension length
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim, dum, LEN = dim_len )
       CALL handle_error( 'get_dimension_length', 526 ) 

#endif
    END SUBROUTINE get_dimension_length

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 1D integer variable from file. 
!------------------------------------------------------------------------------!
     SUBROUTINE get_variable_1d_int( id, variable_name, var )

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name    !< variable name 

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< dimension id

       INTEGER(iwp), DIMENSION(:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
       CALL handle_error( 'get_variable_1d_int', 527 )
!
!--    Inquire dimension length
       nc_stat = NF90_GET_VAR( id, id_var, var )
       CALL handle_error( 'get_variable_1d_int', 527 ) 

#endif
    END SUBROUTINE get_variable_1d_int

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 1D float variable from file. 
!------------------------------------------------------------------------------!
     SUBROUTINE get_variable_1d_real( id, variable_name, var )

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  variable_name    !< variable name 

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< dimension id

       REAL(wp), DIMENSION(:), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    First, inquire variable ID
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
       CALL handle_error( 'get_variable_1d_real', 527 )
!
!--    Inquire dimension length
       nc_stat = NF90_GET_VAR( id, id_var, var )
       CALL handle_error( 'get_variable_1d_real', 527 ) 

#endif
    END SUBROUTINE get_variable_1d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 2D REAL variable from a file. Reading is done processor-wise, 
!> i.e. each core reads its own domain in slices along x. 
!------------------------------------------------------------------------------! 
    SUBROUTINE get_variable_2d_real( id, variable_name, i, var )

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp), INTENT(IN)      ::  i               !< index along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id

       REAL(wp), DIMENSION(0:ny), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, var(0:ny),                          &
                               start = (/ i+1, 1 /),                           &
                               count = (/ 1, ny + 1 /) )

       CALL handle_error( 'get_variable_2d_real', 528 )
#endif
    END SUBROUTINE get_variable_2d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 2D 32-bit INTEGER variable from file. Reading is done processor-wise, 
!> i.e. each core reads its own domain in slices along x. 
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_2d_int32( id, variable_name, i, var )

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp), INTENT(IN)      ::  i               !< index along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp), DIMENSION(0:ny), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, var(0:ny),                          &
                               start = (/ i+1, 1 /),                           &
                               count = (/ 1, ny + 1 /) )

       CALL handle_error( 'get_variable_2d_int32', 529 )
#endif
    END SUBROUTINE get_variable_2d_int32

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 2D 8-bit INTEGER variable from file. Reading is done processor-wise, 
!> i.e. each core reads its own domain in slices along x. 
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_2d_int8( id, variable_name, i, var )

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp), INTENT(IN)      ::  i               !< index along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(KIND=1), DIMENSION(0:ny), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )
!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, var(0:ny),                          &
                               start = (/ i+1, 1 /),                           &
                               count = (/ 1, ny + 1 /) )

       CALL handle_error( 'get_variable_2d_int8', 530 )
#endif
    END SUBROUTINE get_variable_2d_int8

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 3D 8-bit INTEGER variable from file.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_3d_int8( id, variable_name, i, j, var )

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp), INTENT(IN)      ::  i               !< index along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp), INTENT(IN)      ::  j               !< index along y direction
       INTEGER(iwp)                  ::  n_file          !< number of data-points along 3rd dimension

       INTEGER(iwp), DIMENSION(1:3)  ::  id_dim

       INTEGER( KIND = 1 ), DIMENSION(0:buildings_f%nz-1), INTENT(INOUT) ::  var  !< variable to be read
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Get length of first dimension, required for the count parameter.
!--    Therefore, first inquired dimension ids
       nc_stat = NF90_INQUIRE_VARIABLE( id, id_var, DIMIDS = id_dim )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim(3), LEN = n_file )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, var,                                &
                               start = (/ i+1, j+1, 1 /),                      &
                               count = (/ 1, 1, n_file /) )

       CALL handle_error( 'get_variable_3d_int8', 531 )
#endif
    END SUBROUTINE get_variable_3d_int8


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 3D float variable from file.  
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_3d_real( id, variable_name, i, j, var )

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp), INTENT(IN)      ::  i               !< index along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp), INTENT(IN)      ::  j               !< index along y direction
       INTEGER(iwp)                  ::  n3              !< number of data-points along 3rd dimension

       INTEGER(iwp), DIMENSION(3)    ::  id_dim

       REAL(wp), DIMENSION(:), INTENT(INOUT) ::  var     !< variable to be read
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Get length of first dimension, required for the count parameter.
!--    Therefore, first inquired dimension ids
       nc_stat = NF90_INQUIRE_VARIABLE( id, id_var, DIMIDS = id_dim )
       nc_stat = NF90_INQUIRE_DIMENSION( id, id_dim(3), LEN = n3 )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, var,                                &
                               start = (/ i+1, j+1, 1 /),                      &
                               count = (/ 1, 1, n3 /) )

       CALL handle_error( 'get_variable_3d_real', 532 )
#endif
    END SUBROUTINE get_variable_3d_real


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads a 4D float variable from file. Note, in constrast to 3D versions, 
!> dimensions are already inquired and passed so that they are known here.
!------------------------------------------------------------------------------!
    SUBROUTINE get_variable_4d_real( id, variable_name, i, j, var, n3, n4 )

       IMPLICIT NONE

       CHARACTER(LEN=*)              ::  variable_name   !< variable name

       INTEGER(iwp), INTENT(IN)      ::  i               !< index along x direction
       INTEGER(iwp), INTENT(IN)      ::  id              !< file id
       INTEGER(iwp)                  ::  id_var          !< variable id
       INTEGER(iwp), INTENT(IN)      ::  j               !< index along y direction
       INTEGER(iwp), INTENT(IN)      ::  n3              !< number of data-points along 3rd dimension
       INTEGER(iwp), INTENT(IN)      ::  n4              !< number of data-points along 4th dimension

       REAL(wp), DIMENSION(:,:), INTENT(INOUT) ::  var     !< variable to be read
#if defined( __netcdf )

!
!--    Inquire variable id
       nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
!
!--    Get variable
       nc_stat = NF90_GET_VAR( id, id_var, var,                                &
                               start = (/ i+1, j+1, 1, 1 /),                   &
                               count = (/ 1, 1, n3, n4 /) )

       CALL handle_error( 'get_variable_4d_real', 533 )
#endif
    END SUBROUTINE get_variable_4d_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prints out a text message corresponding to the current status.
!------------------------------------------------------------------------------!
    SUBROUTINE handle_error( routine_name, errno )

       IMPLICIT NONE

       CHARACTER(LEN=6) ::  message_identifier
       CHARACTER(LEN=*) ::  routine_name
       CHARACTER(LEN=100) ::  message_string

       INTEGER(iwp) ::  errno
#if defined( __netcdf )

       IF ( nc_stat /= NF90_NOERR )  THEN

          WRITE( message_identifier, '(''NC'',I4.4)' )  errno
          
          message_string = TRIM( NF90_STRERROR( nc_stat ) )

          WRITE(*,*) routine_name,'  ', message_identifier,'  ', TRIM(message_string)
          WRITE(*,*) 'Aborting NavMesh-tool'

       ENDIF

#endif
    END SUBROUTINE handle_error


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inquires the variable names belonging to a file. 
!------------------------------------------------------------------------------!
    SUBROUTINE inquire_variable_names( id, var_names )

       IMPLICIT NONE

       CHARACTER(LEN=*), DIMENSION(:), INTENT(INOUT) ::  var_names   !< return variable - variable names
       INTEGER(iwp)                                  ::  i           !< loop variable
       INTEGER(iwp), INTENT(IN)                      ::  id          !< file id
       INTEGER(iwp)                                  ::  num_vars    !< number of variables (unused return parameter)
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE       ::  varids      !< dummy array to strore variable ids temporarily
#if defined( __netcdf )

       ALLOCATE( varids(1:SIZE(var_names)) )
       nc_stat = NF90_INQ_VARIDS( id, NVARS = num_vars, VARIDS = varids )
       CALL handle_error( 'inquire_variable_names', 535 )

       DO  i = 1, SIZE(var_names)
          nc_stat = NF90_INQUIRE_VARIABLE( id, varids(i), NAME = var_names(i) )
          CALL handle_error( 'inquire_variable_names', 535 )
       ENDDO

       DEALLOCATE( varids )
#endif
    END SUBROUTINE inquire_variable_names

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type INTEGER (32-bit)
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_int32( id, attribute_name, value, global,        &
                                     variable_name )

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL  ::  variable_name    !< variable name 

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id
       INTEGER(iwp), INTENT(INOUT) ::  value            !< read value

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute
#if defined( __netcdf )

!
!--    Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int32 global', 522 )
!
!--    Read attributes referring to a single variable. Therefore, first inquire 
!--    variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_int32', 522 )
          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int32', 522 )       
       ENDIF
#endif
    END SUBROUTINE get_attribute_int32

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type INTEGER (8-bit)
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_int8( id, attribute_name, value, global,         &
                                    variable_name )

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL  ::  variable_name    !< variable name 

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id
       INTEGER(KIND=1), INTENT(INOUT) ::  value         !< read value

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute
#if defined( __netcdf )

!
!--    Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int8 global', 523 )
!
!--    Read attributes referring to a single variable. Therefore, first inquire 
!--    variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_int8', 523 )
          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_int8', 523 )       
       ENDIF
#endif
    END SUBROUTINE get_attribute_int8

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type REAL
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_real( id, attribute_name, value, global,         &
                                    variable_name )

       IMPLICIT NONE

       CHARACTER(LEN=*)            ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL  ::  variable_name    !< variable name 

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute

       REAL(wp), INTENT(INOUT)     ::  value            !< read value
#if defined( __netcdf )


!
!-- Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_real global', 524 )
!
!-- Read attributes referring to a single variable. Therefore, first inquire 
!-- variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_real', 524 )
          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_real', 524 )       
       ENDIF
#endif
    END SUBROUTINE get_attribute_real

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads global or variable-related attributes of type CHARACTER
!> Remark: reading attributes of type NF_STRING return an error code 56 -
!> Attempt to convert between text & numbers.
!------------------------------------------------------------------------------!
     SUBROUTINE get_attribute_string( id, attribute_name, value, global,       &
                                      variable_name )

       IMPLICIT NONE

       CHARACTER(LEN=*)                ::  attribute_name   !< attribute name
       CHARACTER(LEN=*), OPTIONAL      ::  variable_name    !< variable name 
       CHARACTER(LEN=*), INTENT(INOUT) ::  value            !< read value

       INTEGER(iwp), INTENT(IN)    ::  id               !< file id
       INTEGER(iwp)                ::  id_var           !< variable id

       LOGICAL, INTENT(IN) ::  global                   !< flag indicating global attribute
#if defined( __netcdf )

!
!--    Read global attribute
       IF ( global )  THEN
          nc_stat = NF90_GET_ATT( id, NF90_GLOBAL, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_string global', 525 )
!
!--    Read attributes referring to a single variable. Therefore, first inquire 
!--    variable id
       ELSE
          nc_stat = NF90_INQ_VARID( id, TRIM( variable_name ), id_var )
          CALL handle_error( 'get_attribute_string', 525 )

          nc_stat = NF90_GET_ATT( id, id_var, TRIM( attribute_name ), value )
          CALL handle_error( 'get_attribute_string',525 )
       ENDIF
#endif
    END SUBROUTINE get_attribute_string

 END MODULE

 MODULE mod_functions

    USE kinds

    PRIVATE
    PUBLIC dist_point_to_edge, intersect, is_left, is_right

    CONTAINS

!
!-- Calculates distance of point P to edge (A,B). If A = B, calculates
!-- point-to-point distance from A/B to P
    FUNCTION dist_point_to_edge ( a_x, a_y, b_x, b_y, p_x, p_y )

       IMPLICIT NONE

       REAL(wp)  :: ab_x                !< x-coordinate of vector from A to B
       REAL(wp)  :: ab_y                !< y-coordinate of vector from A to B
       REAL(wp)  :: ab_d                !< inverse length of vector from A to B
       REAL(wp)  :: ab_u_x              !< x-coordinate of vector with direction of ab and length 1
       REAL(wp)  :: ab_u_y              !< y-coordinate of vector with direction of ab and length 1
       REAL(wp)  :: ba_x                !< x-coordinate of vector from B to A
       REAL(wp)  :: ba_y                !< y-coordinate of vector from B to A
       REAL(wp)  :: ap_x                !< x-coordinate of vector from A to P
       REAL(wp)  :: ap_y                !< y-coordinate of vector from A to P
       REAL(wp)  :: bp_x                !< x-coordinate of vector from B to P
       REAL(wp)  :: bp_y                !< y-coordinate of vector from B to P
       REAL(wp)  :: a_x                 !< x-coordinate of point A of edge
       REAL(wp)  :: a_y                 !< y-coordinate of point A of edge
       REAL(wp)  :: b_x                 !< x-coordinate of point B of edge
       REAL(wp)  :: b_y                 !< y-coordinate of point B of edge
       REAL(wp)  :: p_x                 !< x-coordinate of point P
       REAL(wp)  :: p_y                 !< y-coordinate of point P
       REAL(wp)  :: dist_x              !< x-coordinate of point P
       REAL(wp)  :: dist_y              !< y-coordinate of point P
       REAL(wp)  :: dist_point_to_edge  !< y-coordinate of point P

       ab_x = - a_x + b_x
       ab_y = - a_y + b_y
       ba_x = - b_x + a_x 
       ba_y = - b_y + a_y 
       ap_x = - a_x + p_x
       ap_y = - a_y + p_y
       bp_x = - b_x + p_x
       bp_y = - b_y + p_y

       IF ( ab_x * ap_x + ab_y * ap_y <= 0. ) THEN
          dist_point_to_edge = SQRT((a_x - p_x)**2 + (a_y - p_y)**2)
       ELSEIF ( ba_x * bp_x + ba_y * bp_y <= 0. ) THEN
          dist_point_to_edge = SQRT((b_x - p_x)**2 + (b_y - p_y)**2)
       ELSE
          ab_d = 1./SQRT((ab_x)**2+(ab_y)**2)
          ab_u_x = ab_x*ab_d
          ab_u_y = ab_y*ab_d
          dist_x = ap_x - (ap_x*ab_u_x+ap_y*ab_u_y)*ab_u_x
          dist_y = ap_y - (ap_x*ab_u_x+ap_y*ab_u_y)*ab_u_y
          dist_point_to_edge = SQRT( dist_x**2 + dist_y**2 )
       ENDIF

       RETURN

    END FUNCTION dist_point_to_edge

!
!-- Returns true if the line segments AB and PQ share an intersection
    FUNCTION intersect ( ax, ay, bx, by, px, py, qx, qy )

       IMPLICIT NONE

       LOGICAL  :: intersect !< return value; TRUE if intersection was found
       LOGICAL  :: la        !< T if a is left of PQ
       LOGICAL  :: lb        !< T if b is left of PQ
       LOGICAL  :: lp        !< T if p is left of AB
       LOGICAL  :: lq        !< T if q is left of AB
       LOGICAL  :: poss      !< flag that indicates if an intersection is still possible
       LOGICAL  :: ra        !< T if a is right of PQ
       LOGICAL  :: rb        !< T if b is right of PQ
       LOGICAL  :: rp        !< T if p is right of AB
       LOGICAL  :: rq        !< T if q is right of AB

       REAL(wp)  :: ax     !< x-coordinate of point A
       REAL(wp)  :: ay     !< y-coordinate of point A
       REAL(wp)  :: bx     !< x-coordinate of point B
       REAL(wp)  :: by     !< y-coordinate of point B
       REAL(wp)  :: px     !< x-coordinate of point P
       REAL(wp)  :: py     !< y-coordinate of point P
       REAL(wp)  :: qx     !< x-coordinate of point Q
       REAL(wp)  :: qy     !< y-coordinate of point Q

       intersect = .FALSE.
       poss      = .FALSE.
!
!--    Intersection is possible only if P and Q are on opposing sides of AB
       lp = is_left(ax,ay,bx,by,px,py)
       rq = is_right(ax,ay,bx,by,qx,qy)
       IF ( lp .AND. rq ) poss = .TRUE.
       IF ( .NOT. poss ) THEN
          lq = is_left(ax,ay,bx,by,qx,qy)
          rp = is_right(ax,ay,bx,by,px,py)
          IF ( lq .AND. rp ) poss = .TRUE.
       ENDIF
!
!--    Intersection occurs only if above test (poss) was true AND
!--    A and B are on opposing sides of PQ
       IF ( poss ) THEN
          la = is_left(px,py,qx,qy,ax,ay)
          rb = is_right(px,py,qx,qy,bx,by)
          IF ( la .AND. rb ) intersect = .TRUE.
          IF ( .NOT. intersect ) THEN
             lb = is_left(px,py,qx,qy,bx,by)
             ra = is_right(px,py,qx,qy,ax,ay)
             IF ( lb .AND. ra ) intersect = .TRUE.
          ENDIF
       ENDIF

       RETURN

    END FUNCTION intersect 

!
!-- Calculates if point P is left of the infinite
!-- line that contains A and B (direction: A to B) 
!-- Concept: 2D rotation of two vectors
    FUNCTION is_left ( ax, ay, bx, by, px, py )

       IMPLICIT NONE

       LOGICAL  :: is_left !< return value; TRUE if P is left of AB

       REAL(wp)  :: ax     !< x-coordinate of point A
       REAL(wp)  :: ay     !< y-coordinate of point A
       REAL(wp)  :: bx     !< x-coordinate of point B
       REAL(wp)  :: by     !< y-coordinate of point B
       REAL(wp)  :: px     !< x-coordinate of point P
       REAL(wp)  :: py     !< y-coordinate of point P
!
!--    2D-rotation
       is_left = (bx-ax)*(py-ay)-(px-ax)*(by-ay) > 0
!
!--    False if the point is on the line (or very close)
       IF ( (ABS(ax-px) < .001 .AND. ABS(ay-py) < .001) .OR.                  &
            (ABS(bx-px) < .001 .AND. ABS(by-py) < .001) )                     &
       THEN
          is_left = .FALSE.
       ENDIF

       RETURN

    END FUNCTION is_left 

!
!-- Calculates if point P is right of the infinite
!-- line that contains A and B (direction: A to B) 
!-- Concept: 2D rotation of two vectors
    FUNCTION is_right ( ax, ay, bx, by, px, py )

       IMPLICIT NONE

       LOGICAL  :: is_right !< return value; TRUE if P is right of AB

       REAL(wp), INTENT(IN)  :: ax     !< x-coordinate of point A
       REAL(wp), INTENT(IN)  :: ay     !< y-coordinate of point A
       REAL(wp), INTENT(IN)  :: bx     !< x-coordinate of point B
       REAL(wp), INTENT(IN)  :: by     !< y-coordinate of point B
       REAL(wp), INTENT(IN)  :: px     !< x-coordinate of point P
       REAL(wp), INTENT(IN)  :: py     !< y-coordinate of point P

!
!--    2D-rotation
       is_right = (bx-ax)*(py-ay)-(px-ax)*(by-ay) < 0
!
!--    False if the point is on the line (or very close)
       IF ( (ABS(ax-px) < .001 .AND. ABS(ay-py) < .001) .OR.                  &
            (ABS(bx-px) < .001 .AND. ABS(by-py) < .001) )                     &
       THEN
          is_right = .FALSE.
       ENDIF

       RETURN

    END FUNCTION is_right 

 END MODULE mod_functions
 
 MODULE polygon_creation

    USE kinds

    USE mod_functions

    USE variables

    USE data_input

    CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialisation, allocation, and reading of some input
!------------------------------------------------------------------------------!
    SUBROUTINE init

       IMPLICIT NONE

       CHARACTER(LEN=255) ::  dirname      !< dummy to read current working directory
       CHARACTER(LEN=255) ::  rundir       !< base run directory
       CHARACTER(LEN=255) ::  input_trunk  !< base filename for run including path
       CHARACTER(LEN=80) ::  line          !< string to identify namelist
       CHARACTER(LEN=80) ::  line_dum      !< line dummy for error output

       CHARACTER(LEN=2),DIMENSION(1:5) ::  run_pars  !< parameters from other namelist

       INTEGER(iwp) ::  ie            !< end index (string manipulation)
       INTEGER(iwp) ::  is            !< start index (string manipulation)
       INTEGER(iwp) ::  line_counter  !< line on which reading error occured

       LOGICAL ::  p3d_flag = .FALSE.  !< indicates whether p3d file was found

       NAMELIST /prepro_par/  flag_2d, internal_buildings, tolerance_dp

       WRITE(*,'(1X,A)')                                                        &
                 "o----------------------------------------------o",           &
                 "| o------------------------------------------o |",           &
                 "| |         __   ____  ____       ____       | |",           &
                 "| |        / _\ (  _ \(_  _) ___ (  _ \      | |",           &
                 "| |       /    \ ) __/  )(  (___) ) __/      | |",           &
                 "| |       \_/\_/(__)   (__)      (__)        | |",           &
                 "| |                                          | |",           &
                 "| |    Agent Preprocessing Tool for PALM     | |",           &
                 "| o------------------------------------------o |",           &
                 "o----------------------------------------------o"

!
!--    Identify run name and Input files
       CALL GET_ENVIRONMENT_VARIABLE('PWD', dirname)
       ie = INDEX(dirname, '/', BACK=.TRUE.)
       is = INDEX(dirname(1:ie-1), '/', BACK=.TRUE.)
       IF ( TRIM(ADJUSTL(dirname(ie+1:))) /= 'INPUT' ) THEN
          WRITE(*,'(3X,A)') 'NavMesh was called from',                         &
                            ' ', TRIM(ADJUSTL(dirname)), ' ',                  &
                            'and is now aborting. Please call this tool',      &
                            'from the INPUT-folder of your job directory.'
          STOP
       ENDIF
       runname = TRIM(ADJUSTL(dirname(is+1:ie-1)))
       rundir  = TRIM(ADJUSTL(dirname(1:ie)))
       input_trunk = TRIM(rundir)//'INPUT/'//TRIM(runname)

!
!--    Check for parameter file
       INQUIRE( FILE = TRIM( input_trunk )//'_p3d', EXIST = p3d_flag )
       IF ( .NOT. p3d_flag ) THEN
          WRITE(*,'(3(3X,A,/))') 'No _p3d file was found. Aborting.',          &
                                 'I was looking for the file',                 &
                                 TRIM( input_trunk )//'_p3d'
          STOP
       ELSE
          WRITE(*,'(2(3X,A,/))') 'The following input file will be used:',     & 
                                  TRIM( input_trunk )//'_p3d'
       ENDIF

!
!--    Read run parameters from run parameter file (_p3d), though not from
!--    namelist.
       run_pars = (/'dx','dy','dz','nx','ny'/)
       OPEN ( 11, FILE=TRIM(input_trunk)//'_p3d', FORM='FORMATTED',    &
                     STATUS='OLD' )
       DO i = 1, SIZE(run_pars)
          REWIND ( 11 )
          line = ' '
          DO   WHILE ( INDEX( line, run_pars(i) ) == 0 )
             READ ( 11, '(A)', END=10 )  line
             IF ( INDEX(line, '!') /= 0 ) THEN
                IF ( INDEX(line, run_pars(i)) > INDEX(line, '!' ) ) THEN
                   line = ' '
                   CYCLE
                ENDIF
             ENDIF
          ENDDO
          line = TRIM(ADJUSTL(line(INDEX(line,'=')+1:INDEX(line,',')-1)))
          SELECT CASE (i)
             CASE(1)
                READ(line,*) dx
             CASE(2) 
                READ(line,*) dy
             CASE(3) 
                READ(line,*) dz
             CASE(4) 
                READ(line,*) nx
             CASE(5) 
                READ(line,*) ny
             CASE DEFAULT
          END SELECT
       ENDDO
10     CONTINUE

!
!--    Try to find prepro package
       REWIND ( 11 )
       line = ' '
       DO   WHILE ( INDEX( line, '&prepro_par' ) == 0 )
          READ ( 11, '(A)', END=40 )  line
       ENDDO
       BACKSPACE ( 11 )

!
!--    Read user-defined namelist
       READ ( 11, prepro_par, ERR = 20, END = 40 )
       GOTO 40

 20    BACKSPACE( 11 )
       READ( 11 , '(A)') line

       line_dum = ' '
       line_counter = 0

       REWIND( 11 )
       DO WHILE ( INDEX( line_dum, TRIM(line) ) == 0 )
          READ ( 11, '(A)', END=30 )  line_dum
          line_counter = line_counter + 1
       ENDDO

 30    WRITE( *, '(A,/,A,I3,A,/,A)' ) 'Error(s) in NAMELIST prepro_par.',      &
                                      'Reading fails on line ', line_counter,  &
                                      ' at ', TRIM(ADJUSTL(line))
       STOP

 40    CONTINUE
       CLOSE( 11 )

!
!--    If tolerance_dp was not set, put in default values
       DO i = 0, 2
          IF ( tolerance_dp(i) > 999998.0_wp ) THEN
             tolerance_dp(i) = SQRT(dx*dy)*1.41/(2**i)
          ELSE
             tolerance_dp(i) = tolerance_dp(i)*SQRT(dx*dy)
          ENDIF
       ENDDO

!
!--    Allocate arrays
       ALLOCATE(obstacle_height(-3:nx+3,-3:ny+3), polygon_id(-3:nx+3,-3:ny+3), &
                wall_flags_0(-3:nx+3,-3:ny+3), grid(-3:nx+3,-3:ny+3))
!
!--    Set null_vertex
       CALL set_vertex(null_vertex,0.0_wp,0.0_wp)
!
!--    Some initializations
       ddx = 1./dx
       ddy = 1./dy

       polygon_id                = 0
       obstacle_height           = 0.

       grid%checked              = .FALSE.
       grid(-3:-1,:)%checked     = .TRUE.
       grid(nx+1:nx+3,:)%checked = .TRUE.
       grid(:,-3:-1)%checked     = .TRUE.
       grid(:,ny+1:ny+3)%checked = .TRUE.
       grid%polygon_id           = 0
!
!--    Open files and topography/building data
       CALL netcdf_data_input_topo ( TRIM(input_trunk) )

    END SUBROUTINE init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Identifies all grid boxes that belong to a building and assigns a building
!> number to each grid box.
!> Method: Scans each grid point. If a grid point was not previously checked
!>         and contains a building, it is added to a new polygon and marked as
!>         checked. Then, all its neighbors that also contain a building are
!>         added to the same polygon and marked as checked. All neighbors of
!>         neighbors are subsequently found, added and checked until none are
!>         left. Then, the initial scan continues, skipping already checked
!>         grid points. Once a grid point with a new building is found, the
!>         polygon_id increases and the grid point and all others that belong
!>         to the same building are added as described above.
!> NOTE:   This procedure will identify grid points that are only connected
!>         diagonally (share only one point with each other) as connected and
!>         have them belonging to the same building. This is necessary, as an
!>         agent will not be able to traverse between these grid points and the
!>         navigation mesh will therefore have to make him circumvent this
!>         point.
!------------------------------------------------------------------------------!
    SUBROUTINE identify_polygons

       IMPLICIT NONE

       INTEGER(iwp) ::  ii    !< local counter
       INTEGER(iwp) ::  il    !< local counter
       INTEGER(iwp) ::  jj    !< local counter
       INTEGER(iwp) ::  jl    !< local counter
       INTEGER(iwp) ::  gpil  !< number of grid points in list
       INTEGER(iwp) ::  gpta  !< number of grid points to add to grid_list

       TYPE(grid_point), DIMENSION(1:7) ::  add_to_grid_list !< grid points to be added to the list

       TYPE(grid_point), DIMENSION(:), ALLOCATABLE ::  dummy_grid_list !< dummy for reallocation of grid_list
       TYPE(grid_point), DIMENSION(:), ALLOCATABLE ::  grid_list       !< list of grid points that belong to the current building but whose neighbors have not been checked yet

!
!--    Initialize wall_flags array: 1 where no buildings, 0 where buildings
       wall_flags_0 = 1
       DO i = 0, nx
          DO j = 0, ny
             IF ( obstacle_height(i,j) > 0 ) THEN
                wall_flags_0(i,j) = 0
             ENDIF
          ENDDO
       ENDDO
       DEALLOCATE(obstacle_height)
       polygon_counter = 0
       gpil = 0
       gpta = 0
       ALLOCATE(grid_list(1:100))
!
!--    Declare all grid points that contain no buildings as already checked.
!--    This way, these points will be skipped in the following calculation and
!--    will have polygon_id = 0
       DO i = 0, nx
          DO j = 0, ny
             IF ( BTEST( wall_flags_0(i,j), 0 ) ) THEN
                grid(i,j)%checked = .TRUE.
             ENDIF
          ENDDO
       ENDDO
!
!--    Check all grid points and process them
       DO i = 0, nx
          DO j = 0, ny
!
!--          If the current grid point has not been checked, mark it as checked.
!--          As it is the first point belonging to a new building, increase the
!--          polygon_id counter and associate the grid point with that id.
             IF ( .NOT. grid(i,j)%checked ) THEN
                polygon_counter = polygon_counter + 1
                grid(i,j)%polygon_id = polygon_counter
                grid(i,j)%checked = .TRUE.
!
!--             Check if any neighbors of the found grid point are part of a 
!--             building too. If so, add them to the list of grid points
!--             that have to be checked and associate them with the same polygon
                gpta = 0
                DO ii = i-1, i+1
                   DO jj = j-1, j+1
                      IF ( ii == i .AND. jj == j ) CYCLE
                      IF ( .NOT. grid(ii,jj)%checked ) THEN
                         gpta = gpta + 1
                         add_to_grid_list(gpta)%i = ii
                         add_to_grid_list(gpta)%j = jj
                      ENDIF
                   ENDDO
                ENDDO

!
!--             Change size of grid_list if it becomes too small
                IF ( gpil + gpta > SIZE(grid_list) ) THEN
                   ALLOCATE(dummy_grid_list(1:gpil))
                   dummy_grid_list = grid_list(1:gpil)
                   DEALLOCATE(grid_list)
                   ALLOCATE(grid_list(1:2*(gpil+gpta)))
                   grid_list(1:gpil) = dummy_grid_list(1:gpil)
                   DEALLOCATE(dummy_grid_list)
                ENDIF
!
!--             If there are grid points to add to grid_list, add them
                IF ( gpta > 0 ) THEN
                   grid_list(gpil+1:gpil+gpta) = add_to_grid_list(1:gpta)
                   gpil = gpil + gpta
                ENDIF
!
!--             Handle all grid points in grid_list until there are none left
                DO WHILE (gpil>0)
                   il = grid_list(gpil)%i
                   jl = grid_list(gpil)%j
                   grid(il,jl)%polygon_id = polygon_counter
                   grid(il,jl)%checked = .TRUE.
!
!--                this grid point in the list is processed, so the number of
!--                grid points in the list can be reduced by one
                   gpil = gpil - 1
                   gpta = 0
!
!--                For the current grid point, check if any unchecked
!--                neighboring grid points also contain a building. All such
!--                grid points are added to the list of grid points to be
!--                handled in this loop
                   DO ii = il-1, il+1
                      DO jj = jl-1, jl+1
                         IF ( jj == jl .AND. ii == il ) CYCLE
                         IF ( .NOT. grid(ii,jj)%checked )                      &
                         THEN
                            gpta = gpta + 1
                            add_to_grid_list(gpta)%i = ii
                            add_to_grid_list(gpta)%j = jj
                         ENDIF
                      ENDDO
                   ENDDO
!
!--                Change size of grid list if it becomes too small
                   IF ( gpil + gpta > SIZE(grid_list) ) THEN
                      ALLOCATE(dummy_grid_list(1:gpil))
                      dummy_grid_list = grid_list(1:gpil)
                      DEALLOCATE(grid_list)
                      ALLOCATE(grid_list(1:2*(gpil+gpta)))
                      grid_list(1:gpil) = dummy_grid_list(1:gpil)
                      DEALLOCATE(dummy_grid_list)
                   ENDIF
!
!--                If there are grid points to add to list, add them
                   IF ( gpta > 0 ) THEN
                      grid_list(gpil+1:gpil+gpta) = add_to_grid_list(1:gpta)
                      gpil = gpil + gpta
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       DEALLOCATE(grid_list)
!
!--    Set size of polygon array and initialize
       ALLOCATE(polygons(1:polygon_counter))
       polygons%nov = 0

    END SUBROUTINE identify_polygons

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Identifies the corners of the PALM building topography and adds them to
!> a specific polygon for each building as vertices. This converts the gridded
!> building data into one polygon per building that contains the coordinates of
!> each inner and outer corner of that building as vertices.
!> A grid point contains an outer corner if it's part of a building and exactly
!>   one of its horizontal and one of its vertical neighbors is also part of a
!>   building (4 cases).
!> A grid point contains an inner corner if it's not part of a building and
!>   exactly one of its horizontal, one of its diagonal and one of its vertical
!>   neighbors are each part of a building and in turn neighbors
!>   to each other (4 cases).
!------------------------------------------------------------------------------!
    SUBROUTINE identify_corners

       IMPLICIT NONE

       INTEGER(iwp) ::  p_id  !< current polygon_id
!
!--    For all grid points, check whether it contains one or more corners
       DO i = 0, nx
          DO j = 0, ny
!
!--          First, check if grid contains topography and has a corner.
             IF ( .NOT. BTEST( wall_flags_0(i,j), 0 ) ) THEN
!
!--             Corner at south left edge of grid cell
                IF ( BTEST( wall_flags_0(i-1,j), 0 ) .AND.                     &
                     BTEST( wall_flags_0(i,j-1), 0 ))                          &
                THEN
                   p_id = grid(i,j)%polygon_id
                   polygons(p_id)%nov = polygons(p_id)%nov + 1
                   nov = polygons(p_id)%nov
                   CALL set_vertex(dummy_vertex, i*dx, j*dy)
                   CALL add_vertex_to_polygon(dummy_vertex, p_id, nov)
                ENDIF
!
!--             Corner at north left edge of grid cell
                IF ( BTEST( wall_flags_0(i-1,j), 0 ) .AND.                     &
                     BTEST( wall_flags_0(i,j+1), 0 ))                          &
                THEN
                   p_id = grid(i,j)%polygon_id
                   polygons(p_id)%nov = polygons(p_id)%nov + 1
                   nov = polygons(p_id)%nov
                   CALL set_vertex(dummy_vertex, i*dx, (j+1)*dy)
                   CALL add_vertex_to_polygon(dummy_vertex, p_id, nov)
                ENDIF
!
!--             Corner at north right edge of grid cell
                IF ( BTEST( wall_flags_0(i+1,j), 0 ) .AND.                     &
                     BTEST( wall_flags_0(i,j+1), 0 ))                          &
                THEN
                   p_id = grid(i,j)%polygon_id
                   polygons(p_id)%nov = polygons(p_id)%nov + 1
                   nov = polygons(p_id)%nov
                   CALL set_vertex(dummy_vertex, (i+1)*dx, (j+1)*dy)
                   CALL add_vertex_to_polygon(dummy_vertex, p_id, nov)
                ENDIF
!
!--             Corner at south right edge of grid cell
                IF ( BTEST( wall_flags_0(i+1,j), 0 ) .AND.                     &
                     BTEST( wall_flags_0(i,j-1), 0 ))                          &
                THEN
                   p_id = grid(i,j)%polygon_id
                   polygons(p_id)%nov = polygons(p_id)%nov + 1
                   nov = polygons(p_id)%nov
                   CALL set_vertex(dummy_vertex, (i+1)*dx, j*dy)
                   CALL add_vertex_to_polygon(dummy_vertex, p_id, nov)
                ENDIF
!
!--          Second, check if grid contains no topography and has a corner.
             ELSE
!
!--             Corner at south left edge of grid cell
                IF ( .NOT. BTEST( wall_flags_0(i-1,j), 0 ) .AND.               &
                     .NOT. BTEST( wall_flags_0(i,j-1), 0 ) .AND.               &
                     .NOT. BTEST( wall_flags_0(i-1,j-1), 0 ) )                 &
                THEN
                   p_id = grid(i-1,j-1)%polygon_id
                   polygons(p_id)%nov = polygons(p_id)%nov + 1
                   nov = polygons(p_id)%nov
                   CALL set_vertex(dummy_vertex, i*dx, j*dy)
                   CALL add_vertex_to_polygon(dummy_vertex, p_id, nov)
                ENDIF
!
!--             Corner at north left edge of grid cell
                IF ( .NOT. BTEST( wall_flags_0(i-1,j), 0 ) .AND.               &
                     .NOT. BTEST( wall_flags_0(i,j+1), 0 ) .AND.               &
                     .NOT. BTEST( wall_flags_0(i-1,j+1), 0 ) )                 &
                THEN
                   p_id = grid(i-1,j+1)%polygon_id
                   polygons(p_id)%nov = polygons(p_id)%nov + 1
                   nov = polygons(p_id)%nov
                   CALL set_vertex(dummy_vertex, i*dx, (j+1)*dy)
                   CALL add_vertex_to_polygon(dummy_vertex, p_id, nov)
                ENDIF
!
!--             Corner at north right edge of grid cell
                IF ( .NOT. BTEST( wall_flags_0(i+1,j), 0 ) .AND.               &
                     .NOT. BTEST( wall_flags_0(i,j+1), 0 ) .AND.               &
                     .NOT. BTEST( wall_flags_0(i+1,j+1), 0 ) )                 &
                THEN
                   p_id = grid(i+1,j+1)%polygon_id
                   polygons(p_id)%nov = polygons(p_id)%nov + 1
                   nov = polygons(p_id)%nov
                   CALL set_vertex(dummy_vertex, (i+1)*dx, (j+1)*dy)
                   CALL add_vertex_to_polygon(dummy_vertex, p_id, nov)
                ENDIF
!
!--             Corner at south right edge of grid cell
                IF ( .NOT. BTEST( wall_flags_0(i+1,j), 0 ) .AND.               &
                     .NOT. BTEST( wall_flags_0(i,j-1), 0 ) .AND.               &
                     .NOT. BTEST( wall_flags_0(i+1,j-1), 0 ) )                 &
                THEN
                   p_id = grid(i+1,j-1)%polygon_id
                   polygons(p_id)%nov = polygons(p_id)%nov + 1
                   nov = polygons(p_id)%nov
                   CALL set_vertex(dummy_vertex, (i+1)*dx, j*dy)
                   CALL add_vertex_to_polygon(dummy_vertex, p_id, nov)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

    END SUBROUTINE identify_corners

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes a vertex
!------------------------------------------------------------------------------!
    SUBROUTINE set_vertex (in_vertex, x, y)

       IMPLICIT NONE

       REAL(wp) ::  x  !< x-coordinate of vertex position
       REAL(wp) ::  y  !< y-coordinate of vertex position

       TYPE(vertex_type) ::  in_vertex  !< vertex to be set

       in_vertex%delete             = .FALSE.
       in_vertex%x                  = x
       in_vertex%y                  = y

    END SUBROUTINE set_vertex

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Adds an existing vertex to the polygon with ID p_id at position in_nov
!------------------------------------------------------------------------------!
    SUBROUTINE add_vertex_to_polygon ( in_vertex, p_id, in_nov)

       IMPLICIT NONE

       INTEGER(iwp) ::  in_nov  !< counter of vertex being added to polygon
       INTEGER(iwp) ::  p_id    !< polygon ID
       INTEGER(iwp) ::  sop     !< size of vertices array

       TYPE(vertex_type) ::  in_vertex  !< vertex to be added

       TYPE(vertex_type), DIMENSION(:), ALLOCATABLE ::  dummy_v_list  !< for reallocation

       polygon => polygons(p_id)
!
!--    Allocate and initialize the vertex array of the polygon, if necessary
       IF ( .NOT. ALLOCATED(polygon%vertices) ) THEN
          ALLOCATE(polygon%vertices(1:100))
          polygon%vertices = null_vertex
       ENDIF
!
!--    Adjust size of polygon, if necessary
       sop = SIZE(polygon%vertices)
       IF ( in_nov > sop ) THEN
          ALLOCATE(dummy_v_list(1:sop))
          dummy_v_list(1:sop) = polygon%vertices(1:sop)
          DEALLOCATE(polygon%vertices)
          ALLOCATE(polygon%vertices(1:2*sop))
          polygon%vertices = null_vertex
          polygon%vertices(1:sop) = dummy_v_list(1:sop)
          DEALLOCATE(dummy_v_list)
       ENDIF
       polygon%vertices(in_nov) = in_vertex
    END SUBROUTINE add_vertex_to_polygon

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Sorts the vertices of a polygon in a counter-clockwise fashion. During this
!> process, all vertices that are not part of the hull of the building
!> (inner courtyards) are deleted.
!------------------------------------------------------------------------------!
    SUBROUTINE sort_polygon(i_p)

       IMPLICIT NONE

       LOGICAL :: starting_vertex_found

       INTEGER(iwp) ::  counter       !< counter for potential starting vertices
       INTEGER(iwp) ::  id_neighbor   !< final ID of neighboring vertex
       INTEGER(iwp) ::  id_neighbor1  !< ID of first potential neighbor
       INTEGER(iwp) ::  id_neighbor2  !< ID of second potential neighbor
       INTEGER(iwp) ::  il            !< local counter
       INTEGER(iwp) ::  i_p           !< index of the current polygon
       INTEGER(iwp) ::  noc           !< number of candidates
       INTEGER(iwp) ::  nosv          !< number of sorted vertices
       INTEGER(iwp) ::  xe            !< x-end-index for building search
       INTEGER(iwp) ::  xs            !< x-start-index for building search
       INTEGER(iwp) ::  ye            !< y-end-index for building search
       INTEGER(iwp) ::  ys            !< y-start-index for building search

       INTEGER, DIMENSION(:), ALLOCATABLE ::  candidate_id  !< ID of the potential neighbors stored in 'candidates'
       INTEGER, DIMENSION(:), ALLOCATABLE ::  dummy_id_arr  !< used for resizing

       REAL(wp) ::  dist  !< distance of one vertex to its neighbor
       REAL(wp) ::  m_x   !< min/max x-value of polygon used for starting vertex
       REAL(wp) ::  m_y   !< min/max y-value of polygon used for starting vertex

       TYPE(vertex_type) ::  current_v     !< current vertex
       TYPE(vertex_type) ::  dummy_vertex  !< dummy vertex for reordering

       TYPE(vertex_type), DIMENSION(:), ALLOCATABLE ::  candidates        !< potential neighbors of the current vertex 
       TYPE(vertex_type), DIMENSION(:), ALLOCATABLE ::  dummy_vertex_arr  !< used for resizing
       TYPE(vertex_type), DIMENSION(:), ALLOCATABLE ::  sorted_p          !< vertices that have been sorted

       starting_vertex_found = .FALSE.
       ALLOCATE(sorted_p(1:nov))
       sorted_p(1:nov) = polygon%vertices(1:nov)
!
!--    Identify a vertex that is certainly a part of the outer hull of the 
!--    current polygon: Get rightmost border of polygon (or if that 
!--    coincides with model border, leftmost) border. Then of those points,
!--    get northmost (or if that coincides with model domain border, southmost).
!--    This identifies exactly one point that is then set to the first index.
       counter = 0
       IF ( MAXVAL(sorted_p%x) < nx*dx ) THEN
          m_x = MAXVAL(sorted_p%x)
       ELSE
          m_x = MINVAL(sorted_p%x)
       ENDIF
       DO il = 1, nov
          IF ( ABS( sorted_p(il)%x - m_x ) < .01 * dx ) THEN
             counter = counter + 1
             dummy_vertex = sorted_p(il)
             sorted_p(il) = sorted_p(counter)
             sorted_p(counter) = dummy_vertex
          ENDIF
       ENDDO
       IF ( MAXVAL(sorted_p(1:counter)%y) < ny*dy ) THEN
          m_y = MAXVAL(sorted_p(1:counter)%y)
       ELSE
          m_y = MINVAL(sorted_p(1:counter)%y)
       ENDIF
       DO il = 1, counter
          IF ( ABS(sorted_p(il)%y - m_y) < .01 * dy ) THEN
             dummy_vertex = sorted_p(il)
             sorted_p(il) = sorted_p(1)
             sorted_p(1) = dummy_vertex
             starting_vertex_found = .TRUE.
             EXIT
          ENDIF
       ENDDO
!
!--    If no starting vertex was found for the current polygon, it will be
!--    deleted and an error message thrown
       IF ( .NOT. starting_vertex_found ) THEN
          WRITE(*,'(A,/,A,1X,I6,/,A)')                                          &
                     'An error occured during polygon sorting:',               &
                     'no starting vertex could be found for polygon',          &
                     i_p, 'This polygon contains the following vertices (x/y)'
          DO il = 1, nov
             WRITE(*,'(4X,F8.1,1X,F8.1)')                                       &
                         polygon%vertices(il)%x, polygon%vertices(il)%x
          ENDDO
          WRITE(*,'(A,/,A)')                                                   &
                     'This polygon will be skipped during sorting and deleted',&
                     'For details on the procedure, see SUBROUTINE sort_polygon.'
          polygon%vertices%delete = .TRUE.
          polygons(i_p)%nov = 0
          CALL delete_empty_polygons
!
!--    Find the unique neighbor of the current vertex. For this, first
!--    determine all possible candidates. Of those, keep only the ones that
!--    are connected to the current vertex along a building edge (the polygon
!--    is sorted counter-clockwise. Therefore, the building is always on the 
!--    left-hand side of the connecting line from the current vertex to its
!--    potential neighbor). This leaves a maximum of two possible neighbors.
!--    This is only the case if the current vertex is the point that diagonally
!--    connects two parts of the same building. In that case, the vertex that
!--    lies to the right of the connecting line between the previous and
!--    current vertex is the neighbor.
       ELSE
          DO nosv = 1, nov
             current_v = sorted_p(nosv)
             noc = 0
             ALLOCATE(candidates(1:100), candidate_id(1:100))
!
!--          Find all candidates that could be neighbors of current vertex:
!--          these are those vertices that share the same x- or y-coordinate
!--          with the current vertex, as the vertices are all inner and outer
!--          corners of the gridded building data
             IF ( nosv < nov ) THEN
                DO il = nosv+1, nov
                   IF ( ABS( current_v%x - sorted_p(il)%x ) < .01 * dx .OR.    &
                        ABS( current_v%y - sorted_p(il)%y ) < .01 * dy )       &
                   THEN
!
!--                   If necessary, resize arrays for candidates
                      IF ( noc >= SIZE(candidates) ) THEN
                         ALLOCATE(dummy_vertex_arr(1:noc), dummy_id_arr(1:noc))
                         dummy_vertex_arr(1:noc) = candidates(1:noc)
                         dummy_id_arr(1:noc)   = candidate_id(1:noc)
                         DEALLOCATE(candidates, candidate_id)
                         ALLOCATE(candidates(1:2*noc), candidate_id(1:2*noc))
                         candidates(1:noc)   = dummy_vertex_arr(1:noc)
                         candidate_id(1:noc) = dummy_id_arr(1:noc)
                         DEALLOCATE(dummy_vertex_arr, dummy_id_arr)
                      ENDIF
                      noc               = noc +1
                      candidates(noc)   = sorted_p(il)
                      candidate_id(noc) = il
                   ENDIF
                ENDDO
             ENDIF
!
!--          Check which one of the candidates is the neighbor of the current
!--          vertex. This is done by several tests that would exclude the
!--          candidate from being the neighbor. Each successful test will
!--          therefore result in a cycle to the next candidate. Only if all
!--          all tests fail, is the candidate one of a maximum of two possible
!--          neighbors.
             id_neighbor  = -999
             id_neighbor1 = -999
             id_neighbor2 = -999
             DO il = 1, noc
!
!--             Exclude the possibility of a vertex with the same coordinates
!--             being chosen as the neighbor. (dist < .9*dx)
!--             NOTE: this could happen, if part of a building is only connected
!--                   to the rest of the building diagonally. In that case, the
!--                   same point is added to the polygon twice. This is necessary
!--                   and not redundant! Two such points can never be neighbors.
!--                   Example: the north right corner of grid point i,j
!--                        AND the south left corner of grid point i+1,j+1.
!--                   See SUBROUTINE identify_corners for the identification
!--                   method.
!--             Also, exclude a connection back to the coordinates of the 
!--             previous vertex.
                dist = SQRT( (candidates(il)%x - current_v%x)**2 +             &
                                     (candidates(il)%y - current_v%y)**2 )
                IF ( nosv > 1 ) THEN
                   IF ( dist < .9 * dx .OR.                                    &
                        ( ( ABS( sorted_p(nosv-1)%x                            &
                        - candidates(il)%x ) < .01 * dx ) .AND.                &
                          ( ABS( sorted_p(nosv-1)%y                            &
                        - candidates(il)%y ) < .01 * dy ) ) )                  &
                   THEN
                      CYCLE
                   ENDIF
                ENDIF
!
!--             Check if there is a building all along only the left-hand side
!--             of the connecting line from current vertex to potential neighbor
!--             (4 cases)
!--             First: for vertical connection
                IF ( ABS( candidates(il)%x - current_v%x ) < .01 * dx ) THEN
                   xs = NINT(current_v%x*ddx)-1
                   xe = xs + 1
 !
!--                Case 1: ys < ye, edge from south to north, building must be
!--                        exclusively in all grid cells left of the edge
                   IF ( current_v%y < candidates(il)%y ) THEN
                      ys = NINT(current_v%y*ddy)
                      ye = NINT(candidates(il)%y*ddy)-1
                      IF ( .NOT.( ALL( .NOT. BTEST( wall_flags_0(xs,ys:ye), 0))&
                           .AND.( ALL( BTEST( wall_flags_0(xe,ys:ye), 0 ) ) )))&
                      THEN
                         CYCLE
                      ENDIF
 !
!--                Case 2: ys > ye, edge from north to south, building must be
!--                        exclusively in all grid cells right of the edge
                   ELSEIF ( current_v%y > candidates(il)%y ) THEN
                      ys = NINT(current_v%y*ddy)-1
                      ye = NINT(candidates(il)%y*ddy)
                      IF ( .NOT.( ALL( .NOT. BTEST( wall_flags_0(xe,ye:ys), 0))&
                           .AND.( ALL( BTEST( wall_flags_0(xs,ye:ys), 0 ) ) )))&
                      THEN
                         CYCLE
                      ENDIF
                   ENDIF
!
!--             Horizontal connection
                ELSEIF ( ABS( candidates(il)%y - current_v%y ) < .01 * dy ) THEN

                   ys = NINT(current_v%y*ddy)-1
                   ye = ys + 1
!
!--                Case 3: xs > xe, edge from right to left, building must be
!--                        exclusively in all grid cells south of the edge
                   IF ( current_v%x > candidates(il)%x ) THEN
                      xs = NINT(current_v%x*ddx)-1
                      xe = NINT(candidates(il)%x*ddx)
                      IF ( .NOT.( ALL( .NOT. BTEST( wall_flags_0(xe:xs,ys), 0))&
                           .AND.( ALL( BTEST( wall_flags_0(xe:xs,ye), 0 ) ) )))&
                      THEN
                         CYCLE
                      ENDIF
!
!--                Case 4: xs < xe, edge from left to right, building must be
!--                        exclusively in all grid cells north of the edge
                   ELSEIF ( current_v%x < candidates(il)%x ) THEN
                      xs = NINT(current_v%x*ddx)
                      xe = NINT(candidates(il)%x*ddx)-1
                      IF ( .NOT.( ALL( .NOT. BTEST( wall_flags_0(xs:xe,ye), 0))&
                           .AND.( ALL( BTEST( wall_flags_0(xs:xe,ys), 0 ) ) )))&
                      THEN
                         CYCLE
                      ENDIF
                   ENDIF
                ENDIF
!
!--             After the tests, only two potential neighbors are possible. The
!--             one found first will get id_neighbor1, the possible 2nd one will
!--             get id_neighbor2
                IF ( id_neighbor1 ==  -999 ) THEN
                   id_neighbor1 = candidate_id(il)
                ELSEIF ( id_neighbor1 /=  -999 .AND.                           &
                       ( ( ABS( sorted_p(id_neighbor1)%x - candidates(il)%x )  &
                                > .01 * dx ) .OR.                              &
                         ( ABS( sorted_p(id_neighbor1)%y - candidates(il)%y )  &
                                > .01 * dy ) ) )                               &
                THEN
                   id_neighbor2 = candidate_id(il)
                ENDIF
             ENDDO
!
!--          If two potential neighbors were found, determine the one that is on 
!--          the right hand side of the line connecting the current and previous 
!--          vertex. It is the real neighbor.
             IF ( id_neighbor2 /= -999 .AND. nosv > 1 ) THEN
                IF ( is_right(sorted_p(nosv-1)%x,sorted_p(nosv-1)%y,           &
                           current_v%x,current_v%y,                            &
                           sorted_p(id_neighbor1)%x,sorted_p(id_neighbor1)%y) )&
                THEN
                   id_neighbor = id_neighbor1
                ELSEIF ( is_right(sorted_p(nosv-1)%x,sorted_p(nosv-1)%y,       &
                           current_v%x,current_v%y,                            &
                           sorted_p(id_neighbor2)%x,sorted_p(id_neighbor2)%y) )&
                THEN
                   id_neighbor = id_neighbor2
                ENDIF
             ELSE
                id_neighbor = id_neighbor1
             ENDIF
!
!--          Put the found neighbor at next index in sorted array and move the
!--          unsorted vertices back one index. This way, only yet unsorted
!--          vertices are eligible to be candidates during the next iteration.
             IF (id_neighbor /= nosv + 1 .AND. id_neighbor /= -999) THEN
                dummy_vertex = sorted_p(id_neighbor)
                sorted_p(nosv+2:id_neighbor) = sorted_p(nosv+1:id_neighbor-1)
                sorted_p(nosv+1) = dummy_vertex
!
!--          If no neighbor was found, sorting is done for this polygon
             ELSEIF ( id_neighbor == -999 ) THEN
                DEALLOCATE(candidates,candidate_id)
                EXIT
             ENDIF
             DEALLOCATE(candidates,candidate_id)
          ENDDO
!
!--       Sorting is done. Reduce size (which means get rid of vertices
!--       that are not part of the outer hull of the building: holes)
!--       of sorted polygon and put it back in polygon%vertices.
!--       Also add first vertex to the end of polygon and last vertex 
!--       before the beginning of polygon.
          DEALLOCATE(polygon%vertices)
          ALLOCATE(polygon%vertices(0:nosv+1))
          polygon%vertices(1:nosv) = sorted_p(1:nosv)
          polygon%vertices(0) = sorted_p(nosv)
          polygon%vertices(nosv+1) = sorted_p(1)
          polygons(i_p)%nov = nosv
          nov = polygons(i_p)%nov
          DEALLOCATE(sorted_p)
       ENDIF

    END SUBROUTINE sort_polygon

!------------------------------------------------------------------------------! 
! Description:
! ------------
!> Reduces the number of vertices in a polygon using the 
!> Douglas-Poiker-Algorithm (1973)
!------------------------------------------------------------------------------! 
    RECURSIVE SUBROUTINE simplify_polygon( id_s, id_e, tol )

       IMPLICIT NONE

       INTEGER(iwp) ::  max_dist_ind  !< Index of vertex with maximum distance
       INTEGER(iwp) ::  il            !< counter
       INTEGER(iwp) ::  id_s          !< End index in polygon
       INTEGER(iwp) ::  id_e          !< End index in polygon

       REAL(wp) ::  max_dist  !< Maximum distance from line
       REAL(wp) ::  dum_dist  !< Distance from line: dummy
       REAL(wp) ::  tol       !< factor that determines how far a vertex can be from the polygon approximation so that the approximation is still accepted

       max_dist = 0.
       max_dist_ind = -999999
!
!--    Find vertex with max distance to id_s and id_e
       DO il = id_s + 1, id_e -1
          dum_dist = dist_point_to_edge(polygon%vertices(id_s)%x,         &
                                        polygon%vertices(id_s)%y,         &
                                        polygon%vertices(id_e)%x,         &
                                        polygon%vertices(id_e)%y,         &
                                        polygon%vertices(il)%x,           &
                                        polygon%vertices(il)%y)
          IF ( dum_dist > max_dist ) THEN
             max_dist     = dum_dist
             max_dist_ind = il
          ENDIF
       ENDDO

       IF ( max_dist > tol ) THEN
          CALL simplify_polygon( id_s, max_dist_ind, tol )
          CALL simplify_polygon( max_dist_ind, id_e, tol )
       ELSE
          polygon%vertices(id_s+1:id_e-1)%delete = .TRUE.
       ENDIF

    END SUBROUTINE simplify_polygon

!------------------------------------------------------------------------------! 
! Description:
! ------------
!> Checks if a vertex of a polygon is inside another polygon and if so, deletes
!> it. The check is done using the crossing number algorithm. If a straight
!> ray starting at a point crosses the borders of one polygon an odd
!> number of times, the point is inside that polygon.
!> This algorithm detects buildings that are completely surrounded by
!> another building. They can be deleted since they can never be navigated.
!> TODO: Maybe add a flag to turn this off and on as it might not be needed.
!>       also, if the domain has buildings at all boundary points, there would
!>       only be one giant building and everything in it deleted. So nothing
!>       could be navigated. relevant?!
!------------------------------------------------------------------------------! 
    SUBROUTINE inside_other_polygon( i_p )

       IMPLICIT NONE

       LOGICAL ::  exit_flag  !< flag to exit loops if an odd crossing number was found for any of a polygons vertices

       INTEGER(iwp) ::  cn         !< number of crossings
       INTEGER(iwp) ::  i_p        !< index of current polygon
       INTEGER(iwp) ::  il         !< index of tested polygon
       INTEGER(iwp) ::  nov_test   !< no. of vertices of test-polygon
       INTEGER(iwp) ::  ref_vert   !< vertex currently being tested if it is inside another polygon
       INTEGER(iwp) ::  test_edge  !< index of edge being tested

       REAL(wp) ::  px  !< x-coord of the point at the crossing of the ray and the vertex
       REAL(wp) ::  xe  !< x-coordinate of end point of edge
       REAL(wp) ::  xr  !< x-coordinate of reference point
       REAL(wp) ::  xs  !< x-coordinate of start point of edge
       REAL(wp) ::  ye  !< y-coordinate of end point of edge
       REAL(wp) ::  yr  !< y-coordinate of reference point
       REAL(wp) ::  ys  !< y-coordinate of start point of edge

       TYPE(polygon_type), POINTER ::  test_pol  !< Polygon to be tested

       exit_flag = .FALSE.
!
!--    Loop over all polygons other than the one being tested
       DO il = 1, polygon_counter
          IF ( il == i_p ) CYCLE
          test_pol => polygons(il)
          nov_test = polygons(il)%nov
!
!--       Inclusion test is done for every vertex of the polygon
          DO ref_vert = 1, nov
             cn = 0
             xr = polygon%vertices(ref_vert)%x
             yr = polygon%vertices(ref_vert)%y
!
!--          All edges of the every polygon il is tested for ray crossing
             DO test_edge = 1, nov_test

!--             It is tested wether the current edge crosses a ray that extends
!--             from the current point to the right indefinitely.
!--             Check if start point of edge is lower than end point. If they
!--             are the same, ignore, since horizontal edges are excluded
                IF ( test_pol%vertices(test_edge)%y <                          &
                     test_pol%vertices(test_edge+1)%y )                        &
                THEN
                   xs = test_pol%vertices(test_edge)%x
                   xe = test_pol%vertices(test_edge+1)%x
                   ys = test_pol%vertices(test_edge)%y
                   ye = test_pol%vertices(test_edge+1)%y
                ELSEIF ( test_pol%vertices(test_edge)%y >                      &
                         test_pol%vertices(test_edge+1)%y )                    &
                THEN
                   xs = test_pol%vertices(test_edge+1)%x
                   xe = test_pol%vertices(test_edge)%x
                   ys = test_pol%vertices(test_edge+1)%y
                   ye = test_pol%vertices(test_edge)%y
                ELSE
                   CYCLE
                ENDIF
!
!--             Only such edges where the starting point of the edge is south of
!--             (or equal to) the reference point and the end point is north of
!--             it are relevant. Note: an edge includes its southern endpoint
!--             and excludes its northern endpoint.
                IF ( .NOT. (ys <= yr .AND. ye > yr )) CYCLE
!
!--             Only edges that are crossed on the right side of the reference 
!--             point are relevant, those on the left are ignored
                IF ( xs <= xr .AND. xe <= xr ) CYCLE
                IF ( ( xs <= xr .AND. xe >= xr ) .OR.                          &
                     ( xs >= xr .AND. xe <= xr ) )                             &
                THEN
                   px = xe - (xe-xs)*(ye-yr)/(ye-ys)
                   IF ( px <= xr ) CYCLE
                ENDIF
!
!--             If none of the previous if clauses were true, a crossing with 
!--             an eligible edge was found and the count increases
                cn = cn + 1
             ENDDO
!
!--          If the number of crossings is odd, the point is inside another 
!--          polyon. The polygon associated with the point will be deleted
             IF (  MOD(cn, 2) /= 0 ) THEN
                exit_flag = .TRUE.
                EXIT
             ENDIF
          ENDDO
          IF ( exit_flag ) EXIT
       ENDDO
       IF ( exit_flag ) polygon%vertices%delete = .TRUE.

    END SUBROUTINE inside_other_polygon

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deletes thoses vertices that are marked for deletion (%delete flag) and
!> resizes the polygon
!------------------------------------------------------------------------------!
    SUBROUTINE delete_extra_vertices (i_p)

       IMPLICIT NONE

       INTEGER(iwp) ::  il        !< local counter
       INTEGER(iwp) ::  vcounter  !< vertex counter
       INTEGER(iwp) ::  i_p       !< polygon ID

       TYPE(vertex_type), DIMENSION(:), ALLOCATABLE ::  dummy_pol !< Temporarily stores non-deleted vertices

       ALLOCATE(dummy_pol(1:nov))
       vcounter = 0
!
!--    Check all vertices and only keep those not marked for deletion
       DO il = 1, nov
          IF ( .NOT. polygon%vertices(il)%delete ) THEN
             vcounter = vcounter + 1
             dummy_pol(vcounter) = polygon%vertices(il)
          ENDIF
       ENDDO
!
!--    Set new number of vertices in the polygon
       nov = vcounter
       polygons(i_p)%nov = nov
!
!--    Resize
       DEALLOCATE(polygon%vertices)
       ALLOCATE(polygon%vertices(0:nov+1))
       polygon%vertices(1:nov) = dummy_pol(1:nov)
       polygon%vertices(0) = polygon%vertices(nov)
       polygon%vertices(nov+1) = polygon%vertices(1)
       DEALLOCATE(dummy_pol)

    END SUBROUTINE delete_extra_vertices

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Deletes polygons that contain no vertices (happens for those polygons that
!> were entirely encompassed by another polygon)
!------------------------------------------------------------------------------!
    SUBROUTINE delete_empty_polygons

       IMPLICIT NONE

       INTEGER(iwp) ::  il  !< local counter
       INTEGER(iwp) ::  pc  !< number of nonempty polygons
       INTEGER(iwp) ::  sv  !< size of vertex array

       TYPE(polygon_type), DIMENSION(:), ALLOCATABLE ::  dummy_polygons  !< temporarily stores non-deletd polygons

       pc = 0
       sv = 0
       ALLOCATE( dummy_polygons(1:polygon_counter) )
!
!--    Keep only those polygons that contain any vertices, skip the rest
       DO il = 1, polygon_counter
          IF ( polygons(il)%nov > 0 ) THEN
             pc = pc + 1
             sv = SIZE(polygons(il)%vertices)
             ALLOCATE(dummy_polygons(pc)%vertices(0:sv-1))
             dummy_polygons(pc) = polygons(il)
          ENDIF
       ENDDO
       polygon_counter = pc
!
!--    Resize polygon array
       DEALLOCATE(polygons)
       ALLOCATE(polygons(1:polygon_counter))
       DO il = 1, polygon_counter
!
!--       give each %vertices array the correct size and information
          sv = SIZE(dummy_polygons(il)%vertices)
          polygons(il)%nov = sv - 2
          ALLOCATE(polygons(il)%vertices(0:sv-1))
          polygons(il) = dummy_polygons(il)
       ENDDO
       DEALLOCATE(dummy_polygons)

    END SUBROUTINE delete_empty_polygons

 END MODULE polygon_creation 

 MODULE mesh_creation

    USE kinds

    USE mod_functions

    USE variables

    CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Creates the navigation mesh:
!>      1) Finds eligible vertices (those that are locally convex)
!>      2) Adds them to the mesh
!>      3) Adds connections between mesh points if they are in line of sight
!>         of each other and the connecting line does not point into either of
!>         the originating polygons (this is known as a visibility graph)
!------------------------------------------------------------------------------!
    SUBROUTINE create_nav_mesh

       IMPLICIT NONE

       LOGICAL ::  add                 !< flag for second cycle of add loop
       LOGICAL ::  intersection_found  !< flag to indicate a found intersection

       INTEGER(iwp) ::  cmp    !< counter: current mesh point
       INTEGER(iwp) ::  il     !< local counter
       INTEGER(iwp) ::  jl     !< local counter
       INTEGER(iwp) ::  pid    !< polygon id of current mesh point
       INTEGER(iwp) ::  pid_t  !< polygon id of tested mesh point
       INTEGER(iwp) ::  pl     !< polygon counter
       INTEGER(iwp) ::  vid    !< vertex id of current mesh point
       INTEGER(iwp) ::  vid_t  !< vertex id of tested mesh point
       INTEGER(iwp) ::  vl     !< vertex counter

       REAL(wp) ::  v1x           !< x-coordinate of test vertex 1 for intersection test
       REAL(wp) ::  v1y           !< y-coordinate of test vertex 1 for intersection test
       REAL(wp) ::  v2x           !< x-coordinate of test vertex 2 for intersection test
       REAL(wp) ::  v2y           !< y-coordinate of test vertex 2 for intersection test
       REAL(wp) ::  x             !< x-coordinate of current mesh point
       REAL(wp) ::  x_t           !< x-coordinate of tested mesh point
       REAL(wp) ::  y             !< y-coordinate of current mesh point
       REAL(wp) ::  y_t           !< y-coordinate of tested mesh point
       REAL(wp) ::  corner_x      !< x-coordinate of shifted corner
       REAL(wp) ::  corner_x_e    !< x-coordinate of end of corner gate
       REAL(wp) ::  corner_y      !< y-coordinate of shifted corner
       REAL(wp) ::  corner_y_e    !< y-coordinate of end of corner gate
       REAL(wp) ::  t_start       !< CPU measure: start
       REAL(wp) ::  t_inter       !< CPU measure: output test time
       REAL(wp) ::  t_inter1      !< CPU measure: output test time
       REAL(wp) ::  t_end         !< CPU measure: end
       REAL(wp) ::  t_left        !< CPU measure: estimate for time left
       REAL(wp) ::  t_done        !< CPU measure: elapsed time
       REAL(wp) ::  percent_done  !< CPU measure: proportion of mesh points checked

!
!--    Add all convex vertices to the mesh.
!--    DO loop will be executed twice. Once to count the mesh points to be
!--    added and allocate the mesh point array, the second time (add == .TRUE.)
!--    to fill the mesh point array.
       WRITE(*,'(1X,A)') 'Adding polygon vertices to mesh ...'
       add = .FALSE.
       DO
          cmp = 0
          DO il = 1, polygon_counter
             polygon => polygons(il)
             nov = polygons(il)%nov
             DO jl = 1, nov
!
!--             In a polygon that is sorted counter-clockwise, if the next vertex
!--             is left of the line connecting the previous and the current vertex,
!--             the current vertex is locally convex.
                IF ( is_left(polygon%vertices(jl-1)%x,polygon%vertices(jl-1)%y,   &
                             polygon%vertices(jl)%x,polygon%vertices(jl)%y,       &
                             polygon%vertices(jl+1)%x,polygon%vertices(jl+1)%y) ) &
                THEN

                   corner_x = polygon%vertices(jl)%x
                   corner_y = polygon%vertices(jl)%y
!
!--                Create end point for corner navigation
                   IF ( add ) THEN
                      CALL shift_corner_outward(                               &
                            polygon%vertices(jl-1)%x, polygon%vertices(jl-1)%y,&
                            polygon%vertices(jl+1)%x, polygon%vertices(jl+1)%y,&
                            polygon%vertices(jl)%x,   polygon%vertices(jl)%y,  &
                            corner_x_e, corner_y_e, 1._wp )
                   ENDIF
!
!--                Disregard corners outside of the domain
                   IF ( corner_x<=(nx+1)*dx .AND. corner_x>=0 .AND.            &
                        corner_y<=(ny+1)*dy .AND. corner_y>=0)                 &
                   THEN
                      cmp = cmp + 1
                      IF ( add ) THEN
                         CALL set_mesh_point( mesh(cmp), il, jl,               &
                                                      corner_x, corner_y,      &
                                                      corner_x_e, corner_y_e )
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          IF ( add ) EXIT
          add = .TRUE.
          ALLOCATE( mesh(1:cmp) )
       ENDDO
       WRITE(*,'(6X,A,1X,I10,1X,A,/)')  'Done. Added',cmp,'vertices to mesh.'
       WRITE(*,'(1X,A)') 'Establishing connections in mesh ...'
!
!--    CPU measurement
       CALL CPU_TIME(t_start)
       CALL CPU_TIME(t_inter)
       DO il = 1, cmp
!--       Output status of processing
          CALL CPU_TIME(t_inter1)
          IF ( t_inter1 - t_inter > 4. ) THEN
             t_done       = (t_inter1-t_start)/60.
             percent_done = REAL(il)/cmp*100.
             t_left       = t_done/percent_done*(100-percent_done)
             WRITE(*,'(3X,2(A,I8),A,F6.2,2(A,F7.1),A,I10)')                    &
                   'Mesh point ',il,' of '                ,cmp,                &
                   ': '                                   ,percent_done,       &
                   ' % || elapsed time : '                ,t_done,             &
                   ' min || ETA: '                        ,t_left,             &
                   ' min || number of connections found: ',number_of_connections
             CALL CPU_TIME(t_inter)
          ENDIF
          x = mesh(il)%x
          y = mesh(il)%y
          pid = mesh(il)%polygon_id
          vid = mesh(il)%vertex_id
          DO jl = 1, cmp
!
!--          No mesh point can be connected to itself
             IF ( il == jl ) CYCLE
             x_t = mesh(jl)%x
             y_t = mesh(jl)%y
             pid_t = mesh(jl)%polygon_id
             vid_t = mesh(jl)%vertex_id
!
!--          Cycle, if a connection had already been established
             IF ( ANY(mesh(jl)%connected_vertices == il) ) CYCLE
!
!--          If the distance between two nodes is larger than 600 m,
!--          no connection will be made since there will typically no be such
!--          long, straight ways in a city that a pedestrian will walk
             IF ( SQRT((x_t-x)**2 +(y_t-y)**2) > 400. ) CYCLE
!
!--          If the connecting line between two mesh points points into either
!--          or both of the corresponding polygons, no connection will be 
!--          established between the two points. This is the case if the 
!--          previous (next) vertex of the polygon is right of the connecting
!--          line and the next (previous) vertex of the polygon is left of the
!--          connecting line. This is checked for both polygons.
             IF ( ((is_left(x_t,y_t,x,y,polygons(pid)%vertices(vid-1)%x,       &
                                       polygons(pid)%vertices(vid-1)%y)        &
                  .AND. is_right(x_t,y_t,x,y,polygons(pid)%vertices(vid+1)%x,  &
                                       polygons(pid)%vertices(vid+1)%y) )      &
                  .OR. (is_right(x_t,y_t,x,y,polygons(pid)%vertices(vid-1)%x,  &
                                       polygons(pid)%vertices(vid-1)%y)        &
                  .AND. is_left(x_t,y_t,x,y,polygons(pid)%vertices(vid+1)%x,   &
                                       polygons(pid)%vertices(vid+1)%y)) )     &
                  .OR. ((is_left(x,y,x_t,y_t,polygons(pid_t)%vertices(vid_t-1)%x, &
                                       polygons(pid_t)%vertices(vid_t-1)%y)        &
                  .AND. is_right(x,y,x_t,y_t,polygons(pid_t)%vertices(vid_t+1)%x,  &
                                       polygons(pid_t)%vertices(vid_t+1)%y) )      &
                  .OR. (is_right(x,y,x_t,y_t,polygons(pid_t)%vertices(vid_t-1)%x,  &
                                       polygons(pid_t)%vertices(vid_t-1)%y)        &
                  .AND. is_left(x,y,x_t,y_t,polygons(pid_t)%vertices(vid_t+1)%x,   &
                                       polygons(pid_t)%vertices(vid_t+1)%y)) ) )   &
             THEN
                CYCLE
             ENDIF
!
!--          For each edge of each polygon, check if it intersects with the 
!--          potential connection. If so, no connection can be made
!--          THIS IS THE BOTTLENECK OF THE PROGRAM
             intersection_found = .FALSE.
             DO pl = pid, polygon_counter
                DO vl = 1, polygons(pl)%nov
                   v1x = polygons(pl)%vertices(vl)%x
                   v1y = polygons(pl)%vertices(vl)%y
                   v2x = polygons(pl)%vertices(vl+1)%x
                   v2y = polygons(pl)%vertices(vl+1)%y
                   intersection_found = intersect(x,y,x_t,y_t,v1x,v1y,v2x,v2y)
                   IF ( intersection_found ) EXIT
                ENDDO
                IF ( intersection_found ) EXIT
             ENDDO
             IF ( intersection_found ) CYCLE
             DO pl = pid, 1, -1
                IF ( pl == pid ) CYCLE
                DO vl = 1, polygons(pl)%nov
                   v1x = polygons(pl)%vertices(vl)%x
                   v1y = polygons(pl)%vertices(vl)%y
                   v2x = polygons(pl)%vertices(vl+1)%x
                   v2y = polygons(pl)%vertices(vl+1)%y
                   intersection_found = intersect(x,y,x_t,y_t,v1x,v1y,v2x,v2y)
                   IF ( intersection_found ) EXIT
                ENDDO
                IF ( intersection_found ) EXIT
             ENDDO
             IF ( intersection_found ) CYCLE
!
!--       If neither of the above two test was true, a connection will be
!--       established between the two mesh points.
          number_of_connections = number_of_connections + 1
          CALL add_connection(mesh(il),jl, mesh(jl))
          CALL add_connection(mesh(jl),il, mesh(il))
          ENDDO
       ENDDO
!
!--    Adapt connected_vertices arrays
       DO il = 1, cmp
          CALL reduce_connections(mesh(il))
       ENDDO
       CALL CPU_TIME(t_end)
!
!--    Output to terminal
       WRITE(*,'(6X,A,I10,A)') 'Done. Established ',number_of_connections,     &
                               ' connections in mesh'
       WRITE(*,'(6X,A,F10.1,A)') 'Time needed for calculation: ',              &
                                 t_end-t_start,' seconds'

    END SUBROUTINE create_nav_mesh

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializes a point of the navigation mesh
!------------------------------------------------------------------------------!
    SUBROUTINE set_mesh_point (in_mp,pid,vid,x,y,x_s,y_s)

       IMPLICIT NONE

       INTEGER(iwp) ::  pid  !< polygon ID
       INTEGER(iwp) ::  vid  !< vertex ID

       REAL(wp) ::  x    !< x-value of mesh point for path calculation
       REAL(wp) ::  x_s  !< x-value shifted outward from corner
       REAL(wp) ::  y    !< y-value of mesh point for path calculation
       REAL(wp) ::  y_s  !< y-value shifted outward from corner

       TYPE(mesh_point) ::  in_mp  !< mesh point to be created

       in_mp%origin_id          = -1
       in_mp%polygon_id         = pid
       in_mp%vertex_id          = vid
       in_mp%cost_so_far        = 1.d12
       in_mp%x                  = x
       in_mp%y                  = y
       in_mp%x_s                = x_s
       in_mp%y_s                = y_s
       in_mp%noc                = 0

       ALLOCATE(in_mp%connected_vertices(1:100),                               &
                in_mp%distance_to_vertex(1:100))

       in_mp%connected_vertices = -999
       in_mp%distance_to_vertex = -999.

    END SUBROUTINE set_mesh_point

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Shifts a corner (middle one of three consecutive points a, b and p) outward
!> by a given length along the angle bisector. Stores the result to res_x/res_y
!------------------------------------------------------------------------------!
    SUBROUTINE shift_corner_outward ( a_x, a_y, b_x, b_y, p_x, p_y, res_x,     &
                                      res_y, shift )

       IMPLICIT NONE

       REAL(wp) ::  a_x     !< x-value of point A
       REAL(wp) ::  a_y     !< y-value of point A
       REAL(wp) ::  abs_ap  !< distance from A to P
       REAL(wp) ::  abs_bp  !< distance from B to P
       REAL(wp) ::  abs_co  !< length of angle bisector
       REAL(wp) ::  b_x     !< x-value of point B
       REAL(wp) ::  b_y     !< y-value of point B
       REAL(wp) ::  eap_x   !< x-value of unit vector from A to P
       REAL(wp) ::  eap_y   !< y-value of unit vector from A to P
       REAL(wp) ::  ebp_x   !< x-value of unit vector from B to P
       REAL(wp) ::  ebp_y   !< y-value of unit vector from B to P
       REAL(wp) ::  p_x     !< x-value of point P
       REAL(wp) ::  p_y     !< y-value of point P
       REAL(wp) ::  res_x   !< x-value of result
       REAL(wp) ::  res_y   !< y-value of result
       REAL(wp) ::  shift   !< distance of shift in meters

!
!--    Get unit vector from previous to current vertex
       eap_x  = p_x - a_x
       eap_y  = p_y - a_y
       abs_ap = SQRT(eap_x**2+eap_y**2)
       eap_x  = eap_x/abs_ap
       eap_y  = eap_y/abs_ap
!
!--    Get unit vector from next to current vertex
       ebp_x  = p_x - b_x
       ebp_y  = p_y - b_y
       abs_bp = SQRT(ebp_x**2+ebp_y**2)
       ebp_x  = ebp_x/abs_bp
       ebp_y  = ebp_y/abs_bp
!
!--    Add previous two vectors to get angle bisector of corner.
!--    Then, set its length to shift and add to original vertex
!--    vector to shift it outward
       res_x   = eap_x + ebp_x
       res_y   = eap_y + ebp_y
       abs_co  = SQRT(res_x**2+res_y**2)
       res_x   = shift*res_x/abs_co + p_x
       res_y   = shift*res_y/abs_co + p_y

    END SUBROUTINE shift_corner_outward

!------------------------------------------------------------------------------! 
! Description:
! ------------
!> Adds a connection between two points of the navigation mesh 
!> (one-way: in_mp1 to in_mp2)
!------------------------------------------------------------------------------! 
    SUBROUTINE add_connection (in_mp1,id2,in_mp2)

       IMPLICIT NONE

       LOGICAL ::  connection_established  !< Flag to indicate if connection has already been established

       INTEGER(iwp) ::  id2  !< ID of in_mp2
       INTEGER(iwp) ::  il   !< local counter
       INTEGER(iwp) ::  noc1 !< number of connections in in_mp1

       INTEGER, DIMENSION(:), ALLOCATABLE ::  dum_cv !< dummy array for connected_vertices

       REAL(wp) ::  dist  !< Distance between the two points

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dum_dtv

       TYPE(mesh_point) ::  in_mp1  !< mesh point that gets a new connection
       TYPE(mesh_point) ::  in_mp2  !< mesh point in_mp1 will be connected to

       connection_established = .FALSE.
!
!--    Check if connection has already been established
       noc1 = SIZE(in_mp1%connected_vertices)
       DO il = 1, in_mp1%noc
          IF ( in_mp1%connected_vertices(il) == id2 ) THEN
             connection_established = .TRUE.
             EXIT
          ENDIF
       ENDDO

       IF ( .NOT. connection_established ) THEN
!
!--       Resize arrays, if necessary
          IF ( in_mp1%noc >= noc1 ) THEN
             ALLOCATE( dum_cv(1:noc1),dum_dtv(1:noc1) )
             dum_cv  = in_mp1%connected_vertices
             dum_dtv = in_mp1%distance_to_vertex
             DEALLOCATE( in_mp1%connected_vertices, in_mp1%distance_to_vertex )
             ALLOCATE( in_mp1%connected_vertices(1:2*noc1),                    &
                       in_mp1%distance_to_vertex(1:2*noc1) )
             in_mp1%connected_vertices         = -999
             in_mp1%distance_to_vertex         = -999.
             in_mp1%connected_vertices(1:noc1) = dum_cv
             in_mp1%distance_to_vertex(1:noc1) = dum_dtv
          ENDIF

!
!--    Add connection
          in_mp1%noc = in_mp1%noc+1
          dist = SQRT( (in_mp1%x - in_mp2%x)**2 + (in_mp1%y - in_mp2%y)**2 )
          in_mp1%connected_vertices(in_mp1%noc) = id2
          in_mp1%distance_to_vertex(in_mp1%noc) = dist
       ENDIF

    END SUBROUTINE add_connection

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reduces the size of connection array to the amount of actual connections
!> after all connetions were added
!------------------------------------------------------------------------------!
    SUBROUTINE reduce_connections (in_mp)

       IMPLICIT NONE

       INTEGER(iwp) ::  noc !< Number of connections

       INTEGER, DIMENSION(:), ALLOCATABLE ::  dum_cv !< dummy: connected_vertices

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dum_dtv !< dummy: distance_to_vertex

       TYPE(mesh_point) ::  in_mp !< Input mesh point

       noc = in_mp%noc
       ALLOCATE( dum_cv(1:noc),dum_dtv(1:noc) )
       dum_cv  = in_mp%connected_vertices(1:noc)
       dum_dtv = in_mp%distance_to_vertex(1:noc)
       DEALLOCATE( in_mp%connected_vertices, in_mp%distance_to_vertex )
       ALLOCATE( in_mp%connected_vertices(1:noc),                              &
                 in_mp%distance_to_vertex(1:noc) )
       in_mp%connected_vertices(1:noc) = dum_cv(1:noc)
       in_mp%distance_to_vertex(1:noc) = dum_dtv(1:noc)

    END SUBROUTINE reduce_connections

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes all NavMesh information into binary file and building data to ASCII
!------------------------------------------------------------------------------!
    SUBROUTINE bin_out_mesh

       IMPLICIT NONE

       INTEGER(iwp) ::  il           !< local counter
       INTEGER(iwp) ::  jl           !< local counter
       INTEGER(iwp) ::  size_of_mesh !< size of mesh
       INTEGER(iwp) ::  size_of_pols !< size of polygon

       WRITE(*,'(1X,A)') 'Writing binary output data ...'

       OPEN ( 14, FILE= TRIM(runname)//'_nav', FORM='UNFORMATTED', STATUS='replace' )
!
!--    Output of mesh data
       size_of_mesh = SIZE(mesh)
       WRITE(14) size_of_mesh
       DO il = 1, size_of_mesh
          WRITE(14) mesh(il)%polygon_id, mesh(il)%vertex_id, mesh(il)%noc,     &
                    mesh(il)%origin_id, mesh(il)%cost_so_far, mesh(il)%x,      &
                    mesh(il)%y, mesh(il)%x_s, mesh(il)%y_s
          DO jl = 1, mesh(il)%noc
             WRITE(14) mesh(il)%connected_vertices(jl),                        &
                       mesh(il)%distance_to_vertex(jl)
          ENDDO
       ENDDO
!
!--    Output of building polygon data
       size_of_pols = SIZE(polygons)
       WRITE(14) size_of_pols
       DO il = 1, size_of_pols
          WRITE(14) polygons(il)%nov
          DO jl = 0, polygons(il)%nov+1
             WRITE(14) polygons(il)%vertices(jl)%delete,                       &
                       polygons(il)%vertices(jl)%x, polygons(il)%vertices(jl)%y
          ENDDO
       ENDDO
       CLOSE(14)
!
!--    Output building data to ASCII file
       OPEN(UNIT=7,FILE='topo.txt',STATUS='replace',ACTION='write')
       DO i = 1, polygon_counter
       IF (polygons(i)%nov == 0) CYCLE
          DO j = 1, polygons(i)%nov
             WRITE(7,150) i,j,polygons(i)%vertices(j)%x,                       &
                              polygons(i)%vertices(j)%y
          ENDDO
       ENDDO
       CLOSE(7)

       WRITE(*,'(6X,A)')  'Done, tool terminating.', '  ',                     &
                          'Before starting your PALM run, please check the',   &
                          'ASCII file topo.txt to see if you are satisfied',   &
                          'with the polygon representation of the building',   &
                          'data. If not, consider adjusting the parameter',    &
                          'tolerance_dp accordingly.', '  ', 'Bye, Bye!', ' '
       CALL CPU_TIME(finish)
       WRITE(*,'(1X,A,F10.4,A)') 'Total runtime: ', finish-start, ' seconds'

       150 FORMAT (2(I7,1X),2(F9.2,1X) )

    END SUBROUTINE bin_out_mesh

 END MODULE mesh_creation

 PROGRAM nav_mesh

    USE mesh_creation
    USE polygon_creation
    USE variables
    IMPLICIT NONE


!
!-- Start CPU mesurement
    CALL CPU_TIME(start)
!
!-- Initialization
    CALL init

    WRITE(*,*) "Converting building data to polygons ..."
!
!-- Convert gridded building data to polygons
    CALL identify_polygons
!
!-- Find corners in topography and add them to polygons
    CALL identify_corners
!
!-- Sort polygons counter-clockwise, then simplify them
    DO i = 1, polygon_counter
       polygon => polygons(i)
       nov = polygons(i)%nov
       CALL sort_polygon(i)
!
!--    Simplify each polygon using douglas-peucker algorithm. If the number
!--    of vertices would fall below 4 due to this procedure, the tolerance
!--    for the algorithm is reduced and it is run again. 
       DO i_sc = 0, 2
          CALL simplify_polygon(1,nov+1,tolerance_dp(i_sc))
          i_cn = 0
          DO j = 1, nov
             IF ( .NOT. polygon%vertices(j)%delete ) i_cn = i_cn + 1
          ENDDO
          IF ( i_cn > 3 ) THEN
             EXIT
          ELSE
             polygon%vertices(:)%delete = .FALSE.
          ENDIF
       ENDDO
       CALL delete_extra_vertices(i)
    ENDDO
!
!-- Remove buildings that are surrounded by another building
    IF ( .NOT. internal_buildings ) THEN
       DO i = 1, polygon_counter
          polygon => polygons(i)
          nov = polygons(i)%nov
          CALL inside_other_polygon(i)
       ENDDO
    ENDIF
!
!-- Delete vertices that are marked for deletion
    DO i = 1, polygon_counter
       polygon => polygons(i)
       nov = polygons(i)%nov
       CALL delete_extra_vertices(i)
    ENDDO
!
!-- Count number of vertices
    vertex_counter = 0
    DO i = 1, polygon_counter
       polygon => polygons(i)
       nov = polygons(i)%nov
       vertex_counter = vertex_counter + nov
    ENDDO
!
!-- Delete polygons with no vertices
    CALL delete_empty_polygons
    WRITE(*,'(2(6X,A,I10,1X,A,/))')                                             &
                  'Done. Created a total of', polygon_counter, 'polygon(s)',   &
                  '         with a total of', vertex_counter, 'vertices'
!
!-- Crate Navigation mesh from polygon data
    CALL create_nav_mesh
!
!-- Binary mesh output
    CALL bin_out_mesh

 END PROGRAM nav_mesh
