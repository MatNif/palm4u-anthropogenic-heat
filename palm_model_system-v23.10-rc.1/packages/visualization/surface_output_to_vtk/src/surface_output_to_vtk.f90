!> @file surface_output_to_vtk.f90
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
! Copyright 1997-2021  Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Matthias Suehring
! @author Klaus Ketelsen
!
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> The following tool combines surface output from PALM-subdomains into a single vtk file. Binary
!> output from each subdomain element is opened, read and written into the vtk file. Output for
!> multiple timesteps is written to multiple vtk files, one vtk file for each output timestep.
!> Output is distinguished between instantaneous and time-averaged data.
!------------------------------------------------------------------------------!
 PROGRAM surface_output_to_vtk

    USE, INTRINSIC ::  ISO_C_BINDING

    USE posix_interface,                                                                           &
        ONLY:  posix_ftell,                                                                        &
               posix_lseek

    IMPLICIT NONE

    CHARACTER(LEN=4)   ::  char_time              !< string indicating simulated time

    CHARACTER(LEN=10)  ::  char_dum               !< dummy string

    CHARACTER(LEN=30)  ::  myid_char              !< combined string indicating binary file


    CHARACTER(LEN=100) ::  path_in                !< path to the binary data
    CHARACTER(LEN=100) ::  path_out = ""          !< path for the output data
    CHARACTER(LEN=100) ::  variable_name          !< name of the processed output variable

    INTEGER(4)   ::  ftell                        !< intrinsic function, get current position in file
    INTEGER(4)   ::  ndum                         !< return parameter of intrinsic function fseek

    INTEGER, PARAMETER ::  iwp = 4                !< integer precision
    INTEGER, PARAMETER ::  wp  = 8                !< float precision
    !INTEGER, PARAMETER ::  OFFSET_KIND = C_SIZE_T !< unsigned integer for the C-interface

    INTEGER(iwp) ::  f                            !< running index over all binary files
    INTEGER(iwp) ::  file_id_in = 18              !< file unit for input binaray file
    INTEGER(iwp) ::  file_id_out = 20             !< file unit for output VTK file
    INTEGER(iwp) ::  file_id_out_header = 19      !< file unit for temporary header file
    INTEGER(iwp) ::  length                       !< length of string on file
    INTEGER(iwp) ::  n                            !< running index over all points and polygons
    INTEGER(iwp) ::  npoints_total                !< total number of points
    INTEGER(iwp) ::  ns_total                     !< total number of polygons
    INTEGER(iwp) ::  num_pe                       !< number of processors of the run

!     INTEGER(OFFSET_KIND),DIMENSION(:), ALLOCATABLE ::  filepos !< current fileposition in binary file
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  filepos !< current fileposition in binary file

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  npoints !< number of points/vertices in a binaray file
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ns      !< number of surface elements in a binaray file

    LOGICAL ::  convert_average_data = .FALSE.          !< namelist parameter to decide whether average or instantaneous data should be converted

    LOGICAL, DIMENSION(:), ALLOCATABLE      ::  eof     !< flag indicating that end of binary file is reached

    REAL(wp)                              ::  simulated_time !< output time

    REAL(wp), DIMENSION(:), ALLOCATABLE   ::  var            !< actual surface data

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  points         !< point / vertex data
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  polygons       !< polygon data

!
!-- Read namelist.
    CALL surface_output_parin
!
!-- Allocate array which contains the file position in each output file,
!-- in order to skip some reading.
    ALLOCATE( eof(0:num_pe-1)     )
    ALLOCATE( filepos(0:num_pe-1) )
    ALLOCATE( npoints(0:num_pe-1) )
    ALLOCATE( ns(0:num_pe-1)      )
!
!-- Initialize file position.
    filepos = 0
!
!-- Open a temporary file which contains all necessary header information for the
!-- VTK format.
    OPEN ( file_id_out_header, FILE = 'HEADER', STATUS = 'REPLACE', FORM = 'FORMATTED' )
!
!-- READ grid setup, i.e. the number and position of vertices and surface elements
!-- and merge all this information into one file. Further, create all the required
!-- header information of the VTK file.
!-- Note, PARAVIEW expects one VTK file for each variable and each timestep.
!-- Hence, header information needs to be duplicated multiple times.
!-- Moreover, Paraview expects consecutive vertex and polygon data, which are
!-- all distributed over the binaray files. Hence, first read vertex data from
!-- binary file, write this to the HEADER file, close the binary file and read
!-- data from the next binary file and so on. This requires several openenings
!-- and closings of the binary files and temporarily storage of the
!-- file-positions.
    DO  f = 0, num_pe - 1
!
!--    Create filename of the treated binary file.
       CALL surface_output_create_file_string
!
!--    Open file with surface output for processor f.
       OPEN ( file_id_in, FILE = TRIM( path_in ) // TRIM( myid_char ), FORM = 'UNFORMATTED' )
!
!--    Read number of vertices / points and surface elements
       READ ( file_id_in )  npoints(f)
       READ ( file_id_in )  npoints_total
       READ ( file_id_in )  ns(f)
       READ ( file_id_in )  ns_total
!
!--    Allocate arrays where all the surfaces and vertices will be stored.
       ALLOCATE( points(3,1:npoints(f))   )
!
!--    Read polygon data and store them in a temporary file.
       READ ( file_id_in )  points
!
!--    Obtain current file position. Will be stored for next file opening.
!        CALL posix_ftell( file_id_in, filepos(f) )
       filepos(f) = INT( ftell( file_id_in ), KIND = iwp )
!
!--    Write header information. Only one time required.
       IF ( f == 0 )  THEN
          WRITE( file_id_out_header,'(A)' ) "# vtk DataFile Version 3.0"
          WRITE( file_id_out_header,'(A,F8.2,A)' ) "legacy vtk File generated by PALM,  simulation time = xxx sec"
          WRITE( file_id_out_header,'(A)') "ASCII"

          WRITE( file_id_out_header,'(A)') "DATASET POLYDATA"
          WRITE( file_id_out_header,'(A,I12,A)') "POINTS ", npoints_total, " float"
       ENDIF
!
!--    Write the vertex data into header file.
       DO  n = 1, npoints(f)
          WRITE( file_id_out_header, '(8F15.4)' )  points(1:3,n)
       ENDDO
!
!--    Deallocate vertex data and close binary file.
       DEALLOCATE( points )

       CLOSE ( file_id_in )
    ENDDO
!
!-- Now, treat polygon data.
    DO  f = 0, num_pe - 1
!
!--    Create filename of the treated binary file .
       CALL surface_output_create_file_string
!
!--    Open file with surface output for processor f.
       OPEN ( file_id_in, FILE = TRIM( path_in ) //  TRIM( myid_char ), FORM = 'UNFORMATTED' )
!
!--    Move to last postion.
!        CALL posix_lseek( file_id_in, filepos(f) )
       CALL FSEEK( file_id_in, filepos(f), 0, ndum )
!
!--    Allocate array for polygon data
       ALLOCATE( polygons(5,1:ns(f)) )
!
!--    Read polygon data and store them in a temporary file.
       READ ( file_id_in )  polygons
!
!--    Obtain current file position after reading the local polygon data.
!--    Will be used for next file opening.
!        CALL posix_ftell( file_id_in, filepos(f) )
       filepos(f) = INT( ftell( file_id_in ), KIND = iwp )
!
!--    Write further header information. Only one time required.
       IF ( f == 0 )  WRITE ( file_id_out_header, '(A,2I10)') "POLYGONS ", ns_total, 5 * ns_total
!
!--    Write the polygons into the header file. Note, polygon data is of type
!--    integer, as it just connects the point-indices which describe the given
!--    surface element.
       DO n = 1, ns(f)
          WRITE ( file_id_out_header, '(5I18)' )  INT( polygons(1:5,n) )
       ENDDO
!
!--    Deallocate array for polygon data and close the file.
       DEALLOCATE( polygons )

       CLOSE ( file_id_in )

    ENDDO

    f = 0
    CALL surface_output_create_file_string
!
!-- Write further header information. Only once required.
    WRITE ( file_id_out_header, '(A,I10)') "CELL_DATA ", ns_total
    WRITE ( file_id_out_header, '(A,I10)') "SCALARS cell_scalars float 1 "
    WRITE ( file_id_out_header, '(A,I10)') "LOOKUP_TABLE default "
!
!-- Header creation has now been completed. File can be closed.
    CLOSE( file_id_out_header )
!
!-- Now, read all the actual data from the surface output. Please note, Paraview
!-- VTK format expects 1 file per variable and time step.
!-- The output file is only created once and includes the variable and the
!-- simulated time.
!-- In the binaray files, several variables and timesteps are stored, data for a
!-- given variable, however, is distributed over all binary files. Hence, read
!-- variable data for a given timestep from the binary file, write this data into
!-- the target output file, remember file position in binary file and close it,
!-- open nex binary file and read the variable, and so on, until all variables
!-- for all timesteps are processed.
    eof = .FALSE.
    DO
       DO  f = 0, num_pe - 1
!
!--       Clean up strings-
          char_time = ''
          variable_name = ''
!
!--       Create filename of the treated binary file.
          CALL surface_output_create_file_string
!
!--       Open binary file with surface output for processor f.
          OPEN ( file_id_in, FILE = TRIM( path_in ) // TRIM( myid_char ), FORM = 'UNFORMATTED' )
!
!--       Move to last postion.
!           CALL posix_lseek( file_id_in, filepos(f) )
          CALL FSEEK( file_id_in, filepos(f), 0, ndum )
!
!--       Read string length and string indicating the output time step.
          READ ( file_id_in ) length
          READ ( file_id_in ) char_time(1:length)
!
!--       If string for the output time indicates the end-of-file, set the eof
!--       flag and skip the read of the loop.
          IF ( char_time(1:length) == 'END' )  THEN
             eof(f) = .TRUE.
             CLOSE ( file_id_in )
             CYCLE
          ENDIF
!
!--       Read output time, and variable name.
          READ ( file_id_in ) simulated_time
          READ ( file_id_in ) length
          READ ( file_id_in ) variable_name(1:length)
!
!--       For first loop index, open the target output file. First create the
!--       filename string. Further, copy HEADER file with the given filename
!--       string. The header information must be given in each VTK file!
          IF ( f == 0 )  THEN
             WRITE( char_dum, '(I9.0)' )  INT( simulated_time )
             print*, char_dum, TRIM(ADJUSTL(char_dum))
!
!--          Copy HEADER file and open VTK file
             IF ( convert_average_data )  THEN
                CALL system('cp HEADER ' // TRIM( path_out ) // TRIM( ADJUSTL(char_dum) ) //       &
                            '_AV_' // 's_' // TRIM( variable_name ) // '.vtk' )

                OPEN ( file_id_out, FILE = TRIM( path_out ) // TRIM( ADJUSTL(char_dum) ) //        &
                       '_AV_' // 's_' // TRIM( variable_name ) // '.vtk',                          &
                       FORM='FORMATTED', POSITION = 'APPEND' )
             ELSE
                CALL system('cp HEADER ' // TRIM( path_out ) // TRIM(ADJUSTL(char_dum)) //         &
                            's_' // TRIM( variable_name ) // '.vtk' )

                OPEN ( file_id_out, FILE = TRIM( path_out ) // TRIM( ADJUSTL(char_dum) ) //        &
                       's_' // TRIM( variable_name ) // '.vtk',                                    &
                       FORM='FORMATTED', POSITION = 'APPEND' )
             ENDIF

          ENDIF
!
!--       Allocate and read array for variable data.
          ALLOCATE( var(1:ns(f)) )

          READ( file_id_in ) var
!
!--       Write variable data into output VTK file.
          DO n = 1, ns(f)
             WRITE( file_id_out, * ) var(n)
          ENDDO
!
!--       Remember file position in binary file and close it.
!           CALL posix_ftell( file_id_in, filepos(f) )
          filepos(f) = INT( ftell( file_id_in ), KIND = iwp )

          CLOSE ( file_id_in )
!
!--       Deallocate temporary array for variable data.
          DEALLOCATE( var )

       ENDDO
!
!--    After data for a variable for one time step is read, close the output
!--    VTK file and go to next variable or timestep.
       CLOSE ( file_id_out )
!
!--    If all files reached the end-of-file, exit the loop.
       IF ( ALL( eof ) )  EXIT

    ENDDO
!
!-- Finally, remove HEADER file
    CALL system( 'rm HEADER' )

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine read the namelist file.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_output_parin

       IMPLICIT NONE

       INTEGER(iwp) ::  file_id_parin = 90

       NAMELIST /surface_output/  convert_average_data, num_pe, path_in, path_out

!
!--    Open namelist file.
       OPEN( file_id_parin, FILE='../namelist', STATUS='OLD', FORM='FORMATTED')
!
!--    Read namelist.
       READ ( file_id_parin, surface_output )
!
!--    Close namelist file.
       CLOSE( file_id_parin )

    END SUBROUTINE surface_output_parin

!------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine creates the filename string of the treated binary file.
!------------------------------------------------------------------------------!
    SUBROUTINE surface_output_create_file_string

       IMPLICIT NONE

       CHARACTER(LEN=3) ::  char_av = ''
!
!--    Create substring for averaged data.
       IF ( convert_average_data )  char_av = '_av'
!
!--    Create substring for the processor id and combine all substrings.
       WRITE( char_dum, '(I6.6)') f
       myid_char = TRIM( char_av ) // '_' // TRIM( char_dum )

    END SUBROUTINE surface_output_create_file_string

 END PROGRAM surface_output_to_vtk
