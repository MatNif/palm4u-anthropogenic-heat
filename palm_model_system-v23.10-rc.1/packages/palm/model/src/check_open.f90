!> @file check_open.f90
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
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Check if file unit is open. If not, open file and, if necessary, write a
!> header or start other initializing actions, respectively.
!--------------------------------------------------------------------------------------------------!
SUBROUTINE check_open( file_id )

    USE control_parameters,                                                                        &
        ONLY:  coupling_char,                                                                      &
               data_output_2d_on_each_pe,                                                          &
               debug_output,                                                                       &
               debug_string,                                                                       &
               max_masks,                                                                          &
               message_string,                                                                     &
               openfile,                                                                           &
               run_description_header

    USE cpulog,                                                                                    &
        ONLY: cpu_log,                                                                             &
              log_point_s

#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  nz_do3d
#endif

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxr,                                                                                &
               ny,                                                                                 &
               nyn,                                                                                &
               nys,                                                                                &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt

    USE kinds

#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                                          &
        ONLY:  id_set_agt,                                                                         &
               id_set_fl,                                                                          &
               id_set_mask,                                                                        &
               id_set_pr,                                                                          &
               id_set_sp,                                                                          &
               id_set_ts,                                                                          &
               id_set_xy,                                                                          &
               id_set_xz,                                                                          &
               id_set_yz,                                                                          &
               id_set_3d,                                                                          &
               nc_stat,                                                                            &
               netcdf_create_file,                                                                 &
               netcdf_data_format,                                                                 &
               netcdf_define_header,                                                               &
               netcdf_handle_error,                                                                &
               netcdf_open_write_file

    USE particle_attributes,                                                                       &
        ONLY:  max_number_of_particle_groups,                                                      &
               number_of_particle_groups,                                                          &
               particle_groups

    USE pegrid

    USE posix_interface,                                                                           &
        ONLY:  fortran_sleep


    IMPLICIT NONE

    CHARACTER (LEN=30)  ::  filename                !<
    CHARACTER (LEN=5)   ::  mask_char               !<
    CHARACTER (LEN=80)  ::  rtext                   !<

    INTEGER(iwp) ::  av          !<
    INTEGER(iwp) ::  file_id     !<
    INTEGER(iwp) ::  ioerr       !< IOSTAT flag for IO-commands ( 0 = no error )
    INTEGER(iwp) ::  mid         !< masked output running index

    LOGICAL ::  netcdf_extend    !<

!
!-- Immediate return if file already open
    IF ( openfile(file_id)%opened )  RETURN

!
!-- Only certain files are allowed to be re-opened
!-- NOTE: some of the other files perhaps also could be re-opened, but it has not been checked so
!-- far, if it works!
    IF ( openfile(file_id)%opened_before )  THEN
       SELECT CASE ( file_id )
          CASE ( 13, 14, 21, 22, 23, 80, 85, 117 )
             IF ( file_id == 14  .AND.  openfile(file_id)%opened_before )  THEN
                message_string = 're-opening of unit ' // '14 not verified'
                CALL message( 'check_open', 'PAC0006', 0, 1, 0, 6, 0 )
             ENDIF

          CASE DEFAULT
             WRITE( message_string, * ) 're-opening of file-id ', file_id, ' is not allowed'
             CALL message( 'check_open', 'PAC0007', 0, 1, 0, 6, 0 )

             RETURN

       END SELECT
    ENDIF

!
!-- Check if file may be opened on the relevant PE
    SELECT CASE ( file_id )

       CASE ( 15, 16, 17, 18, 19, 50:59, 104:105, 107, 109, 117 )

          IF ( myid /= 0 )  THEN
             WRITE( message_string, * ) 'opening file-id ', file_id, ' not allowed for PE ', myid
             CALL message( 'check_open', 'PAC0008', 2, 2, -1, 6, 1 )
          ENDIF

       CASE ( 101:103, 106, 111:113, 116, 201:200+2*max_masks )

          IF ( netcdf_data_format < 5 )  THEN

             IF ( myid /= 0 )  THEN
                WRITE( message_string, * ) 'opening file-id ', file_id, ' not allowed for PE ', myid
                CALL message( 'check_open', 'PAC0008', 2, 2, -1, 6, 1 )
             ENDIF

          ENDIF

       CASE ( 21, 22, 23 )

          IF ( .NOT. data_output_2d_on_each_pe )  THEN
             IF ( myid /= 0 )  THEN
                WRITE( message_string, * ) 'opening file-id ', file_id, ' not allowed for PE ', myid
                CALL message( 'check_open', 'PAC0008', 2, 2, -1, 6, 1 )
             END IF
          ENDIF

       CASE ( 90:99 )

!
!--       File-ids that are used temporarily in other routines
          WRITE( message_string, * ) 'opening file-id ', file_id,                                  &
                                     ' is not allowed since it is used otherwise'
          CALL message( 'check_open', 'PAC0009', 0, 1, 0, 6, 0 )

    END SELECT

!
!-- Open relevant files
    SELECT CASE ( file_id )

       CASE ( 11 )
!
!--       Read the namelist parameter file.
          OPEN ( 11, FILE='PARIN' // TRIM( coupling_char ), FORM='FORMATTED', STATUS='OLD',        &
                     IOSTAT=ioerr )

          IF ( ioerr /= 0 )  THEN
             message_string = 'namelist file "PARIN' // TRIM( coupling_char ) // '" not found'
             CALL message( 'check_open', 'PAC0010', 3, 2, 0, 6, 1 )
          ENDIF

       CASE ( 13 )

          IF ( myid_char == '' )  THEN
             OPEN ( 13, FILE='BININ' // TRIM( coupling_char ) // myid_char, FORM='UNFORMATTED',    &
                        STATUS='OLD' )
          ELSE
!
!--          First opening of unit 13 openes file _000000 on all PEs because only this file contains
!--          the global variables.
             IF ( .NOT. openfile(file_id)%opened_before )  THEN
                OPEN ( 13, FILE='BININ' // TRIM( coupling_char ) // '/_000000', FORM='UNFORMATTED',&
                           STATUS='OLD' )
             ELSE
                OPEN ( 13, FILE='BININ' // TRIM( coupling_char ) // '/' // myid_char,              &
                           FORM='UNFORMATTED', STATUS='OLD' )
             ENDIF
          ENDIF

       CASE ( 14 )

          IF ( myid_char == '' )  THEN
             OPEN ( 14, FILE='BINOUT' // TRIM( coupling_char ) // myid_char, FORM='UNFORMATTED',   &
                        POSITION='APPEND' )
          ELSE
             IF ( myid == 0  .AND. .NOT. openfile(file_id)%opened_before )  THEN
                CALL local_system( 'mkdir  BINOUT' // TRIM( coupling_char ) )
             ENDIF
#if defined( __parallel )
!
!--          Set a barrier in order to allow that all other processors in the directory created by
!--          PE0 can open their file
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             ioerr = 1
             DO WHILE ( ioerr /= 0 )
                OPEN ( 14, FILE='BINOUT' // TRIM( coupling_char )// '/' // myid_char,              &
                           FORM='UNFORMATTED', IOSTAT=ioerr )
                IF ( ioerr /= 0 )  THEN
                   IF ( debug_output )  THEN
                      debug_string = 'could not open "BINOUT' // TRIM(coupling_char) // '/' //     &
                                     myid_char // '"! Trying again in 1 sec.'
                      CALL debug_message( debug_string, 'info' )
                   ENDIF
                   CALL fortran_sleep( 1 )
                ENDIF
             ENDDO

          ENDIF

       CASE ( 15 )

          OPEN ( 15, FILE='RUN_CONTROL' // TRIM( coupling_char ), FORM='FORMATTED' )

       CASE ( 16 )

          OPEN ( 16, FILE='LIST_PROFIL' // TRIM( coupling_char ), FORM='FORMATTED' )

       CASE ( 17 )

          OPEN ( 17, FILE='LIST_PROFIL_1D' // TRIM( coupling_char ), FORM='FORMATTED' )

       CASE ( 18 )

          OPEN ( 18, FILE='CPU_MEASURES' // TRIM( coupling_char ), FORM='FORMATTED' )

       CASE ( 19 )

          OPEN ( 19, FILE='HEADER' // TRIM( coupling_char ), FORM='FORMATTED' )

       CASE ( 20 )

          IF ( myid == 0  .AND. .NOT. openfile(file_id)%opened_before )  THEN
             CALL local_system( 'mkdir  DATA_LOG' // TRIM( coupling_char ) )
          ENDIF
          IF ( myid_char == '' )  THEN
             OPEN ( 20, FILE='DATA_LOG' // TRIM( coupling_char ) // '/_000000', FORM='UNFORMATTED',&
                        POSITION='APPEND' )
          ELSE
#if defined( __parallel )
!
!--          Set a barrier in order to allow that all other processors in the directory created by
!--          PE0 can open their file
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             ioerr = 1
             DO WHILE ( ioerr /= 0 )
                OPEN ( 20, FILE='DATA_LOG' // TRIM( coupling_char ) // '/' // myid_char,           &
                           FORM='UNFORMATTED', POSITION='APPEND', IOSTAT=ioerr )
                IF ( ioerr /= 0 )  THEN
                   IF ( debug_output )  THEN
                      debug_string = 'could not open "DATA_LOG' // TRIM(coupling_char) // '/' //   &
                                     myid_char // '"! Trying again in 1 sec.'
                      CALL debug_message( debug_string, 'info' )
                   ENDIF
                   CALL fortran_sleep( 1 )
                ENDIF
             ENDDO

          ENDIF

       CASE ( 21 )

          IF ( data_output_2d_on_each_pe )  THEN
             OPEN ( 21, FILE='PLOT2D_XY' // TRIM( coupling_char ) // myid_char, FORM='UNFORMATTED',&
                        POSITION='APPEND' )
          ELSE
             OPEN ( 21, FILE='PLOT2D_XY' // TRIM( coupling_char ), FORM='UNFORMATTED',             &
                        POSITION='APPEND' )
          ENDIF

          IF ( myid == 0  .AND.  .NOT. openfile(file_id)%opened_before )  THEN
!
!--          Write index bounds of total domain for combine_plot_fields
             IF ( data_output_2d_on_each_pe  .AND.  myid_char /= '' )  THEN
                WRITE (21)   0, nx,  0, ny
             ENDIF

          ENDIF

       CASE ( 22 )

          IF ( data_output_2d_on_each_pe )  THEN
             OPEN ( 22, FILE='PLOT2D_XZ' // TRIM( coupling_char ) // myid_char, FORM='UNFORMATTED',&
                        POSITION='APPEND' )
          ELSE
             OPEN ( 22, FILE='PLOT2D_XZ' // TRIM( coupling_char ), FORM='UNFORMATTED',             &
                        POSITION='APPEND' )
          ENDIF

          IF ( myid == 0  .AND.  .NOT. openfile(file_id)%opened_before )  THEN
!
!--          Write index bounds of total domain for combine_plot_fields
             IF ( data_output_2d_on_each_pe  .AND.  myid_char /= '' )  THEN
                WRITE (22)   0, nx, 0, nz+1    ! output part
             ENDIF

          ENDIF

       CASE ( 23 )

          IF ( data_output_2d_on_each_pe )  THEN
             OPEN ( 23, FILE='PLOT2D_YZ' // TRIM( coupling_char ) // myid_char, FORM='UNFORMATTED',&
                        POSITION='APPEND' )
          ELSE
             OPEN ( 23, FILE='PLOT2D_YZ' // TRIM( coupling_char ), FORM='UNFORMATTED',             &
                        POSITION='APPEND' )
          ENDIF

          IF ( myid == 0  .AND.  .NOT. openfile(file_id)%opened_before )  THEN
!
!--          Write index bounds of total domain for combine_plot_fields
             IF ( data_output_2d_on_each_pe  .AND.  myid_char /= '' )  THEN
                WRITE (23)   0, ny, 0, nz+1    ! output part
             ENDIF

          ENDIF

       CASE ( 25 )
!
!--       Binary files for surface data
          ! OPEN ( 25, FILE='SURFACE_DATA_BIN' // TRIM( coupling_char ) // myid_char,               &
          !            FORM='UNFORMATTED', POSITION='APPEND' )

          IF ( myid_char == '' )  THEN
             OPEN ( 25, FILE='SURFACE_DATA_BIN' // TRIM( coupling_char ) // myid_char,             &
                        FORM='UNFORMATTED', POSITION='APPEND' )
          ELSE
             IF ( myid == 0  .AND. .NOT. openfile(file_id)%opened_before )  THEN
                CALL local_system( 'mkdir  SURFACE_DATA_BIN' // TRIM( coupling_char ) )
             ENDIF
#if defined( __parallel )
!
!--          Set a barrier in order to allow that all other processors in the directory created by
!--          PE0 can open their file
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             ioerr = 1
             DO WHILE ( ioerr /= 0 )
                OPEN ( 25, FILE='SURFACE_DATA_BIN' // TRIM(coupling_char) //  '/' // myid_char,    &
                           FORM='UNFORMATTED', IOSTAT=ioerr )
                IF ( ioerr /= 0 )  THEN
                   IF ( debug_output )  THEN
                      debug_string = 'could not open "SURFACE_DATA_BIN' // TRIM(coupling_char) //  &
                                     '/' // myid_char // '"! Trying again in 1 sec.'
                      CALL debug_message( debug_string, 'info' )
                   ENDIF
                   CALL fortran_sleep( 1 )
                ENDIF
             ENDDO

          ENDIF

       CASE ( 26 )
!
!--       Binary files for averaged surface data
          ! OPEN ( 26, FILE='SURFACE_DATA_AV_BIN' // TRIM( coupling_char ) // myid_char,            &
          !        FORM='UNFORMATTED', POSITION='APPEND' )

          IF ( myid_char == '' )  THEN
             OPEN ( 26, FILE='SURFACE_DATA_AV_BIN' // TRIM( coupling_char ) // myid_char,          &
                        FORM='UNFORMATTED', POSITION='APPEND' )
          ELSE
             IF ( myid == 0  .AND. .NOT. openfile(file_id)%opened_before )  THEN
                CALL local_system( 'mkdir  SURFACE_DATA_AV_BIN' // TRIM( coupling_char ) )
             ENDIF
#if defined( __parallel )
!
!--          Set a barrier in order to allow that all other processors in the directory created by
!--          PE0 can open their file
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             ioerr = 1
             DO WHILE ( ioerr /= 0 )
                OPEN ( 26, FILE='SURFACE_DATA_AV_BIN' // TRIM( coupling_char ) // '/' // myid_char,&
                           FORM='UNFORMATTED', IOSTAT=ioerr )
                IF ( ioerr /= 0 )  THEN
                   IF ( debug_output )  THEN
                      debug_string = 'could not open "SURFACE_DATA_AV_BIN' // TRIM(coupling_char)  &
                                     // '/' // myid_char // '"! Trying again in 1 sec.'
                      CALL debug_message( debug_string, 'info' )
                   ENDIF
                   CALL fortran_sleep( 1 )
                ENDIF
             ENDDO

          ENDIF

       CASE ( 30 )

          OPEN ( 30, FILE='PLOT3D_DATA' // TRIM( coupling_char ) // myid_char, FORM='UNFORMATTED' )
!
!--       Specifications for combine_plot_fields
          IF ( myid == 0 )  THEN
#if defined( __parallel )
             WRITE ( 30 )  0, nx, 0, ny, 0, nz_do3d
#endif
          ENDIF

       CASE ( 80 )

          IF ( myid_char == '' )  THEN
             OPEN ( 80, FILE='PARTICLE_INFOS'//TRIM(coupling_char)//myid_char, FORM='FORMATTED',   &
                        POSITION='APPEND' )
          ELSE
             IF ( myid == 0  .AND.  .NOT. openfile(80)%opened_before )  THEN
                CALL local_system( 'mkdir  PARTICLE_INFOS' // TRIM( coupling_char ) )
             ENDIF
#if defined( __parallel )
!
!--          Set a barrier in order to allow that thereafter all other processors in the directory
!--          created by PE0 can open their file.
!--          WARNING: The following barrier will lead to hanging jobs, if check_open is first called
!--                   from routine allocate_prt_memory!
             IF ( .NOT. openfile(80)%opened_before )  THEN
                CALL MPI_BARRIER( comm2d, ierr )
             ENDIF
#endif
             OPEN ( 80, FILE='PARTICLE_INFOS' // TRIM( coupling_char ) // '/' // myid_char,        &
                        FORM='FORMATTED', POSITION='APPEND' )
          ENDIF

          IF ( .NOT. openfile(80)%opened_before )  THEN
             WRITE ( 80, 8000 )  TRIM( run_description_header )
          ENDIF

       CASE ( 85 )

          IF ( myid_char == '' )  THEN
             OPEN ( 85, FILE='PARTICLE_DATA' // TRIM(coupling_char) // myid_char,                  &
                        FORM='UNFORMATTED', POSITION='APPEND' )
          ELSE
             IF ( myid == 0  .AND.  .NOT. openfile(85)%opened_before )  THEN
                CALL local_system( 'mkdir  PARTICLE_DATA' // TRIM( coupling_char ) )
             ENDIF
#if defined( __parallel )
!
!--          Set a barrier in order to allow that thereafter all other processors in the directory
!--          created by PE0 can open their file
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             ioerr = 1
             DO WHILE ( ioerr /= 0 )
                OPEN ( 85, FILE='PARTICLE_DATA' // TRIM( coupling_char ) // '/' // myid_char,      &
                           FORM='UNFORMATTED', POSITION='APPEND', IOSTAT=ioerr )
                IF ( ioerr /= 0 )  THEN
                   IF ( debug_output )  THEN
                      debug_string = 'could not open "PARTICLE_DATA' // TRIM(coupling_char) //     &
                                     '/' // myid_char // '"! Trying again in 1 sec.'
                      CALL debug_message( debug_string, 'info' )
                   ENDIF
                   CALL fortran_sleep( 1 )
                ENDIF
             ENDDO

          ENDIF

          IF ( .NOT. openfile(85)%opened_before )  THEN
             WRITE ( 85 )  run_description_header
!
!--          Attention: change version number whenever the output format on unit 85 is changed (see
!--                     also in routine lpm_data_output_particles)
             rtext = 'data format version 3.1'
             WRITE ( 85 )  rtext
             WRITE ( 85 )  number_of_particle_groups, max_number_of_particle_groups
             WRITE ( 85 )  particle_groups
             WRITE ( 85 )  nxl, nxr, nys, nyn, nzb, nzt, nbgp
          ENDIF

!
!--    File where sky-view factors and further required data is stored will be read
       CASE ( 88 )

          IF ( myid_char == '' )  THEN
             OPEN ( 88, FILE='SVFIN' // TRIM( coupling_char ) // myid_char, FORM='UNFORMATTED',    &
                        STATUS='OLD', IOSTAT=ioerr )
          ELSE

             OPEN ( 88, FILE='SVFIN' // TRIM( coupling_char ) // '/' // myid_char,                 &
                        FORM='UNFORMATTED', STATUS='OLD', IOSTAT=ioerr )
          ENDIF

!
!--    File where sky-view factors and further required data is stored will be created
       CASE ( 89 )

          IF ( myid_char == '' )  THEN
             OPEN ( 89, FILE='SVFOUT' // TRIM( coupling_char ) // myid_char, FORM='UNFORMATTED',   &
                        STATUS='NEW' )
          ELSE
             IF ( myid == 0  .AND.  .NOT. openfile(file_id)%opened_before )  THEN
                CALL local_system( 'mkdir  SVFOUT' // TRIM( coupling_char ) )
             ENDIF
#if defined( __parallel )
!
!--          Set a barrier in order to allow that all other processors in the directory created by
!--          PE0 can open their file
             CALL MPI_BARRIER( comm2d, ierr )
#endif
             ioerr = 1
             DO WHILE ( ioerr /= 0 )
                OPEN ( 89, FILE='SVFOUT' // TRIM( coupling_char ) // '/' // myid_char,             &
                           FORM='UNFORMATTED', STATUS='NEW', IOSTAT=ioerr )
                IF ( ioerr /= 0 )  THEN
                   IF ( debug_output )  THEN
                      debug_string = 'could not open "SVFOUT' // TRIM(coupling_char) //            &
                                     '/' // myid_char // '"! Trying again in 1 sec.'
                      CALL debug_message( debug_string, 'info' )
                   ENDIF
                   CALL fortran_sleep( 1 )
                ENDIF
             ENDDO

          ENDIF

!
!--    Progress file that is used by the PALM watchdog
       CASE ( 117 )

          OPEN ( 117, FILE='PROGRESS' // TRIM( coupling_char ), STATUS='REPLACE', FORM='FORMATTED' )

#if defined( __netcdf )
       CASE ( 101, 111 )
!
!--       Set filename depending on unit number
          IF ( file_id == 101 )  THEN
             filename = 'DATA_2D_XY_NETCDF' // TRIM( coupling_char )
             av = 0
          ELSE
             filename = 'DATA_2D_XY_AV_NETCDF' // TRIM( coupling_char )
             av = 1
          ENDIF
!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its dimensions and variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )
          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_xy(av), .TRUE., 20 )
!
!--          Read header information and set all ids. If there is a mismatch between the previous
!--          and the actual run, netcdf_extend is returned as .FALSE.
             CALL netcdf_define_header( 'xy', netcdf_extend, av )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_xy(av) )
                CALL netcdf_handle_error( 'check_open', 21 )
                IF ( myid == 0 )  CALL local_system( 'rm ' // TRIM( filename ) )
#if defined( __parallel )
!
!--             Set a barrier in order to assure that PE0 deleted the old file before any other
!--             processor tries to open a new file.
!--             Barrier is only needed in case of parallel I/O
                IF ( netcdf_data_format > 4 )  CALL MPI_BARRIER( comm2d, ierr )
#endif
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_xy(av), .TRUE., 22 )

!
!--          Define the header
             CALL netcdf_define_header( 'xy', netcdf_extend, av )

!
!--          In case of parallel netCDF output, create flag file which tells combine_plot_fields
!--          that nothing is to do.
             IF ( myid == 0  .AND.  netcdf_data_format > 4 )  THEN
                OPEN( 99, FILE='NO_COMBINE_PLOT_FIELDS_XY' )
                WRITE ( 99, '(A)' )  'no combine_plot_fields.x neccessary'
                CLOSE( 99 )
             ENDIF

          ENDIF

       CASE ( 102, 112 )
!
!--       Set filename depending on unit number
          IF ( file_id == 102 )  THEN
             filename = 'DATA_2D_XZ_NETCDF' // TRIM( coupling_char )
             av = 0
          ELSE
             filename = 'DATA_2D_XZ_AV_NETCDF' // TRIM( coupling_char )
             av = 1
          ENDIF
!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its dimensions and variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_xz(av), .TRUE., 23 )
!
!--          Read header information and set all ids. If there is a mismatch between the previous
!--          and the actual run, netcdf_extend is returned as .FALSE.
             CALL netcdf_define_header( 'xz', netcdf_extend, av )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_xz(av) )
                CALL netcdf_handle_error( 'check_open', 24 )
                IF ( myid == 0 )  CALL local_system( 'rm ' // TRIM( filename ) )
#if defined( __parallel )
!
!--             Set a barrier in order to assure that PE0 deleted the old file before any other
!--             processor tries to open a new file.
!--             Barrier is only needed in case of parallel I/O
                IF ( netcdf_data_format > 4 )  CALL MPI_BARRIER( comm2d, ierr )
#endif
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_xz(av), .TRUE., 25 )

!
!--          Define the header
             CALL netcdf_define_header( 'xz', netcdf_extend, av )

!
!--          In case of parallel netCDF output, create flag file which tells combine_plot_fields
!--          that nothing is to do.
             IF ( myid == 0  .AND.  netcdf_data_format > 4 )  THEN
                OPEN( 99, FILE='NO_COMBINE_PLOT_FIELDS_XZ' )
                WRITE ( 99, '(A)' )  'no combine_plot_fields.x neccessary'
                CLOSE( 99 )
             ENDIF

          ENDIF

       CASE ( 103, 113 )
!
!--       Set filename depending on unit number
          IF ( file_id == 103 )  THEN
             filename = 'DATA_2D_YZ_NETCDF' // TRIM( coupling_char )
             av = 0
          ELSE
             filename = 'DATA_2D_YZ_AV_NETCDF' // TRIM( coupling_char )
             av = 1
          ENDIF
!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its dimensions and variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_yz(av), .TRUE., 26 )
!
!--          Read header information and set all ids. If there is a mismatch between the previous
!--          and the actual run, netcdf_extend is returned as .FALSE.
             CALL netcdf_define_header( 'yz', netcdf_extend, av )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_yz(av) )
                CALL netcdf_handle_error( 'check_open', 27 )
                IF ( myid == 0 )  CALL local_system( 'rm ' // TRIM( filename ) )
#if defined( __parallel )
!
!--             Set a barrier in order to assure that PE0 deleted the old file before any other
!--             processor tries to open a new file.
!--             Barrier is only needed in case of parallel I/O
                IF ( netcdf_data_format > 4 )  CALL MPI_BARRIER( comm2d, ierr )
#endif
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_yz(av), .TRUE., 28 )

!
!--          Define the header
             CALL netcdf_define_header( 'yz', netcdf_extend, av )

!
!--          In case of parallel netCDF output, create flag file which tells combine_plot_fields
!--          that nothing is to do.
             IF ( myid == 0  .AND.  netcdf_data_format > 4 )  THEN
                OPEN( 99, FILE='NO_COMBINE_PLOT_FIELDS_YZ' )
                WRITE ( 99, '(A)' )  'no combine_plot_fields.x neccessary'
                CLOSE( 99 )
             ENDIF

          ENDIF

       CASE ( 104 )
!
!--       Set filename
          filename = 'DATA_1D_PR_NETCDF' // TRIM( coupling_char )

!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_pr, .FALSE., 29 )
!
!--          Read header information and set all ids. If there is a mismatch between the previous
!--          and the actual run, netcdf_extend is returned as .FALSE.
             CALL netcdf_define_header( 'pr', netcdf_extend, 0 )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_pr )
                CALL netcdf_handle_error( 'check_open', 30 )
                CALL local_system( 'rm ' // TRIM( filename ) )
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_pr, .FALSE., 31 )
!
!--          Define the header
             CALL netcdf_define_header( 'pr', netcdf_extend, 0 )

          ENDIF

       CASE ( 105 )
!
!--       Set filename
          filename = 'DATA_1D_TS_NETCDF' // TRIM( coupling_char )

!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_ts, .FALSE., 32 )
!
!--          Read header information and set all ids. If there is a mismatch between the previous
!--          and the actual run, netcdf_extend is returned as .FALSE.
             CALL netcdf_define_header( 'ts', netcdf_extend, 0 )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_ts )
                CALL netcdf_handle_error( 'check_open', 33 )
                CALL local_system( 'rm ' // TRIM( filename ) )
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_ts, .FALSE., 34 )
!
!--          Define the header
             CALL netcdf_define_header( 'ts', netcdf_extend, 0 )

          ENDIF


       CASE ( 106, 116 )

          CALL cpu_log( log_point_s(96), 'output_3d open ', 'start' )
!
!--       Set filename depending on unit number
          IF ( file_id == 106 )  THEN
             filename = 'DATA_3D_NETCDF' // TRIM( coupling_char )
             av = 0
          ELSE
             filename = 'DATA_3D_AV_NETCDF' // TRIM( coupling_char )
             av = 1
          ENDIF
!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its dimensions and variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )
          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_3d(av), .TRUE., 35 )
!
!--          Read header information and set all ids. If there is a mismatch between the previous
!--          and the actual run, netcdf_extend is returned as .FALSE.
             CALL netcdf_define_header( '3d', netcdf_extend, av )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_3d(av) )
                CALL netcdf_handle_error( 'check_open', 36 )
                IF ( myid == 0 )  CALL local_system( 'rm ' // TRIM( filename ) )
#if defined( __parallel )
!
!--             Set a barrier in order to assure that PE0 deleted the old file before any other
!--             processor tries to open a new file.
!--             Barrier is only needed in case of parallel I/O
                IF ( netcdf_data_format > 4 )  CALL MPI_BARRIER( comm2d, ierr )
#endif
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_3d(av), .TRUE., 37 )

!
!--          Define the header
             CALL netcdf_define_header( '3d', netcdf_extend, av )

!
!--          In case of parallel netCDF output, create flag file which tells combine_plot_fields
!--          that nothing is to do.
             IF ( myid == 0  .AND.  netcdf_data_format > 4 )  THEN
                OPEN( 99, FILE='NO_COMBINE_PLOT_FIELDS_3D' )
                WRITE ( 99, '(A)' )  'no combine_plot_fields.x neccessary'
                CLOSE( 99 )
             ENDIF

          ENDIF

          CALL cpu_log( log_point_s(96), 'output_3d open ', 'stop' )

       CASE ( 107 )
!
!--       Set filename
          filename = 'DATA_1D_SP_NETCDF' // TRIM( coupling_char )

!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_sp, .FALSE., 38 )

!
!--          Read header information and set all ids. If there is a mismatch between the previous
!--          and the actual run, netcdf_extend is returned as .FALSE.
             CALL netcdf_define_header( 'sp', netcdf_extend, 0 )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_sp )
                CALL netcdf_handle_error( 'check_open', 39 )
                CALL local_system( 'rm ' // TRIM( filename ) )
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_sp, .FALSE., 40 )
!
!--          Define the header
             CALL netcdf_define_header( 'sp', netcdf_extend, 0 )

          ENDIF

       CASE ( 118 )

          IF ( myid == 0 )  THEN
             filename = 'DATA_AGT_NETCDF'
!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_agt, .FALSE., 43 )

!
!--          Define the header
             CALL netcdf_define_header( 'ag', netcdf_extend, 0 )
          ENDIF

!           IF ( netcdf_extend )  THEN
! !
! !--          Open an existing netCDF file for output
!              CALL netcdf_open_write_file( filename, id_set_agt, .FALSE., 41 )
! !
! !--          Read header information and set all ids. If there is a mismatch
! !--          between the previous and the actual run, netcdf_extend is returned
! !--          as .FALSE.
!              CALL netcdf_define_header( 'ag', netcdf_extend, 0 )
!
! !
! !--          Remove the local file, if it can not be extended
!              IF ( .NOT. netcdf_extend )  THEN
!                 nc_stat = NF90_CLOSE( id_set_agt )
!                 CALL netcdf_handle_error( 'check_open', 42 )
!                 CALL local_system( 'rm ' // TRIM( filename ) )
!              ENDIF
!
!           ENDIF

          IF ( .NOT. netcdf_extend )  THEN

!
! !--          For runs on multiple processors create the subdirectory
!              IF ( myid_char /= '' )  THEN
!                 IF ( myid == 0  .AND. .NOT. openfile(file_id)%opened_before )  &
!                 THEN    ! needs modification in case of non-extendable sets
!                    CALL local_system( 'mkdir  DATA_PRT_NETCDF' //              &
!                                        TRIM( coupling_char ) // '/' )
!                 ENDIF
! #if defined( __parallel )
! !
! !--             Set a barrier in order to allow that all other processors in the
! !--             directory created by PE0 can open their file
!                 CALL MPI_BARRIER( comm2d, ierr )
! #endif
!              ENDIF

          ENDIF


!
!--    nc-file for virtual flight measurements
       CASE ( 199 )
!
!--       Set filename
          filename = 'DATA_1D_FL_NETCDF' // TRIM( coupling_char )

!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_fl, .FALSE., 532 )
!
!--          Read header information and set all ids. If there is a mismatch between the previous
!--          and the actual run, netcdf_extend is returned as .FALSE.
             CALL netcdf_define_header( 'fl', netcdf_extend, 0 )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_fl )
                CALL netcdf_handle_error( 'check_open', 533 )
                CALL local_system( 'rm ' // TRIM( filename ) )
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_fl, .FALSE., 534 )
!
!--          Define the header
             CALL netcdf_define_header( 'fl', netcdf_extend, 0 )

          ENDIF


       CASE ( 201:200+2*max_masks )
!
!--       Set filename depending on unit number
          IF ( file_id <= 200+max_masks )  THEN
             mid = file_id - 200
             WRITE ( mask_char,'(A2,I3.3)')  '_M', mid
             filename = 'DATA_MASK_NETCDF' // TRIM( coupling_char ) // mask_char
             av = 0
          ELSE
             mid = file_id - (200+max_masks)
             WRITE ( mask_char,'(A2,I3.3)')  '_M', mid
             filename = 'DATA_MASK_AV_NETCDF' // TRIM( coupling_char ) // mask_char
             av = 1
          ENDIF
!
!--       Inquire, if there is a netCDF file from a previous run. This should be opened for
!--       extension, if its dimensions and variables match the actual run.
          INQUIRE( FILE=filename, EXIST=netcdf_extend )

          IF ( netcdf_extend )  THEN
!
!--          Open an existing netCDF file for output
             CALL netcdf_open_write_file( filename, id_set_mask(mid,av), .TRUE., 456 )
!
!--          Read header information and set all ids. If there is a mismatch between the previous
!--          and the actual run, netcdf_extend is returned as .FALSE.
             CALL netcdf_define_header( 'ma', netcdf_extend, file_id )

!
!--          Remove the local file, if it can not be extended
             IF ( .NOT. netcdf_extend )  THEN
                nc_stat = NF90_CLOSE( id_set_mask(mid,av) )
                CALL netcdf_handle_error( 'check_open', 457 )
                CALL local_system('rm ' // TRIM( filename ) )
             ENDIF

          ENDIF

          IF ( .NOT. netcdf_extend )  THEN
!
!--          Create a new netCDF output file with requested netCDF format
             CALL netcdf_create_file( filename, id_set_mask(mid,av), .TRUE. , 458 )
!
!--          Define the header
             CALL netcdf_define_header( 'ma', netcdf_extend, file_id )

          ENDIF


#else

       CASE ( 101:107, 111:113, 116, 201:200+2*max_masks )

!
!--       Nothing is done in case of missing netcdf support
          RETURN

#endif

       CASE DEFAULT

          WRITE( message_string, * ) 'no OPEN-statement for file-id ', file_id
          CALL message( 'check_open', 'PAC0011', 2, 2, -1, 6, 1 )

    END SELECT

!
!-- Set open flag
    openfile(file_id)%opened = .TRUE.

!
!-- Formats
8000 FORMAT (A/                                                                                    &
             '  step    time    # of parts     lPE sent/recv  rPE sent/recv  ',                    &
             'sPE sent/recv  nPE sent/recv    max # of parts  '/                                   &
             109('-'))

 END SUBROUTINE check_open
