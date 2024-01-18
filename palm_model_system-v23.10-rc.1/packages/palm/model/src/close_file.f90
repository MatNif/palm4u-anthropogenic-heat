!> @file close_file.f90
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
!> Close specified file or all open files, if "0" has been given as the calling argument. In that
!> case, execute last actions for certain unit numbers, if required.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE close_file( file_id )

    USE control_parameters,                                                                        &
        ONLY:  debug_output,                                                                       &
               debug_string,                                                                       &
               max_masks,                                                                          &
               openfile,                                                                           &
               output_3d_file_size

    USE kinds

#if defined( __parallel )
    USE MPI
#endif

#if defined( __netcdf )
    USE NETCDF
#endif

    USE netcdf_interface,                                                                          &
        ONLY:  id_set_mask,                                                                        &
               id_set_pr,                                                                          &
               id_set_sp,                                                                          &
               id_set_ts,                                                                          &
               id_set_xy,                                                                          &
               id_set_xz,                                                                          &
               id_set_yz,                                                                          &
               id_set_3d,                                                                          &
               id_set_fl,                                                                          &
               nc_stat,                                                                            &
               netcdf_data_format,                                                                 &
               netcdf_handle_error

    USE pegrid,                                                                                    &
        ONLY:  myid
#if defined( __parallel )
    USE pegrid,                                                                                    &
        ONLY:  comm2d,                                                                             &
               ierr
#endif

    IMPLICIT NONE

    CHARACTER (LEN=10)  ::  datform = 'lit_endian' !<
    CHARACTER (LEN=80)  ::  title                  !<

    INTEGER(iwp) ::  av           !<
    INTEGER(iwp) ::  dimx         !<
    INTEGER(iwp) ::  dimy         !<
    INTEGER(iwp) ::  fid          !<
    INTEGER(iwp) ::  file_id      !<
    INTEGER(iwp) ::  mid          !< masked output running index
    INTEGER(iwp) ::  planz        !<

    LOGICAL ::  checkuf = .TRUE.  !<
    LOGICAL ::  datleg = .TRUE.   !<
    LOGICAL ::  dbp = .FALSE.     !<

    NAMELIST /GLOBAL/  checkuf, datform, dbp, dimx, dimy, planz, title
    NAMELIST /RAHMEN/  datleg

!
!-- Close specified unit number (if opened) and set a flag that it has been opened one time at least
    IF ( file_id /= 0 )  THEN
       IF ( openfile(file_id)%opened )  THEN
          CLOSE ( file_id )
          openfile(file_id)%opened        = .FALSE.
          openfile(file_id)%opened_before = .TRUE.
       ENDIF
       RETURN
    ENDIF

!
!-- Close all open unit numbers
    DO  fid = 1, 200+2*max_masks

       IF ( openfile(fid)%opened .OR. openfile(fid)%opened_before )  THEN
!
!--       Last actions for certain unit numbers
          SELECT CASE ( fid )

#if defined( __netcdf )
             CASE ( 101 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_xy(0) )
                   CALL netcdf_handle_error( 'close_file', 44 )
                ENDIF

             CASE ( 102 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_xz(0) )
                   CALL netcdf_handle_error( 'close_file', 45 )
                ENDIF

             CASE ( 103 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_yz(0) )
                   CALL netcdf_handle_error( 'close_file', 46 )
                ENDIF

             CASE ( 104 )

                IF ( myid == 0 )  THEN
                   nc_stat = NF90_CLOSE( id_set_pr )
                   CALL netcdf_handle_error( 'close_file', 47 )
                ENDIF

             CASE ( 105 )

                IF ( myid == 0 )  THEN
                   nc_stat = NF90_CLOSE( id_set_ts )
                   CALL netcdf_handle_error( 'close_file', 48 )
                ENDIF

             CASE ( 106 )
!
!--             Calculate the file size (to be output to the CPU log).
!--             Only makes sense for parallel I/O on all cores. For other cases 3d-data output is
!--             done in Fortran binary format to separate files (one per core) and gathered to
!--             a NetCDF file in post-processing.
                IF ( netcdf_data_format > 4 )  THEN
#if defined( __parallel )
                   CALL MPI_ALLREDUCE( MPI_IN_PLACE, output_3d_file_size, 1, MPI_REAL, MPI_SUM,    &
                                       comm2d, ierr )
#endif
!
!--                Default precision for output is 4 byte. Converted to Mbyte.
                   output_3d_file_size = output_3d_file_size * 4.0_wp / 1.0E6
                ENDIF

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_3d(0) )
                   CALL netcdf_handle_error( 'close_file', 49 )
                ENDIF

             CASE ( 107 )

                IF ( myid == 0 )  THEN
                   nc_stat = NF90_CLOSE( id_set_sp )
                   CALL netcdf_handle_error( 'close_file', 50 )
                ENDIF

             CASE ( 111 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_xy(1) )
                   CALL netcdf_handle_error( 'close_file', 52 )
                ENDIF

             CASE ( 112 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_xz(1) )
                   CALL netcdf_handle_error( 'close_file', 352 )
                ENDIF

             CASE ( 113 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_yz(1) )
                   CALL netcdf_handle_error( 'close_file', 353 )
                ENDIF

             CASE ( 116 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_3d(1) )
                   CALL netcdf_handle_error( 'close_file', 353 )
                ENDIF

             CASE ( 199 )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
                   nc_stat = NF90_CLOSE( id_set_fl )
                   CALL netcdf_handle_error( 'close_file', 353 )
                ENDIF

             CASE ( 201:200+2*max_masks )

                IF ( myid == 0  .OR.  netcdf_data_format > 4 )  THEN
!
!--                Decompose fid into mid and av
                   IF ( fid <= 200+max_masks )  THEN
                      mid = fid - 200
                      av = 0
                   ELSE
                      mid = fid - (200+max_masks)
                      av = 1
                   ENDIF
                   nc_stat = NF90_CLOSE( id_set_mask(mid,av) )
                   CALL netcdf_handle_error( 'close_file', 459 )

                ENDIF

#endif

          END SELECT
!
!--       Close file
          IF ( openfile(fid)%opened )  THEN
             IF ( debug_output )  THEN
                WRITE( debug_string, '(A,1X,I3)' )  'closing file id', fid
                CALL debug_message( debug_string, 'start' )
             ENDIF
             CLOSE ( fid )
             IF ( debug_output )  CALL debug_message( debug_string, 'end' )
          ENDIF

       ENDIF

    ENDDO

 END SUBROUTINE close_file
