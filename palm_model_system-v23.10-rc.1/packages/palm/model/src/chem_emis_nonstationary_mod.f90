!> @file chem_emis_nonstationary_mod.f90
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
! Copyright 2023 pecanode GmbH
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
!> @author Matthias Suehring
!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Modules for nonstationary volume-source emissions.
!--------------------------------------------------------------------------------------------------!
 MODULE chem_emis_nonstationary_mod

    USE arrays_3d,                                                                                 &
        ONLY:  zu,                                                                                 &
               zw

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  pi

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    USE chem_modules,                                                                              &
        ONLY:  chem_species,                                                                       &
               communicator_chem,                                                                  &
               emis_nonstationary_lod,                                                             &
               nc_field_length

    USE control_parameters,                                                                        &
        ONLY:  coupling_char,                                                                      &
               dt_3d,                                                                              &
               message_string,                                                                     &
               time_since_reference_point

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE general_utilities,                                                                         &
        ONLY:  interpolate_linear

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nxl,                                                                                &
               nxr,                                                                                &
               nys,                                                                                &
               nyn,                                                                                &
               nzt,                                                                                &
               topo_top_ind

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY :  get_attribute,                                                                     &
                get_dimension_length,                                                              &
                get_variable,                                                                      &
                init_model,                                                                        &
                open_read_file

    USE palm_date_time_mod,                                                                        &
        ONLY:  diff_in_sec_to_reference_date

    IMPLICIT NONE

    SAVE
    PRIVATE

    CHARACTER(LEN=17), PARAMETER ::  pids_emis_nonstat = 'PIDS_EMIS_NONSTAT' !< input file

    INTEGER(iwp) ::  num_emission_path     !< number of emission pathes
    INTEGER(iwp) ::  pids_emis_nonstat_id  !< NetCDF ID

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  e_utm_tmp ! temporary EUTM coordinate used for rotation
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  n_utm_tmp ! temporary NUTM coordinate used for rotation

    TYPE nonstat_emission

       CHARACTER(LEN=nc_field_length)                            ::  current_time  !< simulation time

       CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE             ::  emis_species  !< emission species
       CHARACTER(LEN=nc_field_length), DIMENSION(:), ALLOCATABLE ::  timestamp     !< individual timestamps of input data

       INTEGER(iwp) ::  num_emis_species !< number of treated species
       INTEGER(iwp) ::  num_timestamp    !< number of timestamps
       INTEGER(iwp) ::  num_vsrc_pos     !< number of volume source points
       INTEGER(iwp) ::  ts               !< index of current timestep
       INTEGER(iwp) ::  tsp              !< index of next timestep

       INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  emis_species_index  !< corresponding index in chem_species structure

       REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  eutm  !< time-dependent EUTM coordinate
       REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  nutm  !< time-dependent NUTM coordinate
       REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::  zag   !< time-dependent z (above ground) coordinate

       REAL(wp), DIMENSION(:), ALLOCATABLE  ::  timediff_ref_point !< time difference with respect to simulation start tme

       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  vsrc_emis_value   !< volume source values

    END TYPE nonstat_emission

    TYPE( nonstat_emission ), DIMENSION(:), ALLOCATABLE ::  emis_ns !< non-stationary emissions


    INTERFACE chem_emis_nonstationary_init
       MODULE PROCEDURE chem_emis_nonstationary_init
    END INTERFACE chem_emis_nonstationary_init

    INTERFACE chem_emis_nonstationary_update
       MODULE PROCEDURE chem_emis_nonstationary_update
    END INTERFACE chem_emis_nonstationary_update

!
!-- Public subroutines.
    PUBLIC :: chem_emis_nonstationary_init,                                                        &
              chem_emis_nonstationary_update

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of nonstationary time-dependent volume-source emissions.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_nonstationary_init

    CHARACTER(LEN=5) :: dum !< dummy string indicating emission path id

    CHARACTER(LEN=512), DIMENSION(:), ALLOCATABLE ::  buf_char !< buffer variable to read strings

    INTEGER(iwp) ::  char_length !< character length
    INTEGER(iwp) ::  k           !< running index over timestamps
    INTEGER(iwp) ::  l           !< running index over chemical species data structure
    INTEGER(iwp) ::  m           !< running index over emission species
    INTEGER(iwp) ::  n           !< running index over emission pathes
    INTEGER(iwp) ::  len_spec    !< character length of species name

    LOGICAL ::  emission_file_exist !< flag to check if emission input file exist

    CALL location_message( 'reading emissions data from ' // TRIM( pids_emis_nonstat ) //          &
                            TRIM( coupling_char ), 'start' )

!
!-- Check for emission input file.
    INQUIRE( FILE = TRIM( pids_emis_nonstat ) // TRIM( coupling_char ),                            &
             EXIST = emission_file_exist )
    IF ( .NOT. emission_file_exist )  THEN
       message_string = 'missing emission input file ( ' // TRIM( pids_emis_nonstat ) //           &
                        TRIM( coupling_char ) // ' for moving emission sources'
       CALL message( 'chem_emis_nonstationary_init', 'CHM0049', 2, 2, 0, 6, 0 )
    ENDIF
!
!-- Open emission file.
    CALL open_read_file( TRIM( pids_emis_nonstat ) // TRIM( coupling_char ), pids_emis_nonstat_id )
!
!-- Read LOD attribute.
    CALL get_attribute( pids_emis_nonstat_id, 'lod', emis_nonstationary_lod, .TRUE. )
!
!-- Read number of emission pathes.
    CALL get_attribute( pids_emis_nonstat_id, 'num_emission_path', num_emission_path, .TRUE. )
!
!-- Read character length - required to process input characters from NetCDF file.
    CALL get_dimension_length( pids_emis_nonstat_id, char_length, 'field_length' )
!
!-- Allocate size of data structure.
    ALLOCATE( emis_ns(1:num_emission_path) )
!
!-- For each emission path, read its dimension sizes (number of timestamps and number of volume
!-- source coordinates), allocate respective arrays sizes and read corresponding initial values.
    DO  n = 1, num_emission_path
!
!--    First of all, determine suffix ID, ordered according to the number of emission pathes.
       IF( n < 10 )  THEN
          WRITE( dum, '(I1)')  n
       ELSEIF( n < 100 )  THEN
          WRITE( dum, '(I2)')  n
       ELSEIF( n < 1000 )  THEN
          WRITE( dum, '(I3)')  n
       ELSEIF( n < 10000 )  THEN
          WRITE( dum, '(I4)')  n
       ELSEIF( n < 100000 )  THEN
          WRITE( dum, '(I5)')  n
       ENDIF

       CALL get_dimension_length( pids_emis_nonstat_id, emis_ns(n)%num_timestamp,                  &
                                  'ntime'    // TRIM( dum ) )
       CALL get_dimension_length( pids_emis_nonstat_id, emis_ns(n)%num_vsrc_pos,                   &
                                  'nvsrc'    // TRIM( dum ) )
       CALL get_dimension_length( pids_emis_nonstat_id, emis_ns(n)%num_emis_species,               &
                                  'nspecies' // TRIM( dum ) )
!
!--    Allocate memory for timestamps, read the data and convert these into model time. For
!--    data input, a buffer array is used, which is trimmed and copied later on. This
!--    is required because of the "special design" of the get_variable routine, which requires
!--    revision. Note, the variable buf_char is allocated in get_variable.
       ALLOCATE( emis_ns(n)%timestamp(0:emis_ns(n)%num_timestamp-1) )
       ALLOCATE( emis_ns(n)%timediff_ref_point(0:emis_ns(n)%num_timestamp-1) )

       CALL get_variable( pids_emis_nonstat_id, 'timestamp' // TRIM( dum ) , buf_char,             &
                          emis_ns(n)%num_timestamp, char_length )
!
!--    Compute the time difference to the actual model time (time_since_reference_point).
       DO  k = 0, emis_ns(n)%num_timestamp-1
          emis_ns(n)%timestamp(k) = buf_char(k+1)(1:char_length)   !buf allocated from 1

          emis_ns(n)%timediff_ref_point(k) = diff_in_sec_to_reference_date(                        &
                                                               TRIM( emis_ns(n)%timestamp(k) )     &
                                                                          )
       ENDDO
!
!--    Deallocate buf_char (allocated in get_variable).
       DEALLOCATE( buf_char )
!
!--    Compute corresponding time index.
       emis_ns(n)%ts = MINLOC( ABS( emis_ns(n)%timediff_ref_point - time_since_reference_point ),  &
                               DIM = 1 ) - 1
       IF ( emis_ns(n)%ts < emis_ns(n)%num_timestamp-1 )  THEN
          emis_ns(n)%tsp = emis_ns(n)%ts + 1
       ELSE
          emis_ns(n)%tsp = emis_ns(n)%ts
       ENDIF
!
!--    Allocate memory for emission species and read the strings.
       ALLOCATE( emis_ns(n)%emis_species(0:emis_ns(n)%num_emis_species-1) )
       CALL get_variable( pids_emis_nonstat_id, 'species' // TRIM( dum ), buf_char,                &
                          emis_ns(n)%num_emis_species, char_length )
!
!--    Store temporary string for species name on data structure.
       DO  k = 0, emis_ns(n)%num_emis_species-1
          emis_ns(n)%emis_species(k) = buf_char(k+1)(1:char_length)   !buf allocated from 1
       ENDDO
!
!--    Deallocate buf_char (allocated in get_variable).
       DEALLOCATE( buf_char )
!
!--    According to the KPP mechanism, assign each emission species a corresponding index with
!--    respect to prognostic variable.
       ALLOCATE( emis_ns(n)%emis_species_index(0:emis_ns(n)%num_emis_species-1) )
       DO  l = 1, nvar
          DO  m = 0, emis_ns(n)%num_emis_species-1

             len_spec = LEN( TRIM( chem_species(l)%name ) )

             IF ( TRIM( chem_species(l)%name ) == emis_ns(n)%emis_species(m)(1:len_spec) )  THEN
                emis_ns(n)%emis_species_index(m) = l
             ENDIF
          ENDDO
       ENDDO
!
!--    Allocate memory for the UTM and height coordinates and read them for the current and
!--    next timestamp.
       ALLOCATE( emis_ns(n)%eutm(0:1,0:emis_ns(n)%num_vsrc_pos-1) )
       ALLOCATE( emis_ns(n)%nutm(0:1,0:emis_ns(n)%num_vsrc_pos-1) )
       ALLOCATE( emis_ns(n)%zag(0:1,0:emis_ns(n)%num_vsrc_pos-1)  )

       emis_ns(n)%eutm = 0.0_wp
       emis_ns(n)%nutm = 0.0_wp
       emis_ns(n)%zag  = 0.0_wp

       CALL get_variable( pids_emis_nonstat_id, 'vsrc' // TRIM( dum ) // '_eutm', emis_ns(n)%eutm, &
                          0, emis_ns(n)%num_vsrc_pos-1, emis_ns(n)%ts, emis_ns(n)%tsp )
       CALL get_variable( pids_emis_nonstat_id, 'vsrc' // TRIM( dum ) // '_nutm', emis_ns(n)%nutm, &
                          0, emis_ns(n)%num_vsrc_pos-1, emis_ns(n)%ts, emis_ns(n)%tsp )
       CALL get_variable( pids_emis_nonstat_id, 'vsrc' // TRIM( dum ) // '_zag', emis_ns(n)%zag,   &
                          0, emis_ns(n)%num_vsrc_pos-1, emis_ns(n)%ts, emis_ns(n)%tsp )
!
!--    If necessary, rotate UTM coordinates.
       IF ( init_model%rotation_angle /= 0.0_wp )  THEN
          ALLOCATE( e_utm_tmp(0:1,0:emis_ns(n)%num_vsrc_pos-1) )
          ALLOCATE( n_utm_tmp(0:1,0:emis_ns(n)%num_vsrc_pos-1) )

          e_utm_tmp = emis_ns(n)%eutm - init_model%origin_x
          n_utm_tmp = emis_ns(n)%nutm - init_model%origin_y

          emis_ns(n)%eutm = COS( init_model%rotation_angle * pi / 180.0_wp ) * e_utm_tmp           &
                          - SIN( init_model%rotation_angle * pi / 180.0_wp ) * n_utm_tmp
          emis_ns(n)%nutm = SIN( init_model%rotation_angle * pi / 180.0_wp ) * e_utm_tmp           &
                          + COS( init_model%rotation_angle * pi / 180.0_wp ) * n_utm_tmp

          DEALLOCATE( e_utm_tmp )
          DEALLOCATE( n_utm_tmp )
       ENDIF
!
!--    Allocate memory for emission values and read data from input file.
       ALLOCATE( emis_ns(n)%vsrc_emis_value(0:1,0:emis_ns(n)%num_vsrc_pos-1,                       &
                                            0:emis_ns(n)%num_emis_species-1) )
       DO  m = 0, emis_ns(n)%num_emis_species-1
          CALL get_variable( pids_emis_nonstat_id, 'vsrc' // TRIM( dum ) // '_' //                 &
                                                  TRIM( emis_ns(n)%emis_species(m) ),              &
                             emis_ns(n)%vsrc_emis_value(:,:,m),                                    &
                             0, emis_ns(n)%num_vsrc_pos-1, emis_ns(n)%ts, emis_ns(n)%tsp )
       ENDDO
    ENDDO

    CALL location_message( 'reading emissions data from ' // TRIM( pids_emis_nonstat ) //          &
                            TRIM( coupling_char ), 'finished' )

 END SUBROUTINE chem_emis_nonstationary_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Updates nonstationary emissions at every time step.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_nonstationary_update

    CHARACTER(LEN=5) :: dum !< dummy string indicating emission path id

    INTEGER(iwp) ::  i        !< grid index in x-direction
    INTEGER(iwp) ::  ispecies !< correspoding index of emitted species in general chem data structure
    INTEGER(iwp) ::  j        !< grid index in y-direction
    INTEGER(iwp) ::  k        !< grid index in z-direction
    INTEGER(iwp) ::  l        !< running index over chemical species
    INTEGER(iwp) ::  n        !< running index over emission pathes
    INTEGER(iwp) ::  m        !< running index over considered emission species
    INTEGER(iwp) ::  p        !< running index over volume source coordinates
    INTEGER(iwp) ::  nts      !< number of available timestamps
    INTEGER(iwp) ::  tl       !< lower time interpolation index (always zero)
    INTEGER(iwp) ::  ts       !< index referring to current timestamp
    INTEGER(iwp) ::  tsp      !< index referring to next timestamp
    INTEGER(iwp) ::  tu       !< upper time interpolation index (zero or one)

    LOGICAL ::  emission_at_timestep !< flag to indicate whether emissions are necessary at current timestep
    LOGICAL ::  read_emissions       !< flag to indicate whether emission data need to be input

    REAL(wp) ::  fac_dt      !< interpolation factor
    REAL(wp) ::  vsrc_term   !< source term
    REAL(wp) ::  zag_tmp     !< height above ground

    CALL cpu_log( log_point_s(103), 'chem. emission nonstat', 'start' )

    DO  n = 1, num_emission_path
!
!--    Check whether volume-source emissions at current emission path are necessary at current
!--    timestep.
       nts = emis_ns(n)%num_timestamp-1
       emission_at_timestep = ( emis_ns(n)%timediff_ref_point(0)   <= time_since_reference_point   &
                         .AND.  emis_ns(n)%timediff_ref_point(nts) >= time_since_reference_point )

       IF ( .NOT. emission_at_timestep )  CYCLE
!
!--    Check whether emission data need to be input.
       ts  = emis_ns(n)%ts
       tsp = emis_ns(n)%tsp
       read_emissions = ( emis_ns(n)%timediff_ref_point(tsp) < time_since_reference_point + dt_3d )
!
!--    Read emission data.
       IF ( read_emissions )  THEN
!
!--       Determine suffix ID, ordered according to the number of emission pathes.
          IF( n < 10 )  THEN
             WRITE( dum, '(I1)')  n
          ELSEIF( n < 100 )  THEN
             WRITE( dum, '(I2)')  n
          ELSEIF( n < 1000 )  THEN
             WRITE( dum, '(I3)')  n
          ELSEIF( n < 10000 )  THEN
             WRITE( dum, '(I4)')  n
          ELSEIF( n < 100000 )  THEN
             WRITE( dum, '(I5)')  n
          ENDIF
!
!--       Update the corresponding time indices.
          emis_ns(n)%ts = MINLOC( ABS( emis_ns(n)%timediff_ref_point                               &
                                     - time_since_reference_point ), DIM = 1 ) - 1

          IF ( emis_ns(n)%ts < emis_ns(n)%num_timestamp-1 )  THEN
             emis_ns(n)%tsp = emis_ns(n)%ts + 1
          ELSE
             emis_ns(n)%tsp = emis_ns(n)%ts
          ENDIF
!
!--       Update coordinates.
          CALL get_variable( pids_emis_nonstat_id, 'vsrc' // TRIM( dum ) // '_eutm',               &
                             emis_ns(n)%eutm, 0, emis_ns(n)%num_vsrc_pos-1,                        &
                             emis_ns(n)%ts, emis_ns(n)%tsp )
          CALL get_variable( pids_emis_nonstat_id, 'vsrc' // TRIM( dum ) // '_nutm',               &
                             emis_ns(n)%nutm, 0, emis_ns(n)%num_vsrc_pos-1,                        &
                             emis_ns(n)%ts, emis_ns(n)%tsp )
          CALL get_variable( pids_emis_nonstat_id, 'vsrc' // TRIM( dum ) // '_zag',                &
                             emis_ns(n)%zag, 0, emis_ns(n)%num_vsrc_pos-1,                         &
                             emis_ns(n)%ts, emis_ns(n)%tsp )
!
!--       If necessary, rotate UTM coordinates.
          IF ( init_model%rotation_angle /= 0.0_wp )  THEN
             ALLOCATE( e_utm_tmp(0:1,0:emis_ns(n)%num_vsrc_pos-1) )
             ALLOCATE( n_utm_tmp(0:1,0:emis_ns(n)%num_vsrc_pos-1) )

             e_utm_tmp = emis_ns(n)%eutm - init_model%origin_x
             n_utm_tmp = emis_ns(n)%nutm - init_model%origin_y

             emis_ns(n)%eutm = COS( init_model%rotation_angle * pi / 180.0_wp ) * e_utm_tmp        &
                             - SIN( init_model%rotation_angle * pi / 180.0_wp ) * n_utm_tmp
             emis_ns(n)%nutm = SIN( init_model%rotation_angle * pi / 180.0_wp ) * e_utm_tmp        &
                             + COS( init_model%rotation_angle * pi / 180.0_wp ) * n_utm_tmp

             DEALLOCATE( e_utm_tmp )
             DEALLOCATE( n_utm_tmp )
          ENDIF
!
!--       Read correspoinding data at time indices ts and tsp.
          DO  m = 0, emis_ns(n)%num_emis_species-1
             CALL get_variable( pids_emis_nonstat_id, 'vsrc' // TRIM( dum ) // '_' //              &
                                                     TRIM( emis_ns(n)%emis_species(m) ),           &
                                emis_ns(n)%vsrc_emis_value(:,:,m),                                 &
                                0, emis_ns(n)%num_vsrc_pos-1, emis_ns(n)%ts, emis_ns(n)%tsp )
          ENDDO

       ENDIF
!
!--    Determine interpolation factor to linearly interpolate in time between two coordinate- and
!--    emission values.
       IF ( emis_ns(n)%timediff_ref_point(emis_ns(n)%tsp)                                          &
          - emis_ns(n)%timediff_ref_point(emis_ns(n)%ts) > 0.0_wp )  THEN

          fac_dt = ( time_since_reference_point + dt_3d                                            &
                   - emis_ns(n)%timediff_ref_point(emis_ns(n)%ts)                                  &
                   ) /                                                                             &
                   ( emis_ns(n)%timediff_ref_point(emis_ns(n)%tsp)                                 &
                   - emis_ns(n)%timediff_ref_point(emis_ns(n)%ts)                                  &
                   )

       ELSE
          fac_dt = 0.0_wp
       ENDIF

!
!--    Impose volume-source term for each species. First loop over all species.
       DO  m = 0, emis_ns(n)%num_emis_species-1
          ispecies = emis_ns(n)%emis_species_index(m)
!
!--       Loop over all volume-source points.
          DO  p = 0, emis_ns(n)%num_vsrc_pos-1
!
!--          Determine lower and upper interpolation indices. In case of a restart, the upper
!--          interpolation case might not be read and could arbitrary values. Hence, determine
!--          the interpolation indices from the time indices. In the described case, the
!--          upper interpolation index will have the same value as the lower one.
             tl = emis_ns(n)%ts  - emis_ns(n)%ts
             tu = emis_ns(n)%tsp - emis_ns(n)%ts
!
!--          Determine horizontal grid indices. Linearly interpolate in between the two emissions
!--          timesteps.
             i = INT( interpolate_linear( emis_ns(n)%eutm(tl,p), emis_ns(n)%eutm(tu,p), fac_dt )   &
                    * ddx )
             j = INT( interpolate_linear( emis_ns(n)%nutm(tl,p), emis_ns(n)%nutm(tu,p), fac_dt )   &
                    * ddy )

             IF ( i >= nxl  .AND.  i <= nxr  .AND.  j >= nys  .AND.  j<= nyn )  THEN
!
!--             Determine vertical grid index.
                zag_tmp = interpolate_linear( emis_ns(n)%zag(0,p), emis_ns(n)%zag(1,p), fac_dt )
                k = MINLOC( ABS( zu - zw(topo_top_ind(j,i,0)) - zag_tmp ), DIM = 1 ) - 1
!
!--             Make sure that emissions are not within topography cells.
                k = MAX( k, topo_top_ind(j,i,0) + 1 )
!
!--             Any conversion required?
                vsrc_term = interpolate_linear( emis_ns(n)%vsrc_emis_value(tl,p,m),                &
                                                emis_ns(n)%vsrc_emis_value(tu,p,m),                &
                                                fac_dt )

                chem_species(ispecies)%conc(k,j,i) = chem_species(ispecies)%conc(k,j,i)            &
                                                   + dt_3d * vsrc_term
             ENDIF
          ENDDO
       ENDDO
    ENDDO
!
!-- Exchange ghost points for the corresponding species. Necessary when not called from prognostic
!-- equations.
    DO  l = 1, nvar
       CALL exchange_horiz( chem_species(l)%conc, nbgp,                                            &
                            alternative_communicator = communicator_chem )
    ENDDO

    CALL cpu_log( log_point_s(103), 'chem. emission nonstat', 'stop' )


 END SUBROUTINE chem_emis_nonstationary_update


 END MODULE chem_emis_nonstationary_mod
