!> @file chem_emis_biogenic_mod.f90
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
! Copyright 2017-2021 Karlsruhe Institute of Technology
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Module for emissions from  biogenic sources i.e. trees and vegetation patchess.
!> Contains interface functions, subroutines and algorithms for calculations of biogenic
!> emissions.Public methods are specific to biogenic emission mode and has
!> prefix 'chem_emis_biogenic'.
!> Currently the bvoc mddule calculates emissions only from resolved (high vegetation with
!> leaf area density) vegetation. Unresolved trees and vegetation patches are treeted as
!> vegetation with 0 emissions.
!--------------------------------------------------------------------------------------------------!
 MODULE chem_emis_biogenic_mod

    USE chem_emis_biogenic_data_mod

    USE chem_emis_generic_mod

    USE chem_emis_vsrc_mod

    USE chem_modules,                                                                              &
        ONLY:  emis_biogenic_lod,                                                                  &
               nc_field_length

    USE control_parameters,                                                                        &
        ONLY:  message_string

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               topo_top_ind

    USE kinds

    USE radiation_model_mod,                                                                       &
        ONLY:  npcbl,                                                                              &
               pcbl,                                                                               &
               pcinsw,                                                                             &
               pcinswdif,                                                                          &
               pcinswdir

    IMPLICIT NONE

    PRIVATE

    SAVE

    CHARACTER(LEN = *), PARAMETER ::  dim_nvegetation_pars  = 'nvegetation_pars'  !< reading from netCDF static file
    CHARACTER(LEN = *), PARAMETER ::  dim_zlad              = 'zlad'              !< reading from netCDF static file
    CHARACTER(LEN = *), PARAMETER ::  str_novalue           = 'novalue'           !< string array initialization
    CHARACTER(LEN = *), PARAMETER ::  var_bad               = 'bad'               !< to read netCDF var fm static file
    CHARACTER(LEN = *), PARAMETER ::  var_lad               = 'lad'               !< to read netCDF var fm static file
    CHARACTER(LEN = *), PARAMETER ::  var_patch_height      = 'patch_height'      !< to read netCDF var fm static file
    CHARACTER(LEN = *), PARAMETER ::  var_soil_type         = 'soil_type'         !< to read netCDF var fm static file
    CHARACTER(LEN = *), PARAMETER ::  var_tree_type         = 'tree_type'         !< to read netCDF var fm static file
    CHARACTER(LEN = *), PARAMETER ::  var_vegetation_height = 'vegetation_height' !< to read netCDF var fm static file
    CHARACTER(LEN = *), PARAMETER ::  var_vegetation_pars   = 'vegetation_pars'   !< to read netCDF var fm static file
    CHARACTER(LEN = *), PARAMETER ::  var_vegetation_type   = 'vegetation_type'   !< to read netCDF var fm static file

    INTEGER(iwp) ::  lod                   !< lod, 0, 1, or 2
    INTEGER(iwp) ::  num_def_matched_spcs  !< number of default spcs from data file matched with mech spcs
    INTEGER(iwp) ::  num_emis_species      !< number of emission species namelist spcs matched with mech spcs
    INTEGER(iwp) ::  num_pft_default       !< number of default pft from data file
    INTEGER(iwp) ::  num_pft_namelist      !< number of pft from namelist
    INTEGER(iwp) ::  num_tree_default      !< number of default tree tyes from data file
    INTEGER(iwp) ::  num_tree_namelist     !< number of tree types from namelist
    INTEGER(iwp) ::  num_vsrc_pos          !< number of volume sources
    INTEGER(iwp) ::  nvegetation_pars      !< number of vegetation parameters,  reading from netCDF static file
    INTEGER(iwp) ::  zlad                  !< number of vertical layers for leaf area density (LAD) from netCDF static file
    INTEGER(iwp) ::  seconds_since_last_update  !<  timestamp control
    INTEGER(iwp) ::  seconds_to_next_update     !< timestamp control

    INTEGER(KIND=1), PARAMETER ::  fill_value_byte = -127  !< fill value for integer kind 1 variables

    INTEGER(iwp), PARAMETER ::  int_novalue = -127  !< fill value for integer variables

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  pft_default    !< list of default pfts from data file
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  pft_namelist   !< list of pft from namelist
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  tree_default   !< list of default tree types from data file
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  tree_namelist  !< list of tree types from namelist
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:) ::  vsrc_pcbl_map  !< pcbl mapping (radiation)

    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  soil_type    !< reading variable from netCDF static file
    INTEGER(iwp), ALLOCATABLE, DIMENSION(:,:) ::  vegetation_type    !< reading variable from netCDF static file

    INTEGER(kind = 1  ), ALLOCATABLE, DIMENSION(:,:,:) ::  tree_type !< reading variable from netCDF static file

    LOGICAL ::  dt_msg_flag = .TRUE.   !< flag for message to appear only once
    LOGICAL ::  is_emis_updated        !< timestamp control
    LOGICAL ::  sm_flag = .TRUE.       !< flat for message to appear only once
    LOGICAL ::  update_emissions       !< timestamp control, emission update
    LOGICAL ::  vsrc_has_pcbl_mapping  !< volume source radiation for the first time step

    REAL(dp), ALLOCATABLE, DIMENSION(:,:) ::  ef_def_pft       !< default emission potentials for pft from data file
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) ::  ef_pft           !< emission potentials for pft from namelist
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) ::  ef_def_tree      !< default emission potentials for tree types from data file
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) ::  ef_tree          !< emission potentials for trees from namelist
    REAL(dp), ALLOCATABLE, DIMENSION(:,:) ::  vsrc_emis_value  !< volume source value

    REAL(wp) ::  current_time       = 0.0_wp !< timestamp control, current time
    REAL(wp) ::  current_second     = 0.0_wp !< timestamp control, current seconds
    REAL(wp) ::  store_current_time = 0.0_wp !< timestamp control, store current time

    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::  patch_height       !< reading variable from netCDF static file
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::  vegetation_height  !< reading variable from netCDF static file
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  bad                !< reading variable from netCDF static file
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  lad                !< reading variable from netCDF static file
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  vegetation_pars    !< reading variable from netCDF static file

    REAL(KIND=wp), PARAMETER ::  fill_value_real = -9999.0_wp
!
!-- data type for saving species and their linkage to KPP chemical mechanism
    TYPE(chem_emis_species),  ALLOCATABLE, DIMENSION(:) ::  def_matched_spcs  !< default emission specs matched with mech spcs
    TYPE(chem_emis_species),  ALLOCATABLE, DIMENSION(:) ::  emis_species      !< namlist emission specs matched with mech spcs
    TYPE(chem_emis_species),  ALLOCATABLE, DIMENSION(:) ::  matched_spcs      !< to save position of def matched spcs in ltab_ebio_nam
    TYPE(chem_emis_vsrc_pos), ALLOCATABLE, DIMENSION(:) ::  vsrc_pos          !< volume source positions


    INTERFACE chem_emis_biogenic_init
       MODULE PROCEDURE chem_emis_biogenic_init
    END INTERFACE chem_emis_biogenic_init

    INTERFACE chem_emis_biogenic_update
       MODULE PROCEDURE chem_emis_biogenic_update
    END INTERFACE chem_emis_biogenic_update
!
!-- public methods
    PUBLIC ::  chem_emis_biogenic_init,                                                            &
               chem_emis_biogenic_update

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Intitializing biogenic emission mode
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_biogenic_init( )

    USE chem_modules,                                                                              &
        ONLY:  ebio_dt,                                                                            &
               emis_biogenic_lod

    IMPLICIT NONE
!
!-- CALL emis_biogenic_netcdf_data_surface
!-- Add code to figure out what LOD it is in
    CALL emis_rename_spcs_ltab2mech
!
!-- initialize emission mode based on LOD
    lod = emis_biogenic_lod
!
!-- select which lod(s) requires timestep control
    IF ( lod < 2 )  THEN
       seconds_to_next_update = ebio_dt
       update_emissions = .TRUE.
       is_emis_updated  = .FALSE.
    ELSE
       update_emissions = .FALSE.
       is_emis_updated  = .FALSE.
    ENDIF
!
!-- initialization of the respective lod
    SELECT CASE( lod )

       CASE( 0 )
          CALL emis_biogenic_init_lod0( )

       CASE( 1 )
          CALL emis_biogenic_init_lod1( )

       CASE( 2 )
          CALL emis_biogenic_init_lod2( )

       CASE DEFAULT
          WRITE( message_string, * ) 'illegal lod = ', lod
          CALL message( 'chem_emis_biogenic_init', 'CHM0022', 1, 2, 0, 6, 0 )

    END SELECT

 END SUBROUTINE chem_emis_biogenic_init

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates biogenic emissions at every tiem step (called in time_integratin)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE chem_emis_biogenic_update( )


    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_get_current_time

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_update_source

    USE chem_modules,                                                                              &
        ONLY:  ebio_dt

    USE control_parameters,                                                                        &
        ONLY:  time_since_reference_point

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    IMPLICIT NONE

    INTEGER(iwp) ::  c_pos  !< counter for gird cell position with vegetation
    INTEGER(iwp) ::  i_spc  !< counter for bvoc
    INTEGER(iwp) ::  kp     !< counter for 1D plant canopy grid cell array
    INTEGER(iwp) ::  ktopo  !< topological offset in vertical direction (for pcbl)

    REAL(wp) ::  dt_timestep, dt_biogenic

    dt_biogenic = ebio_dt
!
!--skip if there are no volume sources
    IF ( num_vsrc_pos == 0 )  RETURN
!
!-- initialize volume source radiation for the first time step only
    IF ( .NOT. vsrc_has_pcbl_mapping )  THEN
       DO  c_pos = 1, num_vsrc_pos
          vsrc_pcbl_map(c_pos) = -1
          ktopo = topo_top_ind(vsrc_pos(c_pos)%j,vsrc_pos(c_pos)%i,0)
          DO  kp = 1, npcbl
             IF ( (vsrc_pos(c_pos)%i == pcbl(4,kp) )  .AND.                                        &
                ( vsrc_pos(c_pos)%j  == pcbl(3,kp) )  .AND.                                        &
                ( vsrc_pos(c_pos)%k  == pcbl(2,kp)+ktopo) )  THEN
                vsrc_pcbl_map(c_pos) = kp
                EXIT
             ENDIF
          ENDDO
       ENDDO
       vsrc_has_pcbl_mapping = .TRUE.
    ENDIF
!
!-- determine whether it is time to update (lod 0/1 only)
    IF ( lod < 2 )  THEN
       CALL get_date_time( time_since_reference_point, second=current_second )
       current_time = time_since_reference_point
       dt_timestep = current_time - store_current_time
!
!--    if ebio_dt not defined then use model time-step
       IF ( dt_biogenic == 0.0 ) dt_biogenic = dt_timestep

       IF ( dt_biogenic > 0.0  .AND.  dt_biogenic < dt_timestep)  THEN
          dt_biogenic = dt_timestep
          IF ( dt_msg_flag )  THEN
             message_string = 'The biogenic emission module time-step is smaller than' //          &
                              'the model time-step.&Reverting to the models time-step.'
             CALL message( 'chem_emis_biogenic_update', 'CHM0023', 0, 1, 0, 6, 0 )
             dt_msg_flag = .FALSE.
          ENDIF
       ENDIF

       store_current_time = current_time
       seconds_since_last_update = MOD( current_time, dt_biogenic )
       seconds_to_next_update = dt_biogenic - seconds_since_last_update

       IF ( (seconds_since_last_update == 0 )  .OR.  ( seconds_to_next_update < 0 ) )  THEN
          IF ( .NOT. is_emis_updated ) update_emissions = .TRUE.
       ELSE
          IF ( is_emis_updated ) is_emis_updated = .FALSE.
       ENDIF
    ENDIF
!
!-- update bvoc emission of he respecitve lod
    SELECT CASE( lod )

       CASE( 0 )
          IF ( update_emissions )  THEN
             CALL emis_biogenic_update_lod0( )
             update_emissions = .FALSE.
             is_emis_updated = .TRUE.
          ENDIF

       CASE( 1 )
          CALL emis_biogenic_update_lod1( )

       CASE( 2 )
          CALL emis_biogenic_update_lod2( )

       CASE DEFAULT
          WRITE( message_string, * ) 'illegal lod = ', lod
          CALL message( 'chem_emis_biogenic_init', 'CHM0022', 1, 2, 0, 6, 0 )

    END SELECT
!
!-- updates volume source values for a given specs
    DO  i_spc = 1, num_emis_species
       CALL chem_emis_vsrc_update_source( emis_species(i_spc)%mech_index,                          &
                                          vsrc_pos, vsrc_emis_value(i_spc,:) )
    ENDDO

 END SUBROUTINE chem_emis_biogenic_update

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> initialize biogenic emission mode module under LOD 0
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_biogenic_init_lod0( )

    IMPLICIT NONE
!
!-- read p3d file, and default data from data file, read static data from _PIDS file
!-- finally get 3d bvoc source locations for namelist pft and trees
    CALL location_message( 'initializing biogenic emissions module for LOD 0', 'start' )
    CALL emis_get_namelist_data( )
    CALL emis_map_default_data( )
    CALL emis_nc_get_vegetation_data( )
    CALL emis_bio_ltab_index( )
    CALL emis_get_vsrc_positions( )
    CALL location_message( 'initializing biogenic emissions module for LOD 0', 'finished' )

 END SUBROUTINE emis_biogenic_init_lod0

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> initialize biogenic emission mode module under LOD 1
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_biogenic_init_lod1( )

    message_string = 'biogenic emission module MEGAN 2.1 has not been implemented yet'
    CALL message( 'emis_biogenic_init_lod1', 'CHM0024', 1, 2, 0, 6, 0 )

 END SUBROUTINE emis_biogenic_init_lod1

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> initialize biogenic emission mode module under LOD 2
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_biogenic_init_lod2( )

    message_string = 'biogenic emission module LOD 2 has not been implemented yet'
    CALL message( 'emis_biogenic_init_lod1', 'CHM0024', 1, 2, 0, 6, 0 )

 END SUBROUTINE emis_biogenic_init_lod2

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates biogenic emission volume source under LOD 0
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_biogenic_update_lod0( )

    USE chem_emis_vsrc_mod,                                                                        &
        ONLY:  chem_emis_vsrc_update_source

    IMPLICIT NONE
!
!-- local variables,
    INTEGER(iwp) ::  c_pos  !< counter for grid cell position
    INTEGER(iwp) ::  i_spc  !< counter for matched bio spcs.

    IF ( num_vsrc_pos == 0 )  RETURN
!
!-- calculate bvoc emissions for all volume source positions
    DO  i_spc = 1, num_def_matched_spcs
       DO  c_pos = 1, num_vsrc_pos
          vsrc_emis_value(i_spc,c_pos) = calc_emis_biogenic ( i_spc, c_pos )
       ENDDO
    ENDDO

 END SUBROUTINE emis_biogenic_update_lod0

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates biogenic emission volume source under LOD 1
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_biogenic_update_lod1( )
!
!-- right now it does nothing

 END SUBROUTINE emis_biogenic_update_lod1

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates biogenic emission volume source under LOD 2
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_biogenic_update_lod2( )
!
!-- right now it does nothing

 END SUBROUTINE emis_biogenic_update_lod2

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> reading data from from the namelist
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_get_namelist_data( )

    USE chem_modules,                                                                              &
        ONLY:  ebio_ef_pft,                                                                        &
               ebio_ef_tree,                                                                       &
               ebio_emis_name,                                                                     &
               ebio_pft,                                                                           &
               ebio_tree,                                                                          &
               nc_field_length



    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_species_match
!
!-- USE chem_gasphase_mod, ONLY : spc_names
    INTEGER(iwp) ::  num_namelist_species  !< total number of namelist spcs.
    INTEGER(iwp) ::  i_spc                 !< index for namelist spcs.
    INTEGER(iwp) ::  i_tree                !< index for tree spcs and pfts
    INTEGER(iwp) ::  j_spc                 !< index for emis spcs

    CHARACTER(LEN=nc_field_length), ALLOCATABLE, DIMENSION(:) ::  namelist_species   !< list of spcs namelist

!
!-- extract species from namelist
    i_spc = 1
    DO WHILE ( ebio_emis_name(i_spc) /= str_novalue )
       i_spc = i_spc + 1
    ENDDO

    num_namelist_species = i_spc - 1

    IF ( num_namelist_species == 0 )  THEN
       message_string = 'no biogenic species defined'
       CALL message( 'emis_get_namelist_data', 'CHM0025', 1, 2, 0, 6, 0 )
    ENDIF

    ALLOCATE( namelist_species(num_namelist_species) )
!
!-- read bvoc's from the namelist file
    DO  i_spc = 1, num_namelist_species
       namelist_species(i_spc) = ebio_emis_name(i_spc)
    ENDDO
!
!-- from namelist species obtain only those appearing in KPP mechanism
    CALL chem_emis_species_match( num_emis_species, emis_species, namelist_species )
!
!-- determine the number of plantf and tree types from the name list
!-- and allocate space for profiles.
    i_tree = 1
    DO WHILE ( ebio_pft(i_tree) /= int_novalue )
       i_tree = i_tree + 1
    ENDDO
    num_pft_namelist = i_tree - 1

    i_tree = 1
    DO WHILE ( ebio_tree(i_tree) /= int_novalue )
       i_tree = i_tree + 1
    ENDDO
    num_tree_namelist = i_tree - 1

    ALLOCATE( pft_namelist(num_pft_namelist) )
    ALLOCATE( tree_namelist(num_tree_namelist) )
!
!-- read pft and single tree infor from nameist,
    pft_namelist(1:num_pft_namelist)   = ebio_pft(1:num_pft_namelist)
    tree_namelist(1:num_tree_namelist) = ebio_tree(1:num_tree_namelist)

    ALLOCATE( ef_pft (num_emis_species,num_pft_namelist) )
    ALLOCATE( ef_tree(num_emis_species,num_tree_namelist  ) )
!
!-- reading EF for pft and trees from namelist.
    j_spc = 1
    DO  i_spc = 1, num_namelist_species
       IF ( TRIM( namelist_species(i_spc) ) == TRIM( emis_species(j_spc)%name_str ) )  THEN
          ef_pft(j_spc,1:num_pft_namelist)   = ebio_ef_pft (i_spc,1:num_pft_namelist)
          ef_tree(j_spc,1:num_tree_namelist) = ebio_ef_tree(i_spc,1:num_tree_namelist )
          j_spc = j_spc + 1
       ENDIF
    ENDDO

    DEALLOCATE( namelist_species )

 END SUBROUTINE emis_get_namelist_data

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine exctracts matched species between chem_mech and ltab (biogeic data file) and then
!> allocate default ef from ltab  to all default pfts given in the  data file )
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_map_default_data( )

    USE chem_modules,                                                                              &
        ONLY:  nc_field_length

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_species_match
!
!-- USE chem_gasphase_mod, ONLY : spc_names
!-- local variables
    CHARACTER(LEN=nc_field_length), ALLOCATABLE, DIMENSION(:) ::  def_spcs_nam

    INTEGER(iwp) ::  i_spc         !< running index for default bvocs
    INTEGER(iwp) ::  j_spc         !< runnig index for the default matched  bvocs
    INTEGER(iwp) ::  num_def_spcs  !< total number of default bvocs


    num_def_spcs = SIZE( ltab_ebio_nam )

    ALLOCATE( def_spcs_nam(num_def_spcs) )

!-- read  bvocs from data file
    DO  i_spc = 1, num_def_spcs
       def_spcs_nam(i_spc) = ltab_ebio_nam(i_spc)
    ENDDO
!
!-- From default species from lookup table obtain only those appearing in KPP mechanism
    CALL chem_emis_species_match( num_def_matched_spcs, def_matched_spcs, def_spcs_nam )
!
!-- determine number of PFTs and tree types from the name list and allocate space for profiles.
!-- get number of default pfts from data file
    num_pft_default = SIZE( ltab_pft )
!
!-- get number of default tree types from data file
    num_tree_default = n_tree_types

    ALLOCATE( pft_default(num_pft_default  )  )
    ALLOCATE( tree_default(1:num_tree_default)  )

!-- read pft and tree types data from look-up tables
    pft_default(1:num_pft_default)   = ltab_pft(1:num_pft_default)
    tree_default(1:num_tree_default) = ltab_tree(1:num_tree_default)

    ALLOCATE( ef_def_pft (num_def_matched_spcs,num_pft_default) )
    ALLOCATE( ef_def_tree(num_def_matched_spcs,num_tree_default  ) )
!
!-- if EF for inidividual trees not available the all single trees are allocated to their
!-- respective pftsfor EF.
    j_spc = 1
    DO  i_spc = 1, num_def_spcs
       IF ( j_spc <= num_def_matched_spcs )  THEN
          IF ( TRIM( def_spcs_nam(i_spc) ) == TRIM( def_matched_spcs(j_spc)%name_str) )  THEN
             ef_def_pft (j_spc, 1:num_pft_default) = ltab_ef_pft(i_spc, 1:num_pft_default)
             ef_def_tree(j_spc, 1:num_pft_default) = ltab_ef_pft(i_spc, 1:num_pft_default)
             j_spc = j_spc + 1
          ENDIF
       ENDIF
    ENDDO

    DEALLOCATE( def_spcs_nam )

 END SUBROUTINE emis_map_default_data

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reads BVOC vegetation Input data from static file.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_nc_get_vegetation_data ( )

    USE control_parameters,                                                                        &
       ONLY:  coupling_char

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  close_input_file,                                                                   &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               input_file_static,                                                                  &
               open_read_file

    IMPLICIT NONE

    INTEGER(iwp) :: ncid  !< netcdf handle for static file

    CALL open_read_file( TRIM(input_file_static) //  TRIM( coupling_char ), ncid )
    CALL get_dimension_length( ncid, nvegetation_pars, dim_nvegetation_pars )
    CALL get_dimension_length( ncid, zlad, dim_zlad )

    ALLOCATE( tree_type(0:zlad-1,nys:nyn,nxl:nxr) )
    ALLOCATE( vegetation_type(nys:nyn,nxl:nxr) )
    ALLOCATE( patch_height(nys:nyn,nxl:nxr) )
    ALLOCATE( vegetation_height(nys:nyn,nxl:nxr) )
    ALLOCATE( vegetation_pars(0:nvegetation_pars-1,nys:nyn,nxl:nxr) )
    ALLOCATE( bad(0:zlad-1,nys:nyn,nxl:nxr) )
    ALLOCATE( lad(0:zlad-1,nys:nyn,nxl:nxr) )
    ALLOCATE( soil_type(nys:nyn,nxl:nxr) )

    CALL get_variable( ncid, var_soil_type, soil_type, nxl, nxr, nys, nyn )
    CALL get_variable( ncid, var_tree_type, tree_type, nxl, nxr, nys, nyn, 0, zlad - 1 )
    CALL get_variable( ncid, var_vegetation_type, vegetation_type, nxl, nxr, nys, nyn )
    CALL get_variable( ncid, var_patch_height, patch_height, nxl, nxr, nys, nyn )
    CALL get_variable( ncid, var_vegetation_height, vegetation_height, nxl, nxr, nys, nyn )
    CALL get_variable( ncid, var_vegetation_pars, vegetation_pars, nxl, nxr, nys, nyn,             &
                                                                        0, nvegetation_pars-1 )
    CALL get_variable( ncid, var_lad, lad, nxl, nxr, nys, nyn, 0, zlad-1 )
    CALL get_variable( ncid, var_bad, bad, nxl, nxr, nys, nyn, 0, zlad-1 )
    CALL close_input_file( ncid )

 END SUBROUTINE emis_nc_get_vegetation_data

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> update volume source positions
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_get_vsrc_positions( )

    USE chem_emis_generic_mod,                                                                     &
        ONLY:  chem_emis_init_volume_source

    IMPLICIT NONE

    INTEGER(iwp) ::  i      !< grid running index in x-direction
    INTEGER(iwp) ::  j      !< grid running index in y-direction
    INTEGER(iwp) ::  k      !< grid running index in z direction
    INTEGER(iwp) ::  kpos   !< counter for vsrc positions.
    INTEGER(iwp) ::  zlad1  !< vertical levels in plant canopy
!
!-- first pass to determine number of volume sources and
!-- allocate memory accordingly
    zlad1 = zlad - 1
    num_vsrc_pos = 0

    DO  j = nys, nyn
       DO  i = nxl, nxr
          DO  k = 0, zlad1
             IF ( tree_type(k,j,i) == fill_value_byte ) CYCLE
             num_vsrc_pos = num_vsrc_pos + 1
          ENDDO
       ENDDO
    ENDDO
!
!-- skip if there are no volume sources
    IF ( num_vsrc_pos == 0 )  RETURN

    ALLOCATE( vsrc_pos(num_vsrc_pos) )
!
!-- second pass to populate volume source positions
    kpos = 0
    DO  j = nys, nyn
       DO  i = nxl, nxr
          DO  k = 0, zlad1
             IF ( tree_type(k,j,i) == fill_value_byte ) CYCLE
                kpos = kpos + 1
                vsrc_pos(kpos)%i = i
                vsrc_pos(kpos)%j = j
                vsrc_pos(kpos)%k = k
          ENDDO
       ENDDO
    ENDDO
!
!-- append volume sources to global vsrc
    CALL chem_emis_init_volume_source( vsrc_pos )
!
!-- allocate space for local volume sources
    ALLOCATE( vsrc_emis_value(num_def_matched_spcs,num_vsrc_pos) )
!
!-- mapping to PCBL
    ALLOCATE( vsrc_pcbl_map(num_vsrc_pos) )
!
!-- vsrc_pcbl_map can only be initialized after initialization of radiation module
    vsrc_has_pcbl_mapping = .FALSE.

 END SUBROUTINE emis_get_vsrc_positions

!--------------------------------------------------------------------------------------------------!
!Description:
!------------
!>  Map emission species to biogenic classes
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_bio_ltab_index( )

    IMPLICIT NONE

    INTEGER(iwp) ::  k_spc  !< running index matched spcs
    INTEGER(iwp) ::  l_spc  !< running index spcs from lookup table

    INTEGER, DIMENSION(1) ::  bvoc_index !< running index for bvoc spcs

    ALLOCATE( matched_spcs(num_def_matched_spcs) )
!
!-- find matched spcs indecies in the default spcs list
    DO  k_spc = 1, num_def_matched_spcs
       matched_spcs(k_spc)%name_str = TRIM( def_matched_spcs(k_spc)%name_str )
!
!--    assign spcs index
       bvoc_index = MINVAL( PACK( [(l_spc,l_spc=1,size(ltab_ebio_nam) )],                          &
                            ltab_ebio_nam == TRIM( matched_spcs(k_spc)%name_str ) ) )
       matched_spcs(k_spc)%mech_index = bvoc_index(1)
    ENDDO

 END SUBROUTINE emis_bio_ltab_index

!--------------------------------------------------------------------------------------------------!
!Description:
!------------
!>  Make respective bvocs in the data file consistent with the bvocs in the chem mechanism.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE emis_rename_spcs_ltab2mech( )

    USE chem_gasphase_mod,                                                                         &
        ONLY:  cs_mech

    SELECT CASE( cs_mech )

       CASE( 'cbm4' )
          ltab_ebio_nam (isop)  = 'ISOP'
          ltab_ebio_nam (xvoc3) = 'HCHO'
          ltab_ebio_nam (ovoc3) = 'ETH'
          ltab_ebio_nam (ovoc5) = 'TOL'

       CASE( 'smog' )
!
!         isop + other vocs
          ltab_ebio_nam (isop)  = 'RH'
!
!--       hcho + other higher aldihydes
          ltab_ebio_nam (xvoc3) = 'RCHO'

       CASE( 'radm2' )
!
!--       right now, do nothing

       CASE( 'mozcart' )
!
!--       right now, do nothing

       CASE( 'cbmz' )
!
!--       right now, do nothing

       CASE DEFAULT
          message_string = 'unknow cs_mech = "' // TRIM( cs_mech ) // '"'
          CALL message( 'chem_emis_biogenic_init', 'CHM0026', 1, 2, 0, 6, 0 )

    END SELECT

 END SUBROUTINE emis_rename_spcs_ltab2mech

!--------------------------------------------------------------------------------------------------!
!Description:
!------------
!>  provides baseline emission factor and apply correction factor to calculate biogenic emissions
!--------------------------------------------------------------------------------------------------!
 FUNCTION calc_emis_biogenic( kspecies, kvsrc ) RESULT( emis_biogenic )

    USE arrays_3d,                                                                                 &
        ONLY:  exner,                                                                              &
               pt

    USE chem_modules,                                                                              &
        ONLY:  ebio_rad_method

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  kspecies        !< species index
    INTEGER(iwp), INTENT(IN) ::  kvsrc           !< volume source index

    INTEGER(iwp) ::  class_pft                   !< pft class for individual trees
    INTEGER(iwp) ::  i                           !< grid running index in  x-direction
    INTEGER(iwp) ::  j                           !< grid running index in y-direction
    INTEGER(iwp) ::  k                           !< grid running index in z direction
    INTEGER(iwp) ::  kp                          !< counter for 1D plant canopy grid cell array
    INTEGER(iwp) ::  this_pft                    !< pft of the grid cell for sm algorithm
    INTEGER(iwp) ::  this_soil_type              !< soil type of current grid cell
    INTEGER(iwp) ::  this_tree_type              !< this_plant_type

    INTEGER(iwp), DIMENSION(1) ::  i_def_pft     !< counter for default pft type
    INTEGER(iwp), DIMENSION(1) ::  i_def_tree    !< counter for default tree type
    INTEGER(iwp), DIMENSION(1) ::  i_pft         !< counter for pft type
    INTEGER(iwp), DIMENSION(1) ::  i_tree        !< counter for tree type

    REAL(KIND=dp) ::  base_emis                  !< emission potentials from namelist or default
    REAL(KIND=dp) ::  correction_factor          !< correction of emissions from various environmental factors
    REAL(KIND=dp) ::  emis_biogenic              !< output emission factor

    REAL(wp) ::  lad_g_per_m3                    !< conversion factor from m2/m3 to g/m3
    REAL(wp) ::  radi_sw                         !< total (dir + diffu) radiation at the current grid cell
    REAL(wp) ::  radi_swdif                      !< diffused radiatoin at the current grid cell
    REAL(wp) ::  radi_swdir                      !< direct radiation at the current grid cell
    REAL(wp) ::  temp_leaf                       !< leaf temperature (sudo) at the current grid cell
    REAL(wp) ::  this_tree_lad                   !< lad at the current grid cell

    radi_sw    = fill_value_real
    radi_swdif = fill_value_real
    radi_swdir = fill_value_real

    kp = vsrc_pcbl_map(kvsrc)
!
!-- calcualte the user defined  radiation method for bvoc emissions
    SELECT CASE( ebio_rad_method )

       CASE( 0 )
          IF ( kp > 0 )  THEN
            radi_sw    = pcinsw (kp)
            radi_swdif = pcinswdif(kp)
            radi_swdir = pcinswdir(kp)
          ENDIF

       CASE( 1 )
          message_string = 'absorption radiation method has not been implemented'
         CALL message( 'gamma_sm', 'CHM0027', 1, 2, 0, 6, 0 )

       CASE( 2 )
          message_string = 'fractional radiation method has not been implemented'
         CALL message( 'gamma_sm', 'CHM0027', 1, 2, 0, 6, 0 )

       CASE DEFAULT
          WRITE( message_string, * ) 'illegal radiation method ebio_rad_method = ', ebio_rad_method
          CALL message( 'gamma_sm', 'CHM0028', 1, 2, 0, 6, 0 )

    END SELECT

    k = vsrc_pos(kvsrc)%k
    j = vsrc_pos(kvsrc)%j
    i = vsrc_pos(kvsrc)%i

    this_tree_type = tree_type(k,j,i)
    this_tree_lad  = lad(k,j,i)

    lad_g_per_m3 = lad_unitconv(this_tree_lad, this_tree_type, 1.0_wp)
!
!-- match the spcs name  an find position of the tree/pft type in the namelist
    i_tree = MINLOC( tree_namelist, MASK = tree_namelist == this_tree_type )
    i_pft = MINLOC( pft_namelist, MASK = pft_namelist == this_tree_type )

    IF ( (i_tree(1) == 1)  .AND.  (tree_namelist(1) /= this_tree_type) ) i_tree(1) = 0
    IF ( (i_pft(1) == 1)  .AND.  (pft_namelist(1) /= this_tree_type) ) i_pft(1) = 0
!
!-- match the spcs name  an find position of the tree/pft type in data file
    i_def_tree = MINLOC( tree_default, MASK = tree_default == this_tree_type )
    i_def_pft = MINLOC( pft_default, MASK = pft_default == this_tree_type )

    IF ( (i_def_tree(1) == 1)  .AND.  (tree_default(1) /= this_tree_type) ) i_def_tree(1) = 0
    IF ( (i_def_pft(1) == 1)  .AND.  (pft_default(1) /= this_tree_type) ) i_def_pft(1) = 0
!
!- -default EF from lookup table allocated to the current grid cell
    base_emis = ltab_def_plant(10)
    IF ( (this_tree_type >= 0)  .AND.  (i_tree(1) > 0)  .AND.  (kspecies <= num_emis_species) )  THEN
       base_emis = ef_tree(kspecies,i_tree(1))
       IF ( k == 0 ) this_pft = ltab_tree_map(i_def_tree(1) )
    ELSEIF ( (this_tree_type <= 0)  .AND.  (i_pft(1) > 0)  .AND.  (kspecies <= num_emis_species) )  THEN
       base_emis = ef_pft(kspecies,i_pft(1))
       IF ( k == 0 ) this_pft = i_pft(1)
    ELSEIF ( (this_tree_type > 0)  .AND.  (i_tree(1) == 0)  .AND.  (i_def_tree(1) > 0) )  THEN
       class_pft = ltab_tree_map(i_def_tree(1) )
       base_emis = ef_def_tree(kspecies, class_pft) * conv_grams2moles(kspecies)
       IF ( k == 0 ) this_pft = ltab_tree_map(i_def_tree(1) )
    ELSEIF ( (this_tree_type < 0)  .AND.  (i_pft(1) == 0)  .AND.  (i_def_pft(1)  >  0) )  THEN
       base_emis = ef_def_pft(kspecies, i_def_pft(1) ) * conv_grams2moles(kspecies)
       IF ( k == 0 ) this_pft = i_def_pft(1)
    ENDIF
!
!-- get air temperature of the grid cell as proxy of leaf temperature.
    temp_leaf = pt(k,j,i) * exner(k)
!
!-- get soil type from the static file
    this_soil_type  = soil_type(j,i)
!
!-- calculate correction factor of the biogenic emissions
    correction_factor = calc_gamma( kspecies, radi_sw, temp_leaf, this_soil_type, this_pft )
!
!-- calculate biogenic emission from tree of in the grid cell
    emis_biogenic = (base_emis * correction_factor * this_tree_lad)

 END FUNCTION calc_emis_biogenic

!--------------------------------------------------------------------------------------------------!
!Description:
!------------
!>  Calcualte correction factor (gamma) to apply on base_emissions
!--------------------------------------------------------------------------------------------------!
 FUNCTION calc_gamma( k_spc, radi_sw1, temp_leaf1, soil_type_id, pft_msoil ) RESULT( gamma )

    USE chem_modules,                                                                              &
        ONLY:  ebio_soilm_method

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  k_spc         !< species index
    INTEGER(iwp), INTENT(IN) ::  pft_msoil     !< pft for soil moisture calculation
    INTEGER(iwp), INTENT(IN) ::  soil_type_id  !< soil type index

    REAL(wp), INTENT(IN) ::  radi_sw1    !< radiation at 'the grid cell
    REAL(wp), INTENT(IN) ::  temp_leaf1  !< should come from 3d temp field in the plant canopy

    REAL(wp) ::  ldf           !< lif dependence factor
    REAL(wp) ::  gam_p_r       !< return result of light dependence algorithm
    REAL(wp) ::  gam_sm_r      !< return result of soil moiture algorithm
    REAL(wp) ::  gam_sn_r      !< return result of seasonal effect algorithm
    REAL(wp) ::  gam_tisop_r   !< return result of temperature  effect algorithm
    REAL(wp) ::  gam_tnisop_r  !< return result of light independent algorithm

    REAL(KIND=dp) ::  gamma  !< output emission correction factor

    gam_p_r      = gamma_p( radi_sw1 )
    gam_tisop_r  = gamma_tisop( temp_leaf1 )
    gam_tnisop_r = gamma_tnisop( k_spc, temp_leaf1 )
    gam_sn_r     = gamma_sn( )

!
!-- calculation of soil moistrue effect based on the user selection
!-- soil moisture factor applied only on isoprene emissions
    IF ( TRIM( ebio_soilm_method ) == 'bulk' )  THEN
       IF ( matched_spcs(k_spc)%mech_index == isop )  THEN
          gam_sm_r = gamma_sm( soil_type_id )
       ELSE
          gam_sm_r = 1.0_wp
       ENDIF
    ELSEIF ( TRIM( ebio_soilm_method ) == 'weighted' )  THEN
       IF ( matched_spcs(k_spc)%mech_index == isop )  THEN
          gam_sm_r = gamma_dstress( soil_type_id, pft_msoil )
       ELSE
             gam_sm_r = 1.0_wp
       ENDIF
    ELSE
       message_string = 'unknown draught stress method ebio_soilm_method = "' //                   &
                        TRIM( ebio_soilm_method ) // '"'
       CALL message( 'gamma_sm', 'CHM0029', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- mapping to predetermined biogenic classes
    ldf = ltab_ldf_spc(matched_spcs(k_spc)%mech_index)
!
!-- calculate correction factor
    gamma = ( ( ( 1.0 - ldf ) * gam_tnisop_r )  +                                                  &
              ( ldf * gam_p_r * gam_tisop_r * gam_sm_r * gam_sn_r ) )

 END FUNCTION calc_gamma

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Light dependence algorithm (gamma_p) for the biogenic species.
!--------------------------------------------------------------------------------------------------!
 FUNCTION gamma_p( radi_sw2 ) RESULT( gam_p )

    USE chem_modules,                                                                              &
        ONLY:  ebio_ppfd_factor

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: radi_sw2  !< radiaton on the particular grid cell

    REAL(wp) :: gam_p  !< result of the gamma_p calculation
    REAL(wp) :: ppfd   !< default ppfd for standard conds under control envi.
    REAL(wp), PARAMETER   :: alpha = 0.0027_wp  !< imperical constant
    REAL(wp), PARAMETER   :: cl_1  = 1.066_wp   !< units J mol-1

    ppfd  = 1000.0_wp
!
!-- calc ppfd from 3-d direct+ diffused radiation in and above the plant canopy
    IF ( .NOT. radi_sw2 == fill_value_real ) ppfd = radi_sw2 * ebio_ppfd_factor
!
!-- calculate light dependence of bvoc emissions
    gam_p = (alpha * cl_1 * ppfd/ SQRT(1.0 + ( alpha * alpha * ppfd * ppfd ) ) )

 END FUNCTION gamma_p

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Temperatrure dependence algorithm (gamma_tisop) for the biogenic species.
!--------------------------------------------------------------------------------------------------!
 FUNCTION gamma_tisop( temp_leaf2 ) RESULT( gam_tisop )

    IMPLICIT NONE

    REAL(wp) ::  gam_tisop        !< return value of gamma_tisop calculation
    REAL(wp) ::  denom_gam_tisop  !< store denominator of the equation
    REAL(wp) ::  numer_gam_tisop  !< store numerator of the equation
    REAL(wp) ::  temp_leaf2       !< units K, to get from 3D radiation field

    REAL(wp), PARAMETER ::  ct_1     = 95000.0_wp   !< imperical constant units J mol-1
    REAL(wp), PARAMETER ::  ct_2     = 230000.0_wp  !< imperical constant units J mol-1
    REAL(wp), PARAMETER ::  ct_3     = 0.961_wp     !< imperical constant
    REAL(wp), PARAMETER ::  r_gas    = 8.314_wp     !< gas constant, units J K-1 mol-1
    REAL(wp), PARAMETER ::  temp_m   = 314.0_wp     !< imperica constant, units K
    REAL(wp), PARAMETER ::  temp_std = 303.0_wp     !< temperature in standard condtions units K

    numer_gam_tisop = EXP( ct_1 * (temp_leaf2 - temp_std) / (r_gas * temp_std * temp_leaf2) )
    denom_gam_tisop = ct_3 + EXP( ct_2 * (temp_leaf2 - temp_m) / (r_gas * temp_std * temp_leaf2) )
!
!-- calcualte temperature dependence of the bvoc emissions
    gam_tisop = numer_gam_tisop / denom_gam_tisop

 END FUNCTION gamma_tisop

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Light independent algorithm (gamma_tnisop) for the biogenic species.
!>  Algorithm based on Guenther et al 1993
!--------------------------------------------------------------------------------------------------!
 FUNCTION gamma_tnisop( k_spc, temp_leaf3 ) RESULT( gam_tnisop )

    IMPLICIT NONE

    INTEGER(iwp) ::  k_spc   !< bvoc spcs

    REAL(wp) ::  temp_leaf3  !< ambient temp in the respective vertical level
    REAL(wp) ::  gam_tnisop  !< light independent emiss_activ_factor
    REAL(wp) ::  beta        !< imperical constant for k_spc for gamma_tnisop

    REAL(wp), PARAMETER ::  temp_std = 303.0_wp  !< standard temp

    beta = ltab_beta_spc(matched_spcs(k_spc)%mech_index)
!
!-- calculate light independent monoterpene emission fator
    gam_tnisop = EXP( beta * (temp_leaf3 - temp_std) )

 END FUNCTION gamma_tnisop

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Soil moisture (gamma_sm)dependence algorithm for the biogenic species.
!> This algorithm is based on Guenther et al. (2006)
!--------------------------------------------------------------------------------------------------!
 FUNCTION gamma_sm( soil_type_id ) RESULT( gam_sm )

    USE surface_mod,                                                                               &
        ONLY:  soil_moisture

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  soil_type_id  !< soil type of the grid cell.

    INTEGER(iwp) ::  i_layer  !< running index for soil layer

    REAL(wp) ::  gam_sm  !< return value of gamma_sm
    REAL(wp) ::  q1      !< imperical parameter
    REAL(wp) ::  qw      !< wilting default value
    REAL(wp) ::  soil_m  !< soil moiture unit m3 m-3

    REAL(wp), PARAMETER ::  del_q1 = 0.04_wp  !< imperical parameter

    i_layer = 0
    IF ( soil_moisture(i_layer) /= fill_value_real )  THEN
       soil_m = soil_moisture(i_layer)
    ELSE
!
!-- default mean value in the top surface layer
       soil_m = 0.3_wp
    ENDIF

    qw = ltab_qw(soil_type_id)
    q1 = qw + del_q1
!
!-- caclulate soil_moisture factor only for isoprene, for other spcs gam_sm = 1.
    IF ( soil_m > q1 )  THEN
       gam_sm = 1.0_wp
       IF ( sm_flag )  THEN
          WRITE(message_string, *) 'Soil is saturated. Draught stress factor is set to 1.0.'
          CALL message( 'gamma_sm', 'CHM0030', 0, 0, 0, 6, 0 )
          sm_flag = .FALSE.
       ENDIF
    ELSEIF ( soil_m > qw  .AND.  soil_m < q1 )  THEN
       gam_sm  = ( soil_m - qw ) / del_q1
    ELSEIF ( soil_m < qw )  THEN
       gam_sm = 0.0_wp
       message_string = 'Soil moisture is below wilting point. Draught stress factor is set to 0.0.'
       CALL message( 'gamma_sm', 'CHM0031', 0, 1, 0, 6, 0 )
    ENDIF

 END FUNCTION gamma_sm

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Soil moisture (gamma_sm)dependence algorithm for the biogenic species.
!> This algorithm is based on Rüdiger,Grote approach.
!--------------------------------------------------------------------------------------------------!
 FUNCTION gamma_dstress( soil_type_id, pft_msoil ) RESULT( gam_dstress )

    USE surface_mod,                                                                               &
        ONLY: soil_moisture

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  soil_type_id  !< soil type from ecmwf-ifs
    INTEGER(iwp), INTENT(IN) ::  pft_msoil     !< PFT for soil moisture calculation

    INTEGER(iwp) ::  i_layer            !< counter to calc soil layers
    INTEGER(iwp) ::  isl                !< running index for soil moisture layers.
    INTEGER(iwp) ::  iscale0            !< to scale root fraction layers to soil moisture layers.
    INTEGER(iwp) ::  iscale1            !< to scale root fraction layers to soil moisture layers.
    INTEGER(iwp) ::  n_msoil_layers     !< total number of soil moisture layers

    REAL(wp) ::  all_soil_layer_depth   !<
    REAL(wp) ::  cumulative_soil_depth  !< cumulative depth of soil moisture layers
    REAL(wp) ::  dstress_this_layer     !< draught stress on 'this' soil layer
    REAL(wp) ::  gam_dstress            !< return value of gamma_dstress
    REAL(wp) ::  qc                     !< field(water holding) capacity
    REAL(wp) ::  qw                     !< wilting default value

    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  relwater_content  !< relative water content
    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  root_dist         !< root distribution (unit m)
    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  root_frac         !< root fraction
    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  soil_frac         !< soil fraction
    REAL(wp), ALLOCATABLE, DIMENSION(:) ::  soil_m            !< soil moiture unit m3 m-3

    i_layer = 0
    DO WHILE ( soil_moisture(i_layer) /= fill_value_real )
       i_layer = i_layer + 1
    ENDDO
    n_msoil_layers = i_layer
!
!-- make sure soil moisture layers and model layers are consistent.
!-- n_msoil_layers must be equal to n_dsoil_layers from data table.
    IF ( n_msoil_layers > n_dsoil_layers )  THEN
       n_msoil_layers = n_dsoil_layers
       WRITE( message_string, * ) 'There are more soil moisture levels (', n_msoil_layers,         &
                                  ') than model soil layers (', n_dsoil_layers, ').&',             &
                                  'Biogenic model defaults to model layers.'
       CALL message( 'gamma_dstress', 'CHM0032', 0, 1, 0, 6, 0 )
    ELSEIF ( n_msoil_layers < n_dsoil_layers )  THEN
       WRITE( message_string, * ) 'There are less soil moisture levels (', n_msoil_layers,         &
                                  ') than model soil layers (', n_dsoil_layers, ').&',             &
                                  'Biogenic model will use soil moisture layers.'
       CALL message( 'gamma_dstress', 'CHM0033', 0, 1, 0, 6, 0 )
    ENDIF

    ALLOCATE( soil_m(n_msoil_layers)           )
    ALLOCATE( relwater_content(n_msoil_layers) )
    ALLOCATE( root_dist(n_msoil_layers)        )
    ALLOCATE( root_frac(n_msoil_layers)        )
    ALLOCATE( soil_frac(n_msoil_layers)        )

    soil_m(1:n_msoil_layers) = soil_moisture(0:n_msoil_layers -1)

    dstress_this_layer = 0.0
    gam_dstress        = 0.0
    root_frac          = 0.0

    qw = ltab_qw(soil_type_id)
    qc = ltab_qc(soil_type_id)

    all_soil_layer_depth = sum(ltab_d_soil)
    cumulative_soil_depth = 0.0

    DO  isl = 1, n_msoil_layers
       cumulative_soil_depth = cumulative_soil_depth + ltab_d_soil(isl)
       root_dist (isl) = get_root_distribution( cumulative_soil_depth, pft_msoil )
    ENDDO
!
!-- scale root fraction from 4 layers to the given soil moisture layers
    iscale0 = 0
    iscale1 = 1

    DO  isl = 1, n_msoil_layers
       iscale0  = iscale0 + 1
       IF ( (isl > 1) )  THEN
          IF ( root_dist(isl) /= root_dist(isl-1) )  THEN
             root_frac(iscale1:isl - 1) = root_dist(isl-1)/(iscale0 - 1)
             iscale1 = isl
             iscale0 = 1
          ELSEIF ( isl == n_msoil_layers )  THEN
             root_frac(iscale1:isl) = root_dist(isl) / (iscale0)
          ENDIF
       ENDIF
    ENDDO
!
!-- calculate draught stress factor (weighted) soil moisture  method
    DO  isl = 1, n_msoil_layers
       soil_frac(isl) = ltab_d_soil(isl)/all_soil_layer_depth
       relwater_content(isl) = (soil_m(isl) - qw) / (qc - qw)
       dstress_this_layer =  (soil_frac(isl) + root_frac(isl)) * 0.5_wp * relwater_content(isl)
       gam_dstress = gam_dstress + dstress_this_layer
    ENDDO

    IF ( gam_dstress > 1.0 )  gam_dstress = 1.0
!
!-- deacllocate local arrays
    DEALLOCATE( root_dist )
    DEALLOCATE( root_frac )
    DEALLOCATE( soil_frac )
    DEALLOCATE( soil_m )
    DEALLOCATE( relwater_content )

 END FUNCTION gamma_dstress

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This function gets the root distribution value of the given PFT for the respective soil layer
!> data : https://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#root_fraction
!--------------------------------------------------------------------------------------------------!
 FUNCTION get_root_distribution( cfreq, pft4root ) RESULT( root_dist )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  pft4root  !< index day of the year

    REAL(wp), INTENT(IN) ::  cfreq  !< day of the year

    REAL(wp) ::  root_dist  !< root distribution (units m)


    IF ( cfreq <= ltab_root_layer_depth(0) )  THEN
          root_dist = ltab_root_dist(0, pft4root)
    ENDIF

    IF ( ( cfreq > ltab_root_layer_depth(0) )  .AND.                                               &
         ( cfreq <= ltab_root_layer_depth(1) ) )  THEN
          root_dist = ltab_root_dist(1, pft4root)
    ENDIF

    IF ( ( cfreq > ltab_root_layer_depth(1) )  .AND.                                               &
         ( cfreq <= ltab_root_layer_depth(2) ) )  THEN
          root_dist = ltab_root_dist(2, pft4root)
    ENDIF

    IF ( cfreq > ltab_root_layer_depth (2) )  THEN
          root_dist = ltab_root_dist(3, pft4root)
    ENDIF

 END FUNCTION get_root_distribution

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Seasonal dependence (gamma_sn) algorithm for the biogenic species.
!> This algorithm is based on Keenan et al. (2009)
!--------------------------------------------------------------------------------------------------!
 FUNCTION gamma_sn( ) RESULT( gam_sn )

    USE control_parameters,                                                                        &
        ONLY:  time_since_reference_point

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    USE chem_modules,                                                                              &
        ONLY:  ebio_max_emis_day

    IMPLICIT NONE

    INTEGER(iwp) ::  day_of_year  !< day of the year
    INTEGER(iwp) ::  index_dd     !< index day of the year

    REAL(wp) ::  gam_sn  !< result of seasonal dependence calc

    REAL(wp), PARAMETER ::  tau = 100.0_wp  !< breadth/kurtosis

    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year )
    index_dd = day_of_year

!-- calcualte seasonal dependence of bvoc emissions
    gam_sn = (EXP(-(index_dd - ebio_max_emis_day) / tau ))**2

 END FUNCTION gamma_sn

! -------------------------------------------------------------------------------------------------!
!> Description:
!> ------------
!>  Algorithm for coversion of lad units from m2/m3 to g/m3.
!--------------------------------------------------------------------------------------------------!
 FUNCTION lad_unitconv( lad_3d, iplant_spc, slm_default ) RESULT( lad_gm3 )

    USE chem_emis_biogenic_data_mod,                                                               &
        ONLY:  ltab_tree_spc_slm

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  iplant_spc  !< plant spcs- (single tree/pft)


    REAL(wp), INTENT(IN) ::  lad_3d       !< leaf area density 3D
    REAL(wp), INTENT(IN) ::  slm_default  !< default specific leaf mass units cm2 g-1

    REAL(wp) ::  lad_gm3          !< Output of matching routine with number
    REAL(wp) ::  slm_to_g_per_m2  !< result of lad_unitconv in units g/m2

    slm_to_g_per_m2 = slm_default

!-- convert m2/m3 to g/m3
    IF (iplant_spc > 0)  THEN
       slm_to_g_per_m2 = 10000.0_wp / ltab_tree_spc_slm(iplant_spc)
    ELSEIF (iplant_spc < 0)  THEN
       slm_to_g_per_m2 = 10000.0_wp / ltab_pft_slm(-iplant_spc)
    ENDIF

    lad_gm3 = lad_3d * slm_to_g_per_m2

 END FUNCTION lad_unitconv

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!>  Function to convert default emission factor in ug m-2 h-1, Guenther (2012) and Henrot (2017)
!> to umoles m-2 s-1.
!--------------------------------------------------------------------------------------------------!
 FUNCTION conv_grams2moles(k_spcs) RESULT( grams2moles )

    IMPLICIT NONE

    INTEGER,  INTENT(IN) ::  k_spcs  !< bvoc spcs

    REAL(KIND=dp) ::  grams2moles    !< conversion from grams moles

    REAL(wp) ::  mwt  !< molecular weight of the bvoc

    REAL(wp), PARAMETER ::  secs_in_one_hour = 3600.0_wp

    mwt = ltab_ebio_mwt(matched_spcs(k_spcs)%mech_index )
    grams2moles = 1.0/(mwt * secs_in_one_hour)

 END FUNCTION  conv_grams2moles


 END MODULE chem_emis_biogenic_mod
