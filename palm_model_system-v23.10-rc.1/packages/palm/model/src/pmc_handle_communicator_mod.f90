!> @file pmc_handle_communicator_mod.f90
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
! Description:
! ------------
! Handle MPI communicator in PALM model coupler
!--------------------------------------------------------------------------------------------------!
 MODULE pmc_handle_communicator

#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  message_string

    USE kinds

    USE MPI

    USE pmc_general,                                                                               &
        ONLY:  pmc_max_models,                                                                     &
               pmc_status_error,                                                                   &
               pmc_status_ok

    IMPLICIT NONE

!
!-- ATTENTION: Do not change the order of variable declarations in the following TYPE definition,
!-- because the order must follow the order in which the layout data is given in the NAMELIST file.
    TYPE pmc_layout

       CHARACTER(LEN=32) ::  name  !<

       INTEGER ::  id            !<
       INTEGER ::  parent_id     !<
       INTEGER ::  npe_total     !<

       REAL(wp) ::  lower_left_x  !<
       REAL(wp) ::  lower_left_y  !<
       REAL(wp) ::  nest_shift_z  !<

    END TYPE pmc_layout

    PUBLIC  pmc_status_ok, pmc_status_error  !<

    INTEGER, PARAMETER, PUBLIC ::  pmc_error_npes        = 1  !< illegal number of processes
    INTEGER, PARAMETER, PUBLIC ::  pmc_namelist_error    = 2  !< error(s) in nesting_parameters namelist
    INTEGER, PARAMETER, PUBLIC ::  pmc_no_namelist_found = 3  !< no couple layout namelist found

    INTEGER ::  m_my_cpl_id   !< coupler id of this model
    INTEGER ::  m_ncpl        !< number of couplers given in nesting_parameters namelist
    INTEGER ::  m_parent_id   !< coupler id of parent of this model
    INTEGER ::  m_world_comm  !< global nesting communicator

    TYPE(pmc_layout), PUBLIC, DIMENSION(pmc_max_models) ::  m_couplers  !< information of all couplers

    INTEGER, PUBLIC ::  m_model_comm          !< communicator of this model
    INTEGER, PUBLIC ::  m_model_npes          !<
    INTEGER, PUBLIC ::  m_model_rank          !<
    INTEGER, PUBLIC ::  m_to_parent_comm      !< communicator to the parent
    INTEGER, PUBLIC ::  m_world_rank          !<
    INTEGER         ::  m_parent_remote_size  !< number of processes in the parent model
    INTEGER         ::  m_world_npes          !<
    INTEGER         ::  peer_comm             !< peer_communicator for inter communicators

    INTEGER, DIMENSION(pmc_max_models), PUBLIC ::  m_to_child_comm   !< communicator to the child(ren)
    INTEGER, DIMENSION(:), POINTER, PUBLIC ::  pmc_parent_for_child  !<


    INTERFACE pmc_get_model_info
       MODULE PROCEDURE pmc_get_model_info
    END INTERFACE pmc_get_model_info

    PUBLIC pmc_get_model_info, pmc_init_model   !<



 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @Todo: Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_init_model( comm, nesting_bounds, nesting_datatransfer_mode, nesting_mode,         &
                            anterpolation_buffer_width, anterpolation_starting_height,             &
                            homogeneous_initialization_child, particle_coupling, pmc_status )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    USE pegrid,                                                                                    &
        ONLY:  myid

    IMPLICIT NONE

    CHARACTER(LEN=14), INTENT(INOUT) ::  nesting_bounds             !<
    CHARACTER(LEN=7),  INTENT(INOUT) ::  nesting_datatransfer_mode  !<
    CHARACTER(LEN=8),  INTENT(INOUT) ::  nesting_mode               !<

    INTEGER, INTENT(INOUT) ::  anterpolation_buffer_width  !< Boundary buffer width for anterpolation
    INTEGER, INTENT(INOUT) ::  comm                        !<
    INTEGER, INTENT(INOUT) ::  pmc_status                  !<

    REAL(wp), INTENT(INOUT) ::  anterpolation_starting_height  !< steering parameter for canopy restricted anterpolation

    LOGICAL, INTENT(INOUT)  ::  homogeneous_initialization_child !< switch to control initialization of child domains (default .FALSE.)
    LOGICAL, INTENT(INOUT)  ::  particle_coupling                !< switch for particle coupling (default .TRUE.)
    
    INTEGER ::  childcount     !<
    INTEGER ::  i              !<
    INTEGER ::  ierr           !<
    INTEGER ::  istat          !<
    INTEGER ::  m_my_cpl_rank  !<
    INTEGER ::  tag            !<

    INTEGER, DIMENSION(pmc_max_models)   ::  activeparent  !< I am active parent for this child ID
    INTEGER, DIMENSION(pmc_max_models+1) ::  start_pe      !<

    pmc_status   = pmc_status_ok
    comm         = -1
    m_world_comm = MPI_COMM_WORLD
    m_my_cpl_id  = -1
    childcount   =  0
    activeparent = -1
    start_pe(:)  =  0

    CALL MPI_COMM_RANK( MPI_COMM_WORLD, m_world_rank, istat )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, m_world_npes, istat )
!
!-- Only process 0 of root model reads
    IF ( m_world_rank == 0 )  THEN

       CALL read_coupling_layout( nesting_bounds, nesting_datatransfer_mode, nesting_mode,         &
                                  anterpolation_buffer_width, anterpolation_starting_height,       &
                                  homogeneous_initialization_child, particle_coupling, pmc_status )

       IF ( pmc_status /= pmc_no_namelist_found  .AND.                                             &
            pmc_status /= pmc_namelist_error )                                                     &
       THEN
!
!--       Calculate start PE of every model
          start_pe(1) = 0
          DO  i = 2, m_ncpl+1
             start_pe(i) = start_pe(i-1) + m_couplers(i-1)%npe_total
          ENDDO

!
!--       The sum of numbers of processes requested by all the domains must be equal to the total
!--       number of processes of the run
          IF ( start_pe(m_ncpl+1) /= m_world_npes )  THEN
             WRITE( message_string, '(2A,I6,2A,I6,A)' )                                            &
                                                'nesting-setup requires different number of ',     &
                                                'MPI procs (', start_pe(m_ncpl+1), ') than ',      &
                                                'provided (', m_world_npes,')'
             CALL message( 'pmc_init_model', 'PMC0001', 3, 2, 0, 6, 0 )
          ENDIF

       ENDIF

    ENDIF
!
!-- Broadcast the read status. This synchronises all other processes with process 0 of the root
!-- model. Without synchronisation, they would not behave in the correct way (e.g. they would not
!-- return in case of a missing NAMELIST).
    CALL MPI_BCAST( pmc_status, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat )

    IF ( pmc_status == pmc_no_namelist_found )  THEN
!
!--    Not a nested run; return the MPI_WORLD communicator
       comm = MPI_COMM_WORLD
       RETURN

    ELSEIF ( pmc_status == pmc_namelist_error )  THEN
!
!--    Only the root model gives the error message. Others are aborted by the message-routine with
!--    MPI_ABORT. Must be done this way since myid and comm2d have not yet been assigned at this
!--    point.
       IF ( m_world_rank == 0 )  THEN
          message_string = 'errors in &nesting_parameters'
          CALL message( 'pmc_init_model', 'PMC0002', 3, 2, 0, 6, 0 )
       ENDIF

    ENDIF

    CALL MPI_BCAST( m_ncpl,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat )
    CALL MPI_BCAST( start_pe, m_ncpl+1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat )
!
!-- Broadcast coupling layout
    DO  i = 1, m_ncpl
       CALL MPI_BCAST( m_couplers(i)%name, LEN( m_couplers(i)%name ),                              &
                       MPI_CHARACTER, 0, MPI_COMM_WORLD, istat )
       CALL MPI_BCAST( m_couplers(i)%id,           1, MPI_INTEGER, 0,                              &
                       MPI_COMM_WORLD, istat )
       CALL MPI_BCAST( m_couplers(i)%Parent_id,    1, MPI_INTEGER, 0,                              &
                       MPI_COMM_WORLD, istat )
       CALL MPI_BCAST( m_couplers(i)%npe_total,    1, MPI_INTEGER, 0,                              &
                       MPI_COMM_WORLD, istat )
       CALL MPI_BCAST( m_couplers(i)%lower_left_x, 1, MPI_REAL,    0,                              &
                       MPI_COMM_WORLD, istat )
       CALL MPI_BCAST( m_couplers(i)%lower_left_y, 1, MPI_REAL,    0,                              &
                       MPI_COMM_WORLD, istat )
       CALL MPI_BCAST( m_couplers(i)%nest_shift_z, 1, MPI_REAL,    0,                              &
                       MPI_COMM_WORLD, istat )
    ENDDO
    CALL MPI_BCAST( nesting_bounds, LEN( nesting_bounds ), MPI_CHARACTER, 0, MPI_COMM_WORLD, istat )
    CALL MPI_BCAST( nesting_mode, LEN( nesting_mode ), MPI_CHARACTER, 0, MPI_COMM_WORLD, istat )
    CALL MPI_BCAST( nesting_datatransfer_mode, LEN(nesting_datatransfer_mode), MPI_CHARACTER, 0,   &
                    MPI_COMM_WORLD, istat )
    CALL MPI_BCAST( anterpolation_buffer_width, 1, MPI_INT, 0, MPI_COMM_WORLD, istat )
    CALL MPI_BCAST( anterpolation_starting_height, 1, MPI_REAL, 0, MPI_COMM_WORLD, istat )
    CALL MPI_BCAST( homogeneous_initialization_child, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat )
    CALL MPI_BCAST( particle_coupling, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat )
!
!-- Assign global MPI processes to individual models by setting the couple id
    DO  i = 1, m_ncpl
       IF ( m_world_rank >= start_pe(i)  .AND.  m_world_rank < start_pe(i+1) )  THEN
          m_my_cpl_id = i
          EXIT
       ENDIF
    ENDDO
    m_my_cpl_rank = m_world_rank - start_pe(i)
!
!-- MPI_COMM_WORLD is the communicator for ALL models (MPI-1 approach). The communictors for the
!-- individual models as created by MPI_COMM_SPLIT. The color of the model is represented by the
!-- coupler id
    CALL MPI_COMM_SPLIT( MPI_COMM_WORLD, m_my_cpl_id, m_my_cpl_rank, comm, istat )
!
!-- Get size and rank of the model running on this process
    CALL  MPI_COMM_RANK( comm, m_model_rank, istat )
    CALL  MPI_COMM_SIZE( comm, m_model_npes, istat )
!
!-- Broadcast (from process 0) the parent id and id of every model
    DO  i = 1, m_ncpl
       CALL MPI_BCAST( m_couplers(i)%parent_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat )
       CALL MPI_BCAST( m_couplers(i)%id,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat )
    ENDDO
!
!-- Save the current model communicator for pmc internal use
    m_model_comm = comm

!
!-- Create intercommunicator between the parent and children.
!-- MPI_INTERCOMM_CREATE creates an intercommunicator between 2 groups of different colors. The
!-- grouping was done above with MPI_COMM_SPLIT. A duplicate of MPI_COMM_WORLD is created and used
!-- as peer communicator (peer_comm) for MPI_INTERCOMM_CREATE.
    CALL MPI_COMM_DUP( MPI_COMM_WORLD, peer_comm, ierr )
    DO  i = 2, m_ncpl
       IF ( m_couplers(i)%parent_id == m_my_cpl_id )  THEN
!
!--       Identify all children models of the current model and create inter-communicators to
!--       connect between the current model and its children models.
          tag = 500 + i
          CALL MPI_INTERCOMM_CREATE( comm, 0, peer_comm, start_pe(i), tag, m_to_child_comm(i),     &
                                     istat )
          childcount = childcount + 1
          activeparent(i) = 1
       ELSEIF ( i == m_my_cpl_id)  THEN
!
!--       Create an inter-communicator to connect between the current model and its parent model.
          tag = 500 + i
          CALL MPI_INTERCOMM_CREATE( comm, 0, peer_comm, start_pe(m_couplers(i)%parent_id), tag,   &
                                     m_to_parent_comm, istat )
       ENDIF
    ENDDO
!
!-- If I am a parent, count the number of children I have. Although this loop is symmetric on all
!-- processes, the "activeparent" flag is true (==1) on the respective individual process only.
    ALLOCATE( pmc_parent_for_child(childcount+1) )

    childcount = 0
    DO  i = 2, m_ncpl
       IF ( activeparent(i) == 1 )  THEN
          childcount = childcount + 1
          pmc_parent_for_child(childcount) = i
       ENDIF
    ENDDO
!
!-- Get the size of the parent model
    IF ( m_my_cpl_id > 1 )  THEN
       CALL MPI_COMM_REMOTE_SIZE( m_to_parent_comm, m_parent_remote_size, istat )
    ELSE
!
!--    The root model does not have a parent
       m_parent_remote_size = -1
    ENDIF
!
!-- Set myid to non-zero value except for the root domain. This is a setting for the message routine
!-- which is called at the end of pmci_init. That routine outputs messages for myid = 0, only.
!-- However, myid has not been assigened so far, so that all processes of the root model would
!-- output a message. To avoid this, set myid to some other value except for process 0 of the root
!-- domain.
    IF ( m_world_rank /= 0 )  myid = 1

 END SUBROUTINE pmc_init_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @Todo: Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmc_get_model_info( comm_world_nesting, cpl_id, cpl_name, cpl_parent_id, lower_left_x, &
                                lower_left_y, nest_shift_z, ncpl, npe_total, request_for_cpl_id,   &
                                atmosphere_ocean_coupled_run, parent_lower_left_x,                 &
                                parent_lower_left_y, parent_nest_shift_z, root_model )
!
!-- Provide module private variables of the pmc for PALM
    USE kinds

    USE, INTRINSIC ::  IEEE_ARITHMETIC,                                                            &
        ONLY:  IEEE_VALUE,                                                                         &
               IEEE_QUIET_NAN

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL ::  cpl_name             !<

    INTEGER, INTENT(IN), OPTIONAL  ::  request_for_cpl_id            !<

    INTEGER, INTENT(OUT), OPTIONAL ::  comm_world_nesting            !<
    INTEGER, INTENT(OUT), OPTIONAL ::  cpl_id                        !<
    INTEGER, INTENT(OUT), OPTIONAL ::  cpl_parent_id                 !<
    INTEGER, INTENT(OUT), OPTIONAL ::  ncpl                          !<
    INTEGER, INTENT(OUT), OPTIONAL ::  npe_total                     !<

    LOGICAL, INTENT(OUT), OPTIONAL ::  atmosphere_ocean_coupled_run  !<
    LOGICAL, INTENT(OUT), OPTIONAL ::  root_model                    !<

    REAL(wp), INTENT(OUT), OPTIONAL ::  lower_left_x         !<
    REAL(wp), INTENT(OUT), OPTIONAL ::  lower_left_y         !<
    REAL(wp), INTENT(OUT), OPTIONAL ::  nest_shift_z         !<
    REAL(wp), INTENT(OUT), OPTIONAL ::  parent_lower_left_x  !<
    REAL(wp), INTENT(OUT), OPTIONAL ::  parent_lower_left_y  !<
    REAL(wp), INTENT(OUT), OPTIONAL ::  parent_nest_shift_z  !<


    INTEGER(iwp) ::  i                                !<    
    INTEGER(iwp) ::  requested_cpl_id                 !<
    INTEGER(iwp) ::  parent_cpl_id                    !<

    REAL(wp) ::  one = 1.0_wp  !<


!
!-- Set the requested coupler id
    IF ( PRESENT( request_for_cpl_id ) )  THEN
       requested_cpl_id = request_for_cpl_id
!
!--    Check for allowed range of values
       IF ( requested_cpl_id < 1  .OR.  requested_cpl_id > m_ncpl )  RETURN
    ELSE
       requested_cpl_id = m_my_cpl_id
    ENDIF
!
!-- Return the requested information
    IF ( PRESENT( comm_world_nesting )  )  THEN
       comm_world_nesting = m_world_comm
    ENDIF
    IF ( PRESENT( cpl_id )        )  THEN
       cpl_id = requested_cpl_id
    ENDIF
    IF ( PRESENT( cpl_parent_id ) )  THEN
       cpl_parent_id = m_couplers(requested_cpl_id)%parent_id
    ENDIF
    IF ( PRESENT( cpl_name )      )  THEN
       cpl_name = m_couplers(requested_cpl_id)%name
    ENDIF
    IF ( PRESENT( ncpl )          )  THEN
       ncpl = m_ncpl
    ENDIF
    IF ( PRESENT( npe_total )     )  THEN
       npe_total = m_couplers(requested_cpl_id)%npe_total
    ENDIF
    IF ( PRESENT( lower_left_x )  )  THEN
       lower_left_x = m_couplers(requested_cpl_id)%lower_left_x
    ENDIF
    IF ( PRESENT( lower_left_y )  )  THEN
       lower_left_y = m_couplers(requested_cpl_id)%lower_left_y
    ENDIF
    IF ( PRESENT( nest_shift_z )  )  THEN
       nest_shift_z = m_couplers(requested_cpl_id)%nest_shift_z
    ENDIF
!
!-- Get lower left edge coordinates of parent model
    parent_cpl_id = m_couplers(requested_cpl_id)%parent_id
    IF ( PRESENT( parent_lower_left_x )  )  THEN
       IF ( parent_cpl_id /= -1 )  THEN
          parent_lower_left_x = m_couplers(parent_cpl_id)%lower_left_x
       ELSE
          parent_lower_left_x = IEEE_VALUE( one, IEEE_QUIET_NAN )
       ENDIF
    ENDIF
    IF ( PRESENT( parent_lower_left_y )  )  THEN
       IF ( parent_cpl_id /= -1 )  THEN
          parent_lower_left_y = m_couplers(parent_cpl_id)%lower_left_y
       ELSE
          parent_lower_left_y = IEEE_VALUE( one, IEEE_QUIET_NAN )
       ENDIF
    ENDIF
    IF ( PRESENT( parent_nest_shift_z )  )  THEN
       IF ( parent_cpl_id /= -1 )  THEN
          parent_nest_shift_z = m_couplers(parent_cpl_id)%nest_shift_z
       ELSE
          parent_nest_shift_z = IEEE_VALUE( one, IEEE_QUIET_NAN )
       ENDIF
    ENDIF

    IF ( PRESENT( atmosphere_ocean_coupled_run ) )  THEN
       atmosphere_ocean_coupled_run = .FALSE.
       DO  i = 1, SIZE( m_couplers )
          IF ( m_couplers(i)%name(1:5) == 'ocean' )  THEN
             atmosphere_ocean_coupled_run = .TRUE.
             EXIT
          ENDIF
       ENDDO
    ENDIF

    IF ( PRESENT( root_model ) )  THEN
       root_model = ( m_my_cpl_id == 1 )
    ENDIF

 END SUBROUTINE pmc_get_model_info


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @Todo: Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE read_coupling_layout( nesting_bounds, nesting_datatransfer_mode, nesting_mode,         &
                                  anterpolation_buffer_width, anterpolation_starting_height,       &
                                  homogeneous_initialization_child, particle_coupling, pmc_status )

    IMPLICIT NONE

    CHARACTER(LEN=14), INTENT(INOUT) ::  nesting_bounds             !<
    CHARACTER(LEN=7),  INTENT(INOUT) ::  nesting_datatransfer_mode  !<
    CHARACTER(LEN=8),  INTENT(INOUT) ::  nesting_mode               !<

    INTEGER, INTENT(INOUT)      ::  anterpolation_buffer_width  !< Boundary buffer width for anterpolation
    INTEGER(iwp), INTENT(INOUT) ::  pmc_status                  !<

    REAL(wp), INTENT(INOUT) ::  anterpolation_starting_height   !< steering parameter for canopy restricted anterpolation

    LOGICAL, INTENT(INOUT)  ::  homogeneous_initialization_child !< switch to control initialization of child domains (default .FALSE.)
    LOGICAL, INTENT(INOUT)  ::  particle_coupling                !< switch for particle coupling (default .TRUE.)

    INTEGER(iwp) ::  bad_llcorner  !<
    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  io_status     !<

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    TYPE(pmc_layout), DIMENSION(pmc_max_models) ::  domain_layouts  !<

    NAMELIST /nesting_parameters/  anterpolation_buffer_width,                                     &
                                   anterpolation_starting_height,                                  &
                                   domain_layouts,                                                 &
                                   homogeneous_initialization_child,                               &
                                   nesting_bounds,                                                 &
                                   nesting_datatransfer_mode,                                      &
                                   nesting_mode,                                                   &
                                   particle_coupling,                                              &
                                   switch_off_module



!
!-- Initialize some coupling variables
    domain_layouts(1:pmc_max_models)%id = -1
    m_ncpl =   0

    pmc_status = pmc_no_namelist_found
!
!-- Open the NAMELIST-file and read the nesting layout
    CALL check_open( 11 )
    READ ( 11, nesting_parameters, IOSTAT = io_status )
!
!-- Set filepointer to the beginning of the file. Otherwise process 0 will later be unable to read
!-- the inipar-NAMELIST
    REWIND ( 11 )

    IF ( io_status == 0 )  THEN
!
!--    nesting_parameters namelist was found and read correctly. Enable the nesting by setting
!--    the palm model coupler status respectively.
       IF ( .NOT. switch_off_module )  THEN
          pmc_status = pmc_status_ok
       ELSE
          RETURN
       ENDIF
    ELSEIF ( io_status < 0 )  THEN
!
!--    No nesting_parameters-NAMELIST found
       RETURN
    ELSEIF ( io_status > 0 )  THEN
!
!--    Errors in reading nesting_parameters-NAMELIST.
!--    Try, if the deprecated domain_layouts format has been used
       CALL read_deprecated_nesting_namelist( nesting_bounds, nesting_datatransfer_mode,           &
                                              nesting_mode, anterpolation_buffer_width,            &
                                              anterpolation_starting_height,                       &
                                              homogeneous_initialization_child,                    &
                                              particle_coupling, io_status, switch_off_module,     &
                                              domain_layouts )
       IF ( io_status == 0 )  THEN
          message_string = 'deprecated format for domain_layouts is used'
          CALL message( 'read_coupling_layout', 'PMC0032', 0, 1, 0, 6, 0 )
          IF ( .NOT. switch_off_module )  THEN
             pmc_status = pmc_status_ok
          ELSE
             RETURN
          ENDIF
       ELSE
          pmc_status = pmc_namelist_error
          RETURN
       ENDIF
    ENDIF

!
!-- Check that at least one child has been defined.
    IF ( domain_layouts(2)%id == -1 )  THEN
       message_string = 'no child defined'
       CALL message( 'read_coupling_layout', 'PMC0030', 3, 2, 0, 6, 0 )
    ENDIF

!
!-- Output location message
    CALL location_message( 'initialize communicators for nesting', 'start' )
!
!-- Assign the layout to the corresponding internally used variable m_couplers
    m_couplers = domain_layouts
!
!-- Get the number of nested models given in the nesting_parameters-NAMELIST
    DO  i = 1, pmc_max_models
!
!--    When id=-1 is found for the first time, the list of domains is finished
       IF ( m_couplers(i)%id == -1  .OR.  i == pmc_max_models )  THEN
          IF ( m_couplers(i)%id == -1 )  THEN
             m_ncpl = i - 1
             EXIT
          ELSE
             m_ncpl = pmc_max_models
          ENDIF
       ENDIF
    ENDDO
!
!-- Make sure that all domains have equal lower left corner in case of vertical nesting
    IF ( TRIM( nesting_bounds ) == 'vertical_only' )  THEN
       bad_llcorner = 0
       DO  i = 1, m_ncpl
          IF ( domain_layouts(i)%lower_left_x /= 0.0_wp .OR.                                       &
               domain_layouts(i)%lower_left_y /= 0.0_wp )  THEN
             bad_llcorner = bad_llcorner + 1
             domain_layouts(i)%lower_left_x = 0.0_wp
             domain_layouts(i)%lower_left_y = 0.0_wp
          ENDIF
       ENDDO
       IF ( bad_llcorner /= 0)  THEN
          WRITE( message_string, *)  'At least one dimension of lower ',                           &
                                     'left corner of one domain is not 0. ',                       &
                                     'All lower left corners were set to (0,0).'
          CALL message( 'read_coupling_layout', 'PMC0003', 0, 0, 0, 6, 0 )
       ENDIF
    ENDIF

    CALL location_message( 'initialize communicators for nesting', 'finished' )

 END SUBROUTINE read_coupling_layout


!
!> @TODO: To be removed in a later revision
!-- This routines reads the nesting_parameters namelist with deprecated domain_layouts format.
 SUBROUTINE read_deprecated_nesting_namelist( nesting_bounds, nesting_datatransfer_mode,           &
                                              nesting_mode, anterpolation_buffer_width,            &
                                              anterpolation_starting_height,                       &
                                              homogeneous_initialization_child,                    &
                                              particle_coupling, io_status, switch_off_module,     &
                                              domain_layouts_new )

    CHARACTER(LEN=14), INTENT(INOUT) ::  nesting_bounds             !<
    CHARACTER(LEN=7),  INTENT(INOUT) ::  nesting_datatransfer_mode  !<
    CHARACTER(LEN=8),  INTENT(INOUT) ::  nesting_mode               !<

    INTEGER, INTENT(INOUT)      ::  anterpolation_buffer_width  !< Boundary buffer width for anterpolation

    REAL(wp), INTENT(INOUT) ::  anterpolation_starting_height   !< steering parameter for canopy restricted anterpolation

    LOGICAL, INTENT(INOUT)  ::  homogeneous_initialization_child !< switch to control initialization of child domains (default .FALSE.)
    LOGICAL, INTENT(INOUT)  ::  particle_coupling                !< switch for particle coupling (default .TRUE.)

    INTEGER(iwp), INTENT(inout) ::  io_status     !<

    LOGICAL, INTENT(inout) ::  switch_off_module  !< local namelist parameter to switch off the module

    TYPE pmc_layout_old

       CHARACTER(LEN=32) ::  name  !<

       INTEGER ::  id            !<
       INTEGER ::  parent_id     !<
       INTEGER ::  npe_total     !<

       REAL(wp) ::  lower_left_x  !<
       REAL(wp) ::  lower_left_y  !<

    END TYPE pmc_layout_old

    TYPE(pmc_layout_old), DIMENSION(pmc_max_models) ::  domain_layouts  !<
    TYPE(pmc_layout), DIMENSION(pmc_max_models), INTENT(inout) ::  domain_layouts_new  !<


    NAMELIST /nesting_parameters/  anterpolation_buffer_width,                                     &
                                   anterpolation_starting_height,                                  &
                                   domain_layouts,                                                 &
                                   homogeneous_initialization_child,                               &
                                   nesting_bounds,                                                 &
                                   nesting_datatransfer_mode,                                      &
                                   nesting_mode,                                                   &
                                   particle_coupling,                                              &
                                   switch_off_module

    READ ( 11, nesting_parameters, IOSTAT = io_status )
!
!-- Set filepointer to the beginning of the file. Otherwise process 0 will later be unable to read
!-- the inipar-NAMELIST
    REWIND ( 11 )

    IF ( io_status == 0 )  THEN
       domain_layouts_new(1:pmc_max_models)%name         = domain_layouts(1:pmc_max_models)%name
       domain_layouts_new(1:pmc_max_models)%id           = domain_layouts(1:pmc_max_models)%id
       domain_layouts_new(1:pmc_max_models)%parent_id    = domain_layouts(1:pmc_max_models)%parent_id
       domain_layouts_new(1:pmc_max_models)%npe_total    = domain_layouts(1:pmc_max_models)%npe_total
       domain_layouts_new(1:pmc_max_models)%lower_left_x = domain_layouts(1:pmc_max_models)%lower_left_x
       domain_layouts_new(1:pmc_max_models)%lower_left_y = domain_layouts(1:pmc_max_models)%lower_left_y
       domain_layouts_new(1:pmc_max_models)%nest_shift_z = 0.0_wp
    ENDIF

 END SUBROUTINE read_deprecated_nesting_namelist

#endif
 END MODULE pmc_handle_communicator
