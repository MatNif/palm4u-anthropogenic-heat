!> @file init_pegrid.f90
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
! Description:
! ------------
!> Determination of the virtual processor topology (if not prescribed by the user) and computation
!> of the grid point number and array bounds of the local domains.
!> @todo: remove MPI-data types for 2D exchange on coarse multigrid level (not used any more)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_pegrid

#if defined( __parallel )
    USE MPI
#endif

    USE control_parameters,                                                                        &
        ONLY:  bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_lr,                                                                              &
               bc_ns,                                                                              &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               grid_level,                                                                         &
               grid_level_count,                                                                   &
               maximum_grid_level,                                                                 &
               message_string,                                                                     &
               mg_switch_to_pe0_level,                                                             &
               pe_grid_prescribed,                                                                 &
               psolver,                                                                            &
               serial_run

#if defined( __parallel )
    USE control_parameters,                                                                        &
        ONLY:  gathered_size,                                                                      &
               momentum_advec,                                                                     &
               nested_run,                                                                         &
               outflow_source_plane,                                                               &
               scalar_advec,                                                                       &
               subdomain_size,                                                                     &
               syn_turb_gen,                                                                       &
               turbulent_inflow,                                                                   &
               turbulent_outflow,                                                                  &
               use_sm_for_poisfft,                                                                 &
               y_shift

    USE grid_variables,                                                                            &
        ONLY:  dx
#endif

    USE indices,                                                                                   &
        ONLY:  nnx,                                                                                &
               nny,                                                                                &
               nnz,                                                                                &
               nx,                                                                                 &
               nxl,                                                                                &
               nxl_mg,                                                                             &
               nxlu,                                                                               &
               nxr,                                                                                &
               nxr_mg,                                                                             &
               ny,                                                                                 &
               nyn,                                                                                &
               nyn_mg,                                                                             &
               nys,                                                                                &
               nys_mg,                                                                             &
               nysv,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt,                                                                                &
               nzt_mg,                                                                             &
               topo_flags_1,                                                                       &
               topo_flags_2,                                                                       &
               topo_flags_3,                                                                       &
               topo_flags_4,                                                                       &
               topo_flags_5,                                                                       &
               topo_flags_6,                                                                       &
               topo_flags_7,                                                                       &
               topo_flags_8,                                                                       &
               topo_flags_9,                                                                       &
               topo_flags_10

#if defined( __parallel )
    USE indices,                                                                                   &
        ONLY:  mg_loc_ind,                                                                         &
               nbgp,                                                                               &
               nnx_pe,                                                                             &
               nny_pe,                                                                             &
               nxl_pe,                                                                             &
               nxr_pe,                                                                             &
               nyn_pe,                                                                             &
               nys_pe
#endif

    USE kinds

    USE pegrid

#if defined( __parallel )
    USE synthetic_turbulence_generator_mod,                                                        &
        ONLY:  id_stg_left,                                                                        &
               id_stg_north,                                                                       &
               id_stg_right,                                                                       &
               id_stg_south
#endif

    USE turbulent_inflow_mod,                                                                      &
        ONLY:  turbulent_inflow_method

    USE transpose_mod,                                                                             &
        ONLY:  nnx_x_max,                                                                          &
               nxl_y,                                                                              &
               nxl_z,                                                                              &
               nxr_x_max,                                                                          &
               nxr_y,                                                                              &
               nxr_z,                                                                              &
               nx_y_max,                                                                           &
               nyn_x,                                                                              &
               nyn_x_max,                                                                          &
               nyn_z,                                                                              &
               nys_x,                                                                              &
               nys_z,                                                                              &
               ny_z_max,                                                                           &
               nzb_x,                                                                              &
               nzb_y,                                                                              &
               nzt_x,                                                                              &
               nzt_y,                                                                              &
               nzt_y_max,                                                                          &
               nz_x_max

#if defined( __parallel )
    USE transpose_mod,                                                                             &
        ONLY:  nnx_y_max,                                                                          &
               nny_yd_max,                                                                         &
               nny_z_max,                                                                          &
               nnz_x_max,                                                                          &
               nnz_yd_max,                                                                         &
               nnz_z_max,                                                                          &
               nxl_y_pe,                                                                           &
               nxl_yd,                                                                             &
               nxr_y_max,                                                                          &
               nxr_y_pe,                                                                           &
               nxr_yd,                                                                             &
               nyn_yd_max,                                                                         &
               nyn_z_max,                                                                          &
               nyn_z_pe,                                                                           &
               nys_z_pe,                                                                           &
               nzb_x_pe,                                                                           &
               nzb_yd,                                                                             &
               nzb_yd_pe,                                                                          &
               nzt_x_max,                                                                          &
               nzt_x_pe,                                                                           &
               nzt_yd,                                                                             &
               nzt_yd_max,                                                                         &
               nzt_yd_pe,                                                                          &
               nz_yd_max

    USE sm_poisfft_mod,                                                                            &
        ONLY:  sm_poisfft
#endif

    IMPLICIT NONE

#if defined( __parallel )
    CHARACTER(LEN=7) ::  myid_char_prel = ''  !< preliminary processor id number
#endif

    INTEGER(iwp) ::  i                        !< running index over number of processors or number of multigrid level
#if defined( __parallel )
    INTEGER(iwp) ::  id_outflow_l             !< local value of id_outflow
    INTEGER(iwp) ::  id_outflow_source_l      !< local value of id_outflow_source
    INTEGER(iwp) ::  id_stg_left_l            !< left lateral boundary local core id in case of turbulence generator
    INTEGER(iwp) ::  id_stg_north_l           !< north lateral boundary local core id in case of turbulence generator
    INTEGER(iwp) ::  id_stg_right_l           !< right lateral boundary local core id in case of turbulence generator
    INTEGER(iwp) ::  id_stg_south_l           !< south lateral boundary local core id in case of turbulence generator
    INTEGER(iwp) ::  ind(5)                   !< array containing the subdomain bounds
#endif
    INTEGER(iwp) ::  j                        !< running index, used for various loops
    INTEGER(iwp) ::  k                        !< number of vertical grid points in different multigrid level
    INTEGER(iwp) ::  maximum_grid_level_l     !< maximum number of grid level without switching to PE 0
    INTEGER(iwp) ::  mg_levels_x              !< maximum number of grid level allowed along x-direction
    INTEGER(iwp) ::  mg_levels_y              !< maximum number of grid level allowed along y-direction
    INTEGER(iwp) ::  mg_levels_z              !< maximum number of grid level allowed along z-direction
    INTEGER(iwp) ::  mg_switch_to_pe0_level_l !< maximum number of grid level with switching to PE 0
#if defined( __parallel )
    INTEGER(iwp) ::  nnx_y                    !< quotient of number of grid points along x-direction and number of PEs used along y-direction
    INTEGER(iwp) ::  nny_x                    !< quotient of number of grid points along y-direction and number of PEs used along x-direction
    INTEGER(iwp) ::  nny_z                    !< quotient of number of grid points along y-direction and number of PEs used along x-direction
    INTEGER(iwp) ::  nnz_x                    !< quotient of number of grid points along z-direction and number of PEs used along x-direction
    INTEGER(iwp) ::  nnz_y                    !< quotient of number of grid points along z-direction and number of PEs used along x-direction
    INTEGER(iwp) ::  numproc_sqr              !< square root of the number of processors
#endif
    INTEGER(iwp) ::  nxl_l                    !< lower index bound along x-direction on subdomain and different multigrid level
    INTEGER(iwp) ::  nxr_l                    !< upper index bound along x-direction on subdomain and different multigrid level
    INTEGER(iwp) ::  nyn_l                    !< lower index bound along y-direction on subdomain and different multigrid level
    INTEGER(iwp) ::  nys_l                    !< upper index bound along y-direction on subdomain and different multigrid level
#if defined( __parallel )
    INTEGER(iwp) ::  nzb_l                    !< lower index bound along z-direction on subdomain and different multigrid level
#endif
    INTEGER(iwp) ::  nzt_l                    !< upper index bound along z-direction on subdomain and different multigrid level
!$  INTEGER(iwp) ::  omp_get_num_threads      !< number of OpenMP threads

#if defined( __parallel )
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ind_all !< dummy array containing index bounds on subdomain, used for gathering
    INTEGER(iwp), DIMENSION(2) ::  npe_xy = 1     !< number of processors along x-y dimension
    INTEGER(iwp)               ::  lcoord(2)      !< PE coordinates of left neighbor along x and y
    INTEGER(iwp)               ::  rcoord(2)      !< PE coordinates of right neighbor along x and y
    INTEGER(iwp)               ::  irest_x        !< remaining grid points along x in case of non-uniform subdomains
    INTEGER(iwp)               ::  irest_y        !< remaining grid points along y in case of non-uniform subdomains

    LOGICAL ::  nnx_too_small              !< global flag to check if subdomain size along x is too small
    LOGICAL ::  nnx_too_small_l = .FALSE.  !< local flag to check if subdomain size along x is too small
    LOGICAL ::  nny_too_small              !< global flag to check if subdomain size along y is too small
    LOGICAL ::  nny_too_small_l = .FALSE.  !< local flag to check if subdomain size along y is too small
#endif
!
!-- Get the number of OpenMP threads
    !$OMP PARALLEL
!$  threads_per_task = omp_get_num_threads()
    !$OMP END PARALLEL


#if defined( __parallel )

    CALL location_message( 'creating virtual PE grids + MPI derived data types', 'start' )

!
!-- Determine the processor topology or check it, if prescribed by the user
    IF ( npex == -1  .AND.  npey == -1 )  THEN

       IF ( psolver /= 'poisfft_sm' )  THEN
!
!--       Automatic determination of the topology.
          numproc_sqr = SQRT( REAL( numprocs, KIND=wp ) )
          npex        = MAX( numproc_sqr , 1 )
          DO  WHILE ( MOD( numprocs , npex ) /= 0 )
             npex = npex - 1
          ENDDO
          npey = numprocs / npex
!
!--       Set internal flag that shared memory Poisson-FFT-solver is not used.
          CALL sm_poisfft%sm_init_pegrid_poisfft( .FALSE. )
       ELSE
!
!--       Topology is defined by the shared memory Poisson-FFT-solver.
!--       First check requirements. Further requirements will be checked in sm_init_pegrid_poisfft.
          IF ( numprocs < 4 )  THEN
             message_string = 'psolver = "poisfft_sm" does not work for less than 4 PEs'
             CALL message( 'init_pegrid', 'PAC0231', 1, 2, 0, 6, 0 )
          ENDIF

!
!--       Set internal flags that shared memory Poisson-FFT-solver used.
          CALL sm_poisfft%sm_init_pegrid_poisfft( .TRUE. )
          use_sm_for_poisfft = .TRUE.
       ENDIF

    ELSEIF ( npex /= -1  .AND.  npey /= -1 )  THEN

       IF ( psolver == 'poisfft_sm' )  THEN
          message_string = 'psolver = "poisfft_sm" does not allow setting of npex and/or npey'
          CALL message( 'init_pegrid', 'PAC0232', 1, 2, 0, 6, 0 )
       ELSE
!
!--       Set internal flag that shared memory Poisson-FFT-solver is not used.
          CALL sm_poisfft%sm_init_pegrid_poisfft( .FALSE. )
       ENDIF

!
!--    Prescribed by user. Number of processors on the prescribed topology must be equal to the
!--    number of PEs available to the job
       IF ( ( npex * npey ) /= numprocs )  THEN
          WRITE( message_string, * ) 'number of PEs of the prescribed ', 'topology (', npex*npey,  &
                                     ') does not match & the number of ',                          &
                                     'PEs available to the job (', numprocs, ')'
          CALL message( 'init_pegrid', 'PAC0233', 1, 2, 0, 6, 0 )
       ENDIF
       pe_grid_prescribed = .TRUE.

    ELSE
!
!--    If the processor topology is prescribed by the user, the number of
!--    PEs must be given in both directions
       message_string = 'both values of "npex" and "npey" must be given in the namelist ' //       &
                        'parameter file'
       CALL message( 'init_pegrid', 'PAC0234', 1, 2, 0, 6, 0 )

    ENDIF
!
!-- Create four default MPI communicators for the 2d virtual PE grid. One of them will be used as
!-- the main communicator for this run, while others might be used for specific quantities like
!-- aerosol, chemical species, or passive scalars), if their horizontal boundary conditions shall
!-- be different from those of the other quantities (e.g. non-cyclic conditions for aerosols, and
!-- cyclic conditions for all others).
    npe_xy(1) = npex
    npe_xy(2) = npey
    DO  i = 1, 4

       IF ( i == 1 )  cyclic = (/  .TRUE., .TRUE.  /)   ! cyclic along x and y
       IF ( i == 2 )  cyclic = (/  .TRUE., .FALSE. /)   ! cyclic along x
       IF ( i == 3 )  cyclic = (/ .FALSE., .TRUE.  /)   ! cyllic along y
       IF ( i == 4 )  cyclic = (/ .FALSE., .FALSE. /)   ! non-cyclic

       CALL MPI_CART_CREATE( comm_palm, ndim, npe_xy, cyclic, reorder,                             &
                             communicator_configurations(i)%mpi_communicator, ierr )

       CALL MPI_CART_SHIFT( communicator_configurations(i)%mpi_communicator, 0, 1,                 &
                            communicator_configurations(i)%pleft,                                  &
                            communicator_configurations(i)%pright, ierr )

       CALL MPI_CART_SHIFT( communicator_configurations(i)%mpi_communicator, 1, 1,                 &
                            communicator_configurations(i)%psouth,                                 &
                            communicator_configurations(i)%pnorth, ierr )

    ENDDO

!
!-- If necessary, set horizontal boundary conditions to non-cyclic
    IF ( bc_lr /= 'cyclic' )  cyclic(1) = .FALSE.
    IF ( bc_ns /= 'cyclic' )  cyclic(2) = .FALSE.


!
!-- Set the main communicator (virtual pe grid) for this run
    IF ( bc_lr == 'cyclic'  .AND.  bc_ns == 'cyclic' )  i = 1
    IF ( bc_lr == 'cyclic'  .AND.  bc_ns /= 'cyclic' )  i = 2
    IF ( bc_lr /= 'cyclic'  .AND.  bc_ns == 'cyclic' )  i = 3
    IF ( bc_lr /= 'cyclic'  .AND.  bc_ns /= 'cyclic' )  i = 4

    comm2d = communicator_configurations(i)%mpi_communicator
    pleft  = communicator_configurations(i)%pleft
    pright = communicator_configurations(i)%pright
    psouth = communicator_configurations(i)%psouth
    pnorth = communicator_configurations(i)%pnorth

!
!-- Set rank and coordinates of the main communicator.
    myid_char_prel = myid_char
    CALL MPI_COMM_RANK( comm2d, myid, ierr )
    WRITE( myid_char, '(''_'',I6.6)' )  myid

    CALL MPI_CART_COORDS( comm2d, myid, ndim, pcoord, ierr )

!
!-- Check, for security reasons, if the PE-id has changed between the preliminary communicator
!-- comm_palm (see main routine palm) and the final communicator comm2d. This should not happen (see
!-- respective comments in main routine).
    IF ( myid_char /= myid_char_prel )  THEN
       message_string = 'mismatch in PE-id: preliminary is "' // myid_char_prel //                 &
                        '"&comm2d-ID is "' // myid_char // '"'
       CALL message( 'init_pegrid', 'PAC0235', 0, 1, myid, 6, 0 )
    ENDIF

!
!-- In case of cyclic boundary conditions, a y-shift at the boundaries in x-direction can be
!-- introduced via parameter y_shift. The shift is done by modifying the processor grid in such a
!-- way that processors located at the x-boundary communicate across it to processors with
!-- y-coordinate shifted by y_shift relative to their own. This feature can not be used in
!-- combination with an fft pressure solver. It has been implemented to counter the effect of streak
!-- structures in case of cyclic boundary conditions. For a description of these see Munters
!-- (2016; dx.doi.org/10.1063/1.4941912)
!--
!-- Get coordinates of left and right neighbor on PE grid
    IF ( y_shift /= 0 ) THEN
       IF ( bc_lr == 'cyclic' ) THEN
          IF ( TRIM( psolver ) /= 'multigrid' .AND.  TRIM( psolver ) /= 'multigrid_noopt')  THEN
             message_string = 'y_shift /= 0 requires a multigrid pressure solver '
             CALL message( 'init_pegrid', 'PAC0236', 1, 2, 0, 6, 0 )
          ENDIF

          CALL MPI_CART_COORDS( comm2d, pright, ndim, rcoord, ierr )
          CALL MPI_CART_COORDS( comm2d, pleft, ndim, lcoord, ierr )

!
!--       If the x(y)-coordinate of the right (left) neighbor is smaller (greater) than that of the
!--       calling process, then the calling process is located on the right (left) boundary of the
!--       processor grid. In that case, the y-coordinate of that neighbor is increased (decreased)
!--       by y_shift.
!--       The rank of the process with that coordinate is then inquired and the neighbor rank for
!--       MPI_SENDRECV, pright (pleft) is set to it.
!--       In this way, the calling process receives a new right (left) neighbor for all future
!--       MPI_SENDRECV calls. That neighbor has a y-coordinate of y+(-)y_shift, where y is the
!--       original right (left) neighbor's y-coordinate. The modulo-operation ensures that if the
!--       neighbor's y-coordinate exceeds the grid-boundary, it will be relocated to the opposite
!--       part of the grid cyclicly.
          IF ( rcoord(1) < pcoord(1) )  THEN
             rcoord(2) = MODULO( rcoord(2) + y_shift, npey )
             CALL MPI_CART_RANK( comm2d, rcoord, pright, ierr )
          ENDIF

          IF ( lcoord(1) > pcoord(1) )  THEN
             lcoord(2) = MODULO( lcoord(2) - y_shift, npey )
             CALL MPI_CART_RANK( comm2d, lcoord, pleft, ierr )
          ENDIF

       ELSE
!
!--       y-shift for non-cyclic boundary conditions is only possible for turbulent inflow using
!--       the recycling method.
          IF ( .NOT. turbulent_inflow  .OR.  TRIM( turbulent_inflow_method ) == 'read_from_file' ) &
          THEN
             message_string = 'y_shift /= 0 is only allowed for cyclic ' //                        &
                              'boundary conditions in both directions '  //                        &
                              '&or with turbulent_inflow_method = "recycle_...".'
             CALL message( 'init_pegrid', 'PAC0237', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       y-shift requires uniform subdomains along y.
          IF ( MOD( ny+1, npey ) /= 0 )  THEN
             message_string = 'y_shift requires that subdomains are uniform along y-direction'
             CALL message( 'init_pegrid', 'PAC0238', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

!
!-- Determine sub-topologies for transpositions
!-- Transposition from z to x:
    remain_dims(1) = .TRUE.
    remain_dims(2) = .FALSE.
    CALL MPI_CART_SUB( comm2d, remain_dims, comm1dx, ierr )
    CALL MPI_COMM_RANK( comm1dx, myidx, ierr )
!
!-- Transposition from x to y
    remain_dims(1) = .FALSE.
    remain_dims(2) = .TRUE.
    CALL MPI_CART_SUB( comm2d, remain_dims, comm1dy, ierr )
    CALL MPI_COMM_RANK( comm1dy, myidy, ierr )
!
!-- Compare, if the virtual shared memory solver layout fits the virtual PE grid.
!-- The PE rank of a shared memory block must fit myidy, the node rank of the shared memory block
!-- must fit myidx.
    IF ( psolver == 'poisfft_sm' )  CALL sm_poisfft%sm_check_layout

!
!-- Check if subdomains are uniform or non-uniform.
    CALL check_subdomain_uniformity

!
!-- The multigrid solver doesn't allow non-uniform subdomains, because then subdomains would have
!-- different numbers of coarsening levels.
    IF ( non_uniform_subdomain  .AND.  psolver(1:9) == 'multigrid' )  THEN
       message_string = 'multigrid-solver does not allow to use non-uniform subdomains'
       CALL message( 'init_pegrid', 'PAC0239', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Calculate array bounds along x-direction for every PE.
    ALLOCATE( nnx_pe(0:npex-1), nny_pe(0:npey-1), nxl_pe(0:npex-1),                    &
              nxr_pe(0:npex-1), nyn_pe(0:npey-1), nys_pe(0:npey-1) )

!
!-- Compute the local index bounds (nxl, nyr, nys, nyn) on each of the subdomains.
!-- Start with array bounds along x. Begin with calculating the minimum grid point
!-- number and the possible rest.
    nnx     = ( nx + 1 ) / npex
    irest_x = MOD( nx+1, npex )
!
!-- Left and right array bounds on each PE.
    DO  i = 0, npex-1
       IF ( i < irest_x )  THEN
!
!--       PEs with low id will get one additional grid point (nnx+1), i.e. the surplus points
!--       irest_x will be distributed among these PEs.
          nxl_pe(i) = i         * ( nnx + 1 )
          nxr_pe(i) = ( i + 1 ) * ( nnx + 1 ) - 1
       ELSE
!
!--       The remaining PEs with higher id get nnx points. Be aware of the larger size
!--       of the low id subdomains when calculating the left index bound.
          nxl_pe(i) = irest_x * ( nnx + 1 ) + ( i - irest_x ) * nnx
          nxr_pe(i) = nxl_pe(i) + nnx - 1
       ENDIF
!
!--    Grid points along x.
       nnx_pe(i) = nxr_pe(i) - nxl_pe(i) + 1
    ENDDO

!
!-- Now calculate array bounds in y-direction. Begin with calculating the minimum grid point
!-- number and the possible rest.
    nny     = ( ny + 1 ) / npey
    irest_y = MOD( ny+1, npey )
!
!-- South and north array bounds on each PE.
    DO  j = 0, npey-1
       IF ( j < irest_y )  THEN
!
!--       PEs with low id will get one additional grid point (nny+1), i.e. the surplus points
!--       irest_y will be distributed among these PEs.
          nys_pe(j) = j         * ( nny + 1 )
          nyn_pe(j) = ( j + 1 ) * ( nny + 1 ) - 1
       ELSE
!
!--       The remaining PEs with higher id get nny points. Be aware of the larger size
!--       of the low id subdomains when calculating the left index bound.
          nys_pe(j) = irest_y * ( nny + 1 ) + ( j - irest_y ) * nny
          nyn_pe(j) = nys_pe(j) + nny - 1
       ENDIF
!
!--    Grid points along y.
       nny_pe(j) = nyn_pe(j) - nys_pe(j) + 1
    ENDDO

!
!-- Local array bounds and grid points on this PE.
    nxl = nxl_pe(pcoord(1))
    nxr = nxr_pe(pcoord(1))
    nys = nys_pe(pcoord(2))
    nyn = nyn_pe(pcoord(2))
    nzb = 0
    nzt = nz
    nnz = nz

    nnx = nxr-nxl+1
    nny = nyn-nys+1

!
!-- Set switches to define if the PE is situated at the border of the virtual processor grid
    IF ( nxl == 0 )   left_border_pe  = .TRUE.
    IF ( nxr == nx )  right_border_pe = .TRUE.
    IF ( nys == 0 )   south_border_pe = .TRUE.
    IF ( nyn == ny )  north_border_pe = .TRUE.

!
!-- Calculate array bounds and gridpoint numbers for the transposed arrays (needed in the pressure
!-- solver and partly for spectra calculation).
    CALL setup_transpose_indices

!
!-- Collect index bounds from other PEs (to be written to Fortran I/O restart file later).
!-- Note: Information in n.._pe arrays above is given for virtual PE grid coordinates, not for
!--       the linear (COMM_WORLD) communicator. Therefore, data is gathered again here.
    ALLOCATE( hor_index_bounds(4,0:numprocs-1) )

    IF ( myid == 0 )  THEN

       hor_index_bounds(1,0) = nxl
       hor_index_bounds(2,0) = nxr
       hor_index_bounds(3,0) = nys
       hor_index_bounds(4,0) = nyn

!
!--    Receive data from all other PEs
       DO  i = 1, numprocs-1
          CALL MPI_RECV( ibuf, 4, MPI_INTEGER, i, MPI_ANY_TAG, comm2d, status, ierr )
          hor_index_bounds(:,i) = ibuf(1:4)
       ENDDO

    ELSE
!
!--    Send index bounds to PE0
       ibuf(1) = nxl
       ibuf(2) = nxr
       ibuf(3) = nys
       ibuf(4) = nyn
       CALL MPI_SEND( ibuf, 4, MPI_INTEGER, 0, myid, comm2d, ierr )

    ENDIF

!
!-- Determine the number of ghost point layers
    IF ( scalar_advec == 'ws-scheme'  .OR.  momentum_advec == 'ws-scheme'  .OR.  nested_run )  THEN
       nbgp = 3
    ELSE
       nbgp = 1
    ENDIF

!
!-- Check that the number of computational grid points is not smaller than the number of ghost
!-- points.
    IF ( nnx < nbgp )  nnx_too_small_l = .TRUE.
    CALL MPI_ALLREDUCE( nnx_too_small_l, nnx_too_small, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr )
    IF ( nnx_too_small )  THEN
       IF ( non_uniform_subdomain )  THEN
          WRITE( message_string, * ) 'number of subdomain grid points along x (', nnx-1, ') is ',  &
                                     'smaller than the number of ghost points (', nbgp, ')'
       ELSE
          WRITE( message_string, * ) 'number of subdomain grid points along x (', nnx, ') is ',    &
                                     'smaller than the number of ghost points (', nbgp, ')'
       ENDIF
       CALL message( 'init_pegrid', 'PAC0240', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( nny < nbgp )  nny_too_small_l = .TRUE.
    CALL MPI_ALLREDUCE( nny_too_small_l, nny_too_small, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr )
    IF ( nny_too_small )  THEN
       IF ( non_uniform_subdomain )  THEN
          WRITE( message_string, * ) 'number of subdomain grid points along y (', nny-1, ') is ',  &
                                     'smaller than the number of ghost points (', nbgp, ')'
       ELSE
          WRITE( message_string, * ) 'number of subdomain grid points along y (', nny, ') is ',    &
                                     'smaller than the number of ghost points (', nbgp, ')'
       ENDIF
       CALL message( 'init_pegrid', 'PAC0240', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Create a new MPI derived datatype for the exchange of surface (xy) data, which is needed for
!-- coupled atmosphere-ocean runs.
!-- First, calculate number of grid points of an xy-plane.
    ngp_xy  = ( nxr - nxl + 1 + 2 * nbgp ) * ( nyn - nys + 1 + 2 * nbgp )
    CALL MPI_TYPE_VECTOR( ngp_xy, 1, nzt-nzb+2, MPI_REAL, type_xy, ierr )
    CALL MPI_TYPE_COMMIT( type_xy, ierr )

    serial_run = .FALSE.
#else

!
!-- Else branch contains settings for serial run.
    serial_run = .TRUE.
    pe_grid_prescribed = .FALSE.
    npex = 1
    npey = 1

!
!-- Array bounds when running on a single PE (respectively a non-parallel machine)
    nxl = 0
    nxr = nx
    nnx = nxr - nxl + 1
    nys = 0
    nyn = ny
    nny = nyn - nys + 1
    nzb = 0
    nzt = nz
    nnz = nz

    ALLOCATE( hor_index_bounds(4,0:0) )
    hor_index_bounds(1,0) = nxl
    hor_index_bounds(2,0) = nxr
    hor_index_bounds(3,0) = nys
    hor_index_bounds(4,0) = nyn

!
!-- Array bounds for the pressure solver (in the parallel code, these bounds are the ones for the
!-- transposed arrays)
    nys_x = nys
    nyn_x = nyn
    nzb_x = nzb + 1
    nzt_x = nzt

    nxl_y = nxl
    nxr_y = nxr
    nzb_y = nzb + 1
    nzt_y = nzt

    nxl_z = nxl
    nxr_z = nxr
    nys_z = nys
    nyn_z = nyn

!
!-- Set dimension variables here to be compatible with non-uniform grid.
    nxr_x_max = nxr
    nnx_x_max = nnx
    nyn_x_max = nyn_x
    nx_y_max  = nx
    nz_x_max  = nz
    nzt_y_max = nzt_y
    ny_z_max  = ny

#endif

!
!-- Calculate number of grid levels necessary for the multigrid poisson solver as well as the
!-- gridpoint indices on each level
    IF ( psolver(1:9) == 'multigrid' )  THEN

!
!--    First calculate number of possible grid levels for the subdomains
       mg_levels_x = 1
       mg_levels_y = 1
       mg_levels_z = 1

       i = nnx
       DO WHILE ( MOD( i, 2 ) == 0  .AND.  i /= 2 )
          i = i / 2
          mg_levels_x = mg_levels_x + 1
       ENDDO

       j = nny
       DO WHILE ( MOD( j, 2 ) == 0  .AND.  j /= 2 )
          j = j / 2
          mg_levels_y = mg_levels_y + 1
       ENDDO

       k = nz    ! do not use nnz because it might be > nz due to transposition
                 ! requirements
       DO WHILE ( MOD( k, 2 ) == 0  .AND.  k /= 2 )
          k = k / 2
          mg_levels_z = mg_levels_z + 1
       ENDDO
!
!--    The optimized MG-solver does not allow odd values for nz at the coarsest grid level
       IF ( TRIM( psolver ) /= 'multigrid_noopt' )  THEN
          IF ( MOD( k, 2 ) /= 0 )  mg_levels_z = mg_levels_z - 1
!
!--       An odd value of nz does not work. The finest level must have an even value.
          IF (  mg_levels_z == 0 )  THEN
             message_string = 'optimized multigrid method requires nz to be even'
             CALL message( 'init_pegrid', 'PAC0241', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

       maximum_grid_level = MIN( mg_levels_x, mg_levels_y, mg_levels_z )
!
!--    Check if subdomain sizes prevents any coarsening.
!--    This case, the maximum number of grid levels is 1, i.e. effectively a Gauss-Seidel scheme is
!--    applied rather than a multigrid approach.
!--    Give a warning in this case.
       IF ( maximum_grid_level == 1  .AND.  mg_switch_to_pe0_level == -1 )  THEN
          message_string = 'no grid coarsening possible, multigrid ' //                            &
                           'approach effectively reduces to a Gauss-Seidel scheme'
          CALL message( 'poismg', 'PAC0242', 0, 1, 0, 6, 0 )
       ENDIF

!
!--    Find out, if the total domain allows more levels. These additional levels are identically
!--    processed on all PEs.
       IF ( numprocs > 1  .AND.  mg_switch_to_pe0_level /= -1 )  THEN

          IF ( mg_levels_z > MIN( mg_levels_x, mg_levels_y ) )  THEN

             mg_switch_to_pe0_level_l = maximum_grid_level

             mg_levels_x = 1
             mg_levels_y = 1

             i = nx+1
             DO WHILE ( MOD( i, 2 ) == 0  .AND.  i /= 2 )
                i = i / 2
                mg_levels_x = mg_levels_x + 1
             ENDDO

             j = ny+1
             DO WHILE ( MOD( j, 2 ) == 0  .AND.  j /= 2 )
                j = j / 2
                mg_levels_y = mg_levels_y + 1
             ENDDO

             maximum_grid_level_l = MIN( mg_levels_x, mg_levels_y, mg_levels_z )

             IF ( maximum_grid_level_l > mg_switch_to_pe0_level_l )  THEN
                mg_switch_to_pe0_level_l = maximum_grid_level_l - mg_switch_to_pe0_level_l + 1
             ELSE
                mg_switch_to_pe0_level_l = 0
             ENDIF

          ELSE

             mg_switch_to_pe0_level_l = 0
             maximum_grid_level_l = maximum_grid_level

          ENDIF

!
!--       Use switch level calculated above only if it is not pre-defined by user
          IF ( mg_switch_to_pe0_level == 0 )  THEN
             IF ( mg_switch_to_pe0_level_l /= 0 )  THEN
                mg_switch_to_pe0_level = mg_switch_to_pe0_level_l
                maximum_grid_level     = maximum_grid_level_l
             ENDIF

          ELSE
!
!--          Check pre-defined value and reset to default, if neccessary
             IF ( mg_switch_to_pe0_level < mg_switch_to_pe0_level_l  .OR.                          &
                  mg_switch_to_pe0_level >= maximum_grid_level_l )  THEN
                message_string = 'mg_switch_to_pe0_level out of range and reset to 0'
                CALL message( 'init_pegrid', 'PAC0243', 0, 1, 0, 6, 0 )
                mg_switch_to_pe0_level = 0
             ELSE
!
!--             Use the largest number of possible levels anyway and recalculate the switch level to
!--             this largest number of possible values
                maximum_grid_level = maximum_grid_level_l

             ENDIF

          ENDIF

       ENDIF

       ALLOCATE( grid_level_count(maximum_grid_level),                                             &
                 nxl_mg(0:maximum_grid_level), nxr_mg(0:maximum_grid_level),                       &
                 nyn_mg(0:maximum_grid_level), nys_mg(0:maximum_grid_level),                       &
                 nzt_mg(0:maximum_grid_level) )

       grid_level_count = 0
!
!--    Index zero required as dummy due to definition of arrays f2 and p2 in recursive subroutine
!--    next_mg_level
       nxl_mg(0) = 0; nxr_mg(0) = 0; nyn_mg(0) = 0; nys_mg(0) = 0; nzt_mg(0) = 0

       nxl_l = nxl; nxr_l = nxr; nys_l = nys; nyn_l = nyn; nzt_l = nzt

       DO  i = maximum_grid_level, 1 , -1

          IF ( i == mg_switch_to_pe0_level )  THEN
#if defined( __parallel )
!
!--          Save the grid size of the subdomain at the switch level, because it is needed in poismg.
             ind(1) = nxl_l; ind(2) = nxr_l
             ind(3) = nys_l; ind(4) = nyn_l
             ind(5) = nzt_l
             ALLOCATE( ind_all(5*numprocs), mg_loc_ind(5,0:numprocs-1) )
             CALL MPI_ALLGATHER( ind, 5, MPI_INTEGER, ind_all, 5, &
                                 MPI_INTEGER, comm2d, ierr )
             DO  j = 0, numprocs-1
                DO  k = 1, 5
                   mg_loc_ind(k,j) = ind_all(k+j*5)
                ENDDO
             ENDDO
             DEALLOCATE( ind_all )
!
!--          Calculate the grid size of the total domain
             nxr_l = ( nxr_l-nxl_l+1 ) * npex - 1
             nxl_l = 0
             nyn_l = ( nyn_l-nys_l+1 ) * npey - 1
             nys_l = 0
!
!--          The size of this gathered array must not be larger than the array tend, which is used
!--          in the multigrid scheme as a temporary array. Therefore the subdomain size of an PE is
!--          calculated and the size of the gathered grid. These values are used in routines pres
!--          and poismg.
             subdomain_size = ( nxr - nxl + 2 * nbgp + 1 ) *                                       &
                              ( nyn - nys + 2 * nbgp + 1 ) * ( nzt - nzb + 2 )
             gathered_size  = ( nxr_l - nxl_l + 3 ) * ( nyn_l - nys_l + 3 ) * ( nzt_l - nzb + 2 )

#else
             message_string = 'multigrid gather/scatter impossible in non parallel mode'
             CALL message( 'init_pegrid', 'PAC0244', 1, 2, 0, 6, 0 )
#endif
          ENDIF

          nxl_mg(i) = nxl_l
          nxr_mg(i) = nxr_l
          nys_mg(i) = nys_l
          nyn_mg(i) = nyn_l
          nzt_mg(i) = nzt_l

          nxl_l = nxl_l / 2
          nxr_l = nxr_l / 2
          nys_l = nys_l / 2
          nyn_l = nyn_l / 2
          nzt_l = nzt_l / 2

       ENDDO

!
!--    Temporary problem: Currently calculation of maxerror in routine poismg crashes if grid data
!--    are collected on PE0 already on the finest grid level.
!--    To be solved later.
       IF ( maximum_grid_level == mg_switch_to_pe0_level )  THEN
          message_string = 'grid coarsening on subdomain level cannot be performed'
          CALL message( 'poismg', 'PAC0245', 1, 2, 0, 6, 0 )
       ENDIF

    ELSE

       maximum_grid_level = 0

    ENDIF

!
!-- Default level 0 tells exchange_horiz that all ghost planes have to be exchanged. grid_level is
!-- adjusted in poismg, where only one ghost plane is required.
    grid_level = 0

#if defined( __parallel )
!
!-- Gridpoint number for the exchange of ghost points (y-line for 2D-arrays)
    ngp_y  = nyn - nys + 1 + 2 * nbgp

!
!-- Define new MPI derived datatypes for the exchange of ghost points in x- and y-direction for
!-- 2D-arrays (line)
    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp, ngp_y, MPI_REAL, type_x, ierr )
    CALL MPI_TYPE_COMMIT( type_x, ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_y, ngp_y, MPI_REAL, type_y, ierr )
    CALL MPI_TYPE_COMMIT( type_y, ierr )
!
!-- Define new MPI derived datatypes for the exchange of ghost points in x- and y-direction for
!-- 2D-INTEGER arrays (line) - on normal grid.
!-- Define types for 32-bit and 8-bit Integer. The 8-bit Integer are only required on normal grid,
!-- while 32-bit Integer may be also required on coarser grid level in case of multigrid solver.
!
!-- 8-bit Integer
    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp, ngp_y, MPI_BYTE, type_x_byte, ierr )
    CALL MPI_TYPE_COMMIT( type_x_byte, ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_y, ngp_y, MPI_BYTE, type_y_byte, ierr )
    CALL MPI_TYPE_COMMIT( type_y_byte, ierr )
!
!-- 32-bit Integer
    ALLOCATE( type_x_int(0:maximum_grid_level), type_y_int(0:maximum_grid_level) )

    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp, ngp_y, MPI_INTEGER, type_x_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_x_int(0), ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_y, ngp_y, MPI_INTEGER, type_y_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_y_int(0), ierr )
!
!-- Calculate gridpoint numbers for the exchange of ghost points along x (yz-plane for 3D-arrays)
!-- and define MPI derived data type(s) for the exchange of ghost points in y-direction (xz-plane).
!-- Do these calculations for the model grid and (if necessary) also for the coarser grid levels
!-- used in the multigrid method
    ALLOCATE ( ngp_xz(0:maximum_grid_level),                                                       &
               ngp_xz_int(0:maximum_grid_level),                                                   &
               ngp_yz(0:maximum_grid_level),                                                       &
               ngp_yz_int(0:maximum_grid_level),                                                   &
               type_xz(0:maximum_grid_level),                                                      &
               type_xz_int(0:maximum_grid_level),                                                  &
               type_yz(0:maximum_grid_level),                                                      &
               type_yz_int(0:maximum_grid_level) )

    nxl_l = nxl; nxr_l = nxr; nys_l = nys; nyn_l = nyn; nzb_l = nzb; nzt_l = nzt

!
!-- Discern between the model grid, which needs nbgp ghost points and grid levels for the multigrid
!-- scheme. In the latter case only one ghost point is necessary.
!-- First definition of MPI-datatypes for exchange of ghost layers on normal grid. The following
!-- loop is needed for data exchange in poismg.f90.
!
!-- Determine number of grid points of yz-layer for exchange
    ngp_yz(0) = (nzt - nzb + 2) * (nyn - nys + 1 + 2 * nbgp)

!
!-- Define an MPI-datatype for the exchange of left/right boundaries.
!-- Although data are contiguous in physical memory (which does not necessarily require an
!-- MPI-derived datatype), the data exchange between left and right PE's using the MPI-derived type
!-- is 10% faster than without.
    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp*(nzt-nzb+2), ngp_yz(0), MPI_REAL, type_xz(0),     &
                          ierr )
    CALL MPI_TYPE_COMMIT( type_xz(0), ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_yz(0), ngp_yz(0), MPI_REAL, type_yz(0), ierr )
    CALL MPI_TYPE_COMMIT( type_yz(0), ierr )

!
!-- Define data types for exchange of 3D Integer arrays.
    ngp_yz_int(0) = (nzt - nzb + 2) * (nyn - nys + 1 + 2 * nbgp)

    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*nbgp, nbgp*(nzt-nzb+2), ngp_yz_int(0), MPI_INTEGER,          &
                          type_xz_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_xz_int(0), ierr )

    CALL MPI_TYPE_VECTOR( nbgp, ngp_yz_int(0), ngp_yz_int(0), MPI_INTEGER, type_yz_int(0), ierr )
    CALL MPI_TYPE_COMMIT( type_yz_int(0), ierr )

!
!-- Definition of MPI-datatypes for multigrid method (coarser level grids)
    IF ( psolver(1:9) == 'multigrid' )  THEN
!
!--    Definition of MPI-datatyoe as above, but only 1 ghost level is used
       DO  i = maximum_grid_level, 1 , -1
!
!--       For 3D-exchange on different multigrid level, one ghost point for REAL arrays, two ghost
!--       points for INTEGER arrays
          ngp_xz(i) = (nzt_l - nzb_l + 2) * (nxr_l - nxl_l + 3)
          ngp_yz(i) = (nzt_l - nzb_l + 2) * (nyn_l - nys_l + 3)

          ngp_xz_int(i) = (nzt_l - nzb_l + 2) * (nxr_l - nxl_l + 3)
          ngp_yz_int(i) = (nzt_l - nzb_l + 2) * (nyn_l - nys_l + 3)
!
!--       MPI data type for REAL arrays, for xz-layers
          CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+3, nzt_l-nzb_l+2, ngp_yz(i), MPI_REAL, type_xz(i),     &
                                ierr )
          CALL MPI_TYPE_COMMIT( type_xz(i), ierr )

!
!--       MPI data type for INTEGER arrays, for xz-layers
          CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+3, nzt_l-nzb_l+2, ngp_yz_int(i), MPI_INTEGER,          &
                                type_xz_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_xz_int(i), ierr )

!
!--       MPI data type for REAL arrays, for yz-layers
          CALL MPI_TYPE_VECTOR( 1, ngp_yz(i), ngp_yz(i), MPI_REAL, type_yz(i), ierr )
          CALL MPI_TYPE_COMMIT( type_yz(i), ierr )
!
!--       MPI data type for INTEGER arrays, for yz-layers
          CALL MPI_TYPE_VECTOR( 1, ngp_yz_int(i), ngp_yz_int(i), MPI_INTEGER, type_yz_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_yz_int(i), ierr )


!--       For 2D-exchange of INTEGER arrays on coarser grid level, where 2 ghost points need to be
!--       exchanged. Only required for 32-bit Integer arrays.
          CALL MPI_TYPE_VECTOR( nxr_l-nxl_l+5, 2, nyn_l-nys_l+5, MPI_INTEGER, type_x_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_x_int(i), ierr )


          CALL MPI_TYPE_VECTOR( 2, nyn_l-nys_l+5, nyn_l-nys_l+5, MPI_INTEGER, type_y_int(i), ierr )
          CALL MPI_TYPE_COMMIT( type_y_int(i), ierr )

          nxl_l = nxl_l / 2
          nxr_l = nxr_l / 2
          nys_l = nys_l / 2
          nyn_l = nyn_l / 2
          nzt_l = nzt_l / 2

       ENDDO

    ENDIF

#endif

#if defined( __parallel )
!
!-- Setting of flags for inflow/outflow/nesting conditions.
    IF ( pleft == MPI_PROC_NULL )  THEN
       IF ( bc_lr == 'dirichlet/radiation'  .OR.  bc_lr == 'nested'  .OR.                          &
            bc_lr == 'nesting_offline' )  THEN
          bc_dirichlet_l  = .TRUE.
       ELSEIF ( bc_lr == 'radiation/dirichlet' )  THEN
          bc_radiation_l = .TRUE.
       ENDIF
    ENDIF

    IF ( pright == MPI_PROC_NULL )  THEN
       IF ( bc_lr == 'dirichlet/radiation' )  THEN
          bc_radiation_r = .TRUE.
       ELSEIF ( bc_lr == 'radiation/dirichlet'  .OR.  bc_lr == 'nested'  .OR.                      &
                bc_lr == 'nesting_offline' )  THEN
          bc_dirichlet_r  = .TRUE.
       ENDIF
    ENDIF

    IF ( psouth == MPI_PROC_NULL )  THEN
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          bc_radiation_s = .TRUE.
       ELSEIF ( bc_ns == 'radiation/dirichlet'  .OR.  bc_ns == 'nested'  .OR.                      &
                bc_ns == 'nesting_offline' )  THEN
          bc_dirichlet_s  = .TRUE.
       ENDIF
    ENDIF

    IF ( pnorth == MPI_PROC_NULL )  THEN
       IF ( bc_ns == 'dirichlet/radiation'  .OR.  bc_ns == 'nested'  .OR.                          &
            bc_ns == 'nesting_offline' )  THEN
          bc_dirichlet_n  = .TRUE.
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          bc_radiation_n = .TRUE.
       ENDIF
    ENDIF
!
!-- In case of synthetic turbulence geneartor determine ids.
!-- Please note, if no forcing or nesting is applied, the generator is applied only at the left
!-- lateral boundary.
    IF ( syn_turb_gen )  THEN
       IF ( bc_dirichlet_l )  THEN
          id_stg_left_l = myidx
       ELSE
          id_stg_left_l = 0
       ENDIF
       IF ( bc_dirichlet_r )  THEN
          id_stg_right_l = myidx
       ELSE
          id_stg_right_l = 0
       ENDIF
       IF ( bc_dirichlet_s )  THEN
          id_stg_south_l = myidy
       ELSE
          id_stg_south_l = 0
       ENDIF
       IF ( bc_dirichlet_n )  THEN
          id_stg_north_l = myidy
       ELSE
          id_stg_north_l = 0
       ENDIF

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_stg_left_l, id_stg_left,   1, MPI_INTEGER, MPI_SUM, comm1dx, ierr )

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_stg_right_l, id_stg_right, 1, MPI_INTEGER, MPI_SUM, comm1dx, ierr )

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_stg_south_l, id_stg_south, 1, MPI_INTEGER, MPI_SUM, comm1dy, ierr )

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_stg_north_l, id_stg_north, 1, MPI_INTEGER, MPI_SUM, comm1dy, ierr )

    ENDIF

!
!-- Broadcast the id of the outflow PE and outflow-source plane
    IF ( turbulent_outflow )  THEN

       IF ( bc_radiation_r )  THEN
          id_outflow_l = myidx
       ELSE
          id_outflow_l = 0
       ENDIF
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_outflow_l, id_outflow, 1, MPI_INTEGER, MPI_SUM, &
                           comm1dx, ierr )

       IF ( NINT( outflow_source_plane / dx ) >= nxl  .AND.                                        &
            NINT( outflow_source_plane / dx ) <= nxr )  THEN
          id_outflow_source_l = myidx
       ELSE
          id_outflow_source_l = 0
       ENDIF
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( id_outflow_source_l, id_outflow_source, 1, MPI_INTEGER, MPI_SUM,        &
                           comm1dx, ierr )

    ENDIF

    CALL location_message( 'creating virtual PE grids + MPI derived data types', 'finished' )

#else
    IF ( bc_lr == 'dirichlet/radiation' )  THEN
       bc_dirichlet_l = .TRUE.
       bc_radiation_r = .TRUE.
    ELSEIF ( bc_lr == 'radiation/dirichlet' )  THEN
       bc_radiation_l = .TRUE.
       bc_dirichlet_r = .TRUE.
    ENDIF

    IF ( bc_ns == 'dirichlet/radiation' )  THEN
       bc_dirichlet_n = .TRUE.
       bc_radiation_s = .TRUE.
    ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
       bc_radiation_n = .TRUE.
       bc_dirichlet_s = .TRUE.
    ENDIF
#endif

!
!-- At the inflow or outflow, u or v, respectively, have to be calculated for one more grid point.
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       nxlu = nxl + 1
    ELSE
       nxlu = nxl
    ENDIF
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       nysv = nys + 1
    ELSE
       nysv = nys
    ENDIF

!
!-- Allocate wall flag arrays used in the multigrid solver
    IF ( psolver(1:9) == 'multigrid' )  THEN

       DO  i = maximum_grid_level, 1, -1

           SELECT CASE ( i )

              CASE ( 1 )
                 ALLOCATE( topo_flags_1(nzb:nzt_mg(i)+1,                                           &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 2 )
                 ALLOCATE( topo_flags_2(nzb:nzt_mg(i)+1,                                           &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 3 )
                 ALLOCATE( topo_flags_3(nzb:nzt_mg(i)+1,                                           &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 4 )
                 ALLOCATE( topo_flags_4(nzb:nzt_mg(i)+1,                                           &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 5 )
                 ALLOCATE( topo_flags_5(nzb:nzt_mg(i)+1,                                           &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 6 )
                 ALLOCATE( topo_flags_6(nzb:nzt_mg(i)+1,                                           &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 7 )
                 ALLOCATE( topo_flags_7(nzb:nzt_mg(i)+1,                                           &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 8 )
                 ALLOCATE( topo_flags_8(nzb:nzt_mg(i)+1,                                           &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 9 )
                 ALLOCATE( topo_flags_9(nzb:nzt_mg(i)+1,                                           &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE ( 10 )
                 ALLOCATE( topo_flags_10(nzb:nzt_mg(i)+1,                                          &
                                        nys_mg(i)-1:nyn_mg(i)+1,                                   &
                                        nxl_mg(i)-1:nxr_mg(i)+1) )

              CASE DEFAULT
                 message_string = 'more than 10 multigrid levels'
                 CALL message( 'init_pegrid', 'PAC0246', 1, 2, 0, 6, 0 )

          END SELECT

       ENDDO

    ENDIF

#if defined( __parallel )
 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks, if the total domain can be split up to uniform subdomains, so that each PE gets the
!> same number of gridpoints.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE check_subdomain_uniformity

!
!-- Check if current grid and PE combination requires non-uniform subdomains.
    IF ( MOD( nx+1, npex ) == 0  .AND.  MOD( ny+1, npey ) == 0 )  THEN

!
!--    Since PE numbers are integral divisors of total grid point numbers,
!--    the subdomains are uniform.
       non_uniform_subdomain = .FALSE.
!
!--    Data arrays to be transposed have additional requirements concerning uniformity.
!--    the following conditions is true.
       IF ( psolver == 'poisfft_sm' )  THEN
          IF ( MOD( nz, npey ) /= 0  .OR.  MOD( nx+1, npey ) /= 0  .OR.  MOD( ny+1, npex ) /= 0 )  &
          THEN
             non_uniform_data_for_transpose = .TRUE.
          ENDIF
       ELSE
          IF ( MOD( nz, npex ) /= 0  .OR.  MOD( nx+1, npey ) /= 0  .OR.  MOD( ny+1, npex ) /= 0 )  &
          THEN
             non_uniform_data_for_transpose = .TRUE.
          ENDIF
       ENDIF

    ELSE

       non_uniform_subdomain = .TRUE.
       non_uniform_data_for_transpose = .TRUE.

    ENDIF

 END SUBROUTINE check_subdomain_uniformity


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates indices and gridpoint numbers required for transposition of 3d-arrays in case
!> of a domain decomposition which has non-uniform subdomains.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE setup_transpose_indices

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  irest_xy
    INTEGER(iwp) ::  irest_yd
    INTEGER(iwp) ::  irest_yz
    INTEGER(iwp) ::  irest_z
    INTEGER(iwp) ::  j
    INTEGER(iwp) ::  k
    INTEGER(iwp) ::  nnz_yd


    ALLOCATE( nxl_y_pe(0:npey-1) )
    ALLOCATE( nxr_y_pe(0:npey-1) )
    ALLOCATE( nyn_z_pe(0:npex-1) )
    ALLOCATE( nys_z_pe(0:npex-1) )
    ALLOCATE( nzb_x_pe(0:npex-1) )
    ALLOCATE( nzt_x_pe(0:npex-1) )

!
!-- 1. transposition  z --> x
    nys_x   = nys
    nyn_x   = nyn
    nny_x   = nny
    nnz_x   = nz / npex
    irest_z = MOD( nz, npex )

    DO  k = 0, npex-1
       IF ( k < irest_z )  THEN
          nzb_x_pe(k) = k * ( nnz_x + 1 ) + 1
          nzt_x_pe(k) = nzb_x_pe(k) + ( nnz_x + 1 ) - 1
       ELSE
          nzb_x_pe(k) = irest_z * ( nnz_x + 1 ) + ( k - irest_z ) * nnz_x + 1
          nzt_x_pe(k) = nzb_x_pe(k) + nnz_x - 1
       ENDIF
    ENDDO

    nzb_x = nzb_x_pe(myidx)
    nzt_x = nzt_x_pe(myidx)
    nnz_x = nzt_x - nzb_x + 1
!
!-- Calculate index bounds and number of grid points that are required for allocating the
!-- 3d-arrays used by the transpose routine. Arrays on all PEs need to have identical size.
!-- Therefore, the dimensions defined on PE0 are used, because in case of non-uniform
!-- subdomains, at least on PE0 the respective arrays contain one more grid point than on the
!-- other PEs.
    nxr_x_max = nxl + nnx_pe(0) - 1
    nnx_x_max = nnx_pe(0)
    nnz_x_max = nzt_x_pe(0) - nzb_x_pe(0) + 1
    nz_x_max  = npex * nnz_x_max
    nzt_x_max = nzb_x + nnz_x_max - 1
    sendrecvcount_zx = nnx_x_max * nny * nnz_x_max

!
!-- 2. transposition  x --> y
    nnz_y    = nnz_x
    nzb_y    = nzb_x
    nzt_y    = nzt_x
    nnx_y    = (nx+1) / npey
    irest_xy = MOD( (nx+1), npey )

    DO  j = 0, npey-1
       IF ( j < irest_xy )  THEN
          nxl_y_pe(j) = j * ( nnx_y + 1 )
          nxr_y_pe(j) = nxl_y_pe(j) + ( nnx_y + 1 ) - 1
       ELSE
          nxl_y_pe(j) = irest_xy * ( nnx_y + 1 ) + ( j - irest_xy ) * nnx_y
          nxr_y_pe(j) = nxl_y_pe(j) + nnx_y - 1
       ENDIF
    ENDDO

    nxl_y = nxl_y_pe(myidy)
    nxr_y = nxr_y_pe(myidy)
    nnx_y = nxr_y - nxl_y + 1
!
!-- Calculate index bounds and number of grid points that are required for allocating the
!-- 3d-arrays used by the transpose routine. Arrays on all PEs need to have identical size.
!-- Therefore, the dimensions defined on PE0 are used, because in case of non-uniform
!-- subdomains, at least on PE0 the respective arrays contain one more grid point than on the
!-- other PEs.
    nyn_x_max = nys_x + nny_pe(0) - 1
    nx_y_max  = npey * ( nxr_y_pe(0) - nxl_y_pe(0) + 1 ) - 1
    nnx_y_max = nxr_y_pe(0) - nxl_y_pe(0) + 1
    nxr_y_max = nxl_y + nxr_y_pe(0) - nxl_y_pe(0)
    sendrecvcount_xy = ( nyn_x_max-nys_x + 1 ) * ( nzt_y-nzb_y + 1 ) * ( nxr_y_max-nxl_y + 1 )

!
!-- 3. transposition  y --> z
!-- (ELSE:  x --> y  in case of 1D-decomposition along x)
    nxl_z    = nxl_y
    nxr_z    = nxr_y
    nny_z    = ( ny+1 ) / npex
    irest_yz = MOD( (ny+1), npex )

    DO  i = 0, npex-1
       IF ( i < irest_yz )  THEN
          nys_z_pe(i) = i * ( nny_z + 1 )
          nyn_z_pe(i) = nys_z_pe(i) + ( nny_z + 1 ) - 1
       ELSE
          nys_z_pe(i) = irest_yz * ( nny_z + 1 ) + ( i - irest_yz ) * nny_z
          nyn_z_pe(i) = nys_z_pe(i) + nny_z - 1
       ENDIF
    ENDDO

    nys_z = nys_z_pe(myidx)
    nyn_z = nyn_z_pe(myidx)

!
!-- Calculate index bounds and number of grid points that are required for allocating the
!-- 3d-arrays used by the transpose routine. Arrays on all PEs need to have identical size.
!-- Therefore, the dimensions defined on PE0 are used, because in case of non-uniform
!-- subdomains, at least on PE0 the respective arrays contain one more grid point than on the
!-- other PEs.
    nzt_y_max = nzt_x_max
    ny_z_max  = npex * ( nyn_z_pe(0) - nys_z_pe(0) + 1 ) - 1
    nnz_z_max = ( nz - 1 ) / npex + 1
    nny_z_max = nyn_z_pe(0) - nys_z_pe(0) + 1
    nyn_z_max = nys_z + nyn_z_pe(0) - nys_z_pe(0)
    sendrecvcount_yz = ( nxr_z-nxl_z + 1 ) * nnz_z_max * ( nyn_z_max-nys_z + 1 )

!
!-- Indices for direct transpositions z --> y (used for calculating spectra)
    ALLOCATE( nzb_yd_pe(0:npey-1) )
    ALLOCATE( nzt_yd_pe(0:npey-1) )

    nxl_yd   = nxl
    nxr_yd   = nxr
    nnz_yd   = nz / npey
    irest_yd = MOD( nz, npey )

    DO  k = 0, npey-1
       IF ( k < irest_yd )  THEN
          nzb_yd_pe(k) = k * ( nnz_yd + 1 ) + 1
          nzt_yd_pe(k) = nzb_yd_pe(k) + ( nnz_yd + 1 ) - 1
       ELSE
          nzb_yd_pe(k) = irest_yd * ( nnz_yd + 1 ) + ( k - irest_yd ) * nnz_yd + 1
          nzt_yd_pe(k) = nzb_yd_pe(k) + nnz_yd - 1
       ENDIF
    ENDDO

    nzb_yd = nzb_yd_pe(myidy)
    nzt_yd = nzt_yd_pe(myidy)
    nnz_yd = nzt_yd - nzb_yd + 1

    nyn_yd_max = nys + nny_pe(0) - 1
    nnz_yd_max = nzt_yd_pe(0) - nzb_yd_pe(0) + 1
    nz_yd_max  = npey * nnz_yd_max
    nny_yd_max = nny_pe(0)
    nzt_yd_max = nzb_yd + nnz_yd_max - 1
    sendrecvcount_zyd = nnx * nnz_yd_max * nny_yd_max

 END SUBROUTINE setup_transpose_indices
#endif

 END SUBROUTINE init_pegrid
