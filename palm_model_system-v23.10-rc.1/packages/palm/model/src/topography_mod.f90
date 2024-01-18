!> @file topography_mod.f90
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
!> Setup of PALM's topography representation
!> @todo: Rearrange topo flag list
!> @todo: reference 3D buildings on top of orography is not tested and may need further improvement
!>        for steep slopes
!> @todo: Use more advanced setting of building type at filled holes
!--------------------------------------------------------------------------------------------------!
 MODULE topography_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  dzw,                                                                                &
               dzu,                                                                                &
               x,                                                                                  &
               y,                                                                                  &
               zu,                                                                                 &
               zw

    USE boundary_settings_mod,                                                                     &
        ONLY:  set_lateral_neumann_bc

    USE control_parameters

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point_s

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz_2d,                                                                  &
               exchange_horiz_2d_byte,                                                             &
               exchange_horiz_2d_int,                                                              &
               exchange_horiz_int

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nys,                                                                                &
               nysg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzb_max,                                                                            &
               nzt,                                                                                &
               topo_min_level,                                                                     &
               topo_top_ind,                                                                       &
               topo_flags

    USE general_utilities,                                                                         &
        ONLY:  gridpoint_id

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy,                                                                                 &
               zu_s_inner,                                                                         &
               zw_w_inner

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  add_ghost_layers,                                                                   &
               buildings_f,                                                                        &
               building_id_f,                                                                      &
               building_type_f,                                                                    &
               char_fill,                                                                          &
               char_lod,                                                                           &
               check_existence,                                                                    &
               close_input_file,                                                                   &
               dims_xy,                                                                            &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               init_model,                                                                         &
               input_file_static,                                                                  &
               input_pids_static,                                                                  &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               open_read_file,                                                                     &
               terrain_height_f

    USE pegrid

#if defined( __parallel )
    USE pmc_handle_communicator,                                                                   &
        ONLY:  pmc_get_model_info
#endif

    USE pmc_interface,                                                                             &
        ONLY:  atmosphere_ocean_coupled_run,                                                       &
               nest_shift_z

    IMPLICIT NONE


    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  topo !< temporary array used to initialize topography

    LOGICAL ::  topo_read_all_domains !< control flag indicating whether topography is read from file in all domains

    SAVE

    PRIVATE
!
!-- Public subroutines
    PUBLIC init_topography

    INTERFACE init_topography
       MODULE PROCEDURE init_topography
    END INTERFACE init_topography

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Routine that controls the basic topography initialization.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_topography

    CALL location_message( 'Setup topography', 'start' )
!
!-- Allocate 3D array to set topography
    ALLOCATE( topo(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    topo = 0
!
!-- Initialize topography by generic topography or read topography from file.
!-- Further process the topography, e.g. filter small cavities on the grid scale, if required, map
!-- buildings onto the underlying terrain, etc..
    CALL define_topography
!
!-- Set flags to mark the topography features on the grid.
    CALL topography_set_flags
!
!-- Determine further topography grid indices, e.g. to output profile data, or to control
!-- the degradation of the advection terms.
    CALL topography_set_indices

    CALL location_message( 'Setup topography', 'finished' )

 END SUBROUTINE init_topography


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Reads topography information from file or sets generic topography. Moreover, all
!> topography-relevant topography arrays are initialized, and grid flags are set.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE define_topography

    INTEGER(iwp) ::  bh                !< temporary vertical index of building height
    INTEGER(iwp) ::  ch                !< temporary vertical index for canyon height
    INTEGER(iwp) ::  hv_in             !< heavyside function to model inner tunnel surface
    INTEGER(iwp) ::  i                 !< index variable along x
    INTEGER(iwp) ::  index_left_bwall  !< index for left building wall
    INTEGER(iwp) ::  index_north_bwall !< index for north building wall
    INTEGER(iwp) ::  index_right_bwall !< index for right building wall
    INTEGER(iwp) ::  index_south_bwall !< index for south building wall
    INTEGER(iwp) ::  index_left_cwall  !< index for left canyon wall
    INTEGER(iwp) ::  index_north_cwall !< index for north canyon wall
    INTEGER(iwp) ::  index_right_cwall !< index for right canyon wall
    INTEGER(iwp) ::  index_south_cwall !< index for south canyon wall
    INTEGER(iwp) ::  j                 !< index variable along y
    INTEGER(iwp) ::  k                 !< index variable along z
    INTEGER(iwp) ::  ngp_bx            !< grid point number of building size along x
    INTEGER(iwp) ::  ngp_by            !< grid point number of building size along y
    INTEGER(iwp) ::  ngp_cx            !< grid point number of canyon size along x
    INTEGER(iwp) ::  ngp_cy            !< grid point number of canyon size along y
    INTEGER(iwp) ::  hv_out            !< heavyside function to model outer tunnel surface
    INTEGER(iwp) ::  td                !< tunnel wall depth
    INTEGER(iwp) ::  th                !< height of outer tunnel wall
    INTEGER(iwp) ::  txe_in            !< end position of inner tunnel wall in x
    INTEGER(iwp) ::  txe_out           !< end position of outer tunnel wall in x
    INTEGER(iwp) ::  txs_in            !< start position of inner tunnel wall in x
    INTEGER(iwp) ::  txs_out           !< start position of outer tunnel wall in x
    INTEGER(iwp) ::  tye_in            !< end position of inner tunnel wall in y
    INTEGER(iwp) ::  tye_out           !< end position of outer tunnel wall in y
    INTEGER(iwp) ::  tys_in            !< start position of inner tunnel wall in y
    INTEGER(iwp) ::  tys_out           !< start position of outer tunnel wall in y

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  nzb_local !< index for topography top at cell-center

    LOGICAL ::  root_model  !< flag to indicate if root or child


!
!-- If a topography input file is available (static input file or ASCII), read the topography
!-- information. Note, this is already done here in order to check for correct parameter settings
!-- and file availability.
    CALL topography_input
!
!-- Check for correct setting of the namelist parameter topography. If topography information is
!-- read from file but topography = 'flat', initialization does not work properly.
    IF ( ( buildings_f%from_file  .OR.  terrain_height_f%from_file )  .AND.                        &
           TRIM( topography ) /= 'read_from_file' )                                                &
    THEN
       message_string = 'wrong setting topography = "' // TRIM( topography ) // '"'
       CALL message( 'topography_mod', 'PAC0319', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check, if 'read_from_file' has been chosen for all domains (parents/childs) and set a
!-- respective flag. Only then lateron global reduction operations (e.g. for getting the lowest
!-- topography throughout all domains) will be allowed.
    topo_read_all_domains = ( topography == 'read_from_file' )
#if defined( __parallel )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, topo_read_all_domains, 1, MPI_LOGICAL, MPI_LAND,             &
                        MPI_COMM_WORLD, ierr )
#endif
!
!-- Define topography and set respective flags. Topography is either flat (only grid point at
!-- k=nzb is flagged), generically defined (block, canyon or tunnel), or it can be read from file.
!-- In case of elevated childs, the lowest grid point at k=nzb may belong to the atmosphere, so
!-- this point is not flagged.
    SELECT CASE ( TRIM( topography ) )

       CASE ( 'flat' )
!
!--       Initialize 3D topography array, used later for initializing flags.
!--       In case of vertically shifted nests, the boundary vertical coordinate does not have its
!--       standard location 0.0, meaning that the respective boundary is "open".
!--       Childs in ocean_mode always have an open bottom boundary.
          root_model = .TRUE.
#if defined( __parallel )
          IF ( nested_run )  CALL pmc_get_model_info( root_model = root_model )
#endif
          IF ( nest_shift_z == 0.0_wp  .AND.  .NOT. ( ocean_mode  .AND.  .NOT. root_model ) )  THEN
             topo(nzb,:,:) = IBSET( topo(nzb,:,:), 0 )
          ENDIF

       CASE ( 'closed_channel' )
!
!--       Initialilize 3D topography array, used later for initializing flags.
          topo(nzb,:,:) = IBSET( topo(nzb,:,:), 0 )

       CASE ( 'single_building' )
!
!--       Single rectangular building, by default centered in the middle of the
!--       total domain.
          ngp_bx = NINT( building_length_x / dx )
          ngp_by = NINT( building_length_y / dy )
          bh  = MINLOC( ABS( zw - building_height ), 1 ) - 1
          IF ( ABS( zw(bh) - building_height ) ==  ABS( zw(bh+1) - building_height ) )  bh = bh + 1
          IF ( building_wall_left == 9999999.9_wp )  THEN
             building_wall_left = ( nx + 1 - ngp_bx ) / 2 * dx
          ENDIF
          index_left_bwall  = NINT( building_wall_left / dx )
          index_right_bwall = index_left_bwall + ngp_bx

          IF ( building_wall_south == 9999999.9_wp )  THEN
              building_wall_south = ( ny + 1 - ngp_by ) / 2 * dy
          ENDIF
          index_south_bwall = NINT( building_wall_south / dy )
          index_north_bwall = index_south_bwall + ngp_by

!
!--       Building size has to meet some requirements.
          IF ( ( index_left_bwall  < 1 )  .OR.  ( index_right_bwall > nx-1 )  .OR.                 &
               ( index_right_bwall < index_left_bwall+3 )  .OR.                                    &
               ( index_south_bwall < 1 )  .OR.  ( index_north_bwall > ny-1 )  .OR.                 &
               ( index_north_bwall < index_south_bwall+3 ) )                                       &
          THEN
             WRITE( message_string, * ) 'inconsistent building parameters:',                       &
                                        '&index_left_bwall=', index_left_bwall,                    &
                                        'index_right_bwall=', index_right_bwall,                   &
                                        'index_south_bwall=', index_south_bwall,                   &
                                        'index_north_bwall=', index_north_bwall,                   &
                                        'nx=', nx, 'ny=', ny
             CALL message( 'topography_mod', 'PAC0320', 1, 2, 0, 6, 0 )
          ENDIF

          ALLOCATE( nzb_local(nysg:nyng,nxlg:nxrg) )
          nzb_local = 0
!
!--       Define the building.
          IF ( index_left_bwall <= nxr  .AND.  index_right_bwall >= nxl  .AND.                     &
               index_south_bwall <= nyn  .AND.  index_north_bwall >= nys )                         &
          THEN
             nzb_local(MAX(nys,index_south_bwall):MIN(nyn,index_north_bwall),                      &
                       MAX(nxl,index_left_bwall):MIN(nxr,index_right_bwall)) = bh
          ENDIF
!
!--       Set bit array on basis of nzb_local.
          DO  i = nxl, nxr
             DO  j = nys, nyn
                topo(nzb:nzb_local(j,i),j,i) = IBSET( topo(nzb:nzb_local(j,i),j,i), 0 )
             ENDDO
          ENDDO

          DEALLOCATE( nzb_local )
!
!--       Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
!--       are set for the topography, i.e. it is assumed that buildings continue with same height
!--       outside the gloabl domain.
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
          CALL topography_set_non_cyc_bc( topo )

       CASE ( 'single_street_canyon' )
!
!--       Single quasi-2D street canyon of infinite length in x or y direction.
!--       The canyon is centered in the other direction by default.
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
!
!--          Street canyon in y direction
             ngp_cx = NINT( canyon_width_x / dx )
             IF ( canyon_wall_left == 9999999.9_wp )  THEN
                canyon_wall_left = ( nx + 1 - ngp_cx ) / 2 * dx
             ENDIF
             index_left_cwall= NINT( canyon_wall_left / dx )
             index_right_cwall= index_left_cwall+ ngp_cx
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
!
!--          Street canyon in x direction
             ngp_cy = NINT( canyon_width_y / dy )
             IF ( canyon_wall_south == 9999999.9_wp )  THEN
                canyon_wall_south = ( ny + 1 - ngp_cy ) / 2 * dy
             ENDIF
             index_south_cwall = NINT( canyon_wall_south / dy )
             index_north_cwall = index_south_cwall + ngp_cy

          ELSE

             message_string = 'no street canyon width given'
             CALL message( 'topography_mod', 'PAC0321', 1, 2, 0, 6, 0 )

          ENDIF

          ch  = MINLOC( ABS( zw - canyon_height ), 1 ) - 1
          IF ( ABS( zw(ch) - canyon_height ) == ABS( zw(ch+1) - canyon_height ) )  ch = ch + 1
          dp_level_ind_b = ch
!
!--       Street canyon size has to meet some requirements.
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
             IF ( ( index_left_cwall< 1 ) .OR. ( index_right_cwall> nx-1 ) .OR. ( ngp_cx < 3 ) )   &
             THEN
                WRITE( message_string, * ) 'inconsistent canyon parameters:',                      &
                                           '&index_left_cwall=', index_left_cwall,                 &
                                           ' index_right_cwall=', index_right_cwall,               &
                                           ' ngp_cx=', ngp_cx, ' ch=', ch, ' nx=', nx, ' ny=', ny
                CALL message( 'topography_mod', 'PAC0322', 1, 2, 0, 6, 0 )
             ENDIF
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
             IF ( ( index_south_cwall < 1 ) .OR. ( index_north_cwall > ny-1 ) .OR. ( ngp_cy < 3 ) )&
             THEN
                WRITE( message_string, * ) 'inconsistent canyon parameters:',                      &
                                           '&index_south_cwall=', index_south_cwall,               &
                                           ' index_north_cwall=', index_north_cwall,               &
                                           ' ngp_cy=', ngp_cy, ' ch=', ch, ' nx=', nx, ' ny=', ny
                CALL message( 'topography_mod', 'PAC0323', 1, 2, 0, 6, 0 )
             ENDIF
          ENDIF
          IF ( canyon_width_x /= 9999999.9_wp  .AND.  canyon_width_y /= 9999999.9_wp )  THEN
             message_string = 'inconsistent canyon parameters:' //                                 &
                              '&street canyon can only be oriented either in x- or in y-direction'
             CALL message( 'topography_mod', 'PAC0324', 1, 2, 0, 6, 0 )
          ENDIF

          ALLOCATE( nzb_local(nysg:nyng,nxlg:nxrg) )
          nzb_local = ch
          IF ( canyon_width_x /= 9999999.9_wp )  THEN
             IF ( index_left_cwall<= nxr  .AND.  index_right_cwall>= nxl )                         &
                nzb_local(:,MAX( nxl, index_left_cwall+1 ):MIN( nxr, index_right_cwall-1) ) = 0
          ELSEIF ( canyon_width_y /= 9999999.9_wp )  THEN
             IF ( index_south_cwall <= nyn  .AND.  index_north_cwall >= nys )                      &
                nzb_local(MAX( nys, index_south_cwall+1 ):MIN( nyn, index_north_cwall-1 ),:) = 0
          ENDIF
!
!--       Set bit array on basis of nzb_local
          DO  i = nxl, nxr
             DO  j = nys, nyn
                topo(nzb:nzb_local(j,i),j,i) = IBSET( topo(nzb:nzb_local(j,i),j,i), 0 )
             ENDDO
          ENDDO
          DEALLOCATE( nzb_local )
!
!--       Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
!--       are set for the topography, i.e. it is assumed that buildings continue with same height
!--       outside the gloabl domain.
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
          CALL topography_set_non_cyc_bc( topo )

       CASE ( 'tunnel' )
!
!--       Initialize surface height with zero and set lowest model grid point to topography.
          topo(nzb,:,:)  = IBSET( topo(nzb,:,:), 0 )
!
!--       Tunnel height.
          IF ( tunnel_height == 9999999.9_wp )  THEN
             th = zw( INT( 0.2 * nz) )
          ELSE
             th = tunnel_height
          ENDIF
!
!--       Tunnel-wall depth.
          IF ( tunnel_wall_depth == 9999999.9_wp )  THEN
             td = MAX( dx, dy, dz(1) )
          ELSE
             td = tunnel_wall_depth
          ENDIF
!
!--       Check for tunnel width
          IF ( tunnel_width_x == 9999999.9_wp  .AND.  tunnel_width_y == 9999999.9_wp  )  THEN
             message_string = 'No tunnel width is given'
             CALL message( 'topography_mod', 'PAC0325', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( tunnel_width_x /= 9999999.9_wp  .AND.  tunnel_width_y /= 9999999.9_wp  )  THEN
             message_string = 'inconsistent tunnel parameters tunnel_width_x / tunnel_width_y'
             CALL message( 'topography_mod', 'PAC0326', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Check for too small tunnel width in x- and y-direction
          IF ( tunnel_width_x /= 9999999.9_wp  .AND.                                               &
               tunnel_width_x - 2.0_wp * td <= 2.0_wp * dx )  THEN
             message_string = 'tunnel_width_x too small'
             CALL message( 'topography_mod', 'PAC0327', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( tunnel_width_y /= 9999999.9_wp  .AND.                                               &
               tunnel_width_y - 2.0_wp * td <= 2.0_wp * dy )  THEN
             message_string = 'tunnel_width_y too small'
             CALL message( 'topography_mod', 'PAC0328', 1, 2, 0, 6, 0 )
          ENDIF
!
!--       Check for too large tunnel width. Tunnel axis along y.
          IF ( tunnel_width_x /= 9999999.9_wp )  THEN
             IF ( tunnel_width_x > ( nx + 1 ) * dx )  THEN
                message_string = 'tunnel_width_x too large'
                CALL message( 'topography_mod', 'PAC0329', 1, 2, 0, 6, 0 )
             ENDIF

             txs_out = INT( ( nx + 1 ) * 0.5_wp * dx - tunnel_width_x * 0.5_wp )
             txe_out = INT( ( nx + 1 ) * 0.5_wp * dx + tunnel_width_x * 0.5_wp )
             txs_in  = INT( ( nx + 1 ) * 0.5_wp * dx - ( tunnel_width_x * 0.5_wp - td ) )
             txe_in  = INT( ( nx + 1 ) * 0.5_wp * dx + ( tunnel_width_x * 0.5_wp - td ) )

             tys_out = INT( ( ny + 1 ) * 0.5_wp * dy - tunnel_length * 0.5_wp )
             tye_out = INT( ( ny + 1 ) * 0.5_wp * dy + tunnel_length * 0.5_wp )
             tys_in  = tys_out
             tye_in  = tye_out
          ENDIF
!
!--       Tunnel axis along x.
          IF ( tunnel_width_y /= 9999999.9_wp )  THEN
             IF ( tunnel_width_y > ( ny + 1 ) * dy )  THEN
                message_string = 'tunnel_width_y too large'
                CALL message( 'topography_mod', 'PAC0330', 1, 2, 0, 6, 0 )
             ENDIF

             txs_out = INT( ( nx + 1 ) * 0.5_wp * dx - tunnel_length * 0.5_wp )
             txe_out = INT( ( nx + 1 ) * 0.5_wp * dx + tunnel_length * 0.5_wp )
             txs_in  = txs_out
             txe_in  = txe_out

             tys_out = INT( ( ny + 1 ) * 0.5_wp * dy - tunnel_width_y * 0.5_wp )
             tye_out = INT( ( ny + 1 ) * 0.5_wp * dy + tunnel_width_y * 0.5_wp )
             tys_in  = INT( ( ny + 1 ) * 0.5_wp * dy - ( tunnel_width_y * 0.5_wp - td ) )
             tye_in  = INT( ( ny + 1 ) * 0.5_wp * dy + ( tunnel_width_y * 0.5_wp - td ) )
          ENDIF

          DO  i = nxl, nxr
             DO  j = nys, nyn
!
!--             Use heaviside function to model outer tunnel surface.
                hv_out = th * 0.5_wp * ( ( SIGN( 1.0_wp, i * dx - txs_out ) + 1.0_wp )             &
                                       - ( SIGN( 1.0_wp, i * dx - txe_out ) + 1.0_wp ) )

                hv_out = hv_out * 0.5_wp * ( ( SIGN( 1.0_wp, j * dy - tys_out ) + 1.0_wp )         &
                                           - ( SIGN( 1.0_wp, j * dy - tye_out ) + 1.0_wp ) )
!
!--             Use heaviside function to model inner tunnel surface.
                hv_in  = ( th - td ) * 0.5_wp * ( ( SIGN( 1.0_wp, i * dx - txs_in ) + 1.0_wp )     &
                                                - ( SIGN( 1.0_wp, i * dx - txe_in ) + 1.0_wp ) )

                hv_in = hv_in * 0.5_wp * ( ( SIGN( 1.0_wp, j * dy - tys_in ) + 1.0_wp )            &
                                         - ( SIGN( 1.0_wp, j * dy - tye_in ) + 1.0_wp ) )

                IF ( hv_out - hv_in == 0.0_wp )  THEN
!
!--                Set flags at x-y-positions without any tunnel surface.
                   topo(nzb+1:nzt+1,j,i) = IBCLR( topo(nzb+1:nzt+1,j,i), 0 )
                ELSE
!
!--                Set flags at x-y-positions with tunnel surfaces.
                   DO  k = nzb + 1, nzt + 1
!
!--                   Inner tunnel.
                      IF ( hv_out - hv_in == th )  THEN
                         IF ( zw(k) <= hv_out )  THEN
                            topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         ELSE
                            topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                         ENDIF
                      ENDIF
!
!--                   Lateral tunnel walls
                      IF ( hv_out - hv_in == td )  THEN
                         IF ( zw(k) <= hv_in )  THEN
                            topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                         ELSEIF ( zw(k) > hv_in  .AND.  zw(k) <= hv_out )  THEN
                            topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         ELSEIF ( zw(k) > hv_out )  THEN
                            topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
!
!--       Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
!--       are set for the topography, i.e. it is assumed that buildings continue with same height
!--       outside the gloabl domain.
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
          CALL topography_set_non_cyc_bc( topo )

       CASE ( 'read_from_file' )
!
!--       Note, topography information has already been read.
!--       If required, further process topography, i.e. reference buildings on top of orography and
!--       set temporary 3D topography array, which is used later to set grid flags. Calling of this
!--       routine is also required in case of ASCII input, even though no distinction between
!--       terrain- and building height is made in this case.
          CALL process_topography
!
!--       Filter holes resolved by only one grid-point.
          CALL filter_topography
!
!--       Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
!--       are set for the topography, i.e. it is assumed that buildings continue with same height
!--       outside the gloabl domain.
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
          CALL topography_set_non_cyc_bc( topo )

       CASE DEFAULT
!
!--       The DEFAULT case is reached either if the parameter topography contains a wrong character
!--       string or if the user has defined a special case in the user interface. There, the
!--       subroutine user_init_grid checks which of these two conditions applies.
          CALL user_init_grid( topo )
          CALL filter_topography

    END SELECT

!
!-- Consistency checks and index array initialization are only required for non-flat topography.
    IF ( TRIM( topography ) /= 'flat' )  THEN
!
!--    In case of non-flat topography, check whether the convention how to define the topography
!--    grid has been set correctly, or whether the default is applicable. If this is not possible,
!--    abort.
       IF ( TRIM( topography_grid_convention ) == ' ' )  THEN
          IF ( TRIM( topography ) /= 'closed_channel'        .AND.                                 &
               TRIM( topography ) /= 'single_building'       .AND.                                 &
               TRIM( topography ) /= 'single_street_canyon'  .AND.                                 &
               TRIM( topography ) /= 'tunnel'                .AND.                                 &
               TRIM( topography ) /= 'read_from_file')                                             &
          THEN
!
!--          The default value is not applicable here, because it is only valid for the four
!--          standard cases 'single_building', 'single_street_canyon', 'tunnel' and 'read_from_file'.
             message_string = 'missing value for topography_grid_convention'
             CALL message( 'topography_mod', 'PAC0331', 1, 2, 0, 6, 0 )
          ELSE
!
!--          The default value is applicable here. Set convention according to topography.
             IF ( TRIM( topography ) == 'single_building'  .OR.                                    &
                  TRIM( topography ) == 'single_street_canyon' )                                   &
             THEN
                topography_grid_convention = 'cell_edge'
             ELSEIF ( TRIM( topography ) == 'read_from_file'  .OR.  TRIM( topography ) == 'tunnel')&
             THEN
                topography_grid_convention = 'cell_center'
             ENDIF
          ENDIF
       ELSEIF ( TRIM( topography_grid_convention ) /= 'cell_edge'  .AND.                           &
                TRIM( topography_grid_convention ) /= 'cell_center' )                              &
       THEN
          WRITE( message_string, * )  'illegal value for topography_grid_convention: "' //         &
                                      TRIM( topography_grid_convention ) // '"'
          CALL message( 'topography_mod', 'PAC0332', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( topography_grid_convention == 'cell_edge' )  THEN
!
!--       The topography as defined using the 'cell_edge' convention describes the actual total size
!--       of topography which is defined at the cell edges where u=0 on the topography walls in
!--       x-direction and v=0 on the topography walls in y-direction.
!--       Therefore, the existence of topography in the grid center is now reduced b1dx at the east
!--       (left) topography walls and by 1dy at the north topography walls to form the basis for
!--       the grid-center flag.
!--       Note, the reverse memory access (first j loop, then i loop) is absolutely required at
!--       this point.
          DO  j = nys+1, nyn+1
             DO  i = nxl-1, nxr
                DO  k = nzb, nzt+1
                   IF ( .NOT. BTEST( topo(k,j,i), 0 )  .OR.  .NOT. BTEST( topo(k,j,i+1), 0 ) )     &
                      topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                ENDDO
             ENDDO
          ENDDO
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )

          DO  i = nxl, nxr+1
             DO  j = nys-1, nyn
                DO  k = nzb, nzt+1
                   IF ( .NOT. BTEST( topo(k,j,i), 0 )  .OR.  .NOT. BTEST( topo(k,j+1,i), 0 ) )     &
                      topo(k,j,i) = IBCLR( topo(k,j,i), 0 )
                ENDDO
             ENDDO
          ENDDO
          CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )

       ENDIF
    ENDIF

 END SUBROUTINE define_topography


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Reads orography and building information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE topography_input

    CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE ::  var_names  !< variable names in static input file

    INTEGER(iwp) ::  i             !< running index along x-direction
    INTEGER(iwp) ::  id_topo       !< NetCDF id of topograhy input file
    INTEGER(iwp) ::  ii            !< running index for IO blocks
    INTEGER(iwp) ::  io_status     !< status after reading the ascii topo file
    INTEGER(iwp) ::  j             !< running index along y-direction
    INTEGER(iwp) ::  k             !< running index along z-direction
    INTEGER(iwp) ::  num_vars      !< number of variables in netcdf input file
    INTEGER(iwp) ::  skip_n_rows   !< counting variable to skip rows while reading topography file

    REAL(wp) ::  dum           !< dummy variable to skip columns while reading topography file

    TYPE(dims_xy)   ::  dim_static     !< data structure for x, y-dimension in static input file

!
!-- CPU measurement.
    CALL cpu_log( log_point_s(83), 'NetCDF/ASCII input topo', 'start' )
!
!-- Input via palm-input data standard.
    IF ( input_pids_static )  THEN
#if defined ( __netcdf )
!
!--    Open file in read-only mode.
       CALL open_read_file( TRIM( input_file_static ) // TRIM( coupling_char ), id_topo )
!
!--    At first, inquire all variable names.
!--    This will be used to check whether an input variable exists or not.
       CALL inquire_num_variables( id_topo, num_vars )
!
!--    Allocate memory to store variable names and inquire them.
       ALLOCATE( var_names(1:num_vars) )
       CALL inquire_variable_names( id_topo, var_names )
!
!--    Read x, y - dimensions. Only required for consistency checks.
       CALL get_dimension_length( id_topo, dim_static%nx, 'x' )
       CALL get_dimension_length( id_topo, dim_static%ny, 'y' )
       ALLOCATE( dim_static%x(0:dim_static%nx-1) )
       ALLOCATE( dim_static%y(0:dim_static%ny-1) )
       CALL get_variable( id_topo, 'x', dim_static%x )
       CALL get_variable( id_topo, 'y', dim_static%y )
!
!--    Check whether dimension size in input file matches the model dimensions.
       IF ( dim_static%nx-1 /= nx )  THEN
          WRITE( message_string, * )  'static driver: horizontal dimension in x-direction (=',     &
                                      dim_static%nx-1, ') does not match the respective model ',   &
                                      'dimension (=', nx, ')'
          CALL message( 'topography_mod', 'PAC0333', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( dim_static%ny-1 /= ny )  THEN
          WRITE( message_string, * )  'static driver: horizontal dimension in y-direction (=',     &
                                      dim_static%ny-1, ') does not match the respective model ',   &
                                      'dimension (=', ny, ')'
          CALL message( 'topography_mod', 'PAC0333', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Check if grid spacing of provided input data matches the respective grid spacing in the
!--    model. The allowed tolerance is 0.1% of the respective model grid spacing.
       IF ( ABS( dim_static%x(1) - dim_static%x(0) - dx ) > 0.001_wp * dx )  THEN
          WRITE( message_string, * )  'static driver: horizontal grid spacing in x-direction (=',  &
                                      dim_static%x(1) - dim_static%x(0), ') does not match the ',  &
                                      'respective model model grid spacing (=', dx, ')'
          CALL message( 'topography_mod', 'PAC0334', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( ABS( dim_static%y(1) - dim_static%y(0) - dy ) > 0.001_wp * dy )  THEN
          WRITE( message_string, * )  'static driver: horizontal grid spacing in y-direction (=',  &
                                      dim_static%y(1) - dim_static%y(0), ') does not match the ',  &
                                      'respective model model grid spacing (=', dy, ')'
          CALL message( 'topography_mod', 'PAC0334', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Terrain height. First, get variable-related _fillvalue attribute.
       IF ( check_existence( var_names, 'zt' ) )  THEN
          terrain_height_f%from_file = .TRUE.
          CALL get_attribute( id_topo, char_fill, terrain_height_f%fill, .FALSE., 'zt' )
!
!--       Input 2D terrain height.
          ALLOCATE( terrain_height_f%var(nys:nyn,nxl:nxr)  )

          CALL get_variable( id_topo, 'zt', terrain_height_f%var, nxl, nxr, nys, nyn )
       ELSE
          terrain_height_f%from_file = .FALSE.
       ENDIF

!
!--    Read building height. First, read its _fillvalue attribute, as well as lod attribute.
       buildings_f%from_file = .FALSE.
       IF ( check_existence( var_names, 'buildings_2d' ) )  THEN
          buildings_f%from_file = .TRUE.
          CALL get_attribute( id_topo, char_lod, buildings_f%lod, .FALSE., 'buildings_2d' )
          CALL get_attribute( id_topo, char_fill, buildings_f%fill1, .FALSE., 'buildings_2d' )

!
!--       Read 2D buildings.
          IF ( buildings_f%lod == 1 )  THEN
             ALLOCATE( buildings_f%var_2d(nys:nyn,nxl:nxr) )
             CALL get_variable( id_topo, 'buildings_2d', buildings_f%var_2d, nxl, nxr, nys, nyn )
          ELSE
             WRITE( message_string, * )  'static driver: wrong netCDF attribute lod (=',           &
                                         buildings_f%lod, ') for buildings_2d'
             CALL message( 'topography_mod', 'PAC0335', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    If available, also read 3D building information. If both are available, use 3D information.
       IF ( check_existence( var_names, 'buildings_3d' ) )  THEN
          buildings_f%from_file = .TRUE.
          CALL get_attribute( id_topo, char_lod, buildings_f%lod, .FALSE., 'buildings_3d' )
          CALL get_attribute( id_topo, char_fill, buildings_f%fill2, .FALSE., 'buildings_3d' )
          CALL get_dimension_length( id_topo, buildings_f%nz, 'z' )
!
!--       Read 3D buildings
          IF ( buildings_f%lod == 2 )  THEN
             ALLOCATE( buildings_f%z(nzb:buildings_f%nz-1) )
             CALL get_variable( id_topo, 'z', buildings_f%z )
!
!--          Check if building information is consistent to numeric grid.
             IF ( buildings_f%nz > SIZE( zu ) )  THEN
                WRITE( message_string, * ) 'static driver: too much data points (=',               &
                                           buildings_f%nz, ') along the vertical coordinate for',  &
                                           ' 3d building data (maximum allowed=', SIZE( zu ), ')'
                CALL message( 'topography_mod', 'PAC0336', 2, 2, 0, 6, 0 )
             ENDIF
             IF ( ANY( ABS( buildings_f%z(0:buildings_f%nz-1) - zu(0:buildings_f%nz-1) ) >         &
                       0.001_wp * MINVAL( dz(1:number_dz) ) ) )  THEN
                k = 0
                DO  WHILE ( k <= buildings_f%nz-1 )
                   IF ( ABS( buildings_f%z(k) - zu(k) ) >                                          &
                        0.001_wp * MINVAL( dz(1:number_dz) ) )  EXIT
                   k = k + 1
                ENDDO
                WRITE( message_string, * ) 'static driver: vertical coordinate do not match ',     &
                                           'numeric grid at z(', k, ') for 3d building data'
                CALL message( 'topography_mod', 'PAC0337', 2, 2, 0, 6, 0 )
             ENDIF

             ALLOCATE( buildings_f%var_3d(nzb:buildings_f%nz-1, nys:nyn,nxl:nxr) )
             buildings_f%var_3d = 0
             CALL get_variable( id_topo, 'buildings_3d', buildings_f%var_3d, nxl, nxr, nys, nyn, 0,&
                                buildings_f%nz-1 )
          ELSE
             WRITE( message_string, * )  'static driver: wrong netCDF attribute lod (=',           &
                                         buildings_f%lod, ') for buildings_3d'
             CALL message( 'topography_mod', 'PAC0338', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
!
!--    Read building IDs and its FillValue attribute. Further required for mapping buildings on top
!--    of orography.
       IF ( check_existence( var_names, 'building_id' ) )  THEN
          building_id_f%from_file = .TRUE.
          CALL get_attribute( id_topo, char_fill, building_id_f%fill, .FALSE., 'building_id' )
          ALLOCATE( building_id_f%var(nys:nyn,nxl:nxr) )
          CALL get_variable( id_topo, 'building_id', building_id_f%var, nxl, nxr, nys, nyn )
       ELSE
          building_id_f%from_file = .FALSE.
       ENDIF
!
!--    Read building_type and required attributes.
       IF ( check_existence( var_names, 'building_type' ) )  THEN
          building_type_f%from_file = .TRUE.
          CALL get_attribute( id_topo, char_fill, building_type_f%fill, .FALSE., 'building_type' )
          ALLOCATE( building_type_f%var(nys:nyn,nxl:nxr) )
          CALL get_variable( id_topo, 'building_type', building_type_f%var, nxl, nxr, nys, nyn )
       ELSE
          building_type_f%from_file = .FALSE.
       ENDIF
!
!--    Close topography input file.
       CALL close_input_file( id_topo )
#else
       CONTINUE
#endif
!
!-- ASCII input
    ELSEIF ( TRIM( topography ) == 'read_from_file' )  THEN

       DO  ii = 0, io_blocks-1
          IF ( ii == io_group )  THEN

             OPEN( 90, FILE='TOPOGRAPHY_DATA' // TRIM( coupling_char ), STATUS='OLD',              &
                       FORM='FORMATTED', IOSTAT=io_status )

             IF ( io_status > 0 )  THEN
                message_string = 'file TOPOGRAPHY_DATA' // TRIM(coupling_char) // ' does not exist'
                CALL message( 'topography_mod', 'PAC0339', 1, 2, 0, 6, 0 )
             ENDIF

!
!--          Read topography PE-wise. Rows are read from nyn to nys, columns are read from nxl to
!--          nxr. At first, ny-nyn rows need to be skipped.
             skip_n_rows = 0
             DO  WHILE ( skip_n_rows < ny - nyn )
                READ( 90, * )
                skip_n_rows = skip_n_rows + 1
             ENDDO
!
!--          Read data from nyn to nys and nxl to nxr. Therefore, skip column until nxl-1 is reached
             ALLOCATE( buildings_f%var_2d(nys:nyn,nxl:nxr) )
             DO  j = nyn, nys, -1

                READ( 90, *, IOSTAT=io_status )  ( dum, i = 0, nxl-1 ),                            &
                                                 ( buildings_f%var_2d(j,i), i = nxl, nxr )

                IF ( io_status > 0 )  THEN
                   WRITE( message_string, '(A,1X,I5,1X,A)' ) 'error reading line', ny-j+1,         &
                                          'of file TOPOGRAPHY_DATA' // TRIM( coupling_char )
                   CALL message( 'topography_mod', 'PAC0340', 2, 2, myid, 6, 0 )
                ELSEIF ( io_status < 0 )  THEN
                   WRITE( message_string, '(A,1X,I5)' ) 'end of line or file detected for '//      &
                               'file TOPOGRAPHY_DATA' // TRIM( coupling_char ) // ' at line', ny-j+1
                   CALL message( 'topography_mod', 'PAC0341', 2, 2, myid, 6, 0 )
                ENDIF

             ENDDO

             CLOSE( 90 )
             buildings_f%from_file = .TRUE.

          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO

    ENDIF
!
!-- End of CPU measurement.
    CALL cpu_log( log_point_s(83), 'NetCDF/ASCII input topo', 'stop' )

!
!-- Check for minimum requirement to setup building topography. If buildings are provided, also an
!-- ID and a type are required.
!-- Note that performing this check in check_parameters will be too late (data will be used for grid
!-- inititialization before).
    IF ( input_pids_static )  THEN
       IF ( buildings_f%from_file  .AND.  .NOT. building_id_f%from_file )  THEN
          message_string = 'static driver: building ID is missing'
          CALL message( 'topography_mod', 'PAC0342', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( terrain_height_f%from_file )  THEN
!
!--    Check orography for fill-values.
!--    For the moment, give an error message. More advanced methods, e.g. a nearest neighbor
!--    algorithm as used in GIS systems might be implemented later.
!--    Note: This check must be placed here as terrain_height_f is altered later.
       IF ( ANY( terrain_height_f%var == terrain_height_f%fill ) )  THEN
          message_string = 'static driver: fill value for variable zt found'
          CALL message( 'topography_mod', 'PAC0343', 2, 2, myid, 6, 0 )
       ENDIF
    ELSE
!
!--    In case no terrain height is provided by static input file, allocate array nevertheless and
!--    set terrain height to 0, which simplifies topography initialization.
       ALLOCATE( terrain_height_f%var(nys:nyn,nxl:nxr) )
       terrain_height_f%var = 0.0_wp
    ENDIF
!
!-- Finally, exchange 1 ghost point for building ID and type.
!-- In case of non-cyclic boundary conditions set Neumann conditions at the lateral boundaries.
    IF ( building_id_f%from_file )  THEN
       CALL add_ghost_layers( building_id_f%var )
       CALL exchange_horiz_2d_int( building_id_f%var, nys, nyn, nxl, nxr, nbgp )
       CALL set_lateral_neumann_bc( building_id_f%var )
    ENDIF

    IF ( building_type_f%from_file )  THEN
       CALL add_ghost_layers( building_type_f%var )
       CALL exchange_horiz_2d_byte( building_type_f%var, nys, nyn, nxl, nxr, nbgp )
       CALL set_lateral_neumann_bc( building_type_f%var )
    ENDIF

 END SUBROUTINE topography_input


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Filter topography. This subroutine contains two filter steps. First, one-grid point wide
!> structures are filled. More precisely, all fluid grid points that are surrounded by more than 4
!> topography grid points in the x-,y-,z-direction are filled. This filtering is applied
!> iteratively until no grid-point wide holes exit any more.
!> Such holes are suspected to lead to velocity blow-ups as the continuity equation on a discrete
!> grid cannot be fulfilled in such case.
!> In a second step, enclosed narrow cavities are filled. These are also suspected to lead to
!> numerical instabilities, in particular when surface scalar fluxes are not zero. This case, scalar
!> values can increase to unrealistic levels since almost no mixing between the narrow cavity and
!> its surroundings takes place. At the moment, enclosed cavities of up to 9 grid points are
!> filtered.
!> Attention: The first filter step removes narrow elongated structures like streets, while the
!> second step only filters small fluid cavities that are completely surrounded by topography/
!> buildings.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE filter_topography

    INTEGER(iwp) ::  bl              !< search bound at left subdomain boundary
    INTEGER(iwp) ::  bn              !< search bound at north subdomain boundary
    INTEGER(iwp) ::  br              !< search bound at right subdomain boundary
    INTEGER(iwp) ::  bs              !< search bound at south subdomain boundary
    INTEGER(iwp) ::  f               !< running index over all fluid indices flagged to be filtered
    INTEGER(iwp) ::  i               !< running index along x-direction
    INTEGER(iwp) ::  i_f             !< grid index of filtered grid point in x-direction
    INTEGER(iwp) ::  j               !< running index along y-direction
    INTEGER(iwp) ::  j_f             !< grid index of filtered grid point in y-direction
    INTEGER(iwp) ::  k               !< running index along z-direction
#if defined( __parallel )
    INTEGER(iwp) ::  ngp_yz          !< number of extended ghost points
#endif
    INTEGER(iwp) ::  num_cavity      !< number of narrow cavities that have been filled
    INTEGER(iwp) ::  num_hole        !< number of holes (in topography) resolved by only one grid point
    INTEGER(iwp) ::  num_wall        !< number of surrounding vertical walls for a single grid point
    INTEGER(iwp) ::  type_xz_ext_int !< derived MPI datatype for 3-D integer ghost-point exchange with extended number of ghost points - left / right
    INTEGER(iwp) ::  type_yz_ext_int !< derived MPI datatype for 3-D integer ghost-point exchange with extended number of ghost points - south / north

    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  cav_filled  !< flag indiating whether a cavity has been filled or not

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  topo_tmp !< temporary 3D-topography used to fill holes

    LOGICAL ::  filled = .FALSE. !< flag indicating if holes were filled

    TYPE filter_type
       INTEGER(iwp) ::  num_gp         !< number of fluid grid point in current trace
       INTEGER(iwp) ::  num_thresh = 9 !< maximum number of grid points forming a cavity that will be filtered

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  i     !< grid index in x-direction indicating fluid grid point in current trace
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  j     !< grid index in y-direction indicating fluid grid point in current trace
       INTEGER(idp), DIMENSION(:), ALLOCATABLE ::  ij_id !< unique ID for each (ji)-pair
    END TYPE filter_type

    TYPE(filter_type) ::  filter !< derived structure summarizing all fluid-grid point in current trace

!
!-- Before checking for holes, set lateral boundary conditions for topography. After hole-filling,
!-- boundary conditions must be set again. Several iterations are performed, in order to fill
!-- holes which might emerge by applying the filling-algorithm itself.
!-- Attention: Beside small holes, also narrow elongated structures without a length limit (like
!-- streets with a width of only one gridpoint) will be filtered!
    ALLOCATE( topo_tmp(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    topo_tmp = 0

    num_hole = 99999
    DO WHILE ( num_hole > 0 )

       num_hole = 0
       CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
!
!--    Exchange also building ID and type. Note, building_type is an one-byte variable.
       IF ( building_id_f%from_file )  THEN
          CALL exchange_horiz_2d_int( building_id_f%var, nys, nyn, nxl, nxr, nbgp )
          CALL set_lateral_neumann_bc( building_id_f%var )
       ENDIF
       IF ( building_type_f%from_file )  THEN
          CALL exchange_horiz_2d_byte( building_type_f%var, nys, nyn, nxl, nxr, nbgp )
          CALL set_lateral_neumann_bc( building_type_f%var )
       ENDIF

       topo_tmp = topo
!
!--    In case of non-cyclic lateral boundaries, assume lateral boundary to be a solid wall. Thus,
!--    intermediate spaces of one grid point between boundary and some topographic structure will be
!--    filled.
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( nys == 0  )  topo_tmp(:,-1,:)   = IBSET( topo_tmp(:,0,:),  0 )
          IF ( nyn == ny )  topo_tmp(:,ny+1,:) = IBSET( topo_tmp(:,ny,:), 0 )
       ENDIF

       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( nxl == 0  )  topo_tmp(:,:,-1)   = IBSET( topo_tmp(:,:,0),  0 )
          IF ( nxr == nx )  topo_tmp(:,:,nx+1) = IBSET( topo_tmp(:,:,nx), 0 )
       ENDIF

       num_hole = 0
       DO i = nxl, nxr
          DO j = nys, nyn
             DO  k = nzb+1, nzt
                IF ( .NOT. BTEST( topo_tmp(k,j,i), 0 ) )  THEN
                   num_wall = 0
                   IF ( BTEST( topo_tmp(k,j-1,i), 0 ) )  num_wall = num_wall + 1
                   IF ( BTEST( topo_tmp(k,j+1,i), 0 ) )  num_wall = num_wall + 1
                   IF ( BTEST( topo_tmp(k,j,i-1), 0 ) )  num_wall = num_wall + 1
                   IF ( BTEST( topo_tmp(k,j,i+1), 0 ) )  num_wall = num_wall + 1
                   IF ( BTEST( topo_tmp(k-1,j,i), 0 ) )  num_wall = num_wall + 1
                   IF ( BTEST( topo_tmp(k+1,j,i), 0 ) )  num_wall = num_wall + 1

                   IF ( num_wall >= 4 )  THEN
                      num_hole = num_hole + 1
!
!--                   Clear flag 0 and set special flag ( bit 4) to indicate that new topography
!--                   point is a result of filtering process.
                      topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                      topo(k,j,i) = IBSET( topo(k,j,i), 4 )
!
!--                   If filled grid point is occupied by a building, classify it as building grid
!--                   point.
                      IF ( building_type_f%from_file )  THEN
                         IF ( building_type_f%var(j,i)   /= building_type_f%fill  .OR.             &
                              building_type_f%var(j+1,i) /= building_type_f%fill  .OR.             &
                              building_type_f%var(j-1,i) /= building_type_f%fill  .OR.             &
                              building_type_f%var(j,i+1) /= building_type_f%fill  .OR.             &
                              building_type_f%var(j,i-1) /= building_type_f%fill )                 &
                         THEN
!
!--                         Set flag indicating building surfaces
                            topo(k,j,i) = IBSET( topo(k,j,i), 2 )
!
!--                         Set building_type and ID at this position if not already set. This is
!--                         required for proper initialization of urban-surface energy balance
!--                         solver.
                            IF ( building_type_f%var(j,i) == building_type_f%fill )  THEN

                               IF ( building_type_f%var(j+1,i) /= building_type_f%fill )  THEN
                                  building_type_f%var(j,i) = building_type_f%var(j+1,i)
                                  building_id_f%var(j,i)   = building_id_f%var(j+1,i)
                               ELSEIF ( building_type_f%var(j-1,i) /= building_type_f%fill )  THEN
                                  building_type_f%var(j,i) = building_type_f%var(j-1,i)
                                  building_id_f%var(j,i)   = building_id_f%var(j-1,i)
                               ELSEIF ( building_type_f%var(j,i+1) /= building_type_f%fill )  THEN
                                  building_type_f%var(j,i) = building_type_f%var(j,i+1)
                                  building_id_f%var(j,i)   = building_id_f%var(j,i+1)
                               ELSEIF ( building_type_f%var(j,i-1) /= building_type_f%fill )  THEN
                                  building_type_f%var(j,i) = building_type_f%var(j,i-1)
                                  building_id_f%var(j,i)   = building_id_f%var(j,i-1)
                               ENDIF
                            ENDIF
                         ENDIF
                      ENDIF
!
!--                   If filled grid point is already classified as building everything is fine,
!--                   else classify this grid point as natural type grid point. This case, values
!--                   for the surface type are already set.
                      IF ( .NOT. BTEST( topo(k,j,i), 2 ) )  THEN
                         topo(k,j,i) = IBSET( topo(k,j,i), 1 )
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!--    Count the total number of holes, required for informative message.
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, num_hole, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#endif
       IF ( num_hole > 0  .AND.  .NOT. filled )  filled = .TRUE.

    ENDDO
!
!-- Create an informative message if any 1-grid point wide holes were filled.
    IF ( filled )  THEN
       WRITE( message_string, '(A,I6,A)' )                                                         &
              'topography was filtered: &', num_hole, ' hole(s) resolved by only one grid ' //     &
              'point were filled during initialization'
       CALL message( 'topography_mod', 'PAC0344', 0, 0, 0, 6, 0 )
    ENDIF

    DEALLOCATE( topo_tmp )
!
!-- Finally, exchange topo array again and if necessary set Neumann boundary condition in case of
!-- non-cyclic lateral boundaries.
    CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
    CALL topography_set_non_cyc_bc( topo )
!
!-- Exchange building ID and type. Note, building_type is an one-byte variable.
    IF ( building_id_f%from_file )  THEN
       CALL exchange_horiz_2d_int( building_id_f%var, nys, nyn, nxl, nxr, nbgp )
       CALL set_lateral_neumann_bc( building_id_f%var )
    ENDIF
    IF ( building_type_f%from_file )  THEN
       CALL exchange_horiz_2d_byte( building_type_f%var, nys, nyn, nxl, nxr, nbgp )
       CALL set_lateral_neumann_bc( building_type_f%var )
    ENDIF


!
!-- On top of this 1-gridpoint cavity filtering, also larger structures are filtered, i.e.
!-- enclosed cavities of less than 9 gridpoints.
!-- Note, at the moment only courtyard-like cavities in the xy-plance are filtered, while
!-- cavities in the yz- and xz-plane are not treated at the moment.
!-- In a first step, search for all grid points that are somehow topography bounded in their
!-- vicinity and store their indices. Do this layer by layer, which is the reason for the reverse
!-- loop structure.
!-- First of all, determine the search bounds at the subdomain boundaries. Note, for the filter
!-- algorithm extended ghost layers are requied.
    bl = MERGE( nxl-filter%num_thresh, 0,    .NOT. ( bc_dirichlet_l  .OR.  bc_radiation_l ) )
    br = MERGE( nxr+filter%num_thresh, nx+1, .NOT. ( bc_dirichlet_r  .OR.  bc_radiation_r ) )
    bn = MERGE( nyn+filter%num_thresh, ny+1, .NOT. ( bc_dirichlet_n  .OR.  bc_radiation_n ) )
    bs = MERGE( nys-filter%num_thresh, 0,    .NOT. ( bc_dirichlet_s  .OR.  bc_radiation_s ) )

!
!-- Define MPI-datatypes for extended ghost-point exchange. This is necessary to detect also
!-- elongated cavities with length of up to 9 grid points.
!-- However, extended ghost layers require the subdomains in x- and y-direction to hold at least
!-- the same number of grid points. Else (at least at the moment), ghost point exchange will not
!-- work. Hence, in case of smaller subdomains, the number of ghost points and the filter
!-- threshold must be reduced. In unlucky cases this can result in the situation that elongated
!-- cavities are not fully filtered.
    filter%num_thresh = MIN( filter%num_thresh, nyn-nys+1, nxr-nxl+1 )

#if defined( __parallel )
    ngp_yz = ( nzt - nzb + 2 ) * ( nyn - nys + 1 + 2 * filter%num_thresh )

    CALL MPI_TYPE_VECTOR( nxr-nxl+1+2*filter%num_thresh, filter%num_thresh*(nzt-nzb+2), ngp_yz,    &
                          MPI_INTEGER, type_xz_ext_int, ierr )
    CALL MPI_TYPE_COMMIT( type_xz_ext_int, ierr )

    CALL MPI_TYPE_VECTOR( filter%num_thresh, ngp_yz, ngp_yz, MPI_INTEGER, type_yz_ext_int, ierr )
    CALL MPI_TYPE_COMMIT( type_yz_ext_int, ierr )
!
!-- Set a barrier so that all MPI datatypes are defined before the ghost-point exchange starts.
!-- Without such a barrier, a MPI error concerning insufficiently large buffer size may occur.
    CALL MPI_BARRIER( comm2d, ierr )
#endif

    ALLOCATE( topo_tmp(nzb:nzt+1,nys-filter%num_thresh:nyn+filter%num_thresh,                      &
                       nxl-filter%num_thresh:nxr+filter%num_thresh) )


    topo_tmp = 0
    topo_tmp(:,nys:nyn,nxl:nxr) = topo(:,nys:nyn,nxl:nxr)
    CALL exchange_horiz_int( topo_tmp, nys, nyn, nxl, nxr, nzt, filter%num_thresh,                 &
                             type_xz_ext_int, type_yz_ext_int )
!
!-- Allocate arrays containing the flagged grid-indices and their corresponding grid-point ID
!-- that are flagged to be potentially filtered.
    ALLOCATE( filter%i(filter%num_thresh)     )
    ALLOCATE( filter%j(filter%num_thresh)     )
    ALLOCATE( filter%ij_id(filter%num_thresh) )

    cav_filled = 0
    DO  k = nzb+1, nzt
       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Only employ filter algorithm to non-topography grid points..
             IF ( .NOT. BTEST( topo_tmp(k,j,i), 0 ) )  THEN
!
!--             Starting from each grid point, search
                filter%num_gp = 1
                filter%i(:)     = -HUGE( 1 )
                filter%j(:)     = -HUGE( 1 )
                filter%ij_id(:) = -HUGE( 1 )

                filter%i(1)     = i
                filter%j(1)     = j
                filter%ij_id(1) = gridpoint_id( j, i )
!
!--             Search for fluid grid points in the surroundings of (j,i). Note, the search is
!--             done recursively. If a surrounding grid point is flagged as fluid, further extend
!--             the search in the surrounding of this flagged grid point, and so on. The search
!--             for fluid grid points will be stopped if a treshold value of fluid grid point
!--             is reached. This case, no topography filter is employed. However, if the number
!--             of continuously connected fluid grid points in the farther surrounding of (j,i) is
!--             less than this threshold number, topography will be filtered and the cavity is
!--             filled.
                CALL trace_surrounding_gridpoints( k, j, i )

                IF ( filter%num_gp < filter%num_thresh )  THEN
                   cav_filled(j,i) = 1

                   DO  f = 1, filter%num_gp
                      i_f = filter%i(f)
                      j_f = filter%j(f)
!
!--                   Set topography bit at flagged grid points. Only set topography if the
!--                   k-1 level is also topography. By this condition, filling
!--                   of courtyards with lateral openings, urn-shaped cavities, or elevated
!--                   "bottlenecks", which are narrower in their upper levels should be avoided.
!--                   Same as for the 1 grid point hole filling algorithm above, set bit 4 to
!--                   indicate filtered grid points.
                      IF ( BTEST( topo_tmp(k-1,j_f,i_f), 0 ) )  THEN
                         topo_tmp(k,j_f,i_f) = IBSET( topo_tmp(k,j_f,i_f), 0 )
                         topo_tmp(k,j_f,i_f) = IBSET( topo_tmp(k,j_f,i_f), 4 )

!
!--                      If filled grid point is occupied by a building, classify it as building
!--                      grid point.
                         IF ( building_type_f%from_file )  THEN
                            IF ( j_f >= nysg  .AND.  j_f <= nyng  .AND.                            &
                                 i_f >= nxlg  .AND.  i_f <= nxrg )  THEN
                               IF ( building_type_f%var(j_f,i_f) /= building_type_f%fill )  THEN
!
!--                               Set flag indicating building surfaces.
                                  topo_tmp(k,j_f,i_f) = IBSET( topo_tmp(k,j_f,i_f), 2 )
                               ENDIF
                            ENDIF
                         ENDIF
!
!--                      If filled grid point is already classified as building everything is fine,
!--                      else classify this grid point as natural type grid point. This case,
!--                      values for the surface type are already set.
                         IF ( .NOT. BTEST( topo_tmp(k,j_f,i_f), 2 ) )  THEN
                            topo_tmp(k,j_f,i_f) = IBSET( topo_tmp(k,j_f,i_f), 1 )
                         ENDIF

                      ENDIF
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE ( filter%i     )
    DEALLOCATE ( filter%j     )
    DEALLOCATE ( filter%ij_id )
!
!-- Count the total number of filled cavities, required for informative message.
    num_cavity = SUM( cav_filled )
#if defined( __parallel )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, num_cavity, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#endif
!
!-- Create an informative message if any narrow cavities were filled.
    IF ( num_cavity > 0 )  THEN
       WRITE( message_string, * ) num_cavity, 'narrow cavities of less than ', filter%num_thresh,  &
                                  'horizontal grid points were filled during initialization'
       CALL message( 'topography_mod', 'PAC0345', 0, 0, 0, 6, 0 )
    ENDIF

    topo(:,nys:nyn,nxl:nxr) = topo_tmp(:,nys:nyn,nxl:nxr)
!
!-- Finally, exchange topo array again and if necessary set Neumann boundary condition in case of
!-- non-cyclic lateral boundaries.
    CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
    CALL topography_set_non_cyc_bc( topo )
!
!-- Exchange building ID and type. Note, building_type is an one-byte variable.
    IF ( building_id_f%from_file )  THEN
       CALL exchange_horiz_2d_int( building_id_f%var, nys, nyn, nxl, nxr, nbgp )
       CALL set_lateral_neumann_bc( building_id_f%var )
    ENDIF
    IF ( building_type_f%from_file )  THEN
       CALL exchange_horiz_2d_byte( building_type_f%var, nys, nyn, nxl, nxr, nbgp )
       CALL set_lateral_neumann_bc( building_type_f%var )
    ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Core of the cavity filter algorithm. This routine checks if horizontally surrounding grid points
!> of (j,i) belong to fluid or topography. If they belong to fluid, recursively extend the search
!> until a threshold value is reached (meaning that the cavity is wide enough), or until all farther
!> surrounding fluid grid points are counted.
!--------------------------------------------------------------------------------------------------!
    RECURSIVE SUBROUTINE trace_surrounding_gridpoints( k, j, i )

       INTEGER(iwp) ::  i       !< input grid index in x-direction
       INTEGER(iwp) ::  i_trace !< grid index in x-direction for next trace grid cell
       INTEGER(iwp) ::  j       !< input grid index in y-direction
       INTEGER(iwp) ::  j_trace !< grid index in y-direction for next trace grid cell
       INTEGER(iwp) ::  k       !< input grid index in z-direction
       INTEGER(iwp) ::  n       !< loop variable for the 4 grid-line directions

       INTEGER(iwp), DIMENSION(4) ::  off_x = (/ -1, 1,  0, 0 /)
       INTEGER(iwp), DIMENSION(4) ::  off_y = (/  0, 0, -1, 1 /)

!
!--    From the input grid point, start to look around the 4 grid-line surrounding grid points.
       DO  n = 1, SIZE( off_x )
!
!--       Exit subroutine if more than filter%num_thresh (currently 9) fluid grid points are found.
          IF ( filter%num_gp >= filter%num_thresh )  RETURN

          i_trace = i + off_x(n)
          j_trace = j + off_y(n)
!
!--       Skip trace coordinates that are out of the search- or subdomain boundaries.
          IF ( j_trace < bs  .OR.  j_trace > bn  .OR.  i_trace < bl  .OR.  i_trace > br )  CYCLE
!
!--       Check if the grid point has been already counted. This is identified by a unique ID each
!--       (j,i) is given.
          IF ( .NOT. ANY( gridpoint_id( j_trace, i_trace ) == filter%ij_id ) )  THEN
!
!--          If the trace grid point (j_trace,i_trace) is a fluid grid point, add its
!--          coordinates to the already existing list of coordinates.
             IF ( .NOT. BTEST( topo_tmp(k,j_trace,i_trace), 0 ) )  THEN
!
!--             Increment the fluid grid point counter.
                filter%num_gp = filter%num_gp + 1
!
!--             Add new index pair and a corresponding unique ID.
                filter%i(filter%num_gp)     = i_trace
                filter%j(filter%num_gp)     = j_trace
                filter%ij_id(filter%num_gp) = gridpoint_id( j_trace, i_trace )
!
!--             Based on the updated (j_trace,i_trace) coordinate, further search for
!--             fluid grid points in its vicinity.
                IF ( filter%num_gp < filter%num_thresh )                                           &
                   CALL trace_surrounding_gridpoints( k, j_trace, i_trace )
             ENDIF
          ENDIF
       ENDDO

    END SUBROUTINE trace_surrounding_gridpoints

 END SUBROUTINE filter_topography


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Set static topography flags
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE topography_set_flags

    INTEGER(iwp) ::  i    !< index variable along x
    INTEGER(iwp) ::  ibit !< integer bit position of topgraphy masking array
    INTEGER(iwp) ::  j    !< index variable along y
    INTEGER(iwp) ::  k    !< index variable along z


    ALLOCATE( topo_flags(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    topo_flags = 0
!
!-- Set-up topography flags. First, set flags only for s, u, v and w-grid.
!-- Further special flags will be set in following loops. Note, topography is defined when the
!-- temporary array "topo", bit 0 is one, else it is atmosphere.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt+1
!
!--          scalar grid
             IF ( .NOT. BTEST( topo(k,j,i), 0 ) )  topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 0 )
          ENDDO
!
!--       Cartesian topography.
          DO  k = nzb, nzt+1
!
!--          u grid.
             IF ( .NOT. BTEST( topo(k,j,i), 0 )  .AND.  .NOT. BTEST( topo(k,j,i-1), 0 ) )          &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 1 )
!
!--          v grid.
             IF ( .NOT. BTEST( topo(k,j,i), 0 )  .AND.  .NOT. BTEST( topo(k,j-1,i), 0 ) )          &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 2 )
          ENDDO

          DO k = nzb, nzt
!
!--          w grid.
             IF ( .NOT. BTEST( topo(k,j,i), 0 )  .AND.  .NOT. BTEST( topo(k+1,j,i), 0 ) )          &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 3 )
          ENDDO

          IF ( topography /= 'closed_channel' )  THEN
             topo_flags(nzt+1,j,i) = IBSET( topo_flags(nzt+1,j,i), 3 )
          ENDIF

       ENDDO
    ENDDO

    CALL exchange_horiz_int( topo_flags, nys, nyn, nxl, nxr, nzt, nbgp )

!
!-- Set outer array for scalars to mask near-surface grid points. Note, on basis of flag 24 futher
!-- flags will be derived which are used to control production of subgrid TKE production near walls.
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
             IF ( BTEST( topo_flags(k,j-1,i), 0 )    .AND.                                         &
                  BTEST( topo_flags(k,j+1,i), 0 )    .AND.                                         &
                  BTEST( topo_flags(k,j,i-1), 0 )    .AND.                                         &
                  BTEST( topo_flags(k,j,i+1), 0 )    .AND.                                         &
                  BTEST( topo_flags(k,j-1,i-1), 0 )  .AND.                                         &
                  BTEST( topo_flags(k,j+1,i-1), 0 )  .AND.                                         &
                  BTEST( topo_flags(k,j-1,i+1), 0 )  .AND.                                         &
                  BTEST( topo_flags(k,j+1,i+1), 0 ) )                                              &
             THEN
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 24 )
             ENDIF
          ENDDO
       ENDDO
    ENDDO
!
!-- Set further special flags.
    DO i = nxl, nxr
       DO j = nys, nyn
          DO k = nzb, nzt+1
!
!--          Scalar grid, former nzb_diff_s_inner.
!--          Note, use this flag also to mask topography in diffusion_u and diffusion_v along the
!--          vertical direction. In case of use_surface_fluxes, fluxes are calculated via MOST,
!--          else, simple gradient approach is applied. Please note, in case of u- and v-diffuison,
!--          a small error is made at edges (on the east side for u, at the north side for v), since
!--          topography on scalar grid point is used instead of topography on u/v-grid. As number of
!--          topography grid points on uv-grid is different than s-grid, different number of surface
!--          elements would be required. In order to avoid this, treat edges (u(k,j,i+1)) simply by
!--          a gradient approach, i.e. these points are not masked within diffusion_u. Tests had
!--          shown that the effect on the flow is negligible.
             IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
                IF ( BTEST( topo_flags(k,j,i), 0 ) )                                               &
                   topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 8 )
             ELSE
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 8 )
             ENDIF

          ENDDO
!
!--       Special flag to control vertical diffusion at model top - former nzt_diff.
          topo_flags(:,j,i) = IBSET( topo_flags(:,j,i), 9 )
          IF ( use_top_fluxes )  topo_flags(nzt+1,j,i) = IBCLR( topo_flags(nzt+1,j,i), 9 )

          DO k = nzb+1, nzt
!
!--          Special flag on u grid, former nzb_u_inner + 1, required for disturb_field and
!--          initialization. Do not disturb directly at topography, as well as initialize u with
!--          zero one grid point outside of topography.
             IF ( BTEST( topo_flags(k-1,j,i), 1 )  .AND.                                           &
                  BTEST( topo_flags(k,j,i),   1 )  .AND.                                           &
                  BTEST( topo_flags(k+1,j,i), 1 ) )                                                &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 20 )
!
!--          Special flag on v grid, former nzb_v_inner + 1, required for disturb_field and
!--          initialization. Do not disturb directly at topography, as well as initialize v with
!--          zero one grid point outside of topography.
             IF ( BTEST( topo_flags(k-1,j,i), 2 )  .AND.                                           &
                  BTEST( topo_flags(k,j,i),   2 )  .AND.                                           &
                  BTEST( topo_flags(k+1,j,i), 2 ) )                                                &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 21 )
!
!--          Special flag on scalar grid, former nzb_s_inner+1. Used for lpm_sgs_tke.
             IF ( BTEST( topo_flags(k,j,i),   0 )  .AND.                                           &
                  BTEST( topo_flags(k-1,j,i), 0 )  .AND.                                           &
                  BTEST( topo_flags(k+1,j,i), 0 ) )                                                &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 25 )
!
!--          Special flag on scalar grid, nzb_diff_s_outer - 1, required in in production_e.
             IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
                IF ( BTEST( topo_flags(k,j,i),   24 )  .AND.                                       &
                     BTEST( topo_flags(k-1,j,i), 24 )  .AND.                                       &
                     BTEST( topo_flags(k+1,j,i), 0 ) )                                             &
                   topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 29 )
             ELSE
                IF ( BTEST( topo_flags(k,j,i), 0 ) )                                               &
                   topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 29 )
             ENDIF
!
!--          Special flag on scalar grid, nzb_diff_s_outer - 1, required in
!--          in production_e.
             IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
                IF ( BTEST( topo_flags(k,j,i),   0 )  .AND.                                        &
                     BTEST( topo_flags(k-1,j,i), 0 )  .AND.                                        &
                     BTEST( topo_flags(k+1,j,i), 0 ) )                                             &
                   topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 30 )
             ELSE
                IF ( BTEST( topo_flags(k,j,i), 0 ) )                                               &
                   topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 30 )
             ENDIF
          ENDDO
!
!--       Flags indicating downward facing walls.
          DO k = nzb+1, nzt+1
!
!--          Scalar grid.
             IF (       BTEST( topo_flags(k-1,j,i), 0 )  .AND.                                     &
                  .NOT. BTEST( topo_flags(k,j,i), 0   ) )                                          &
                 topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 13 )
!
!--          Downward facing wall on u grid.
             IF (       BTEST( topo_flags(k-1,j,i), 1 )  .AND.                                     &
                  .NOT. BTEST( topo_flags(k,j,i), 1   ) )                                          &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 15 )
!
!--          Downward facing wall on v grid.
             IF (       BTEST( topo_flags(k-1,j,i), 2 )  .AND.                                     &
                  .NOT. BTEST( topo_flags(k,j,i), 2   ) )                                          &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 17 )
!
!--          Downward facing wall on w grid.
             IF (       BTEST( topo_flags(k-1,j,i), 3 )  .AND.                                     &
                  .NOT. BTEST( topo_flags(k,j,i), 3 ) )                                            &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 19 )
          ENDDO
!
!--       Flags indicating upward facing walls.
          DO k = nzb, nzt
!
!--          Upward facing wall on scalar grid.
             IF ( .NOT. BTEST( topo_flags(k,j,i),   0 )  .AND.                                     &
                        BTEST( topo_flags(k+1,j,i), 0 ) )                                          &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 12 )
!
!--          Upward facing wall on u grid.
             IF ( .NOT. BTEST( topo_flags(k,j,i),   1 )  .AND.                                     &
                        BTEST( topo_flags(k+1,j,i), 1 ) )                                          &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 14 )

!
!--          Upward facing wall on v grid.
             IF ( .NOT. BTEST( topo_flags(k,j,i),   2 )  .AND.                                     &
                        BTEST( topo_flags(k+1,j,i), 2 ) )                                          &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 16 )

!
!--          Upward facing wall on w grid.
             IF ( .NOT. BTEST( topo_flags(k,j,i),   3 )  .AND.                                     &
                        BTEST( topo_flags(k+1,j,i), 3 ) )                                          &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 18 )
!
!--          Special flag on scalar grid, former nzb_s_inner.
             IF ( BTEST( topo_flags(k,j,i), 0 )  .OR.                                              &
                  BTEST( topo_flags(k,j,i), 12 ) .OR.                                              &
                  BTEST( topo_flags(k,j,i), 13 ) )                                                 &
                   topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 22 )
!
!--          Special flag on scalar grid, nzb_diff_s_inner - 1, required for flow_statistics.
             IF ( constant_flux_layer  .OR.  use_surface_fluxes )  THEN
                IF ( BTEST( topo_flags(k,j,i),   0 )  .AND.                                        &
                     BTEST( topo_flags(k+1,j,i), 0 ) )                                             &
                  topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 23 )
             ELSE
                IF ( BTEST( topo_flags(k,j,i), 22 ) )                                              &
                   topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 23 )
             ENDIF


          ENDDO
          topo_flags(nzt+1,j,i) = IBSET( topo_flags(nzt+1,j,i), 22 )
          topo_flags(nzt+1,j,i) = IBSET( topo_flags(nzt+1,j,i), 23 )
!
!--       Set flags indicating that topography is close by in horizontal direction, i.e. flags that
!--       infold the topography. These will be used to set advection flags for passive scalars,
!--       where due to large gradients near buildings stationary numerical oscillations can produce
!--       unrealistically high concentrations. This is only necessary if WS-scheme is applied for
!--       scalar advection. Note, these flags will be only used for passive scalars such as chemical
!--       species or aerosols.
          IF ( scalar_advec == 'ws-scheme' )  THEN
             DO  k = nzb, nzt
                IF ( BTEST( topo_flags(k,j,i), 0 )  .AND. (                                        &
                     ANY( .NOT. BTEST( topo_flags(k,j-3:j+3,i-1), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j-3:j+3,i-2), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j-3:j+3,i-3), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j-3:j+3,i+1), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j-3:j+3,i+2), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j-3:j+3,i+3), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j-1,i-3:i+3), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j-2,i-3:i+3), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j-3,i-3:i+3), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j+1,i-3:i+3), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j+2,i-3:i+3), 0 ) )  .OR.                      &
                     ANY( .NOT. BTEST( topo_flags(k,j+3,i-3:i+3), 0 ) )                            &
                                                          ) )                                      &
                THEN
                   topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 31 )
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDDO
!
!-- Finally, set identification flags indicating natural terrain or buildings.
!-- Natural terrain grid points. Information on the type of the surface is stored in bit 1 of
!-- 3D Integer array topo. However, this bit is only set when topography is read from file. In order
!-- to run the land-surface model also without topography information, set bit 1 explicitely in this
!-- case.
!-- Natural terrain grid points.
!-- If no topography is initialized, the land-surface is at k = nzb.
    IF ( TRIM( topography ) /= 'read_from_file' )  THEN
       topo_flags(nzb,:,:) = IBSET( topo_flags(nzb,:,:), 5 )
    ELSE
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt+1
!
!--             Natural terrain grid point.
                IF ( BTEST( topo(k,j,i), 1 ) )                                                     &
                   topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 5 )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Building grid points.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt+1
             IF ( BTEST( topo(k,j,i), 2 ) )                                                        &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 6 )
          ENDDO
       ENDDO
    ENDDO
!
!-- Set flag 4, indicating new topography grid points due to filtering.
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt+1
             IF ( BTEST( topo(k,j,i), 4 ) )                                                        &
                topo_flags(k,j,i) = IBSET( topo_flags(k,j,i), 4 )
          ENDDO
       ENDDO
    ENDDO

!
!-- Exchange ghost points for wall flags
    CALL exchange_horiz_int( topo_flags, nys, nyn, nxl, nxr, nzt, nbgp )
!
!-- Set boundary conditions also for flags. Can be interpreted as Neumann boundary conditions for
!-- topography.
    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( nys == 0  )  THEN
          DO  i = 1, nbgp
             topo_flags(:,nys-i,:) = topo_flags(:,nys,:)
          ENDDO
       ENDIF
       IF ( nyn == ny )  THEN
          DO  i = 1, nbgp
             topo_flags(:,nyn+i,:) = topo_flags(:,nyn,:)
          ENDDO
       ENDIF
    ENDIF
    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( nxl == 0  )  THEN
          DO  i = 1, nbgp
             topo_flags(:,:,nxl-i) = topo_flags(:,:,nxl)
          ENDDO
       ENDIF
       IF ( nxr == nx )  THEN
          DO  i = 1, nbgp
             topo_flags(:,:,nxr+i) = topo_flags(:,:,nxr)
          ENDDO
       ENDIF
    ENDIF
!
!-- Pre-calculate topography top indices.
    ALLOCATE( topo_top_ind(nysg:nyng,nxlg:nxrg,0:6) )
!
!-- Uppermost topography index on scalar grid.
    ibit = 12
    topo_top_ind(:,:,0) = MAXLOC( MERGE( 1, 0, BTEST( topo_flags(:,:,:), ibit ) ), DIM=1 ) - 1
!
!-- Uppermost topography index on u grid.
    ibit = 14
    topo_top_ind(:,:,1) = MAXLOC( MERGE( 1, 0, BTEST( topo_flags(:,:,:), ibit ) ), DIM=1 ) - 1
!
!-- Uppermost topography index on v grid.
    ibit = 16
    topo_top_ind(:,:,2) = MAXLOC( MERGE( 1, 0, BTEST( topo_flags(:,:,:), ibit ) ), DIM=1 ) - 1
!
!-- Uppermost topography index on w grid.
    ibit = 18
    topo_top_ind(:,:,3) = MAXLOC( MERGE( 1, 0, BTEST( topo_flags(:,:,:), ibit ) ), DIM=1 ) - 1
!
!-- Uppermost topography index on scalar outer grid.
    ibit = 24
    topo_top_ind(:,:,4) = MAXLOC( MERGE( 1, 0, BTEST( topo_flags(:,:,:), ibit ) ), DIM=1 ) - 1
!
!-- Uppermost topography index including full-3D geometry.
    ibit = 12
    DO  k = nzb, nzt+1
       WHERE( BTEST( topo_flags(k,:,:), ibit ) )  topo_top_ind(:,:,5) = k
    ENDDO
!
!-- Pre-calculate top index of uppermost terrain grid point. This is used for terrain-following
!-- masked data output.
    topo_top_ind(:,:,6) = MINLOC( MERGE( 1, 0, BTEST( topo_flags, 5 ) ), DIM = 1 ) - 1

 END SUBROUTINE topography_set_flags


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Set temporary topography flags and reference buildings on top of underlying orography.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE process_topography

    INTEGER(iwp) ::  dim_builds       !< total number of buildings within the model domain
    INTEGER(iwp) ::  i                !< running index along x-direction
    INTEGER(iwp) ::  j                !< running index along y-direction
    INTEGER(iwp) ::  k                !< running index along z-direction with respect to numeric grid
    INTEGER(iwp) ::  k2               !< running index along z-direction with respect to netcdf grid
    INTEGER(iwp) ::  nr               !< index variable indication maximum terrain height for respective building ID
    INTEGER(iwp) ::  num_build        !< counter for number of buildings
    INTEGER(iwp) ::  topo_top_index   !< orography top index, used to map 3D buildings onto terrain

#if defined( __parallel )
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  displace_dum        !< displacements of start addresses, used for MPI_ALLGATHERV
#endif
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids           !< building IDs on entire model domain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_final     !< building IDs on entire model domain, multiple occurences are
                                                                    !< sorted out
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_final_tmp !< temporary array used for resizing
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_l         !< building IDs on local subdomain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_l_tmp     !< temporary array used to resize array of building IDs

    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_buildings     !< number of buildings with different ID on entire model domain
    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_buildings_l   !< number of buildings with different ID on local subdomain

    REAL(wp)                            ::  ocean_offset        !< offset to consider inverse vertical coordinate at topography
                                                                !< definition
    REAL(wp)                            ::  oro_min = 0.0_wp    !< minimum terrain height in entire model domain, used to reference
                                                                !< terrain to zero

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  oro_max             !< maximum terrain height occupied by an building with certain id
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  oro_max_l           !< maximum terrain height occupied by an building with certain id,

                                                                !< on local subdomain
!
!-- Reference lowest terrain height to zero. This ensures that first, non-required gird levels
!-- (those which lie entirely below the minimum orography) are avoided, and second, that also
!-- negative orography can be used within the input file.
!-- Please note, in case of a nested run, the global minimum from all parent and childs needs to be
!-- removed to avoid steep edges at the child-domain boundaries.
!-- Moreover, please note, the global minimum topography is only substracted if topography is
!-- defined in all coupled models. If it is defined only in some of the models, this step is
!-- skipped, same as in coupled atmosphere-ocean runs where topography height relates to
!-- different reference levels.
    IF ( topo_read_all_domains  .AND.  .NOT. atmosphere_ocean_coupled_run )  THEN

       IF ( input_pids_static )  THEN
#if defined( __parallel )
          CALL MPI_ALLREDUCE( MINVAL( terrain_height_f%var ), oro_min, 1, MPI_REAL, MPI_MIN,       &
                              MPI_COMM_WORLD, ierr )
#else
          oro_min = MINVAL( terrain_height_f%var )
#endif
          terrain_height_f%var = terrain_height_f%var - oro_min
!
!--       Update reference height used within output files
          init_model%origin_z = init_model%origin_z + oro_min
!
!--    ASCII topography branch. In this case, in contrast to the static driver input, topography is
!--    input via the variable buildings_2d, which is used as a dummy here. This case, the minimum
!--    building height is subtracted from the building array.
       ELSE
#if defined( __parallel )
          CALL MPI_ALLREDUCE( MINVAL( buildings_f%var_2d ), oro_min, 1, MPI_REAL, MPI_MIN,         &
                              MPI_COMM_WORLD, ierr )
#else
          oro_min = MINVAL( buildings_f%var_2d )
#endif
          buildings_f%var_2d = buildings_f%var_2d - oro_min
!
!--       Update reference height used within output files.
          init_model%origin_z = init_model%origin_z + oro_min
       ENDIF

    ENDIF

!
!-- In the following, buildings and orography are further preprocessed before they are mapped on the
!-- LES grid.
!-- Buildings are mapped on top of the orography by maintaining the roof shape of the building. This
!-- can be achieved by referencing building on top of the maximum terrain height within the area
!-- occupied by the respective building. As buildings and terrain height are defined only locally,
!-- parallelization of this referencing is required (a building can be distributed between different
!-- cores).
!-- In a first step, determine the number of buildings with different building id on each PE. In a
!-- next step, all building ids are gathered into one array which is present to all PEs. For each
!-- building ID, the maximum terrain height occupied by the respective building is computed and
!-- distributed to each PE.
!-- Finally, for each building id and its respective reference orography, builidings are mapped on
!-- top.
!--
!-- First, set topography flags. Bit 1 indicates orography, bit 2 buildings.
!-- Grid point nzb is topography on the staggered grid, but
!-- in case of vertically shifted nests, the boundary vertical coordinate does not have its standard
!-- location 0.0, meaning that the respective boundary is "open".
    topo          = IBCLR( topo, 0 )
    IF ( nest_shift_z == 0.0_wp )  topo(nzb,:,:) = IBSET( topo(nzb,:,:), 0 )
    topo(nzb,:,:) = IBSET( topo(nzb,:,:), 11 )
    topo(nzb,:,:) = IBSET( topo(nzb,:,:), 12 )
!
!-- In order to map topography on PALM grid also in case of ocean simulations, pre-calculate an
!-- offset value.
    ocean_offset = MERGE( zw(0), 0.0_wp, ocean_mode )
!
!-- Reference buildings on top of orography. This is not necessary if topography is read from ASCII
!-- file as no distinction between buildings and terrain height can be made. Moreover, this is also
!-- not necessary if urban-surface and land-surface model are used at the same time.
    IF ( input_pids_static )  THEN

       IF ( buildings_f%from_file )  THEN
          num_buildings_l = 0
          num_buildings   = 0
!
!--       Allocate at least one element for building ids and give it an inital negative value that
!--       will be overwritten later. This, however, is necessary in case there all IDs in the model
!--       domain are fill values.
          ALLOCATE( build_ids_l(1) )
          build_ids_l = -1
          DO  i = nxl, nxr
             DO  j = nys, nyn
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
                   IF ( num_buildings_l(myid) > 0 )  THEN
                      IF ( ANY( building_id_f%var(j,i) ==  build_ids_l ) )  THEN
                         CYCLE
                      ELSE
                         num_buildings_l(myid) = num_buildings_l(myid) + 1
!
!--                      Resize array with different local building ids
                         ALLOCATE( build_ids_l_tmp(1:SIZE(build_ids_l)) )
                         build_ids_l_tmp = build_ids_l
                         DEALLOCATE( build_ids_l )
                         ALLOCATE( build_ids_l(1:num_buildings_l(myid)) )
                         build_ids_l(1:num_buildings_l(myid)-1) =                                  &
                                                          build_ids_l_tmp(1:num_buildings_l(myid)-1)
                         build_ids_l(num_buildings_l(myid)) = building_id_f%var(j,i)
                         DEALLOCATE( build_ids_l_tmp )
                      ENDIF
!
!--                First occuring building id on PE.
                   ELSE
                      num_buildings_l(myid) = num_buildings_l(myid) + 1
                      build_ids_l(1) = building_id_f%var(j,i)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
!
!--       Determine number of different building ids for the entire domain.
#if defined( __parallel )
          CALL MPI_ALLREDUCE( num_buildings_l, num_buildings, numprocs, MPI_INTEGER, MPI_SUM,      &
                              comm2d, ierr )
#else
          num_buildings = num_buildings_l
#endif
!
!--       Gather all buildings ids on each PEs.
!--       First, allocate array encompassing all building ids in model domain.
          ALLOCATE( build_ids(1:SUM( num_buildings )) )
#if defined( __parallel )
!
!--       Allocate array for displacements.
!--       As each PE may has a different number of buildings, so that the block sizes send by each
!--       PE may not be equal. Hence,  information about the respective displacement is required,
!--       indicating the respective adress where each MPI-task writes into the receive buffer array.
          ALLOCATE( displace_dum(0:numprocs-1) )
          displace_dum(0) = 0
          DO i = 1, numprocs-1
             displace_dum(i) = displace_dum(i-1) + num_buildings(i-1)
          ENDDO

          CALL MPI_ALLGATHERV( build_ids_l(1:num_buildings_l(myid)), num_buildings(myid),          &
                               MPI_INTEGER, build_ids, num_buildings, displace_dum, MPI_INTEGER,   &
                               comm2d, ierr )
          DEALLOCATE( displace_dum )

#else
          build_ids = build_ids_l
#endif
!
!--       Note, in parallel mode building ids can occure mutliple times, as each PE has send its own
!--       ids. Therefore, sort out building ids which appear more than one time.
          num_build = 0
          DO  nr = 1, SIZE( build_ids )

             IF ( ALLOCATED( build_ids_final ) )  THEN
                IF ( ANY( build_ids(nr) == build_ids_final ) )  THEN
                   CYCLE
                ELSE
                   num_build = num_build + 1
!
!--                Resize.
                   ALLOCATE( build_ids_final_tmp(1:num_build) )
                   build_ids_final_tmp(1:num_build-1) = build_ids_final(1:num_build-1)
                   DEALLOCATE( build_ids_final )
                   ALLOCATE( build_ids_final(1:num_build) )
                   build_ids_final(1:num_build-1) = build_ids_final_tmp(1:num_build-1)
                   build_ids_final(num_build) = build_ids(nr)
                   DEALLOCATE( build_ids_final_tmp )
                ENDIF
             ELSE
                num_build = num_build + 1
                ALLOCATE( build_ids_final(1:num_build) )
                build_ids_final(num_build) = build_ids(nr)
             ENDIF
          ENDDO
!
!--       Determine maximumum terrain height occupied by the respective building and temporalily
!--       store on oro_max. Before, check whether any buildings are defined within the domain.
          IF ( ALLOCATED( build_ids_final ) )  THEN
             dim_builds = SIZE(build_ids_final)
          ELSE
             dim_builds = 0
          ENDIF

          ALLOCATE( oro_max_l(1:dim_builds) )
          ALLOCATE( oro_max(1:dim_builds) )
          oro_max_l = 0.0_wp

          DO  nr = 1, dim_builds
             oro_max_l(nr) = MAXVAL( MERGE( terrain_height_f%var(nys:nyn,nxl:nxr), 0.0_wp,         &
                                            building_id_f%var(nys:nyn,nxl:nxr) ==                  &
                                            build_ids_final(nr) ) )
          ENDDO

#if defined( __parallel )
          IF ( dim_builds >= 1 )  THEN
             CALL MPI_ALLREDUCE( oro_max_l, oro_max, SIZE( oro_max ), MPI_REAL, MPI_MAX, comm2d,   &
                                 ierr )
          ENDIF
#else
          oro_max = oro_max_l
#endif
!
!--       Finally, determine discrete grid height of maximum orography occupied by a building. Use
!--       all-or-nothing approach, i.e. if terrain exceeds the scalar level the grid box is fully
!--       terrain and the maximum terrain is set to the zw level.
          oro_max_l = 0.0
          DO  nr = 1, dim_builds
             DO  k = nzb, nzt
                IF ( zu(k) - ocean_offset <= oro_max(nr) )  oro_max_l(nr) = zw(k) - ocean_offset
             ENDDO
             oro_max(nr) = oro_max_l(nr)
          ENDDO
       ENDIF
!
!--    Allocate array for storing terrain height under buildings.
       IF ( buildings_f%from_file )  THEN
          ALLOCATE( buildings_f%oro_max(nysg:nyng,nxlg:nxrg) )
          buildings_f%oro_max = buildings_f%fill1
       ENDIF
!
!--    Map orography as well as buildings onto grid.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             topo_top_index = 0
!
!--          Obtain index in global building_id array.
             IF ( buildings_f%from_file )  THEN
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
!
!--                Determine index where maximum terrain height occupied by the respective building
!--                height is stored.
                   nr = MINLOC( ABS( build_ids_final - building_id_f%var(j,i) ), DIM=1 )
!
!--                Save grid-indexed oro_max.
                   buildings_f%oro_max(j,i) = oro_max(nr)
                ENDIF
             ENDIF
             DO  k = nzb, nzt
!
!--             In a first step, if grid point is below or equal the given terrain height, grid
!--             point is flagged to be of type natural.
!--             Please note, in case there is also a building which is lower than the vertical grid
!--             spacing, initialization of surface attributes will not be correct as given surface
!--             information will not be in accordance to the classified grid points.
!--             Hence, in this case, also a building flag.
                IF ( zu(k) - ocean_offset <= terrain_height_f%var(j,i) )  THEN
                   topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                   topo(k,j,i) = IBSET( topo(k,j,i), 1 )
                   topo_top_index = k ! topo_top_index + 1
                ENDIF
!
!--             Set building grid points. Here, only consider 2D buildings.
!--             3D buildings require separate treatment.
                IF ( buildings_f%from_file  .AND.  buildings_f%lod == 1 )  THEN
!
!--                Fill-up the terrain to the level of maximum orography within the building-covered
!--                area.
                   IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
!
!--                   Note, oro_max is always on zw level.
                      IF ( zu(k) - ocean_offset < oro_max(nr) )  THEN
                         topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         topo(k,j,i) = IBSET( topo(k,j,i), 1 )
                      ELSEIF ( zu(k) - ocean_offset <= oro_max(nr) + buildings_f%var_2d(j,i) )  THEN
                         topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                         topo(k,j,i) = IBSET( topo(k,j,i), 2 )
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
!
!--          Special treatment for non grid-resolved buildings. This case, the uppermost terrain
!--          grid point is flagged as building as well, even though no building exists at all.
!--          However, the surface element will be identified as urban-surface and the input data
!--          provided by the drivers is consistent to the surface classification. Else, all non
!--          grid-resolved buildings would vanish and identified as terrain grid points, which,
!--          however, won't be consistent with the input data.
             IF ( buildings_f%from_file  .AND.  buildings_f%lod == 1 )  THEN
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
                   DO  k = nzb, nzt
                      IF( zw(k) - ocean_offset == oro_max(nr) )  THEN
                         IF ( buildings_f%var_2d(j,i) <= zu(k+1) - zw(k) )  THEN
                            topo(k,j,i) = IBSET( topo(k,j,i), 2 )
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF
!
!--          Map 3D buildings onto terrain height.
!--          In case of any slopes, map building on top of maximum terrain height covered by the
!--          building. In other words, extend building down to the respective local terrain-surface
!--          height.
             IF ( buildings_f%from_file  .AND.  buildings_f%lod == 2 )  THEN
                IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
!
!--                Extend building down to the terrain surface, i.e. fill-up surface irregularities
!--                below a building. Note, oro_max is already a discrete height according to the
!--                all-or-nothing approach, i.e. grid box is either topography or atmosphere,
!--                terrain top is defined at upper bound of the grid box.
!--                Hence, check for zw in this case.
!--                Note, do this only for buildings which are surface mounted, i.e. building types
!--                1-6. Below bridges, which are represented exclusively by building type 7, terrain
!--                shape should be maintained.
                   IF ( building_type_f%from_file )  THEN
                      IF ( building_type_f%var(j,i) /= 7 )  THEN
                         DO k = topo_top_index + 1, nzt + 1
                            IF ( zu(k) - ocean_offset <= oro_max(nr) )  THEN
                               topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                               topo(k,j,i) = IBSET( topo(k,j,i), 1 )
                            ENDIF
                         ENDDO
!
!--                      After surface irregularities are smoothen, determine lower start index
!--                      where building starts.
                         DO  k = nzb, nzt
                            IF ( zu(k) - ocean_offset <= oro_max(nr) )  topo_top_index = k
                         ENDDO
                      ENDIF
                   ENDIF
!
!--                Finally, map building on top.
                   k2 = 0
                   DO k = topo_top_index, nzt + 1
                      IF ( k2 <= buildings_f%nz-1 )  THEN
                         IF ( buildings_f%var_3d(k2,j,i) == 1 )  THEN
                            topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                            topo(k,j,i) = IBSET( topo(k,j,i), 2 )
                         ENDIF
                      ENDIF
                      k2 = k2 + 1
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDDO
!
!--    Horizontal exchange the oro_max array, which is required to for initialization of
!--    building-surface properties.
       IF ( ALLOCATED( buildings_f%oro_max ) )  THEN
          CALL exchange_horiz_2d( buildings_f%oro_max )
          CALL set_lateral_neumann_bc( buildings_f%oro_max )
       ENDIF
!
!--    Deallocate temporary arrays required for processing and reading data
       IF ( ALLOCATED( oro_max         ) )  DEALLOCATE( oro_max )
       IF ( ALLOCATED( oro_max_l       ) )  DEALLOCATE( oro_max_l )
       IF ( ALLOCATED( build_ids_final ) )  DEALLOCATE( build_ids_final )
!
!-- Topography input via ASCII format.
    ELSE
       ocean_offset = MERGE( zw(0), 0.0_wp, ocean_mode )
!
!--    Initialize topography bit 0 (indicates obstacle) everywhere to zero and clear all grid points
!--    at nzb, where always a surface is defined.
!--    Further, set also bit 1 (indicates terrain) at nzb, which is further used for masked data
!--    output and further processing. Note, in the ASCII case no distinction is made between
!--    buildings and terrain, so that setting of bit 1 and 2 at the same time has no effect.
       topo          = IBCLR( topo, 0 )
       topo(nzb,:,:) = IBSET( topo(nzb,:,:), 0 )
       topo(nzb,:,:) = IBSET( topo(nzb,:,:), 1 )
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb, nzt
!
!--             Flag topography for all grid points which are below the local topography height.
!--             Note, each topography is flagged as building (bit 2) as well as terrain (bit 1) in
!--             order to employ urban-surface as well as land-surface model.
                IF ( ( zu(k) - ocean_offset ) <= buildings_f%var_2d(j,i) )  THEN
                   topo(k,j,i) = IBSET( topo(k,j,i), 0 )
                   topo(k,j,i) = IBSET( topo(k,j,i), 1 )
                   topo(k,j,i) = IBSET( topo(k,j,i), 2 )
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Exchange ghost points. Further, in case of non-cyclic boundary conditions Neumann BC
!-- are set for the topography.
    CALL exchange_horiz_int( topo, nys, nyn, nxl, nxr, nzt, nbgp )
    CALL topography_set_non_cyc_bc( topo )

 END SUBROUTINE process_topography


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Setting of further topography inidices on the grid.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE topography_set_indices

    INTEGER(iwp) ::  i             !< grid index in x-direction
    INTEGER(iwp) ::  j             !< grid index in y-direction
    INTEGER(iwp) ::  k             !< grid index in z-direction
    INTEGER(iwp) ::  k_top         !< topography top index on local PE
    INTEGER(iwp) ::  nzb_local_max !< vertical grid index of maximum topography height
    INTEGER(iwp) ::  nzb_local_min !< vertical grid index of minimum topography height


!
!-- Determine the maximum level of topography. It is used for steering the degradation of order of
!-- the applied advection scheme, as well in the lpm.
    k_top = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt + 1
             k_top = MAX( k_top, MERGE( k, 0, BTEST( topo(k,j,i), 0 ) ) )
          ENDDO
       ENDDO
    ENDDO
#if defined( __parallel )
    CALL MPI_ALLREDUCE( k_top, nzb_max, 1, MPI_INTEGER, MPI_MAX, comm2d, ierr )
#else
    nzb_max = k_top
#endif
!
!-- Increment nzb_max by 1 in order to allow for proper diverengence correction.
!-- Further, in case topography extents up to the model top, limit to nzt.
    nzb_max = MIN( nzb_max+1, nzt )
!
!-- Determine minimum index of topography. Usually, this will be nzb. In case there is elevated
!-- topography, however, the lowest topography will be higher.
!-- This index is e.g. used to calculate mean first-grid point atmosphere temperature, surface
!-- pressure and density, etc. .
    topo_min_level   = 0
#if defined( __parallel )
    CALL MPI_ALLREDUCE( MINVAL( topo_top_ind(nys:nyn,nxl:nxr,0) ), topo_min_level, 1, MPI_INTEGER, &
                        MPI_MIN, comm2d, ierr )
#else
    topo_min_level = MINVAL( topo_top_ind(nys:nyn,nxl:nxr,0) )
#endif

!
!-- Check topography for consistency with model domain. Therefore, use maximum and minium
!-- topography-top indices. Note, minimum topography top index is already calculated.
    IF ( TRIM( topography ) /= 'flat' )  THEN
#if defined( __parallel )
       CALL MPI_ALLREDUCE( MAXVAL( topo_top_ind(nys:nyn,nxl:nxr,0) ), nzb_local_max, 1,            &
                           MPI_INTEGER, MPI_MAX, comm2d, ierr )
#else
       nzb_local_max = MAXVAL( topo_top_ind(nys:nyn,nxl:nxr,0) )
#endif
       nzb_local_min = topo_min_level
!
!--    Consistency checks
       IF ( nzb_local_min < 0  .OR.  nzb_local_max  > ( nz + 1 ) )  THEN
          WRITE( message_string, * ) 'nzb_local values are outside the model domain',              &
                                     '&MINVAL( nzb_local ) = ', nzb_local_min,                     &
                                     '&MAXVAL( nzb_local ) = ', nzb_local_max
          CALL message( 'topography_mod', 'PAC0346', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- Allocate and set the arrays containing the topography height (for output reasons only).
    IF ( TRIM( topography ) /= 'flat' )  THEN

       IF ( nxr == nx  .AND.  nyn /= ny )  THEN
          ALLOCATE( zu_s_inner(nxl:nxr+1,nys:nyn), zw_w_inner(nxl:nxr+1,nys:nyn) )
       ELSEIF ( nxr /= nx  .AND.  nyn == ny )  THEN
          ALLOCATE( zu_s_inner(nxl:nxr,nys:nyn+1), zw_w_inner(nxl:nxr,nys:nyn+1) )
       ELSEIF ( nxr == nx  .AND.  nyn == ny )  THEN
          ALLOCATE( zu_s_inner(nxl:nxr+1,nys:nyn+1), zw_w_inner(nxl:nxr+1,nys:nyn+1) )
       ELSE
          ALLOCATE( zu_s_inner(nxl:nxr,nys:nyn), zw_w_inner(nxl:nxr,nys:nyn) )
       ENDIF

       zu_s_inner = 0.0_wp
       zw_w_inner = 0.0_wp
!
!--    Determine local topography height on scalar and w-grid. Note, setting lateral boundary values
!--    is not necessary, realized via topo_flags array. Further, please note that loop
!--    bounds are different from nxl to nxr and nys to nyn on south and right model boundary, hence,
!--    use intrinsic lbound and ubound functions to infer array bounds.
       DO  i = LBOUND( zu_s_inner, 1 ), UBOUND( zu_s_inner, 1 )
          DO  j = LBOUND( zu_s_inner, 2 ), UBOUND( zu_s_inner, 2 )
!
!--          Topography height on scalar grid. Therefore, determine index of upward-facing surface
!--          element on scalar grid.
             zu_s_inner(i,j) = zu(topo_top_ind(j,i,0))
!
!--          Topography height on w grid. Therefore, determine index of upward-facing surface
!--          element on w grid.
             zw_w_inner(i,j) = zw(topo_top_ind(j,i,3))
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE topography_set_indices


!--------------------------------------------------------------------------------------------------!
! Description:
! -------------------------------------------------------------------------------------------------!
!> Routine to set boundary conditions which is repeatedly necessary during topography
!> initialization. Note, this routine is just to avoid lengthy code.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE topography_set_non_cyc_bc( topo_array )

    INTEGER(iwp) ::  i  !< running index over ghost layers

    INTEGER(iwp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  topo_array !< input array for 3D topography


!
!-- Set boundary conditions also for flags. Can be interpreted as Neumann boundary conditions
!-- for topography.
    IF ( .NOT. bc_ns_cyc )  THEN
       IF ( nys == 0  )  THEN
          DO  i = 1, nbgp
             topo_array(:,nys-i,:) = topo_array(:,nys,:)
          ENDDO
       ENDIF
       IF ( nyn == ny )  THEN
          DO  i = 1, nbgp
             topo_array(:,nyn+i,:) = topo_array(:,nyn,:)
          ENDDO
       ENDIF
    ENDIF
    IF ( .NOT. bc_lr_cyc )  THEN
       IF ( nxl == 0  )  THEN
          DO  i = 1, nbgp
             topo_array(:,:,nxl-i) = topo_array(:,:,nxl)
          ENDDO
       ENDIF
       IF ( nxr == nx )  THEN
          DO  i = 1, nbgp
             topo_array(:,:,nxr+i) = topo_array(:,:,nxr)
          ENDDO
       ENDIF
    ENDIF

 END SUBROUTINE topography_set_non_cyc_bc

 END MODULE topography_mod
