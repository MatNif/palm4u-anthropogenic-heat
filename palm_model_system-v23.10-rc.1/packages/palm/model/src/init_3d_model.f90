!> @file init_3d_model.f90
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
!> Allocation of arrays and initialization of the 3D model via
!> a) pre-run the 1D model
!> or
!> b) pre-set constant linear profiles
!> or
!> c) read values of a previous run
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_3d_model

#if defined( __parallel )
    USE MPI
#endif

    USE advec_ws

    USE arrays_3d

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  barometric_formula,                                                                 &
               c_p,                                                                                &
               exner_function,                                                                     &
               exner_function_invers,                                                              &
               g,                                                                                  &
               ideal_gas_law_rho,                                                                  &
               ideal_gas_law_rho_pt,                                                               &
               l_v,                                                                                &
               pi

    USE boundary_settings_mod,                                                                     &
        ONLY:  set_lateral_neumann_bc

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model

    USE control_parameters

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz_2d

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy,                                                                                 &
               ddx2_mg,                                                                            &
               ddy2_mg

    USE indices

    USE kinds

    USE lsf_nudging_mod,                                                                           &
        ONLY:  ls_forcing_surf

    USE model_1d_mod,                                                                              &
        ONLY:  init_1d_model,                                                                      &
               l1d,                                                                                &
               u1d,                                                                                &
               v1d

    USE module_interface,                                                                          &
        ONLY:  module_interface_check_data_output_ts,                                              &
               module_interface_init_before_pressure_solver,                                       &
               module_interface_init_after_pressure_solver,                                        &
               module_interface_init_arrays,                                                       &
               module_interface_init_checks

    USE multi_agent_system_mod,                                                                    &
        ONLY:  agents_active, mas_init

    USE netcdf_interface,                                                                          &
        ONLY:  dots_label,                                                                         &
               dots_max,                                                                           &
               dots_num,                                                                           &
               dots_unit


    USE netcdf_data_input_mod,                                                                     &
        ONLY:  add_ghost_layers,                                                                   &
               char_fill,                                                                          &
               check_existence,                                                                    &
               close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_variable,                                                                       &
               init_3d,                                                                            &
               input_pids_dynamic,                                                                 &
               input_pids_static,                                                                  &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               input_file_dynamic,                                                                 &
               input_file_static,                                                                  &
               netcdf_data_input_1d,                                                               &
               netcdf_data_input_init_3d,                                                          &
               num_var_pids,                                                                       &
               open_read_file,                                                                     &
               pids_id,                                                                            &
               real_1d_3d,                                                                         &
               real_2d,                                                                            &
               vars_pids

    USE nesting_offl_mod,                                                                          &
        ONLY:  nesting_offl_init_modules

    USE palm_date_time_mod,                                                                        &
        ONLY:  init_date_time

    USE pegrid

    USE pmc_interface,                                                                             &
        ONLY:  atmosphere_ocean_coupled_run,                                                       &
               nesting_bounds

    USE random_function_mod

    USE random_generator_parallel,                                                                 &
        ONLY:  init_parallel_random_generator

    USE read_restart_data_mod,                                                                     &
        ONLY:  rrd_global_spinup,                                                                  &
               rrd_local,                                                                          &
               rrd_local_spinup,                                                                   &
               rrd_read_parts_of_global

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               hom_sum,                                                                            &
               mean_surface_level_height,                                                          &
               pr_max,                                                                             &
               pr_palm,                                                                            &
               rmask,                                                                              &
               statistic_regions,                                                                  &
               sums,                                                                               &
               sums_divnew_l,                                                                      &
               sums_divold_l,                                                                      &
               sums_l,                                                                             &
               sums_l_l,                                                                           &
               sums_wsts_bc_l,                                                                     &
               ts_value,                                                                           &
               weight_pres,                                                                        &
               weight_substep

    USE surface_layer_fluxes_mod,                                                                  &
        ONLY:  init_surface_layer_fluxes

    USE surface_mod,                                                                               &
        ONLY:  init_single_surface_properties,                                                     &
               init_surface_arrays,                                                                &
               init_surfaces,                                                                      &
               surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_usm

    USE surface_data_output_mod,                                                                   &
        ONLY:  surface_data_output_init,                                                           &
               surface_data_output_init_arrays

    IMPLICIT NONE


    INTEGER(iwp) ::  i                    !< grid index in x direction
    INTEGER(iwp) ::  ind_array(1)         !< dummy used to determine start index for external pressure forcing
    INTEGER(iwp) ::  j                    !< grid index in y direction
    INTEGER(iwp) ::  k                    !< grid index in z direction
    INTEGER(iwp) ::  l                    !< running index over multigrid levels
    INTEGER(iwp) ::  nz_s_shift           !< topography-top index on scalar-grid, used to vertically shift initial profiles
    INTEGER(iwp) ::  nz_u_shift           !< topography-top index on u-grid, used to vertically shift initial profiles
    INTEGER(iwp) ::  nz_v_shift           !< topography-top index on v-grid, used to vertically shift initial profiles
    INTEGER(iwp) ::  nz_w_shift           !< topography-top index on w-grid, used to vertically shift initial profiles
    INTEGER(iwp) ::  nzt_l                !< index of top PE boundary for multigrid level
    INTEGER(iwp) ::  sr                   !< index of statistic region

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE   ::  ngp_2dh_l  !< toal number of horizontal grid points in statistical region on
                                                             !< subdomain

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  ngp_2dh_outer_l    !< number of horizontal non-wall bounded grid points on
                                                                     !< subdomain
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  ngp_2dh_s_inner_l  !< number of horizontal non-topography grid points on
                                                                     !< subdomain

    REAL(wp) ::  dx_l !< grid spacing along x on different multigrid level
    REAL(wp) ::  dy_l !< grid spacing along y on different multigrid level

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  init_l                       !< dummy array used for averaging 3D data to obtain
                                                                         !< inital profiles
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mean_surface_level_height_l  !< mean surface level height on subdomain
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ngp_3d_inner_l               !< total number of non-topography grid points on subdomain
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ngp_3d_inner_tmp             !< total number of non-topography grid points
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  p_hydrostatic                !< hydrostatic pressure

    REAL(wp), DIMENSION(1:3) ::  volume_flow_area_l     !< area of lateral and top model domain surface on local subdomain
    REAL(wp), DIMENSION(1:3) ::  volume_flow_initial_l  !< initial volume flow into model domain

    TYPE(real_1d_3d) ::  tmp_1d !< temporary variable to input additional data from static/dynamic file
    TYPE(real_2d)    ::  tmp_2d !< temporary variable to input additional surface-data from static file

    CALL location_message( 'model initialization', 'start' )
!
!-- Set reference date-time
    CALL init_date_time( date_time_str=origin_date_time,                                           &
                         use_fixed_date=use_fixed_date,                                            &
                         use_fixed_time=use_fixed_time )

    IF ( debug_output )  CALL debug_message( 'allocating arrays', 'start' )
!
!-- Allocate arrays
    ALLOCATE( sums_divnew_l(0:statistic_regions),                                                  &
              sums_divold_l(0:statistic_regions) )
    ALLOCATE( dp_smooth_factor(nzb:nzt), rdf(nzb+1:nzt), rdf_sc(nzb+1:nzt) )
    ALLOCATE( rmask(nysg:nyng,nxlg:nxrg,0:statistic_regions),                                      &
              sums(nzb:nzt+1,pr_max),                                                              &
              sums_l(nzb:nzt+1,pr_max,0:threads_per_task-1),                                       &
              sums_l_l(nzb:nzt+1,0:statistic_regions,0:threads_per_task-1),                        &
              sums_wsts_bc_l(nzb:nzt+1,0:statistic_regions) )
    ALLOCATE( ts_value(dots_max,0:statistic_regions) )
    ALLOCATE( odf_x(nxlg:nxrg), odf_y(nysg:nyng) )
    ALLOCATE( ptdf_x(nxlg:nxrg), ptdf_y(nysg:nyng) )

    ALLOCATE( d(nzb+1:nzt,nys:nyn,nxl:nxr),                                                        &
              p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                    &
              tend(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    ALLOCATE( pt_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                 &
              pt_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                 &
              u_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                  &
              u_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                  &
              u_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                  &
              v_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                  &
              v_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                  &
              v_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                  &
              w_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                  &
              w_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                                  &
              w_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    IF (  .NOT.  neutral )  THEN
       ALLOCATE( pt_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF
!
!-- Pre-set masks for regional statistics. Default is the total model domain.
!-- Ghost points are excluded because counting values at the ghost boundaries would bias the
!-- statistics.
    rmask = 1.0_wp
    rmask(:,nxlg:nxl-1,:) = 0.0_wp;  rmask(:,nxr+1:nxrg,:) = 0.0_wp
    rmask(nysg:nys-1,:,:) = 0.0_wp;  rmask(nyn+1:nyng,:,:) = 0.0_wp
!
!-- Following array is required for perturbation pressure within the iterative pressure solvers. For
!-- the multistep schemes (Runge-Kutta), array p holds the weighted average of the substeps and
!-- cannot be used in the Poisson solver.
    IF ( psolver == 'sor' )  THEN
       ALLOCATE( p_loc(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ELSEIF ( psolver(1:9) == 'multigrid' )  THEN
!
!--    For performance reasons, multigrid is using one ghost layer only
       ALLOCATE( p_loc(nzb:nzt+1,nys-1:nyn+1,nxl-1:nxr+1) )
    ENDIF

    IF ( humidity )  THEN
!
!--    3D-humidity
       ALLOCATE( q_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                               &
                 q_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                               &
                 q_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                               &
                 vpt_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF

    IF ( passive_scalar )  THEN

!
!--    3D scalar arrays
       ALLOCATE( s_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                               &
                 s_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                               &
                 s_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    ENDIF

!
!-- Allocate and set 1d-profiles for Stokes drift velocity. It may be set to non-zero values later
!-- in ocean_init.
    ALLOCATE( u_stokes_zu(nzb:nzt+1), u_stokes_zw(nzb:nzt+1),                                      &
              v_stokes_zu(nzb:nzt+1), v_stokes_zw(nzb:nzt+1) )
    u_stokes_zu(:) = 0.0_wp
    u_stokes_zw(:) = 0.0_wp
    v_stokes_zu(:) = 0.0_wp
    v_stokes_zw(:) = 0.0_wp

!
!-- Read surface pressure from dynamic file (if present).
    IF ( input_pids_dynamic  .AND.  INDEX( initializing_actions, 'read_from_file' ) /= 0 )  THEN
       CALL netcdf_data_input_1d( TRIM( input_file_dynamic ) // TRIM( coupling_char ),             &
                                  'surface_forcing_surface_pressure', 'time', tmp_1d%var1d )
!
!--    Assign value to surface_pressure variable. Be aware that surface_pressure is given in hPa,
!--    while the unit in the dynamic file is Pa. After this has been done, deallocate memory.
       surface_pressure = tmp_1d%var1d(0) * 0.01_wp
       DEALLOCATE( tmp_1d%var1d )
    ENDIF
!
!-- Allocation of anelastic and Boussinesq approximation specific arrays
    ALLOCATE( p_hydrostatic(nzb:nzt+1) )
    ALLOCATE( rho_air(nzb:nzt+1) )
    ALLOCATE( rho_air_zw(nzb:nzt+1) )
    ALLOCATE( drho_air(nzb:nzt+1) )
    ALLOCATE( drho_air_zw(nzb:nzt+1) )

!-- Density profile calculation for anelastic and Boussinesq approximation.
!-- In case of a Boussinesq approximation, a constant density is calculated mainly for output
!-- purposes. This density does not need to be considered in the model's system of equations.
    IF ( TRIM( approximation ) == 'anelastic' )  THEN

       DO  k = nzb, nzt+1
          p_hydrostatic(k) = barometric_formula( zu(k), pt_surface *                               &
                                                 exner_function( surface_pressure * 100.0_wp ),    &
                                                 surface_pressure * 100.0_wp )

          rho_air(k) = ideal_gas_law_rho_pt( p_hydrostatic(k), pt_init(k) )
       ENDDO

       DO  k = nzb, nzt
          rho_air_zw(k) = 0.5_wp * ( rho_air(k) + rho_air(k+1) )
       ENDDO
       rho_air_zw(nzt+1)  = rho_air_zw(nzt) + 2.0_wp * ( rho_air(nzt+1) - rho_air_zw(nzt)  )

    ELSE
!
!--    Boussinesq-Approximation: density is assumed constant. The actual value of density does not
!--    effect the fluid simulation. Only some output quantities as dynamic pressure, divergence,
!--    or flux output in dynamic units rely on it. We use the surface density as the reference
!--    value here.
       IF ( ocean_mode )  THEN
!
!--       Set density to water density near the ocean surface.
          rho_air(:) = 1027.62_wp
       ELSE
          p_hydrostatic(:) = barometric_formula( zu(nzb), pt_surface *                             &
                                                 exner_function( surface_pressure * 100.0_wp ),    &
                                                 surface_pressure * 100.0_wp )

          rho_air(:) = ideal_gas_law_rho_pt( p_hydrostatic(nzb), pt_init(nzb) )
       ENDIF
       rho_air_zw(:) = rho_air(:)

    ENDIF

!
!-- Compute the inverse density array in order to avoid expencive divisions
    drho_air    = 1.0_wp / rho_air
    drho_air_zw = 1.0_wp / rho_air_zw

!
!-- Allocation of flux conversion arrays
    ALLOCATE( heatflux_input_conversion(nzb:nzt+1)      )
    ALLOCATE( heatflux_output_conversion(nzb:nzt+1)     )
    ALLOCATE( momentumflux_input_conversion(nzb:nzt+1)  )
    ALLOCATE( momentumflux_output_conversion(nzb:nzt+1) )
    ALLOCATE( scalarflux_input_conversion(nzb:nzt+1)    )
    ALLOCATE( scalarflux_output_conversion(nzb:nzt+1)   )
    ALLOCATE( waterflux_input_conversion(nzb:nzt+1)     )
    ALLOCATE( waterflux_output_conversion(nzb:nzt+1)    )


!
!-- Calculate flux conversion factors according to approximation and in-/output mode
    DO  k = nzb, nzt+1

        IF ( TRIM( flux_input_mode ) == 'kinematic' )  THEN
            heatflux_input_conversion(k)      = rho_air_zw(k)
            momentumflux_input_conversion(k)  = rho_air_zw(k)
            scalarflux_input_conversion(k)    = rho_air_zw(k)
            waterflux_input_conversion(k)     = rho_air_zw(k)
        ELSEIF ( TRIM( flux_input_mode ) == 'dynamic' ) THEN
            heatflux_input_conversion(k)      = 1.0_wp / c_p
            momentumflux_input_conversion(k)  = 1.0_wp
            scalarflux_input_conversion(k)    = 1.0_wp
            waterflux_input_conversion(k)     = 1.0_wp / l_v
        ENDIF

        IF ( TRIM( flux_output_mode ) == 'kinematic' )  THEN
            heatflux_output_conversion(k)     = drho_air_zw(k)
            momentumflux_output_conversion(k) = drho_air_zw(k)
            scalarflux_output_conversion(k)   = drho_air_zw(k)
            waterflux_output_conversion(k)    = drho_air_zw(k)
        ELSEIF ( TRIM( flux_output_mode ) == 'dynamic' ) THEN
            heatflux_output_conversion(k)     = c_p
            momentumflux_output_conversion(k) = 1.0_wp
            scalarflux_output_conversion(k)   = 1.0_wp
            waterflux_output_conversion(k)    = l_v
        ENDIF

        IF ( .NOT. humidity ) THEN
            waterflux_input_conversion(k)  = 1.0_wp
            waterflux_output_conversion(k) = 1.0_wp
        ENDIF

    ENDDO

!
!-- In case of multigrid method, compute grid lengths and grid factors for the grid levels with
!-- respective density on each grid.
    IF ( psolver(1:9) == 'multigrid' )  THEN

       ALLOCATE( ddx2_mg(maximum_grid_level) )
       ALLOCATE( ddy2_mg(maximum_grid_level) )
       ALLOCATE( dzu_mg(nzb+1:nzt+1,maximum_grid_level) )
       ALLOCATE( dzw_mg(nzb+1:nzt+1,maximum_grid_level) )
       ALLOCATE( f1_mg(nzb+1:nzt,maximum_grid_level) )
       ALLOCATE( f2_mg(nzb+1:nzt,maximum_grid_level) )
       ALLOCATE( f3_mg(nzb+1:nzt,maximum_grid_level) )
       ALLOCATE( rho_air_mg(nzb:nzt+1,maximum_grid_level) )
       ALLOCATE( rho_air_zw_mg(nzb:nzt+1,maximum_grid_level) )

       dzu_mg(:,maximum_grid_level) = dzu
       rho_air_mg(:,maximum_grid_level) = rho_air
!
!--    Next line to ensure an equally spaced grid.
       dzu_mg(1,maximum_grid_level) = dzu(2)
       rho_air_mg(nzb,maximum_grid_level) = rho_air(nzb) + (rho_air(nzb) - rho_air(nzb+1))

       dzw_mg(:,maximum_grid_level) = dzw
       rho_air_zw_mg(:,maximum_grid_level) = rho_air_zw
       nzt_l = nzt
       DO  l = maximum_grid_level-1, 1, -1
           dzu_mg(nzb+1,l) = 2.0_wp * dzu_mg(nzb+1,l+1)
           dzw_mg(nzb+1,l) = 2.0_wp * dzw_mg(nzb+1,l+1)
           rho_air_mg(nzb,l)    = rho_air_mg(nzb,l+1)    + ( rho_air_mg(nzb,l+1)    -              &
                                                             rho_air_mg(nzb+1,l+1)    )
           rho_air_zw_mg(nzb,l) = rho_air_zw_mg(nzb,l+1) + ( rho_air_zw_mg(nzb,l+1) -              &
                                                             rho_air_zw_mg(nzb+1,l+1) )
           rho_air_mg(nzb+1,l)    = rho_air_mg(nzb+1,l+1)
           rho_air_zw_mg(nzb+1,l) = rho_air_zw_mg(nzb+1,l+1)
           nzt_l = nzt_l / 2
           DO  k = 2, nzt_l+1
              dzu_mg(k,l) = dzu_mg(2*k-2,l+1) + dzu_mg(2*k-1,l+1)
              dzw_mg(k,l) = dzw_mg(2*k-2,l+1) + dzw_mg(2*k-1,l+1)
              rho_air_mg(k,l)    = rho_air_mg(2*k-1,l+1)
              rho_air_zw_mg(k,l) = rho_air_zw_mg(2*k-1,l+1)
           ENDDO
       ENDDO

       nzt_l = nzt
       dx_l  = dx
       dy_l  = dy
       DO  l = maximum_grid_level, 1, -1
          ddx2_mg(l) = 1.0_wp / dx_l**2
          ddy2_mg(l) = 1.0_wp / dy_l**2
          DO  k = nzb+1, nzt_l
             f2_mg(k,l) = rho_air_zw_mg(k,l) / ( dzu_mg(k+1,l) * dzw_mg(k,l) )
             f3_mg(k,l) = rho_air_zw_mg(k-1,l) / ( dzu_mg(k,l)   * dzw_mg(k,l) )
             f1_mg(k,l) = 2.0_wp * ( ddx2_mg(l) + ddy2_mg(l) )                                     &
                          * rho_air_mg(k,l) + f2_mg(k,l) + f3_mg(k,l)
          ENDDO
          nzt_l = nzt_l / 2
          dx_l  = dx_l * 2.0_wp
          dy_l  = dy_l * 2.0_wp
       ENDDO

    ENDIF

!
!-- 1D-array for large scale subsidence velocity
    IF ( .NOT. ALLOCATED( w_subs ) )  THEN
       ALLOCATE ( w_subs(nzb:nzt+1) )
       w_subs = 0.0_wp
    ENDIF

!
!-- Initial assignment of the pointers
    IF ( .NOT. neutral )  THEN
       pt => pt_1;  pt_p => pt_2;  tpt_m => pt_3
    ELSE
       pt => pt_1;  pt_p => pt_1;  tpt_m => pt_3
    ENDIF
    u  => u_1;   u_p  => u_2;   tu_m  => u_3
    v  => v_1;   v_p  => v_2;   tv_m  => v_3
    w  => w_1;   w_p  => w_2;   tw_m  => w_3

    IF ( humidity )  THEN
       q => q_1;  q_p => q_2;  tq_m => q_3
       vpt  => vpt_1
    ENDIF

    IF ( passive_scalar )  THEN
       s => s_1;  s_p => s_2;  ts_m => s_3
    ENDIF

!
!-- Initialize potential temperature in case of neutral runs. Although in such a case the prognostic
!-- equation for temperature isn't calculated, a non-initialized pt might cause aborts at other
!-- locations, e.g. in TKE production terms, where the temperature gradient is used.
    IF ( neutral )  pt = pt_surface
!
!-- Initialize surface arrays
    CALL init_surface_arrays
!
!-- Allocate arrays for surface data output
    IF ( surface_output )  THEN
       CALL surface_data_output_init_arrays
    ENDIF
!
!-- Allocate arrays for other modules
    CALL module_interface_init_arrays


!
!-- Allocate arrays containing the RK coefficient for calculation of perturbation pressure and
!-- turbulent fluxes. At this point values are set for pressure calculation during initialization
!-- (where no timestep is done). Further below the values needed within the timestep scheme will be
!-- set.
    ALLOCATE( weight_substep(1:intermediate_timestep_count_max),                                   &
              weight_pres(1:intermediate_timestep_count_max) )
    weight_substep = 1.0_wp
    weight_pres    = 1.0_wp
    intermediate_timestep_count = 0  ! needed when simulated_time = 0.0

    IF ( debug_output )  CALL debug_message( 'allocating arrays', 'end' )

!
!-- Initialize time series
    ts_value = 0.0_wp

!
!-- Initialize local summation arrays for routine flow_statistics.
!-- This is necessary because they may not yet have been initialized when they are called from
!-- flow_statistics (or - depending on the chosen model run - are never initialized)
    sums_divnew_l  = 0.0_wp
    sums_divold_l  = 0.0_wp
    sums_l_l       = 0.0_wp
    sums_wsts_bc_l = 0.0_wp
!
!-- Initialize model variables.
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.                                &
         .NOT. cyclic_fill_initialization )                                                        &
    THEN
!
!--    Initialization with provided input data derived from larger-scale model.
       IF ( INDEX( initializing_actions, 'read_from_file' ) /= 0 )  THEN
          IF ( debug_output )  CALL debug_message( 'initializing with external data', 'start' )
!
!--       Read initial 1D profiles or 3D data from NetCDF file, depending on the provided
!--       level-of-detail.
!--       At the moment, only u, v, w, pt and q are provided.
          CALL netcdf_data_input_init_3d
!
!--       Please note, data from dynamic input file is from nzb+1 to nzt.
!--       Bottom and top boundary conditions for profiles are already set (just after
!--       reading), so that this is not necessary here.
!--       Depending on the provided level-of-detail, initial data is either stored on data
!--       type (lod=1), or directly on 3D arrays (lod=2).
!--       In order to obtain also initial profiles in case of lod=2 (which is required for e.g.
!--       damping), average over 3D data.
          IF( init_3d%lod_u == 1 )  THEN
             u_init = init_3d%u_init
          ELSEIF( init_3d%lod_u == 2 )  THEN
             ALLOCATE( init_l(nzb:nzt+1) )
             DO  k = nzb, nzt+1
                init_l(k) = SUM( u(k,nys:nyn,nxl:nxr) )
             ENDDO
             init_l = init_l / REAL( ( nx + 1 ) * ( ny + 1 ), KIND = wp )

#if defined( __parallel )
             CALL MPI_ALLREDUCE( init_l, u_init, nzt+1-nzb+1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
             u_init = init_l
#endif
             DEALLOCATE( init_l )

          ENDIF

          IF( init_3d%lod_v == 1 )  THEN
             v_init = init_3d%v_init
          ELSEIF( init_3d%lod_v == 2 )  THEN
             ALLOCATE( init_l(nzb:nzt+1) )
             DO  k = nzb, nzt+1
                init_l(k) = SUM( v(k,nys:nyn,nxl:nxr) )
             ENDDO
             init_l = init_l / REAL( ( nx + 1 ) * ( ny + 1 ), KIND = wp )

#if defined( __parallel )
             CALL MPI_ALLREDUCE( init_l, v_init, nzt+1-nzb+1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
             v_init = init_l
#endif
             DEALLOCATE( init_l )
          ENDIF
          IF( .NOT. neutral )  THEN
             IF( init_3d%lod_pt == 1 )  THEN
                pt_init = init_3d%pt_init
             ELSEIF( init_3d%lod_pt == 2 )  THEN
                ALLOCATE( init_l(nzb:nzt+1) )
                DO  k = nzb, nzt+1
                   init_l(k) = SUM( pt(k,nys:nyn,nxl:nxr) )
                ENDDO
                init_l = init_l / REAL( ( nx + 1 ) * ( ny + 1 ), KIND = wp )

#if defined( __parallel )
                CALL MPI_ALLREDUCE( init_l, pt_init, nzt+1-nzb+1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
                pt_init = init_l
#endif
                DEALLOCATE( init_l )
             ENDIF
          ENDIF


          IF( humidity )  THEN
             IF( init_3d%lod_q == 1 )  THEN
                q_init = init_3d%q_init
             ELSEIF( init_3d%lod_q == 2 )  THEN
                ALLOCATE( init_l(nzb:nzt+1) )
                DO  k = nzb, nzt+1
                   init_l(k) = SUM( q(k,nys:nyn,nxl:nxr) )
                ENDDO
                init_l = init_l / REAL( ( nx + 1 ) * ( ny + 1 ), KIND = wp )

#if defined( __parallel )
                CALL MPI_ALLREDUCE( init_l, q_init, nzt+1-nzb+1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
                q_init = init_l
#endif
                DEALLOCATE( init_l )
             ENDIF
          ENDIF

!
!--       Write initial profiles onto 3D arrays.
!--       Work-around, 3D initialization of u,v,w creates artificial structures which correlate with
!--       the processor grid. The reason for this is still unknown. To work-around this, 3D
!--       initialization will be effectively reduce to a 1D initialization where no such artificial
!--       structures appear.
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                IF( init_3d%lod_u == 1  .OR.  init_3d%lod_u == 2 )  u(:,j,i) = u_init(:)
                IF( init_3d%lod_v == 1  .OR.  init_3d%lod_u == 2 )  v(:,j,i) = v_init(:)
                IF( .NOT. neutral  .AND.  ( init_3d%lod_pt == 1  .OR.  init_3d%lod_pt == 2 ) )     &
                   pt(:,j,i) = pt_init(:)
                IF( humidity  .AND.  ( init_3d%lod_q == 1  .OR.  init_3d%lod_q == 2 ) )            &
                   q(:,j,i) = q_init(:)
             ENDDO
          ENDDO
!
!--       Set geostrophic wind components.
          IF ( init_3d%from_file_ug )  THEN
             ug(:) = init_3d%ug_init(:)
          ENDIF
          IF ( init_3d%from_file_vg )  THEN
             vg(:) = init_3d%vg_init(:)
          ENDIF
!
!--       Set bottom and top boundary condition for geostrophic wind
          ug(nzt+1) = ug(nzt)
          vg(nzt+1) = vg(nzt)
          ug(nzb)   = ug(nzb+1)
          vg(nzb)   = vg(nzb+1)
!
!--       Set inital w to 0.
          w = 0.0_wp

          IF ( passive_scalar )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   s(:,j,i) = s_init
                ENDDO
             ENDDO
          ENDIF

!
!--       Set velocity components at non-atmospheric / oceanic grid points to zero.
          u = MERGE( u, 0.0_wp, BTEST( topo_flags, 1 ) )
          v = MERGE( v, 0.0_wp, BTEST( topo_flags, 2 ) )
          w = MERGE( w, 0.0_wp, BTEST( topo_flags, 3 ) )
!
!--       Initialize surface variables, e.g. friction velocity, momentum fluxes, etc.
          CALL  init_surfaces

          IF ( debug_output )  CALL  debug_message( 'initializing with external data', 'end' )
!
!--    Initialization via computed 1D-model profiles.
       ELSEIF ( INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN

          IF ( debug_output )  CALL  debug_message( 'initializing with 1D model profiles', 'start' )
!
!--       Use solutions of the 1D model as initial profiles.
!--       Start 1D model.
          CALL init_1d_model
!
!--       Transfer initial profiles to the arrays of the 3D model.
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                pt(:,j,i) = pt_init
                u(:,j,i)  = u1d
                v(:,j,i)  = v1d
             ENDDO
          ENDDO

          IF ( humidity )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   q(:,j,i) = q_init
                ENDDO
             ENDDO
          ENDIF

          IF ( passive_scalar )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   s(:,j,i) = s_init
                ENDDO
             ENDDO
          ENDIF
!
!--          Store initial profiles for output purposes etc.
          IF ( .NOT. constant_diffusion )  THEN
             hom(:,1,25,:) = SPREAD( l1d, 2, statistic_regions+1 )
          ENDIF
!
!--       Set velocities back to zero.
          u = MERGE( u, 0.0_wp, BTEST( topo_flags, 1 ) )
          v = MERGE( v, 0.0_wp, BTEST( topo_flags, 2 ) )
!
!--       WARNING: The extra boundary conditions set after running the 1D model impose an error on
!--       -------- the divergence one layer below the topography; need to correct later
!--       ATTENTION: Provisional correction for Piacsek & Williams advection scheme: keep u and v
!--       ---------- zero one layer below the topography.
          IF ( ibc_uv_b == 1 )  THEN
!
!--          Neumann condition
             DO  i = nxl-1, nxr+1
                DO  j = nys-1, nyn+1
                   u(nzb,j,i) = u(nzb+1,j,i)
                   v(nzb,j,i) = v(nzb+1,j,i)
                ENDDO
             ENDDO

          ENDIF
!
!--       Initialize surface variables, e.g. friction velocity, momentum fluxes, etc.
          CALL init_surfaces

          IF ( debug_output )  CALL  debug_message( 'initializing with 1D model profiles', 'end' )

       ELSEIF ( ( INDEX(initializing_actions, 'set_constant_profiles') /= 0 )  .OR.                 &
                ( INDEX(initializing_actions, 'interpolate_from_parent') /= 0 ) )  THEN

          IF ( debug_output )  CALL  debug_message( 'initializing with constant profiles', 'start' )

!
!--       Use constructed initial profiles (velocity constant with height, temperature profile with
!--       constant gradient)
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                pt(:,j,i) = pt_init
                u(:,j,i)  = u_init
                v(:,j,i)  = v_init
             ENDDO
          ENDDO
!
!--       Mask topography
          u = MERGE( u, 0.0_wp, BTEST( topo_flags, 1 ) )
          v = MERGE( v, 0.0_wp, BTEST( topo_flags, 2 ) )
!
!--       Set initial horizontal velocities at the lowest computational grid levels to zero in order
!--       to avoid too small time steps caused by the diffusion limit in the initial phase of a run
!--       (at k=1, dz/2 occurs in the limiting formula!).
!--       Please note, in case land- or urban-surface model is used and a spinup is applied, masking
!--       the lowest computational level is not possible as MOST as well as energy-balance
!--       parametrizations will not work with zero wind velocity. For sake of comparison
          IF ( ibc_uv_b /= 1  .AND.  .NOT. spinup  .AND.  .NOT. read_spinup_data )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt
                      u(k,j,i) = MERGE( u(k,j,i), 0.0_wp, BTEST( topo_flags(k,j,i), 20 ) )
                      v(k,j,i) = MERGE( v(k,j,i), 0.0_wp, BTEST( topo_flags(k,j,i), 21 ) )
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          IF ( humidity )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   q(:,j,i) = q_init
                ENDDO
             ENDDO
          ENDIF

          IF ( passive_scalar )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   s(:,j,i) = s_init
                ENDDO
             ENDDO
          ENDIF

!
!--       Compute initial temperature field and other constants used in case of a sloping surface.
          IF ( sloping_surface )  CALL init_slope
!
!--       Initialize surface variables, e.g. friction velocity, momentum fluxes, etc.
          CALL init_surfaces

          IF ( debug_output )  CALL debug_message( 'initializing with constant profiles', 'end' )

       ELSEIF ( INDEX(initializing_actions, 'by_user') /= 0 )  THEN

          IF ( debug_output )  CALL debug_message( 'initializing by user', 'start' )
!
!--       Pre-initialize surface variables, i.e. setting start- and end-indices at each
!--       (j,i)-location. Please note, this does not supersede user-defined initialization of
!--       surface quantities.
          CALL init_surfaces
!
!--       Initialization will completely be done by the user
          CALL user_init_3d_model

          IF ( debug_output )  CALL debug_message( 'initializing by user', 'end' )

       ENDIF

       IF ( debug_output )  THEN
          CALL debug_message( 'initializing statistics, boundary conditions, etc.', 'start' )
       ENDIF

!
!--    Bottom boundary
       IF ( ibc_uv_b == 0  .OR.  ibc_uv_b == 2  )  THEN
          u(nzb,:,:) = 0.0_wp
          v(nzb,:,:) = 0.0_wp
       ENDIF

!
!--    Apply channel flow boundary condition
       IF ( TRIM( bc_uv_t ) == 'dirichlet_0' )  THEN
          u(nzt+1,:,:) = 0.0_wp
          v(nzt+1,:,:) = 0.0_wp
       ENDIF

!
!--    Calculate virtual potential temperature
       IF ( humidity )  vpt = pt * ( 1.0_wp + 0.61_wp * q )

!
!--    Store initial profiles for output purposes etc.. Please note, in case of initialization of u,
!--    v, w, pt, and q via output data derived from larger scale models, data will not be
!--    horizontally homogeneous. Actually, a mean profile should be calculated before.
       hom(:,1,5,:) = SPREAD( u(:,nys,nxl), 2, statistic_regions+1 )
       hom(:,1,6,:) = SPREAD( v(:,nys,nxl), 2, statistic_regions+1 )
       IF ( ibc_uv_b == 0 .OR. ibc_uv_b == 2)  THEN
          hom(nzb,1,5,:) = 0.0_wp
          hom(nzb,1,6,:) = 0.0_wp
       ENDIF
       hom(:,1,7,:)  = SPREAD( pt(:,nys,nxl), 2, statistic_regions+1 )

       IF ( humidity )  THEN
!
!--       Store initial profile of total water content, virtual potential temperature
          hom(:,1,26,:) = SPREAD(   q(:,nys,nxl), 2, statistic_regions+1 )
          hom(:,1,29,:) = SPREAD( vpt(:,nys,nxl), 2, statistic_regions+1 )
!
!--       Store initial profile of mixing ratio and potential temperature
          IF ( bulk_cloud_model  .OR.  cloud_droplets ) THEN
             hom(:,1,27,:) = SPREAD(  q(:,nys,nxl), 2, statistic_regions+1 )
             hom(:,1,28,:) = SPREAD( pt(:,nys,nxl), 2, statistic_regions+1 )
          ENDIF
       ENDIF

!
!--    Store initial scalar profile
       IF ( passive_scalar )  THEN
          hom(:,1,121,:) = SPREAD(  s(:,nys,nxl), 2, statistic_regions+1 )
       ENDIF

!
!--    Initialize the random number generators (from numerical recipes)
       CALL random_function_ini

       IF ( random_generator == 'random-parallel' )  THEN
          CALL init_parallel_random_generator( nx, nys, nyn, nxl, nxr )
       ENDIF
!
!--    Set the reference state to be used in the buoyancy terms (for ocean runs the reference state
!--    will be set (overwritten) in init_ocean).
       IF ( use_single_reference_value )  THEN
          IF ( .NOT. humidity )  THEN
             ref_state(:) = pt_reference
          ELSE
             ref_state(:) = vpt_reference
          ENDIF
       ELSE
          IF ( .NOT. humidity )  THEN
             ref_state(:) = pt_init(:)
          ELSE
             ref_state(:) = vpt(:,nys,nxl)
          ENDIF
       ENDIF

!
!--    For the moment, vertical velocity is zero
       w = 0.0_wp

!
!--    Initialize array sums (must be defined in first call of pres)
       sums = 0.0_wp

!
!--    In case of iterative solvers, p must get an initial value
       IF ( psolver(1:9) == 'multigrid'  .OR.  psolver == 'sor' )  p = 0.0_wp
!
!--    Impose vortex with vertical axis on the initial velocity profile
       IF ( INDEX( initializing_actions, 'initialize_vortex' ) /= 0 )  THEN
          CALL init_rankine
       ENDIF

!
!--    Impose temperature anomaly (advection test only) or warm air bubble close to surface.
       IF ( INDEX( initializing_actions, 'initialize_ptanom' ) /= 0  .OR.                          &
            INDEX( initializing_actions, 'initialize_bubble' ) /= 0  )  THEN
          CALL init_pt_anomaly
       ENDIF

!
!--    Initialize old and new time levels.
       tpt_m = 0.0_wp; tu_m = 0.0_wp; tv_m = 0.0_wp; tw_m = 0.0_wp
       pt_p = pt; u_p = u; v_p = v; w_p = w

       IF ( humidity  )  THEN
          tq_m = 0.0_wp
          q_p = q
       ENDIF

       IF ( passive_scalar )  THEN
          ts_m = 0.0_wp
          s_p  = s
       ENDIF

       IF ( debug_output )  THEN
          CALL debug_message( 'initializing statistics, boundary conditions, etc.', 'end' )
       ENDIF

    ELSEIF ( TRIM( initializing_actions ) == 'read_restart_data'  .OR.                             &
             cyclic_fill_initialization )                                                          &
    THEN

       IF ( debug_output )  THEN
          CALL debug_message( 'initialization in case of restart / cyclic_fill', 'start' )
       ENDIF
!
!--    Initialize surface elements and its attributes, e.g. heat- and momentumfluxes, roughness,
!--    scaling parameters. As number of surface elements might be different between runs, e.g. in
!--    case of cyclic fill, and not all surface elements are read, surface elements need to be
!--    initialized before.
!--    Please note, in case of cyclic fill, surfaces should be initialized after restart data is
!--    read, else, individual settings of surface parameters will be overwritten from data of
!--    precursor run, hence, init_surfaces is called a second time after reading the restart data.
       CALL init_surfaces
!
!--    When reading prerun data for cyclic fill, read some of the global variables from the restart
!--    data file which are required for initializing the inflow.
       IF ( cyclic_fill_initialization )  THEN

!
!--       Blockwise I/O does not work together with MPI-I/O
          IF ( restart_data_format_input(1:3) == 'mpi' )  THEN
             CALL rrd_read_parts_of_global
          ELSE
             DO  i = 0, io_blocks-1
                IF ( i == io_group )  THEN
                   CALL rrd_read_parts_of_global
                ENDIF
#if defined( __parallel )
                CALL MPI_BARRIER( comm2d, ierr )
#endif
             ENDDO
          ENDIF

       ENDIF

!
!--    Read processor specific binary data from restart file.
!--    Blockwise I/O does not work together with MPI-I/O
       IF ( restart_data_format_input(1:3) == 'mpi' )  THEN
          CALL rrd_local
       ELSE
          DO  i = 0, io_blocks-1
             IF ( i == io_group )  THEN
                CALL rrd_local
             ENDIF
#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
          ENDDO
       ENDIF


       IF ( cyclic_fill_initialization )  THEN

!
!--       In case of cyclic fill, call init_surfaces a second time, so that surface properties such
!--       as heat fluxes are initialized as prescribed.
          CALL init_surfaces

!
!--       Overwrite u_init, v_init, pt_init, q_init and s_init with the horizontally mean (hom)
!--       vertical profiles from the end of the prerun, because these profiles shall be used as the
!--       reference state for the rayleigh damping and the pt_damping. This is especially important
!--       for the use of large_scale_subsidence, because the reference temperature in the free
!--       atmosphere changes in time.
          u_init(:) = hom_sum(:,1,0)
          v_init(:) = hom_sum(:,2,0)
          pt_init(:) = hom_sum(:,4,0)
          IF ( humidity )  q_init(:) = hom_sum(:,41,0)
          IF ( passive_scalar )  s_init(:) = hom_sum(:,115,0)
       ENDIF
!
!--    In case of complex terrain and cyclic fill method as initialization, shift initial data in
!--    the vertical direction for each point in the x-y-plane depending on local surface height.
       IF ( terrain_following_mapping  .AND.  cyclic_fill_initialization )  THEN
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                nz_u_shift = topo_top_ind(j,i,1)
                nz_v_shift = topo_top_ind(j,i,2)
                nz_w_shift = topo_top_ind(j,i,3)
                nz_s_shift = topo_top_ind(j,i,0)

                u(nz_u_shift:nzt+1,j,i)  = u(0:nzt+1-nz_u_shift,j,i)

                v(nz_v_shift:nzt+1,j,i)  = v(0:nzt+1-nz_v_shift,j,i)

                w(nz_w_shift:nzt+1,j,i)  = w(0:nzt+1-nz_w_shift,j,i)

                p(nz_s_shift:nzt+1,j,i)  =  p(0:nzt+1-nz_s_shift,j,i)
                pt(nz_s_shift:nzt+1,j,i) = pt(0:nzt+1-nz_s_shift,j,i)
             ENDDO
          ENDDO
       ENDIF
!
!--    Inside buildings set velocities back to zero
       IF ( cyclic_fill_initialization  .AND.  topography /= 'flat' )  THEN
!
!--       Inside buildings set velocities back to zero.
!--       Other scalars (pt, q, s, p, sa, ...) are ignored at present,
!--       maybe revise later.
          DO  i = nxlg, nxrg
             DO  j = nysg, nyng
                DO  k = nzb, nzt
                   u(k,j,i) = MERGE( u(k,j,i), 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
                   v(k,j,i) = MERGE( v(k,j,i), 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                   w(k,j,i) = MERGE( w(k,j,i), 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
                ENDDO
             ENDDO
          ENDDO

       ENDIF
!
!--    Calculate initial temperature field and other constants used in case of a sloping surface.
       IF ( sloping_surface )  CALL init_slope
!
!--    Initialize new time levels (only done in order to set boundary values including ghost
!--    points).
       pt_p = pt; u_p = u; v_p = v; w_p = w
       IF ( humidity )  THEN
          q_p = q
       ENDIF
       IF ( passive_scalar )  s_p = s
!
!--    Allthough tendency arrays are set in prognostic_equations, they have have to be predefined
!--    here because they are used (but multiplied with 0) there before they are set.
       tpt_m = 0.0_wp; tu_m = 0.0_wp; tv_m = 0.0_wp; tw_m = 0.0_wp
       IF ( humidity )  THEN
          tq_m = 0.0_wp
       ENDIF
       IF ( passive_scalar )  ts_m  = 0.0_wp

       IF ( debug_output )  THEN
          CALL debug_message( 'initialization in case of restart / cyclic_fill', 'end' )
       ENDIF

    ELSE
!
!--    Actually this part of the programm should not be reached
       message_string = 'unknown initialization problem'
       CALL message( 'init_3d_model', 'PAC0208', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- If required, read surface spinup data from a previous run to initialize surfaces.
!-- Please note that the surface spinup data will overwrite the previously initialized surface
!-- Moreover, please note that these action needs to be done before the LSM and the USM are
!-- initialized (i.e. before module_interface_init is invoked), else also the prognostic
!-- time levels (_p variables) would need to be initialized.
    IF ( read_spinup_data )  THEN
       CALL location_message( 'Reading spinup data', 'start' )
       CALL rrd_global_spinup
       CALL rrd_local_spinup
       CALL location_message( 'Reading spinup data', 'finished' )
    ENDIF
!
!-- Calculate the initial volume flow at the right and north boundary
    IF ( conserve_volume_flow )  THEN

       IF ( use_prescribed_profile_data )  THEN

          volume_flow_initial_l = 0.0_wp
          volume_flow_area_l    = 0.0_wp

          IF ( nxr == nx )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   volume_flow_initial_l(1) = volume_flow_initial_l(1) +                           &
                                              u_init(k) * dzw(k)                                   &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,j,nxr), 1 )             &
                                                     )

                   volume_flow_area_l(1)    = volume_flow_area_l(1) + dzw(k)                       &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,j,nxr), 1 )             &
                                                     )
                ENDDO
             ENDDO
          ENDIF

          IF ( nyn == ny )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   volume_flow_initial_l(2) = volume_flow_initial_l(2) +                           &
                                              v_init(k) * dzw(k)                                   &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,nyn,i), 2 )             &
                                                     )
                   volume_flow_area_l(2)    = volume_flow_area_l(2) + dzw(k)                       &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,nyn,i), 2 )             &
                                                     )
                ENDDO
             ENDDO
          ENDIF

#if defined( __parallel )
          CALL MPI_ALLREDUCE( volume_flow_initial_l(1), volume_flow_initial(1), 2, MPI_REAL,       &
                              MPI_SUM, comm2d, ierr )
          CALL MPI_ALLREDUCE( volume_flow_area_l(1), volume_flow_area(1), 2, MPI_REAL, MPI_SUM,    &
                              comm2d, ierr )

#else
          volume_flow_initial = volume_flow_initial_l
          volume_flow_area    = volume_flow_area_l
#endif

       ELSEIF ( cyclic_fill_initialization )  THEN

          volume_flow_initial_l = 0.0_wp
          volume_flow_area_l    = 0.0_wp

          IF ( nxr == nx )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   volume_flow_initial_l(1) = volume_flow_initial_l(1) +                           &
                                              hom_sum(k,1,0) * dzw(k)                              &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,j,nx), 1 )              &
                                                     )
                   volume_flow_area_l(1)    = volume_flow_area_l(1) + dzw(k)                       &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,j,nx), 1 )              &
                                                     )
                ENDDO
             ENDDO
          ENDIF

          IF ( nyn == ny )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   volume_flow_initial_l(2) = volume_flow_initial_l(2) +                           &
                                              hom_sum(k,2,0) * dzw(k)                              &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,ny,i), 2 )              &
                                                     )
                   volume_flow_area_l(2)    = volume_flow_area_l(2) + dzw(k)                       &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,ny,i), 2 )              &
                                                     )
                ENDDO
             ENDDO
          ENDIF

#if defined( __parallel )
          CALL MPI_ALLREDUCE( volume_flow_initial_l(1), volume_flow_initial(1), 2, MPI_REAL,       &
                              MPI_SUM, comm2d, ierr )
          CALL MPI_ALLREDUCE( volume_flow_area_l(1), volume_flow_area(1), 2, MPI_REAL, MPI_SUM,    &
                              comm2d, ierr )

#else
          volume_flow_initial = volume_flow_initial_l
          volume_flow_area    = volume_flow_area_l
#endif

       ELSEIF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

          volume_flow_initial_l = 0.0_wp
          volume_flow_area_l    = 0.0_wp

          IF ( nxr == nx )  THEN
             DO  j = nys, nyn
                DO  k = nzb+1, nzt
                   volume_flow_initial_l(1) = volume_flow_initial_l(1) +                           &
                                              u(k,j,nx) * dzw(k)                                   &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,j,nx), 1 )              &
                                                     )
                   volume_flow_area_l(1)    = volume_flow_area_l(1) + dzw(k)                       &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,j,nx), 1 )              &
                                                     )
                ENDDO
             ENDDO
          ENDIF

          IF ( nyn == ny )  THEN
             DO  i = nxl, nxr
                DO  k = nzb+1, nzt
                   volume_flow_initial_l(2) = volume_flow_initial_l(2) +                           &
                                              v(k,ny,i) * dzw(k)                                   &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,ny,i), 2 )              &
                                                     )
                   volume_flow_area_l(2)    = volume_flow_area_l(2) + dzw(k)                       &
                                              * MERGE( 1.0_wp, 0.0_wp,                             &
                                                       BTEST( topo_flags(k,ny,i), 2 )              &
                                                     )
                ENDDO
             ENDDO
          ENDIF

#if defined( __parallel )
          CALL MPI_ALLREDUCE( volume_flow_initial_l(1), volume_flow_initial(1), 2, MPI_REAL,       &
                              MPI_SUM, comm2d, ierr )
          CALL MPI_ALLREDUCE( volume_flow_area_l(1), volume_flow_area(1), 2, MPI_REAL, MPI_SUM,    &
                              comm2d, ierr )

#else
          volume_flow_initial = volume_flow_initial_l
          volume_flow_area    = volume_flow_area_l
#endif

       ENDIF

!
!--    In case of 'bulk_velocity' mode, volume_flow_initial is calculated from u|v_bulk instead
       IF ( TRIM( conserve_volume_flow_mode ) == 'bulk_velocity' )  THEN
          volume_flow_initial(1) = u_bulk * volume_flow_area(1)
          volume_flow_initial(2) = v_bulk * volume_flow_area(2)
       ENDIF

    ENDIF
!
!-- In the following, surface properties can be further initialized with input from static driver
!-- file.
!-- At the moment this affects only default surfaces. For example, roughness length or sensible /
!-- latent heat fluxes can be initialized heterogeneously for default surfaces. Therefore, a generic
!-- routine from netcdf_data_input_mod is called to read a 2D array.
    IF ( input_pids_static )  THEN
!
!--    Open the static input file
#if defined( __netcdf )
       CALL open_read_file( TRIM( input_file_static ) //                                           &
                            TRIM( coupling_char ), pids_id )

       CALL inquire_num_variables( pids_id, num_var_pids )
!
!--    Allocate memory to store variable names and read them
       ALLOCATE( vars_pids(1:num_var_pids) )
       CALL inquire_variable_names( pids_id, vars_pids )
!
!--    Allocate memory for possible static input
       ALLOCATE( tmp_2d%var(nys:nyn,nxl:nxr) )
!
!--    Input roughness length.
       IF ( check_existence( vars_pids, 'z0' ) )  THEN

          tmp_2d%var = 0.0_wp
!
!--       Read _FillValue attribute
          CALL get_attribute( pids_id, char_fill, tmp_2d%fill, .FALSE., 'z0' )
!
!--       Read variable
          CALL get_variable( pids_id, 'z0', tmp_2d%var, nxl, nxr, nys, nyn )
          CALL add_ghost_layers( tmp_2d%var )
          CALL exchange_horiz_2d( tmp_2d%var )
          CALL set_lateral_neumann_bc( tmp_2d%var )
!
!--       Initialize roughness length. Note, z0 will be only initialized at default-type surfaces.
!--       At natural or urban z0 is implicitly initialized by the respective parameter lists.
!--       Initialize horizontal surface elements.
!--       Note, the actual 2D input arrays are only defined on the subdomain. Therefore, pass the
!--       index arrays with their respective offset values.
          CALL init_single_surface_properties( surf_def%z0, tmp_2d%var, surf_def%ns, tmp_2d%fill,  &
                                               surf_def%i+surf_def%ioff, surf_def%j+surf_def%joff, &
                                               surf_def%k+surf_def%koff)
       ENDIF
!
!--    Input surface sensible heat flux.
       IF ( check_existence( vars_pids, 'shf' ) )  THEN

          tmp_2d%var = 0.0_wp
!
!--       Read _FillValue attribute
          CALL get_attribute( pids_id, char_fill, tmp_2d%fill, .FALSE., 'shf' )
!
!--       Read variable
          CALL get_variable( pids_id, 'shf', tmp_2d%var, nxl, nxr, nys, nyn )
          CALL add_ghost_layers( tmp_2d%var )
          CALL exchange_horiz_2d( tmp_2d%var )
          CALL set_lateral_neumann_bc( tmp_2d%var )

          CALL init_single_surface_properties( surf_def%shf, tmp_2d%var, surf_def%ns, tmp_2d%fill, &
                                               surf_def%i+surf_def%ioff, surf_def%j+surf_def%joff, &
                                               surf_def%k+surf_def%koff, heatflux_input_conversion )
       ENDIF
!
!--    Input surface latent heat flux.
       IF ( humidity )  THEN
          IF ( check_existence( vars_pids, 'qsws' ) )  THEN

             tmp_2d%var = 0.0_wp
!
!--          Read _FillValue attribute
             CALL get_attribute( pids_id, char_fill, tmp_2d%fill, .FALSE., 'qsws' )
!
!--          Read variable
             CALL get_variable( pids_id, 'qsws', tmp_2d%var, nxl, nxr, nys, nyn )
             CALL add_ghost_layers( tmp_2d%var )
             CALL exchange_horiz_2d( tmp_2d%var )
             CALL set_lateral_neumann_bc( tmp_2d%var )

             CALL init_single_surface_properties( surf_def%qsws, tmp_2d%var, surf_def%ns,          &
                                                  tmp_2d%fill, surf_def%i+surf_def%ioff,           &
                                                  surf_def%j+surf_def%joff,                        &
                                                  surf_def%k+surf_def%koff,                        &
                                                  waterflux_input_conversion )
          ENDIF
       ENDIF
!
!--    Input passive scalar flux.
       IF ( passive_scalar )  THEN
          IF ( check_existence( vars_pids, 'ssws' ) )  THEN

             tmp_2d%var = 0.0_wp
!
!--          Read _FillValue attribute
             CALL get_attribute( pids_id, char_fill, tmp_2d%fill, .FALSE., 'ssws' )
!
!--          Read variable
             CALL get_variable( pids_id, 'ssws', tmp_2d%var, nxl, nxr, nys, nyn )
             CALL add_ghost_layers( tmp_2d%var )
             CALL exchange_horiz_2d( tmp_2d%var )
             CALL set_lateral_neumann_bc( tmp_2d%var )

             CALL init_single_surface_properties( surf_def%ssws, tmp_2d%var, surf_def%ns,          &
                                                  tmp_2d%fill, surf_def%i+surf_def%ioff,           &
                                                  surf_def%j+surf_def%joff,                        &
                                                  surf_def%k+surf_def%koff,                        &
                                                  scalarflux_input_conversion )
          ENDIF
       ENDIF
!
!--    Additional variables, can be initialized the
!--    same way.

!
!--    Finally, close the input file and deallocate temporary arrays
       DEALLOCATE( tmp_2d%var )
       DEALLOCATE( vars_pids )

       CALL close_input_file( pids_id )
#endif
    ENDIF
!
!-- Finally, if random_heatflux is set, disturb shf at horizontal surfaces. Actually, this should be
!-- done in surface_mod, where all other initializations of surface quantities are done. However,
!-- this would create a ring dependency, hence, it is done here. Maybe delete disturb_heatflux and
!-- tranfer the respective code directly into the initialization in surface_mod.
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.                                &
         .NOT. cyclic_fill_initialization )                                                        &
    THEN
       IF ( use_surface_fluxes  .AND.  constant_heatflux  .AND.  random_heatflux )  THEN
          IF ( surf_def%ns >= 1 )  CALL disturb_heatflux( surf_def )
          IF ( surf_lsm%ns >= 1 )  CALL disturb_heatflux( surf_lsm )
          IF ( surf_usm%ns >= 1 )  CALL disturb_heatflux( surf_usm )
       ENDIF
    ENDIF
!
!-- Initialize and count number of grid points used to calculate domain-averages
!-- including/excluding topography.
    ALLOCATE( mean_surface_level_height(0:statistic_regions),                                      &
              mean_surface_level_height_l(0:statistic_regions),                                    &
              ngp_2dh(0:statistic_regions),                                                        &
              ngp_2dh_l(0:statistic_regions),                                                      &
              ngp_3d(0:statistic_regions),                                                         &
              ngp_3d_inner(0:statistic_regions),                                                   &
              ngp_3d_inner_l(0:statistic_regions),                                                 &
              ngp_3d_inner_tmp(0:statistic_regions) )

    ALLOCATE( ngp_2dh_outer(nzb:nzt+1,0:statistic_regions),                                        &
              ngp_2dh_outer_l(nzb:nzt+1,0:statistic_regions),                                      &
              ngp_2dh_s_inner(nzb:nzt+1,0:statistic_regions),                                      &
              ngp_2dh_s_inner_l(nzb:nzt+1,0:statistic_regions) )
!
!-- Compute total sum of grid points and the mean surface level height for each statistic region.
!-- These are mainly used for horizontal averaging of turbulence statistics.
!-- ngp_2dh: number of grid points of a horizontal cross section through the respective statistic
!--          region
!-- ngp_3d:  number of grid points of the respective statistic region
    ngp_2dh_outer_l   = 0
    ngp_2dh_outer     = 0
    ngp_2dh_s_inner_l = 0
    ngp_2dh_s_inner   = 0
    ngp_2dh_l         = 0
    ngp_2dh           = 0
    ngp_3d_inner_l    = 0.0_wp
    ngp_3d_inner      = 0
    ngp_3d            = 0
    ngp_sums          = ( nz + 2 ) * ( pr_palm + max_pr_user )

    mean_surface_level_height   = 0.0_wp
    mean_surface_level_height_l = 0.0_wp
!
!-- To do: New concept for these non-topography grid points!
    DO  sr = 0, statistic_regions
       DO  i = nxl, nxr
          DO  j = nys, nyn
             IF ( rmask(j,i,sr) == 1.0_wp )  THEN
!
!--             All xy-grid points
                ngp_2dh_l(sr) = ngp_2dh_l(sr) + 1
!
!--             Determine mean surface-level height. In case of downward-facing walls are present,
!--             more than one surface level exist.
!--             In this case, use the lowest surface-level height.
                mean_surface_level_height_l(sr) = mean_surface_level_height_l(sr) +                &
                                                  zw(topo_top_ind(j,i,3))

                DO  k = nzb, nzt+1
!
!--                xy-grid points above topography.
!--                Calculate the number of atmosphere grid points not bounded by any walls
!--                (ngp_2dh_outer) indicated by bit 24, as well as the number of atmosphere grid
!--                points (ngp_2dh_s_inner) including boundary values, indicated by bit 22.
                   ngp_2dh_outer_l(k,sr) = ngp_2dh_outer_l(k,sr)     +                             &
                                           MERGE( 1, 0, BTEST( topo_flags(k,j,i), 24 ) )

                   ngp_2dh_s_inner_l(k,sr) = ngp_2dh_s_inner_l(k,sr) +                             &
                                             MERGE( 1, 0, BTEST( topo_flags(k,j,i), 22 ) )

                ENDDO
!
!--             All grid points of the total domain above topography
                ngp_3d_inner_l(sr) = ngp_3d_inner_l(sr) + ( nz - topo_top_ind(j,i,0) + 2 )

             ENDIF
          ENDDO
       ENDDO
    ENDDO

    sr = statistic_regions + 1
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_2dh_l(0), ngp_2dh(0), sr, MPI_INTEGER, MPI_SUM, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_2dh_outer_l(0,0), ngp_2dh_outer(0,0), (nz+2)*sr, MPI_INTEGER, MPI_SUM, &
                        comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_2dh_s_inner_l(0,0), ngp_2dh_s_inner(0,0), (nz+2)*sr, MPI_INTEGER,      &
                        MPI_SUM, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_3d_inner_l(0), ngp_3d_inner_tmp(0), sr, MPI_REAL, MPI_SUM, comm2d,     &
                        ierr )
    ngp_3d_inner = INT( ngp_3d_inner_tmp, KIND = SELECTED_INT_KIND( 18 ) )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( mean_surface_level_height_l(0), mean_surface_level_height(0), sr, MPI_REAL,&
                        MPI_SUM, comm2d, ierr )
    mean_surface_level_height = mean_surface_level_height / REAL( ngp_2dh )
#else
    ngp_2dh         = ngp_2dh_l
    ngp_2dh_outer   = ngp_2dh_outer_l
    ngp_2dh_s_inner = ngp_2dh_s_inner_l
    ngp_3d_inner    = INT( ngp_3d_inner_l, KIND = SELECTED_INT_KIND( 18 ) )
    mean_surface_level_height = mean_surface_level_height_l / REAL( ngp_2dh_l )
#endif

    ngp_3d = INT ( ngp_2dh, KIND = SELECTED_INT_KIND( 18 ) ) *                                     &
             INT ( (nz + 2 ), KIND = SELECTED_INT_KIND( 18 ) )

!
!-- Set a lower limit of 1 in order to avoid zero divisions in flow_statistics, buoyancy, etc. A
!-- zero value will occur for cases where all grid points of the respective subdomain lie below the
!-- surface topography
    ngp_2dh_outer   = MAX( 1, ngp_2dh_outer(:,:)   )
    ngp_3d_inner    = MAX( INT(1, KIND = SELECTED_INT_KIND( 18 )), ngp_3d_inner(:) )
    ngp_2dh_s_inner = MAX( 1, ngp_2dh_s_inner(:,:) )

    DEALLOCATE( mean_surface_level_height_l, ngp_2dh_l, ngp_2dh_outer_l, ngp_3d_inner_l,           &
                ngp_3d_inner_tmp )
!
!-- Compute number of prognostic w-grid points. This is only required for mean vertical velocity
!-- removal in case of bottom and top Neumann boundary conditions before the pressure solver is
!-- invoked. Note, the removal will not be done for offline nested simulations, and for only nested
!-- simulations, where childs do not cover the complete root domain.
    IF ( ibc_p_b == 1  .AND.  ibc_p_t == 1  .AND.  .NOT. nesting_offline  .AND.                    &
         .NOT. ( child_domain  .AND.  nesting_bounds /= 'vertical_only' ) )  THEN
       ALLOCATE( ngp_2dh_wgrid(nzb+1:nzt) )
       ngp_2dh_wgrid = 0
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                ngp_2dh_wgrid(k) = ngp_2dh_wgrid(k) + MERGE( 1, 0, BTEST( topo_flags(k,j,i), 3 ) )
             ENDDO
          ENDDO
       ENDDO
#if defined( __parallel )
       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, ngp_2dh_wgrid(1), nzt-nzb, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#endif
!
!--    To avoid divisions by zero (this may happen if an entire prognostic level is occupied by
!--    topography) but also to avoid recurrent checks on this, set a minimum value of 1
!--    (at these levels there won't be any correction at all).
       ngp_2dh_wgrid = MERGE( 1, ngp_2dh_wgrid, ngp_2dh_wgrid == 0 )
    ENDIF

!
!-- Initializing actions for all modules which impact the boundary conditions and need to be
!-- initialized before the pressure solver is called.
    CALL module_interface_init_before_pressure_solver
!
!-- Initialize quantities for special advections schemes
    CALL init_advec

!
!-- Impose random perturbation on the horizontal velocity field and then
!-- remove the divergences from the velocity field at the initial stage
    IF ( create_disturbances  .AND.  disturbance_energy_limit /= 0.0_wp  .AND.                     &
         TRIM( initializing_actions ) /= 'read_restart_data'  .AND.                                &
         .NOT. cyclic_fill_initialization )                                                        &
    THEN

       IF ( debug_output )  THEN
          CALL debug_message( 'creating disturbances + applying pressure solver', 'start' )
       ENDIF
!
!--    Needed for both disturb_field and pres
!$ACC DATA &
!$ACC CREATE(tend(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(u(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(v(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

       CALL disturb_field( 'u', tend, u )
       CALL disturb_field( 'v', tend, v )

!$ACC DATA &
!$ACC CREATE(d(nzb+1:nzt,nys:nyn,nxl:nxr)) &
!$ACC COPY(w(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPY(p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)) &
!$ACC COPYIN(rho_air(nzb:nzt+1), rho_air_zw(nzb:nzt+1)) &
!$ACC COPYIN(ddzu(1:nzt+1), ddzw(1:nzt+1)) &
!$ACC COPYIN(topo_flags(nzb:nzt+1,nysg:nyng,nxlg:nxrg))

       n_sor = nsor_ini
       CALL pres
       n_sor = nsor

!$ACC END DATA
!$ACC END DATA

       IF ( debug_output )  THEN
          CALL debug_message( 'creating disturbances + applying pressure solver', 'end' )
       ENDIF

    ENDIF

    IF ( .NOT. ocean_mode )  THEN

       ALLOCATE( hyp(nzb:nzt+1) )
       ALLOCATE( d_exner(nzb:nzt+1) )
       ALLOCATE( exner(nzb:nzt+1) )
       ALLOCATE( hyrho(nzb:nzt+1) )
!
!--    Check temperature in case of too large domain height
       DO  k = nzb, nzt+1
          IF ( ( pt_surface * exner_function( surface_pressure * 100.0_wp ) - g/c_p * zu(k) )      &
                 < 0.0_wp )  THEN
             WRITE( message_string, * )  'absolute temperature < 0.0 at zu(', k, ') = ', zu(k)
             CALL message( 'init_3d_model', 'PAC0209', 1, 2, 0, 6, 0 )
          ENDIF
       ENDDO

!
!--    Calculate vertical profile of the hydrostatic pressure (hyp)
       hyp    = barometric_formula( zu, pt_surface * exner_function( surface_pressure * 100.0_wp ),&
                                    surface_pressure * 100.0_wp )
       d_exner = exner_function_invers( hyp )
       exner = 1.0_wp / exner_function_invers( hyp )
       hyrho  = ideal_gas_law_rho_pt( hyp, pt_init )
!
!--    Compute reference density
       rho_surface = ideal_gas_law_rho( surface_pressure * 100.0_wp,                               &
                                        pt_surface * exner_function( surface_pressure * 100.0_wp ) )

    ENDIF

!
!-- If required, initialize particles
    IF ( agents_active )  CALL mas_init
!
!-- Initializing actions for all other modules which can be initialized after the pressure
!-- solver has been called.
    CALL module_interface_init_after_pressure_solver
!
!-- Initialize mesoscale offline nesting for other modules that do not modify the flow field.
    CALL nesting_offl_init_modules
!
!-- Initialize surface layer (done after LSM as roughness length are required for initialization
    IF ( constant_flux_layer )  CALL init_surface_layer_fluxes
!
!-- Initialize surface data output
    IF ( surface_output )  CALL surface_data_output_init
!
!-- Initialize the ws-scheme.
    IF ( ws_scheme_sca  .OR.  ws_scheme_mom )  CALL ws_init
!
!-- Perform post-initializing checks for all other modules
    CALL module_interface_init_checks
!
!-- Add other module specific timeseries
    CALL module_interface_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )

!
!-- Initialize surface forcing corresponding to large-scale forcing. Therein,
!-- initialize heat-fluxes, etc. via datatype. Revise it later!
    IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
       IF ( use_surface_fluxes  .AND.  constant_heatflux )  THEN
          CALL ls_forcing_surf( simulated_time )
       ENDIF
    ENDIF
!
!-- Setting weighting factors for calculation of perturbation pressure and turbulent quantities from
!-- the RK substeps.
    IF ( TRIM( timestep_scheme ) == 'runge-kutta-3' )  THEN      ! for RK3-method

       weight_substep(1) = 1._wp/6._wp
       weight_substep(2) = 3._wp/10._wp
       weight_substep(3) = 8._wp/15._wp

       weight_pres(1)    = 1._wp/3._wp
       weight_pres(2)    = 5._wp/12._wp
       weight_pres(3)    = 1._wp/4._wp

    ELSEIF ( TRIM( timestep_scheme ) == 'runge-kutta-2' )  THEN  ! for RK2-method

       weight_substep(1) = 1._wp/2._wp
       weight_substep(2) = 1._wp/2._wp

       weight_pres(1)    = 1._wp/2._wp
       weight_pres(2)    = 1._wp/2._wp

    ELSE                                     ! for Euler-method

       weight_substep(1) = 1.0_wp
       weight_pres(1)    = 1.0_wp

    ENDIF

!
!-- Initialize Rayleigh damping factors
    rdf    = 0.0_wp
    rdf_sc = 0.0_wp
    IF ( rayleigh_damping_factor /= 0.0_wp )  THEN

       IF (  .NOT.  ocean_mode )  THEN
          DO  k = nzb+1, nzt
             IF ( zu(k) >= rayleigh_damping_height )  THEN
                rdf(k) = rayleigh_damping_factor *                                                 &
                         ( SIN( pi * 0.5_wp * ( zu(k) - rayleigh_damping_height )                  &
                                / ( zu(nzt) - rayleigh_damping_height ) )                          &
                         )**2
             ENDIF
          ENDDO
       ELSE
!
!--       In ocean mode, rayleigh damping is applied in the lower part of the model domain
          DO  k = nzt, nzb+1, -1
             IF ( zu(k) <= rayleigh_damping_height )  THEN
                rdf(k) = rayleigh_damping_factor *                                                 &
                         ( SIN( pi * 0.5_wp * ( rayleigh_damping_height - zu(k) )                  &
                                / ( rayleigh_damping_height - zu(nzb+1) ) )                        &
                         )**2
             ENDIF
          ENDDO
       ENDIF

    ENDIF
    IF ( scalar_rayleigh_damping )  rdf_sc = rdf

!
!-- Initialize the starting level and the vertical smoothing factor used for the external pressure
!-- gradient
    dp_smooth_factor = 1.0_wp
    IF ( dp_external )  THEN
!
!--    Set the starting level dp_level_ind_b only if it has not been set before (e.g. in init_grid).
       IF ( dp_level_ind_b == 0 )  THEN
          ind_array = MINLOC( ABS( dp_level_b - zu ) )
          dp_level_ind_b = ind_array(1) - 1 + nzb
                                        ! MINLOC uses lower array bound 1
       ENDIF
       IF ( dp_smooth )  THEN
          dp_smooth_factor(:dp_level_ind_b) = 0.0_wp
          DO  k = dp_level_ind_b+1, nzt
             dp_smooth_factor(k) = 0.5_wp * ( 1.0_wp + SIN( pi *                                   &
                                             ( REAL( k - dp_level_ind_b, KIND=wp ) /               &
                                               REAL( nzt - dp_level_ind_b, KIND=wp ) - 0.5_wp ) ) )
          ENDDO
       ENDIF
    ENDIF

!
!-- Initialize damping zone for the potential temperature in case of non-cyclic lateral boundaries.
!-- The damping zone has the maximum value at the inflow boundary and decreases to zero at
!-- pt_damping_width.
    ptdf_x = 0.0_wp
    ptdf_y = 0.0_wp
    IF ( bc_lr_dirrad )  THEN
       DO  i = nxl, nxr
          IF ( ( i * dx ) < pt_damping_width )  THEN
             ptdf_x(i) = pt_damping_factor * ( SIN( pi * 0.5_wp *                                  &
                                                    REAL( pt_damping_width - i * dx, KIND=wp ) /   &
                                                    REAL( pt_damping_width, KIND=wp ) ) )**2
                                  ENDIF
       ENDDO
    ELSEIF ( bc_lr_raddir )  THEN
       DO  i = nxl, nxr
          IF ( ( i * dx ) > ( nx * dx - pt_damping_width ) )  THEN
             ptdf_x(i) = pt_damping_factor * SIN( pi * 0.5_wp *                                    &
                                                  ( ( i - nx ) * dx + pt_damping_width ) /         &
                                                  REAL( pt_damping_width, KIND=wp ) )**2
          ENDIF
       ENDDO
    ELSEIF ( bc_ns_dirrad )  THEN
       DO  j = nys, nyn
          IF ( ( j * dy ) > ( ny * dy - pt_damping_width ) )  THEN
             ptdf_y(j) = pt_damping_factor * SIN( pi * 0.5_wp *                                    &
                                                  ( ( j - ny ) * dy + pt_damping_width ) /         &
                                                  REAL( pt_damping_width, KIND=wp ) )**2
          ENDIF
       ENDDO
    ELSEIF ( bc_ns_raddir )  THEN
       DO  j = nys, nyn
          IF ( ( j * dy ) < pt_damping_width )  THEN
             ptdf_y(j) = pt_damping_factor * SIN( pi * 0.5_wp *                                    &
                                                  ( pt_damping_width - j * dy ) /                  &
                                                  REAL( pt_damping_width, KIND=wp ) )**2
          ENDIF
       ENDDO
    ENDIF

!
!-- Initialize outflow damping zone for u, v, w and pt in case of non-cyclic lateral boundaries.
!-- The damping zone has the maximum value at the outflow boundary and decreases to zero at a
!-- distance of outflow_damping_width to the outflow boundary.
    odf_x = 0.0_wp
    odf_y = 0.0_wp
    IF ( bc_lr_dirrad )  THEN
       DO  i = nxl, nxr
          IF ( ( i * dx ) > ( nx * dx - outflow_damping_width ) )  THEN
             odf_x(i) = outflow_damping_factor *                                                   &
                        SIN( pi * 0.5_wp * ( ( i - nx ) * dx + outflow_damping_width ) /           &
                             outflow_damping_width )**2
          ENDIF
       ENDDO
    ELSEIF ( bc_lr_raddir )  THEN
       DO  i = nxl, nxr
          IF ( ( i * dx ) < outflow_damping_width )  THEN
             odf_x(i) = outflow_damping_factor *                                                   &
                        SIN( pi * 0.5_wp * ( outflow_damping_width - i * dx ) /                    &
                             outflow_damping_width )**2
          ENDIF
       ENDDO
    ELSEIF ( bc_ns_dirrad )  THEN
       DO  j = nys, nyn
          IF ( ( j * dy ) < outflow_damping_width )  THEN
             odf_y(j) = outflow_damping_factor *                                                   &
                        SIN( pi * 0.5_wp * ( outflow_damping_width - j * dy ) /                    &
                             outflow_damping_width )**2
          ENDIF
       ENDDO
    ELSEIF ( bc_ns_raddir )  THEN
       DO  j = nys, nyn
          IF ( ( j * dy ) > ( ny * dy - outflow_damping_width ) )  THEN
             odf_y(j) = outflow_damping_factor *                                                   &
                        SIN( pi * 0.5_wp * ( ( j - ny ) * dy + outflow_damping_width ) /           &
                             outflow_damping_width )**2
          ENDIF
       ENDDO
    ENDIF

!
!-- Input binary data file is not needed anymore. This line must be placed after call of user_init!
    CALL close_file( 13 )
!
!-- Finally, initialize new time levels again. This is to guarantee that boundary values
!-- are set adequately.
    pt_p = pt; u_p = u; v_p = v; w_p = w
    IF ( humidity )  THEN
       q_p = q
    ENDIF
    IF ( passive_scalar )  s_p = s
!
!-- In case of nesting/coupling, put a barrier to assure that all parent and child domains finished
!-- initialization.
    IF ( nested_run  .OR.  atmosphere_ocean_coupled_run )  THEN
#if defined( __parallel )
       CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
#else
       CONTINUE
#endif
    ENDIF

    CALL location_message( 'model initialization', 'finished' )

 END SUBROUTINE init_3d_model
