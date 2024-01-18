!> @file pmc_interface_mod.f90
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
! Domain nesting interface routines. The low-level inter-domain communication is conducted by the
! PMC-library routines.
!
! @todo Improve child initializing_actions steering
! @todo Remove array_3d variables from USE statements thate not used in the routine
! @todo Data transfer of qc and nc is prepared but not activated
!--------------------------------------------------------------------------------------------------!
 MODULE pmc_interface

#if ! defined( __parallel )
!
!-- Serial mode does not allow nesting, but requires the following variables as steering quantities
    USE kinds

    IMPLICIT NONE

    PUBLIC

    CHARACTER(LEN=14), SAVE ::  nesting_bounds = 'none'  !< steering parameter for outer boundary conditions of child domains
    CHARACTER(LEN=8),  SAVE ::  nesting_mode   = 'none'  !< steering parameter for 1- or 2-way nesting

    INTEGER(iwp), SAVE ::  comm_world_nesting  !< Global nesting communicator
    INTEGER(iwp), SAVE ::  cpl_id  = 1         !<

    LOGICAL, SAVE ::  atmosphere_ocean_coupled_run = .FALSE.     !< general switch
    LOGICAL, SAVE ::  homogeneous_initialization_child = .FALSE. !< switch to control initialization of child domains (default .FALSE.)
    LOGICAL, SAVE ::  nesting_bounds_vertical_only = .FALSE.     !< general switch
    LOGICAL, SAVE ::  particle_coupling = .FALSE.                !< switch for particle coupling (meaningful only when lpm is used)
    LOGICAL, SAVE ::  rans_mode_parent = .FALSE.                 !< parent model mode (.F.-LES mode, .T.-RANS mode)
    LOGICAL, SAVE ::  root_model = .TRUE.                        !< root model flag

    REAL(wp), SAVE ::  lower_left_coord_x = 0.0_wp  !< x-coordinate of the lower left corner of the domain
    REAL(wp), SAVE ::  lower_left_coord_y = 0.0_wp  !< y-coordinate of the lower left corner of the domain
    REAL(wp), SAVE ::  nest_shift_z = 0.0_wp        !< z-coordinate of the bottom of the domain

#else

    USE ISO_C_BINDING


    USE arrays_3d,                                                                                 &
        ONLY:  diss,                                                                               &
               diss_2,                                                                             &
               dzu,                                                                                &
               dzw,                                                                                &
               e,                                                                                  &
               e_p,                                                                                &
               e_2,                                                                                &
               nc,                                                                                 &
               nc_2,                                                                               &
               nc_p,                                                                               &
               nr,                                                                                 &
               nr_2,                                                                               &
               pt,                                                                                 &
               pt_2,                                                                               &
               q,                                                                                  &
               q_2,                                                                                &
               qc,                                                                                 &
               qc_2,                                                                               &
               qr,                                                                                 &
               qr_2,                                                                               &
               s,                                                                                  &
               s_2,                                                                                &
               sa,                                                                                 &
               sa_2,                                                                               &
               u,                                                                                  &
               u_p,                                                                                &
               u_2,                                                                                &
               v,                                                                                  &
               v_p,                                                                                &
               v_2,                                                                                &
               w,                                                                                  &
               w_p,                                                                                &
               w_2,                                                                                &
               x,                                                                                  &
               y,                                                                                  &
               zu,                                                                                 &
               zw

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model,                                                                   &
               microphysics_morrison,                                                              &
               microphysics_seifert

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nspec

    USE chem_modules,                                                                              &
        ONLY:  chem_species,                                                                       &
               ibc_cs_b,                                                                           &
               nesting_chem

    USE chemistry_model_mod,                                                                       &
        ONLY:  spec_conc_2

    USE control_parameters,                                                                        &
        ONLY:  air_chemistry,                                                                      &
               atmosphere_run_coupled_to_ocean,                                                    &
               bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_lr_cyc,                                                                          &
               bc_ns_cyc,                                                                          &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               child_domain,                                                                       &
               constant_diffusion,                                                                 &
               constant_flux_layer,                                                                &
               coupling_char,                                                                      &
               coupling_start_time,                                                                &
               debug_output,                                                                       &
               debug_output_timestep,                                                              &
               debug_string,                                                                       &
               dt_coupling,                                                                        &
               dt_restart,                                                                         &
               dt_3d,                                                                              &
               dz,                                                                                 &
               end_time,                                                                           &
               homogenize_surface_temperature,                                                     &
               humidity,                                                                           &
               humidity_remote,                                                                    &
               ibc_pt_b,                                                                           &
               ibc_q_b,                                                                            &
               ibc_s_b,                                                                            &
               ibc_uv_b,                                                                           &
               message_string,                                                                     &
               nested_run,                                                                         &
               neutral,                                                                            &
               ocean_mode,                                                                         &
               ocean_run_coupled_to_atmosphere,                                                    &
               passive_scalar,                                                                     &
               rans_mode,                                                                          &
               rans_tke_e,                                                                         &
               restart_time,                                                                       &
               roughness_length,                                                                   &
               salinity,                                                                           &
               salsa,                                                                              &
               time_restart,                                                                       &
               time_since_reference_point,                                                         &
               terminate_coupled,                                                                  &
               topography,                                                                         &
               volume_flow

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxlu,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nysv,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzb_max,                                                                            &
               nzt,                                                                                &
               topo_top_ind,                                                                       &
               topo_flags

    USE particle_attributes,                                                                       &
        ONLY:  particle_advection

    USE kinds

    USE MPI

    USE pegrid,                                                                                    &
        ONLY:  collective_wait,                                                                    &
               comm1dx,                                                                            &
               comm1dy,                                                                            &
               comm2d,                                                                             &
               myid,                                                                               &
               myidx,                                                                              &
               myidy,                                                                              &
               npex,                                                                               &
               npey,                                                                               &
               numprocs,                                                                           &
               pleft,                                                                              &
               pnorth,                                                                             &
               pright,                                                                             &
               psouth,                                                                             &
               status

    USE pmc_child,                                                                                 &
        ONLY:  pmc_childinit,                                                                      &
               pmc_c_clear_next_array_list,                                                        &
               pmc_c_finalize,                                                                     &
               pmc_c_getnextarray,                                                                 &
               pmc_c_get_2d_index_list,                                                            &
               pmc_c_getbuffer,                                                                    &
               pmc_c_putbuffer,                                                                    &
               pmc_c_setind_and_allocmem,                                                          &
               pmc_c_set_dataarray,                                                                &
               pmc_set_dataarray_name

    USE pmc_general,                                                                               &
        ONLY:  da_namelen,                                                                         &
               pmc_max_array

    USE pmc_handle_communicator,                                                                   &
        ONLY:  pmc_get_model_info,                                                                 &
               pmc_init_model,                                                                     &
               pmc_no_namelist_found,                                                              &
               pmc_parent_for_child,                                                               &
               m_couplers

    USE pmc_mpi_wrapper,                                                                           &
        ONLY:  pmc_bcast,                                                                          &
               pmc_recv_from_child,                                                                &
               pmc_recv_from_parent,                                                               &
               pmc_send_to_child,                                                                  &
               pmc_send_to_parent

    USE pmc_parent,                                                                                &
        ONLY:  pmc_parentinit,                                                                     &
               pmc_s_clear_next_array_list,                                                        &
               pmc_s_fillbuffer,                                                                   &
               pmc_s_finalize,                                                                     &
               pmc_s_getdata_from_buffer,                                                          &
               pmc_s_getnextarray,                                                                 &
               pmc_s_setind_and_allocmem,                                                          &
               pmc_s_set_active_data_array,                                                        &
               pmc_s_set_dataarray,                                                                &
               pmc_s_set_2d_index_list

    USE salsa_mod,                                                                                 &
        ONLY:  aerosol_mass,                                                                       &
               aerosol_number,                                                                     &
               gconc_2,                                                                            &
               ibc_aer_b,                                                                          &
               mconc_2,                                                                            &
               nbins_aerosol,                                                                      &
               ncomponents_mass,                                                                   &
               nconc_2,                                                                            &
               nesting_salsa,                                                                      &
               ngases_salsa,                                                                       &
               salsa_gas,                                                                          &
               salsa_gases_from_chem

    USE surface_coupler_mod,                                                                       &
        ONLY:  child_recv,                                                                         &
               child_send,                                                                         &
               parent_recv,                                                                        &
               parent_send,                                                                        &
               surface_coupler_buffer_handling,                                                    &
               surface_coupler_alloc_mem,                                                          &
               surface_coupler_exchange_array_1,                                                   &
               surface_coupler_exchange_array_2,                                                   &
               surface_coupler_exchange_array_3,                                                   &
               surface_coupler_exchange_array_4

    USE surface_mod,                                                                               &
        ONLY:  bc_hv

    IMPLICIT NONE


    PRIVATE
!
!-- Constants
    INTEGER(iwp), PARAMETER ::  child_to_parent = 2            !< Parameter for pmci_parent_datatrans indicating the direction of
                                                               !< transfer
    INTEGER(iwp), PARAMETER ::  interpolation_scheme_lrsn = 2  !< Interpolation scheme to be used on lateral boundaries
    INTEGER(iwp), PARAMETER ::  interpolation_scheme_t = 3     !< Interpolation scheme to be used on top boundary
    INTEGER(iwp), PARAMETER ::  parent_to_child = 1            !< Parameter for pmci_parent_datatrans indicating the direction of
                                                               !< transfer

    REAL(wp), PARAMETER ::  tolefac = 1.0E-6_wp  !< Relative tolerence for grid-line matching tests and comparisons
!
!-- Coupler setup
    CHARACTER(LEN=32), SAVE ::  cpl_name  !<

    INTEGER(iwp), SAVE ::  comm_world_nesting  !< Global nesting communicator
    INTEGER(iwp), SAVE ::  cpl_id  = 1         !< Model (domain) id (1 for root, 2,... for nested domains).
    INTEGER(iwp), SAVE ::  cpl_npe_total       !<
    INTEGER(iwp), SAVE ::  cpl_parent_id       !<

!
!-- Control parameters
    CHARACTER(LEN=14), SAVE ::  nesting_bounds = '3d_nested'  !< steering parameter for outer boundary conditions of child domains
    CHARACTER(LEN=7),  SAVE ::  nesting_datatransfer_mode = 'mixed'  !< steering parameter for data-transfer mode
    CHARACTER(LEN=8),  SAVE ::  nesting_mode = 'two-way'             !< steering parameter for 1- or 2-way nesting

    INTEGER(iwp), SAVE ::  anterpolation_buffer_width = 2  !< Boundary buffer width for anterpolation

    REAL(wp), SAVE ::  anterpolation_starting_height = 9999999.9_wp  !< steering parameter for canopy restricted anterpolation

    LOGICAL, SAVE ::  homogeneous_initialization_child = .FALSE. !< switch to control initialization of child domains (default .FALSE.)
    LOGICAL, SAVE ::  atmosphere_ocean_coupled_run = .FALSE. !< general switch
    LOGICAL, SAVE ::  nested_bc_at_bottom = .FALSE.          !< general switch
    LOGICAL, SAVE ::  nested_bc_at_top = .FALSE.             !< general switch
    LOGICAL, SAVE ::  nesting_bounds_vertical_only = .FALSE. !< general switch
    LOGICAL, SAVE ::  particle_coupling = .TRUE.  !< switch for particle coupling (meaningful only when lpm is used)
    LOGICAL, SAVE ::  rans_mode_parent = .FALSE.  !< mode of parent model (.F. - LES mode, .T. - RANS mode)
    LOGICAL, SAVE ::  root_model = .TRUE.         !< indicates, if model is root model
!
!-- Geometry
    REAL(wp), SAVE, DIMENSION(:), ALLOCATABLE, PUBLIC ::  coord_x                !< Array for the absolute x-coordinates
    REAL(wp), SAVE, DIMENSION(:), ALLOCATABLE, PUBLIC ::  coord_y                !< Array for the absolute y-coordinates
    REAL(wp), SAVE, PUBLIC                            ::  lower_left_coord_x     !< x-coordinate of the lower left corner of the domain
    REAL(wp), SAVE, PUBLIC                            ::  lower_left_coord_y     !< y-coordinate of the lower left corner of the domain
    REAL(wp), SAVE                                    ::  nest_shift_z = 0.0_wp  !< z-coordinate of the bottom of the domain
    REAL(wp), SAVE                                    ::  parent_nest_shift_z    !< Parent z-coordinate of the bottom of the domain
!
!-- Children's parent-grid arrays
    INTEGER(iwp), SAVE, DIMENSION(:,:), ALLOCATABLE :: kpb_anterp  !< Lower limit of kp in anterpolation
    INTEGER(iwp), SAVE, DIMENSION(5), PUBLIC ::  parent_bound  !< subdomain index bounds for children's parent-grid arrays

    INTEGER(idp), SAVE, DIMENSION(:,:), ALLOCATABLE, TARGET, PUBLIC ::  nr_partc   !<
    INTEGER(idp), SAVE, DIMENSION(:,:), ALLOCATABLE, TARGET, PUBLIC ::  part_adrc  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_p_init             !< Parent-grid 1D profile on child domain - u-component
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_p_init             !< Parent-grid 1D profile on child domain - v-component
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  e_p_init             !< Parent-grid 1D profile on child domain - SGS-TKE
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  diss_p_init          !< Parent-grid 1D profile on child domain - dissipation rate
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_p_init            !< Parent-grid 1D profile on child domain - pot. temperature
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  q_p_init             !< Parent-grid 1D profile on child domain - mixing ratio
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qc_p_init            !< Parent-grid 1D profile on child domain - cloud water
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  nc_p_init            !< Parent-grid 1D profile on child domain - cloud water number conc.
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  qr_p_init            !< Parent-grid 1D profile on child domain - rain water
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  nr_p_init            !< Parent-grid 1D profile on child domain - rain water number conc.
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  s_p_init             !< Parent-grid 1D profile on child domain - passive scalar

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  chem_p_init          !< Parent-grid 1D profile on child domain - chemistry
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  aerosol_number_p_init!< Parent-grid 1D profile on child domain - aerosol number
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  aerosol_mass_p_init  !< Parent-grid 1D profile on child domain - aerosol mass
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  salsa_gas_p_init     !< Parent-grid 1D profile on child domain - aerosol gas

    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  dissc  !< Parent-grid array on child domain - dissipation rate
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ec     !< Parent-grid array on child domain - SGS TKE
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  nrc    !< Parent-grid array on child domain -
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ncc    !< Parent-grid array on child domain -
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  ptc    !< Parent-grid array on child domain - potential temperature
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  q_c    !< Parent-grid array on child domain -
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  qcc    !< Parent-grid array on child domain -
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  qrc    !< Parent-grid array on child domain -
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sc     !< Parent-grid array on child domain - scalar
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sac    !< Parent-grid array on child domain - salinity
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  uc     !< Parent-grid array on child domain - velocity component u
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  vc     !< Parent-grid array on child domain - velocity component v
    REAL(wp), SAVE, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  wc     !< Parent-grid array on child domain - velocity component w

    REAL(wp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  aerosol_mass_c    !< Aerosol mass
    REAL(wp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  aerosol_number_c  !< Aerosol number
    REAL(wp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  chem_spec_c       !< Parent-grid array on child domain
                                                                                  !< - chemical species
    REAL(wp), SAVE, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::  salsa_gas_c       !< SALSA gases

    REAL(wp), SAVE, DIMENSION(:,:), ALLOCATABLE, TARGET, PUBLIC ::  surface_coupler_exchange_array_1c  !< for atmosphere-ocean data exchange
    REAL(wp), SAVE, DIMENSION(:,:), ALLOCATABLE, TARGET, PUBLIC ::  surface_coupler_exchange_array_2c  !< for atmosphere-ocean data exchange
    REAL(wp), SAVE, DIMENSION(:,:), ALLOCATABLE, TARGET, PUBLIC ::  surface_coupler_exchange_array_3c  !< for atmosphere-ocean data exchange
    REAL(wp), SAVE, DIMENSION(:,:), ALLOCATABLE, TARGET, PUBLIC ::  surface_coupler_exchange_array_4c  !< for atmosphere-ocean data exchange

!
!-- Grid-spacing ratios.
    INTEGER(iwp), SAVE ::  igsr  !< Integer grid-spacing ratio in i-direction
    INTEGER(iwp), SAVE ::  jgsr  !< Integer grid-spacing ratio in j-direction
    INTEGER(iwp), SAVE ::  kgsr  !< Integer grid-spacing ratio in k-direction
!
!-- Global parent-grid index bounds
    INTEGER(iwp), SAVE ::  iplg  !< Leftmost parent-grid array ip index of the whole child domain
    INTEGER(iwp), SAVE ::  iprg  !< Rightmost parent-grid array ip index of the whole child domain
    INTEGER(iwp), SAVE ::  jpsg  !< Southmost parent-grid array jp index of the whole child domain
    INTEGER(iwp), SAVE ::  jpng  !< Northmost parent-grid array jp index of the whole child domain
!
!-- Local parent-grid index bounds. Different sets of index bounds are needed for parent-grid arrays
!-- (uc, etc), for index mapping arrays (iflu, etc) and for work arrays (workarr_lr, etc). This is
!-- because these arrays have different dimensions depending on the location of the subdomain
!-- relative to boundaries and corners.
    INTEGER(iwp), SAVE ::  ipl   !< Left index limit for children's parent-grid arrays
    INTEGER(iwp), SAVE ::  ipla  !< Left index limit for allocation of index-mapping and other auxiliary arrays
    INTEGER(iwp), SAVE ::  iplw  !< Left index limit for children's parent-grid work arrays
    INTEGER(iwp), SAVE ::  ipr   !< Right index limit for children's parent-grid arrays
    INTEGER(iwp), SAVE ::  ipra  !< Right index limit for allocation of index-mapping and other auxiliary arrays
    INTEGER(iwp), SAVE ::  iprw  !< Right index limit for children's parent-grid work arrays
    INTEGER(iwp), SAVE ::  jpn   !< North index limit for children's parent-grid arrays
    INTEGER(iwp), SAVE ::  jpna  !< North index limit for allocation of index-mapping and other auxiliary arrays
    INTEGER(iwp), SAVE ::  jpnw  !< North index limit for children's parent-grid work arrays
    INTEGER(iwp), SAVE ::  jps   !< South index limit for children's parent-grid arrays
    INTEGER(iwp), SAVE ::  jpsa  !< South index limit for allocation of index-mapping and other auxiliary arrays
    INTEGER(iwp), SAVE ::  jpsw  !< South index limit for children's parent-grid work arrays
!
!-- Highest prognostic parent-grid k-indices.
    INTEGER(iwp), SAVE ::  kcto  !< Upper bound for k in anterpolation of variables other than w.
    INTEGER(iwp), SAVE ::  kctw  !< Upper bound for k in anterpolation of w.
!
!-- Child-array indices to be precomputed and stored for anterpolation.
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  iflu   !< child index indicating left bound of parent grid box on u-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  ifuu   !< child index indicating right bound of parent grid box on u-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  iflo   !< child index indicating left bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  ifuo   !< child index indicating right bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jflv   !< child index indicating south bound of parent grid box on v-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jfuv   !< child index indicating north bound of parent grid box on v-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jflo   !< child index indicating south bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  jfuo   !< child index indicating north bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kflw   !< child index indicating lower bound of parent grid box on w-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kfuw   !< child index indicating upper bound of parent grid box on w-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kflo   !< child index indicating lower bound of parent grid box on scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:) ::  kfuo   !< child index indicating upper bound of parent grid box on scalar-grid
!
!-- Number of child-grid nodes within anterpolation cells to be precomputed for anterpolation.
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  ijkfc_s  !< number of child grid points contributing to a parent grid
                                                                   !< node in anterpolation, scalar-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  ijkfc_u  !< number of child grid points contributing to a parent grid
                                                                   !< node in anterpolation, u-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  ijkfc_v  !< number of child grid points contributing to a parent grid
                                                                   !< node in anterpolation, v-grid
    INTEGER(iwp), SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::  ijkfc_w  !< number of child grid points contributing to a parent grid
                                                                   !< node in anterpolation, w-grid
!
!-- Work arrays for interpolation and user-defined type definitions for horizontal work-array exchange
    INTEGER(iwp) ::  workarr_lr_exchange_type    !< type definition for work-array exchange on left and right boundaries
    INTEGER(iwp) ::  workarr_sn_exchange_type    !< type definition for work-array exchange on south and north boundaries
    INTEGER(iwp) ::  workarr_bt_exchange_type_x  !< type definition for work-array exchange on bottom and top boundary between left-right
                                                 !< neighbouring subdomains
    INTEGER(iwp) ::  workarr_bt_exchange_type_y  !< type definition for work-array exchange on bottom and top boundary between south-north
                                                 !< neighbouring subdomains

    INTEGER(iwp), DIMENSION(3) ::  parent_grid_info_int  !< array for communicating the parent-array dimensions to its children

    REAL(wp), DIMENSION(7) ::  face_area              !< Surface area of each boundary face
    REAL(wp), DIMENSION(7) ::  parent_grid_info_real  !< Array for communicating the real-type parent-grid parameters to its
                                                      !< children.

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  workarr_lr  !< work array for interpolation on left and right boundaries
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  workarr_sn  !< work array for interpolation on south and north boundaries
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  workarr_bt  !< work array for interpolation on bottom and top boundary

    TYPE parentgrid_def
       INTEGER(iwp) ::  nx  !<
       INTEGER(iwp) ::  ny  !<
       INTEGER(iwp) ::  nz  !<
       REAL(wp) ::  dx                  !<
       REAL(wp) ::  dy                  !<
       REAL(wp) ::  dz                  !<
       REAL(wp) ::  lower_left_coord_x  !<
       REAL(wp) ::  lower_left_coord_y  !<
       REAL(wp) ::  xend                !<
       REAL(wp) ::  yend                !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  coord_x  !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  coord_y  !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dzu      !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  dzw      !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu       !<
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw       !<
    END TYPE parentgrid_def

    TYPE(parentgrid_def), SAVE, PUBLIC ::  pg  !< Parent-grid information package of type parentgrid_def
!
!-- Variables for particle coupling
    TYPE, PUBLIC :: childgrid_def
       INTEGER(iwp) ::  nx          !<
       INTEGER(iwp) ::  ny          !<
       INTEGER(iwp) ::  nz          !<
       INTEGER(iwp) ::  ks          !<
       INTEGER(iwp) ::  ke          !<
       REAL(wp)     ::  dx          !<
       REAL(wp)     ::  dy          !<
       REAL(wp)     ::  dz          !<
       REAL(wp)     ::  lx_coord    !<
       REAL(wp)     ::  rx_coord    !<
       REAL(wp)     ::  sy_coord    !<
       REAL(wp)     ::  ny_coord    !<
       REAL(wp)     ::  uz_coord    !<
       REAL(wp)     ::  lx_coord_b  !<
       REAL(wp)     ::  rx_coord_b  !<
       REAL(wp)     ::  sy_coord_b  !<
       REAL(wp)     ::  ny_coord_b  !<
       REAL(wp)     ::  uz_coord_b  !<
    END TYPE childgrid_def

    TYPE(childgrid_def), SAVE, DIMENSION(:), ALLOCATABLE ::  childgrid  !<

    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE, TARGET ::  nr_part   !<
    INTEGER(idp), DIMENSION(:,:), ALLOCATABLE, TARGET ::  part_adr  !<


    INTERFACE get_childid
       MODULE PROCEDURE get_childid
    END  INTERFACE get_childid

    INTERFACE get_child_edges
       MODULE PROCEDURE get_child_edges
    END  INTERFACE get_child_edges

    INTERFACE get_child_gridspacing
       MODULE PROCEDURE get_child_gridspacing
    END  INTERFACE get_child_gridspacing

    INTERFACE get_number_of_children
       MODULE PROCEDURE get_number_of_children
    END  INTERFACE get_number_of_children

    INTERFACE pmci_adjust_dt_coupling
       MODULE PROCEDURE pmci_adjust_dt_coupling
    END INTERFACE pmci_adjust_dt_coupling

    INTERFACE pmci_atmos_ocean
       MODULE PROCEDURE pmci_atmos_ocean
    END INTERFACE pmci_atmos_ocean

    INTERFACE pmci_boundary_conds
       MODULE PROCEDURE pmci_boundary_conds
    END INTERFACE pmci_boundary_conds

    INTERFACE pmci_check_setting_mismatches
       MODULE PROCEDURE pmci_check_setting_mismatches
    END INTERFACE

    INTERFACE pmci_child_initialize
       MODULE PROCEDURE pmci_child_initialize
    END INTERFACE

    INTERFACE pmci_datatrans
       MODULE PROCEDURE pmci_datatrans
    END INTERFACE pmci_datatrans

    INTERFACE pmci_ensure_nest_mass_conservation
       MODULE PROCEDURE pmci_ensure_nest_mass_conservation
    END INTERFACE pmci_ensure_nest_mass_conservation

    INTERFACE pmci_finalize
       MODULE PROCEDURE pmci_finalize
    END INTERFACE pmci_finalize

    INTERFACE pmci_init
       MODULE PROCEDURE pmci_init
    END INTERFACE

    INTERFACE pmci_modelconfiguration
       MODULE PROCEDURE pmci_modelconfiguration
    END INTERFACE

    INTERFACE pmci_parent_initialize
       MODULE PROCEDURE pmci_parent_initialize
    END INTERFACE

    INTERFACE pmci_set_swaplevel
       MODULE PROCEDURE pmci_set_swaplevel
    END INTERFACE pmci_set_swaplevel

    INTERFACE pmci_synchronize
       MODULE PROCEDURE pmci_synchronize
    END INTERFACE

    PUBLIC atmosphere_ocean_coupled_run,                                                           &
           childgrid,                                                                              &
           child_to_parent,                                                                        &
           comm_world_nesting,                                                                     &
           cpl_id,                                                                                 &
           homogeneous_initialization_child,                                                       &
           nesting_bounds,                                                                         &
           nesting_datatransfer_mode,                                                              &
           nesting_mode,                                                                           &
           nest_shift_z,                                                                           &
           nr_part,                                                                                &
           part_adr,                                                                               &
           parent_to_child,                                                                        &
           particle_coupling,                                                                      &
           rans_mode_parent,                                                                       &
           root_model

    PUBLIC get_childid
    PUBLIC get_child_edges
    PUBLIC get_child_gridspacing
    PUBLIC get_number_of_children
    PUBLIC pmci_adjust_dt_coupling
    PUBLIC pmci_atmos_ocean
    PUBLIC pmci_boundary_conds
    PUBLIC pmci_child_initialize
    PUBLIC pmci_datatrans
    PUBLIC pmci_ensure_nest_mass_conservation
    PUBLIC pmci_finalize
    PUBLIC pmci_init
    PUBLIC pmci_modelconfiguration
    PUBLIC pmci_parent_initialize
    PUBLIC pmci_synchronize
    PUBLIC pmci_set_swaplevel
    PUBLIC pmc_get_model_info

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Find out if this is a nested run and if so, read and broadcast the nesting parameters and set
!> the communicators accordingly.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_init( world_comm )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(OUT) ::  world_comm  !< global communicator


#if defined( __parallel )

    INTEGER(iwp) ::  pmc_status  !< status parameter indicating if the nesting_parameters namelist
                                 !< was succesfully input or not


    CALL pmc_init_model( world_comm, nesting_bounds, nesting_datatransfer_mode, nesting_mode,      &
                         anterpolation_buffer_width, anterpolation_starting_height,                &
                         homogeneous_initialization_child, particle_coupling, pmc_status )
    IF ( pmc_status == pmc_no_namelist_found )  THEN
!
!--    This is not a nested run
       world_comm = MPI_COMM_WORLD
       cpl_id     = 1
       cpl_name   = ''

       RETURN

    ENDIF
!
!-- Check steering parameter values
    IF ( TRIM( nesting_bounds ) /= '3d_nested'  .AND.                                              &
         TRIM( nesting_bounds ) /= 'cyclic_along_x'  .AND.                                         &
         TRIM( nesting_bounds ) /= 'cyclic_along_y'  .AND.                                         &
         TRIM( nesting_bounds ) /= 'vertical_only' )                                               &
    THEN
       message_string = 'illegal nesting boundary condition: ' // TRIM( nesting_bounds )
       CALL message( 'pmci_init', 'PMC0005', 3, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( nesting_mode ) /= 'one-way'  .AND.  TRIM( nesting_mode ) /= 'two-way' )  THEN
       message_string = 'illegal nesting mode: ' // TRIM( nesting_mode )
       CALL message( 'pmci_init', 'PMC0006', 3, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( nesting_datatransfer_mode ) /= 'cascade'  .AND.                                     &
         TRIM( nesting_datatransfer_mode ) /= 'mixed'    .AND.                                     &
         TRIM( nesting_datatransfer_mode ) /= 'overlap' )                                          &
    THEN
       message_string = 'illegal nesting datatransfer mode: ' // TRIM( nesting_datatransfer_mode )
       CALL message( 'pmci_init', 'PMC0007', 3, 2, 0, 6, 0 )
    ENDIF
!
!-- Set the steering switch for a pure vertical nesting.
    IF ( TRIM( nesting_bounds ) == 'vertical_only' )  nesting_bounds_vertical_only = .TRUE.
!
!-- Get some variables required by the pmc-interface (and in some cases in the PALM code out of the
!-- pmci) out of the pmc-core
    CALL pmc_get_model_info( comm_world_nesting = comm_world_nesting, cpl_id = cpl_id,             &
                             cpl_parent_id = cpl_parent_id, cpl_name = cpl_name,                   &
                             npe_total = cpl_npe_total, lower_left_x = lower_left_coord_x,         &
                             lower_left_y = lower_left_coord_y, nest_shift_z = nest_shift_z,       &
                             parent_nest_shift_z = parent_nest_shift_z,                            &
                             atmosphere_ocean_coupled_run = atmosphere_ocean_coupled_run,          &
                             root_model = root_model )
!
!-- Set the steering switch which tells the models that they are nested (of course the root domain
!-- is not nested)
    IF ( .NOT. root_model )  THEN
       IF ( cpl_name(1:5) == 'ocean' )  THEN
           coupling_char = '_O'
       ELSE
          child_domain = .TRUE.
          WRITE( coupling_char, '(A2,I2.2)') '_N', cpl_id
       ENDIF
    ENDIF

!
!-- Set the general steering switch which tells PALM that it is a nested or a coupled run.
    IF ( atmosphere_ocean_coupled_run )  THEN
       IF ( root_model )  THEN
          atmosphere_run_coupled_to_ocean = .TRUE.
       ELSE
          ocean_run_coupled_to_atmosphere = .TRUE.
       ENDIF
    ELSE
       nested_run = .TRUE.
    ENDIF

!
!-- Message that communicators for nesting are initialized.
!-- Attention: myid has been set at the end of pmc_init_model in order to guarantee that only PE0 of
!-- the root domain does the output.
    CALL location_message( 'initialize model nesting', 'finished' )
!
!-- Reset myid to its default value
    myid = 0
#else
!
!-- Nesting cannot be used in serial mode. cpl_id is set to root domain (1) because no location
!-- messages would be generated otherwise. world_comm is given a dummy value to avoid compiler
!-- warnings (INTENT(OUT) must get an explicit value).
!-- Note that this branch is only to avoid compiler warnings. The actual execution never reaches
!-- here because the call of this subroutine is already enclosed by  #if defined( __parallel ).
    cpl_id     = 1
    nested_run = .FALSE.
    world_comm = 1
#endif

 END SUBROUTINE pmci_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define the nesting setup.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_modelconfiguration

    IMPLICIT NONE

    INTEGER(iwp) ::  ncpl  !< number of nest domains


#if defined( __parallel )
    CALL location_message( 'setup the nested model configuration', 'start' )
    CALL cpu_log( log_point_s(79), 'pmci_model_config', 'start' )
!
!-- Compute absolute coordinates for all models
    CALL pmci_setup_coordinates
!
!-- Allocate memory for 2-D surface arrays to be coupled
    CALL surface_coupler_alloc_mem
!
!-- Determine the number of coupled arrays
    CALL pmci_num_arrays
!
!-- Initialize the child (must be called before pmc_setup_parent)
!-- Klaus, extend this comment to explain why it must be called before
    CALL pmci_setup_child
!
!-- Initialize PMC parent
    CALL pmci_setup_parent
!
!-- Check for mismatches between settings of master and child variables
!-- (e.g., all children have to follow the end_time settings of the root master)
    CALL pmci_check_setting_mismatches
!
!-- Set flag file for combine_plot_fields for processing the nest / coupling output data
    IF ( myid == 0 )  THEN
       IF ( nested_run )  THEN
          OPEN( 90, FILE = '3DNESTING', FORM = 'FORMATTED' )
          CALL pmc_get_model_info( ncpl = ncpl )
          WRITE( 90, '(I2)' )  ncpl
          CLOSE( 90 )
       ELSEIF( atmosphere_ocean_coupled_run )  THEN
          OPEN( 90, FILE = 'COUPLED_ATMOSPHERE_OCEAN', FORM = 'FORMATTED' )
          WRITE( 90, '(A)' )  'this is a coupled atmosphere-ocean-run'
          CLOSE( 90 )
       ENDIF
    ENDIF

    CALL cpu_log( log_point_s(79), 'pmci_model_config', 'stop' )
    CALL location_message( 'setup the nested model configuration', 'finished' )
!
!-- Check for invalid combinations
    IF ( nested_run  .AND.  atmosphere_ocean_coupled_run )  THEN
       message_string = 'combination of nesting and atmosphere-ocean coupling is not allowed'
       CALL message( 'pmci_modelconfiguration', 'PMC0008', 1, 2, 0, 6, 0 )
    ENDIF
#endif

 END SUBROUTINE pmci_modelconfiguration


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Prepare the coupling environment for the current parent domain.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_setup_parent

#if defined( __parallel )
    IMPLICIT NONE

    CHARACTER(LEN=32) ::  myname  !< String for variable name such as 'u'

    INTEGER(iwp) ::  child_id    !< Child id-number for the child m
    INTEGER(iwp) ::  ierr        !< MPI-error code
    INTEGER(iwp) ::  kp          !< Parent-grid index n the z-direction
    INTEGER(iwp) ::  lb = 1      !< Running index for aerosol size bins
    INTEGER(iwp) ::  lc = 1      !< Running index for aerosol mass bins
    INTEGER(iwp) ::  lg = 1      !< Running index for SALSA gases
    INTEGER(iwp) ::  m           !< Loop index over all children of the current parent
    INTEGER(iwp) ::  msib        !< Loop index over all other children than m in case of siblings (parallel children)
    INTEGER(iwp) ::  n = 1       !< Running index for chemical species
    INTEGER(iwp) ::  nx_child    !< Number of child-grid points in the x-direction
    INTEGER(iwp) ::  ny_child    !< Number of child-grid points in the y-direction
    INTEGER(iwp) ::  kz_child_start  !< parent k value of the lower level of the child domain  ! AH: change to kp_child_bottom?
    INTEGER(iwp) ::  kz_child_end    !< parent k value of the upper level of the child domain  ! AH: change to kp_child_top?
    INTEGER(iwp) ::  nz_child    !< Number of child-grid points in the z-direction
    INTEGER(iwp) ::  sibling_id  !< Child id-number for the child msib (sibling of child m)

    INTEGER(iwp), DIMENSION(3) ::  child_grid_dim  !< Array for receiving the child-grid dimensions from the children

    LOGICAL ::  child_matches_parent_along_x  !< switch that tells if child matches parent along x
    LOGICAL ::  child_matches_parent_along_y  !< switch that tells if child matches parent along y
    LOGICAL ::  m_left_in_msib   !< Logical auxiliary parameter for the overlap test: true if the left border
                                 !< of the child m is within the x-range of the child msib
    LOGICAL ::  m_right_in_msib  !< Logical auxiliary parameter for the overlap test: true if the right border
                                 !< of the child m is within the x-range of the child msib
    LOGICAL ::  msib_left_in_m   !< Logical auxiliary parameter for the overlap test: true if the left border
                                 !< of the child msib is within the x-range of the child m
    LOGICAL ::  msib_right_in_m  !< Logical auxiliary parameter for the overlap test: true if the right border
                                 !< of the child msib is within the x-range of the child m
    LOGICAL ::  m_south_in_msib  !< Logical auxiliary parameter for the overlap test: true if the south border
                                 !< of the child m is within the y-range of the child msib
    LOGICAL ::  m_north_in_msib  !< Logical auxiliary parameter for the overlap test: true if the north border
                                 !< of the child m is within the y-range of the child msib
    LOGICAL ::  msib_south_in_m  !< Logical auxiliary parameter for the overlap test: true if the south border
                                 !< of the child msib is within the y-range of the child m
    LOGICAL ::  msib_north_in_m  !< Logical auxiliary parameter for the overlap test: true if the north border
                                 !< of the child msib is within the y-range of the child m

    REAL(wp) ::  child_height         !< Height of the child domain defined on the child side as zw(nzt+1)
    REAL(wp) ::  child_low            !< Height of the child domain at the bottom defined on the child side as zw(nzb)
    REAL(wp) ::  dx_child             !< Child-grid spacing in the x-direction
    REAL(wp) ::  dy_child             !< Child-grid spacing in the y-direction
    REAL(wp) ::  dz_child             !< Child-grid spacing in the z-direction
    REAL(wp) ::  left_limit           !< Left limit for the absolute x-coordinate of the child left boundary
    REAL(wp) ::  north_limit          !< North limit for the absolute y-coordinate of the child north boundary
    REAL(wp) ::  right_limit          !< Right limit for the absolute x-coordinate of the child right boundary
    REAL(wp) ::  south_limit          !< South limit for the absolute y-coordinate of the child south boundary
    REAL(wp) ::  upper_right_coord_x  !< Absolute x-coordinate of the upper right corner of the child domain
    REAL(wp) ::  upper_right_coord_y  !< Absolute y-coordinate of the upper right corner of the child domain
    REAL(wp) ::  xez                  !< Minimum separation in the x-direction required between the child and
                                      !< parent boundaries (left or right)
    REAL(wp) ::  yez                  !< Minimum separation in the y-direction required between the child and
                                      !< parent boundaries (south or north)
    REAL(wp) ::  tolex                !< Tolerance for grid-line matching in x-direction
    REAL(wp) ::  toley                !< Tolerance for grid-line matching in y-direction
    REAL(wp) ::  tolez                !< Tolerance for grid-line matching in z-direction

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_coord_x  !< Child domain x-coordinate array
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_coord_y  !< Child domain y-coordinate array
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_x_left   !< Minimum x-coordinate of the child domain including the ghost
                                                           !< point layers
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_x_right  !< Maximum x-coordinate of the child domain including the ghost
                                                           !< point layers
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_y_north  !< Maximum y-coordinate of the child domain including the ghost
                                                           !< point layers
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  child_y_south  !< Minimum y-coordinate of the child domain including the ghost
                                                           !< point layers

    REAL(wp), DIMENSION(6) ::  child_grid_info             !< Array for receiving the child-grid spacings etc from the children

!
!-- Grid-line tolerances.
    tolex = tolefac * dx
    toley = tolefac * dy
    tolez = tolefac * dz(1)
!
!-- Initialize the current pmc parent.
    CALL pmc_parentinit
!
!-- Corners of all children of the present parent. Note that SIZE( pmc_parent_for_child ) = 1 if we
!-- have no children.
    IF ( ( SIZE( pmc_parent_for_child ) - 1 > 0 )  .AND.  myid == 0 )  THEN
       ALLOCATE( child_x_left(1:SIZE(  pmc_parent_for_child ) - 1) )
       ALLOCATE( child_x_right(1:SIZE( pmc_parent_for_child ) - 1) )
       ALLOCATE( child_y_south(1:SIZE( pmc_parent_for_child ) - 1) )
       ALLOCATE( child_y_north(1:SIZE( pmc_parent_for_child ) - 1) )
    ENDIF
    IF ( ( SIZE( pmc_parent_for_child ) - 1 > 0 ) )  THEN
       ALLOCATE( childgrid(1:SIZE( pmc_parent_for_child ) - 1) )
    ENDIF
!
!-- Get coordinates from all children and check that the children match the parent domain and each
!-- others. Note that SIZE( pmc_parent_for_child ) = 1 if we have no children, hence the loop is
!-- not executed at all.
    DO  m = 1, SIZE( pmc_parent_for_child ) - 1

       child_id = pmc_parent_for_child(m)
!
!--    Set counter variables for chemical and aerosol species back to one for each child domain
       n  = 1
       lb = 1
       lc = 1
       lg = 1

       IF ( myid == 0 )  THEN

          CALL pmc_recv_from_child( child_id, child_grid_dim,  SIZE( child_grid_dim ), 0, 123,     &
                                    ierr )
          CALL pmc_recv_from_child( child_id, child_grid_info, SIZE( child_grid_info ), 0, 124,    &
                                    ierr )

          nx_child     = child_grid_dim(1)
          ny_child     = child_grid_dim(2)
          dx_child     = child_grid_info(3)
          dy_child     = child_grid_info(4)
          dz_child     = child_grid_info(5)
          child_height = child_grid_info(1)
          child_low    = child_grid_info(6)

          IF ( .NOT. ocean_mode )  THEN
!
!--          Find the lowest child-domain level in the parent grid for the reduced z transfer
             DO  kp = 1, nzt
                IF ( zw(kp) - child_low > tolez )  THEN
                   kz_child_start = kp
                   EXIT
                ENDIF
             ENDDO
!
!--          Find the highest child-domain level in the parent grid for the reduced z transfer
             DO  kp = 1, nzt
                IF ( zw(kp) - child_height > tolez )  THEN
                   kz_child_end = kp
                   EXIT
                ENDIF
             ENDDO

          ELSE
!
!--          Find the lowest child-domain level in the parent grid for the reduced z transfer
             DO  kp = 1, nzt
                IF ( zw(kp) - child_height > tolez )  THEN
                   kz_child_start = kp
                   EXIT
                ENDIF
             ENDDO
!
!--          In ocean mode the top of the child is always positioned at the top of the parent (the
!--          ocean surface).
             kz_child_end = nzt

          ENDIF

          nz_child = kz_child_end - kz_child_start + 1

!
!--       Get absolute coordinates from the child
          ALLOCATE( child_coord_x(-nbgp:nx_child+nbgp) )
          ALLOCATE( child_coord_y(-nbgp:ny_child+nbgp) )

          CALL pmc_recv_from_child( child_id, child_coord_x, SIZE( child_coord_x ), 0, 11, ierr )
          CALL pmc_recv_from_child( child_id, child_coord_y, SIZE( child_coord_y ), 0, 12, ierr )

          parent_grid_info_real(1) = lower_left_coord_x
          parent_grid_info_real(2) = lower_left_coord_y
          parent_grid_info_real(3) = dx
          parent_grid_info_real(4) = dy

          upper_right_coord_x      = lower_left_coord_x + ( nx + 1 ) * dx
          upper_right_coord_y      = lower_left_coord_y + ( ny + 1 ) * dy
          parent_grid_info_real(5) = upper_right_coord_x
          parent_grid_info_real(6) = upper_right_coord_y
          parent_grid_info_real(7) = dz(1)

          parent_grid_info_int(1)  = nx
          parent_grid_info_int(2)  = ny
          parent_grid_info_int(3)  = nz_child

!
!--       Check, if lateral boundaries of parent and child match exactly
          IF ( ( ABS( child_coord_x(0)          - lower_left_coord_x  ) > tolex )  .OR.            &
               ( ABS( child_coord_x(nx_child+1) - upper_right_coord_x ) > tolex ) )                &
          THEN
             child_matches_parent_along_x = .FALSE.
          ELSE
             child_matches_parent_along_x = .TRUE.
          ENDIF
          IF ( ( ABS( child_coord_y(0)          - lower_left_coord_y  ) > toley )  .OR.            &
               ( ABS( child_coord_y(ny_child+1) - upper_right_coord_y ) > toley ) )                &
          THEN
             child_matches_parent_along_y = .FALSE.
          ELSE
             child_matches_parent_along_y = .TRUE.
          ENDIF
!
!--       Depending on the nesting-/coupling mode, check that the child domain matches its parent
!--       domain respectively.
          IF ( nesting_bounds_vertical_only  .OR.  atmosphere_ocean_coupled_run )  THEN
!
!--          In case of vertical nesting or atmosphere-ocean coupling, the lateral boundaries must
!--          match exactly.
             IF ( .NOT. child_matches_parent_along_x )  THEN
                WRITE( message_string, '(A,I2,A)' ) 'nested child (id: ', child_id,                &
                       ') domain boundaries along x do not match the respective parent boundaries'
                CALL message( 'pmci_setup_parent', 'PMC0009', 3, 2, 0, 6, 0 )
             ENDIF
             IF ( .NOT. child_matches_parent_along_y )  THEN
                WRITE( message_string, '(A,I2,A)' ) 'nested child (id: ', child_id,                &
                       ') domain boundaries along y do not match the respective parent boundaries'
                CALL message( 'pmci_setup_parent', 'PMC0009', 3, 2, 0, 6, 0 )
             ENDIF

!
!--          In case of coupling, due to interpolation requirements, the total number of ocean grid
!--          points along x or y must be a multiple of the total number of atmosphere grid points.
             IF ( atmosphere_ocean_coupled_run ) THEN

                IF ( MOD( nx_child+1, nx+1 ) /= 0 )  THEN
                   message_string = 'nx+1 in ocean is not divisible by nx+1 in atmosphere ' //     &
                                    'without remainder'
                   CALL message( 'pmci_setup_parent', 'PMC0010', 1, 2, 0, 6, 0 )
                ENDIF

                IF ( MOD( ny_child+1, ny+1 ) /= 0 )  THEN
                   message_string = 'ny+1 in ocean is not divisible by ny+1 in atmosphere ' //     &
                                    'without remainder'
                   CALL message( 'pmci_setup_parent', 'PMC0011', 1, 2, 0, 6, 0 )
                ENDIF

             ENDIF

          ELSE
!
!--          In case of 3-D nesting, check that the child domain is completely inside its parent
!--          domain. This check  also cares for the case, where the user set cyclic conditions for
!--          the child, but didn't care for that the child has the same extension as the parent
!--          along that direction.
             IF ( nesting_bounds == 'cyclic_along_x'  .AND. .NOT. child_matches_parent_along_x )   &
             THEN
!
!--             Lateral boundaries along x must match exactly.
                IF ( .NOT. child_matches_parent_along_x )  THEN
                   WRITE( message_string, '(A,I2,A)' ) 'nested child (id: ', child_id,             &
                         ') domain boundaries along x do not match the respective parent boundaries'
                   CALL message( 'pmci_setup_parent', 'PMC0009', 3, 2, 0, 6, 0 )
                ENDIF

             ELSEIF ( nesting_bounds /= 'cyclic_along_x' )  THEN

                xez = ( nbgp + 1 ) * dx
                left_limit  = lower_left_coord_x + xez
                right_limit = upper_right_coord_x - xez
                IF ( left_limit - child_coord_x(0) > tolex )  THEN
                   WRITE( message_string, '(A,I2,A)' ) 'nested child (id: ',child_id,              &
                          ') domain does not fit in its parent domain, left edge is either ' //    &
                          'too close or outside its parent left edge'
                   CALL message( 'pmci_setup_parent', 'PMC0012', 3, 2, 0, 6, 0 )
                ENDIF
                IF ( child_coord_x(nx_child+1) - right_limit > tolex )  THEN
                   WRITE( message_string, '(A,I2,A)' ) 'nested child (id: ',child_id,              &
                          ') domain does not fit in its parent domain, right edge is either ' //   &
                          'too close or outside its parent right edge'
                   CALL message( 'pmci_setup_parent', 'PMC0012', 3, 2, 0, 6, 0 )
                ENDIF

             ENDIF

             IF ( nesting_bounds == 'cyclic_along_y'  .AND. .NOT. child_matches_parent_along_y )   &
             THEN
!
!--             Lateral boundaries along y must match exactly.
                IF ( .NOT. child_matches_parent_along_y )  THEN
                   WRITE( message_string, '(A,I2,A)' ) 'nested child (id: ', child_id,             &
                         ') domain boundaries along y do not match the respective parent boundaries'
                   CALL message( 'pmci_setup_parent', 'PMC0009', 3, 2, 0, 6, 0 )
                ENDIF

             ELSEIF ( nesting_bounds /= 'cyclic_along_y' )  THEN

                yez = ( nbgp + 1 ) * dy
                south_limit = lower_left_coord_y + yez
                north_limit = upper_right_coord_y - yez
                IF ( south_limit - child_coord_y(0) > toley )  THEN
                   WRITE( message_string, '(A,I2,A)' ) 'nested child (id: ',child_id,              &
                          ') domain does not fit in its parent domain, south edge is either ' //   &
                          'too close or outside its parent south edge'
                   CALL message( 'pmci_setup_parent', 'PMC0012', 3, 2, 0, 6, 0 )
                ENDIF
                IF ( child_coord_y(ny_child+1) - north_limit > toley )  THEN
                   WRITE( message_string, '(A,I2,A)' ) 'nested child (id: ',child_id,              &
                          ') domain does not fit in its parent domain, north edge is either ' //   &
                          'too close or outside its parent north edge'
                   CALL message( 'pmci_setup_parent', 'PMC0012', 3, 2, 0, 6, 0 )
                ENDIF
             ENDIF

          ENDIF
!
!--       Child domain must be lower than the parent domain such that the top ghost layer of the
!--       child grid does not exceed the parent domain top boundary.
          IF ( nested_run .AND. child_height - zw(nzt) > tolez )  THEN
             WRITE( message_string, '(A,i2,a)' ) 'nested child (id: ',child_id,                    &
                    ') domain does not fit in its parent domain, top edge is either too ' //       &
                    'close or above its parent top edge'
             CALL message( 'pmci_setup_parent', 'PMC0012', 3, 2, 0, 6, 0 )
          ENDIF
!
!--       If parallel child domains (siblings) do exist ( m > 1 ), check that they do not overlap.
          child_x_left(m)  = child_coord_x(-nbgp)
          child_x_right(m) = child_coord_x(nx_child+nbgp)
          child_y_south(m) = child_coord_y(-nbgp)
          child_y_north(m) = child_coord_y(ny_child+nbgp)

          IF ( .NOT. nesting_bounds_vertical_only  .AND.  .NOT. atmosphere_ocean_coupled_run )  THEN
!
!--          Note that the msib-loop is executed only if ( m > 1 ).
!--          Also note that the tests have to be done both ways (m vs msib and msib vs m) in order
!--          to detect all the possible overlap situations.
             DO  msib = 1, m - 1
!
!--             Set some logical auxiliary parameters to simplify the IF-condition.
                m_left_in_msib  = ( child_x_left(m)  >= child_x_left(msib)  - tolex )  .AND.       &
                                  ( child_x_left(m)  <= child_x_right(msib) + tolex )
                m_right_in_msib = ( child_x_right(m) >= child_x_left(msib)  - tolex )  .AND.       &
                                  ( child_x_right(m) <= child_x_right(msib) + tolex )
                msib_left_in_m  = ( child_x_left(msib)  >= child_x_left(m)  - tolex )  .AND.       &
                                  ( child_x_left(msib)  <= child_x_right(m) + tolex )
                msib_right_in_m = ( child_x_right(msib) >= child_x_left(m)  - tolex )  .AND.       &
                                  ( child_x_right(msib) <= child_x_right(m) + tolex )
                m_south_in_msib = ( child_y_south(m) >= child_y_south(msib) - toley )  .AND.       &
                                  ( child_y_south(m) <= child_y_north(msib) + toley )
                m_north_in_msib = ( child_y_north(m) >= child_y_south(msib) - toley )  .AND.       &
                                  ( child_y_north(m) <= child_y_north(msib) + toley )
                msib_south_in_m = ( child_y_south(msib) >= child_y_south(m) - toley )  .AND.       &
                                  ( child_y_south(msib) <= child_y_north(m) + toley )
                msib_north_in_m = ( child_y_north(msib) >= child_y_south(m) - toley )  .AND.       &
                                  ( child_y_north(msib) <= child_y_north(m) + toley )

                IF ( ( m_left_in_msib  .OR.  m_right_in_msib  .OR.                                 &
                       msib_left_in_m  .OR.  msib_right_in_m )  .AND.                              &
                     ( m_south_in_msib  .OR.  m_north_in_msib  .OR.                                &
                       msib_south_in_m  .OR.  msib_north_in_m ) )  THEN
                   sibling_id = pmc_parent_for_child(msib)
                   WRITE( message_string, '(A,I2,A,I2,A)' ) 'nested parallel child domains (ids: ',&
                          child_id, ' and ', sibling_id, ') overlap'
                   CALL message( 'pmci_setup_parent', 'PMC0013', 3, 2, 0, 6, 0 )
                ENDIF

             ENDDO
          ENDIF

          CALL pmci_set_child_edge_coords

          DEALLOCATE( child_coord_x )
          DEALLOCATE( child_coord_y )
!
!--       Send information about operating mode (LES or RANS) to child. This will be used to
!--       control TKE nesting and setting boundary conditions properly.
          CALL pmc_send_to_child( child_id, rans_mode, 1, 0, 19, ierr )
!
!--       Send parent grid information to child
          CALL pmc_send_to_child( child_id, parent_grid_info_real, SIZE( parent_grid_info_real ),  &
                                  0, 21, ierr )
          CALL pmc_send_to_child( child_id, parent_grid_info_int,  SIZE( parent_grid_info_int ),   &
                                  0, 22, ierr )
!
!--       Send local grid to child
          CALL pmc_send_to_child( child_id, coord_x, nx+1+2*nbgp, 0, 24, ierr )
          CALL pmc_send_to_child( child_id, coord_y, ny+1+2*nbgp, 0, 25, ierr )
!
!--       Also send the dzu-, dzw-, zu- and zw-arrays here
          CALL pmc_send_to_child( child_id, dzu, nz_child + 1, 0, 26, ierr )
          CALL pmc_send_to_child( child_id, dzw, nz_child + 1, 0, 27, ierr )
          CALL pmc_send_to_child( child_id, zu(kz_child_start-1:kz_child_end+1),  nz_child + 2, 0, 28, ierr )
          CALL pmc_send_to_child( child_id, zw(kz_child_start-1:kz_child_end+1),  nz_child + 2, 0, 29, ierr )
!
!--       Required in case of atmosphere-ocean coupling.
          CALL pmc_send_to_child( child_id, humidity,  1, 0, 30, ierr )

       ENDIF  ! ( myid == 0 )

       CALL MPI_BCAST( kz_child_start, 1, MPI_INTEGER, 0, comm2d, ierr )
       CALL MPI_BCAST( kz_child_end, 1, MPI_INTEGER, 0, comm2d, ierr )
       CALL MPI_BCAST( nz_child, 1, MPI_INTEGER, 0, comm2d, ierr )

       CALL MPI_BCAST( childgrid(m), STORAGE_SIZE( childgrid( 1 ) ) / 8, MPI_BYTE, 0, comm2d, ierr )

!
!--    Set up the index-list which is an integer array that maps the child index space on the parent
!--    index- and subdomain spaces.
       CALL pmci_create_index_list
!
!--    Include couple arrays into parent content.
!--    The adresses of the PALM 2D or 3D array (here parent grid) which are candidates for coupling
!--    are stored once into the pmc context. While data transfer, the arrays do not have to be
!--    specified again
       CALL pmc_s_clear_next_array_list
       DO WHILE ( pmc_s_getnextarray( child_id, myname ) )
          IF ( INDEX( TRIM( myname ), 'chem_' ) /= 0 )  THEN
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child,        &
                                          ks_child = kz_child_start, ke_child = kz_child_end,      &
                                          n = n )
             n = n + 1
          ELSEIF ( INDEX( TRIM( myname ), 'an_' ) /= 0 )  THEN
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child,        &
                                          ks_child = kz_child_start, ke_child = kz_child_end,      &
                                          n = lb )
             lb = lb + 1
          ELSEIF ( INDEX( TRIM( myname ), 'am_' ) /= 0 )  THEN
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child,        &
                                          ks_child = kz_child_start, ke_child = kz_child_end,      &
                                          n = lc )
             lc = lc + 1
          ELSEIF ( INDEX( TRIM( myname ), 'sg_' ) /= 0  .AND.  .NOT. salsa_gases_from_chem )  THEN
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child,        &
                                          ks_child = kz_child_start, ke_child = kz_child_end,      &
                                          n = lg )
             lg = lg + 1
          ELSE
             CALL pmci_set_array_pointer( myname, child_id = child_id, nz_child = nz_child,        &
                                          ks_child = kz_child_start, ke_child = kz_child_end )
          ENDIF
       ENDDO

       CALL pmc_s_setind_and_allocmem( child_id )

    ENDDO  ! m

    IF ( ( SIZE( pmc_parent_for_child ) - 1 > 0 )  .AND.  myid == 0 )  THEN
       DEALLOCATE( child_x_left )
       DEALLOCATE( child_x_right )
       DEALLOCATE( child_y_south )
       DEALLOCATE( child_y_north )
    ENDIF

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create the index list
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_create_index_list

    IMPLICIT NONE

    INTEGER(iwp) ::  ilist            !< Index-list index running over the child's parent-grid jc,ic-space
    INTEGER(iwp) ::  index_list_size  !< Dimension 2 of the array index_list
    INTEGER(iwp) ::  ierr             !< MPI error code
    INTEGER(iwp) ::  ip               !< Running parent-grid index on the child domain in the x-direction
    INTEGER(iwp) ::  jp               !< Running parent-grid index on the child domain in the y-direction
    INTEGER(iwp) ::  n                !< Running index over child subdomains
    INTEGER(iwp) ::  nrx              !< Parent subdomain dimension in the x-direction
    INTEGER(iwp) ::  nry              !< Parent subdomain dimension in the y-direction
    INTEGER(iwp) ::  pex              !< Two-dimensional subdomain (pe) index in the x-direction
    INTEGER(iwp) ::  pey              !< Two-dimensional subdomain (pe) index in the y-direction
    INTEGER(iwp) ::  parent_pe        !< Parent subdomain index (one-dimensional)

    INTEGER(iwp), DIMENSION(2) ::  pe_indices_2d                          !< Array for two-dimensional subdomain (pe)
                                                                          !< indices needed for MPI_CART_RANK
    INTEGER(iwp), DIMENSION(2) ::  size_of_childs_parent_grid_bounds_all  !< Dimensions of childs_parent_grid_bounds_all

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE  ::  childs_parent_grid_bounds_all  !< Array that contains the child's
                                                                                  !< parent-grid index
                                                                                  !< bounds for all its subdomains (pes)
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE  ::  index_list                     !< Array that maps the child index space on
                                                                                  !< the parent index- and subdomain spaces

    IF ( myid == 0 )  THEN

       CALL pmc_recv_from_child( child_id, size_of_childs_parent_grid_bounds_all, 2, 0, 40, ierr )
       ALLOCATE( childs_parent_grid_bounds_all(size_of_childs_parent_grid_bounds_all(1),           &
                                               size_of_childs_parent_grid_bounds_all(2)) )
       CALL pmc_recv_from_child( child_id, childs_parent_grid_bounds_all,                          &
                                 SIZE( childs_parent_grid_bounds_all ), 0, 41, ierr )
!
!--    Compute size (dimension) of the index_list.
       index_list_size = 0
       DO  n = 1, size_of_childs_parent_grid_bounds_all(2)
          index_list_size = index_list_size +                                                      &
               ( childs_parent_grid_bounds_all(4,n) - childs_parent_grid_bounds_all(3,n) + 1 ) *   &
               ( childs_parent_grid_bounds_all(2,n) - childs_parent_grid_bounds_all(1,n) + 1 )
       ENDDO

       ALLOCATE( index_list(6,index_list_size) )

       nrx = nxr - nxl + 1
       nry = nyn - nys + 1
       ilist = 0
!
!--    Loop over all children PEs
       DO  n = 1, size_of_childs_parent_grid_bounds_all(2)
!
!--       Subspace along y required by actual child PE
          DO  jp = childs_parent_grid_bounds_all(3,n), childs_parent_grid_bounds_all(4,n)  ! jp = jps, jpn of child PE# n
!
!--          Subspace along x required by actual child PE
             DO  ip = childs_parent_grid_bounds_all(1,n), childs_parent_grid_bounds_all(2,n)  ! ip = ipl, ipr of child PE# n

                pex = ip / nrx
                pey = jp / nry
                pe_indices_2d(1) = pex
                pe_indices_2d(2) = pey
                CALL MPI_CART_RANK( comm2d, pe_indices_2d, parent_pe, ierr )

                ilist = ilist + 1
!
!--             First index in parent array  ! TO_DO: Klaus, please explain better
                index_list(1,ilist) = ip - ( pex * nrx ) + 1 + nbgp
!
!--             Second index in parent array  ! TO_DO: Klaus, please explain better
                index_list(2,ilist) = jp - ( pey * nry ) + 1 + nbgp
!
!--             x index of child's parent grid
                index_list(3,ilist) = ip - childs_parent_grid_bounds_all(1,n) + 1
!
!--             y index of child's parent grid
                index_list(4,ilist) = jp - childs_parent_grid_bounds_all(3,n) + 1
!
!--             PE number of child
                index_list(5,ilist) = n - 1
!
!--             PE number of parent
                index_list(6,ilist) = parent_pe

             ENDDO
          ENDDO
       ENDDO
!
!--    TO_DO: Klaus: comment what is done here
       CALL pmc_s_set_2d_index_list( child_id, index_list(:,1:ilist) )

    ELSE
!
!--    TO_DO: Klaus: comment why this dummy allocation is required
       ALLOCATE( index_list(6,1) )
       CALL pmc_s_set_2d_index_list( child_id, index_list )
    ENDIF

    DEALLOCATE( index_list )

 END SUBROUTINE pmci_create_index_list


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store the child-edge coordinates.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_set_child_edge_coords
    IMPLICIT  NONE

    INTEGER(iwp) ::  nbgp_lpm = 1  !< Number of ghost-point layers used for lpm (Klaus, is this correct?)


    nbgp_lpm = MIN( nbgp_lpm, nbgp )

    childgrid(m)%nx = nx_child
    childgrid(m)%ny = ny_child
    childgrid(m)%ks = kz_child_start
    childgrid(m)%ke = kz_child_end
    childgrid(m)%nz = nz_child
    childgrid(m)%dx = dx_child
    childgrid(m)%dy = dy_child
    childgrid(m)%dz = dz_child

    childgrid(m)%lx_coord   = child_coord_x(0)
    childgrid(m)%lx_coord_b = child_coord_x(-nbgp_lpm)
    childgrid(m)%rx_coord   = child_coord_x(nx_child) + dx_child
    childgrid(m)%rx_coord_b = child_coord_x(nx_child+nbgp_lpm) + dx_child
    childgrid(m)%sy_coord   = child_coord_y(0)
    childgrid(m)%sy_coord_b = child_coord_y(-nbgp_lpm)
    childgrid(m)%ny_coord   = child_coord_y(ny_child) + dy_child
    childgrid(m)%ny_coord_b = child_coord_y(ny_child+nbgp_lpm) + dy_child
    childgrid(m)%uz_coord   = child_grid_info(2)
    childgrid(m)%uz_coord_b = child_grid_info(1)

 END SUBROUTINE pmci_set_child_edge_coords

#endif
 END SUBROUTINE pmci_setup_parent


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Prepare the coupling environment for the current child domain.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_setup_child

#if defined( __parallel )
    IMPLICIT NONE

    CHARACTER(LEN=da_namelen) ::  myname      !< Name of the variable to be coupled
    CHARACTER(LEN=5)          ::  salsa_char  !< Name extension for the variable name in case of SALSA variable

    INTEGER(iwp) ::  ierr  !< MPI error code
    INTEGER(iwp) ::  lb    !< Running index for aerosol size bins
    INTEGER(iwp) ::  lc    !< Running index for aerosol mass bins
    INTEGER(iwp) ::  lg    !< Running index for SALSA gases
    INTEGER(iwp) ::  n     !< Running index for number of chemical species

    INTEGER(iwp), DIMENSION(3) ::  child_grid_dim  !< Array for sending the child-grid dimensions to parent

    REAL(wp), DIMENSION(6) ::  child_grid_info  !< Array for sending the child-grid spacings etc to parent

!
!-- Child setup
!-- Root model does not have a parent and is not a child, therefore no child setup on root model
    IF ( .NOT. root_model )  THEN
!
!--    KLaus, add a description here what pmc_childinit does
       CALL pmc_childinit

       IF ( atmosphere_ocean_coupled_run )  THEN

          CALL pmc_set_dataarray_name( 'parent', 'exchange_array_1', 'child', 'exchange_array_1',  &
                                       ierr )
          CALL pmc_set_dataarray_name( 'parent', 'exchange_array_2', 'child', 'exchange_array_2',  &
                                       ierr )
          CALL pmc_set_dataarray_name( 'parent', 'exchange_array_3', 'child', 'exchange_array_3',  &
                                       ierr )
          CALL pmc_set_dataarray_name( 'parent', 'exchange_array_4', 'child', 'exchange_array_4',  &
                                       ierr )

       ELSE

!
!--       The arrays, which actually will be exchanged between child and parent are defined Here AND
!--       ONLY HERE. If a variable is removed, it only has to be removed from here. Please check, if
!--       the arrays are in the list of POSSIBLE exchange arrays in subroutines:
!--       pmci_set_array_pointer (for parent arrays)
!--       pmci_create_childs_parent_grid_arrays (for child's parent-grid arrays)
          CALL pmc_set_dataarray_name( 'parent', 'u', 'child', 'u', ierr )
          CALL pmc_set_dataarray_name( 'parent', 'v', 'child', 'v', ierr )
          CALL pmc_set_dataarray_name( 'parent', 'w', 'child', 'w', ierr )
!
!--       Set data array name for TKE. Please note, nesting of TKE is actually only done if both
!--       parent and child are in LES or in RANS mode. Due to design of model coupler, however, data
!--       array names must be already available at this point, though the control flag whether data
!--       shall be interpolated / anterpolated is not available yet.
          CALL pmc_set_dataarray_name( 'parent', 'e', 'child', 'e', ierr )
!
!--       Nesting of dissipation rate only if both parent and child are in RANS mode and TKE-epsilon
!--       closure is applied. Please see also comment for TKE above.
          CALL pmc_set_dataarray_name( 'parent', 'diss', 'child', 'diss', ierr )

          IF ( .NOT. neutral )  THEN
             CALL pmc_set_dataarray_name( 'parent', 'pt' ,'child', 'pt', ierr )
          ENDIF

          IF ( humidity )  THEN

             CALL pmc_set_dataarray_name( 'parent', 'q', 'child', 'q', ierr )

             IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
                CALL pmc_set_dataarray_name( 'parent', 'qc', 'child', 'qc', ierr )
                CALL pmc_set_dataarray_name( 'parent', 'nc', 'child', 'nc', ierr )
             ENDIF

             IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
                CALL pmc_set_dataarray_name( 'parent', 'qr', 'child', 'qr', ierr )
                CALL pmc_set_dataarray_name( 'parent', 'nr', 'child', 'nr', ierr )
             ENDIF

          ENDIF

       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmc_set_dataarray_name( 'parent', 's', 'child', 's', ierr )
       ENDIF

       IF ( ocean_mode  .AND.  salinity )  THEN
          CALL pmc_set_dataarray_name( 'parent', 'sa', 'child', 'sa', ierr )
       ENDIF

       IF ( particle_advection )  THEN
          CALL pmc_set_dataarray_name( 'parent', 'nr_part', 'child', 'nr_part', ierr )
          CALL pmc_set_dataarray_name( 'parent', 'part_adr', 'child', 'part_adr', ierr )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO n = 1, nspec
             CALL pmc_set_dataarray_name( 'parent', 'chem_' // TRIM( chem_species(n)%name ),       &
                                          'child',  'chem_' // TRIM( chem_species(n)%name ), ierr )
          ENDDO
       ENDIF

       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             WRITE( salsa_char,'(i0)' ) lb
             CALL pmc_set_dataarray_name( 'parent', 'an_' // TRIM( salsa_char ),                   &
                                          'child',  'an_' // TRIM( salsa_char ), ierr )
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             WRITE( salsa_char,'(i0)' ) lc
             CALL pmc_set_dataarray_name( 'parent', 'am_' // TRIM( salsa_char ),                   &
                                          'child',  'am_' // TRIM( salsa_char ), ierr )
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                WRITE( salsa_char,'(i0)' ) lg
                CALL pmc_set_dataarray_name( 'parent', 'sg_' // TRIM( salsa_char ),                &
                                             'child',  'sg_' // TRIM( salsa_char ), ierr )
             ENDDO
          ENDIF
       ENDIF

       CALL pmc_set_dataarray_name( lastentry = .TRUE. )

!
!--    Set switches to control interpolation at bottom and top boundaries.
       IF ( .NOT. ocean_mode )  THEN
!
!--       Bottom condition for elevated childs. Then the bottom boundary is "open".
          IF ( nest_shift_z /= 0.0_wp )  THEN
             nested_bc_at_bottom = .TRUE.
          ENDIF
!
!--       Child top is always "open".
!> TODO:  Should be adjusted for the case if the top of the child matches the top of the root.
          nested_bc_at_top = .TRUE.
       ELSE
!
!--       By default the ocean child is elevated, since its top boundary matches the sea surface.
!--       Here nest_shift_z is interpreted as the distance of the top child boundary from the sea
!--       surface.
!--       ATTENTION: The current code does not allow nest_shift_z /= 0.0 for the ocean mode.
          IF ( nest_shift_z /= 0.0_wp )  THEN
             nested_bc_at_top = .TRUE.
          ENDIF
!
!--       Child bottom is always "open".
          nested_bc_at_bottom = .TRUE.

       ENDIF

!
!--    Send grid to parent
       child_grid_dim(1)  = nx
       child_grid_dim(2)  = ny
       child_grid_dim(3)  = nz
       IF ( .NOT. ocean_mode )  THEN
          child_grid_info(1) = zw(nzt+1)
          child_grid_info(2) = zw(nzt)
       ELSE
          child_grid_info(1) = zw(nzb)
          child_grid_info(2) = zw(nzb)
       ENDIF
       child_grid_info(3) = dx
       child_grid_info(4) = dy
       child_grid_info(5) = dz(1)
       child_grid_info(6) = nest_shift_z

       IF ( myid == 0 )  THEN

          CALL pmc_send_to_parent( child_grid_dim, SIZE( child_grid_dim ), 0, 123, ierr )
          CALL pmc_send_to_parent( child_grid_info, SIZE( child_grid_info ), 0, 124, ierr )
          CALL pmc_send_to_parent( coord_x, nx + 1 + 2 * nbgp, 0, 11, ierr )
          CALL pmc_send_to_parent( coord_y, ny + 1 + 2 * nbgp, 0, 12, ierr )

          CALL pmc_recv_from_parent( rans_mode_parent, 1, 0, 19, ierr )
!
!--       Receive parent-grid information.
          CALL pmc_recv_from_parent( parent_grid_info_real, SIZE( parent_grid_info_real ), 0, 21,  &
                                     ierr )
          CALL pmc_recv_from_parent( parent_grid_info_int, SIZE( parent_grid_info_int ), 0, 22,    &
                                     ierr )

       ENDIF

       CALL MPI_BCAST( parent_grid_info_real, SIZE( parent_grid_info_real ), MPI_REAL, 0, comm2d,  &
                       ierr )
       CALL MPI_BCAST( parent_grid_info_int, SIZE( parent_grid_info_int ), MPI_INTEGER, 0, comm2d, &
                       ierr )

       pg%dx = parent_grid_info_real(3)
       pg%dy = parent_grid_info_real(4)
       pg%dz = parent_grid_info_real(7)
       pg%nx = parent_grid_info_int(1)
       pg%ny = parent_grid_info_int(2)
       pg%nz = parent_grid_info_int(3)
!
!--    Allocate 1-D arrays for parent-grid coordinates and grid-spacings in the z-direction
       ALLOCATE( pg%coord_x(-nbgp:pg%nx+nbgp) )
       ALLOCATE( pg%coord_y(-nbgp:pg%ny+nbgp) )
       ALLOCATE( pg%dzu(1:pg%nz+1) )
       ALLOCATE( pg%dzw(1:pg%nz+1) )
       ALLOCATE( pg%zu(0:pg%nz+1) )
       ALLOCATE( pg%zw(0:pg%nz+1) )
!
!--    Get parent-grid coordinates and grid-spacings in the z-direction from the parent
       IF ( myid == 0)  THEN
          CALL pmc_recv_from_parent( pg%coord_x, pg%nx+1+2*nbgp, 0, 24, ierr )
          CALL pmc_recv_from_parent( pg%coord_y, pg%ny+1+2*nbgp, 0, 25, ierr )
          CALL pmc_recv_from_parent( pg%dzu, pg%nz+1, 0, 26, ierr )
          CALL pmc_recv_from_parent( pg%dzw, pg%nz+1, 0, 27, ierr )
          CALL pmc_recv_from_parent( pg%zu, pg%nz+2, 0, 28, ierr )
          CALL pmc_recv_from_parent( pg%zw, pg%nz+2, 0, 29, ierr )
          CALL pmc_recv_from_parent( humidity_remote, 1, 0, 30, ierr )
       ENDIF
!
!--    Broadcast this information
       CALL MPI_BCAST( pg%coord_x, pg%nx+1+2*nbgp, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%coord_y, pg%ny+1+2*nbgp, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%dzu, pg%nz+1, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%dzw, pg%nz+1, MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%zu, pg%nz+2,  MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( pg%zw, pg%nz+2,  MPI_REAL, 0, comm2d, ierr )
       CALL MPI_BCAST( rans_mode_parent, 1, MPI_LOGICAL, 0, comm2d, ierr )
       CALL MPI_BCAST( humidity_remote, 1, MPI_LOGICAL, 0, comm2d, ierr )
!
!--    Find the index bounds for the nest domain in the parent-grid index space
       CALL pmci_map_child_grid_to_parent_grid
!
!--    TO_DO: Klaus give a comment what is happening here
       CALL pmc_c_get_2d_index_list
!
!--    Include couple arrays into child content
!--    TO_DO: Klaus: better explain the above comment (what is child content?)
       CALL  pmc_c_clear_next_array_list

       n  = 1
       lb = 1
       lc = 1
       lg = 1

       DO  WHILE ( pmc_c_getnextarray( myname ) )
!
!--       Note that pg%nz is not the original nz of parent, but the highest parent-grid level needed
!--       for nesting. Note that in case of chemical species or SALSA variables an additional
!--       parameter needs to be passed. The parameter is required to set the pointer correctly to
!--       the chemical-species or SALSA data structure. Hence, first check if the current variable
!--       is a chemical species or a SALSA variable. If so, pass index id of respective sub-variable
!--       (species or bin) and increment this subsequently.
          IF ( INDEX( TRIM( myname ), 'chem_' ) /= 0 )  THEN
             CALL pmci_create_childs_parent_grid_arrays( myname, ipl, ipr, jps, jpn, pg%nz, n )
             n = n + 1
          ELSEIF ( INDEX( TRIM( myname ), 'an_' ) /= 0 )  THEN
             CALL pmci_create_childs_parent_grid_arrays( myname, ipl, ipr, jps, jpn, pg%nz, lb )
             lb = lb + 1
          ELSEIF ( INDEX( TRIM( myname ), 'am_' ) /= 0 )  THEN
             CALL pmci_create_childs_parent_grid_arrays( myname, ipl, ipr, jps, jpn, pg%nz, lc )
             lc = lc + 1
          ELSEIF ( INDEX( TRIM( myname ), 'sg_' ) /= 0  .AND.  .NOT.  salsa_gases_from_chem )  THEN
             CALL pmci_create_childs_parent_grid_arrays( myname, ipl, ipr, jps, jpn, pg%nz, lg )
             lg = lg + 1
          ELSE
             CALL pmci_create_childs_parent_grid_arrays( myname, ipl, ipr, jps, jpn, pg%nz )
          ENDIF
       ENDDO
       CALL pmc_c_setind_and_allocmem
!
!--    Precompute the index-mapping arrays
       CALL pmci_define_index_mapping
!
!--    Check that the child and parent grid lines do match.
!--    Atmosphere-ocean coupled runs use same size model domains and don't require checking.
       IF ( .NOT. atmosphere_ocean_coupled_run )  THEN
          CALL pmci_check_grid_matching
!
!--       Compute surface areas of the nest-boundary faces
          CALL pmci_compute_face_areas
!
!--       Compute anterpolation lower kp-index bounds if necessary
          CALL pmci_compute_kpb_anterp
       ENDIF
    ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Determine index bounds of interpolation/anterpolation area in the parent-grid index space.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_map_child_grid_to_parent_grid

    IMPLICIT NONE

    INTEGER(iwp) ::  ip     !< Running parent-grid index in the x-direction
    INTEGER(iwp) ::  iauxl  !< Offset between the index bound ipl and the auxiliary index bound ipla
    INTEGER(iwp) ::  iauxr  !< Offset between the index bound ipr and the auxiliary index bound ipra
    INTEGER(iwp) ::  ijaux  !< Temporary variable for receiving the index bound from the neighbouring subdomain
    INTEGER(iwp) ::  jp     !< Running parent-grid index in the y-direction
    INTEGER(iwp) ::  jauxs  !< Offset between the index bound jps and the auxiliary index bound jpsa
    INTEGER(iwp) ::  jauxn  !< Offset between the index bound jpn and the auxiliary index bound jpna

    INTEGER(iwp), DIMENSION(4) ::  parent_bound_global  !< Transfer array for global parent-grid index bounds
    INTEGER(iwp), DIMENSION(2) ::  size_of_array        !< For sending the dimensions of parent_bound_all to parent

    INTEGER(iwp), DIMENSION(5,numprocs) ::  parent_bound_all  !< Transfer array for parent-grid index bounds

    REAL(wp) ::  tolex  !< Tolerance for grid-line matching in x-direction
    REAL(wp) ::  toley  !< Tolerance for grid-line matching in y-direction
    REAL(wp) ::  xexl   !< Parent-grid array exceedance behind the left edge of the child PE subdomain
    REAL(wp) ::  xexr   !< Parent-grid array exceedance behind the right edge of the child PE subdomain
    REAL(wp) ::  xpl    !< Requested left-edge x-coordinate of the parent-grid array domain (at the internal boundaries
                        !< the real edge may differ from this in some cases as explained in the comment block below)
    REAL(wp) ::  xpr    !< Requested right-edge x-coordinate of the parent-grid array domain (at the internal boundaries
                        !< the real edge may differ from this in some cases as explained in the comment block below)
    REAL(wp) ::  yexs   !< Parent-grid array exceedance behind the south edge of the child PE subdomain
    REAL(wp) ::  yexn   !< Parent-grid array exceedance behind the north edge of the child PE subdomain
    REAL(wp) ::  yps    !< Requested south-edge y-coordinate of the parent-grid array domain (at the internal boundaries
                        !< the real edge may differ from this in some cases as explained in the comment block below)
    REAL(wp) ::  ypn    !< Requested south-edge y-coordinate of the parent-grid array domain (at the internal boundaries
                        !< the real edge may differ from this in some cases as explained in the comment block below)

!
!-- Determine the index limits for the child's parent-grid arrays (such as uc for example).
!-- Note that at the outer edges of the child domain (nest boundaries) these arrays exceed the
!-- boundary by two parent-grid cells. At the internal boundaries, there are no exceedances and
!-- thus no overlaps with the neighbouring subdomain. If at least half of the parent-grid cell is
!-- within the current child sub-domain, then it is included in the current sub-domain's
!-- parent-grid array. Else the parent-grid cell is included in the neighbouring subdomain's
!-- parent-grid array, or not included at all if we are at the outer edge of the child domain.
!-- This may occur especially when a large grid-spacing ratio is used.
!
!-- Tolerances for grid-line matching.
    tolex = tolefac * dx
    toley = tolefac * dy
!
!-- Left boundary.
!-- Extension by two parent-grid cells behind the boundary, see the comment block above.
    IF ( bc_dirichlet_l  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       xexl  = 2.0_wp * pg%dx
       iauxl = 0
    ELSE
       xexl  = 0.0_wp
       iauxl = 1
    ENDIF
    xpl     = coord_x(nxl) - xexl
    DO  ip = 0, pg%nx
       IF ( pg%coord_x(ip) + 0.5_wp * pg%dx >= xpl - tolex )  THEN
          ipl = MAX( 0, ip )
          EXIT
       ENDIF
    ENDDO
!
!-- Right boundary.
!-- Extension by two parent-grid cells behind the boundary, see the comment block above.
    IF ( bc_dirichlet_r  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       xexr  = 2.0_wp * pg%dx
       iauxr = 0
    ELSE
       xexr  = 0.0_wp
       iauxr = 1
    ENDIF
    xpr  = coord_x(nxr+1) + xexr
    DO  ip = pg%nx, 0 , -1
       IF ( pg%coord_x(ip) + 0.5_wp * pg%dx <= xpr + tolex )  THEN
          ipr = MIN( pg%nx, MAX( ipl, ip ) )
          EXIT
       ENDIF
    ENDDO
!
!-- South boundary.
!-- Extension by two parent-grid cells behind the boundary, see the comment block above.
    IF ( bc_dirichlet_s  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       yexs  = 2.0_wp * pg%dy
       jauxs = 0
    ELSE
       yexs  = 0.0_wp
       jauxs = 1
    ENDIF
    yps  = coord_y(nys) - yexs
    DO  jp = 0, pg%ny
       IF ( pg%coord_y(jp) + 0.5_wp * pg%dy >= yps - toley )  THEN
          jps = MAX( 0, jp )
          EXIT
       ENDIF
    ENDDO
!
!-- North boundary.
!-- Extension by two parent-grid cells behind the boundary, see the comment block above.
    IF  ( bc_dirichlet_n  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       yexn  = 2.0_wp * pg%dy
       jauxn = 0
    ELSE
       yexn  = 0.0_wp
       jauxn = 1
    ENDIF
    ypn  = coord_y(nyn+1) + yexn
    DO  jp = pg%ny, 0 , -1
       IF ( pg%coord_y(jp) + 0.5_wp * pg%dy <= ypn + toley )  THEN
          jpn = MIN( pg%ny, MAX( jps, jp ) )
          EXIT
       ENDIF
    ENDDO

!
!-- For ocean coupling, one more parent cell is required at right and north boundary
!-- to enable bilinear interpolation.
    IF( atmosphere_ocean_coupled_run )  THEN

       ipr = ipr + 1
       jpn = jpn + 1

    ELSE

!
!--    Make sure that the indexing is contiguous (no gaps, no overlaps). This is a safety measure
!--    mainly for cases with high grid-spacing ratio and narrow child subdomains.
       IF ( npex > 1 )  THEN
          IF ( nxl == 0 )  THEN
             CALL MPI_SEND( ipr, 1, MPI_INTEGER, pright, 717, comm2d, ierr )
          ELSEIF ( nxr == nx )  THEN
             CALL MPI_RECV( ijaux, 1, MPI_INTEGER, pleft, 717, comm2d, status, ierr )
             ipl = ijaux + 1
          ELSE
             CALL MPI_SEND( ipr, 1, MPI_INTEGER, pright, 717, comm2d, ierr )
             CALL MPI_RECV( ijaux, 1, MPI_INTEGER, pleft, 717, comm2d, status, ierr )
             ipl = ijaux + 1
          ENDIF
       ENDIF

       IF ( npey > 1 )  THEN
          IF ( nys == 0 )  THEN
             CALL MPI_SEND( jpn, 1, MPI_INTEGER, pnorth, 719, comm2d, ierr )
          ELSEIF ( nyn == ny )  THEN
             CALL MPI_RECV( ijaux, 1, MPI_INTEGER, psouth, 719, comm2d, status, ierr )
             jps = ijaux + 1
          ELSE
             CALL MPI_SEND( jpn, 1, MPI_INTEGER, pnorth, 719, comm2d, ierr )
             CALL MPI_RECV( ijaux, 1, MPI_INTEGER, psouth, 719, comm2d, status, ierr )
             jps = ijaux + 1
          ENDIF
       ENDIF
    ENDIF

    IF ( debug_output )  THEN
       WRITE( debug_string,                                                                        &
              '("pmci_map_child_grid_to_parent_grid. Parent-grid array bounds: ",4(I4,2X))' )      &
            ipl, ipr, jps, jpn
       CALL debug_message( debug_string, 'info' )
    ENDIF
    parent_bound(1) = ipl
    parent_bound(2) = ipr
    parent_bound(3) = jps
    parent_bound(4) = jpn
    parent_bound(5) = myid
!
!-- The following auxiliary index bounds are used for allocating index mapping and some other
!-- auxiliary arrays.
    ipla = ipl - iauxl
    ipra = ipr + iauxr
    jpsa = jps - jauxs
    jpna = jpn + jauxn
!
!-- The index-bounds parent_bound of all subdomains of the current child domain must be sent to the
!-- parent in order for the parent to create the index list. For this reason, the parent_bound
!-- arrays are packed together in single array parent_bound_all using MPI_GATHER. Note that
!-- MPI_Gather receives data from all processes in the rank order This fact is exploited in creating
!-- the index list in pmci_create_index_list.
    CALL MPI_GATHER( parent_bound, 5, MPI_INTEGER, parent_bound_all, 5, MPI_INTEGER, 0, comm2d,    &
                     ierr )

    IF ( myid == 0 )  THEN
       size_of_array(1) = SIZE( parent_bound_all, 1 )
       size_of_array(2) = SIZE( parent_bound_all, 2 )
       CALL pmc_send_to_parent( size_of_array, 2, 0, 40, ierr )
       CALL pmc_send_to_parent( parent_bound_all, SIZE( parent_bound_all ), 0, 41, ierr )
!
!--    Determine the global parent-grid index bounds
       parent_bound_global(1) = MINVAL( parent_bound_all(1,:) )
       parent_bound_global(2) = MAXVAL( parent_bound_all(2,:) )
       parent_bound_global(3) = MINVAL( parent_bound_all(3,:) )
       parent_bound_global(4) = MAXVAL( parent_bound_all(4,:) )
    ENDIF
!
!-- Broadcast the global parent-grid index bounds to all current child processes
    CALL MPI_BCAST( parent_bound_global, 4, MPI_INTEGER, 0, comm2d, ierr )
    iplg = parent_bound_global(1)
    iprg = parent_bound_global(2)
    jpsg = parent_bound_global(3)
    jpng = parent_bound_global(4)

    IF ( debug_output )  THEN
       WRITE( debug_string,                                                                        &
              '("pmci_map_child_grid_to_parent_grid: global parent-grid index bounds ",            &
              & "iplg, iprg, jpsg, jpng: ",4(I4,2X))' )  iplg, iprg, jpsg, jpng
       CALL debug_message( debug_string, 'info' )
    ENDIF

 END SUBROUTINE pmci_map_child_grid_to_parent_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Precomputation of the mapping between the child- and parent-grid indices.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_define_index_mapping

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< child-grid index in the x-direction
    INTEGER(iwp) ::  ii      !< parent-grid index in the x-direction
    INTEGER(iwp) ::  istart  !< starting index for the index-mapping search loop in the x-direction
    INTEGER(iwp) ::  ir      !< search-loop running index in the x-direction
    INTEGER(iwp) ::  iw      !< child-grid index limited to -1 <= iw <= nx+1 for topo_flags
    INTEGER(iwp) ::  j       !< child-grid index in the y-direction
    INTEGER(iwp) ::  jj      !< parent-grid index in the y-direction
    INTEGER(iwp) ::  jstart  !< starting index for the index-mapping search loop in the y-direction
    INTEGER(iwp) ::  jr      !< search-loop running index in the y-direction
    INTEGER(iwp) ::  jw      !< child-grid index limited to -1 <= jw <= ny+1 for topo_flags
    INTEGER(iwp) ::  k       !< child-grid index in the z-direction
    INTEGER(iwp) ::  kk      !< parent-grid index in the z-direction
    INTEGER(iwp) ::  kstart  !< starting index for the index-mapping search loop in the z-direction
    INTEGER(iwp) ::  kw      !< child-grid index limited to kw <= nzt+1 for topo_flags

    REAL(wp) ::  tolex  !< tolerance for grid-line matching in x-direction
    REAL(wp) ::  toley  !< tolerance for grid-line matching in y-direction
    REAL(wp) ::  tolez  !< tolerance for grid-line matching in z-direction

!
!-- Grid-line tolerances.
    tolex = tolefac * dx
    toley = tolefac * dy
    tolez = tolefac * dz(1)
!
!-- Allocate child-grid work arrays for interpolation.
    igsr = NINT( pg%dx / dx, iwp )
    jgsr = NINT( pg%dy / dy, iwp )
    kgsr = NINT( pg%dzw(1) / dzw(1), iwp )
!
!-- Determine index bounds for the parent-grid work arrays for interpolation and allocate them.
    CALL pmci_allocate_workarrays
!
!-- Define the MPI-datatypes for parent-grid work array exchange between the PE-subdomains.
    CALL pmci_create_workarray_exchange_datatypes
!
!-- The setups further below are for nesting only, therefore return here in case of coupling.
    IF ( atmosphere_ocean_coupled_run )  RETURN

!
!-- First determine kcto and kctw which refer to the uppermost parent-grid levels below the child
!-- top-boundary level. Note that these comparison tests are not round-off-error sensitive and
!-- therefore tolerance buffering is not needed here.
    kk = 0
    DO WHILE ( pg%zu(kk) <= zu(nzt) )
       kk = kk + 1
    ENDDO
    kcto = kk - 1

    kk = 0
    DO WHILE ( pg%zw(kk) <= zw(nzt-1) )
       kk = kk + 1
    ENDDO
    kctw = kk - 1

    IF ( debug_output )  THEN
       WRITE( debug_string, '("kcto, kctw = ", 2(I3,2X))' )  kcto, kctw
       CALL debug_message( debug_string, 'info' )
    ENDIF
!
!-- In case of two-way coupling, check that the child domain is sufficiently large in terms of the
!-- number of parent-grid cells covered. Otherwise anterpolation is not possible.
    IF ( nesting_mode == 'two-way' )  THEN
       CALL pmci_check_child_domain_size
    ENDIF

    ALLOCATE( iflu(ipla:ipra) )
    ALLOCATE( iflo(ipla:ipra) )
    ALLOCATE( ifuu(ipla:ipra) )
    ALLOCATE( ifuo(ipla:ipra) )
    ALLOCATE( jflv(jpsa:jpna) )
    ALLOCATE( jflo(jpsa:jpna) )
    ALLOCATE( jfuv(jpsa:jpna) )
    ALLOCATE( jfuo(jpsa:jpna) )
    ALLOCATE( kflw(0:pg%nz+1) )
    ALLOCATE( kflo(0:pg%nz+1) )
    ALLOCATE( kfuw(0:pg%nz+1) )
    ALLOCATE( kfuo(0:pg%nz+1) )
    ALLOCATE( ijkfc_u(0:pg%nz+1,jpsa:jpna,ipla:ipra) )
    ALLOCATE( ijkfc_v(0:pg%nz+1,jpsa:jpna,ipla:ipra) )
    ALLOCATE( ijkfc_w(0:pg%nz+1,jpsa:jpna,ipla:ipra) )
    ALLOCATE( ijkfc_s(0:pg%nz+1,jpsa:jpna,ipla:ipra) )

    ijkfc_u = 0
    ijkfc_v = 0
    ijkfc_w = 0
    ijkfc_s = 0
!
!-- i-indices of u for each ii-index value
    istart = nxlg
    DO  ii = ipla, ipra
!
!--    The parent and child grid lines do always match in x, hence we use only the local
!--    k,j-child-grid plane for the anterpolation. However, iflu still has to be stored separately
!--    as these index bounds are passed as arguments to the interpolation and anterpolation
!--    subroutines. Note that this comparison test is round-off-error sensitive and therefore
!--    tolerance buffering is needed here.
       i = istart
       DO WHILE ( pg%coord_x(ii) - coord_x(i) > tolex  .AND. i < nxrg )
          i = i + 1
       ENDDO
       iflu(ii) = MIN( MAX( i, nxlg ), nxrg )
       ifuu(ii) = iflu(ii)
       istart   = iflu(ii)
!
!--    Print out the index bounds for checking and debugging purposes
       IF ( debug_output )  THEN
          WRITE( debug_string, '("pmci_define_index_mapping, ii, iflu, ifuu: ", 3(I4,2X))' )       &
                 ii, iflu(ii), ifuu(ii)
          CALL debug_message( debug_string, 'info' )
       ENDIF
    ENDDO
!
!-- i-indices of others for each ii-index value. Note that these comparison tests are not
!-- round-off-error sensitive and therefore tolerance buffering is not needed here.
    istart = nxlg
    DO  ii = ipla, ipra
       i = istart
       DO WHILE ( ( coord_x(i) + 0.5_wp * dx < pg%coord_x(ii) )  .AND.  ( i < nxrg ) )
          i  = i + 1
       ENDDO
       iflo(ii) = MIN( MAX( i, nxlg ), nxrg )
       ir = i
       DO WHILE ( ( coord_x(ir) + 0.5_wp * dx < pg%coord_x(ii) + pg%dx )  .AND.  ( i < nxrg+1 ) )
          i  = i + 1
          ir = MIN( i, nxrg )
       ENDDO
       ifuo(ii) = MIN( MAX( i-1, iflo(ii) ), nxrg )
       istart = iflo(ii)
!
!--    Print out the index bounds for checking and debugging purposes
       IF ( debug_output )  THEN
          WRITE( debug_string, '("pmci_define_index_mapping, ii, iflo, ifuo: ", 3(I4,2X))' )       &
                 ii, iflo(ii), ifuo(ii)
          CALL debug_message( debug_string, 'info' )
       ENDIF
    ENDDO
!
!-- j-indices of v for each jj-index value
    jstart = nysg
    DO  jj = jpsa, jpna
!
!--    The parent and child grid lines do always match in y, hence we use only the local
!--    k,i-child-grid plane for the anterpolation. However, jcnv still has to be stored separately
!--    as these index bounds are passed as arguments to the interpolation and anterpolation
!--    subroutines. Note that this comparison test is round-off-error sensitive and therefore
!--    tolerance buffering is needed here.
       j = jstart
       DO WHILE ( pg%coord_y(jj) - coord_y(j) > toley  .AND.  j < nyng )
          j = j + 1
       ENDDO
       jflv(jj) = MIN( MAX( j, nysg ), nyng )
       jfuv(jj) = jflv(jj)
       jstart   = jflv(jj)
!
!--    Print out the index bounds for checking and debugging purposes
       IF ( debug_output )  THEN
          WRITE( debug_string, '("pmci_define_index_mapping, jj, jflv, jfuv: ", 3(I4,2X))' )       &
                 jj, jflv(jj), jfuv(jj)
          CALL debug_message( debug_string, 'info' )
       ENDIF
    ENDDO
!
!-- j-indices of others for each jj-index value
!-- Note that these comparison tests are not round-off-error sensitive and therefore tolerance
!-- buffering is not needed here.
    jstart = nysg
    DO  jj = jpsa, jpna
       j = jstart
       DO WHILE ( ( coord_y(j) + 0.5_wp * dy < pg%coord_y(jj) ) .AND. ( j < nyng ) )
          j  = j + 1
       ENDDO
       jflo(jj) = MIN( MAX( j, nysg ), nyng )
       jr = j
       DO WHILE ( ( coord_y(jr) + 0.5_wp * dy < pg%coord_y(jj) + pg%dy ) .AND. ( j < nyng+1 ) )
          j  = j + 1
          jr = MIN( j, nyng )
       ENDDO
       jfuo(jj) = MIN( MAX( j-1, jflo(jj) ), nyng )
       jstart = jflo(jj)
!
!--    Print out the index bounds for checking and debugging purposes
       IF ( debug_output )  THEN
          WRITE( debug_string, '("pmci_define_index_mapping, jj, jflo, jfuo: ", 3(I4,2X))' )       &
                 jj, jflo(jj), jfuo(jj)
          CALL debug_message( debug_string, 'info' )
       ENDIF
    ENDDO
!
!-- k-indices of w for each kk-index value
!-- Note that anterpolation index limits are needed also for the top boundary ghost cell level
!-- because they are used also in the interpolation.
    kstart  = 0
    kflw(0) = 0
    kfuw(0) = 0
    DO  kk = 1, pg%nz+1
!
!--    The parent and child grid lines do always match in z, hence we use only the local
!--    j,i-child-grid plane for the anterpolation. However, kctw still has to be stored separately
!--    as these index bounds are passed as arguments to the interpolation and anterpolation
!--    subroutines. Note that this comparison test is round-off-error sensitive and therefore
!--    tolerance buffering is needed here.
       k = kstart
       DO WHILE ( ( pg%zw(kk) - zw(k) > tolez )  .AND.  ( k < nzt+1 ) )
          k = k + 1
       ENDDO
       kflw(kk) = MIN( MAX( k, 1 ), nzt + 1 )
       kfuw(kk) = kflw(kk)
       kstart   = kflw(kk)
!
!--    Print out the index bounds for checking and debugging purposes
       IF ( debug_output )  THEN
          WRITE( debug_string, '("pmci_define_index_mapping, kk, kflw, kfuw: ", 4(I4,2X),          &
                 & 2(E12.5,2X))' )  kk, kflw(kk), kfuw(kk), nzt,  pg%zu(kk), pg%zw(kk)
          CALL debug_message( debug_string, 'info' )
       ENDIF
    ENDDO
!
!-- k-indices of others for each kk-index value
    kstart  = 0
    kflo(0) = 0
    kfuo(0) = 0
!
!-- Note that anterpolation index limits are needed also for the top boundary ghost cell level
!-- because they are used also in the interpolation. Note that these comparison tests are not
!-- round-off-error sensitive and therefore tolerance buffering is not needed here.
    DO  kk = 1, pg%nz+1
       k = kstart
       DO WHILE ( ( zu(k) < pg%zw(kk-1) )  .AND.  ( k <= nzt ) )
          k = k + 1
       ENDDO
       kflo(kk) = MIN( MAX( k, 1 ), nzt + 1 )
       DO WHILE ( ( zu(k) < pg%zw(kk) )  .AND.  ( k <= nzt+1 ) )
          k = k + 1
          IF ( k > nzt + 1 ) EXIT  ! This EXIT is to prevent zu(k) from flowing over.
       ENDDO
       kfuo(kk) = MIN( MAX( k-1, kflo(kk) ), nzt + 1 )
       kstart = kflo(kk)
    ENDDO
!
!-- Print out the index bounds for checking and debugging purposes
    IF ( debug_output )  THEN
       DO  kk = 1, pg%nz+1
          WRITE( debug_string, '("pmci_define_index_mapping, kk, kflo, kfuo: ", 4(I4,2X),          &
                 & 2(E12.5,2X))' )  kk, kflo(kk), kfuo(kk), nzt,  pg%zu(kk), pg%zw(kk)
          CALL debug_message( debug_string, 'info' )
       ENDDO
    ENDIF
!
!-- Precomputation of number of child-grid nodes inside parent-grid cells. Note that ii, jj, and kk
!-- are parent-grid indices. This information is needed in the anterpolation. The indices for
!-- topo_flags (kw,jw,iw) must be limited to the range [-1,...,nx/ny/nzt+1] in order to
!-- avoid zero values on the outer ghost nodes.
    DO  ii = ipla, ipra
       DO  jj = jpsa, jpna
          DO  kk = 0, pg%nz+1
!
!--          u-component
             DO  i = iflu(ii), ifuu(ii)
                iw = MAX( MIN( i, nx+1 ), -1 )
                DO  j = jflo(jj), jfuo(jj)
                   jw = MAX( MIN( j, ny+1 ), -1 )
                   DO  k = kflo(kk), kfuo(kk)
                      kw = MIN( k, nzt+1 )
                      ijkfc_u(kk,jj,ii) = ijkfc_u(kk,jj,ii)                                        &
                                          + MERGE( 1, 0, BTEST( topo_flags(kw,jw,iw), 1 ) )
                   ENDDO
                ENDDO
             ENDDO
!
!--          v-component
             DO  i = iflo(ii), ifuo(ii)
                iw = MAX( MIN( i, nx+1 ), -1 )
                DO  j = jflv(jj), jfuv(jj)
                   jw = MAX( MIN( j, ny+1 ), -1 )
                   DO  k = kflo(kk), kfuo(kk)
                      kw = MIN( k, nzt+1 )
                      ijkfc_v(kk,jj,ii) = ijkfc_v(kk,jj,ii)                                        &
                                          + MERGE( 1, 0, BTEST( topo_flags(kw,jw,iw), 2 ) )
                   ENDDO
                ENDDO
             ENDDO
!
!--          Scalars
             DO  i = iflo(ii), ifuo(ii)
                iw = MAX( MIN( i, nx+1 ), -1 )
                DO  j = jflo(jj), jfuo(jj)
                   jw = MAX( MIN( j, ny+1 ), -1 )
                   DO  k = kflo(kk), kfuo(kk)
                      kw = MIN( k, nzt+1 )
                      ijkfc_s(kk,jj,ii) = ijkfc_s(kk,jj,ii)                                        &
                                          + MERGE( 1, 0, BTEST( topo_flags(kw,jw,iw), 0 ) )
                   ENDDO
                ENDDO
             ENDDO
!
!--          w-component
             DO  i = iflo(ii), ifuo(ii)
                iw = MAX( MIN( i, nx+1 ), -1 )
                DO  j = jflo(jj), jfuo(jj)
                   jw = MAX( MIN( j, ny+1 ), -1 )
                   DO  k = kflw(kk), kfuw(kk)
                      kw = MIN( k, nzt+1 )
                      ijkfc_w(kk,jj,ii) = ijkfc_w(kk,jj,ii)                                        &
                                          + MERGE( 1, 0, BTEST( topo_flags(kw,jw,iw), 3 ) )
                   ENDDO
                ENDDO
             ENDDO

          ENDDO  ! kk
       ENDDO  ! jj
    ENDDO  ! ii

 END SUBROUTINE pmci_define_index_mapping



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check if the child domain is too small in terms of number of parent-grid cells covered so that
!> anterpolation buffers fill the whole domain so that anterpolation is not possible. Also, check
!> that anterpolation_buffer_width is not too large to prevent anterpolation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_check_child_domain_size

    IMPLICIT NONE

!
!-- First x-direction
    IF ( iplg + 3 + anterpolation_buffer_width > iprg - 3 - anterpolation_buffer_width )  THEN
       IF ( iprg - iplg + 1 < 7 )  THEN
!
!--       Error
          message_string = 'child domain too narrow for anterpolation in x-direction'
          CALL message( 'pmci_check_child_domain_size', 'PMC0014', 3, 2, 0, 6, 0 )
       ELSEIF ( iprg - iplg + 1 < 11 )  THEN
!
!--       Warning
          message_string = 'anterpolation_buffer_width value too high, reset to 0'
          CALL message( 'pmci_check_child_domain_size', 'PMC0015', 0, 1, 0, 6, 0 )
          anterpolation_buffer_width = 0
       ELSE
!
!--       Informative message
          message_string = 'anterpolation_buffer_width value too high, reset to default value 2'
          CALL message( 'pmci_check_child_domain_size', 'PMC0016', 0, 0, 0, 6, 0 )
          anterpolation_buffer_width = 2
       ENDIF
    ENDIF
!
!-- Then y-direction
    IF ( jpsg + 3 + anterpolation_buffer_width > jpng - 3 - anterpolation_buffer_width )  THEN
       IF ( jpng - jpsg + 1 < 7 )  THEN
!
!--       Error
          message_string = 'child domain too narrow for anterpolation in y-direction'
          CALL message( 'pmci_check_child_domain_size', 'PMC0014', 3, 2, 0, 6, 0 )
       ELSEIF ( jpng - jpsg + 1 < 11 )  THEN
!
!--       Warning
          message_string = 'anterpolation_buffer_width value too high, reset to 0'
          CALL message( 'pmci_check_child_domain_size', 'PMC0015', 0, 1, 0, 6, 0 )
          anterpolation_buffer_width = 0
       ELSE
!
!--       Informative message
          message_string = 'anterpolation_buffer_width value too high, reset to default value 2'
          CALL message( 'pmci_check_child_domain_size', 'PMC0016', 0, 0, 0, 6, 0 )
          anterpolation_buffer_width = 2
       ENDIF
    ENDIF
!
!-- Finally z-direction
    IF ( nested_bc_at_bottom )  THEN
       IF ( 1 + anterpolation_buffer_width > kctw - 1 - anterpolation_buffer_width )  THEN
          IF ( kctw - 1 < 1 )  THEN
!
!--          Error
             message_string = 'child domain too shallow for anterpolation in z-direction'
             CALL message( 'pmci_check_child_domain_size', 'PMC0014', 3, 2, 0, 6, 0 )
          ELSEIF ( kctw - 3 < 1 )  THEN
!
!--          Warning
             message_string = 'anterpolation_buffer_width value too high, reset to 0'
             CALL message( 'pmci_check_child_domain_size', 'PMC0015', 0, 1, 0, 6, 0 )
             anterpolation_buffer_width = 0
          ELSE
!
!--          Informative message
             message_string = 'anterpolation_buffer_width value too high, reset to default value 2'
             CALL message( 'pmci_check_child_domain_size', 'PMC0016', 0, 0, 0, 6, 0 )
             anterpolation_buffer_width = 2
          ENDIF
       ENDIF
    ELSE
       IF ( kctw - 1 - anterpolation_buffer_width < 1 )  THEN
          IF ( kctw - 1 < 1 )  THEN
!
!--          Error
             message_string = 'child domain too shallow for anterpolation in z-direction'
             CALL message( 'pmci_check_child_domain_size', 'PMC0014', 3, 2, 0, 6, 0 )
          ELSEIF ( kctw - 3 < 1 )  THEN
!
!--          Warning
             message_string = 'anterpolation_buffer_width value too high, reset to 0'
             CALL message( 'pmci_check_child_domain_size', 'PMC0015', 0, 1, 0, 6, 0 )
             anterpolation_buffer_width = 0
          ELSE
!
!--          Informative message
             message_string = 'anterpolation_buffer_width value too high, reset to default value 2'
             CALL message( 'pmci_check_child_domain_size', 'PMC0016', 0, 0, 0, 6, 0 )
             anterpolation_buffer_width = 2
          ENDIF
       ENDIF
    ENDIF

 END SUBROUTINE pmci_check_child_domain_size


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate parent-grid work-arrays for interpolation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_allocate_workarrays

    IMPLICIT NONE

!
!-- Determine and store the PE-subdomain dependent index bounds
    IF ( bc_dirichlet_l  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       iplw = ipl + 1
    ELSE
       iplw = ipl - 1
    ENDIF

    IF ( bc_dirichlet_r  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       iprw = ipr - 1
    ELSE
       iprw = ipr + 1
    ENDIF

    IF ( bc_dirichlet_s  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       jpsw = jps + 1
    ELSE
       jpsw = jps - 1
    ENDIF

    IF ( bc_dirichlet_n  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       jpnw = jpn - 1
    ELSE
       jpnw = jpn + 1
    ENDIF
!
!-- Left and right boundaries.
    ALLOCATE( workarr_lr(0:pg%nz+1,jpsw:jpnw,0:2) )
!
!-- South and north boundaries.
    ALLOCATE( workarr_sn(0:pg%nz+1,0:2,iplw:iprw) )
!
!-- Bottom and top boundary.
    ALLOCATE( workarr_bt(0:2,jpsw:jpnw,iplw:iprw) )

 END SUBROUTINE pmci_allocate_workarrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define and create specific MPI-datatypes for the interpolation work-array exchange.
!> Exhcanges are needed to make these arrays contiguous over the horizontal direction on the
!> plane of the boundary.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_create_workarray_exchange_datatypes

    IMPLICIT NONE

!
!-- For the left and right boundaries
    CALL MPI_TYPE_VECTOR( 3, pg%nz+2, (jpnw-jpsw+1)*(pg%nz+2), MPI_REAL, workarr_lr_exchange_type, &
                          ierr )
    CALL MPI_TYPE_COMMIT( workarr_lr_exchange_type, ierr )
!
!-- For the south and north boundaries
    CALL MPI_TYPE_VECTOR( 1, 3*(pg%nz+2), 3*(pg%nz+2), MPI_REAL, workarr_sn_exchange_type, ierr )
    CALL MPI_TYPE_COMMIT( workarr_sn_exchange_type, ierr )
!
!-- For the bottom and top boundary x-slices
    CALL MPI_TYPE_VECTOR( iprw-iplw+1, 3, 3*(jpnw-jpsw+1), MPI_REAL, workarr_bt_exchange_type_x,   &
                          ierr )
    CALL MPI_TYPE_COMMIT( workarr_bt_exchange_type_x, ierr )
!
!-- For the bottom and top boundary y-slices
    CALL MPI_TYPE_VECTOR( 1, 3*(jpnw-jpsw+1), 3*(jpnw-jpsw+1), MPI_REAL,                           &
                          workarr_bt_exchange_type_y, ierr )
    CALL MPI_TYPE_COMMIT( workarr_bt_exchange_type_y, ierr )

 END SUBROUTINE pmci_create_workarray_exchange_datatypes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check that the grid lines of child and parent do match. Also check that the child subdomain
!> width is not smaller than the parent grid spacing in the respective direction.
!> In case of elevated childs, their bottom boundaries must not lie below the bottom boundary of
!> their respective parent.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_check_grid_matching

    IMPLICIT NONE

    INTEGER(iwp) ::  non_int_gsr_x = 0                    !< Flag for non-integer grid-spacing ration in x-direction
    INTEGER(iwp) ::  non_int_gsr_y = 0                    !< Flag for non-integer grid-spacing ration in y-direction
    INTEGER(iwp) ::  non_int_gsr_z = 0                    !< Flag for non-integer grid-spacing ration in z-direction
    INTEGER(iwp) ::  non_matching_height = 0              !< Flag for non-matching child-domain height
    INTEGER(iwp) ::  non_matching_lower_left_corner = 0   !< Flag for non-matching lower left corner
    INTEGER(iwp) ::  non_matching_upper_right_corner = 0  !< Flag for non-matching upper right corner
    INTEGER(iwp) ::  too_narrow_pesd_x = 0                !< Flag for too narrow pe-subdomain in x-direction
    INTEGER(iwp) ::  too_narrow_pesd_y = 0                !< Flag for too narrow pe-subdomain in y-direction

    REAL(wp) ::  child_ngp_x_l                            !< Number of gridpoints in child subdomain in x-direction
                                                          !< converted to REAL(wp)
    REAL(wp) ::  child_ngp_y_l                            !< Number of gridpoints in child subdomain in y-direction
                                                          !< converted to REAL(wp)
    REAL(wp) ::  gridline_mismatch_x                      !< Mismatch between the parent and child gridlines in the x-direction
    REAL(wp) ::  gridline_mismatch_y                      !< Mismatch between the parent and child gridlines in the y-direction
    REAL(wp) ::  gridline_mismatch_z                      !< Mismatch between the parent and child gridlines in the z-direction
    REAL(wp) ::  gsr_mismatch_x                           !< Deviation of the grid-spacing ratio from the nearest integer value,
                                                          !< the x-direction
    REAL(wp) ::  gsr_mismatch_y                           !< Deviation of the grid-spacing ratio from the nearest integer value, the
                                                          !< y-direction
    REAL(wp) ::  tolex                                    !< Tolerance for grid-line matching in x-direction
    REAL(wp) ::  toley                                    !< Tolerance for grid-line matching in y-direction
    REAL(wp) ::  tolez                                    !< Tolerance for grid-line matching in z-direction
    REAL(wp) ::  upper_right_coord_x                      !< X-coordinate of the upper right corner of the child domain
    REAL(wp) ::  upper_right_coord_y                      !< Y-coordinate of the upper right corner of the child domain


    IF ( myid == 0 )  THEN

       tolex = tolefac * dx
       toley = tolefac * dy
       tolez = tolefac * dz(1)
!
!--    First check that the child domain lower left corner matches the parent grid lines.
       gridline_mismatch_x = ABS( NINT( lower_left_coord_x / pg%dx ) * pg%dx - lower_left_coord_x )
       gridline_mismatch_y = ABS( NINT( lower_left_coord_y / pg%dy ) * pg%dy - lower_left_coord_y )
       IF ( gridline_mismatch_x > tolex ) non_matching_lower_left_corner = 1
       IF ( gridline_mismatch_y > toley ) non_matching_lower_left_corner = 1
!
!--    Then check that the child doman upper right corner matches the parent grid lines.
       upper_right_coord_x = lower_left_coord_x + ( nx + 1 ) * dx
       upper_right_coord_y = lower_left_coord_y + ( ny + 1 ) * dy
       gridline_mismatch_x = ABS( NINT( upper_right_coord_x / pg%dx ) * pg%dx - upper_right_coord_x )
       gridline_mismatch_y = ABS( NINT( upper_right_coord_y / pg%dy ) * pg%dy - upper_right_coord_y )
       IF ( gridline_mismatch_x > tolex )  non_matching_upper_right_corner = 1
       IF ( gridline_mismatch_y > toley )  non_matching_upper_right_corner = 1
!
!--    Also check that the child domain height matches the parent grid lines.
       gridline_mismatch_z = ABS( NINT( zw(nzt) / pg%dz ) * pg%dz - zw(nzt) )
       IF ( gridline_mismatch_z > tolez )  non_matching_height = 1
!
!--    Check that the grid-spacing ratios in each direction are integer valued.
       gsr_mismatch_x = ABS( NINT( pg%dx / dx ) * dx - pg%dx )
       gsr_mismatch_y = ABS( NINT( pg%dy / dy ) * dy - pg%dy )
       IF ( gsr_mismatch_x > tolex )  non_int_gsr_x = 1
       IF ( gsr_mismatch_y > toley )  non_int_gsr_y = 1
!
!--    In the z-direction, all levels need to be checked separately against grid stretching which is
!--    not allowed.
       DO  n = 0, kctw+1
          IF ( ABS( pg%zw(n) - zw(kflw(n)) ) > tolez )  non_int_gsr_z = 1
       ENDDO

       child_ngp_x_l = REAL( nxr - nxl + 1, KIND=wp )
       IF ( child_ngp_x_l / REAL( igsr, KIND=wp ) < 1.0_wp )  too_narrow_pesd_x = 1
       child_ngp_y_l = REAL( nyn - nys + 1, KIND=wp )
       IF ( child_ngp_y_l / REAL( jgsr, KIND=wp ) < 1.0_wp )  too_narrow_pesd_y = 1

       IF ( non_matching_height > 0 )  THEN
          message_string = 'nested child domain height must match its parent grid lines'
          CALL message( 'pmci_check_grid_matching', 'PMC0017', 3, 2, 0, 6, 0 )
       ENDIF

       IF ( non_matching_lower_left_corner > 0 )  THEN
          message_string = 'nested child domain lower left corner must match its parent grid lines'
          CALL message( 'pmci_check_grid_matching', 'PMC0017', 3, 2, 0, 6, 0 )
       ENDIF

       IF ( non_matching_upper_right_corner > 0 )  THEN
          message_string = 'nested child domain upper right corner must match its parent grid lines'
          CALL message( 'pmci_check_grid_matching', 'PMC0017', 3, 2, 0, 6, 0 )
       ENDIF

       IF ( non_int_gsr_x > 0 )  THEN
          message_string = 'nesting grid-spacing ratio ( parent dx / child dx ) ' //               &
                           'must have an integer value'
          CALL message( 'pmci_check_grid_matching', 'PMC0018', 3, 2, 0, 6, 0 )
       ENDIF

       IF ( non_int_gsr_y > 0 )  THEN
          message_string = 'nesting grid-spacing ratio ( parent dy / child dy ) ' //               &
                           'must have an integer value'
          CALL message( 'pmci_check_grid_matching', 'PMC0018', 3, 2, 0, 6, 0 )
       ENDIF

       IF ( non_int_gsr_z > 0 )  THEN
          message_string = 'nesting grid-spacing ratio ( parent dz / child dz ) ' //               &
                           'must have an integer value for each z-level'
          CALL message( 'pmci_check_grid_matching', 'PMC0018', 3, 2, 0, 6, 0 )
       ENDIF

       IF ( too_narrow_pesd_x > 0 )  THEN
          message_string = 'child subdomain width in x-direction must not be smaller than its ' // &
                           'parent grid dx. &Change the PE-grid setting (npex, npey) to satisfy' //&
                           ' this requirement.'
          CALL message( 'pmci_check_grid_matching', 'PMC0019', 3, 2, 0, 6, 0 )
       ENDIF

       IF ( too_narrow_pesd_y > 0 )  THEN
          message_string = 'child subdomain width in y-direction must not be smaller than its ' // &
                           'parent grid dy. &Change the PE-grid setting (npex, npey) to satisfy' //&
                           ' this requirement.'
          CALL message( 'pmci_check_grid_matching', 'PMC0019', 3, 2, 0, 6, 0 )
       ENDIF

       IF ( parent_nest_shift_z > nest_shift_z )  THEN
          WRITE( message_string, '(A,F8.1,A,F8.1,A)' )  'heigth of child bottom boundary (',       &
             nest_shift_z, ' m) must not lie below the bottom boundary of the respective ' //      &
             'parent (', parent_nest_shift_z, ' m)'
          CALL message( 'pmci_check_grid_matching', 'PMC0031', 3, 2, 0, 6, 0 )
       ENDIF

    ENDIF  !  ( myid == 0 )

 END SUBROUTINE pmci_check_grid_matching


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the areas of the domain boundary faces: left, right, south, north and top.
!> For cyclic boundaries face areas are zero along the respective direction.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_compute_face_areas

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< Running index in the x-direction
    INTEGER(iwp) ::  j       !< Running index in the y-direction
    INTEGER(iwp) ::  k       !< Running index in the z-direction
    INTEGER(iwp) ::  n       !< Running index over boundary faces

    REAL(wp) ::  face_area_local  !< Local (for the current pe) integral face area of the left boundary
    REAL(wp) ::  sub_sum          !< Intermediate sum in order to improve the accuracy of the summation

!
!-- Sum up the volume flow face area of the left boundary
    face_area(1) = 0.0_wp
    face_area_local = 0.0_wp
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       i = 0
       DO  j = nys, nyn
          sub_sum = 0.0_wp
          DO  k = nzb + 1, nzt
             sub_sum = sub_sum + dzw(k) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
          face_area_local =  face_area_local + dy * sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( face_area_local, face_area(1), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    face_area(1) = face_area_local
#endif
!
!-- Sum up the volume flow face area of the right boundary
    face_area(2) = 0.0_wp
    face_area_local = 0.0_wp
    IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
       i = nx
       DO  j = nys, nyn
          sub_sum = 0.0_wp
          DO  k = nzb + 1, nzt
             sub_sum = sub_sum + dzw(k) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
          face_area_local =  face_area_local + dy * sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( face_area_local, face_area(2), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    face_area(2) = face_area_local
#endif
!
!-- Sum up the volume flow face area of the south boundary
    face_area(3) = 0.0_wp
    face_area_local = 0.0_wp
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       j = 0
       DO  i = nxl, nxr
          sub_sum = 0.0_wp
          DO  k = nzb + 1, nzt
             sub_sum = sub_sum + dzw(k) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO
          face_area_local = face_area_local + dx * sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( face_area_local, face_area(3), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    face_area(3) = face_area_local
#endif
!
!-- Sum up the volume flow face area of the north boundary
    face_area(4) = 0.0_wp
    face_area_local = 0.0_wp
    IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
       j = ny
       DO  i = nxl, nxr
          sub_sum = 0.0_wp
          DO  k = nzb + 1, nzt
             sub_sum = sub_sum + dzw(k) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO
          face_area_local = face_area_local + dx * sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( face_area_local, face_area(4), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    face_area(4) = face_area_local
#endif
!
!-- Sum up the volume flow face area of the top boundary
    face_area(5) = 0.0_wp
    face_area_local = 0.0_wp
    k = nzt
    DO  i = nxl, nxr
       DO  j = nys, nyn
          face_area_local = face_area_local +                                                      &
                            dx * dy * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
       ENDDO
    ENDDO

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( face_area_local, face_area(5), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    face_area(5) = face_area_local
#endif

!
!-- Sum up the volume flow face area of the bottom boundary
    face_area(6) = 0.0_wp
    face_area_local = 0.0_wp
    k = 0
    DO  i = nxl, nxr
       DO  j = nys, nyn
          face_area_local = face_area_local +                                                      &
                            dx * dy * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
       ENDDO
    ENDDO

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( face_area_local, face_area(6), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    face_area(6) = face_area_local
#endif

!
!-- The 7th element is used for the total area
    face_area(7) = 0.0_wp
    DO  n = 1, 6
       face_area(7) = face_area(7) + face_area(n)
    ENDDO

 END SUBROUTINE pmci_compute_face_areas


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define the anterpolation starting kp-indices for all jp,ip based on the obstacle-canopy
!> topograhy for canopy-restricted anterpolation. Note that this is based on the child terrain
!> topography information since it is difficult to access the parent topography information
!> from the child. This means that these topographies are assumed to be close to each other.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_compute_kpb_anterp

    INTEGER(iwp) ::  building_node          !< Shorthand for building flag: inside building 1 otherwise 0
    INTEGER(iwp) ::  ic                     !< Child-grid index in the x-direction
    INTEGER(iwp) ::  ip                     !< Parent-grid index in the x-direction
    INTEGER(iwp) ::  jc                     !< Child-grid index in the y-direction
    INTEGER(iwp) ::  jp                     !< Parent-grid index in the y-direction
    INTEGER(iwp) ::  kc                     !< Child-grid index in the z-direction
    INTEGER(iwp) ::  num_building_nodes     !< Number of building nodes in the whole child domain
    INTEGER(iwp) ::  num_building_nodes_l   !< Number of building nodes in the current subdomain
    INTEGER(iwp) ::  terrain_surf_k_child   !< Local terrain height k index in the child grid
    INTEGER(iwp) ::  terrain_surf_k_parent  !< Local terrain height k index in the parent grid


    ALLOCATE( kpb_anterp(jps:jpn,ipl:ipr) )
    IF ( nested_bc_at_bottom )  THEN
       kpb_anterp(jps:jpn,ipl:ipr) = 1 + anterpolation_buffer_width
    ELSE
       kpb_anterp(jps:jpn,ipl:ipr) = 0
    ENDIF

    IF ( nesting_mode /= 'one-way'  .AND.  topography /= 'flat' )  THEN
!
!--    Check if the value was given by user or not
       IF ( anterpolation_starting_height > 9000000.0_wp )  THEN
          anterpolation_starting_height = 0.0_wp
          RETURN
       ENDIF
!
!--    Determine the number of building-canopy grid nodes in this child domain. 
!--    If zero, canopy-restriction of anterpolation makes no sense and kpb_anterp
!--    is left to its initial values set above.
       num_building_nodes_l = 0
       num_building_nodes   = 0
       DO  ic = nxl, nxr
          DO  jc = nys, nyn
             DO  kc = nzb, nzt
                building_node = MERGE( 1, 0, BTEST( topo_flags(kc,jc,ic), 6 ) )
                num_building_nodes_l = num_building_nodes_l + building_node
             ENDDO
          ENDDO
       ENDDO

       CALL MPI_ALLREDUCE( num_building_nodes_l, num_building_nodes, 1, MPI_INT, MPI_SUM, comm2d,  &
                           ierr )
!  
!--    Test if there is any building canopy in this child domain and only compute 
!--    kpb_anterp if building canopy exists
       IF ( num_building_nodes > 0 )  THEN
          DO  ip = ipl, ipr
             DO  jp = jps, jpn
                terrain_surf_k_parent = 0
                DO  ic = iflo(ip), ifuo(ip)
                   DO  jc = jflo(jp), jfuo(jp)
                      terrain_surf_k_child = MINLOC(                                               &
                           MERGE( 1, 0, BTEST( topo_flags(:,jc,ic), 5 ) ), DIM = 1 ) - 1
                      terrain_surf_k_parent = MAX( terrain_surf_k_child, terrain_surf_k_parent )
                   ENDDO
                ENDDO
                terrain_surf_k_parent = NINT( terrain_surf_k_parent / REAL( kgsr, wp ) )
                kpb_anterp(jp,ip) = terrain_surf_k_parent +                                        &
                     NINT( anterpolation_starting_height / pg%dz )
             ENDDO
          ENDDO
       ENDIF

    ENDIF

 END SUBROUTINE pmci_compute_kpb_anterp



#endif

 END SUBROUTINE pmci_setup_child


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Translate the coordinates of the current domain into the root coordinate system and store them
!> (z remains the same).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_setup_coordinates

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index in the x-direction
    INTEGER(iwp) ::  j  !< grid index in the y-direction

!
!-- Create coordinate arrays.
    ALLOCATE( coord_x(-nbgp:nx+nbgp) )
    ALLOCATE( coord_y(-nbgp:ny+nbgp) )

    DO  i = -nbgp, nx + nbgp
       coord_x(i) = lower_left_coord_x + i * dx
    ENDDO

    DO  j = -nbgp, ny + nbgp
       coord_y(j) = lower_left_coord_y + j * dy
    ENDDO

#endif

 END SUBROUTINE pmci_setup_coordinates


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> In this subroutine the number of coupled arrays is determined.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_num_arrays

#if defined( __parallel )
    IMPLICIT NONE


    IF ( atmosphere_ocean_coupled_run )  THEN
       pmc_max_array = 4
       RETURN
    ENDIF

!
!-- The number of coupled arrays depends on the model settings. At least 5 arrays need to be
!-- coupled (u, v, w, e, diss).  Please note, actually e and diss (TKE and dissipation rate) are
!-- only required if RANS-RANS nesting is applied, but memory is allocated nevertheless. This is
!-- because the information whether they are needed or not is retrieved at a later point in time.
!-- In case e and diss are not needed, they are also not exchanged between parent and child.
    pmc_max_array = 5
!
!-- pt
    IF ( .NOT. neutral )  pmc_max_array = pmc_max_array + 1

    IF ( humidity )  THEN
!
!--    q
       pmc_max_array = pmc_max_array + 1
!
!--    qc, nc
       IF ( bulk_cloud_model  .AND.  microphysics_morrison )                                       &
          pmc_max_array = pmc_max_array + 2
!
!--    qr, nr
       IF ( bulk_cloud_model  .AND.  microphysics_seifert )                                        &
          pmc_max_array = pmc_max_array + 2
    ENDIF
!
!-- s
    IF ( passive_scalar )  pmc_max_array = pmc_max_array + 1
!
!-- sa
    IF ( ocean_mode  .AND.  salinity )  pmc_max_array = pmc_max_array + 1
!
!-- nr_part, part_adr
    IF ( particle_advection )  pmc_max_array = pmc_max_array + 2
!
!-- Chemistry, depends on number of species
    IF ( air_chemistry  .AND.  nesting_chem )  pmc_max_array = pmc_max_array + nspec
!
!-- SALSA, depens on the number aerosol size bins and chemical components + the number of default
!-- gases
    IF ( salsa  .AND.  nesting_salsa )  pmc_max_array = pmc_max_array + nbins_aerosol +            &
                                                        nbins_aerosol * ncomponents_mass
    IF ( .NOT. salsa_gases_from_chem )  pmc_max_array = pmc_max_array + ngases_salsa

#endif

 END SUBROUTINE pmci_num_arrays



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Assigns the pointer to the array to be coupled.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_set_array_pointer( name, child_id, nz_child, ks_child, ke_child, n)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name  !<

    INTEGER(iwp), INTENT(IN) ::  child_id  !<
    INTEGER(iwp), INTENT(IN) ::  nz_child  !<

    INTEGER(iwp), INTENT(IN), OPTIONAL ::  n  !< index of chemical species
    INTEGER(iwp), INTENT(IN), OPTIONAL ::  ks_child  !<
    INTEGER(iwp), INTENT(IN), OPTIONAL ::  ke_child  !<

#if defined( __parallel )
!
!-- Local variables:
    INTEGER(iwp) ::  ierr  !< MPI error code

    INTEGER(idp), POINTER, DIMENSION(:,:) ::  i_2d  !<

    REAL(wp), POINTER, DIMENSION(:,:)   ::  p_2d      !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  p_3d      !<
    REAL(wp), POINTER, DIMENSION(:,:,:) ::  p_3d_sec  !<


    NULLIFY( p_3d )
    NULLIFY( p_2d )
    NULLIFY( i_2d )
!
!-- List of array names, which can be coupled.
!-- In case of 3D please change also the second array for the pointer version
    IF ( TRIM(name) == 'u'          )  p_3d => u
    IF ( TRIM(name) == 'v'          )  p_3d => v
    IF ( TRIM(name) == 'w'          )  p_3d => w
    IF ( TRIM(name) == 'e'          )  p_3d => e
    IF ( TRIM(name) == 'pt'         )  p_3d => pt
    IF ( TRIM(name) == 'q'          )  p_3d => q
    IF ( TRIM(name) == 'qc'         )  p_3d => qc
    IF ( TRIM(name) == 'qr'         )  p_3d => qr
    IF ( TRIM(name) == 'nr'         )  p_3d => nr
    IF ( TRIM(name) == 'nc'         )  p_3d => nc
    IF ( TRIM(name) == 's'          )  p_3d => s
    IF ( TRIM(name) == 'sa'         )  p_3d => sa
    IF ( TRIM(name) == 'diss'       )  p_3d => diss
    IF ( TRIM(name) == 'nr_part'    )  i_2d => nr_part
    IF ( TRIM(name) == 'part_adr'   )  i_2d => part_adr
    IF ( INDEX( TRIM(name), 'chem_' ) /= 0 )  p_3d => chem_species(n)%conc
    IF ( INDEX( TRIM(name), 'an_' ) /= 0 )  p_3d => aerosol_number(n)%conc
    IF ( INDEX( TRIM(name), 'am_' ) /= 0 )  p_3d => aerosol_mass(n)%conc
    IF ( INDEX( TRIM(name), 'sg_' ) /= 0  .AND.  .NOT. salsa_gases_from_chem )                     &
       p_3d => salsa_gas(n)%conc

!
!-- Ocean coupling
    IF ( TRIM(name) == 'exchange_array_1' )  p_2d => surface_coupler_exchange_array_1
    IF ( TRIM(name) == 'exchange_array_2' )  p_2d => surface_coupler_exchange_array_2
    IF ( TRIM(name) == 'exchange_array_3' )  p_2d => surface_coupler_exchange_array_3
    IF ( TRIM(name) == 'exchange_array_4' )  p_2d => surface_coupler_exchange_array_4

!
!-- Next line is just an example for a 2D array (not active for coupling!)
!-- Please note, that z0 has to be declared as TARGET array in modules.f90.
!    IF ( TRIM(name) == 'z0' )    p_2d => z0
    IF ( TRIM(name) == 'u'    )  p_3d_sec => u_2
    IF ( TRIM(name) == 'v'    )  p_3d_sec => v_2
    IF ( TRIM(name) == 'w'    )  p_3d_sec => w_2
    IF ( TRIM(name) == 'e'    )  p_3d_sec => e_2
    IF ( TRIM(name) == 'pt'   )  p_3d_sec => pt_2
    IF ( TRIM(name) == 'q'    )  p_3d_sec => q_2
    IF ( TRIM(name) == 'qc'   )  p_3d_sec => qc_2
    IF ( TRIM(name) == 'qr'   )  p_3d_sec => qr_2
    IF ( TRIM(name) == 'nr'   )  p_3d_sec => nr_2
    IF ( TRIM(name) == 'nc'   )  p_3d_sec => nc_2
    IF ( TRIM(name) == 's'    )  p_3d_sec => s_2
    IF ( TRIM(name) == 'sa'   )  p_3d_sec => sa_2
    IF ( TRIM(name) == 'diss' )  p_3d_sec => diss_2
    IF ( INDEX( TRIM(name), 'chem_' ) /= 0 )  p_3d_sec => spec_conc_2(:,:,:,n)
    IF ( INDEX( TRIM(name), 'an_' )   /= 0 )  p_3d_sec => nconc_2(:,:,:,n)
    IF ( INDEX( TRIM(name), 'am_' )   /= 0 )  p_3d_sec => mconc_2(:,:,:,n)
    IF ( INDEX( TRIM(name), 'sg_' )   /= 0  .AND.  .NOT.  salsa_gases_from_chem )                  &
       p_3d_sec => gconc_2(:,:,:,n)

    IF ( ASSOCIATED( p_3d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, p_3d, nz_child, nz, array_2 = p_3d_sec, ks_cl=ks_child, ke_cl = ke_child )
    ELSEIF  ( ASSOCIATED( p_2d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, p_2d )
    ELSEIF  ( ASSOCIATED( i_2d ) )  THEN
       CALL pmc_s_set_dataarray( child_id, i_2d )
    ELSE
!
!--    Give only one message for the root domain
       IF ( root_model  .AND.  myid == 0 )  THEN
          message_string = 'pointer for array "' // TRIM( name ) // '" can''t be associated ' //   &
                           'for root domain'
          CALL message( 'pmci_set_array_pointer', 'PMC0020', 3, 2, 0, 6, 0 )
       ELSE
!
!--       Avoid others to continue
          CALL MPI_BARRIER( comm2d, ierr )
       ENDIF

    ENDIF

#endif

 END SUBROUTINE pmci_set_array_pointer



 INTEGER FUNCTION get_number_of_children()

    IMPLICIT NONE


#if defined( __parallel )
    get_number_of_children = SIZE( pmc_parent_for_child ) - 1
#else
    get_number_of_children = 0
#endif

    RETURN

 END FUNCTION get_number_of_children



 INTEGER FUNCTION get_childid( id_index )

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  id_index   !<


#if defined( __parallel )
    get_childid = pmc_parent_for_child(id_index)
#else
    get_childid = 0
#endif

    RETURN

 END FUNCTION get_childid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Provides the child grid edge coordinates for pmc_particle_interface.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE get_child_edges( m, lx_coord, lx_coord_b, rx_coord, rx_coord_b, sy_coord, sy_coord_b,  &
                             ny_coord, ny_coord_b, uz_coord, uz_coord_b )

    IMPLICIT NONE

    INTEGER,INTENT(IN) ::  m  !<

    REAL(wp),INTENT(OUT) ::  lx_coord, lx_coord_b  !<
    REAL(wp),INTENT(OUT) ::  ny_coord, ny_coord_b  !<
    REAL(wp),INTENT(OUT) ::  rx_coord, rx_coord_b  !<
    REAL(wp),INTENT(OUT) ::  sy_coord, sy_coord_b  !<
    REAL(wp),INTENT(OUT) ::  uz_coord, uz_coord_b  !<


#if defined( __parallel )

    lx_coord = childgrid(m)%lx_coord
    rx_coord = childgrid(m)%rx_coord
    sy_coord = childgrid(m)%sy_coord
    ny_coord = childgrid(m)%ny_coord
    uz_coord = childgrid(m)%uz_coord

    lx_coord_b = childgrid(m)%lx_coord_b
    rx_coord_b = childgrid(m)%rx_coord_b
    sy_coord_b = childgrid(m)%sy_coord_b
    ny_coord_b = childgrid(m)%ny_coord_b
    uz_coord_b = childgrid(m)%uz_coord_b

#endif

 END SUBROUTINE get_child_edges



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Provides the child grid spacings for pmc_particle_interface.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE  get_child_gridspacing( m, dx, dy,dz )

    IMPLICIT NONE

    INTEGER, INTENT(IN) ::  m  !<

    REAL(wp), INTENT(OUT) ::  dx,dy  !<

    REAL(wp), INTENT(OUT), OPTIONAL ::  dz  !<


#if defined( __parallel )

    dx = childgrid(m)%dx
    dy = childgrid(m)%dy
    IF ( PRESENT( dz ) )  THEN
       dz = childgrid(m)%dz
    ENDIF

#endif

 END SUBROUTINE get_child_gridspacing


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate child's parent-grid arrays for data transfer between parent and child.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_create_childs_parent_grid_arrays( name, is, ie, js, je, nzc, n  )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  name  !<

    INTEGER(iwp), INTENT(IN) ::  ie   !<  RENAME ie, is, je, js?
    INTEGER(iwp), INTENT(IN) ::  is   !<
    INTEGER(iwp), INTENT(IN) ::  je   !<
    INTEGER(iwp), INTENT(IN) ::  js   !<
    INTEGER(iwp), INTENT(IN) ::  nzc  !<  nzc is pg%nz, but note that pg%nz is not the original nz of parent,
                                      !<  but the highest parent-grid level needed for nesting.
    INTEGER(iwp), INTENT(IN), OPTIONAL ::  n  !< number of chemical species / salsa variables

#if defined( __parallel )
!
!-- Local variables:
    INTEGER(iwp) ::  ierr  !<

    INTEGER(idp), POINTER,DIMENSION(:,:) ::  i_2d  !<

    REAL(wp), POINTER,DIMENSION(:,:) ::  p_2d  !<

    REAL(wp), POINTER,DIMENSION(:,:,:) ::  p_3d  !<

    NULLIFY( p_3d )
    NULLIFY( p_2d )
    NULLIFY( i_2d )
!
!-- List of array names, which can be coupled
    IF ( TRIM( name ) == 'u' )  THEN
       IF ( .NOT. ALLOCATED( uc ) )  ALLOCATE( uc(0:nzc+1,js:je,is:ie) )
       p_3d => uc
    ELSEIF ( TRIM( name ) == 'v' )  THEN
       IF ( .NOT. ALLOCATED( vc ) )  ALLOCATE( vc(0:nzc+1,js:je,is:ie) )
       p_3d => vc
    ELSEIF ( TRIM( name ) == 'w' )  THEN
       IF ( .NOT. ALLOCATED( wc ) )  ALLOCATE( wc(0:nzc+1,js:je,is:ie) )
       p_3d => wc
    ELSEIF ( TRIM( name ) == 'e' )  THEN
       IF ( .NOT. ALLOCATED( ec ) )  ALLOCATE( ec(0:nzc+1,js:je,is:ie) )
       p_3d => ec
    ELSEIF ( TRIM( name ) == 'diss' )  THEN
       IF ( .NOT. ALLOCATED( dissc ) )  ALLOCATE( dissc(0:nzc+1,js:je,is:ie) )
       p_3d => dissc
    ELSEIF ( TRIM( name ) == 'pt')  THEN
       IF ( .NOT. ALLOCATED( ptc ) )  ALLOCATE( ptc(0:nzc+1,js:je,is:ie) )
       p_3d => ptc
    ELSEIF ( TRIM( name ) == 'q')  THEN
       IF ( .NOT. ALLOCATED( q_c ) )  ALLOCATE( q_c(0:nzc+1,js:je,is:ie) )
       p_3d => q_c
    ELSEIF ( TRIM( name ) == 'qc')  THEN
       IF ( .NOT. ALLOCATED( qcc ) )  ALLOCATE( qcc(0:nzc+1,js:je,is:ie) )
       p_3d => qcc
    ELSEIF ( TRIM( name ) == 'qr')  THEN
       IF ( .NOT. ALLOCATED( qrc ) )  ALLOCATE( qrc(0:nzc+1,js:je,is:ie) )
       p_3d => qrc
    ELSEIF ( TRIM( name ) == 'nr')  THEN
       IF ( .NOT. ALLOCATED( nrc ) )  ALLOCATE( nrc(0:nzc+1,js:je,is:ie) )
       p_3d => nrc
    ELSEIF ( TRIM( name ) == 'nc')  THEN
       IF ( .NOT. ALLOCATED( ncc ) )  ALLOCATE( ncc(0:nzc+1,js:je,is:ie) )
       p_3d => ncc
    ELSEIF ( TRIM( name ) == 's')  THEN
       IF ( .NOT. ALLOCATED( sc ) )  ALLOCATE( sc(0:nzc+1,js:je,is:ie) )
       p_3d => sc
    ELSEIF ( TRIM( name ) == 'sa')  THEN
       IF ( .NOT. ALLOCATED( sac ) )  ALLOCATE( sac(0:nzc+1,js:je,is:ie) )
       p_3d => sac
    ELSEIF ( TRIM( name ) == 'nr_part') THEN
       IF ( .NOT. ALLOCATED( nr_partc ) )  ALLOCATE( nr_partc(js:je,is:ie) )
       i_2d => nr_partc
    ELSEIF ( TRIM( name ) == 'part_adr') THEN
       IF ( .NOT. ALLOCATED( part_adrc ) )  ALLOCATE( part_adrc(js:je,is:ie) )
       i_2d => part_adrc
    ELSEIF ( TRIM( name ) == 'exchange_array_1' )  THEN
       IF ( .NOT. ALLOCATED( surface_coupler_exchange_array_1c ) )  THEN
          ALLOCATE( surface_coupler_exchange_array_1c(js:je,is:ie) )
       ENDIF
       p_2d => surface_coupler_exchange_array_1c
    ELSEIF ( TRIM( name ) == 'exchange_array_2' )  THEN
       IF ( .NOT. ALLOCATED( surface_coupler_exchange_array_2c ) )  THEN
          ALLOCATE( surface_coupler_exchange_array_2c(js:je,is:ie) )
       ENDIF
       p_2d => surface_coupler_exchange_array_2c
    ELSEIF ( TRIM( name ) == 'exchange_array_3' )  THEN
       IF ( .NOT. ALLOCATED( surface_coupler_exchange_array_3c ) )  THEN
          ALLOCATE( surface_coupler_exchange_array_3c(js:je,is:ie) )
       ENDIF
       p_2d => surface_coupler_exchange_array_3c
    ELSEIF ( TRIM( name ) == 'exchange_array_4' )  THEN
       IF ( .NOT. ALLOCATED( surface_coupler_exchange_array_4c ) )  THEN
          ALLOCATE( surface_coupler_exchange_array_4c(js:je,is:ie) )
       ENDIF
       p_2d => surface_coupler_exchange_array_4c
    ELSEIF ( TRIM( name(1:5) ) == 'chem_' )  THEN
       IF ( .NOT. ALLOCATED( chem_spec_c ) ) ALLOCATE( chem_spec_c(0:nzc+1,js:je,is:ie,1:nspec) )
       p_3d => chem_spec_c(:,:,:,n)
    ELSEIF ( TRIM( name(1:3) ) == 'an_' )  THEN
       IF ( .NOT. ALLOCATED( aerosol_number_c ) )                                                  &
          ALLOCATE( aerosol_number_c(0:nzc+1,js:je,is:ie,1:nbins_aerosol) )
       p_3d => aerosol_number_c(:,:,:,n)
    ELSEIF ( TRIM( name(1:3) ) == 'am_' )  THEN
       IF ( .NOT. ALLOCATED( aerosol_mass_c ) )                                                    &
          ALLOCATE( aerosol_mass_c(0:nzc+1,js:je,is:ie,1:(nbins_aerosol*ncomponents_mass) ) )
       p_3d => aerosol_mass_c(:,:,:,n)
    ELSEIF ( TRIM( name(1:3) ) == 'sg_'  .AND.  .NOT. salsa_gases_from_chem )  THEN
       IF ( .NOT. ALLOCATED( salsa_gas_c ) )                                                       &
          ALLOCATE( salsa_gas_c(0:nzc+1,js:je,is:ie,1:ngases_salsa) )
       p_3d => salsa_gas_c(:,:,:,n)
    !ELSEIF (trim(name) == 'z0') then
       !IF (.not.allocated(z0c))  allocate(z0c(js:je, is:ie))
       !p_2d => z0c
    ENDIF

    IF ( ASSOCIATED( p_3d ) )  THEN
       CALL pmc_c_set_dataarray( p_3d )
    ELSEIF ( ASSOCIATED( p_2d ) )  THEN
       CALL pmc_c_set_dataarray( p_2d )
    ELSEIF ( ASSOCIATED( i_2d ) )  THEN
       CALL pmc_c_set_dataarray( i_2d )
    ELSE
!
!--    Give only one message for the first child domain.
       IF ( cpl_id == 2  .AND.  myid == 0 )  THEN
          message_string = 'pointer for array "' // TRIM( name ) // '" can''t be associated ' //   &
                           'for child domain'
          CALL message( 'pmci_create_childs_parent_grid_arrays', 'PMC0020', 3, 2, 0, 6, 0 )
       ELSE
!
!--          Prevent others from continuing in case the abort is to come.
          CALL MPI_BARRIER( comm2d, ierr )
       ENDIF

    ENDIF

#endif
 END SUBROUTINE pmci_create_childs_parent_grid_arrays



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Send data for the children in order to let them create initial conditions by interpolating the
!> parent-domain fields.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_parent_initialize

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  child_id  !< interal ID of the corresponding child domain
    INTEGER(iwp) ::  m         !< running index over all child domains

    REAL(wp) ::  waittime  !<


    IF ( atmosphere_ocean_coupled_run )  RETURN

    DO  m = 1, SIZE( pmc_parent_for_child ) - 1
       IF ( homogeneous_initialization_child )  THEN
!
!--       Compute domain-averaged profiles in the parent and send them to the child.
!--       These are later used to initialize the corresponding variable in the child-domain.
          child_id = pmc_parent_for_child(m)
          CALL pmci_send_domain_averaged_profiles( child_id )
       ELSE
          child_id = pmc_parent_for_child(m)
          CALL pmc_s_fillbuffer( child_id, waittime=waittime )
       ENDIF
    ENDDO

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute and send domain-averaged profiles of prognostic variables within the parent domain and
!> send these to the respective child domain.
!> Attention: If the order of the send operations is altered it needs to be altered in
!>            pmci_recv_domain_averaged_profiles the same way!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_send_domain_averaged_profiles( child_id )

    INTEGER(iwp) ::  child_id !< interal ID of the corresponding child domain
    INTEGER(iwp) ::  ierr     !< MPI error code
    INTEGER(iwp) ::  n        !< running index for species (chemistry or aerosols)
    INTEGER(iwp) ::  nz_child !< number of vertical grid points in parent covering vertical depth of child domain
    INTEGER(iwp) ::  tag_nr   !< tag number to identify send/receive operations

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mean_profile !< domain-averaged profile in the parent domain

    nz_child = childgrid(m)%nz
!
!-- Compute domain-average u/v-profile in the parent. These are later used to initialize u/v in the
!-- child-domain.
!-- u-component
    ALLOCATE( mean_profile(nzb:nz_child+1) )
!
!-- Start with tag number 50 (far beyond other used tag numbers).
    tag_nr = 50
    CALL pmci_compute_average( u, m, 1, nz_child, mean_profile, 'u' )
    IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
    tag_nr = tag_nr + 1
!
!-- v-component
    CALL pmci_compute_average( v, m, 2, nz_child, mean_profile, 'v' )
    IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
    tag_nr = tag_nr + 1
!
!-- SGS-TKE
    IF ( ( rans_mode_parent  .AND.         rans_mode )  .OR.                                       &
         ( .NOT. rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                                  &
           .NOT. constant_diffusion ) )  THEN
       CALL pmci_compute_average( e, m, 0, nz_child, mean_profile, 's' )
       IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
       tag_nr = tag_nr + 1
    ENDIF
!
!-- dissipation rate
    IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
       CALL pmci_compute_average( diss, m, 0, nz_child, mean_profile, 's' )
       IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
       tag_nr = tag_nr + 1
    ENDIF
!
!-- potential temperature
    IF ( .NOT. neutral )  THEN
       CALL pmci_compute_average( pt, m, 0, nz_child, mean_profile, 'pt' )
       IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
       tag_nr = tag_nr + 1
    ENDIF

    IF ( humidity )  THEN
!
!--    mixing ratio
       CALL pmci_compute_average( q, m, 0, nz_child, mean_profile, 'q' )
       IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
       tag_nr = tag_nr + 1
!
!--    qc and nc
       IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
!
!--       qc
          CALL pmci_compute_average( qc, m, 0, nz_child, mean_profile, 's' )
          IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
!
!--       nc
          CALL pmci_compute_average( nc, m, 0, nz_child, mean_profile, 's' )
          IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDIF
!
!--    qr and nr
       IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
!
!--       qr
          CALL pmci_compute_average( qr, m, 0, nz_child, mean_profile, 's' )
          IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
!
!--       nr
          CALL pmci_compute_average( nr, m, 0, nz_child, mean_profile, 's' )
          IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDIF
    ENDIF
!
!-- passive scalar
    IF ( passive_scalar )  THEN
       CALL pmci_compute_average( s, m, 0, nz_child, mean_profile, 's' )
       IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
       tag_nr = tag_nr + 1
    ENDIF
!
!-- chemistry
    IF ( air_chemistry  .AND.  nesting_chem )  THEN
       DO  n = 1, nspec
          CALL pmci_compute_average( chem_species(n)%conc, m, 0, nz_child, mean_profile, 's' )
          IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDDO
    ENDIF
!
!-- aerosols
    IF ( salsa  .AND.  nesting_salsa )  THEN
       DO  n = 1, nbins_aerosol
          CALL pmci_compute_average( aerosol_number(n)%conc, m, 0, nz_child, mean_profile, 's' )
          IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDDO
       DO  n = 1, nbins_aerosol * ncomponents_mass
          CALL pmci_compute_average( aerosol_mass(n)%conc, m, 0, nz_child, mean_profile, 's' )
          IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDDO
       IF ( .NOT. salsa_gases_from_chem )  THEN
          DO  n = 1, ngases_salsa
             CALL pmci_compute_average( salsa_gas(n)%conc, m, 0, nz_child, mean_profile, 's' )
             IF ( myid == 0 )  CALL pmc_send_to_child( child_id, mean_profile, nz_child+2, 0, tag_nr, ierr )
             tag_nr = tag_nr + 1
          ENDDO
       ENDIF

    ENDIF

    DEALLOCATE( mean_profile )

    END SUBROUTINE pmci_send_domain_averaged_profiles

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute domain-averaged profile of of the given array.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE pmci_compute_average( ar, m, bit_nr, nz_child, mean_profile, var )

       CHARACTER(LEN=*) ::  var  !< variable string.

       INTEGER(iwp) ::  bit_nr   !< bit position used to masked topography for the respective variable
       INTEGER(iwp) ::  i        !< running index in x-direction
       INTEGER(iwp) ::  ierr     !< MPI error code
       INTEGER(iwp) ::  j        !< running index in y-direction
       INTEGER(iwp) ::  k        !< running index in z-direction
       INTEGER(iwp) ::  m        !< index for current child in childgrid
       INTEGER(iwp) ::  nz_child !< number of vertical grid points in parent covering vertical depth of child domain

       INTEGER(iwp), DIMENSION(nzb:nz_child+1) ::  num_gp !< domain-averaged profile in the parent domain

       LOGICAL ::  within_child_area !< flag indicating that current grid point location also belongs to a child domain

       REAL(wp), DIMENSION(nzb:nz_child+1), INTENT(INOUT)             ::  mean_profile !< domain-averaged profile in the parent domain

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  ar           !< array that will be domain-averaged

!
!--    Branch with masking of topography for all quantities except for potential temperature and
!--    mixing ratio.
       mean_profile = 0.0_wp
       num_gp = 0_iwp
       IF ( var /= 'pt'  .AND.  var /= 'q' )  THEN
          DO  i = nxl, nxr
             DO  j = nys, nyn
                within_child_area = x(i) + lower_left_coord_x >= childgrid(m)%lx_coord  .AND.      &
                                    x(i) + lower_left_coord_x <= childgrid(m)%rx_coord  .AND.      &
                                    y(j) + lower_left_coord_y >= childgrid(m)%sy_coord  .AND.      &
                                    y(j) + lower_left_coord_y <= childgrid(m)%ny_coord

                DO  k = nzb, nz_child+1
                   mean_profile(k) = mean_profile(k) + ar(k,j,i)                                   &
                           * MERGE( 1.0_wp, 0.0_wp,  BTEST( topo_flags(k,j,i), bit_nr ) )          &
                           * MERGE( 1.0_wp, 0.0_wp, within_child_area )

                   num_gp(k) =  num_gp(k)                                                          &
                           + MERGE( 1_iwp, 0_iwp, BTEST( topo_flags(k,j,i), bit_nr ) )             &
                           * MERGE( 1_iwp, 0_iwp, within_child_area )
                ENDDO
             ENDDO
          ENDDO
!
!--    No masking of topography for potential temperature and mixing ratio. This is to avoid zero
!--    values in the profiles when the parent's domain lowest prognostic level is higher than nzb+1.
!--    Especially for potential temperature and humidity this can become critical since zero
!--    values will cause divisions by zero in the buoyancy term. Anyhow, since both quantities
!--    have been defined also within topography, this will have no consequence at all.
       ELSE
          DO  i = nxl, nxr
             DO  j = nys, nyn
                within_child_area = x(i) + lower_left_coord_x >= childgrid(m)%lx_coord  .AND.      &
                                    x(i) + lower_left_coord_x <= childgrid(m)%rx_coord  .AND.      &
                                    y(j) + lower_left_coord_y >= childgrid(m)%sy_coord  .AND.      &
                                    y(j) + lower_left_coord_y <= childgrid(m)%ny_coord

                DO  k = nzb, nz_child+1
                   mean_profile(k) = mean_profile(k) + ar(k,j,i)                                   &
                           * MERGE( 1.0_wp, 0.0_wp, within_child_area )

                   num_gp(k) =  num_gp(k)                                                          &
                           + MERGE( 1_iwp, 0_iwp, within_child_area )
                ENDDO
             ENDDO
          ENDDO
       ENDIF

       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, mean_profile, nz_child+1-nzb+1, MPI_REAL,    MPI_SUM,     &
                           comm2d, ierr )
       CALL MPI_ALLREDUCE( MPI_IN_PLACE, num_gp,       nz_child+1-nzb+1, MPI_INTEGER, MPI_SUM,     &
                           comm2d, ierr )
!
!--    Divide by number of grid points. Avoid division by zero.
       num_gp = MERGE( num_gp, 1_iwp, num_gp >= 1_iwp )
       mean_profile = mean_profile / REAL( num_gp, KIND = wp )

    END SUBROUTINE pmci_compute_average

#endif
 END SUBROUTINE pmci_parent_initialize


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Create initial conditions for the current child domain by interpolating the parent-domain fields.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_child_initialize

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  ic  !< Child-grid index in x-direction
    INTEGER(iwp) ::  jc  !< Child-grid index in y-direction
    INTEGER(iwp) ::  kc  !< Child-grid index in z-direction
    INTEGER(iwp) ::  lb  !< Running index for aerosol size bins
    INTEGER(iwp) ::  lc  !< Running index for aerosol mass bins
    INTEGER(iwp) ::  lg  !< Running index for salsa gases
    INTEGER(iwp) ::  n   !< Running index for chemical species

    REAL(wp) ::  waittime  !< Waiting time

!
!-- Root model is never anyone's child.
    IF ( .NOT.  root_model .AND. .NOT. atmosphere_ocean_coupled_run )  THEN
!
!--    1d initialization with domain-averaged profiles from the parent domain.
       IF ( homogeneous_initialization_child )  THEN
!
!--       Receive domain-adveraged u-/v-profiles from the parent domain and distribute them
!--       among comm2d.
          CALL pmci_recv_domain_averaged_profiles
!
!--       Initialize variables with 1d profiles.
!--       u- and v-component.
          CALL pmci_interp_1d( u, u_p_init, kcto, kflo, kfuo )
          CALL pmci_interp_1d( v, v_p_init, kcto, kflo, kfuo )
!
!--       With homogeneous initialization set w = 0.0
          w = 0.0_wp
!
!--       SGS-TKE
          IF ( (    rans_mode_parent  .AND.         rans_mode )  .OR.                              &
            ( .NOT. rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                               &
              .NOT. constant_diffusion ) )  THEN
             CALL pmci_interp_1d( e, e_p_init, kcto, kflo, kfuo )
          ENDIF
!
!--       dissipation rate
          IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
             CALL pmci_interp_1d( diss, diss_p_init, kcto, kflo, kfuo )
          ENDIF
!
!--       potential temperature
          IF ( .NOT. neutral )  THEN
             CALL pmci_interp_1d( pt, pt_p_init, kcto, kflo, kfuo )
          ENDIF

          IF ( humidity )  THEN
!
!--          mixing ratio
             CALL pmci_interp_1d( q, q_p_init, kcto, kflo, kfuo )
!
!--          qc, nc
             IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
                 CALL pmci_interp_1d( qc, qc_p_init, kcto, kflo, kfuo )
                 CALL pmci_interp_1d( nc, nc_p_init, kcto, kflo, kfuo )
             ENDIF
!
!--          qr, nr
             IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
                 CALL pmci_interp_1d( qr, qr_p_init, kcto, kflo, kfuo )
                 CALL pmci_interp_1d( nr, nr_p_init, kcto, kflo, kfuo )
             ENDIF
          ENDIF
!
!--       passive scalar
          IF ( passive_scalar )  THEN
             CALL pmci_interp_1d( s, s_p_init, kcto, kflo, kfuo )
          ENDIF
!
!--       chemistry
          IF ( air_chemistry  .AND.  nesting_chem )  THEN
             DO  n = 1, nspec
                CALL pmci_interp_1d( chem_species(n)%conc, chem_p_init(:,n), kcto, kflo, kfuo )
             ENDDO
          ENDIF
!
!--       aerosols
          IF ( salsa  .AND.  nesting_salsa )  THEN
             DO  lb = 1, nbins_aerosol
                CALL pmci_interp_1d( aerosol_number(lb)%conc, aerosol_number_p_init(:,lb), kcto,   &
                                     kflo, kfuo )
             ENDDO
             DO  lc = 1, nbins_aerosol * ncomponents_mass
                CALL pmci_interp_1d( aerosol_mass(lc)%conc, aerosol_mass_p_init(:,lc), kcto, kflo, &
                                     kfuo )
             ENDDO
             IF ( .NOT. salsa_gases_from_chem )  THEN
                DO  lg = 1, ngases_salsa
                   CALL pmci_interp_1d( salsa_gas(lg)%conc, salsa_gas_p_init(:,lg), kcto, kflo,    &
                                        kfuo )
                ENDDO
             ENDIF
          ENDIF
!
!--    3d initialization from parent.
       ELSE
!
!--       Get 3d data from the parent.
          CALL pmc_c_getbuffer( waittime = waittime )
!
!--       Initialize child with 3d parent arrays.
          CALL pmci_interp_all( u, uc, kcto, iflu, ifuu, jflo, jfuo, kflo, kfuo, 'u' )
          CALL pmci_interp_all( v, vc, kcto, iflo, ifuo, jflv, jfuv, kflo, kfuo, 'v' )

          CALL pmci_interp_all( w, wc, kctw, iflo, ifuo, jflo, jfuo, kflw, kfuw, 'w' )
          CALL exchange_horiz( u, nbgp )
          CALL exchange_horiz( v, nbgp )
          CALL exchange_horiz( w, nbgp )

          IF ( (       rans_mode_parent  .AND.         rans_mode )  .OR.                           &
               ( .NOT. rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                            &
                 .NOT. constant_diffusion ) )  THEN
             CALL pmci_interp_all( e, ec, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 'e' )
             CALL exchange_horiz( e, nbgp )
          ENDIF

          IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
             CALL pmci_interp_all( diss, dissc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
             CALL exchange_horiz( diss, nbgp )
          ENDIF

          IF ( .NOT. neutral )  THEN
             CALL pmci_interp_all( pt, ptc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
             CALL exchange_horiz( pt, nbgp )
          ENDIF

          IF ( humidity )  THEN

             CALL pmci_interp_all( q, q_c, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
             CALL exchange_horiz( q, nbgp )

             IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
                CALL pmci_interp_all( qc, qcc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
                CALL pmci_interp_all( nc, ncc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
                CALL exchange_horiz( qc, nbgp )
                CALL exchange_horiz( nc, nbgp )
             ENDIF

             IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
                CALL pmci_interp_all( qr, qrc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
                CALL pmci_interp_all( nr, nrc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
                CALL exchange_horiz( qr, nbgp )
                CALL exchange_horiz( nr, nbgp )
             ENDIF

          ENDIF

          IF ( passive_scalar )  THEN
             CALL pmci_interp_all( s, sc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
             CALL exchange_horiz( s, nbgp )
          ENDIF

          IF ( ocean_mode  .AND.  salinity )  THEN
             CALL pmci_interp_all( sa, sac, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
             CALL exchange_horiz( sa, nbgp )
          ENDIF

          IF ( air_chemistry  .AND.  nesting_chem )  THEN
             DO  n = 1, nspec
                CALL pmci_interp_all ( chem_species(n)%conc, chem_spec_c(:,:,:,n), kcto, iflo,     &
                                       ifuo, jflo, jfuo, kflo, kfuo, 's' )
                CALL exchange_horiz( chem_species(n)%conc, nbgp )
             ENDDO
          ENDIF

          IF ( salsa  .AND.  nesting_salsa )  THEN
             DO  lb = 1, nbins_aerosol
                CALL pmci_interp_all ( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb), kcto,  &
                                       iflo, ifuo, jflo, jfuo, kflo, kfuo, 's' )
                CALL exchange_horiz( aerosol_number(lb)%conc, nbgp )
             ENDDO
             DO  lc = 1, nbins_aerosol * ncomponents_mass
                CALL pmci_interp_all ( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc), kcto, iflo,&
                                       ifuo, jflo, jfuo, kflo, kfuo, 's' )
                CALL exchange_horiz( aerosol_mass(lc)%conc, nbgp )
             ENDDO
             IF ( .NOT. salsa_gases_from_chem )  THEN
                DO  lg = 1, ngases_salsa
                   CALL pmci_interp_all ( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg), kcto, iflo,   &
                                          ifuo, jflo, jfuo, kflo, kfuo, 's' )
                   CALL exchange_horiz( salsa_gas(lg)%conc, nbgp )
                ENDDO
             ENDIF
          ENDIF
       ENDIF

       IF ( topography /= 'flat' .AND. .NOT. atmosphere_ocean_coupled_run )  THEN
!
!--       Inside buildings set velocities back to zero.
          DO  ic = nxlg, nxrg
             DO  jc = nysg, nyng
                DO  kc = nzb, nzt
                   u(kc,jc,ic)   = MERGE( u(kc,jc,ic), 0.0_wp, BTEST( topo_flags(kc,jc,ic), 1 ) )
                   v(kc,jc,ic)   = MERGE( v(kc,jc,ic), 0.0_wp, BTEST( topo_flags(kc,jc,ic), 2 ) )
                   w(kc,jc,ic)   = MERGE( w(kc,jc,ic), 0.0_wp, BTEST( topo_flags(kc,jc,ic), 3 ) )
                   u_p(kc,jc,ic) = MERGE( u_p(kc,jc,ic), 0.0_wp, BTEST( topo_flags(kc,jc,ic), 1 ) )
                   v_p(kc,jc,ic) = MERGE( v_p(kc,jc,ic), 0.0_wp, BTEST( topo_flags(kc,jc,ic), 2 ) )
                   w_p(kc,jc,ic) = MERGE( w_p(kc,jc,ic), 0.0_wp, BTEST( topo_flags(kc,jc,ic), 3 ) )
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF


 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation of the internal values for the child-domain initialization.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_interp_1d( child_array, parent_profile, kct, kfl, kfu )

    INTEGER(iwp) ::  kc               !< Running child-grid index in the z-direction
    INTEGER(iwp) ::  kp               !< Running parent-grid index in the z-direction
    INTEGER(iwp), INTENT(IN) ::  kct  !< The parent-grid index in z-direction just below the boundary value node

    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfl  !<  Indicates start index of child cells belonging to certain
                                                            !<  parent cell - z direction
    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfu  !<  Indicates end index of child cells belonging to certain
                                                            !<  parent cell - z direction
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) ::  child_array  !<  Child-grid 3D array
    REAL(wp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  parent_profile                      !<  Parent-grid 1D array



    child_array(:,:,:) = 0.0_wp
    DO  kp = 0, kct + 1
       DO  kc = kfl(kp), MIN( kfu(kp), nzt+1 )
          child_array(kc,:,:) = parent_profile(kp)
       ENDDO
    ENDDO

 END SUBROUTINE pmci_interp_1d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation of the internal values for the child-domain initialization.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_interp_all( child_array, parent_array, kct, ifl, ifu, jfl, jfu, kfl, kfu, var )

    IMPLICIT NONE

    CHARACTER(LEN=1), INTENT(IN) ::  var  !<  Variable symbol: 'u', 'v', 'w' or 's'

    INTEGER(iwp), INTENT(IN) ::  kct  !< The parent-grid index in z-direction just below the boundary value node

    INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifl  !<  Indicates start index of child cells belonging to certain
                                                            !<  parent cell - x direction
    INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifu  !<  Indicates end index of child cells belonging to certain
                                                            !<  parent cell - x direction
    INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfl  !<  Indicates start index of child cells belonging to certain
                                                            !<  parent cell - y direction
    INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfu  !<  Indicates end index of child cells belonging to certain
                                                            !<  parent cell - y direction
    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfl  !<  Indicates start index of child cells belonging to certain
                                                            !<  parent cell - z direction
    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfu  !<  Indicates end index of child cells belonging to certain
                                                            !<  parent cell - z direction
    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) ::  child_array  !<  Child-grid array

    REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(IN) ::  parent_array  !<  Parent-grid array

!
!-- Local variables:
    INTEGER(iwp) ::  ic        !< Running child-grid index in the x-direction
    INTEGER(iwp) ::  icb       !< Index pointing to the first redundant ghost point layer behind the actual boundary
                               !< ghost point layer in the x-direction
    INTEGER(iwp) ::  icbc      !< Index pointing to the boundary ghost point layer in the x-direction
    INTEGER(iwp) ::  icfirst   !< Leftmost child-grid index initialized by the main loops (usually icfirst == icl_init)
    INTEGER(iwp) ::  iclast    !< Rightmost child-grid index initialized by the main loops (usually iclast == icr_init)
    INTEGER(iwp) ::  icl_init  !< Left child-grid index bound for initialization in the x-direction
    INTEGER(iwp) ::  icr_init  !< Right child-grid index bound for initialization in the x-direction
    INTEGER(iwp) ::  ip        !< Running parent-grid index in the x-direction
    INTEGER(iwp) ::  ipl_init  !< Left parent-grid index bound for initialization in the x-direction
    INTEGER(iwp) ::  ipr_init  !< Right parent-grid index bound for initialization in the x-direction
    INTEGER(iwp) ::  jc        !< Running child-grid index in the y-direction
    INTEGER(iwp) ::  jcb       !< Index pointing to the first redundant ghost point layer behind the actual boundary
                               !< ghost point layer in the y-direction
    INTEGER(iwp) ::  jcbc      !< Index pointing to the boundary ghost point layer in the y-direction
    INTEGER(iwp) ::  jcfirst   !< Southmost child-grid index initialized by the main loops (usually jcfirst == jcs_init)
    INTEGER(iwp) ::  jclast    !< Northmost child-grid index initialized by the main loops (usually jclast == jcn_init)
    INTEGER(iwp) ::  jcs_init  !< South child-grid index bound for initialization in the y-direction
    INTEGER(iwp) ::  jcn_init  !< North child-grid index bound for initialization in the y-direction
    INTEGER(iwp) ::  jp        !< Running parent-grid index in the y-direction
    INTEGER(iwp) ::  jps_init  !< South parent-grid index bound for initialization in the y-direction
    INTEGER(iwp) ::  jpn_init  !< North parent-grid index bound for initialization in the y-direction
    INTEGER(iwp) ::  kc        !< Running child-grid index in the z-direction
    INTEGER(iwp) ::  kp        !< Running parent-grid index in the z-direction


    ipl_init = ipl
    ipr_init = ipr
    jps_init = jps
    jpn_init = jpn
    icl_init = nxl
    icr_init = nxr
    jcs_init = nys
    jcn_init = nyn

    icbc = -1
    icb  = -2
    jcbc = -1
    jcb  = -2
    IF ( var == 'u' )  THEN
       icbc =  0
       icb  = -1
    ELSEIF ( var == 'v' )  THEN
       jcbc =  0
       jcb  = -1
    ENDIF

    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       ipl_init = ipl + 1
       icl_init = nxl - 1
!
!--    For u, nxl is a ghost node, but not for the other variables
       IF ( var == 'u' )  THEN
          ipl_init = ipl + 2
          icl_init = nxl
       ENDIF
    ENDIF
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       jps_init = jps + 1
       jcs_init = nys - 1
!
!--    For v, nys is a ghost node, but not for the other variables
       IF ( var == 'v' )  THEN
          jps_init = jps + 2
          jcs_init = nys
       ENDIF
    ENDIF
    IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
       ipr_init = ipr - 1
       icr_init = nxr + 1
    ENDIF
    IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
       jpn_init = jpn - 1
       jcn_init = nyn + 1
    ENDIF

    child_array(:,:,:) = 0.0_wp

    IF ( var == 'u' )  THEN

       icfirst = ifl(ipl_init)
       iclast  = ifl(ipr_init+1) - 1
       jcfirst = jfl(jps_init)
       jclast  = jfu(jpn_init)
       DO  ip = ipl_init, ipr_init
          DO  jp = jps_init, jpn_init
             DO  kp = 0, kct + 1

                DO  ic = ifl(ip), ifl(ip+1)-1
                   DO  jc = jfl(jp), jfu(jp)
                      DO  kc = kfl(kp), MIN( kfu(kp), nzt+1 )
                         child_array(kc,jc,ic) = parent_array(kp,jp,ip)
                      ENDDO
                   ENDDO
                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( var == 'v' )  THEN

       icfirst = ifl(ipl_init)
       iclast  = ifu(ipr_init)
       jcfirst = jfl(jps_init)
       jclast  = jfl(jpn_init+1) - 1
       DO  ip = ipl_init, ipr_init
          DO  jp = jps_init, jpn_init
             DO  kp = 0, kct + 1

                DO  ic = ifl(ip), ifu(ip)
                   DO  jc = jfl(jp), jfl(jp+1)-1
                      DO  kc = kfl(kp), MIN( kfu(kp), nzt+1 )
                         child_array(kc,jc,ic) = parent_array(kp,jp,ip)
                      ENDDO
                   ENDDO
                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ELSEIF ( var == 'w' )  THEN

       icfirst = ifl(ipl_init)
       iclast  = ifu(ipr_init)
       jcfirst = jfl(jps_init)
       jclast  = jfu(jpn_init)
       DO  ip = ipl_init, ipr_init
          DO  jp = jps_init, jpn_init
             DO  kp = 1, kct + 1

                DO  ic = ifl(ip), ifu(ip)
                   DO  jc = jfl(jp), jfu(jp)
!
!--                   Because the kp-loop for w starts from kp=1 instead of 0
                      child_array(nzb,jc,ic) = 0.0_wp
                      DO  kc = kfu(kp-1)+1, kfu(kp)
                         child_array(kc,jc,ic) = parent_array(kp,jp,ip)
                      ENDDO
                   ENDDO
                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ELSE   ! Scalars

       icfirst = ifl(ipl_init)
       iclast  = ifu(ipr_init)
       jcfirst = jfl(jps_init)
       jclast  = jfu(jpn_init)
       DO  ip = ipl_init, ipr_init
          DO  jp = jps_init, jpn_init
             DO  kp = 0, kct + 1

                DO  ic = ifl(ip), ifu(ip)
                   DO  jc = jfl(jp), jfu(jp)
                      DO  kc = kfl(kp), MIN( kfu(kp), nzt+1 )
                         child_array(kc,jc,ic) = parent_array(kp,jp,ip)
                      ENDDO
                   ENDDO
                ENDDO

             ENDDO
          ENDDO
       ENDDO

    ENDIF  ! var
!
!-- If the number of grid points in child subdomain in x- or y-direction
!-- (nxr - nxl + 1 and/or nyn - nys + 1) is not integer divisible by the grid spacing ratio in its
!-- direction (igsr and/or jgsr), the above loops will return with unfilled gaps in the initial
!-- fields. These gaps, if present, are filled here.
    IF ( icfirst > icl_init )  THEN
       DO  ic = icl_init, icfirst - 1
          child_array(:,:,ic) = child_array(:,:,icfirst)
       ENDDO
    ENDIF
    IF ( iclast < icr_init )  THEN
       DO  ic = iclast + 1, icr_init
          child_array(:,:,ic) = child_array(:,:,iclast)
       ENDDO
    ENDIF
    IF ( jcfirst > jcs_init )  THEN
       DO  jc = jcs_init, jcfirst - 1
          child_array(:,jc,:) = child_array(:,jcfirst,:)
       ENDDO
    ENDIF
    IF ( jclast < jcn_init )  THEN
       DO  jc = jclast + 1, jcn_init
          child_array(:,jc,:) = child_array(:,jclast,:)
       ENDDO
    ENDIF
!
!-- Finally, make sure that also the redundant 2nd and 3rd ghost-node layers including the corners
!-- are properly filled up.
    IF ( nys == 0 )  THEN
       DO  jc = -nbgp, jcb  ! jcb = -2 if var == v, else jcb = -1
          child_array(0:nzt+1,jc,nxlg:nxrg) = child_array(0:nzt+1,jcbc,nxlg:nxrg)
       ENDDO
    ENDIF
    IF ( nyn == ny )  THEN
       DO  jc = ny+2, ny+nbgp
          child_array(0:nzt+1,jc,nxlg:nxrg) = child_array(0:nzt+1,ny+1,nxlg:nxrg)
       ENDDO
    ENDIF
    IF ( nxl == 0 )  THEN
       DO  ic = -nbgp, icb  ! icb = -2 if var == u, else icb = -1
          child_array(0:nzt+1,nysg:nyng,ic) = child_array(0:nzt+1,nysg:nyng,icbc)
       ENDDO
    ENDIF
    IF ( nxr == nx )  THEN
       DO  ic = nx+2, nx+nbgp
          child_array(0:nzt+1,nysg:nyng,ic) = child_array(0:nzt+1,nysg:nyng,nx+1)
       ENDDO
    ENDIF

 END SUBROUTINE pmci_interp_all

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Receive domain-averaged profiles of the parent domain and broadcast them among all child cores.
!> Attention: If the order of the receive operations is altered it needs to be altered in
!>            pmci_send_domain_averaged_profiles the same way!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_recv_domain_averaged_profiles

    INTEGER(iwp) ::  ierr   !< MPI error code
    INTEGER(iwp) ::  tag_nr !< tag number to identify send/receive operations
!
!-- u- and v-component
    ALLOCATE( u_p_init(0:pg%nz+1) )
    ALLOCATE( v_p_init(0:pg%nz+1) )
    IF ( myid == 0 )  THEN
       tag_nr = 50
       CALL pmc_recv_from_parent( u_p_init, pg%nz+2, 0, tag_nr, ierr )
       tag_nr = tag_nr + 1

       CALL pmc_recv_from_parent( v_p_init, pg%nz+2, 0, tag_nr, ierr )
       tag_nr = tag_nr + 1
    ENDIF
!
!-- SGS-TKE
    IF ( ( rans_mode_parent  .AND.         rans_mode )  .OR.                                    &
         ( .NOT. rans_mode_parent  .AND.  .NOT.  rans_mode  .AND.                               &
           .NOT. constant_diffusion ) )  THEN
       ALLOCATE( e_p_init(0:pg%nz+1) )
       IF ( myid == 0 )  THEN
          CALL pmc_recv_from_parent( e_p_init, pg%nz+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDIF
    ENDIF
!
!-- dissipation rate
    IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
       ALLOCATE( diss_p_init(0:pg%nz+1) )
       IF ( myid == 0 )  THEN
          CALL pmc_recv_from_parent( diss_p_init, pg%nz+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDIF
    ENDIF
!
!-- potential temperature
    IF ( .NOT. neutral )  THEN
       ALLOCATE( pt_p_init(0:pg%nz+1) )
       IF ( myid == 0 )  THEN
          CALL pmc_recv_from_parent( pt_p_init, pg%nz+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDIF
    ENDIF

    IF ( humidity )  THEN
!
!--    mixing ratio
       ALLOCATE( q_p_init(0:pg%nz+1) )
       IF ( myid == 0 )  THEN
          CALL pmc_recv_from_parent( q_p_init, pg%nz+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDIF
!
!--    qc and nc
       IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
          ALLOCATE( qc_p_init(0:pg%nz+1) )
          ALLOCATE( nc_p_init(0:pg%nz+1) )
          IF ( myid == 0 )  THEN
             CALL pmc_recv_from_parent( qc_p_init, pg%nz+2, 0, tag_nr, ierr )
             tag_nr = tag_nr + 1

             CALL pmc_recv_from_parent( nc_p_init, pg%nz+2, 0, tag_nr, ierr )
             tag_nr = tag_nr + 1
          ENDIF
       ENDIF
!
!--    qr and nr
       IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
          ALLOCATE( qr_p_init(0:pg%nz+1) )
          ALLOCATE( nr_p_init(0:pg%nz+1) )
          IF ( myid == 0 )  THEN
             CALL pmc_recv_from_parent( qr_p_init, pg%nz+2, 0, tag_nr, ierr )
             tag_nr = tag_nr + 1

             CALL pmc_recv_from_parent( nr_p_init, pg%nz+2, 0, tag_nr, ierr )
             tag_nr = tag_nr + 1
          ENDIF
       ENDIF
    ENDIF
!
!-- passive scalar
    IF ( passive_scalar )  THEN
       ALLOCATE( s_p_init(0:pg%nz+1) )
       IF ( myid == 0 )  THEN
          CALL pmc_recv_from_parent( s_p_init, pg%nz+2, 0, tag_nr, ierr )
          tag_nr = tag_nr + 1
       ENDIF
    ENDIF
!
!-- chemistry
    IF ( air_chemistry  .AND.  nesting_chem )  THEN
       ALLOCATE( chem_p_init(0:pg%nz+1,1:nspec) )
       IF ( myid == 0 )  THEN
          DO  n = 1, nspec
             CALL pmc_recv_from_parent( chem_p_init(:,n), pg%nz+2, 0, tag_nr, ierr )
             tag_nr = tag_nr + 1
          ENDDO
       ENDIF
    ENDIF
!
!-- aerosols
    IF ( salsa  .AND.  nesting_salsa )  THEN
       ALLOCATE( aerosol_number_p_init(0:pg%nz+1,1:nbins_aerosol)                  )
       ALLOCATE( aerosol_mass_p_init(0:pg%nz+1,1:nbins_aerosol * ncomponents_mass) )
       ALLOCATE( salsa_gas_p_init(0:pg%nz+1,1:ngases_salsa)                        )
       IF ( myid == 0 )  THEN
          DO  n = 1, nbins_aerosol
             CALL pmc_recv_from_parent( aerosol_number_p_init(:,n), pg%nz+2, 0, tag_nr, ierr )
             tag_nr = tag_nr + 1
          ENDDO
          DO  n = 1, nbins_aerosol * ncomponents_mass
             CALL pmc_recv_from_parent( aerosol_mass_p_init(:,n), pg%nz+2, 0, tag_nr, ierr )
             tag_nr = tag_nr + 1
          ENDDO
          DO  n = 1, ngases_salsa
             CALL pmc_recv_from_parent( salsa_gas_p_init(:,n), pg%nz+2, 0, tag_nr, ierr )
             tag_nr = tag_nr + 1
          ENDDO
       ENDIF
    ENDIF
!
!-- Broadcast the initial profiles among all child cores.
    IF ( ALLOCATED( u_p_init    ) )  CALL MPI_BCAST( u_p_init,    pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( v_p_init    ) )  CALL MPI_BCAST( v_p_init,    pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( e_p_init    ) )  CALL MPI_BCAST( e_p_init,    pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( diss_p_init ) )  CALL MPI_BCAST( diss_p_init, pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( pt_p_init   ) )  CALL MPI_BCAST( pt_p_init,   pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( q_p_init    ) )  CALL MPI_BCAST( q_p_init,    pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( qc_p_init   ) )  CALL MPI_BCAST( qc_p_init,   pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( nc_p_init   ) )  CALL MPI_BCAST( nc_p_init,   pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( qr_p_init   ) )  CALL MPI_BCAST( qr_p_init,   pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( nr_p_init   ) )  CALL MPI_BCAST( nr_p_init,   pg%nz+2, MPI_REAL, 0, comm2d, ierr )
    IF ( ALLOCATED( s_p_init    ) )  CALL MPI_BCAST( s_p_init,    pg%nz+2, MPI_REAL, 0, comm2d, ierr )

    IF ( ALLOCATED( chem_p_init ) )                                                                &
       CALL MPI_BCAST( chem_p_init, (pg%nz+2)*nspec, MPI_REAL, 0, comm2d, ierr )

    IF ( ALLOCATED( aerosol_number_p_init ) )                                                      &
       CALL MPI_BCAST( aerosol_number_p_init, (pg%nz+2)*nbins_aerosol, MPI_REAL, 0, comm2d, ierr )

    IF ( ALLOCATED( aerosol_mass_p_init ) )                                                        &
       CALL MPI_BCAST( aerosol_mass_p_init, (pg%nz+2)*nbins_aerosol*ncomponents_mass, MPI_REAL, 0, &
                       comm2d, ierr )

    IF ( ALLOCATED( salsa_gas_p_init ) )                                                           &
       CALL MPI_BCAST( salsa_gas_p_init, (pg%nz+2)*ngases_salsa, MPI_REAL, 0, comm2d, ierr )

 END SUBROUTINE pmci_recv_domain_averaged_profiles

#endif
 END SUBROUTINE pmci_child_initialize


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check for mismatches between settings of root and child variables (e.g., all children have to
!> follow the end_time settings of the root model). The root model overwrites variables in the
!> other models, so these variables only need to be set once in file PARIN.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_check_setting_mismatches

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER ::  ierr  !<  MPI error code

    LOGICAL ::  homogenize_surface_temperature_root  !< root value of variable
    LOGICAL ::  salinity_root                        !< root value of variable

    REAL(wp) ::  dt_coupling_root          !< root value of variable
    REAL(wp) ::  dt_restart_root           !< root value of variable
    REAL(wp) ::  dx_root                   !< root value of variable
    REAL(wp) ::  dy_root                   !< root value of variable
    REAL(wp) ::  end_time_root             !< root value of variable
    REAL(wp) ::  restart_time_root         !< root value of variable
    REAL(wp) ::  time_from_ref_point       !< root value of variable
    REAL(wp) ::  time_from_ref_point_root  !< root value of variable
    REAL(wp) ::  time_restart_root         !< root value of variable


!
!-- Check the time to be simulated. Here, and in the following, the root process communicates the
!-- respective variable to all others, and its value will then be compared with the local values.
    IF ( root_model )  end_time_root = end_time
    CALL MPI_BCAST( end_time_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. root_model )  THEN
       IF ( end_time /= end_time_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and child settings:& ' //       &
                                      'end_time(root) = ', end_time_root,                          &
                                      '& end_time(child) = ', end_time, '& child value is set',    &
                                      ' to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0021', 0, 1, 0, 6, 0 )
          end_time = end_time_root
       ENDIF
    ENDIF
!
!-- Same for restart time
    IF ( root_model )  restart_time_root = restart_time
    CALL MPI_BCAST( restart_time_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. root_model )  THEN
       IF ( restart_time /= restart_time_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and child settings: & ' //      &
                                      'restart_time(root) = ', restart_time_root,                  &
                                      '& restart_time(child) = ', restart_time, '& child ',        &
                                      'value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0021', 0, 1, 0, 6, 0 )
          restart_time = restart_time_root
       ENDIF
    ENDIF
!
!-- Same for dt_restart
    IF ( root_model )  dt_restart_root = dt_restart
    CALL MPI_BCAST( dt_restart_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. root_model )  THEN
       IF ( dt_restart /= dt_restart_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and ',                          &
                                      'child settings: & dt_restart(root) = ', dt_restart_root,    &
                                      '& dt_restart(child) = ', dt_restart, '& child ',            &
                                      'value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0021', 0, 1, 0, 6, 0 )
          dt_restart = dt_restart_root
       ENDIF
    ENDIF
!
!-- Same for time_restart
    IF ( root_model )  time_restart_root = time_restart
    CALL MPI_BCAST( time_restart_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. root_model )  THEN
       IF ( time_restart /= time_restart_root )  THEN
          WRITE( message_string, * )  'mismatch between root model and child settings: & ' //      &
                                      'time_restart(root) = ', time_restart_root,                  &
                                      '& time_restart(child) = ', time_restart, '& child ',        &
                                      'value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0021', 0, 1, 0, 6, 0 )
          time_restart = time_restart_root
       ENDIF
    ENDIF

!
!-- Homogenize_surface_temperature must be set in all domains.
    IF ( root_model )  homogenize_surface_temperature_root = homogenize_surface_temperature
    CALL MPI_BCAST( homogenize_surface_temperature_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. root_model )  THEN
       IF ( homogenize_surface_temperature .NEQV. homogenize_surface_temperature_root )  THEN
          message_string = 'mismatch between root model and child settings: & ' //                 &
                           'homogenize_surface_temperature must be set the same for all models' // &
                           '. &child value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0022', 0, 1, 0, 6, 0 )
          homogenize_surface_temperature = homogenize_surface_temperature_root
       ENDIF
    ENDIF

!
!-- In ocean mode salinity requires to be set in all models in the same way.
    IF ( root_model )  salinity_root = salinity
    CALL MPI_BCAST( salinity_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

    IF ( .NOT. root_model )  THEN
       IF ( salinity .NEQV. salinity_root )  THEN
          message_string = 'mismatch between root model and child settings: & ' //                 &
                           'salinity must be set the same for all models' //                       &
                           '. &child value is set to root value'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0021', 0, 1, 0, 6, 0 )
          salinity = salinity_root
       ENDIF
    ENDIF

!
!-- Further checks in case of atmosphere-ocean coupled runs.
    IF ( atmosphere_ocean_coupled_run )  THEN

       IF ( root_model  .AND.  dt_coupling == 9999999.9_wp )  THEN
          message_string = 'dt_coupling is required to be set in atmosphere runtime parameters' // &
                           ' namelist'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0023', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( root_model  .AND.  dt_coupling < 0.0_wp )  THEN
          message_string = 'dt_coupling < 0.0 is not allowed'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0024', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( root_model )  dt_coupling_root = dt_coupling
       CALL MPI_BCAST( dt_coupling_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

       IF ( dt_coupling /= dt_coupling_root )  THEN
          WRITE( message_string, * )  'mismatch between atmosphere and ocean model: & ' //         &
                                      'dt_coupling (atmosphere) = ', dt_coupling_root,             &
                                      '& dt_coupling (ocean) = ', dt_coupling,                     &
                                      '& ocean value is set to atmosphere value'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0025', 0, 1, 0, 6, 0 )
          dt_coupling = dt_coupling_root
       ENDIF

       time_from_ref_point = end_time - coupling_start_time
       IF ( root_model )  time_from_ref_point_root = time_from_ref_point
       CALL MPI_BCAST( time_from_ref_point_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

       IF ( time_from_ref_point /= time_from_ref_point_root )  THEN
          WRITE( message_string, * )  'mismatch between atmosphere and ocean model: & ' //         &
                            'time from reference point (atmosphere) = ', time_from_ref_point_root, &
                            '& time_from reference point (ocean) = ', time_from_ref_point
          CALL message( 'pmci_check_setting_mismatches', 'PMC0026', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( root_model )  dx_root = dx
       CALL MPI_BCAST( dx_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

       IF ( dx_root < dx )  THEN
          message_string = 'dx in atmosphere is smaller than dx in ocean'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0027', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( root_model )  dy_root = dy
       CALL MPI_BCAST( dy_root, 1, MPI_REAL, 0, comm_world_nesting, ierr )

       IF ( dy_root < dy )  THEN
          message_string = 'dy in atmosphere is smaller than dy in ocean'
          CALL message( 'pmci_check_setting_mismatches', 'PMC0027', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF
#endif

 END SUBROUTINE pmci_check_setting_mismatches


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Unify the time steps for each model and synchronize using MPI_ALLREDUCE with the MPI_MIN
!> operator over all processes using the global communicator MPI_COMM_WORLD.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_synchronize

#if defined( __parallel )

   IMPLICIT NONE

   INTEGER(iwp) ::  ierr  !< MPI error code

   REAL(wp) ::  dtl  !< Local time step of the current process
   REAL(wp) ::  dtg  !< Global time step defined as the global minimum of dtl of all processes


   IF ( debug_output_timestep )  CALL debug_message( 'pmci_synchronize', 'start' )

   dtl = dt_3d
   CALL MPI_ALLREDUCE( dtl, dtg, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr )
   dt_3d  = dtg

   IF ( debug_output_timestep )  CALL debug_message( 'pmci_synchronize', 'end' )

#endif
 END SUBROUTINE pmci_synchronize


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> After each Runge-Kutta sub-timestep, alternately set buffer one or buffer two active
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_set_swaplevel( swaplevel )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  swaplevel  !< swaplevel (1 or 2) of PALM's timestep

    INTEGER(iwp) ::  child_id  !<  Child id of the child number m
    INTEGER(iwp) ::  m         !<  Loop index over all children of the current parent

#if defined( __parallel )
    DO  m = 1, SIZE( pmc_parent_for_child ) - 1
       child_id = pmc_parent_for_child(m)
       CALL pmc_s_set_active_data_array( child_id, swaplevel )
    ENDDO
#endif
 END SUBROUTINE pmci_set_swaplevel


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine ...
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_atmos_ocean( )

    IMPLICIT NONE

#if defined( __parallel )
    CALL cpu_log( log_point(39), 'atmos_ocean_coupler', 'start' )

!
!-- Copy 1d-surface data to 2d-xy-arrays (sensible/latent heat- and momentum-fluxes).
    IF ( root_model )  CALL surface_coupler_buffer_handling( parent_send )
!
!-- Transfer data from atmosphere to ocean.
    CALL pmci_child_datatrans( parent_to_child )
    CALL pmci_parent_datatrans( parent_to_child )
!
!-- Interpolation from atmosphere to ocean grid (and backwards) is done in the ocean model.
    IF ( .NOT.  root_model )  THEN
!
!--    Interpolate the fluxes.
       CALL pmc_interpolate_to_ocean( surface_coupler_exchange_array_1c,                           &
                                      surface_coupler_exchange_array_1 )
       CALL pmc_interpolate_to_ocean( surface_coupler_exchange_array_2c,                           &
                                      surface_coupler_exchange_array_2 )
       CALL pmc_interpolate_to_ocean( surface_coupler_exchange_array_3c,                           &
                                      surface_coupler_exchange_array_3 )
       CALL pmc_interpolate_to_ocean( surface_coupler_exchange_array_4c,                           &
                                      surface_coupler_exchange_array_4 )
!
!--    Copy interpolated 2d-xy-arrays (sensible/latent heat- and momentum-fluxes) to 1d-surface-data
       CALL surface_coupler_buffer_handling( child_recv )
!
!--    Copy pt, u, v at ocean surface to 2d-xy-arrays.
       CALL surface_coupler_buffer_handling( child_send )
!
!--    Anterpolate pt, u, and v
       CALL pmc_interpolate_to_atmos( surface_coupler_exchange_array_1,                            &
                                      surface_coupler_exchange_array_1c )
       CALL pmc_interpolate_to_atmos( surface_coupler_exchange_array_2,                            &
                                      surface_coupler_exchange_array_2c )
       CALL pmc_interpolate_to_atmos( surface_coupler_exchange_array_3,                            &
                                      surface_coupler_exchange_array_3c )
    ENDIF

    CALL pmci_parent_datatrans( child_to_parent )
    CALL pmci_child_datatrans( child_to_parent )
    IF ( root_model )  CALL surface_coupler_buffer_handling( parent_recv )

    CALL cpu_log( log_point(39), 'atmos_ocean_coupler', 'stop' )

 CONTAINS

 SUBROUTINE pmc_interpolate_to_ocean( surface_coupler_exchange_array_parent,                       &
                                      surface_coupler_exchange_array )

    IMPLICIT NONE

    INTEGER(iwp) ::  dx_nr  !<
    INTEGER(iwp) ::  dy_nr  !<
    INTEGER(iwp) ::  i      !<
    INTEGER(iwp) ::  ii     !<
    INTEGER(iwp) ::  ix     !<
    INTEGER(iwp) ::  j      !<
    INTEGER(iwp) ::  jj     !<
    INTEGER(iwp) ::  jy     !<

    REAL(wp), PARAMETER ::  tolerance = 1.0E-10  !<

    REAL(wp) :: dx_fact  !<
    REAL(wp) :: dy_fact  !<
    REAL(wp) :: fyx      !<
    REAL(wp) :: f00      !<
    REAL(wp) :: f01      !<
    REAL(wp) :: f10      !<
    REAL(wp) :: f11      !<
    REAL(wp) :: x        !<
    REAL(wp) :: y        !<

    REAL(wp), INTENT(IN), DIMENSION(jps:jpn,ipl:ipr)      ::  surface_coupler_exchange_array_parent  !<
    REAL(wp), INTENT(OUT), DIMENSION(nysg:nyng,nxlg:nxrg) ::  surface_coupler_exchange_array         !<


    surface_coupler_exchange_array = 0.0_wp
    dx_fact = pg%dx / dx
    dy_fact = pg%dy / dy
    dx_nr   = INT( dx_fact + 0.01_wp )
    dy_nr   = INT( dy_fact + 0.01_wp )

    IF ( ABS( pg%dx - dx ) <= tolerance  .AND.  ABS( pg%dy - dy ) <= tolerance )  THEN
!
!--    No interpolation.
       surface_coupler_exchange_array(jps:jpn,ipl:ipr) =                                           &
                                             surface_coupler_exchange_array_parent(jps:jpn,ipl:ipr)
    ELSE
!
!--    Interpolation from atmosphere-grid to ocean-grid.
       DO  i = ipl, ipr-1
          DO  j = jps, jpn-1
              f00 = surface_coupler_exchange_array_parent(j,i)
              f01 = surface_coupler_exchange_array_parent(j,i+1)
              f10 = surface_coupler_exchange_array_parent(j+1,i)
              f11 = surface_coupler_exchange_array_parent(j+1,i+1)
              DO  ix =0, dx_nr-1
                 ii  = i * dx_nr + ix
                 DO  jy = 0, dy_nr-1
                    x = ix / dx_fact
                    y = jy / dy_fact
                    fyx = f00 * (1-y) * (1-x) + f10 * y * (1-x) + f01 * (1-y) * x + f11 * y * x
                    jj = j * dy_nr + jy
                    surface_coupler_exchange_array(jj,ii) = fyx
                 ENDDO
             ENDDO
          ENDDO
       ENDDO

    ENDIF

 END SUBROUTINE pmc_interpolate_to_ocean


 SUBROUTINE pmc_interpolate_to_atmos( surface_coupler_exchange_array,                              &
                                      surface_coupler_exchange_array_parent )

    IMPLICIT NONE

    INTEGER(iwp) :: dx_nr  !<
    INTEGER(iwp) :: dy_nr  !<
    INTEGER(iwp) :: i      !<
    INTEGER(iwp) :: ii     !<
    INTEGER(iwp) :: j      !<
    INTEGER(iwp) :: jj     !<

    REAL(wp), PARAMETER :: tolerance = 1.0E-10  !<

    REAL(wp) ::  dx_fact  !<
    REAL(wp) ::  dy_fact  !<

    REAL(wp), INTENT(IN), DIMENSION(nysg:nyng,nxlg:nxrg) ::  surface_coupler_exchange_array         !<
    REAL(wp), INTENT(OUT), DIMENSION(jps:jpn,ipl:ipr)    ::  surface_coupler_exchange_array_parent  !<


    dx_fact = pg%dx / dx
    dy_fact = pg%dy / dy
    dx_nr   = INT( dx_fact + 0.01_wp )
    dy_nr   = INT( dy_fact + 0.01_wp )

    surface_coupler_exchange_array_parent = 0.0_wp

    IF ( ABS( pg%dx - dx ) <= tolerance  .AND.  ABS( pg%dy - dy ) <= tolerance )  THEN
!
!--    No Interpolation.
       surface_coupler_exchange_array_parent(jps:jpn,ipl:ipr) =                                    &
                                                    surface_coupler_exchange_array(jps:jpn,ipl:ipr)
    ELSE
!
!--    Interpolation from ocean-grid to atmosphere-grid.
       DO  i = ipl, ipr
          DO  j = jps, jpn
             DO  ii = i*dx_nr, (i+1)*dx_nr-1
                DO  jj = j*dy_nr, (j+1)*dy_nr-1
                   surface_coupler_exchange_array_parent(j,i) =                                    &
                                                      surface_coupler_exchange_array_parent(j,i) + &
                                                      surface_coupler_exchange_array(jj,ii)
                ENDDO
             ENDDO
             surface_coupler_exchange_array_parent(j,i) =                                          &
                                                      surface_coupler_exchange_array_parent(j,i) / &
                                                      ( dx_fact * dy_fact )
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE pmc_interpolate_to_atmos
#endif

 END SUBROUTINE pmci_atmos_ocean


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust the time interval after which atmosphere and ocean are coupled, to be at least
!> as large as the maximum timestep of the atmosphere and the ocean model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_adjust_dt_coupling

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  ierr  !< MPI error code

    REAL(wp), SAVE ::  dt_coupling_org = -1.0_wp  !< original value of dt_coupling as specified in namelist

    REAL(wp), DIMENSION(2) ::  maxvalue  !< local variable to store global maximum values

!
!-- At first call, save the coupling interval that has been set by the user. If the maximum of the
!-- atmosphere/ocean timestep is smaller, then this interval will be used. Otherwise, the coupling
!-- interval will be adjusted (decreased).
    IF ( dt_coupling_org < 0.0_wp )  dt_coupling_org = dt_coupling

    maxvalue(1) = dt_3d
!
!-- Set flag to terminate the run, if either the atmosphere or the ocean model has reached its end.
    IF ( time_since_reference_point >= end_time )  THEN
       maxvalue(2) = 1.0_wp
    ELSE
       maxvalue(2) = 0.0_wp
    ENDIF
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, maxvalue, 2, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )
    terminate_coupled = maxvalue(2)

!
!-- The new coupling interval is taken as the maximum of the original interval (as set by the user)
!-- and the maximum timestep of the atmosphere/ocean model.
    dt_coupling = MAX( dt_coupling_org, maxvalue(1) )
#endif

 END SUBROUTINE pmci_adjust_dt_coupling


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine controls the nesting according to the nestpar parameter nesting_mode (two-way
!> (default) or one-way) and the order of anterpolations according to the nestpar parameter
!> nesting_datatransfer_mode (cascade, overlap or mixed (default)). Although nesting_mode is a
!> variable of this module, pass it as an argument to allow for example to force one-way
!> initialization phase. Note that interpolation ( parent_to_child ) must always be carried out
!> before anterpolation ( child_to_parent ).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_datatrans( local_nesting_mode )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  local_nesting_mode  !<  Nesting mode: 'one-way', 'two-way' or 'vertical'

#if defined( __parallel )

    IF ( debug_output_timestep )  CALL debug_message( 'pmci_datatrans', 'start' )

    IF ( TRIM( local_nesting_mode ) == 'one-way' )  THEN

       CALL pmci_child_datatrans( parent_to_child )
       CALL pmci_parent_datatrans( parent_to_child )

    ELSE

       IF ( nesting_datatransfer_mode == 'cascade' )  THEN

          CALL pmci_child_datatrans( parent_to_child )
          CALL pmci_parent_datatrans( parent_to_child )

          CALL pmci_parent_datatrans( child_to_parent )
          CALL pmci_child_datatrans( child_to_parent )

       ELSEIF ( nesting_datatransfer_mode == 'overlap')  THEN

          CALL pmci_parent_datatrans( parent_to_child )
          CALL pmci_child_datatrans( parent_to_child )

          CALL pmci_child_datatrans( child_to_parent )
          CALL pmci_parent_datatrans( child_to_parent )

       ELSEIF ( TRIM( nesting_datatransfer_mode ) == 'mixed' )  THEN

          CALL pmci_parent_datatrans( parent_to_child )
          CALL pmci_child_datatrans( parent_to_child )

          CALL pmci_parent_datatrans( child_to_parent )
          CALL pmci_child_datatrans( child_to_parent )

       ENDIF

    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'pmci_datatrans', 'end' )

#endif
 END SUBROUTINE pmci_datatrans


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Last action in pmc: currently only freeing of MPI windows
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_finalize

#if defined( __parallel )
    INTEGER(iwp) ::  child_id
    INTEGER(iwp) ::  m

!
!-- Free windows for transfer to/from parent
    IF ( .NOT. root_model )  THEN
       CALL pmc_c_finalize()
    ENDIF
!
!-- Free windows for transfer to/from child(s)
    DO  m = 1, SIZE( pmc_parent_for_child ) - 1
       child_id = pmc_parent_for_child(m)
       CALL pmc_s_finalize( child_id )
    ENDDO
#endif

 END SUBROUTINE pmci_finalize

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Wrapper routine including the parent-side calls to all the data-transfer routines.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_parent_datatrans( direction )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  direction  !<  Direction of the data transfer: 'parent_to_child' or 'child_to_parent'

#if defined( __parallel )
    INTEGER(iwp) ::  child_id  !<  Child id of the child number m
    INTEGER(iwp) ::  i         !<  Parent-grid index in x-direction
    INTEGER(iwp) ::  j         !<  Parent-grid index in y-direction
    INTEGER(iwp) ::  k         !<  Parent-grid index in z-direction
    INTEGER(iwp) ::  m         !<  Loop index over all children of the current parent


    DO  m = 1, SIZE( pmc_parent_for_child ) - 1
       child_id = pmc_parent_for_child(m)
       IF ( direction == parent_to_child )  THEN
          CALL cpu_log( log_point_s(71), 'pmc parent send', 'start' )
          CALL pmc_s_fillbuffer( child_id )
          CALL cpu_log( log_point_s(71), 'pmc parent send', 'stop' )
       ELSE
!
!--       Communication from child to parent
          CALL cpu_log( log_point_s(72), 'pmc parent recv', 'start' )
          child_id = pmc_parent_for_child(m)
          CALL pmc_s_getdata_from_buffer( child_id )
          CALL cpu_log( log_point_s(72), 'pmc parent recv', 'stop' )
!
!--       The anterpolated data is now available in u etc
          IF ( topography /= 'flat' .AND. .NOT. atmosphere_ocean_coupled_run )  THEN
!
!--          Inside buildings/topography reset velocities back to zero.
!--          Scalars (pt, q, s, km, kh, p, sa, ...) are ignored at present, maybe revise later.
!--          Resetting of e is removed as unnecessary since e is not interpolated, and as incorrect
!--          since it overran the default Neumann condition (bc_e_b).
             DO   i = nxlg, nxrg
                DO   j = nysg, nyng
                   DO  k = nzb, nzt+1
                      u(k,j,i) = MERGE( u(k,j,i), 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
                      v(k,j,i) = MERGE( v(k,j,i), 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
                      w(k,j,i) = MERGE( w(k,j,i), 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
!
!--                 TO_DO: zero setting of temperature within topography creates wrong results
!                   pt(nzb:nzb_s_inner(j,i),j,i) = 0.0_wp
!                   IF ( humidity  .OR.  passive_scalar )  THEN
!                      q(nzb:nzb_s_inner(j,i),j,i) = 0.0_wp
!                   ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDDO  ! m

#endif
 END SUBROUTINE pmci_parent_datatrans


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Wrapper routine including the child-side calls to all the data-transfer and processing routines.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_child_datatrans( direction )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  direction  !< transfer direction: parent_to_child or child_to_parent

#if defined( __parallel )

    REAL(wp), DIMENSION(1) ::  dtl  !< time step size


    dtl = dt_3d
    IF ( .NOT.  root_model )  THEN

       IF ( direction == parent_to_child )  THEN

          CALL cpu_log( log_point_s(73), 'pmc child recv', 'start' )
          CALL pmc_c_getbuffer( )
          CALL cpu_log( log_point_s(73), 'pmc child recv', 'stop' )

          IF ( .NOT. atmosphere_ocean_coupled_run )  THEN
             CALL cpu_log( log_point_s(75), 'pmc interpolation', 'start' )
             CALL pmci_interpolation
             CALL cpu_log( log_point_s(75), 'pmc interpolation', 'stop' )
          ENDIF

       ELSE
!
!--       direction == child_to_parent
          IF ( .NOT. atmosphere_ocean_coupled_run )  THEN
             CALL cpu_log( log_point_s(76), 'pmc anterpolation', 'start' )
             CALL pmci_anterpolation
             CALL cpu_log( log_point_s(76), 'pmc anterpolation', 'stop' )
          ENDIF

          CALL cpu_log( log_point_s(74), 'pmc child send', 'start' )
          CALL pmc_c_putbuffer( )
          CALL cpu_log( log_point_s(74), 'pmc child send', 'stop' )

       ENDIF
    ENDIF

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> A wrapper routine for all interpolation actions.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_interpolation

    IMPLICIT NONE

    INTEGER(iwp) ::  ibgp  !< Index running over the nbgp boundary ghost points in i-direction
    INTEGER(iwp) ::  jbgp  !< Index running over the nbgp boundary ghost points in j-direction
    INTEGER(iwp) ::  lb    !< Running index for aerosol size bins
    INTEGER(iwp) ::  lc    !< Running index for aerosol mass bins
    INTEGER(iwp) ::  lg    !< Running index for salsa gases
    INTEGER(iwp) ::  n     !< Running index for number of chemical species

!
!-- Interpolation is only needed for the lateral child boundaries, if they are non-cyclic.
!-- Start with boundary on the left-border PE:
    IF ( bc_dirichlet_l  .AND.  .NOT. nesting_bounds_vertical_only )  THEN

       CALL pmci_interp_lr( u, uc, kcto, jflo, jfuo, kflo, kfuo, 'l', 'u' )
       CALL pmci_interp_lr( v, vc, kcto, jflv, jfuv, kflo, kfuo, 'l', 'v' )
       CALL pmci_interp_lr( w, wc, kctw, jflo, jfuo, kflw, kfuw, 'l', 'w' )
!
!--    Treatment of TKE. Interpolation is only required if parent and child operate in RANS mode,
!--    else, interpolation is replaced by a Neumann condition.
       IF ( rans_mode_parent  .AND.  rans_mode )  THEN
           CALL pmci_interp_lr( e, ec, kcto, jflo, jfuo, kflo, kfuo, 'l', 'e' )
       ELSE
          DO  ibgp = -nbgp, -1
             e(nzb:nzt,nys:nyn,ibgp) = e(nzb:nzt,nys:nyn,0)
          ENDDO
       ENDIF

       IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
          CALL pmci_interp_lr( diss, dissc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
       ENDIF

       IF ( .NOT. neutral )  THEN
          CALL pmci_interp_lr( pt, ptc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
       ENDIF

       IF ( humidity )  THEN

          CALL pmci_interp_lr( q, q_c, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )

          IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
             CALL pmci_interp_lr( qc, qcc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
             CALL pmci_interp_lr( nc, ncc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
          ENDIF

          IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
             CALL pmci_interp_lr( qr, qrc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
             CALL pmci_interp_lr( nr, nrc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
          ENDIF

       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmci_interp_lr( s, sc, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
       ENDIF

       IF ( ocean_mode  .AND.  salinity )  THEN
          CALL pmci_interp_lr( sa, sac, kcto, jflo, jfuo, kflo, kfuo, 'l', 's' )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO  n = 1, nspec
             CALL pmci_interp_lr( chem_species(n)%conc, chem_spec_c(:,:,:,n), kcto, jflo, jfuo,    &
                                  kflo, kfuo, 'l', 's' )
          ENDDO
       ENDIF

       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             CALL pmci_interp_lr( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb), kcto, jflo, &
                                  jfuo, kflo, kfuo, 'l', 's')
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             CALL pmci_interp_lr( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc), kcto, jflo,     &
                                  jfuo, kflo, kfuo, 'l', 's')
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                CALL pmci_interp_lr( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg), kcto, jflo, jfuo,  &
                                     kflo, kfuo, 'l', 's')
             ENDDO
          ENDIF
       ENDIF

    ENDIF
!
!-- Right-border PE:
    IF ( bc_dirichlet_r  .AND.  .NOT. nesting_bounds_vertical_only )  THEN

       CALL pmci_interp_lr( u, uc, kcto, jflo, jfuo, kflo, kfuo, 'r', 'u' )
       CALL pmci_interp_lr( v, vc, kcto, jflv, jfuv, kflo, kfuo, 'r', 'v' )
       CALL pmci_interp_lr( w, wc, kctw, jflo, jfuo, kflw, kfuw, 'r', 'w' )
!
!--    Treatment of TKE. Interpolation is only required if parent and child operate in RANS mode,
!--    else, interpolation is replaced by a Neumann condition.
       IF ( rans_mode_parent  .AND.  rans_mode )  THEN
          CALL pmci_interp_lr( e, ec, kcto, jflo, jfuo, kflo, kfuo, 'r', 'e' )
       ELSE
          DO  ibgp = nx+1, nx+nbgp
             e(nzb:nzt,nys:nyn,ibgp) = e(nzb:nzt,nys:nyn,nx)
          ENDDO
       ENDIF

       IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
          CALL pmci_interp_lr( diss, dissc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
       ENDIF

       IF (  .NOT.  neutral )  THEN
          CALL pmci_interp_lr( pt, ptc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
       ENDIF

       IF ( humidity )  THEN

          CALL pmci_interp_lr( q, q_c, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )

          IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
             CALL pmci_interp_lr( qc, qcc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
             CALL pmci_interp_lr( nc, ncc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
          ENDIF

          IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
             CALL pmci_interp_lr( qr, qrc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
             CALL pmci_interp_lr( nr, nrc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
          ENDIF

       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmci_interp_lr( s, sc, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
       ENDIF

       IF ( ocean_mode  .AND.  salinity )  THEN
          CALL pmci_interp_lr( sa, sac, kcto, jflo, jfuo, kflo, kfuo, 'r', 's' )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO  n = 1, nspec
             CALL pmci_interp_lr( chem_species(n)%conc, chem_spec_c(:,:,:,n), kcto, jflo, jfuo,    &
                                  kflo, kfuo, 'r', 's' )
          ENDDO
       ENDIF

       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             CALL pmci_interp_lr( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb), kcto, jflo, &
                                  jfuo, kflo, kfuo, 'r', 's' )
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             CALL pmci_interp_lr( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc), kcto, jflo,     &
                                  jfuo, kflo, kfuo, 'r', 's' )
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                CALL pmci_interp_lr( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg), kcto, jflo, jfuo,  &
                                     kflo, kfuo, 'r', 's' )
             ENDDO
          ENDIF
       ENDIF

    ENDIF
!
!-- South-border PE:
    IF ( bc_dirichlet_s  .AND.  .NOT. nesting_bounds_vertical_only )  THEN

       CALL pmci_interp_sn( v, vc, kcto, iflo, ifuo, kflo, kfuo, 's', 'v' )
       CALL pmci_interp_sn( w, wc, kctw, iflo, ifuo, kflw, kfuw, 's', 'w' )
       CALL pmci_interp_sn( u, uc, kcto, iflu, ifuu, kflo, kfuo, 's', 'u' )
!
!--    Treatment of TKE. Interpolation is only required if parent and child operate in RANS mode,
!--    else, interpolation is replaced by a Neumann condition.
       IF ( rans_mode_parent  .AND.  rans_mode )  THEN
          CALL pmci_interp_sn( e, ec, kcto, iflo, ifuo, kflo, kfuo, 's', 'e' )
       ELSE
          DO  jbgp = -nbgp, -1
             e(nzb:nzt,jbgp,nxl:nxr) = e(nzb:nzt,0,nxl:nxr)
          ENDDO
       ENDIF

       IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
          CALL pmci_interp_sn( diss, dissc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
       ENDIF

       IF (  .NOT.  neutral )  THEN
          CALL pmci_interp_sn( pt, ptc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
       ENDIF

       IF ( humidity )  THEN

          CALL pmci_interp_sn( q, q_c, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )

          IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
             CALL pmci_interp_sn( qc, qcc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
             CALL pmci_interp_sn( nc, ncc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
          ENDIF

          IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
             CALL pmci_interp_sn( qr, qrc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
             CALL pmci_interp_sn( nr, nrc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
          ENDIF

       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmci_interp_sn( s, sc, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
       ENDIF

       IF ( ocean_mode  .AND.  salinity )  THEN
          CALL pmci_interp_sn( sa, sac, kcto, iflo, ifuo, kflo, kfuo, 's', 's' )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO  n = 1, nspec
             CALL pmci_interp_sn( chem_species(n)%conc, chem_spec_c(:,:,:,n), kcto, iflo, ifuo,    &
                                  kflo, kfuo, 's', 's' )
          ENDDO
       ENDIF

       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             CALL pmci_interp_sn( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb), kcto, iflo, &
                                  ifuo, kflo, kfuo, 's', 's' )
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             CALL pmci_interp_sn( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc), kcto, iflo,     &
                                  ifuo, kflo, kfuo, 's', 's' )
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                CALL pmci_interp_sn( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg), kcto, iflo, ifuo,  &
                                     kflo, kfuo, 's', 's' )
             ENDDO
          ENDIF
       ENDIF

    ENDIF
!
!-- North-border PE:
    IF ( bc_dirichlet_n  .AND.  .NOT. nesting_bounds_vertical_only )  THEN

       CALL pmci_interp_sn( v, vc, kcto, iflo, ifuo, kflo, kfuo, 'n', 'v' )
       CALL pmci_interp_sn( w, wc, kctw, iflo, ifuo, kflw, kfuw, 'n', 'w' )
       CALL pmci_interp_sn( u, uc, kcto, iflu, ifuu, kflo, kfuo, 'n', 'u' )
!
!--    Treatment of TKE. Interpolation is only required if parent and child operate in RANS mode,
!--    else, interpolation is replaced by a Neumann condition.
       IF ( rans_mode_parent  .AND.  rans_mode )  THEN
          CALL pmci_interp_sn( e, ec, kcto, iflo, ifuo, kflo, kfuo, 'n', 'e' )
       ELSE
          DO  jbgp = ny+1, ny+nbgp
             e(nzb:nzt,jbgp,nxl:nxr) = e(nzb:nzt,ny,nxl:nxr)
          ENDDO
       ENDIF

       IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
          CALL pmci_interp_sn( diss, dissc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
       ENDIF

       IF (  .NOT.  neutral )  THEN
          CALL pmci_interp_sn( pt, ptc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
       ENDIF

       IF ( humidity )  THEN

          CALL pmci_interp_sn( q, q_c, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )

          IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
             CALL pmci_interp_sn( qc, qcc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
             CALL pmci_interp_sn( nc, ncc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
          ENDIF

          IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
             CALL pmci_interp_sn( qr, qrc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
             CALL pmci_interp_sn( nr, nrc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
          ENDIF

       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmci_interp_sn( s, sc, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
       ENDIF

       IF ( ocean_mode  .AND.  salinity )  THEN
          CALL pmci_interp_sn( sa, sac, kcto, iflo, ifuo, kflo, kfuo, 'n', 's' )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO  n = 1, nspec
             CALL pmci_interp_sn( chem_species(n)%conc, chem_spec_c(:,:,:,n), kcto, iflo, ifuo,    &
                                  kflo, kfuo, 'n', 's' )
          ENDDO
       ENDIF

       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             CALL pmci_interp_sn( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb), kcto, iflo, &
                                  ifuo, kflo, kfuo, 'n', 's' )
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             CALL pmci_interp_sn( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc), kcto, iflo,     &
                                  ifuo, kflo, kfuo, 'n', 's' )
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                CALL pmci_interp_sn( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg), kcto, iflo, ifuo,  &
                                     kflo, kfuo, 'n', 's' )
             ENDDO
          ENDIF
       ENDIF

    ENDIF

!
!-- Top-border PEs:
    IF ( nested_bc_at_top )  THEN

       CALL pmci_interp_bt( w, wc, kctw, iflo, ifuo, jflo, jfuo, 't', 'w' )
       CALL pmci_interp_bt( u, uc, kcto, iflu, ifuu, jflo, jfuo, 't', 'u' )
       CALL pmci_interp_bt( v, vc, kcto, iflo, ifuo, jflv, jfuv, 't', 'v' )
!
!--    Treatment of TKE. Interpolation is only required if parent and child operate in RANS mode,
!--    else, interpolation is replaced by a Neumann condition.
       IF ( rans_mode_parent  .AND.  rans_mode )  THEN
          CALL pmci_interp_bt( e, ec, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
       ELSE
          e(nzt+1,nys:nyn,nxl:nxr) = e(nzt,nys:nyn,nxl:nxr)
       ENDIF

       IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
          CALL pmci_interp_bt( diss, dissc, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
       ENDIF

       IF ( .NOT. neutral )  THEN
          CALL pmci_interp_bt( pt, ptc, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
       ENDIF

       IF ( humidity )  THEN
          CALL pmci_interp_bt( q, q_c, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
          IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
             CALL pmci_interp_bt( qc, qcc, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
             CALL pmci_interp_bt( nc, ncc, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
          ENDIF
          IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
             CALL pmci_interp_bt( qr, qrc, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
             CALL pmci_interp_bt( nr, nrc, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
          ENDIF
       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmci_interp_bt( s, sc, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
       ENDIF

       IF ( ocean_mode  .AND.  salinity )  THEN
          CALL pmci_interp_bt( sa, sac, kcto, iflo, ifuo, jflo, jfuo, 't', 's' )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO  n = 1, nspec
             CALL pmci_interp_bt( chem_species(n)%conc, chem_spec_c(:,:,:,n), kcto, iflo, ifuo,    &
                                  jflo, jfuo, 't', 's' )
          ENDDO
       ENDIF

       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             CALL pmci_interp_bt( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb), kcto, iflo, &
                                  ifuo, jflo, jfuo, 't', 's' )
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             CALL pmci_interp_bt( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc), kcto, iflo,     &
                                  ifuo, jflo, jfuo, 't', 's' )
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                CALL pmci_interp_bt( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg), kcto, iflo, ifuo,  &
                                     jflo, jfuo, 't', 's' )
             ENDDO
          ENDIF
       ENDIF

    ENDIF

!
!-- Bottom-border PEs:
    IF ( nested_bc_at_bottom )  THEN

       CALL pmci_interp_bt( w, wc, 0, iflo, ifuo, jflo, jfuo, 'b', 'w' )
       CALL pmci_interp_bt( u, uc, 0, iflu, ifuu, jflo, jfuo, 'b', 'u' )
       CALL pmci_interp_bt( v, vc, 0, iflo, ifuo, jflv, jfuv, 'b', 'v' )
!
!--    Treatment of TKE. Interpolation is only required if parent and child operate in RANS mode,
!--    else, interpolation is replaced by a Neumann condition.
       IF ( rans_mode_parent  .AND.  rans_mode )  THEN
          CALL pmci_interp_bt( e, ec, 0, iflo, ifuo, jflo, jfuo, 'b', 's' )
       ELSE
          e(nzt+1,nys:nyn,nxl:nxr) = e(nzt,nys:nyn,nxl:nxr)
       ENDIF

       IF ( rans_mode_parent  .AND.  rans_mode  .AND.  rans_tke_e )  THEN
          CALL pmci_interp_bt( diss, dissc, 0, iflo, ifuo, jflo, jfuo, 'b', 's' )
       ENDIF

       IF ( .NOT. neutral )  THEN
          CALL pmci_interp_bt( pt, ptc, 0, iflo, ifuo, jflo, jfuo, 'b', 's' )
       ENDIF

       IF ( humidity )  THEN
          CALL pmci_interp_bt( q, q_c, 0, iflo, ifuo, jflo, jfuo, 'b', 's' )
          IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
             CALL pmci_interp_bt( qc, qcc, 0, iflo, ifuo, jflo, jfuo, 'b', 's' )
             CALL pmci_interp_bt( nc, ncc, 0, iflo, ifuo, jflo, jfuo, 'b', 's' )
          ENDIF
          IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
             CALL pmci_interp_bt( qr, qrc, 0, iflo, ifuo, jflo, jfuo, 'b', 's' )
             CALL pmci_interp_bt( nr, nrc, 0, iflo, ifuo, jflo, jfuo, 'b', 's' )
          ENDIF
       ENDIF

       IF ( passive_scalar )  THEN
          CALL pmci_interp_bt( s, sc, 0, iflo, ifuo, jflo, jfuo, 'b', 's' )
       ENDIF

       IF ( ocean_mode  .AND.  salinity )  THEN
          CALL pmci_interp_bt( sa, sac, kcto, iflo, ifuo, jflo, jfuo, 'b', 's' )
       ENDIF

       IF ( air_chemistry  .AND.  nesting_chem )  THEN
          DO  n = 1, nspec
             CALL pmci_interp_bt( chem_species(n)%conc, chem_spec_c(:,:,:,n), 0, iflo, ifuo,       &
                                  jflo, jfuo, 'b', 's' )
          ENDDO
       ENDIF

       IF ( salsa  .AND.  nesting_salsa )  THEN
          DO  lb = 1, nbins_aerosol
             CALL pmci_interp_bt( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb), 0,          &
                                  iflo, ifuo, jflo, jfuo, 'b', 's' )
          ENDDO
          DO  lc = 1, nbins_aerosol * ncomponents_mass
             CALL pmci_interp_bt( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc), 0,              &
                                  iflo, ifuo, jflo, jfuo, 'b', 's' )
          ENDDO
          IF ( .NOT. salsa_gases_from_chem )  THEN
             DO  lg = 1, ngases_salsa
                CALL pmci_interp_bt( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg), 0, iflo, ifuo,     &
                                     jflo, jfuo, 'b', 's' )
             ENDDO
          ENDIF
       ENDIF

    ENDIF

 END SUBROUTINE pmci_interpolation



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> A wrapper routine for all anterpolation actions. Note that TKE is not anterpolated.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_anterpolation

    IMPLICIT NONE
    INTEGER(iwp) ::  lb  !< running index for aerosol size bins
    INTEGER(iwp) ::  lc  !< running index for aerosol mass bins
    INTEGER(iwp) ::  lg  !< running index for salsa gases
    INTEGER(iwp) ::  n   !< running index for number of chemical species


    CALL pmci_anterp_var( u,  uc,  kcto, iflu, ifuu, jflo, jfuo, kflo, kfuo, ijkfc_u, 'u' )
    CALL pmci_anterp_var( v,  vc,  kcto, iflo, ifuo, jflv, jfuv, kflo, kfuo, ijkfc_v, 'v' )
    CALL pmci_anterp_var( w,  wc,  kctw, iflo, ifuo, jflo, jfuo, kflw, kfuw, ijkfc_w, 'w' )
!
!-- Anterpolation of TKE and dissipation rate if parent and child are in
!-- RANS mode.
    IF ( rans_mode_parent  .AND.  rans_mode )  THEN
       CALL pmci_anterp_var( e, ec, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'e' )
!
!--    Anterpolation of dissipation rate only if TKE-e closure is applied.
       IF ( rans_tke_e )  THEN
          CALL pmci_anterp_var( diss, dissc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s,    &
                                'diss' )
       ENDIF

    ENDIF

    IF ( .NOT. neutral )  THEN
       CALL pmci_anterp_var( pt, ptc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'pt' )
    ENDIF

    IF ( humidity )  THEN

       CALL pmci_anterp_var( q, q_c, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'q' )

       IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
          CALL pmci_anterp_var( qc, qcc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'qc' )
          CALL pmci_anterp_var( nc, ncc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'nc' )
       ENDIF

       IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
          CALL pmci_anterp_var( qr, qrc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'qr' )
          CALL pmci_anterp_var( nr, nrc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'nr' )
       ENDIF

    ENDIF

    IF ( passive_scalar )  THEN
       CALL pmci_anterp_var( s, sc, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
    ENDIF

    IF ( ocean_mode  .AND.  salinity )  THEN
       CALL pmci_anterp_var( sa, sac, kcto, iflo, ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 'sa' )
    ENDIF

    IF ( air_chemistry  .AND.  nesting_chem )  THEN
       DO  n = 1, nspec
          CALL pmci_anterp_var( chem_species(n)%conc, chem_spec_c(:,:,:,n), kcto, iflo, ifuo,      &
                                jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
       ENDDO
    ENDIF

    IF ( salsa  .AND.  nesting_salsa )  THEN
       DO  lb = 1, nbins_aerosol
          CALL pmci_anterp_var( aerosol_number(lb)%conc, aerosol_number_c(:,:,:,lb), kcto, iflo,   &
                                ifuo, jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
       ENDDO
       DO  lc = 1, nbins_aerosol * ncomponents_mass
          CALL pmci_anterp_var( aerosol_mass(lc)%conc, aerosol_mass_c(:,:,:,lc), kcto, iflo, ifuo, &
                                jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
       ENDDO
       IF ( .NOT. salsa_gases_from_chem )  THEN
          DO  lg = 1, ngases_salsa
             CALL pmci_anterp_var( salsa_gas(lg)%conc, salsa_gas_c(:,:,:,lg), kcto, iflo, ifuo,    &
                                   jflo, jfuo, kflo, kfuo, ijkfc_s, 's' )
          ENDDO
       ENDIF
    ENDIF

 END SUBROUTINE pmci_anterpolation


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation of ghost-node values used as the child-domain boundary conditions. This subroutine
!> handles the left and right boundaries.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_interp_lr( child_array, parent_array, kct, jfl, jfu, kfl, kfu, edge, var )

    IMPLICIT NONE

    CHARACTER(LEN=1), INTENT(IN) ::  edge  !< edge symbol: 'l' or 'r'
    CHARACTER(LEN=1), INTENT(IN) ::  var   !< variable symbol: 'u', 'v', 'w' or 's'

    INTEGER(iwp), INTENT(IN) ::  kct  !< the parent-grid index in z-direction just below the boundary value node

    INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfl  !< indicates start index of child cells belonging to certain
                                                            !< parent cell - y direction
    INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfu  !< indicates end index of child cells belonging to certain
                                                            !< parent cell - y direction
    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfl  !< indicates start index of child cells belonging to certain
                                                            !< parent cell - z direction
    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfu  !< indicates end index of child cells belonging to certain
                                                            !< parent cell - z direction

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) ::  child_array  !< Child-grid array

    REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(IN) ::  parent_array  !< Parent-grid array

!
!-- Local variables:
    INTEGER(iwp) ::  icb    !< fixed child-grid index in the x-direction pointing to the node just behind the
                            !< boundary-value node
    INTEGER(iwp) ::  icbc   !< fixed child-grid index in the x-direction pointing to the boundary-value nodes
    INTEGER(iwp) ::  icbgp  !< index running over the redundant boundary ghost points in the x-direction
    INTEGER(iwp) ::  ierr   !< MPI error code
    INTEGER(iwp) ::  ipbeg  !< parent-grid index in the x-direction pointing to the starting point of workarr_lr
                            !< in the parent-grid array
    INTEGER(iwp) ::  ipw    !< reduced parent-grid index in the x-direction for workarr_lr pointing to
                            !< the boundary ghost node
    INTEGER(iwp) ::  ipwp   !< reduced parent-grid index in the x-direction for workarr_lr pointing to
                            !< the first prognostic node
    INTEGER(iwp) ::  jc     !< running child-grid index in the y-direction
    INTEGER(iwp) ::  jp     !< running parent-grid index in the y-direction
    INTEGER(iwp) ::  kc     !< running child-grid index in the z-direction
    INTEGER(iwp) ::  kp     !< running parent-grid index in the z-direction

    REAL(wp) ::  cb          !< interpolation coefficient for the boundary ghost node
    REAL(wp) ::  cp          !< interpolation coefficient for the first prognostic node
    REAL(wp) ::  c_interp_1  !< value interpolated to the flux point in x direction from the parent-grid data
    REAL(wp) ::  c_interp_2  !< auxiliary value interpolated  to the flux point in x direction from the parent-grid data

!
!-- Check which edge is to be handled
    IF ( edge == 'l' )  THEN
!
!--    For u, nxl is a ghost node, but not for the other variables
       IF ( var == 'u' )  THEN
          icbc  = nxl
          icb   = icbc - 1
          ipw   = 2
          ipwp  = ipw        ! This is redundant when var == 'u'
          ipbeg = ipl
       ELSE
          icbc  = nxl - 1
          icb   = icbc - 1
          ipw   = 1
          ipwp  = 2
          ipbeg = ipl
       ENDIF
    ELSEIF ( edge == 'r' )  THEN
       IF ( var == 'u' )  THEN
          icbc  = nxr + 1
          icb   = icbc + 1
          ipw   = 1
          ipwp  = ipw        ! This is redundant when var == 'u'
          ipbeg = ipr - 2
       ELSE
          icbc  = nxr + 1
          icb   = icbc + 1
          ipw   = 1
          ipwp  = 0
          ipbeg = ipr - 2
       ENDIF
    ENDIF
!
!-- Interpolation coefficients
    IF ( interpolation_scheme_lrsn == 1 )  THEN
       cb = 1.0_wp  ! 1st-order upwind
    ELSEIF ( interpolation_scheme_lrsn == 2 )  THEN
       cb = 0.5_wp  ! 2nd-order central
    ELSE
       cb = 0.5_wp  ! 2nd-order central (default)
    ENDIF
    cp = 1.0_wp - cb
!
!-- Substitute the necessary parent-grid data to the work array workarr_lr.
    workarr_lr = 0.0_wp
    IF ( npey > 1 )  THEN

       IF ( bc_dirichlet_s )  THEN
          workarr_lr(0:pg%nz+1,jpsw:jpnw-1,0:2) = parent_array(0:pg%nz+1,jpsw:jpnw-1,ipbeg:ipbeg+2)
       ELSEIF ( bc_dirichlet_n )  THEN
          workarr_lr(0:pg%nz+1,jpsw+1:jpnw,0:2) = parent_array(0:pg%nz+1,jpsw+1:jpnw,ipbeg:ipbeg+2)
       ELSE
          workarr_lr(0:pg%nz+1,jpsw+1:jpnw-1,0:2) =                                                &
                                                parent_array(0:pg%nz+1,jpsw+1:jpnw-1,ipbeg:ipbeg+2)
       ENDIF
!
!--    South-north exchange if more than one PE subdomain in the y-direction. Note that in case of
!--    3-D nesting the south (psouth == MPI_PROC_NULL) and north (pnorth == MPI_PROC_NULL)
!--    boundaries are not exchanged because the nest domain is not cyclic.
!--    From south to north
       CALL MPI_SENDRECV( workarr_lr(0,jpsw+1,0), 1, workarr_lr_exchange_type, psouth,  0,         &
                          workarr_lr(0,jpnw,0), 1, workarr_lr_exchange_type, pnorth,  0, comm2d,   &
                          status, ierr )
!
!--    From north to south
       CALL MPI_SENDRECV( workarr_lr(0,jpnw-1,0), 1, workarr_lr_exchange_type, pnorth,  1,         &
                          workarr_lr(0,jpsw,0), 1, workarr_lr_exchange_type, psouth,  1, comm2d,   &
                          status, ierr )

    ELSE
       workarr_lr(0:pg%nz+1,jpsw:jpnw,0:2) = parent_array(0:pg%nz+1,jpsw:jpnw,ipbeg:ipbeg+2)
    ENDIF

    IF ( var == 'u' )  THEN

       DO  jp = jpsw, jpnw
          DO  kp = 0, kct

             DO  jc = jfl(jp), jfu(jp)
                DO  kc = kfl(kp), kfu(kp)
                   child_array(kc,jc,icbc) = workarr_lr(kp,jp,ipw)
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    ELSEIF ( var == 'v' )  THEN

       DO  jp = jpsw, jpnw-1
          DO  kp = 0, kct
!
!--          First interpolate to the flux point
             c_interp_1 = cb * workarr_lr(kp,jp,ipw)   + cp * workarr_lr(kp,jp,ipwp)
             c_interp_2 = cb * workarr_lr(kp,jp+1,ipw) + cp * workarr_lr(kp,jp+1,ipwp)
!
!--          Use averages of the neighbouring matching grid-line values
             DO  jc = jfl(jp), jfl(jp+1)
                child_array(kfl(kp):kfu(kp),jc,icbc) = 0.5_wp * ( c_interp_1 + c_interp_2 )
             ENDDO
!
!--          Then set the values along the matching grid-lines
             IF  ( MOD( jfl(jp), jgsr ) == 0 )  THEN
                child_array(kfl(kp):kfu(kp),jfl(jp),icbc) = c_interp_1
             ENDIF
          ENDDO
       ENDDO
!
!--    Finally, set the values along the last matching grid-line
       IF ( MOD( jfl(jpnw), jgsr ) == 0 )  THEN
          DO  kp = 0, kct
             c_interp_1 = cb * workarr_lr(kp,jpnw,ipw) + cp * workarr_lr(kp,jpnw,ipwp)
             child_array(kfl(kp):kfu(kp),jfl(jpnw),icbc) = c_interp_1
          ENDDO
       ENDIF
!
!--    A gap may still remain in some cases if the subdomain size is not divisible by the
!--    grid-spacing ratio. In such a case, fill the gap. Note however, this operation may produce
!--    some additional momentum conservation error.
       IF ( jfl(jpnw) < nyn )  THEN
          DO  kp = 0, kct
             DO  jc = jfl(jpnw) + 1, nyn
                child_array(kfl(kp):kfu(kp),jc,icbc) = child_array(kfl(kp):kfu(kp),jfl(jpnw),icbc)
             ENDDO
          ENDDO
       ENDIF

    ELSEIF ( var == 'w' )  THEN

       DO  jp = jpsw, jpnw
          DO  kp = 0, kct + 1   ! It is important to go up to kct+1
!
!--          Interpolate to the flux point
             c_interp_1 = cb * workarr_lr(kp,jp,ipw) + cp * workarr_lr(kp,jp,ipwp)
!
!--          First substitute only the matching-node values
             child_array(kfu(kp),jfl(jp):jfu(jp),icbc) = c_interp_1

          ENDDO
       ENDDO

       DO  jp = jpsw, jpnw
          DO  kp = 1, kct + 1   ! It is important to go up to kct+1
!
!--          Then fill up the nodes in between with the averages
             DO  kc = kfu(kp-1) + 1, kfu(kp) - 1
                child_array(kc,jfl(jp):jfu(jp),icbc) = 0.5_wp * (                                  &
                                                       child_array(kfu(kp-1),jfl(jp):jfu(jp),icbc) &
                                                     + child_array(kfu(kp),jfl(jp):jfu(jp),icbc)   &
                                                                )
             ENDDO

          ENDDO
       ENDDO

    ELSE   ! Any scalar

       DO  jp = jpsw, jpnw
          DO  kp = 0, kct
!
!--          Interpolate to the flux point
             c_interp_1 = cb * workarr_lr(kp,jp,ipw) + cp * workarr_lr(kp,jp,ipwp)
             DO  jc = jfl(jp), jfu(jp)
                DO  kc = kfl(kp), kfu(kp)
                   child_array(kc,jc,icbc) = c_interp_1
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    ENDIF  ! var
!
!-- Fill up also the redundant 2nd and 3rd ghost-node layers
    IF ( edge == 'l' )  THEN
       DO  icbgp = -nbgp, icb
          child_array(0:nzt+1,nysg:nyng,icbgp) = child_array(0:nzt+1,nysg:nyng,icbc)
       ENDDO
    ELSEIF ( edge == 'r' )  THEN
       DO  icbgp = icb, nx+nbgp
          child_array(0:nzt+1,nysg:nyng,icbgp) = child_array(0:nzt+1,nysg:nyng,icbc)
       ENDDO
    ENDIF

 END SUBROUTINE pmci_interp_lr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation of ghost-node values used as the child-domain boundary conditions. This subroutine
!> handles the south and north boundaries.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_interp_sn( child_array, parent_array, kct, ifl, ifu, kfl, kfu, edge, var )

    IMPLICIT NONE

    CHARACTER(LEN=1), INTENT(IN) ::  edge  !< edge symbol: 's' or 'n'
    CHARACTER(LEN=1), INTENT(IN) ::  var   !< variable symbol: 'u', 'v', 'w' or 's'

    INTEGER(iwp), INTENT(IN) ::  kct  !< the parent-grid index in z-direction just below the boundary value node

    INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifl  !< indicates start index of child cells belonging to certain
                                                            !< parent cell - x direction
    INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifu  !< indicates end index of child cells belonging to certain
                                                            !< parent cell - x direction
    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfl  !< indicates start index of child cells belonging to certain
                                                            !< parent cell - z direction
    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfu  !< indicates end index of child cells belonging to certain
                                                            !< parent cell - z direction

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) ::  child_array  !< Child-grid array

    REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(IN) ::  parent_array  !< Parent-grid array
!
!-- Local variables:
    INTEGER(iwp) ::  ic     !< running child-grid index in the x-direction
    INTEGER(iwp) ::  ierr   !< MPI error code
    INTEGER(iwp) ::  ip     !< running parent-grid index in the x-direction
    INTEGER(iwp) ::  jcb    !< fixed child-grid index in the y-direction pointing to the node just behind the
                            !< boundary-value node
    INTEGER(iwp) ::  jcbc   !< fixed child-grid index in the y-direction pointing to the boundary-value nodes
    INTEGER(iwp) ::  jcbgp  !< index running over the redundant boundary ghost points in y-direction
    INTEGER(iwp) ::  jpbeg  !< parent-grid index in the y-direction pointing to the starting point of workarr_sn
                            !< in the parent-grid array
    INTEGER(iwp) ::  jpw    !< reduced parent-grid index in the y-direction for workarr_sn pointing to
                            !< the boundary ghost node
    INTEGER(iwp) ::  jpwp   !< reduced parent-grid index in the y-direction for workarr_sn pointing to
                            !< the first prognostic node
    INTEGER(iwp) ::  kc     !< running child-grid index in the z-direction
    INTEGER(iwp) ::  kp     !< running parent-grid index in the z-direction

    REAL(wp) ::  cb          !< interpolation coefficient for the boundary ghost node
    REAL(wp) ::  cp          !< interpolation coefficient for the first prognostic node
    REAL(wp) ::  c_interp_1  !< value interpolated to the flux point in x direction from the parent-grid data
    REAL(wp) ::  c_interp_2  !< auxiliary value interpolated  to the flux point in x direction from the parent-grid data

!
!-- Check which edge is to be handled: south or north
    IF ( edge == 's' )  THEN
!
!--    For v, nys is a ghost node, but not for the other variables
       IF ( var == 'v' )  THEN
          jcbc  = nys
          jcb   = jcbc - 1
          jpw   = 2
          jpwp  = 2         ! This is redundant when var == 'v'
          jpbeg = jps
       ELSE
          jcbc  = nys - 1
          jcb   = jcbc - 1
          jpw   = 1
          jpwp  = 2
          jpbeg = jps
       ENDIF
    ELSEIF ( edge == 'n' )  THEN
       IF ( var == 'v' )  THEN
          jcbc  = nyn + 1
          jcb   = jcbc + 1
          jpw   = 1
          jpwp  = 0         ! This is redundant when var == 'v'
          jpbeg = jpn - 2
       ELSE
          jcbc  = nyn + 1
          jcb   = jcbc + 1
          jpw   = 1
          jpwp  = 0
          jpbeg = jpn - 2
       ENDIF
    ENDIF
!
!-- Interpolation coefficients
    IF ( interpolation_scheme_lrsn == 1 )  THEN
       cb = 1.0_wp  ! 1st-order upwind
    ELSEIF ( interpolation_scheme_lrsn == 2 )  THEN
       cb = 0.5_wp  ! 2nd-order central
    ELSE
       cb = 0.5_wp  ! 2nd-order central (default)
    ENDIF
    cp = 1.0_wp - cb
!
!-- Substitute the necessary parent-grid data to the work array workarr_sn.
    workarr_sn = 0.0_wp
    IF ( npex > 1 )  THEN

       IF ( bc_dirichlet_l )  THEN
          workarr_sn(0:pg%nz+1,0:2,iplw:iprw-1) = parent_array(0:pg%nz+1,jpbeg:jpbeg+2,iplw:iprw-1)
       ELSEIF ( bc_dirichlet_r )  THEN
          workarr_sn(0:pg%nz+1,0:2,iplw+1:iprw) = parent_array(0:pg%nz+1,jpbeg:jpbeg+2,iplw+1:iprw)
       ELSE
          workarr_sn(0:pg%nz+1,0:2,iplw+1:iprw-1) =                                                &
                                                parent_array(0:pg%nz+1,jpbeg:jpbeg+2,iplw+1:iprw-1)
       ENDIF
!
!--    Left-right exchange if more than one PE subdomain in the x-direction. Note that in case of
!--    3-D nesting the left (pleft == MPI_PROC_NULL) and right (pright == MPI_PROC_NULL) boundaries
!--    are not exchanged because the nest domain is not cyclic.
!--    From left to right
       CALL MPI_SENDRECV( workarr_sn(0,0,iplw+1), 1, workarr_sn_exchange_type, pleft,   0,         &
                          workarr_sn(0,0,iprw), 1, workarr_sn_exchange_type, pright, 0, comm2d,    &
                          status, ierr )
!
!--    From right to left
       CALL MPI_SENDRECV( workarr_sn(0,0,iprw-1), 1, workarr_sn_exchange_type, pright,  1,         &
                          workarr_sn(0,0,iplw), 1, workarr_sn_exchange_type, pleft, 1, comm2d,     &
                          status, ierr )

    ELSE
       workarr_sn(0:pg%nz+1,0:2,iplw:iprw) = parent_array(0:pg%nz+1,jpbeg:jpbeg+2,iplw:iprw)
    ENDIF

    IF ( var == 'v' )  THEN

       DO  ip = iplw, iprw
          DO  kp = 0, kct

             DO  ic = ifl(ip), ifu(ip)
                DO  kc = kfl(kp), kfu(kp)
                   child_array(kc,jcbc,ic) = workarr_sn(kp,jpw,ip)
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    ELSEIF ( var == 'u' )  THEN

       DO  ip = iplw, iprw - 1
          DO  kp = 0, kct
!
!--          First interpolate to the flux point
             c_interp_1 = cb * workarr_sn(kp,jpw,ip)   + cp * workarr_sn(kp,jpwp,ip)
             c_interp_2 = cb * workarr_sn(kp,jpw,ip+1) + cp * workarr_sn(kp,jpwp,ip+1)
!
!--          Use averages of the neighbouring matching grid-line values
             DO  ic = ifl(ip), ifl(ip+1)
                child_array(kfl(kp):kfu(kp),jcbc,ic) = 0.5_wp * ( c_interp_1 + c_interp_2 )
             ENDDO
!
!--          Then set the values along the matching grid-lines
             IF ( MOD( ifl(ip), igsr ) == 0 )  THEN
                child_array(kfl(kp):kfu(kp),jcbc,ifl(ip)) = c_interp_1
             ENDIF

          ENDDO
       ENDDO
!
!--    Finally, set the values along the last matching grid-line
       IF ( MOD( ifl(iprw), igsr ) == 0 )  THEN
          DO  kp = 0, kct
             c_interp_1 = cb * workarr_sn(kp,jpw,iprw) + cp * workarr_sn(kp,jpwp,iprw)
             child_array(kfl(kp):kfu(kp),jcbc,ifl(iprw)) = c_interp_1
          ENDDO
       ENDIF
!
!--    A gap may still remain in some cases if the subdomain size is not divisible by the
!--    grid-spacing ratio. In such a case, fill the gap. Note however, this operation may produce
!--    some additional momentum conservation error.
       IF ( ifl(iprw) < nxr )  THEN
          DO  kp = 0, kct
             DO  ic = ifl(iprw) + 1, nxr
                child_array(kfl(kp):kfu(kp),jcbc,ic) = child_array(kfl(kp):kfu(kp),jcbc,ifl(iprw))
             ENDDO
          ENDDO
       ENDIF

    ELSEIF ( var == 'w' )  THEN

       DO  ip = iplw, iprw
          DO  kp = 0, kct + 1   ! It is important to go up to kct+1
!
!--          Interpolate to the flux point
             c_interp_1 = cb * workarr_sn(kp,jpw,ip) + cp * workarr_sn(kp,jpwp,ip)
!
!--          First substitute only the matching-node values
             child_array(kfu(kp),jcbc,ifl(ip):ifu(ip)) = c_interp_1

          ENDDO
       ENDDO

       DO  ip = iplw, iprw
          DO  kp = 1, kct + 1   ! It is important to go up to kct + 1
!
!--          Then fill up the nodes in between with the averages
             DO  kc = kfu(kp-1) + 1, kfu(kp) - 1
                child_array(kc,jcbc,ifl(ip):ifu(ip)) = 0.5_wp * (                                  &
                                                       child_array(kfu(kp-1),jcbc,ifl(ip):ifu(ip)) &
                                                     + child_array(kfu(kp),jcbc,ifl(ip):ifu(ip))   &
                                                                )
             ENDDO

          ENDDO
       ENDDO

    ELSE   ! Any scalar

       DO  ip = iplw, iprw
          DO  kp = 0, kct
!
!--          Interpolate to the flux point
             c_interp_1 = cb * workarr_sn(kp,jpw,ip) + cp * workarr_sn(kp,jpwp,ip)
             DO  ic = ifl(ip), ifu(ip)
                DO  kc = kfl(kp), kfu(kp)
                   child_array(kc,jcbc,ic) = c_interp_1
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    ENDIF  ! var
!
!-- Fill up also the redundant 2nd and 3rd ghost-node layers
    IF ( edge == 's' )  THEN
       DO  jcbgp = -nbgp, jcb
          child_array(0:nzt+1,jcbgp,nxlg:nxrg) = child_array(0:nzt+1,jcbc,nxlg:nxrg)
       ENDDO
    ELSEIF ( edge == 'n' )  THEN
       DO  jcbgp = jcb, ny+nbgp
          child_array(0:nzt+1,jcbgp,nxlg:nxrg) = child_array(0:nzt+1,jcbc,nxlg:nxrg)
       ENDDO
    ENDIF

 END SUBROUTINE pmci_interp_sn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Interpolation of ghost-node values used as the child-domain boundary conditions. This subroutine
!> handles the bottom and top boundary.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_interp_bt( child_array, parent_array, kct, ifl, ifu, jfl, jfu, edge, var )

    IMPLICIT NONE

    CHARACTER(LEN=1), INTENT(IN) ::  edge  !< edge symbol: 'b' or 't'
    CHARACTER(LEN=1), INTENT(IN) ::  var   !< variable grid symbol: 'u', 'v', 'w' or 's'

    INTEGER(iwp), INTENT(IN) ::  kct  !< the parent-grid index in z-direction just below the top-boundary-value node
                                      !< however for the bottom boundary kct = 0 i.e. the boundary-value node

    INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifl  !< indicates start index of child cells belonging to certain
                                                            !< parent cell - x direction
    INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifu  !< indicates end index of child cells belonging to certain
                                                            !< parent cell - x direction
    INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfl  !< indicates start index of child cells belonging to certain
                                                            !< parent cell - y direction
    INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfu  !< indicates end index of child cells belonging to certain
                                                            !< parent cell - y direction

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(INOUT) ::  child_array  !< Child-grid array

    REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(IN) ::  parent_array  !< Parent-grid array

!
!-- Local variables:
    INTEGER(iwp) ::  ic          !< running child-grid index in the x-direction
    INTEGER(iwp) ::  ierr        !< MPI error code
    INTEGER(iwp) ::  iplc        !< lower parent-grid index limit in the x-direction for copying parent-grid
                                 !< array data to workarr_bt
    INTEGER(iwp) ::  iprc        !< upper parent-grid index limit in the x-direction for copying parent-grid
                                 !< array data to workarr_bt
    INTEGER(iwp) ::  jc          !< running child-grid index in the y-direction
    INTEGER(iwp) ::  jpsc        !< lower parent-grid index limit in the y-direction for copying parent-grid
                                 !< array data to workarr_bt
    INTEGER(iwp) ::  jpnc        !< upper parent-grid-index limit in the y-direction for copying parent-grid
                                 !< array data to workarr_bt
    INTEGER(iwp) ::  kc          !< vertical child-grid index fixed to the boundary-value level
    INTEGER(iwp) ::  ip          !< running parent-grid index in the x-direction
    INTEGER(iwp) ::  jp          !< running parent-grid index in the y-direction
    INTEGER(iwp) ::  kpw         !< reduced parent-grid index in the z-direction for workarr_bt pointing to
                                 !< the boundary ghost node
    REAL(wp) ::  c31         !< interpolation coefficient for the 3rd-order WS scheme
    REAL(wp) ::  c32         !< interpolation coefficient for the 3rd-order WS scheme
    REAL(wp) ::  c33         !< interpolation coefficient for the 3rd-order WS scheme
    REAL(wp) ::  c_interp_1  !< value interpolated to the flux point in z direction from the parent-grid data
    REAL(wp) ::  c_interp_2  !< auxiliary value interpolated to the flux point in z direction from the parent-grid data


    IF ( edge == 't' )  THEN
       IF ( var == 'w' )  THEN
          kc = nzt
       ELSE
          kc = nzt + 1
       ENDIF
    ELSEIF (edge == 'b' )  THEN
       kc = 0
    ENDIF
    kpw = 1
!
!-- Interpolation coefficients. The 3rd-order WS scheme cannot be used for the bottom boundary.
!-- Therefore the 2nd-order central is used for the bottom boundary.    
    IF ( interpolation_scheme_t == 1 )  THEN
       c31 =  0.0_wp           ! 1st-order upwind
       c32 =  1.0_wp
       c33 =  0.0_wp
    ELSEIF ( interpolation_scheme_t == 2  .OR.  edge == 'b' )  THEN
       c31 =  0.5_wp           ! 2nd-order central (default for the bottom boundary) 
       c32 =  0.5_wp
       c33 =  0.0_wp
    ELSE
       c31 =  2.0_wp / 6.0_wp  ! 3rd-order WS upwind biased (default for the top boundary)
       c32 =  5.0_wp / 6.0_wp
       c33 = -1.0_wp / 6.0_wp
    ENDIF
!
!-- Substitute the necessary parent-grid data to the work array. Note that the dimension of
!-- workarr_bt is (0:2,jpsw:jpnw,iplw:iprw) and the jc?w and ic?w-index bounds depend on the
!-- location of the PE-subdomain relative to the side boundaries.
    iplc = iplw + 1
    iprc = iprw - 1
    jpsc = jpsw + 1
    jpnc = jpnw - 1
    IF ( bc_dirichlet_l  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       iplc = iplw
    ENDIF
    IF ( bc_dirichlet_r  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       iprc = iprw
    ENDIF
    IF ( bc_dirichlet_s  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       jpsc = jpsw
    ENDIF
    IF ( bc_dirichlet_n  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       jpnc = jpnw
    ENDIF
    workarr_bt = 0.0_wp
!    
!-- Note that kct = 0 for the bottom boundary
    workarr_bt(0:2,jpsc:jpnc,iplc:iprc) = parent_array(kct:kct+2,jpsc:jpnc,iplc:iprc)
!
!-- Left-right exchange if more than one PE subdomain in the x-direction. Note that in case of
!-- 3-D nesting the left and right boundaries are not exchanged because the nest domain is not
!-- cyclic.
    IF ( npex > 1 )  THEN
!
!--    From left to right
       CALL MPI_SENDRECV( workarr_bt(0,jpsw,iplw+1), 1, workarr_bt_exchange_type_y, pleft, 0,      &
                          workarr_bt(0,jpsw,iprw), 1, workarr_bt_exchange_type_y, pright, 0,       &
                          comm2d, status, ierr )
!
!--    From right to left
       CALL MPI_SENDRECV( workarr_bt(0,jpsw,iprw-1), 1, workarr_bt_exchange_type_y, pright, 1,     &
                          workarr_bt(0,jpsw,iplw), 1, workarr_bt_exchange_type_y, pleft,  1,       &
                          comm2d, status, ierr )
    ENDIF
!
!-- South-north exchange if more than one PE subdomain in the y-direction.
!-- Note that in case of 3-D nesting the south and north boundaries are not exchanged because the
!-- nest domain is not cyclic.
    IF ( npey > 1 )  THEN
!
!--    From south to north
       CALL MPI_SENDRECV( workarr_bt(0,jpsw+1,iplw), 1, workarr_bt_exchange_type_x, psouth, 2,     &
                          workarr_bt(0,jpnw,iplw), 1, workarr_bt_exchange_type_x, pnorth, 2,       &
                          comm2d, status, ierr )
!
!--    From north to south
       CALL MPI_SENDRECV( workarr_bt(0,jpnw-1,iplw), 1, workarr_bt_exchange_type_x, pnorth, 3,     &
                          workarr_bt(0,jpsw,iplw), 1, workarr_bt_exchange_type_x, psouth, 3,       &
                          comm2d, status, ierr )
    ENDIF

    IF  ( var == 'w' )  THEN
       DO  ip = iplw, iprw
          DO  jp = jpsw, jpnw
             
             DO  ic = ifl(ip), ifu(ip)
                DO  jc = jfl(jp), jfu(jp)
                   child_array(kc,jc,ic) = workarr_bt(kpw,jp,ip)
                ENDDO
             ENDDO
             
          ENDDO
       ENDDO

    ELSEIF  ( var == 'u' )  THEN

       DO  ip = iplw, iprw - 1
          DO  jp = jpsw, jpnw
!
!--          First interpolate to the flux point
             c_interp_1 = c31 * workarr_bt(kpw-1,jp,ip)   + c32 * workarr_bt(kpw,jp,ip) +          &
                          c33 * workarr_bt(kpw+1,jp,ip)
             c_interp_2 = c31 * workarr_bt(kpw-1,jp,ip+1) + c32 * workarr_bt(kpw,jp,ip+1) +        &
                          c33 * workarr_bt(kpw+1,jp,ip+1)
!
!--          Use averages of the neighbouring matching grid-line values
             DO  ic = ifl(ip), ifl(ip+1)
                child_array(kc,jfl(jp):jfu(jp),ic) = 0.5_wp * ( c_interp_1 + c_interp_2 )
             ENDDO
!
!--          Then set the values along the matching grid-lines
             IF ( MOD( ifl(ip), igsr ) == 0 )  THEN
!
!--             First interpolate to the flux point
                c_interp_1 = c31 * workarr_bt(kpw-1,jp,ip) + c32 * workarr_bt(kpw,jp,ip) +         &
                             c33 * workarr_bt(kpw+1,jp,ip)
                child_array(kc,jfl(jp):jfu(jp),ifl(ip)) = c_interp_1
             ENDIF

          ENDDO
       ENDDO
!
!--    Finally, set the values along the last matching grid-line
       IF  ( MOD( ifl(iprw), igsr ) == 0 )  THEN
          DO  jp = jpsw, jpnw
!
!--          First interpolate to the flux point using the 3rd-order WS scheme
             c_interp_1 = c31 * workarr_bt(kpw-1,jp,iprw) + c32 * workarr_bt(kpw,jp,iprw) +        &
                          c33 * workarr_bt(kpw+1,jp,iprw)
             child_array(kc,jfl(jp):jfu(jp),ifl(iprw)) = c_interp_1
          ENDDO
       ENDIF
!
!--    A gap may still remain in some cases if the subdomain size is not divisible by the
!--    grid-spacing ratio. In such a case, fill the gap. Note however, this operation may produce
!--    some additional momentum conservation error.
       IF ( ifl(iprw) < nxr )  THEN
          DO  jp = jpsw, jpnw
             DO  ic = ifl(iprw) + 1, nxr
                child_array(kc,jfl(jp):jfu(jp),ic) = child_array(kc,jfl(jp):jfu(jp),ifl(iprw))
             ENDDO
          ENDDO
       ENDIF

    ELSEIF  ( var == 'v' )  THEN

       DO  ip = iplw, iprw
          DO  jp = jpsw, jpnw-1
!
!--          First interpolate to the flux point using the 3rd-order WS scheme
             c_interp_1 = c31 * workarr_bt(kpw-1,jp,ip)   + c32 * workarr_bt(kpw,jp,ip) +          &
                          c33 * workarr_bt(kpw+1,jp,ip)
             c_interp_2 = c31 * workarr_bt(kpw-1,jp+1,ip) + c32 * workarr_bt(kpw,jp+1,ip) +        &
                          c33 * workarr_bt(kpw+1,jp+1,ip)
!
!--          Use averages of the neighbouring matching grid-line values
             DO  jc = jfl(jp), jfl(jp+1)
                child_array(kc,jc,ifl(ip):ifu(ip)) = 0.5_wp * ( c_interp_1 + c_interp_2 )
             ENDDO
!
!--          Then set the values along the matching grid-lines
             IF ( MOD( jfl(jp), jgsr ) == 0 )  THEN
                c_interp_1 = c31 * workarr_bt(kpw-1,jp,ip) + c32 * workarr_bt(kpw,jp,ip) +         &
                             c33 * workarr_bt(kpw+1,jp,ip)
                child_array(kc,jfl(jp),ifl(ip):ifu(ip)) = c_interp_1
             ENDIF
             
          ENDDO

       ENDDO
!
!--    Finally, set the values along the last matching grid-line
       IF ( MOD( jfl(jpnw), jgsr ) == 0 )  THEN
          DO  ip = iplw, iprw
!
!--          First interpolate to the flux point using the 3rd-order WS scheme
             c_interp_1 = c31 * workarr_bt(kpw-1,jpnw,ip) + c32 * workarr_bt(kpw,jpnw,ip) +        &
                          c33 * workarr_bt(kpw+1,jpnw,ip)
             child_array(kc,jfl(jpnw),ifl(ip):ifu(ip)) = c_interp_1
          ENDDO
       ENDIF
!
!--    A gap may still remain in some cases if the subdomain size is not divisible by the
!--    grid-spacing ratio. In such a case, fill the gap. Note however, this operation may produce
!--    some additional momentum conservation error.
       IF  ( jfl(jpnw) < nyn )  THEN
          DO  ip = iplw, iprw
             DO  jc = jfl(jpnw)+1, nyn
                child_array(kc,jc,ifl(ip):ifu(ip)) = child_array(kc,jfl(jpnw),ifl(ip):ifu(ip))
             ENDDO
          ENDDO
       ENDIF

    ELSE  ! Any scalar variable

       DO  ip = iplw, iprw
          DO  jp = jpsw, jpnw
!
!--          First interpolate to the flux point using the 3rd-order WS scheme
             c_interp_1 = c31 * workarr_bt(kpw-1,jp,ip) + c32 * workarr_bt(kpw,jp,ip) +            &
                          c33 * workarr_bt(kpw+1,jp,ip)
             DO  ic = ifl(ip), ifu(ip)
                DO  jc = jfl(jp), jfu(jp)
                   child_array(kc,jc,ic) = c_interp_1
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    ENDIF  ! var
!
!-- Just fill up the redundant second ghost-node layer in case of edge == t and var == w.
    IF ( edge == 't'  .AND.  var == 'w' )  THEN
       child_array(nzt+1,:,:) = child_array(nzt,:,:)
    ENDIF

 END SUBROUTINE pmci_interp_bt


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Anterpolation of internal-node values of one variable to be used as the parent-domain values.
!> This subroutine is based on the first-order numerical integration of the child-grid values
!> contained within the anterpolation cell (Clark & Farley, Journal of the Atmospheric
!> Sciences 41(3), 1984).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_anterp_var( child_array, parent_array, kct, ifl, ifu, jfl, jfu, kfl, kfu,         &
                             ijkfc, var )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  var  !< variable symbol: 'u', 'v', 'w' or 's'

    INTEGER(iwp), INTENT(IN) ::  kct  !< top boundary index for anterpolation along z

    INTEGER(iwp), DIMENSION(0:pg%nz+1,jpsa:jpna,ipla:ipra), INTENT(IN) ::  ijkfc  !< number of child grid points contributing
                                                                                  !< to a parent grid box
    INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifl  !< indicates start index of child cells belonging to certain
                                                            !< parent cell - x direction
    INTEGER(iwp), DIMENSION(ipla:ipra), INTENT(IN) ::  ifu  !< indicates end index of child cells belonging to certain
                                                            !< parent cell - x direction
    INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfl  !< indicates start index of child cells belonging to certain
                                                            !< parent cell - y direction
    INTEGER(iwp), DIMENSION(jpsa:jpna), INTENT(IN) ::  jfu  !< indicates end index of child cells belonging to certain
                                                            !< parent cell - y direction
    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfl  !< indicates start index of child cells belonging to certain
                                                            !< parent cell - z direction
    INTEGER(iwp), DIMENSION(0:pg%nz+1), INTENT(IN) ::  kfu  !< indicates end index of child cells belonging to certain
                                                            !< parent cell - z direction

    REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg), INTENT(IN) ::  child_array  !< child-grid array

    REAL(wp), DIMENSION(0:pg%nz+1,jps:jpn,ipl:ipr), INTENT(INOUT) ::  parent_array  !< parent-grid array

!
!-- Local variables:
    INTEGER(iwp) ::  ic              !< running index x-direction - child grid
    INTEGER(iwp) ::  ip              !< running index x-direction - parent grid
    INTEGER(iwp) ::  ipl_anterp      !< left boundary index for anterpolation along x
    INTEGER(iwp) ::  ipr_anterp      !< right boundary index for anterpolation along x
    INTEGER(iwp) ::  jc              !< running index y-direction - child grid
    INTEGER(iwp) ::  jp              !< running index y-direction - parent grid
    INTEGER(iwp) ::  jpn_anterp      !< north boundary index for anterpolation along y
    INTEGER(iwp) ::  jps_anterp      !< south boundary index for anterpolation along y
    INTEGER(iwp) ::  kc              !< running index z-direction - child grid
    INTEGER(iwp) ::  kp              !< running index z-direction - parent grid
    INTEGER(iwp) ::  kpt_anterp      !< top boundary index for anterpolation along z
    INTEGER(iwp) ::  var_flag        !< bit number used to flag topography on respective grid

    REAL(wp) ::  cellsum  !< sum of respective child cells belonging to parent cell

!
!-- Define the index bounds ipl_anterp, ipr_anterp, jps_anterp and jpn_anterp.
!-- Note that kcb_anterp is simply zero and kct_anterp depends on kct which enters here as a
!-- parameter and it is determined in pmci_define_index_mapping. Note that the grid points directly
!-- used also for interpolation (from parent to child) are always excluded from anterpolation, e.g.
!-- anterpolation is maximally only from 0:kct-1, since kct is directly used for interpolation.
!-- Similar restriction is applied to the lateral boundaries as well. An additional buffer is also
!-- applied (default value for anterpolation_buffer_width = 2) in order to avoid unphysical
!-- accumulation of kinetic energy.
    ipl_anterp = ipl
    ipr_anterp = ipr
    jps_anterp = jps
    jpn_anterp = jpn
!
!-- kpb_anterp is a function of jp and ip and it is set in the initialization phase in
!-- pmci_compute_kpb_anterp.
    kpt_anterp = kct - 1 - anterpolation_buffer_width

!
!-- Set the anterpolation buffers on the lateral boundaries (not required, if the child has
!-- cyclic lateral boundaries and if pure vertical nesting is used).
    IF ( .NOT. bc_lr_cyc  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       ipl_anterp = MAX( ipl, iplg + 3 + anterpolation_buffer_width )
       ipr_anterp = MIN( ipr, iprg - 3 - anterpolation_buffer_width )
    ENDIF
    IF ( .NOT. bc_ns_cyc  .AND.  .NOT. nesting_bounds_vertical_only )  THEN
       jps_anterp = MAX( jps, jpsg + 3 + anterpolation_buffer_width )
       jpn_anterp = MIN( jpn, jpng - 3 - anterpolation_buffer_width )
    ENDIF
!
!-- Set masking bit for topography flags
    IF ( var == 'u' )  THEN
       var_flag = 1
    ELSEIF ( var == 'v' )  THEN
       var_flag = 2
    ELSEIF ( var == 'w' )  THEN
       var_flag = 3
    ELSE
       var_flag = 0
    ENDIF
!
!-- Note that ip, jp, and kp are parent-grid indices and ic,jc, and kc are child-grid indices.
    DO  ip = ipl_anterp, ipr_anterp
       DO  jp = jps_anterp, jpn_anterp
!
!--       If the user has set anterpolation_starting_height less than the canopy height, the
!--       anterpolation is made also within buildings for simplicity, and even under elevated
!--       terrain if anterpolation_starting_height is set smaller than terrain height.
          DO  kp = kpb_anterp(jp,ip), kpt_anterp
             cellsum = 0.0_wp
             DO  ic = ifl(ip), ifu(ip)
                DO  jc = jfl(jp), jfu(jp)
                   DO  kc = kfl(kp), kfu(kp)
                      cellsum = cellsum + MERGE( child_array(kc,jc,ic), 0.0_wp,                    &
                                                 BTEST( topo_flags(kc,jc,ic), var_flag ) )
                   ENDDO
                ENDDO
             ENDDO
!
!--          In case all child grid points are inside topography, i.e. ijkfc and cellsum are zero,
!--          also parent solution would have zero values at that grid point, which may cause
!--          problems in particular for the temperature. Therefore, in case cellsum is zero, keep
!--          the parent solution at this point.
             IF ( ijkfc(kp,jp,ip) /= 0 )  THEN
                parent_array(kp,jp,ip) = cellsum / REAL( ijkfc(kp,jp,ip), KIND = wp )
             ENDIF

          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE pmci_anterp_var

#endif

 END SUBROUTINE pmci_child_datatrans


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set boundary conditions for the prognostic quantities after interpolation and anterpolation at
!> upward- and downward facing surfaces.
!> @todo: add Dirichlet boundary conditions for pot. temperature, humdidity and passive scalar.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_boundary_conds

#if defined( __parallel )
    IMPLICIT NONE

    INTEGER(iwp) ::  ic  !< index along x-direction
    INTEGER(iwp) ::  jc  !< index along y-direction
    INTEGER(iwp) ::  kc  !< index along z-direction
    INTEGER(iwp) ::  lb  !< running index for aerosol size bins
    INTEGER(iwp) ::  lc  !< running index for aerosol mass bins
    INTEGER(iwp) ::  lg  !< running index for salsa gases
    INTEGER(iwp) ::  m   !< running index for surface type
    INTEGER(iwp) ::  n   !< running index for number of chemical species


    IF ( debug_output_timestep )  CALL debug_message( 'pmci_boundary_conds', 'start' )
!
!-- Set Dirichlet boundary conditions for horizontal velocity components at topography grid points.
    IF ( ibc_uv_b == 0 )  THEN
       DO  m = 1, bc_hv%ns
!
!--       Only consider horizontal surfaces. This is to maintain comparability to prior versions.
!--       Later on, treat also vertical surfaces.
          IF ( bc_hv%koff(m) /= 0 )  THEN
             ic = bc_hv%i(m)
             jc = bc_hv%j(m)
             kc = bc_hv%k(m)
             u(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = 0.0_wp
             v(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = 0.0_wp
          ENDIF
       ENDDO
    ENDIF
!
!-- Set Dirichlet boundary conditions for vertical velocity component
    DO  m = 1, bc_hv%ns
       IF ( bc_hv%koff(m) /= 0  )  THEN
          ic = bc_hv%i(m)
          jc = bc_hv%j(m)
          kc = bc_hv%k(m)
          w(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = 0.0_wp
       ENDIF
    ENDDO
!
!-- Set Neumann boundary conditions for potential temperature
    IF ( .NOT. neutral )  THEN
       IF ( ibc_pt_b == 1 )  THEN
          DO  m = 1, bc_hv%ns
             IF ( bc_hv%koff(m) /= 0  )  THEN
                ic = bc_hv%i(m)
                jc = bc_hv%j(m)
                kc = bc_hv%k(m)
                pt(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = pt(kc,jc,ic)
             ENDIF
          ENDDO
       ENDIF
    ENDIF
!
!-- Set Neumann boundary conditions for humidity and cloud-physical quantities
    IF ( humidity )  THEN
       IF ( ibc_q_b == 1 )  THEN
          DO  m = 1, bc_hv%ns
             IF ( bc_hv%koff(m) /= 0  )  THEN
                ic = bc_hv%i(m)
                jc = bc_hv%j(m)
                kc = bc_hv%k(m)
                q(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = q(kc,jc,ic)
             ENDIF
          ENDDO
       ENDIF
       IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
          DO  m = 1, bc_hv%ns
             IF ( bc_hv%koff(m) /= 0  )  THEN
                ic = bc_hv%i(m)
                jc = bc_hv%j(m)
                kc = bc_hv%k(m)
                nc(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = 0.0_wp
                qc(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = 0.0_wp
             ENDIF
          ENDDO
       ENDIF

       IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
          DO  m = 1, bc_hv%ns
             IF ( bc_hv%koff(m) /= 0  )  THEN
                ic = bc_hv%i(m)
                jc = bc_hv%j(m)
                kc = bc_hv%k(m)
                nr(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = 0.0_wp
                qr(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = 0.0_wp
             ENDIF
          ENDDO
       ENDIF

    ENDIF
!
!-- Set Neumann boundary conditions for passive scalar
    IF ( passive_scalar )  THEN
       IF ( ibc_s_b == 1 )  THEN
          DO  m = 1, bc_hv%ns
             IF ( bc_hv%koff(m) /= 0  )  THEN
                ic = bc_hv%i(m)
                jc = bc_hv%j(m)
                kc = bc_hv%k(m)
                s(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = s(kc,jc,ic)
             ENDIF
          ENDDO
       ENDIF
    ENDIF
!
!-- Set Neumann boundary conditions for salinity.
    IF ( ocean_mode  .AND.  salinity )  THEN
       IF ( .NOT. nested_bc_at_bottom )  THEN
          DO  m = 1, bc_hv%ns
             IF ( bc_hv%koff(m) /= 0  )  THEN
                ic = bc_hv%i(m)
                jc = bc_hv%j(m)
                kc = bc_hv%k(m)
                sa(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) = sa(kc,jc,ic)
             ENDIF
          ENDDO
       ENDIF
    ENDIF
!
!-- Set Neumann boundary conditions for chemical species
    IF ( air_chemistry  .AND.  nesting_chem )  THEN
       IF ( ibc_cs_b == 1 )  THEN
          DO  n = 1, nspec
             DO  m = 1, bc_hv%ns
                IF ( bc_hv%koff(m) /= 0  )  THEN
                   ic = bc_hv%i(m)
                   jc = bc_hv%j(m)
                   kc = bc_hv%k(m)
                   chem_species(n)%conc(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) =      &
                                                                     chem_species(n)%conc(kc,jc,ic)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDIF
!
!-- Set Neumann boundary conditions for aerosols and salsa gases
    IF ( salsa  .AND.  nesting_salsa )  THEN
       IF ( ibc_aer_b == 1 )  THEN
          DO  m = 1, bc_hv%ns
             IF ( bc_hv%koff(m) /= 0  )  THEN
                ic = bc_hv%i(m)
                jc = bc_hv%j(m)
                kc = bc_hv%k(m)
                DO  lb = 1, nbins_aerosol
                   aerosol_number(lb)%conc(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) =   &
                                                                   aerosol_number(lb)%conc(kc,jc,ic)
                ENDDO
                DO  lc = 1, nbins_aerosol * ncomponents_mass
                   aerosol_mass(lc)%conc(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) =     &
                                                                   aerosol_mass(lc)%conc(kc,jc,ic)
                ENDDO
                IF ( .NOT. salsa_gases_from_chem )  THEN
                   DO  lg = 1, ngases_salsa
                      salsa_gas(lg)%conc(kc+bc_hv%koff(m),jc+bc_hv%joff(m),ic+bc_hv%ioff(m)) =     &
                                                                   salsa_gas(lg)%conc(kc,jc,ic)
                   ENDDO
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'pmci_boundary_conds', 'end' )

#endif
 END SUBROUTINE pmci_boundary_conds



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Adjust the volume-flow rate through the nested boundaries so that the net volume flow through
!> all boundaries of the current nest domain becomes zero.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pmci_ensure_nest_mass_conservation

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !< running index in the x-direction
    INTEGER(iwp) ::  ierr  !< MPI error code
    INTEGER(iwp) ::  j     !< running index in the y-direction
    INTEGER(iwp) ::  k     !< running index in the z-direction
    INTEGER(iwp) ::  n     !< running index over the boundary faces: l, r, s, n and t

    REAL(wp) ::  dxdy                  !< surface area of grid cell top face
    REAL(wp) ::  sub_sum               !< intermediate sum for reducing the loss of signifigant digits in 2-D summations
    REAL(wp) ::  u_corr_left           !< correction added to the left boundary value of u
    REAL(wp) ::  u_corr_right          !< correction added to the right boundary value of u
    REAL(wp) ::  v_corr_south          !< correction added to the south boundary value of v
    REAL(wp) ::  v_corr_north          !< correction added to the north boundary value of v
    REAL(wp) ::  volume_flux_integral  !< surface integral of volume flux over the domain boundaries
    REAL(wp) ::  volume_flux_local     !< surface integral of volume flux over the subdomain boundary face
    REAL(wp) ::  w_corr_top            !< correction added to the top boundary value of w

    REAL(wp), DIMENSION(6) ::  volume_flux  !< surface integral of volume flux over each boundary face of the domain

!
!-- Sum up the volume flow through the left boundary
    volume_flux(1) = 0.0_wp
    volume_flux_local = 0.0_wp
    IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  THEN
       i = 0
       DO   j = nys, nyn
          sub_sum = 0.0_wp
          DO   k = nzb+1, nzt
             sub_sum = sub_sum + dy * dzw(k) * u(k,j,i) *                                          &
                       MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
          volume_flux_local = volume_flux_local + sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(1), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(1) = volume_flux_local
#endif
!
!-- Sum up the volume flow through the right boundary
    volume_flux(2) = 0.0_wp
    volume_flux_local = 0.0_wp
    IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  THEN
       i = nx + 1
       DO   j = nys, nyn
          sub_sum = 0.0_wp
          DO   k = nzb+1, nzt
             sub_sum = sub_sum - dy * dzw(k) * u(k,j,i) *                                          &
                       MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
          volume_flux_local = volume_flux_local + sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(2), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(2) = volume_flux_local
#endif
!
!-- Sum up the volume flow through the south boundary
    volume_flux(3) = 0.0_wp
    volume_flux_local = 0.0_wp
    IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  THEN
       j = 0
       DO   i = nxl, nxr
          sub_sum = 0.0_wp
          DO   k = nzb+1, nzt
             sub_sum = sub_sum + dx * dzw(k) * v(k,j,i) *                                          &
                       MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO
          volume_flux_local = volume_flux_local + sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(3), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(3) = volume_flux_local
#endif
!
!-- Sum up the volume flow through the north boundary
    volume_flux(4) = 0.0_wp
    volume_flux_local = 0.0_wp
    IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  THEN
       j = ny + 1
       DO  i = nxl, nxr
          sub_sum = 0.0_wp
          DO  k = nzb+1, nzt
             sub_sum = sub_sum - dx * dzw(k) * v(k,j,i) *                                          &
                       MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
          ENDDO
          volume_flux_local = volume_flux_local + sub_sum
       ENDDO
    ENDIF

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(4), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(4) = volume_flux_local
#endif
!
!-- Sum up the volume flow through the top boundary
    volume_flux(5) = 0.0_wp
    volume_flux_local = 0.0_wp
    dxdy = dx * dy
    k = nzt
    DO  i = nxl, nxr
       sub_sum = 0.0_wp
       DO   j = nys, nyn
          sub_sum = sub_sum - dxdy * w(k,j,i) *                                                    &
                    MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
       ENDDO
       volume_flux_local = volume_flux_local + sub_sum
    ENDDO

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(5), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(5) = volume_flux_local
#endif

!
!-- Sum up the volume flow through the bottom boundary
    volume_flux(6) = 0.0_wp
    volume_flux_local = 0.0_wp
    dxdy = dx * dy
    k = 0
    DO  i = nxl, nxr
       sub_sum = 0.0_wp
       DO   j = nys, nyn
          sub_sum = sub_sum + dxdy * w(k,j,i) *                                                    &
                    MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 3 ) )
       ENDDO
       volume_flux_local = volume_flux_local + sub_sum
    ENDDO

#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( volume_flux_local, volume_flux(6), 1, MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    volume_flux(6) = volume_flux_local
#endif

    volume_flux_integral = 0.0_wp
    DO  n = 1, 6
       volume_flux_integral = volume_flux_integral + volume_flux(n)
    ENDDO
!
!-- Adjust the velocities.
    IF ( nesting_bounds_vertical_only )  THEN
!
!--    For pure vertical nesting and non-cyclic conditions along x or y only the top boundary must
!--    compensate, because velocity values at boundaries along x or y are prescribed elsewhere
!--    (via turbulence recycling, Dirichlet- or radiation boundary conditions). Therefore,
!--    normal velocities (u,v) at these boundaries must not be changed here. They would be
!--    overwritten anyhow at a later time (e.g. when turbulence recycling is carried out).
!--    For the case with cyclic conditions along x and y no correction is required at all
!--    (horizontally averaged vertical velocity may be removed in the pressure solver).
       IF ( .NOT. ( bc_lr_cyc  .AND.  bc_ns_cyc ) )  THEN
          w_corr_top = volume_flux_integral / ( face_area(5) + face_area(6) )
       ELSE
          w_corr_top = 0.0_wp
       ENDIF
       u_corr_left  = 0.0_wp
       u_corr_right = 0.0_wp
       v_corr_south = 0.0_wp
       v_corr_north = 0.0_wp
    ELSE
!
!--    Correction equally distributed to all nest boundaries, area_total must be used as area.
!--    Note that face_area(7) is the total area (=sum from 1 to 6)
       w_corr_top   = volume_flux_integral / face_area(7)
       IF ( bc_dirichlet_l )  u_corr_left  = -w_corr_top
       IF ( bc_dirichlet_r )  u_corr_right =  w_corr_top
       IF ( bc_dirichlet_s )  v_corr_south = -w_corr_top
       IF ( bc_dirichlet_n )  v_corr_north =  w_corr_top
   ENDIF
!
!-- Correct the top-boundary value of w
    DO   i = nxl, nxr
       DO   j = nys, nyn
          DO  k = nzt, nzt + 1
             w(k,j,i) = w(k,j,i) + w_corr_top *                                                    &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
       ENDDO
    ENDDO
!
!-- Correct the bottom-boundary value of w
    DO   i = nxl, nxr
       DO   j = nys, nyn
          DO  k = nzb, nzb
             w(k,j,i) = w(k,j,i) - w_corr_top *                                                    &
                                   MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
          ENDDO
       ENDDO
    ENDDO
!
!-- Correct the left-boundary value of u
    IF ( bc_dirichlet_l )  THEN
       DO  i = nxlg, nxl
          DO  j = nys, nyn
             DO  k = nzb + 1, nzt
                u(k,j,i) = u(k,j,i) + u_corr_left                                                  &
                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Correct the right-boundary value of u
    IF ( bc_dirichlet_r )  THEN
       DO  i = nxr+1, nxrg
          DO  j = nys, nyn
             DO  k = nzb + 1, nzt
                u(k,j,i) = u(k,j,i) + u_corr_right                                                 &
                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Correct the south-boundary value of v
    IF ( bc_dirichlet_s )  THEN
       DO  i = nxl, nxr
          DO  j = nysg, nys
             DO  k = nzb + 1, nzt
                v(k,j,i) = v(k,j,i) + v_corr_south                                                 &
                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
!-- Correct the north-boundary value of v
    IF ( bc_dirichlet_n )  THEN
       DO  i = nxl, nxr
          DO  j = nyn+1, nyng
             DO  k = nzb + 1, nzt
                v(k,j,i) = v(k,j,i) + v_corr_north                                                 &
                           * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
             ENDDO
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE pmci_ensure_nest_mass_conservation

#endif
END MODULE pmc_interface
