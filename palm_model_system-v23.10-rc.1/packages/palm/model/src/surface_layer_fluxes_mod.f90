!> @file surface_layer_fluxes_mod.f90
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
!> Diagnostic computation of vertical fluxes in the constant flux layer from the values of the
!> variables at grid point k=1 based on Newton iteration.
!>
!> @todo (Re)move large_scale_forcing actions
!> @todo Check/optimize OpenMP directives
!> @todo Simplify if conditions (which flux need to be computed in which case)
!--------------------------------------------------------------------------------------------------!
 MODULE surface_layer_fluxes_mod

    USE arrays_3d,                                                                                 &
        ONLY:  d_exner,                                                                            &
               drho_air_zw,                                                                        &
               e,                                                                                  &
               nc,                                                                                 &
               nr,                                                                                 &
               pt,                                                                                 &
               q,                                                                                  &
               ql,                                                                                 &
               qc,                                                                                 &
               qr,                                                                                 &
               s,                                                                                  &
               u,                                                                                  &
               v,                                                                                  &
               vpt,                                                                                &
               w,                                                                                  &
               rho_air_zw


    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  g,                                                                                  &
               kappa,                                                                              &
               lv_d_cp,                                                                            &
               pi,                                                                                 &
               rd_d_rv

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    USE chem_modules,                                                                              &
        ONLY:  constant_csflux

    USE cpulog

    USE control_parameters,                                                                        &
        ONLY:  air_chemistry,                                                                      &
               atmosphere_run_coupled_to_ocean,                                                    &
               cloud_droplets,                                                                     &
               constant_heatflux,                                                                  &
               constant_scalarflux,                                                                &
               constant_waterflux,                                                                 &
               debug_output_timestep,                                                              &
               humidity,                                                                           &
               ibc_e_b,                                                                            &
               ibc_pt_b,                                                                           &
               indoor_model,                                                                       &
               land_surface,                                                                       &
               large_scale_forcing,                                                                &
               loop_optimization,                                                                  &
               lsf_surf,                                                                           &
               neutral,                                                                            &
               passive_scalar,                                                                     &
               pt_surface,                                                                         &
               q_surface,                                                                          &
               surface_pressure,                                                                   &
               simulated_time,                                                                     &
               time_since_reference_point,                                                         &
               urban_surface,                                                                      &
               use_free_convection_scaling

    USE kinds

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model,                                                                   &
               microphysics_morrison,                                                              &
               microphysics_seifert

    USE pegrid

    USE land_surface_model_mod,                                                                    &
        ONLY:  aero_resist_kray,                                                                   &
               skip_time_do_lsm

    USE surface_mod,                                                                               &
        ONLY:  surf_type,                                                                          &
               surf_def,                                                                           &
               surf_lsm,                                                                           &
               surf_usm


    IMPLICIT NONE

    INTEGER(iwp) ::  i      !< loop index x direction
    INTEGER(iwp) ::  i_off  !< offset index between surface and reference grid point in x direction
    INTEGER(iwp) ::  j      !< loop index y direction
    INTEGER(iwp) ::  j_off  !< offset index between surface and reference grid point in y direction
    INTEGER(iwp) ::  k      !< loop index z direction
    INTEGER(iwp) ::  k_off  !< offset index between surface and reference grid point in z direction
    INTEGER(iwp) ::  m      !< running index surface elements

    REAL(wp) :: e_s                !< Saturation water vapor pressure
    REAL(wp) :: ol_max = 1.0E6_wp  !< Maximum Obukhov length
    REAL(wp) :: z_mo               !< Height of the constant flux layer where MOST is assumed

    TYPE(surf_type), POINTER ::  surf  !< surf-type array, used to generalize subroutines


    SAVE

    PRIVATE

    PUBLIC init_surface_layer_fluxes,                                                              &
           phi_m,                                                                                  &
           psi_h,                                                                                  &
           psi_m,                                                                                  &
           surface_layer_fluxes

    INTERFACE init_surface_layer_fluxes
       MODULE PROCEDURE init_surface_layer_fluxes
    END INTERFACE init_surface_layer_fluxes

    INTERFACE phi_m
       MODULE PROCEDURE phi_m
    END INTERFACE phi_m

    INTERFACE psi_h
       MODULE PROCEDURE psi_h
    END INTERFACE psi_h

    INTERFACE psi_m
       MODULE PROCEDURE psi_m
    END INTERFACE psi_m

    INTERFACE surface_layer_fluxes
       MODULE PROCEDURE surface_layer_fluxes
    END INTERFACE surface_layer_fluxes


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Main routine to compute the surface fluxes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_layer_fluxes

    IMPLICIT NONE

    IF ( debug_output_timestep )  CALL debug_message( 'surface_layer_fluxes', 'start' )
!
!-- First, calculate the new Obukhov length from precalculated values of log(z/z0) and wind speeds.
!-- As a second step, then calculate new friction velocity, followed by the new scaling
!-- parameters (th*, q*, etc.), and the new surface fluxes, if required. Note, each routine is called
!-- for different surface types. First call for default-type surfaces, then for natural- and
!-- urban-type surfaces.
!-- Start with default-type surfaces
    IF ( surf_def%ns >= 1 )  THEN
       surf => surf_def
!
!--    First, precalculate ln(z/z0) for all surfaces. This is done each timestep, in order
!--    to account for time-dependent roughness or user-modifications.
       CALL calc_ln
!
!--    Calculate temperatures and specific humidity at the first computational grid level
!--    above the surface and store values in the 1d surface data arrays.
       CALL calc_pt_q
!
!--    Store surface temperatures and specific humidity of the Eularian fields in the
!--    1d surface data arrays.
       IF ( .NOT. neutral )  THEN
          CALL store_pt_surface
          IF ( humidity )  THEN
             CALL store_q_surface
             CALL store_vpt_surface
          ENDIF
       ENDIF
!
!--    Calculate surface-parallel absolute velocity on different positions of the staggered grid.
       CALL calc_uvw_abs
!
!--    Calculate Obukhov length.
       IF ( .NOT. neutral )  CALL calc_ol
!
!--    Calculate friction velocity representative for different positions of the staggered grid.
!--    Note, stability functions will be only employed at horizontal upward-facing surfaces.
       CALL calc_us
!
!--    Calculate scaling parameters.
       CALL calc_scaling_parameters
!
!--    Calculate surface fluxes for scalars.
       CALL calc_surface_fluxes
!
!--    Calculate surface momentum fluxes. Note, not all fluxes become effective at all surface
!--    orientations. For example, u'w'_0 is zero at vertical walls.
       CALL calc_usws
       CALL calc_vsws
       CALL calc_usvs
       CALL calc_vsus
       CALL calc_wsus_wsvs
!
!--    Calculate surface momentum fluxes on scalar grid. This is required for TKE production.
       CALL calc_usws_vsws_for_tke
    ENDIF
!
!-- Natural land surfaces.
    IF ( surf_lsm%ns >= 1 )  THEN
       surf => surf_lsm
!
!--    First, precalculate ln(z/z0) for all surfaces. This is done each timestep, in order
!--    to account for time-dependent roughness or user-modifications.
       CALL calc_ln
!
!--    Derive potential temperature and specific humidity at first grid level from the fields
!--    pt and q
       CALL calc_pt_q
!
!--    Calculate surface-parallel absolute velocity on different positions of the staggered grid.
       CALL calc_uvw_abs
!
!--    Calculate Obukhov length.
       IF ( .NOT. neutral )  CALL calc_ol
!
!--    Calculate friction velocity.
!--    Note, stability functions will be only employed at horizontal upward-facing surfaces.
       CALL calc_us
!
!--    Calculate scaling parameters.
       CALL calc_scaling_parameters
!
!--    Calculate surface fluxes for scalars.
       CALL calc_surface_fluxes
!
!--    Calculate surface momentum fluxes. Note, not all fluxes become effective at all surface
!--    orientations. For example, u'w'_0 is zero at vertical walls.
       CALL calc_usws
       CALL calc_vsws
       CALL calc_usvs
       CALL calc_vsus
       CALL calc_wsus_wsvs
!
!--    Calculate surface momentum fluxes on scalar grid. This is required for TKE production.
       CALL calc_usws_vsws_for_tke
    ENDIF
!
!-- Building surfaces.
    IF ( surf_usm%ns >= 1 )  THEN
       surf => surf_usm
!
!--    First, precalculate ln(z/z0) for all surfaces. This is done each timestep, in order
!--    to account for time-dependent roughness or user-modifications.
       CALL calc_ln
!
!--    Derive potential temperature and specific humidity at first grid level from the fields
!--    pt and q
       CALL calc_pt_q
!
!--    Calculate surface-parallel absolute velocity on different positions of the staggered grid.
       CALL calc_uvw_abs
!
!--    Calculate Obukhov length.
       IF ( .NOT. neutral )  CALL calc_ol
!
!--    Calculate friction velocity  on different positions of the staggered grid.
!--    Note, stability functions will be only employed at horizontal upward-facing surfaces.
       CALL calc_us
!
!--    Calculate scaling parameters.
       CALL calc_scaling_parameters
!
!--    Calculate surface fluxes for scalars.
       CALL calc_surface_fluxes
!
!--    Calculate surface momentum fluxes. Note, not all fluxes become effective at all surface
!--    orientations. For example, u'w'_0 is zero at vertical walls.
       CALL calc_usws
       CALL calc_vsws
       CALL calc_usvs
       CALL calc_vsus
       CALL calc_wsus_wsvs
!
!--    Calculate surface momentum fluxes on scalar grid. This is required for TKE production.
       CALL calc_usws_vsws_for_tke
    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'surface_layer_fluxes', 'end' )

 END SUBROUTINE surface_layer_fluxes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initializing actions for the surface layer routine.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_surface_layer_fluxes

    IMPLICIT NONE

    CALL location_message( 'initializing surface layer', 'start' )
!
!-- In case of runs with neutral statification, set Obukhov length to a large value
    IF ( neutral )  THEN
       IF ( surf_def%ns >= 1  .AND.  ALLOCATED( surf_def%ol ) )  surf_def%ol = 1.0E10_wp
       IF ( surf_lsm%ns >= 1  .AND.  ALLOCATED( surf_lsm%ol ) )  surf_lsm%ol = 1.0E10_wp
       IF ( surf_usm%ns >= 1  .AND.  ALLOCATED( surf_usm%ol ) )  surf_usm%ol = 1.0E10_wp
    ENDIF

    CALL location_message( 'initializing surface layer', 'finished' )

 END SUBROUTINE init_surface_layer_fluxes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute ln(z/z0).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_ln

!
!-- Note, ln(z/z0h) and ln(z/z0q) is also calculated in neutral situations.
!-- This is because the roughness for scalars are also used for other quantities such as passive
!-- scalar, chemistry and aerosols.
    !$OMP PARALLEL DO PRIVATE( z_mo )
    !$ACC PARALLEL LOOP PRIVATE(z_mo) &
    !$ACC PRESENT(surf)
    DO  m = 1, surf%ns
       z_mo = surf%z_mo(m)
       surf%ln_z_z0(m)  = LOG( z_mo / surf%z0(m) )
       surf%ln_z_z0h(m) = LOG( z_mo / surf%z0h(m) )
       surf%ln_z_z0q(m) = LOG( z_mo / surf%z0q(m) )
    ENDDO

 END SUBROUTINE calc_ln

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute the absolute value of the surface-parallel velocity (relative to the surface).
!> This is required for all surfaces.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_uvw_abs

    IMPLICIT NONE

    INTEGER(iwp) ::  ibit  !< flag to mask computation of relative velocity in case of downward-facing surfaces

    LOGICAL  ::  u_grid   !< flag indicating that absolute velocity is required on u-grid needed for north/soutward-facing surfaces

    REAL(wp) ::  u_comp   !< u-component of the velocity vector interpolated onto respective staggered grid position
    REAL(wp) ::  v_comp   !< v-component of the velocity vector interpolated onto respective staggered grid position
    REAL(wp) ::  vel_dot  !< dot product of velocity and surface normal vector
    REAL(wp) ::  w_comp   !< w-component of the velocity vector interpolated onto respective staggered grid position

    REAL(wp), DIMENSION(1:surf%ns) ::  w_lfc   !< local free convection velocity scale
    !$ACC DECLARE CREATE( w_lfc )
!
!-- Pre-calculate local free-convection scaling parameter if required. This
!-- will maintain a horizontal velocity even for very weak wind convective conditions. SIGN is
!-- used to set w_lfc to zero under stable conditions. Note, free-convection scaling velocity
!-- is only used at horizontal surfaces (multiplication with normal-vector component).
    IF ( use_free_convection_scaling )  THEN
       !$OMP PARALLEL DO
       !$ACC PARALLEL LOOP &
       !$ACC PRESENT(surf)
       DO  m = 1, surf%ns
          w_lfc(m) = ABS( g / surf%pt1(m) * surf%z_mo(m) * surf%shf(m) * surf%n_s(m,3) )
          w_lfc(m) = ( 0.5_wp * ( w_lfc(m)                                                         &
                                + SIGN( w_lfc(m) , surf%shf(m) * surf%n_s(m,3) ) ) )**(0.33333_wp) &
                     * MERGE( 1.0_wp, 0.0_wp, surf%upward(m) )
       ENDDO
    ELSE
       !$ACC PARALLEL LOOP
       DO  m = 1, surf%ns
          w_lfc(m) = 0.0_wp
       ENDDO
    ENDIF

    !$OMP PARALLEL DO PRIVATE(i, ibit, j, k, u_comp, u_grid, v_comp, w_comp, vel_dot)
    !$ACC PARALLEL LOOP PRIVATE(i, ibit, j, k, u_comp, u_grid, v_comp, w_comp, vel_dot) &
    !$ACC PRESENT(surf, u, v, w)
    DO  m = 1, surf%ns
       i   = surf%i(m)
       j   = surf%j(m)
       k   = surf%k(m)
!
!--    ibit is 1 for upward-facing surfaces, zero for downward-facing surfaces.
       ibit = MERGE( 1, 0, surf%upward(m) )
!
!--    Compute surface-parallel velocity. Therefore, first compute dot product between the
!--    the wind velocity and surface normal vector components. Note, horizontal velocity
!--    components are considered as relative velocities, which takes coupled atmosphere ocean
!--    surfaces into account (see the k-1 values). Relative velocities, however, are only
!--    considered in case of horizontal upward-facing surfaces (see ibit).
       u_comp = 0.5 * ( u(k,j,i) + u(k,j,i+1) - ( u(k-1,j,i) + u(k-1,j,i+1) ) * ibit )
       v_comp = 0.5 * ( v(k,j,i) + v(k,j+1,i) - ( v(k-1,j,i) + v(k-1,j+1,i) ) * ibit )
       w_comp = 0.5 * ( w(k,j,i) + w(k-1,j,i) )

       vel_dot = u_comp * surf%n_s(m,1) + v_comp * surf%n_s(m,2) + w_comp * surf%n_s(m,3)
       surf%uvw_abs(m) = SQRT( ( u_comp - vel_dot * surf%n_s(m,1) )**2                             &
                             + ( v_comp - vel_dot * surf%n_s(m,2) )**2                             &
                             + ( w_comp - vel_dot * surf%n_s(m,3) )**2 + w_lfc(m)**2 )
!
!--    Now compute the absolute velocity on the u- or v-grid, depending on the respective
!--    surface grid point orientation.
       u_grid = surf%joff(m) /= 0
!
!--    Compute the surface-parallel absolute velocity on u- and v-grid. Note,
!--    at north/south facing walls, the respective wall-parallel velocity actually does not
!--    include any contribution from the v-component. Even in case of straight walls this is
!--    maintained by the following code (which yields to identical results as before), even though
!--    the respective v-component on the u-grid is not necessarily zero here. However, later on,
!--    this portion cancels out by the normal vector component.
       u_comp = MERGE( u(k,j,i),                                                                   & ! u-grid
                       0.25_wp * ( u(k,j,i) + u(k,j,i+1) + u(k,j-1,i) + u(k,j-1,i+1) ),            & ! v-grid
                       u_grid )
       v_comp = MERGE( 0.25_wp * ( v(k,j,i-1) + v(k,j,i) + v(k,j+1,i-1) + v(k,j+1,i) ),            & ! u-grid
                       v(k,j,i),                                                                   & ! v-grid
                       u_grid )
       w_comp = MERGE( 0.25_wp * ( w(k-1,j,i-1) + w(k-1,j,i) + w(k,j,i-1) + w(k,j,i) ),            & ! u-grid
                       0.25_wp * ( w(k-1,j-1,i) + w(k-1,j,i) + w(k,j-1,i) + w(k,j,i) ),            & ! v-grid
                       u_grid )

       vel_dot = u_comp * surf%n_s(m,1) + v_comp * surf%n_s(m,2) + w_comp * surf%n_s(m,3)
       surf%uvw_abs_uv(m) = SQRT( ( u_comp - vel_dot * surf%n_s(m,1) )**2                          &
                                + ( v_comp - vel_dot * surf%n_s(m,2) )**2                          &
                                + ( w_comp - vel_dot * surf%n_s(m,3) )**2 )
!
!--    Now compute the surface-parallel absolute velocity on w-grid.
       u_comp = 0.25_wp * ( u(k+1,j,i+1) + u(k+1,j,i) + u(k,j,i+1) + u(k,j,i) )
       v_comp = 0.25_wp * ( v(k+1,j+1,i) + v(k+1,j,i) + v(k,j+1,i) + v(k,j,i) )
       w_comp = w(k,j,i)

       vel_dot = u_comp * surf%n_s(m,1) + v_comp * surf%n_s(m,2) + w_comp * surf%n_s(m,3)
       surf%uvw_abs_w(m) = SQRT( ( u_comp - vel_dot * surf%n_s(m,1) )**2                           &
                               + ( v_comp - vel_dot * surf%n_s(m,2) )**2                           &
                               + ( w_comp - vel_dot * surf%n_s(m,3) )**2 )
    ENDDO

 END SUBROUTINE calc_uvw_abs


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the Obukhov length (L) and Richardson flux number (z/L)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_ol

    IMPLICIT NONE

    INTEGER(iwp) ::  iter  !< Newton iteration step

    LOGICAL, DIMENSION(surf%ns) ::  convergence_reached  !< convergence switch for vectorization
    !$ACC DECLARE CREATE( convergence_reached )

    REAL(wp) ::  f       !< Function for Newton iteration: f = Ri - [...]/[...]^2 = 0
    REAL(wp) ::  f_d_ol  !< Derivative of f
    REAL(wp) ::  ol_l    !< Lower bound of L for Newton iteration
    REAL(wp) ::  ol_m    !< Previous value of L for Newton iteration
    REAL(wp) ::  ol_old  !< Previous time step value of L
    REAL(wp) ::  ol_u    !< Upper bound of L for Newton iteration

    REAL(wp), DIMENSION(surf%ns) ::  ol_old_vec  !< temporary array required for vectorization
    REAL(wp), DIMENSION(surf%ns) ::  z_mo_vec    !< temporary array required for vectorization
    !$ACC DECLARE CREATE( ol_old_vec, z_mo_vec )

!
!-- Evaluate bulk Richardson number (calculation depends on definition based on setting of boundary
!-- conditions)
    IF ( ibc_pt_b /= 1 )  THEN
       IF ( humidity )  THEN
          !$OMP PARALLEL DO PRIVATE( z_mo )
          DO  m = 1, surf%ns
             z_mo = surf%z_mo(m)
             surf%rib(m) = g * z_mo * ( surf%vpt1(m) - surf%vpt_surface(m) ) /                     &
                           ( surf%uvw_abs(m)**2 * surf%vpt1(m) + 1.0E-20_wp )
          ENDDO
       ELSE
          !$OMP PARALLEL DO PRIVATE( z_mo )
          DO  m = 1, surf%ns
             z_mo = surf%z_mo(m)
             surf%rib(m) = g * z_mo * ( surf%pt1(m) - surf%pt_surface(m) ) /                       &
                           ( surf%uvw_abs(m)**2 * surf%pt1(m) + 1.0E-20_wp )
          ENDDO
       ENDIF
    ELSE
       IF ( humidity )  THEN
          !$OMP PARALLEL DO PRIVATE( k, k_off, z_mo )
          DO  m = 1, surf%ns
             k = surf%k(m)
             k_off = surf%koff(m)
             z_mo = surf%z_mo(m)
             surf%rib(m) = - g * z_mo * ( ( 1.0_wp + 0.61_wp * surf%qv1(m) ) * surf%shf(m)         &
                                          + 0.61_wp * surf%pt1(m) * surf%qsws(m) )                 &
                           / ( surf%uvw_abs(m)**3 * surf%vpt1(m) * kappa**2 + 1.0E-20_wp )
!
!--          Note, at upward or downward-facing surfaces the sensible and latent fluxes include
!--          density. To make Rib dimensionless, multiply with 1/density.
             surf%rib(m) = MERGE( surf%rib(m) * drho_air_zw(k+k_off), surf%rib(m),                 &
                                  surf%upward(m)  .OR.  surf%downward(m) )
          ENDDO
       ELSE
          !$OMP PARALLEL DO PRIVATE( k, k_off, z_mo )
          !$ACC PARALLEL LOOP PRIVATE(k, k_off, z_mo) &
          !$ACC PRESENT(surf, drho_air_zw)
          DO  m = 1, surf%ns
             k = surf%k(m)
             k_off = surf%koff(m)
             z_mo = surf%z_mo(m)
             surf%rib(m) = - g * z_mo * surf%shf(m)                                                &
                           / ( surf%uvw_abs(m)**3 * surf%pt1(m) * kappa**2 + 1.0E-20_wp )
!
!--          Note, at upward or downward-facing surfaces the sensible flux includes density.
!--          To make Rib dimensionless, multiply with 1/density.
             surf%rib(m) = MERGE( surf%rib(m) * drho_air_zw(k+k_off), surf%rib(m),                 &
                                  surf%upward(m)  .OR.  surf%downward(m) )
          ENDDO
       ENDIF
    ENDIF

    IF ( loop_optimization == 'cache' )  THEN
!
!--    Calculate the Obukhov length using Newton iteration
       !$OMP PARALLEL DO PRIVATE(i, j, z_mo, ol_old, iter, ol_m, ol_l, ol_u, f, f_d_ol)
       !$ACC PARALLEL LOOP PRIVATE(i, j, z_mo) &
       !$ACC PRIVATE(ol_old, ol_m, ol_l, ol_u, f, f_d_ol) &
       !$ACC PRESENT(surf)
       DO  m = 1, surf%ns
          i   = surf%i(m)
          j   = surf%j(m)

          z_mo = surf%z_mo(m)

!
!--       Store current value in case the Newton iteration fails
          ol_old = surf%ol(m)

!
!--       Ensure that the bulk Richardson number and the Obukhov length have the same sign
          IF ( surf%rib(m) * surf%ol(m) < 0.0_wp  .OR.  ABS( surf%ol(m) ) == ol_max )  THEN
             IF ( surf%rib(m) > 1.0_wp ) surf%ol(m) =  0.01_wp
             IF ( surf%rib(m) < 0.0_wp ) surf%ol(m) = -0.01_wp
          ENDIF
!
!--       Iteration to find Obukhov length
          iter = 0
          DO
             iter = iter + 1
!
!--          In case of divergence, use the value of the previous time step
             IF ( iter > 1000 )  THEN
                surf%ol(m) = ol_old
                EXIT
             ENDIF

             ol_m = surf%ol(m)
             ol_l = ol_m - 0.001_wp * ol_m
             ol_u = ol_m + 0.001_wp * ol_m


             IF ( ibc_pt_b /= 1 )  THEN
!
!--             Calculate f = Ri - [...]/[...]^2 = 0
                f = surf%rib(m) - ( z_mo / ol_m ) * ( surf%ln_z_z0h(m)                             &
                                                      - psi_h( z_mo / ol_m )                       &
                                                      + psi_h( surf%z0h(m) / ol_m ) ) /            &
                    ( surf%ln_z_z0(m) - psi_m( z_mo / ol_m )                                       &
                                      + psi_m( surf%z0(m) /  ol_m ) )**2

!
!--             Calculate df/dL
                f_d_ol = ( - ( z_mo / ol_u ) * ( surf%ln_z_z0h(m)                                  &
                                                 - psi_h( z_mo / ol_u )                            &
                                                 + psi_h( surf%z0h(m) / ol_u ) ) /                 &
                           ( surf%ln_z_z0(m) - psi_m( z_mo / ol_u )                                &
                                             + psi_m( surf%z0(m) / ol_u ) )**2                     &
                           + ( z_mo / ol_l ) * ( surf%ln_z_z0h(m) - psi_h( z_mo / ol_l )           &
                                                                  + psi_h( surf%z0h(m) / ol_l ) ) /&
                           ( surf%ln_z_z0(m) - psi_m( z_mo / ol_l )                                &
                                             + psi_m( surf%z0(m) / ol_l ) )**2 ) / ( ol_u - ol_l )
             ELSE
!
!--             Calculate f = Ri - 1 /[...]^3 = 0
                f = surf%rib(m) - ( z_mo / ol_m ) /                                                &
                    ( surf%ln_z_z0(m) - psi_m( z_mo / ol_m ) + psi_m( surf%z0(m) / ol_m ) )**3

!
!--             Calculate df/dL
                f_d_ol = ( - ( z_mo / ol_u ) / ( surf%ln_z_z0(m)                                   &
                                                 - psi_m( z_mo / ol_u )                            &
                                                 + psi_m( surf%z0(m) / ol_u )                      &
                                               )**3                                                &
                           + ( z_mo / ol_l ) / ( surf%ln_z_z0(m)                                   &
                                                 - psi_m( z_mo / ol_l )                            &
                                                 + psi_m( surf%z0(m) / ol_l )                      &
                                               )**3                                                &
                          ) / ( ol_u - ol_l )
             ENDIF
!
!--          Calculate new L
             surf%ol(m) = ol_m - f / f_d_ol

!
!--          Ensure that the bulk Richardson number and the Obukhov length have the same sign and
!--          ensure convergence.
             IF ( surf%ol(m) * ol_m < 0.0_wp )  surf%ol(m) = ol_m * 0.5_wp
!
!--          If unrealistic value occurs, set L to the maximum value that is allowed
             IF ( ABS( surf%ol(m) ) > ol_max )  THEN
                surf%ol(m) = ol_max
                EXIT
             ENDIF
!
!--          Assure that Obukhov length does not become zero. If the limit is reached, exit the loop.
             IF ( ABS( surf%ol(m) ) < 1E-5_wp )  THEN
                surf%ol(m) = SIGN( 1E-5_wp, surf%ol(m) )
                EXIT
             ENDIF
!
!--          Check for convergence
             IF ( ABS( ( surf%ol(m) - ol_m ) /  surf%ol(m) ) < 1.0E-4_wp )  EXIT
          ENDDO
       ENDDO

!
!-- Vector Version
    ELSE
!
!--    Calculate the Obukhov length using Newton iteration
!--    First set arrays required for vectorization
       !$ACC PARALLEL LOOP &
       !$ACC PRESENT(surf)
       DO  m = 1, surf%ns
          z_mo_vec(m) = surf%z_mo(m)
!
!--       Store current value in case the Newton iteration fails
          ol_old_vec(m) = surf%ol(m)
!
!--       Ensure that the bulk Richardson number and the Obukhov length have the same sign
          IF ( surf%rib(m) * surf%ol(m) < 0.0_wp  .OR.  ABS( surf%ol(m) ) == ol_max )  THEN
             IF ( surf%rib(m) > 1.0_wp ) surf%ol(m) =  0.01_wp
             IF ( surf%rib(m) < 0.0_wp ) surf%ol(m) = -0.01_wp
          ENDIF
!
!--       Initialize convergence flag
          convergence_reached(m) = .FALSE.
       ENDDO

!
!--    Iteration to find Obukhov length
       iter = 0
       DO
          iter = iter + 1
!
!--       In case of divergence, use the value(s) of the previous time step
          IF ( iter > 1000 )  THEN
             !$ACC PARALLEL LOOP &
             !$ACC PRESENT(surf)
             DO  m = 1, surf%ns
                IF ( .NOT. convergence_reached(m) )  surf%ol(m) = ol_old_vec(m)
             ENDDO
             EXIT
          ENDIF

          !$ACC PARALLEL LOOP PRIVATE(ol_m, ol_l, ol_u, f, f_d_ol) &
          !$ACC PRESENT(surf)
          DO  m = 1, surf%ns
             IF ( convergence_reached(m) )  CYCLE

             ol_m = surf%ol(m)
             ol_l = ol_m - 0.001_wp * ol_m
             ol_u = ol_m + 0.001_wp * ol_m


             IF ( ibc_pt_b /= 1 )  THEN
!
!--             Calculate f = Ri - [...]/[...]^2 = 0
                f = surf%rib(m) - ( z_mo_vec(m) / ol_m ) * ( surf%ln_z_z0h(m)                      &
                                                           - psi_h( z_mo_vec(m) / ol_m )           &
                                                           + psi_h( surf%z0h(m) / ol_m )           &
                                                           ) /                                     &
                                                           ( surf%ln_z_z0(m)                       &
                                                          - psi_m( z_mo_vec(m) / ol_m )            &
                                                          + psi_m( surf%z0(m) /  ol_m )            &
                                                           )**2

!
!--             Calculate df/dL
                f_d_ol = ( - ( z_mo_vec(m) / ol_u ) * ( surf%ln_z_z0h(m)                           &
                                                        - psi_h( z_mo_vec(m) / ol_u )              &
                                                        + psi_h( surf%z0h(m) / ol_u )              &
                                                      ) /                                          &
                                                      ( surf%ln_z_z0(m)                            &
                                                        - psi_m( z_mo_vec(m) / ol_u )              &
                                                        + psi_m( surf%z0(m) / ol_u )               &
                                                      )**2                                         &
                           + ( z_mo_vec(m) / ol_l ) * ( surf%ln_z_z0h(m)                           &
                                                        - psi_h( z_mo_vec(m) / ol_l )              &
                                                        + psi_h( surf%z0h(m) / ol_l )              &
                                                      ) /                                          &
                                                      ( surf%ln_z_z0(m)                            &
                                                        - psi_m( z_mo_vec(m) / ol_l )              &
                                                        + psi_m( surf%z0(m) / ol_l )               &
                                                      )**2                                         &
                         ) / ( ol_u - ol_l )
             ELSE
!
!--             Calculate f = Ri - 1 /[...]^3 = 0
                f = surf%rib(m) - ( z_mo_vec(m) / ol_m ) / ( surf%ln_z_z0(m)                       &
                                                             - psi_m( z_mo_vec(m) / ol_m )         &
                                                             + psi_m( surf%z0(m) / ol_m )          &
                                                           )**3

!
!--             Calculate df/dL
                f_d_ol = ( - ( z_mo_vec(m) / ol_u ) / ( surf%ln_z_z0(m)                            &
                                                        - psi_m( z_mo_vec(m) / ol_u )              &
                                                        + psi_m( surf%z0(m) / ol_u )               &
                                                      )**3                                         &
                           + ( z_mo_vec(m) / ol_l ) / ( surf%ln_z_z0(m)                            &
                                                        - psi_m( z_mo_vec(m) / ol_l )              &
                                                        + psi_m( surf%z0(m) / ol_l )               &
                                                      )**3                                         &
                         ) / ( ol_u - ol_l )
             ENDIF
!
!--          Calculate new L
             surf%ol(m) = ol_m - f / f_d_ol

!
!--          Ensure that the bulk Richardson number and the Obukhov length have the same sign and
!--          ensure convergence.
             IF ( surf%ol(m) * ol_m < 0.0_wp )  surf%ol(m) = ol_m * 0.5_wp

!
!--          Check for convergence
!--          This check does not modify surf%ol, therefore this is done first
             IF ( ABS( ( surf%ol(m) - ol_m ) /  surf%ol(m) ) < 1.0E-4_wp )  THEN
                convergence_reached(m) = .TRUE.
             ENDIF
!
!--          If unrealistic value occurs, set L to the maximum allowed value
             IF ( ABS( surf%ol(m) ) > ol_max )  THEN
                surf%ol(m) = ol_max
                convergence_reached(m) = .TRUE.
             ENDIF
          ENDDO
!
!--       Assure that Obukhov length does not become zero
          !$ACC PARALLEL LOOP &
          !$ACC PRESENT(surf)
          DO  m = 1, surf%ns
             IF ( convergence_reached(m) )  CYCLE
             IF ( ABS( surf%ol(m) ) < 1E-5_wp )  THEN
                surf%ol(m) = SIGN( 10E-6_wp, surf%ol(m) )
                convergence_reached(m) = .TRUE.
             ENDIF
          ENDDO

          IF ( ALL( convergence_reached ) )  EXIT

       ENDDO  ! End of iteration loop

    ENDIF  ! End of vector branch

 END SUBROUTINE calc_ol


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate friction velocity u* representative for grid center.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_us

    IMPLICIT NONE

    !$OMP PARALLEL  DO PRIVATE( z_mo, k_off )
    !$ACC PARALLEL LOOP PRIVATE( z_mo, k_off ) &
    !$ACC PRESENT(surf)
    DO  m = 1, surf%ns
       k_off = surf%koff(m)
       z_mo = surf%z_mo(m)
!
!--    Compute u* at the scalars' grid points. At horizonally upward-facing surfaces use
!--    stability correction, else take purely neutral solution.
       surf%us(m) = MERGE( kappa * surf%uvw_abs(m) / ( surf%ln_z_z0(m)                             &
                           - psi_m( z_mo / surf%ol(m) ) + psi_m( surf%z0(m) / surf%ol(m) ) ),      &
                           kappa * surf%uvw_abs(m) / surf%ln_z_z0(m),                              &
                           surf%upward(m) )
!
!--    Compute u* on u-/v-grid. No stability corrections are employed here. This is because
!--    these are always vertical walls where stability corrections are not valid.
       surf%us_uvgrid(m) = kappa * surf%uvw_abs_uv(m) / surf%ln_z_z0(m)
!
!--    Compute u* on w-grid.
       surf%us_wgrid(m) = kappa * surf%uvw_abs_w(m) / surf%ln_z_z0(m)

    ENDDO

 END SUBROUTINE calc_us

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate potential temperature, specific humidity, and virtual potential temperature at first
!> grid level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_pt_q

    IMPLICIT NONE
!
!-- @todo: Following loop need to be split-up.

    !$OMP PARALLEL DO PRIVATE( i, j, k )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k) &
    !$ACC PRESENT(surf, pt)
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k = surf%k(m)

#ifndef _OPENACC
       IF ( bulk_cloud_model ) THEN
          surf%pt1(m) = pt(k,j,i) + lv_d_cp * d_exner(k) * ql(k,j,i)
          surf%qv1(m) = q(k,j,i) - ql(k,j,i)
       ELSEIF( cloud_droplets ) THEN
          surf%pt1(m) = pt(k,j,i) + lv_d_cp * d_exner(k) * ql(k,j,i)
          surf%qv1(m) = q(k,j,i)
       ELSE
#endif
          surf%pt1(m) = pt(k,j,i)
#ifndef _OPENACC
          IF ( humidity )  THEN
             surf%qv1(m) = q(k,j,i)
          ELSE
#endif
             surf%qv1(m) = 0.0_wp
#ifndef _OPENACC
          ENDIF
       ENDIF

       IF ( humidity )  THEN
          surf%vpt1(m) = pt(k,j,i) * ( 1.0_wp + 0.61_wp * q(k,j,i) )
       ENDIF
#endif
    ENDDO

 END SUBROUTINE calc_pt_q


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store potential temperature at surface grid level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE store_pt_surface

    IMPLICIT NONE

    !$OMP PARALLEL DO PRIVATE( i, j, k )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k ) &
    !$ACC PRESENT(surf, pt)
    DO  m = 1, surf%ns
       i = surf%i(m) + surf%ioff(m)
       j = surf%j(m) + surf%joff(m)
       k = surf%k(m) + surf%koff(m)
       surf%pt_surface(m) = pt(k,j,i)
    ENDDO

 END SUBROUTINE store_pt_surface


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store mixing ratio at surface grid level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE store_q_surface

    IMPLICIT NONE

    !$OMP PARALLEL DO PRIVATE( i, j, k )
    DO  m = 1, surf%ns
       i = surf%i(m) + surf%ioff(m)
       j = surf%j(m) + surf%joff(m)
       k = surf%k(m) + surf%koff(m)
       surf%q_surface(m) = q(k,j,i)
    ENDDO

 END SUBROUTINE store_q_surface


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Store virtual potential temperature at surface grid level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE store_vpt_surface

    IMPLICIT NONE

    !$OMP PARALLEL DO PRIVATE( i, j, k )
    DO  m = 1, surf%ns
       i = surf%i(m) + surf%ioff(m)
       j = surf%j(m) + surf%joff(m)
       k = surf%k(m) + surf%koff(m)
       surf%vpt_surface(m) = vpt(k,j,i)
    ENDDO

 END SUBROUTINE store_vpt_surface


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate the other MOST scaling parameters theta*, q*, (qc*, qr*, nc*, nr*)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_scaling_parameters

    IMPLICIT NONE

    INTEGER(iwp)  ::  lsp   !< running index for chemical species

    LOGICAL  ::  lsm_switch !<

!
!-- Compute theta* at horizontal surfaces
    IF ( constant_heatflux )  THEN
!
!--    For a given heat flux in the surface layer:

       !$OMP PARALLEL DO PRIVATE( k, k_off )
       !$ACC PARALLEL LOOP PRIVATE( k, k_off ) &
       !$ACC PRESENT( surf, drho_air_zw)
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
!
!--       Compute ts for horizontally upward-facing surfaces.
          surf%ts(m) = MERGE( -surf%shf(m) * drho_air_zw(k+k_off) / ( surf%us(m) + 1E-30_wp ),     &
                              surf%ts(m),                                                          &
                              surf%upward(m) )
!
!--       ts must be limited, because otherwise overflow may occur in case of us=0 when computing
!--       ol further below
          IF ( surf%ts(m) < -1.05E5_wp )  surf%ts(m) = -1.0E5_wp
          IF ( surf%ts(m) >  1.0E5_wp  )  surf%ts(m) =  1.0E5_wp
       ENDDO

    ELSE
!
!--    For a given surface temperature:
       IF ( large_scale_forcing  .AND.  lsf_surf )  THEN

          !$OMP PARALLEL DO PRIVATE( i, i_off, j, j_off, k, k_off )
          DO  m = 1, surf%ns
             i   = surf%i(m)
             j   = surf%j(m)
             k   = surf%k(m)
             i_off = surf%ioff(m)
             j_off = surf%joff(m)
             k_off = surf%koff(m)
! MS: This need to be changed to surf%pt_surface later
             pt(k+k_off,j+j_off,i+i_off) = pt_surface
          ENDDO
       ENDIF

       !$OMP PARALLEL DO PRIVATE( z_mo )
       DO  m = 1, surf%ns
          z_mo = surf%z_mo(m)
!
!--       Compute ts for horizontally upward-facing surfaces.
          surf%ts(m) = MERGE( kappa * ( surf%pt1(m) - surf%pt_surface(m) )                         &
                              / ( surf%ln_z_z0h(m) - psi_h( z_mo / surf%ol(m) )                    &
                                                   + psi_h( surf%z0h(m) / surf%ol(m) ) ),          &
                              surf%ts(m),                                                          &
                              surf%upward(m) )
       ENDDO

    ENDIF
!
!-- Compute theta* again at vertical surfaces. This is only required for natural surfaces when
!-- aerodynamical resistance are computed via MOST relations.
    lsm_switch = land_surface  .AND.  .NOT. aero_resist_kray  .AND.                                &
                 ( ALLOCATED( surf%pavement_surface )  .OR.                                        &
                   ALLOCATED( surf%vegetation_surface )  .OR.                                      &
                   ALLOCATED( surf%water_surface ) )
    !$OMP PARALLEL DO
    !$ACC PARALLEL LOOP &
    !$ACC PRESENT( surf)
    DO  m = 1, surf%ns
!
!--    Save already computed values at horizontal surfaces.
       surf%ts(m) = MERGE( -surf%shf(m) / ( surf%us(m) + 1E-30_wp ),                               &
                           surf%ts(m),                                                             &
                           .NOT. ( surf%upward(m)  .OR.  surf%downward(m) )  .AND.  lsm_switch )
!
!--    ts must be limited, because otherwise overflow may occur in case of us=0 when computing ol
!--    further below
       IF ( surf%ts(m) < -1.05E5_wp )  surf%ts(m) = -1.0E5_wp
       IF ( surf%ts(m) >  1.0E5_wp  )  surf%ts(m) =  1.0E5_wp
    ENDDO

!
!-- If required compute q* at horizontal surfaces
    IF ( humidity )  THEN
       IF ( constant_waterflux )  THEN
!
!--       For a given water flux in the surface layer
          !$OMP PARALLEL DO PRIVATE( i, j, k, k_off )
          DO  m = 1, surf%ns
             i = surf%i(m)
             j = surf%j(m)
             k = surf%k(m)
             k_off = surf%koff(m)
             surf%qs(m) = MERGE( -surf%qsws(m) * drho_air_zw(k+k_off) / ( surf%us(m) + 1E-30_wp ), &
                                 surf%qs(m),                                                       &
                                 surf%upward(m) )
          ENDDO

       ELSE

          IF ( large_scale_forcing  .AND.  lsf_surf )  THEN
             !$OMP PARALLEL DO PRIVATE( i, i_off, j, j_off, k, k_off )
             DO  m = 1, surf%ns
                i = surf%i(m)
                j = surf%j(m)
                k = surf%k(m)
                i_off = surf%ioff(m)
                j_off = surf%joff(m)
                k_off = surf%koff(m)
                q(k+k_off,j+j_off,i+i_off) = q_surface
             ENDDO
          ENDIF

!
!--       Assume saturation for atmosphere coupled to ocean (but not in case of precursor runs).
          IF ( atmosphere_run_coupled_to_ocean )  THEN
             !$OMP PARALLEL DO PRIVATE( i, j, k, e_s )
             DO  m = 1, surf%ns
                i   = surf%i(m)
                j   = surf%j(m)
                k   = surf%k(m)
                e_s = 6.1_wp * EXP( 0.07_wp * ( MIN( pt(k-1,j,i), pt(k,j,i) ) - 273.15_wp ) )
                q(k-1,j,i) = rd_d_rv * e_s / ( surface_pressure - e_s )
             ENDDO
          ENDIF

          !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, z_mo )
          DO  m = 1, surf%ns
             i = surf%i(m)
             j = surf%j(m)
             k = surf%k(m)
             k_off = surf%koff(m)
             z_mo = surf%z_mo(m)
             surf%qs(m) = MERGE( kappa * ( surf%qv1(m) - surf%q_surface(m) )                       &
                                / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                  &
                                                     + psi_h( surf%z0q(m) / surf%ol(m) ) ),        &
                                 surf%qs(m),                                                       &
                                 surf%upward(m) )
          ENDDO
       ENDIF
!
!--    Compute q* at vertical surfaces
       !$OMP PARALLEL DO
       DO  m = 1, surf%ns
!
!--       Save already computed values at horizontal surfaces.
          surf%qs(m) = MERGE( -surf%qsws(m) / ( surf%us(m) + 1E-30_wp ),                           &
                              surf%qs(m),                                                          &
                              .NOT. ( surf%upward(m)  .OR.  surf%downward(m) ) )
       ENDDO
    ENDIF

!
!-- If required compute s*
    IF ( passive_scalar )  THEN
!
!--    At horizontal surfaces
       IF ( constant_scalarflux  )  THEN
!
!--       For a given scalar flux in the surface layer
          !$OMP PARALLEL DO
          DO  m = 1, surf%ns
             surf%ss(m) = MERGE( -surf%ssws(m) / ( surf%us(m) + 1E-30_wp ), surf%ss(m), surf%upward(m) )
          ENDDO
       ELSE

          !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, z_mo )
          DO  m = 1, surf%ns
             i = surf%i(m)
             j = surf%j(m)
             k = surf%k(m)
             k_off = surf%koff(m)
             z_mo = surf%z_mo(m)

             surf%ss(m) = MERGE( kappa * ( s(k,j,i) - s(k+k_off,j,i) )                             &
                                 / ( surf%ln_z_z0h(m) - psi_h( z_mo / surf%ol(m) )                 &
                                               + psi_h( surf%z0h(m) / surf%ol(m) ) ),              &
                                 surf%ss(m),                                                       &
                                 surf%upward(m) )
          ENDDO
       ENDIF
!
!--    At vertical surfaces
       !$OMP PARALLEL DO
       DO  m = 1, surf%ns
!
!--       Save already computed values at horizontal surfaces.
          surf%ss(m) = MERGE( -surf%ssws(m) / ( surf%us(m) + 1E-30_wp ),                           &
                              surf%ss(m),                                                          &
                              .NOT. ( surf%upward(m)  .OR.  surf%downward(m) ) )
       ENDDO
    ENDIF

!
!-- If required compute cs* (chemical species)
    IF ( air_chemistry  )  THEN
!
!--    At horizontal surfaces
       DO  lsp = 1, nvar
          IF ( constant_csflux(lsp) )  THEN
!--          For a given chemical species' flux in the surface layer
             !$OMP PARALLEL DO
             DO  m = 1, surf%ns
                surf%css(lsp,m) = MERGE( -surf%cssws(lsp,m) / ( surf%us(m) + 1E-30_wp ),           &
                                         surf%css(lsp,m),                                          &
                                         surf%upward(m) )
             ENDDO
          ENDIF
       ENDDO
!
!--    At vertical surfaces
       DO  lsp = 1, nvar
          !$OMP PARALLEL DO
          DO  m = 1, surf%ns
!
!--          Save already computed values at horizontal surfaces.
             surf%css(lsp,m) = MERGE( -surf%cssws(lsp,m) / ( surf%us(m) + 1E-30_wp ),              &
                                      surf%css(lsp,m),                                             &
                                      .NOT. ( surf%upward(m)  .OR.  surf%downward(m) ) )
          ENDDO
       ENDDO
    ENDIF

!
!-- If required compute qc* and nc*
    IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, z_mo )
       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)
          k_off = surf%koff(m)

          z_mo = surf%z_mo(m)

          surf%qcs(m) = MERGE( kappa * ( qc(k,j,i) - qc(k+k_off,j,i) )                             &
                               / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                   &
                                                    + psi_h( surf%z0q(m) / surf%ol(m) ) ),         &
                               surf%qcs(m),                                                        &
                               surf%upward(m) )

          surf%ncs(m) = MERGE( kappa * ( nc(k,j,i) - nc(k+k_off,j,i) )                             &
                               / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                   &
                                                    + psi_h( surf%z0q(m) / surf%ol(m) ) ),         &
                               surf%ncs(m),                                                        &
                               surf%upward(m) )
       ENDDO

    ENDIF

!
!-- If required compute qr* and nr*
    IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, z_mo )
       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)
          k_off = surf%koff(m)

          z_mo = surf%z_mo(m)

          surf%qrs(m) = MERGE( kappa * ( qr(k,j,i) - qr(k+k_off,j,i) )                             &
                               / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                   &
                                                    + psi_h( surf%z0q(m) / surf%ol(m) ) ),         &
                               surf%qrs(m),                                                        &
                               surf%upward(m) )

          surf%nrs(m) = MERGE( kappa * ( nr(k,j,i) - nr(k+k_off,j,i) )                             &
                               / ( surf%ln_z_z0q(m) - psi_h( z_mo / surf%ol(m) )                   &
                                                    + psi_h( surf%z0q(m) / surf%ol(m) ) ),         &
                               surf%nrs(m),                                                        &
                               surf%upward(m) )
       ENDDO

    ENDIF

 END SUBROUTINE calc_scaling_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux usws which later enters the u-equation. Note, usws is actually only
!> required at horizontal surfaces, though the loops runs over all surfaces.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_usws

    IMPLICIT NONE

    REAL(wp)      ::  denom   !< denominator including ln(z/z0) + integrated profile functions to account for stability
    REAL(wp)      ::  u_comp  !< u-component


!
!-- Compute u'w'
    !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, z_mo, u_comp, denom )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k, k_off, z_mo, u_comp, denom ) &
    !$ACC PRESENT(surf, u, rho_air_zw)
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k = surf%k(m)
       k_off = surf%koff(m)
!
!--    Depending on the orientation of the surface the resulting u-component is either defined
!--    as a relative velocity to the surface (in case of coupled atmosphere-ocean runs), or a
!--    absolute value (at downward-facing surfaces). The same is done for the stability
!--    correction which is only considered at upward-facing walls.
       z_mo    = surf%z_mo(m)
       u_comp  = MERGE( u(k,j,i) - u(k-1,j,i), u(k,j,i), surf%upward(m) )
       denom   = MERGE( surf%ln_z_z0(m) - psi_m( z_mo / surf%ol(m) )                               &
                                        + psi_m( surf%z0(m) / surf%ol(m) ),                        &
                        surf%ln_z_z0(m),                                                           &
                        surf%upward(m) )
!
!--    Please note, the computation of usws is not fully accurate. Actually a further
!--    interpolation of ol onto the u-grid, where usws is defined, is required. However, this
!--    is not done here since this would require several data transfers between the surface-data
!--    structures. Moreover, please note, usws is calculated relative to the surface speed, which
!--    is necessary for coupled atmosphere-ocean runs. In case of no-slip conditions this does not
!--    harm. Further, to account for different facings (up/downward), multiply with the respective
!--    normal-vector component.
       surf%usws(m) = kappa * u_comp / denom
       surf%usws(m) = -surf%usws(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_eff(m)
    ENDDO
!
!-- Mask usws at vertically bounded grid points.
    !$ACC PARALLEL LOOP &
    !$ACC PRESENT(surf)
    DO  m = 1, surf%ns
       surf%usws(m) = MERGE( surf%usws(m), 0.0_wp, surf%upward(m)  .OR.  surf%downward(m) )
    ENDDO

 END SUBROUTINE calc_usws


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux vsws which later enters the v-equation. Note, vsws is actually only
!> required at horizontal surfaces, though the loops runs over all surfaces.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_vsws

    IMPLICIT NONE

    REAL(wp)      ::  denom   !< denominator including ln(z/z0) + integrated profile functions to account for stability
    REAL(wp)      ::  v_comp  !< v-component


!
!-- Compute v'w'
    !$OMP PARALLEL DO PRIVATE( i, j, k, k_off, z_mo, v_comp, denom )
    !$ACC PARALLEL LOOP PRIVATE(i, j, k, k_off, z_mo, v_comp, denom ) &
    !$ACC PRESENT(surf, v, rho_air_zw)
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k = surf%k(m)
       k_off = surf%koff(m)
!
!--    Depending on the orientation of the surface the resulting v-component is either defined
!--    as a relative velocity to the surface (in case of coupled atmosphere-ocean runs), or a
!--    absolute value (at downward-facing surfaces). The same is done for the stability
!--    correction which is only considered at upward-facing walls.
       z_mo    = surf%z_mo(m)
       v_comp  = MERGE( v(k,j,i) - v(k-1,j,i), v(k,j,i), surf%upward(m) )
       denom   = MERGE( surf%ln_z_z0(m) - psi_m( z_mo / surf%ol(m) )                               &
                                        + psi_m( surf%z0(m) / surf%ol(m) ),                        &
                        surf%ln_z_z0(m),                                                           &
                        surf%upward(m) )
!
!--    Please note, the computation of vsws is not fully accurate. Actually a further
!--    interpolation of ol onto the v-grid, where vsws is defined, is required. However, this
!--    is not done here since this would require several data transfers between the surface-data
!--    structures. Moreover, please note, vsws is calculated relative to the surface speed, which
!--    is necessary for coupled atmosphere-ocean runs. In case of no-slip conditions this does not
!--    harm. Further, to account for different facings (up/downward), multiply with the respective
!--    normal-vector component.
       surf%vsws(m) = kappa * v_comp / denom
       surf%vsws(m) = -surf%vsws(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_eff(m)
    ENDDO
!
!-- Mask vsws at vertically bounded grid points.
    !$ACC PARALLEL LOOP &
    !$ACC PRESENT(surf)
    DO  m = 1, surf%ns
       surf%vsws(m) = MERGE( surf%vsws(m), 0.0_wp, surf%upward(m)  .OR.  surf%downward(m) )
    ENDDO

 END SUBROUTINE calc_vsws


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux vsus at vertical surfaces which later enter the u-equation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_usvs

    IMPLICIT NONE

    REAL(wp) ::  flag_u  !< flag indicating u-grid, used for calculation of horizontal momentum fluxes at vertical surfaces


!
!-- Generalize computation by introducing flags. At north- and south-facing surfaces
!-- u-component is used, at east- and west-facing surfaces v-component is used.
    !$OMP PARALLEL  DO PRIVATE( i, j, k, flag_u )
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k = surf%k(m)

       flag_u = MERGE( 1.0_wp, 0.0_wp, surf%northward(m)  .OR.  surf%southward(m) )
!
!--    Compute the fluxes for the horizontal transport of u at vertical surfaces.
!--    Mask fluxes at non-relevant, horizontal surface-orientations by multiplication with 0.
       surf%usvs(m) = kappa * ( flag_u * u(k,j,i) ) / surf%ln_z_z0(m)
       surf%usvs(m) = -surf%usvs(m) * surf%us_uvgrid(m) * surf%n_eff(m)
    ENDDO

 END SUBROUTINE calc_usvs


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux usvs and vsus at vertical surfaces which later enter the u- and
!> v-equation, respectively.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_vsus

    IMPLICIT NONE

    REAL(wp) ::  flag_v  !< flag indicating v-grid, used for calculation of horizontal momentum fluxes at vertical surfaces


!
!-- Generalize computation by introducing flags. At north- and south-facing surfaces
!-- u-component is used, at east- and west-facing surfaces v-component is used.
    !$OMP PARALLEL  DO PRIVATE( i, j, k, flag_v )
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k = surf%k(m)

       flag_v = MERGE( 1.0_wp, 0.0_wp, surf%eastward(m)  .OR.  surf%westward(m) )
!
!--    Compute the fluxes for the horizontal transport of v at vertical surfaces.
!--    Mask fluxes at non-relevant, horizontal surface-orientations by multiplication with 0.
       surf%vsus(m) = kappa * ( flag_v * v(k,j,i) ) / surf%ln_z_z0(m)
       surf%vsus(m) = -surf%vsus(m) * surf%us_uvgrid(m) * surf%n_eff(m)
    ENDDO

 END SUBROUTINE calc_vsus


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux wsus and wsvs at vertical surfaces which later enter the w-equation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_wsus_wsvs

    IMPLICIT NONE



    !$OMP PARALLEL DO PRIVATE( i, j, k )
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k = surf%k(m)
!
!--    Compute the fluxes for the horizontal transport of w at vertical surfaces. Mask horizontal
!--    surfaces by multiplication with 0.
       surf%wsus_wsvs(m) = kappa * w(k,j,i) / surf%ln_z_z0(m)
       surf%wsus_wsvs(m) = -surf%wsus_wsvs(m) * surf%us_wgrid(m) * surf%n_eff(m) *                 &
                           MERGE( 1.0_wp, 0.0_wp, surf%northward(m)  .OR.  surf%southward(m)  .OR. &
                                                  surf%eastward(m)   .OR.  surf%westward(m) )
    ENDDO

 END SUBROUTINE calc_wsus_wsvs


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface flux usws and vsws on scalar grid. These fluxes later enter the
!> SGS-TKE-equation.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_usws_vsws_for_tke

    IMPLICIT NONE

    REAL(wp) ::  dum     !< dummy to precalculate logarithm
    REAL(wp) ::  u_comp  !< u-component interpolated onto scalar grid point
    REAL(wp) ::  v_comp  !< v-component interpolated onto scalar grid point
    REAL(wp) ::  w_comp  !< w-component interpolated onto scalar grid point


!
!-- Compute momentum fluxes used for subgrid-scale TKE production at vertical surfaces.
!-- Fluxes required for the TKE prognostic are located at the grid center.
!-- Please note, the note runs over all surfaces but the resulting fluxes will effectively computed
!-- only at vertical surfaces.
    !$OMP PARALLEL DO PRIVATE( i, j, k, dum, u_comp, v_comp, w_comp )
    DO  m = 1, surf%ns
       i = surf%i(m)
       j = surf%j(m)
       k = surf%k(m)

       u_comp = MERGE( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ), 0.0_wp,                                 &
                       surf%northward(m)  .OR.  surf%southward(m) )
       v_comp = MERGE( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ), 0.0_wp,                                 &
                       surf%eastward(m)  .OR.  surf%westward(m) )
       w_comp = MERGE( 0.5_wp * ( w(k,j,i) + w(k-1,j,i) ), 0.0_wp,                                 &
                       .NOT. ( surf%upward(m)  .OR.  surf%downward(m) ) )

       dum = kappa / surf%ln_z_z0(m)
!
!--    usvs  at north/southward-facing walls (joff/=0) or vsus at
!--    east/westward-facing walls (ioff/=0).
       surf%mom_flux_tke(0,m) = dum * ( u_comp + v_comp )
!
!--    wsvs at north/southward-facing walls (joff/=0) or wsus at
!--    east/westward-facing walls (ioff/=0).
       surf%mom_flux_tke(1,m) = dum * w_comp
!
!--    Finally, the momentum fluxes needs to be multiplied with u*. Note, multiplication with
!--    normal vector is done directly in the calculation of the shear-production term in
!--    turbulence_closure_mod
       surf%mom_flux_tke(0:1,m) = -surf%mom_flux_tke(0:1,m) * surf%us(m)
    ENDDO

 END SUBROUTINE calc_usws_vsws_for_tke


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate surface fluxes usws, vsws, shf, qsws, (qcsws, qrsws, ncsws, nrsws)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_surface_fluxes

    IMPLICIT NONE

    INTEGER(iwp)  ::  lsp   !< running index for chemical species

    LOGICAL       ::  compute_shf    !< control flag for computation of surface sensible heat flux
    LOGICAL       ::  compute_qsws   !< control flag for computation of surface latent heat flux
    LOGICAL       ::  compute_ssws   !< control flag for computation of surface passive scalar flux
    LOGICAL       ::  compute_qcsws  !< control flag for computation of surface cloud-water content flux (is this necessary?)
    LOGICAL       ::  compute_qrsws  !< control flag for computation of surface cloud-water content flux (is this necessary?)


!
!-- Set control flags to decide whether fluxes need to be computed or not.
    compute_shf = ( .NOT.  constant_heatflux  .AND.  .NOT. neutral  .AND.                          &
                    ( ( time_since_reference_point <=  skip_time_do_lsm  .AND.                     &
                        simulated_time > 0.0_wp )  .OR.  .NOT.  land_surface )  .AND.              &
                    .NOT.  urban_surface )
    compute_qsws = ( .NOT.  constant_heatflux  .AND.                                               &
                     ( ( time_since_reference_point <=  skip_time_do_lsm  .AND.                    &
                         simulated_time > 0.0_wp )  .OR.  .NOT.  land_surface )  .AND.             &
                     .NOT.  urban_surface  .AND.  humidity )
    compute_ssws = ( .NOT.  constant_scalarflux  .AND.  passive_scalar )
    compute_qcsws = ( bulk_cloud_model  .AND.  microphysics_morrison  .AND.  humidity )
    compute_qrsws = ( bulk_cloud_model  .AND.  microphysics_seifert   .AND.  humidity )

    IF ( compute_shf )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
          surf%shf(m) = MERGE( -surf%ts(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,3),     &
                               0.0_wp,                                                             &
                               surf%upward(m) )
       ENDDO
    ENDIF
    IF ( compute_qsws )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
          surf%qsws(m) = MERGE( -surf%qs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,3),    &
                                0.0_wp,                                                            &
                                surf%upward(m) )
       ENDDO
    ENDIF
    IF ( compute_ssws )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
          surf%ssws(m) = MERGE( -surf%ss(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,3),    &
                                0.0_wp,                                                            &
                                surf%upward(m) )
       ENDDO
    ENDIF
!
!-- Compute (turbulent) fluxes of cloud water content and cloud drop conc.
    IF ( compute_qcsws )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)
!
!--       Compute (turbulent) fluxes of cloud water content and cloud drop conc.
          surf%qcsws(m) = MERGE( -surf%qcs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,3),  &
                                0.0_wp,                                                            &
                                surf%upward(m) )
          surf%ncsws(m) = MERGE( -surf%ncs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,3),  &
                                0.0_wp,                                                            &
                                surf%upward(m) )
       ENDDO
    ENDIF
    IF ( compute_qrsws )  THEN
       !$OMP PARALLEL DO PRIVATE( k, k_off )
       DO  m = 1, surf%ns
          k = surf%k(m)
          k_off = surf%koff(m)

          surf%qrsws(m) = MERGE( -surf%qrs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,3),  &
                                0.0_wp,                                                            &
                                surf%upward(m) )
          surf%nrsws(m) = MERGE( -surf%nrs(m) * surf%us(m) * rho_air_zw(k+k_off) * surf%n_s(m,3),  &
                                0.0_wp,                                                            &
                                surf%upward(m) )
       ENDDO
    ENDIF
!
!-- Compute the vertical chemical species' flux
    DO  lsp = 1, nvar
       IF (  .NOT.  constant_csflux(lsp)  .AND.  air_chemistry )  THEN
          !$OMP PARALLEL DO PRIVATE( k, k_off )
          DO  m = 1, surf%ns
             k = surf%k(m)
             k_off = surf%koff(m)
             surf%cssws(lsp,m) = MERGE( -surf%css(lsp,m) * surf%us(m) * rho_air_zw(k+k_off)        &
                                                                      * surf%n_s(m,3),             &
                                        0.0_wp,                                                    &
                                        surf%upward(m) )
          ENDDO
       ENDIF
    ENDDO

!
!-- Boundary condition for the TKE.
    IF ( ibc_e_b == 2 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j, k, i_off, j_off, k_off )
       DO  m = 1, surf%ns
          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)
          i_off = surf%ioff(m)
          j_off = surf%joff(m)
          k_off = surf%koff(m)

          e(k,j,i) = ( surf%us(m) / 0.1_wp )**2
          e(k+k_off,j+j_off,i+i_off) = e(k,j,i)
       ENDDO
    ENDIF

 END SUBROUTINE calc_surface_fluxes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Integrated stability function for momentum.
!--------------------------------------------------------------------------------------------------!
 FUNCTION psi_m( zeta )
    !$ACC ROUTINE SEQ

    USE kinds

    IMPLICIT NONE

    REAL(wp) ::  psi_m  !< Integrated similarity function result
    REAL(wp) ::  zeta   !< Stability parameter z/L
    REAL(wp) ::  x      !< dummy variable

    REAL(wp), PARAMETER ::  a = 1.0_wp            !< constant
    REAL(wp), PARAMETER ::  b = 0.66666666666_wp  !< constant
    REAL(wp), PARAMETER ::  c = 5.0_wp            !< constant
    REAL(wp), PARAMETER ::  d = 0.35_wp           !< constant
    REAL(wp), PARAMETER ::  c_d_d = c / d         !< constant
    REAL(wp), PARAMETER ::  bc_d_d = b * c / d    !< constant


    IF ( zeta < 0.0_wp )  THEN
       x = SQRT( SQRT( 1.0_wp  - 16.0_wp * zeta ) )
       psi_m = pi * 0.5_wp - 2.0_wp * ATAN( x ) + LOG( ( 1.0_wp + x )**2                           &
               * ( 1.0_wp + x**2 ) * 0.125_wp )
    ELSE

       psi_m = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - a * zeta - bc_d_d
!
!--    Old version for stable conditions (only valid for z/L < 0.5) psi_m = - 5.0_wp * zeta

    ENDIF

 END FUNCTION psi_m


!--------------------------------------------------------------------------------------------------!
! Description:
!------------
!> Integrated stability function for heat and moisture.
!--------------------------------------------------------------------------------------------------!
 FUNCTION psi_h( zeta )
    !$ACC ROUTINE SEQ

    USE kinds

    IMPLICIT NONE

    REAL(wp) ::  psi_h  !< Integrated similarity function result
    REAL(wp) ::  zeta   !< Stability parameter z/L
    REAL(wp) ::  x      !< dummy variable

    REAL(wp), PARAMETER ::  a = 1.0_wp            !< constant
    REAL(wp), PARAMETER ::  b = 0.66666666666_wp  !< constant
    REAL(wp), PARAMETER ::  c = 5.0_wp            !< constant
    REAL(wp), PARAMETER ::  d = 0.35_wp           !< constant
    REAL(wp), PARAMETER ::  c_d_d = c / d         !< constant
    REAL(wp), PARAMETER ::  bc_d_d = b * c / d    !< constant


    IF ( zeta < 0.0_wp )  THEN
       x = SQRT( 1.0_wp  - 16.0_wp * zeta )
       psi_h = 2.0_wp * LOG( (1.0_wp + x ) / 2.0_wp )
    ELSE
       psi_h = - b * ( zeta - c_d_d ) * EXP( -d * zeta ) - (1.0_wp                                 &
               + 0.66666666666_wp * a * zeta )**1.5_wp - bc_d_d + 1.0_wp
!
!--    Old version for stable conditions (only valid for z/L < 0.5)
!--    psi_h = - 5.0_wp * zeta
    ENDIF

 END FUNCTION psi_h


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates stability function for momentum
!>
!> @author Hauke Wurps
!--------------------------------------------------------------------------------------------------!
 FUNCTION phi_m( zeta )
    !$ACC ROUTINE SEQ

    IMPLICIT NONE

    REAL(wp) ::  phi_m  !< Value of the function
    REAL(wp) ::  zeta   !< Stability parameter z/L

    REAL(wp), PARAMETER ::  a = 16.0_wp  !< constant
    REAL(wp), PARAMETER ::  c = 5.0_wp   !< constant

    IF ( zeta < 0.0_wp )  THEN
       phi_m = 1.0_wp / SQRT( SQRT( 1.0_wp - a * zeta ) )
    ELSE
       phi_m = 1.0_wp + c * zeta
    ENDIF

 END FUNCTION phi_m

 END MODULE surface_layer_fluxes_mod
