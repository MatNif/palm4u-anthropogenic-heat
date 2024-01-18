!> @file surface_coupler_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> PMC utilities for atmosphere-ocean coupling
!--------------------------------------------------------------------------------------------------!
 MODULE surface_coupler_mod

   USE arrays_3d,                                                                                  &
       ONLY:  pt,                                                                                  &
              pt_p,                                                                                &
              rho_ocean,                                                                           &
              sa,                                                                                  &
              u,                                                                                   &
              u_p,                                                                                 &
              v,                                                                                   &
              v_p

   USE basic_constants_and_equations_mod,                                                          &
       ONLY:  c_p,                                                                                 &
              l_v

   USE boundary_settings_mod,                                                                      &
       ONLY:  set_lateral_neumann_bc

   USE control_parameters,                                                                         &
       ONLY:  humidity_remote,                                                                     &
              land_surface,                                                                        &
              message_string,                                                                      &
              urban_surface

   USE exchange_horiz_mod,                                                                         &
       ONLY:  exchange_horiz,                                                                      &
              exchange_horiz_2d

   USE indices,                                                                                    &
       ONLY:  nbgp,                                                                                &
              nx,                                                                                  &
              nxl,                                                                                 &
              nxlg,                                                                                &
              nxr,                                                                                 &
              nxrg,                                                                                &
              ny,                                                                                  &
              nyn,                                                                                 &
              nyng,                                                                                &
              nys,                                                                                 &
              nysg,                                                                                &
              nzt

   USE kinds

   USE pegrid

   USE surface_mod,                                                                                &
       ONLY:  surf_def,                                                                            &
              surf_top

   IMPLICIT NONE
   PRIVATE

   INTEGER(iwp), PARAMETER, PUBLIC ::  parent_send = 1               !<
   INTEGER(iwp), PARAMETER, PUBLIC ::  parent_recv = 2               !<
   INTEGER(iwp), PARAMETER, PUBLIC ::  child_send  = 3               !<
   INTEGER(iwp), PARAMETER, PUBLIC ::  child_recv  = 4               !<

   REAL(wp), PARAMETER ::  cpw = 4218.0_wp  !< heat capacity of water at constant pressure

   REAL(wp), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:,:) ::  surface_coupler_exchange_array_1  !< 2d-array for exchange data
   REAL(wp), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:,:) ::  surface_coupler_exchange_array_2  !< 2d-array for exchange data
   REAL(wp), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:,:) ::  surface_coupler_exchange_array_3  !< 2d-array for exchange data
   REAL(wp), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:,:) ::  surface_coupler_exchange_array_4  !< 2d-array for exchange data

   INTERFACE surface_coupler_alloc_mem
      MODULE PROCEDURE surface_coupler_alloc_mem
   END INTERFACE surface_coupler_alloc_mem

   INTERFACE surface_coupler_buffer_handling
      MODULE PROCEDURE surface_coupler_buffer_handling
   END INTERFACE surface_coupler_buffer_handling

   PUBLIC surface_coupler_alloc_mem, surface_coupler_buffer_handling

 CONTAINS


 SUBROUTINE surface_coupler_alloc_mem

    IMPLICIT NONE


    ALLOCATE( surface_coupler_exchange_array_1(nysg:nyng,nxlg:nxrg) )
    ALLOCATE( surface_coupler_exchange_array_2(nysg:nyng,nxlg:nxrg) )
    ALLOCATE( surface_coupler_exchange_array_3(nysg:nyng,nxlg:nxrg) )
    ALLOCATE( surface_coupler_exchange_array_4(nysg:nyng,nxlg:nxrg) )

 END SUBROUTINE surface_coupler_alloc_mem


 SUBROUTINE surface_coupler_buffer_handling( action )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN) ::  action  !<

    INTEGER(iwp) :: i  !<
    INTEGER(iwp) :: j  !<
    INTEGER(iwp) :: m  !<


    SELECT CASE ( action )

       CASE ( parent_send )

          CALL transfer_1d_to_2d( surface_coupler_exchange_array_1, surf_def%shf  )
          CALL transfer_1d_to_2d( surface_coupler_exchange_array_2, surf_def%qsws )
          CALL transfer_1d_to_2d( surface_coupler_exchange_array_3, surf_def%usws )
          CALL transfer_1d_to_2d( surface_coupler_exchange_array_4, surf_def%vsws )

       CASE ( child_recv )

          CALL transfer_2d_to_1d( surface_coupler_exchange_array_1, surf_top%shf  )
          CALL transfer_2d_to_1d( surface_coupler_exchange_array_2, surf_top%qsws )
          CALL transfer_2d_to_1d( surface_coupler_exchange_array_3, surf_top%usws )
          CALL transfer_2d_to_1d( surface_coupler_exchange_array_4, surf_top%vsws )
!
!--       Conversions of fluxes received from atmosphere.
          IF ( humidity_remote )  THEN
!
!--          Here top heat flux is still the sum of atmospheric bottom heat fluxes,
!--          * latent heat of vaporization in m2/s2, or 540 cal/g, or 40.65 kJ/mol
!--          /(rho_atm(=1.0)*c_p)
             DO  m = 1, surf_top%ns
                i = surf_top%i(m)
                j = surf_top%j(m)

                surf_top%shf(m) = surf_top%shf(m) + surf_top%qsws(m) * l_v / c_p
!
!--             ...and convert it to a salinity flux at the sea surface (top)
!--             following Steinhorn (1991), JPO 21, pp. 1681-1683:
!--             S'w' = -S * evaporation / ( rho_water * ( 1 - S ) )
                surf_top%sasws(m) = -1.0_wp * sa(nzt,j,i) * 0.001_wp * surf_top%qsws(m) /          &
                                    ( rho_ocean(nzt,j,i) * ( 1.0_wp - sa(nzt,j,i) * 0.001_wp )  )
             ENDDO
          ENDIF
!
!--       Adjust the kinematic heat flux with respect to ocean density
!--       (constants are the specific heat capacities for air and water), as well
!--       as momentum fluxes
          DO  m = 1, surf_top%ns
             i = surf_top%i(m)
             j = surf_top%j(m)
             surf_top%shf(m)  = surf_top%shf(m)  / rho_ocean(nzt,j,i) * c_p / cpw
             surf_top%usws(m) = surf_top%usws(m) / rho_ocean(nzt,j,i)
             surf_top%vsws(m) = surf_top%vsws(m) / rho_ocean(nzt,j,i)
          ENDDO

       CASE ( child_send )

          surface_coupler_exchange_array_1 = pt(nzt,:,:)
          surface_coupler_exchange_array_2 = u(nzt,:,:)
          surface_coupler_exchange_array_3 = v(nzt,:,:)

       CASE ( parent_recv )

          CALL exchange_horiz_2d( surface_coupler_exchange_array_1 )
          CALL exchange_horiz_2d( surface_coupler_exchange_array_2 )
          CALL exchange_horiz_2d( surface_coupler_exchange_array_3 )

          CALL set_lateral_neumann_bc( surface_coupler_exchange_array_1 )
          CALL set_lateral_neumann_bc( surface_coupler_exchange_array_2 )
          CALL set_lateral_neumann_bc( surface_coupler_exchange_array_3 )

!
!--       Timelevel t+dt must be set, too, since it is used in the Runge-Kutta sub-steps due to the
!--       swapping of time levels.
          pt(0,:,:)   = surface_coupler_exchange_array_1
          pt_p(0,:,:) = surface_coupler_exchange_array_1
          u(0,:,:)   = surface_coupler_exchange_array_2
          u_p(0,:,:) = surface_coupler_exchange_array_2
          v(0,:,:)   = surface_coupler_exchange_array_3
          v_p(0,:,:) = surface_coupler_exchange_array_3

!
!--       Storage on surface data type is only done to guarantee correct output and averaging of
!--       surface data, since the assignment of pt to pt_surface is regularly done within routine
!--       surface_layer_fluxes, which is called before the surface coupling.
          DO  m = 1, surf_def%ns
             i = surf_def%i(m)
             j = surf_def%j(m)
             surf_def%pt_surface(m) = pt(0,j,i)
          ENDDO

       CASE DEFAULT

          WRITE( message_string, * ) 'surface_coupler_buffer_handling: illegal action', action
          CALL message( 'surface_coupler_buffer_handling', 'PAC0315', 3, 2, 0, 6, 1 )

    END SELECT


 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Data transfer from 1d surface-data type to 2d dummy array for equal grids in atmosphere and
!> ocean.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE transfer_1d_to_2d( surface_coupler_exchange_array, surface_data_1d )

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< running index x
       INTEGER(iwp) ::  j   !< running index y
       INTEGER(iwp) ::  m   !< running index surface type

       REAL(wp), DIMENSION(nysg:nyng,nxlg:nxrg) ::  surface_coupler_exchange_array  !< 1D surface flux, default surfaces
       REAL(wp), DIMENSION(1:surf_def%ns) ::  surface_data_1d  !< 1d surface data


!
!--    Transfer the surface flux to 2d grid. Only fluxes from upward-facing surfaces are considered.
       DO  m = 1, surf_def%ns
          IF ( surf_def%upward(m) )  THEN
             i = surf_def%i(m)
             j = surf_def%j(m)
             surface_coupler_exchange_array(j,i) = surface_data_1d(m)
          ENDIF
       ENDDO

    END SUBROUTINE transfer_1d_to_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Data transfer from 2d surface-data array to 1d surface data type for equal grids in atmosphere
!> and ocean.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE transfer_2d_to_1d( surface_coupler_exchange_array, surface_data_1d )

       IMPLICIT NONE

       INTEGER(iwp) ::  i   !< running index x
       INTEGER(iwp) ::  j   !< running index y
       INTEGER(iwp) ::  m   !< running index surface type

       REAL(wp), INTENT(IN), DIMENSION(nysg:nyng,nxlg:nxrg) ::  surface_coupler_exchange_array  !< 2d surface flux
       REAL(wp), INTENT(OUT), DIMENSION(1:surf_top%ns) ::  surface_data_1d  !< 1d surface data


!
!--    Transfer surface flux to 1d surface type.
       surface_data_1d = 0.0_wp
       DO  m = 1, surf_top%ns
          i = surf_top%i(m)
          j = surf_top%j(m)
          IF ( i >= nxl  .AND.  i <= nxr )  THEN
             IF ( j >= nys  .AND.  j <= nyn )  THEN
                surface_data_1d(m) = surface_coupler_exchange_array(j,i)
             ENDIF
          ENDIF
       ENDDO

    END SUBROUTINE transfer_2d_to_1d

 END SUBROUTINE surface_coupler_buffer_handling

 END MODULE surface_coupler_mod

