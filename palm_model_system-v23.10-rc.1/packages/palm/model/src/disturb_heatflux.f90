!> @file disturb_heatflux.f90
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
! Description:
! ------------
!> Generate random, normally distributed heatflux values and store them as the near-surface
!> heatflux.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE disturb_heatflux( surf )

    USE arrays_3d,                                                                                 &
        ONLY:  heatflux_input_conversion

    USE control_parameters,                                                                        &
        ONLY:  iran,                                                                               &
               surface_heatflux,                                                                   &
               random_generator,                                                                   &
               wall_heatflux

    USE indices,                                                                                   &
        ONLY:  nzb

    USE kinds

    USE random_generator_parallel,                                                                 &
        ONLY:  random_dummy,                                                                       &
               random_number_parallel,                                                             &
               random_seed_parallel,                                                               &
               seq_random_array

    USE surface_mod,                                                                               &
        ONLY:  surf_type

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index, x direction
    INTEGER(iwp) ::  j  !< grid index, y direction
    INTEGER(iwp) ::  k  !< grid index, z direction
    INTEGER(iwp) ::  m  !< loop variables over surface elements

    REAL(wp) ::  random_gauss  !<
    REAL(wp) ::  randomnumber  !<

    TYPE(surf_type) ::  surf   !< surface-type variable


!
!-- Generate random disturbances and store them. Note, if
!-- random_generator /= 'random-parallel' it is not guaranteed to obtain
!-- the same random distribution if the number of processors is changed.
    IF ( random_generator /= 'random-parallel' )  THEN

       DO  m = 1, surf%ns

          k = surf%k(m)

          randomnumber = random_gauss( iran, 5.0_wp )
!
!--       Only disturb the heat flux at upward-facing surfaces. Surfaces located at k=0 are set to
!--       surface_heatflux, all other surfaces are replaced by wall_heatflux(0).
          IF ( surf%upward(m) )  THEN
             IF ( k+surf%koff(m) == 0 )  THEN
                surf%shf(m) = randomnumber * surface_heatflux * heatflux_input_conversion(nzb)
             ELSE
                surf%shf(m) = randomnumber * wall_heatflux(0) * heatflux_input_conversion(k-1)
             ENDIF
          ENDIF
       ENDDO
    ELSE

       DO  m = 1, surf%ns

          i = surf%i(m)
          j = surf%j(m)
          k = surf%k(m)

          CALL random_seed_parallel( put=seq_random_array(:,j,i) )
          CALL random_number_parallel( random_dummy )
!
!--       Only disturb the heat flux at upward-facing surfaces. Surfaces located at k=0 are set to
!--       surface_heatflux, all other surfaces are replaced by wall_heatflux(0).
          IF ( surf%upward(m) )  THEN
             IF ( k+surf%koff(m) == 0 )  THEN
                surf%shf(m) = ( random_dummy - 0.5_wp ) * surface_heatflux * heatflux_input_conversion(nzb)
             ELSE
                surf%shf(m) = ( random_dummy - 0.5_wp ) * wall_heatflux(0) * heatflux_input_conversion(k-1)
             ENDIF
          ENDIF

          CALL random_seed_parallel( get=seq_random_array(:,j,i) )

       ENDDO

    ENDIF



 END SUBROUTINE disturb_heatflux
