!> @file calc_mean_profile.f90
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
!
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Calculate the horizontally averaged vertical temperature profile (pr=4 in case
!> of potential temperature, 44 in case of virtual potential temperature, and 64
!> in case of density (ocean runs)).
!------------------------------------------------------------------------------!
 MODULE calc_mean_profile_mod

#if defined( __parallel )
    USE MPI
#endif

    PRIVATE
    PUBLIC calc_mean_profile

    INTERFACE calc_mean_profile
       MODULE PROCEDURE calc_mean_profile
    END INTERFACE calc_mean_profile

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE calc_mean_profile( var, pr, force_calc )

       USE control_parameters,                                                                     &
           ONLY:  intermediate_timestep_count

       USE indices,                                                                                &
           ONLY:  ngp_2dh_s_inner,                                                                 &
                  nxl,                                                                             &
                  nxr,                                                                             &
                  nyn,                                                                             &
                  nys,                                                                             &
                  nzb,                                                                             &
                  nzb,                                                                             &
                  nzt,                                                                             &
                  topo_flags

       USE kinds

       USE pegrid

       USE statistics,                                                                             &
           ONLY:  flow_statistics_called,                                                          &
                  hom,                                                                             &
                  sums,                                                                            &
                  sums_l


       IMPLICIT NONE

       LOGICAL, OPTIONAL ::  force_calc             !< passed flag used to force calculation of mean profiles independend on data output
       LOGICAL           ::  force_calc_l = .FALSE. !< control flag

       INTEGER(iwp) ::  i                  !<
       INTEGER(iwp) ::  j                  !<
       INTEGER(iwp) ::  k                  !<
       INTEGER(iwp) ::  pr                 !<
!$     INTEGER(iwp) ::  omp_get_thread_num !<
       INTEGER(iwp) ::  tn                 !<

       REAL(wp), DIMENSION(:,:,:), POINTER ::  var

!
!--    Computation of the horizontally averaged profile of variable var, unless already done by the
!--    relevant call from flow_statistics. The calculation is done only for the first respective
!--    intermediate timestep in order to spare communication time and to produce identical model
!--    results with jobs which are calling flow_statistics at different time intervals. At
!--    initialization, intermediate_timestep_count = 0 is considered as well.
!--    Note, calc_mean_profile is also called from the radiation. Especially during the spinup when
!--    no data output is called, the values in hom do not necessarily reflect horizontal mean.
!--    The same is also for nested simulations during the spinup, where hom is initialized on a
!--    ealier point in time than the initialization via the parent domain takes place, meaning
!--    that hom does not necessarily reflect the horizontal mean during the spinup, too.
!--    In order to force the computation of horizontal mean profiles independent on data output,
!--    i.e. independent on flow_statistics, add a special control flag.
       IF ( PRESENT( force_calc ) )  THEN
          force_calc_l = force_calc
       ELSE
          force_calc_l = .FALSE.
       ENDIF

       IF ( ( .NOT. flow_statistics_called  .AND.  intermediate_timestep_count <= 1 ) .OR.         &
            force_calc_l )  THEN

!
!--       Horizontal average of variable var
          tn           =   0  ! Default thread number in case of one thread
          !$OMP PARALLEL PRIVATE( i, j, k, tn )
!$        tn = omp_get_thread_num()
          sums_l(:,pr,tn) = 0.0_wp
          !$OMP DO
          DO  i = nxl, nxr
             DO  j =  nys, nyn
                DO  k = nzb, nzt+1
                   sums_l(k,pr,tn) = sums_l(k,pr,tn) + var(k,j,i) * MERGE( 1.0_wp, 0.0_wp,         &
                                                            BTEST( topo_flags(k,j,i), 22 ) )
                ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL

          DO  i = 1, threads_per_task-1
             sums_l(:,pr,0) = sums_l(:,pr,0) + sums_l(:,pr,i)
          ENDDO

#if defined( __parallel )

          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( sums_l(nzb,pr,0), sums(nzb,pr), nzt+2-nzb, MPI_REAL, MPI_SUM, comm2d,&
                              ierr )

#else

          sums(:,pr) = sums_l(:,pr,0)

#endif

          DO  k = nzb, nzt+1
             IF ( ngp_2dh_s_inner(k,0) /= 0 )  THEN
                hom(k,1,pr,0) = sums(k,pr) / ngp_2dh_s_inner(k,0)
             ENDIF
          ENDDO

       ENDIF


    END SUBROUTINE calc_mean_profile

 END MODULE calc_mean_profile_mod
