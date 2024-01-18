!> @file biometeorology_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of PALM-4U.
!
! PALM-4U is free software: You can redistribute it and/or modify it under the terms of the GNU
! General Public License as published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! PALM-4U is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
! the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 2018-2021 Deutscher Wetterdienst (DWD)
! Copyright 2018-2021 Institute of Computer Science, Academy of Sciences, Prague
! Copyright 2018-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Dominik Froehlich <dominik.froehlich@dwd.de>, thermal indices
! @author Jaroslav Resler <resler@cs.cas.cz>, mean radiant temperature
! @author Michael Schrempf <schrempf@muk.uni-hannover.de>, uv exposure
!
!
! Description:
! ------------
!> Biometeorology module consisting of two parts:
!> 1.: Human thermal comfort module calculating thermal perception of a sample human being under the
!> current meteorological conditions.
!> 2.: Calculation of vitamin-D weighted UV exposure
!>
!> @todo Alphabetical sorting of "USE ..." lists, "ONLY" list, variable declarations
!>       (per subroutine: first all CHARACTERs, then INTEGERs, LOGICALs, REALs, )
!> @todo Comments start with capital letter --> "!-- Include..."
!> @todo uv_vitd3dose-->new output type necessary (cumulative)
!> @todo consider upwelling radiation in UV
!> @todo re-design module to work with PALM's module interface and reduce number of workarounds
!> @todo some queries for initial values (e.g. <= -998.0) should be adjusted to the new initial
!>       value -999999.0
!>
!> @note nothing now
!>
!> @bug  checks for proper parameter settings and required input data are missing. Currently
!<       implemented only by a workaround!
!--------------------------------------------------------------------------------------------------!
 MODULE biometeorology_mod

    USE arrays_3d,                                                                                 &
        ONLY:  p,                                                                                  &
               pt,                                                                                 &
               q,                                                                                  &
               u,                                                                                  &
               v,                                                                                  &
               w

    USE averaging,                                                                                 &
        ONLY:  pt_av,                                                                              &
               q_av,                                                                               &
               u_av,                                                                               &
               v_av,                                                                               &
               w_av

    USE basic_constants_and_equations_mod,                                                         &
        ONLY: c_p,                                                                                 &
              degc_to_k,                                                                           &
              l_v,                                                                                 &
              magnus,                                                                              &
              pi,                                                                                  &
              sigma_sb

    USE control_parameters,                                                                        &
        ONLY:  average_count_3d,                                                                   &
               biometeorology,                                                                     &
               cyclic_fill_initialization,                                                         &
               debug_output,                                                                       &
               dz,                                                                                 &
               dz_stretch_factor,                                                                  &
               dz_stretch_level,                                                                   &
               humidity,                                                                           &
               initializing_actions, message_string,                                               &
               nz_do3d,                                                                            &
               output_fill_value,                                                                  &
               restart_data_format_output,                                                         &
               surface_pressure

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddy,                                                                                &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_top_ind

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY: building_obstruction_f,                                                              &
              input_file_uvem,                                                                     &
              input_pids_uvem,                                                                     &
              netcdf_data_input_uvem,                                                              &
              uvem_integration_f,                                                                  &
              uvem_irradiance_f,                                                                   &
              uvem_projarea_f,                                                                     &
              uvem_radiance_f

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time

    USE radiation_model_mod,                                                                       &
        ONLY:  id,                                                                                 &
               ix,                                                                                 &
               iy,                                                                                 &
               iz,                                                                                 &
               mrt_include_sw,                                                                     &
               mrt_minlevel,                                                                       &
               mrt_nlevels,                                                                        &
               mrtbl,                                                                              &
               mrtinlw,                                                                            &
               mrtinsw,                                                                            &
               nmrtbl,                                                                             &
               radiation,                                                                          &
               rad_lw_in,                                                                          &
               rad_lw_out,                                                                         &
               rad_sw_in,                                                                          &
               rad_sw_out,                                                                         &
               radiation_interactions

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rrd_mpi_io,                                                                         &
               rd_mpi_io_check_array,                                                              &
               wrd_mpi_io


    IMPLICIT NONE

    REAL(wp), PARAMETER ::  bio_initial_value = -999999.0_wp  !< initial value for some quantities
    REAL(wp), PARAMETER ::  human_absorb = 0.7_wp             !< SW absorbtivity of a human body (Fanger 1972)
    REAL(wp), PARAMETER ::  human_emiss = 0.97_wp             !< LW emissivity of a human body after (Fanger 1972)

    INTEGER(iwp) ::  bio_cell_level  !< cell level biom calculates for

    LOGICAL ::  thermal_comfort        = .FALSE.  !< Enables or disables the entire thermal comfort part
    LOGICAL ::  do_average_theta       = .FALSE.  !< switch: do theta averaging in this module? (if .FALSE. this is done globally)
    LOGICAL ::  do_average_q           = .FALSE.  !< switch: do e averaging in this module?
    LOGICAL ::  do_average_u           = .FALSE.  !< switch: do u averaging in this module?
    LOGICAL ::  do_average_v           = .FALSE.  !< switch: do v averaging in this module?
    LOGICAL ::  do_average_w           = .FALSE.  !< switch: do w averaging in this module?
    LOGICAL ::  do_average_mrt         = .FALSE.  !< switch: do mrt averaging in this module?
    LOGICAL ::  average_trigger_perct  = .FALSE.  !< update averaged input on call to bio_perct?
    LOGICAL ::  average_trigger_utci   = .FALSE.  !< update averaged input on call to bio_utci?
    LOGICAL ::  average_trigger_pet    = .FALSE.  !< update averaged input on call to bio_pet?
    LOGICAL ::  average_trigger_mrt    = .FALSE.  !< update averaged input on call to bio_pet?
    LOGICAL ::  do_calculate_perct     = .FALSE.  !< Turn index PT (instant. input) on or off
    LOGICAL ::  do_calculate_perct_av  = .FALSE.  !< Turn index PT (averaged input) on or off
    LOGICAL ::  do_calculate_pet       = .FALSE.  !< Turn index PET (instant. input) on or off
    LOGICAL ::  do_calculate_pet_av    = .FALSE.  !< Turn index PET (averaged input)                                                                        & on or off
    LOGICAL ::  do_calculate_utci      = .FALSE.  !< Turn index UTCI (instant. input) on or off
    LOGICAL ::  do_calculate_utci_av   = .FALSE.  !< Turn index UTCI (averaged input) on or off
    LOGICAL ::  do_calculate_mrt2d     = .FALSE.  !< Turn index MRT 2D (averaged or inst) on or off

    REAL(wp)    ::  bio_output_height  !< height output is calculated in m

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  perct         !< PT results   (degree_C)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  perct_av      !< PT results (aver. input)   (degree_C)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pet           !< PET results  (degree_C)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  pet_av        !< PET results (aver. input)  (degree_C)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tmrt_grid     !< tmrt results (degree_C)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  tmrt_av_grid  !< tmrt results (degree_C)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  utci          !< UTCI results (degree_C)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  utci_av       !< UTCI results (aver. input) (degree_C)

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  mrt_av_grid   !< time average mean radiant temperature

!
!-- UVEM parameters from here
    INTEGER(iwp) ::  ai                      = 0  !< loop index in azimuth direction
    INTEGER(iwp) ::  bi                      = 0  !< loop index of bit location within an 8bit-integer (one Byte)
    INTEGER(iwp) ::  bio_nmrtbl
    INTEGER(iwp) ::  clothing                = 1  !< clothing (0=unclothed, 1=Arms,Hands,Face free, 3=Hand,Face free)
    INTEGER(iwp) ::  iq                      = 0  !< loop index of irradiance quantity
    INTEGER(iwp) ::  pobi                    = 0  !< loop index of the position of corresponding byte within ibset byte vektor
    INTEGER(iwp) ::  obstruction_direct_beam = 0  !< Obstruction information for direct beam
    INTEGER(iwp) ::  zi                      = 0  !< loop index in zenith direction

    INTEGER(ibp), DIMENSION(0:44)  ::  obstruction_temp1 = 0  !< temporary obstruction information stored with ibset
    INTEGER(iwp), DIMENSION(0:359) ::  obstruction_temp2 = 0  !< restored temporary obstruction information from ibset file

    INTEGER(iwp), DIMENSION(0:35,0:9) ::  obstruction       = 1  !< final 2D obstruction information array

    LOGICAL ::  consider_obstructions = .TRUE.   !< namelist parameter (see documentation)
    LOGICAL ::  sun_in_south          = .FALSE.  !< namelist parameter (see documentation)
    LOGICAL ::  turn_to_sun           = .TRUE.   !< namelist parameter (see documentation)
    LOGICAL ::  uv_exposure           = .FALSE.  !< namelist parameter (see documentation)

    REAL(wp) ::  diffuse_exposure            =   0.0_wp  !< calculated exposure by diffuse radiation
    REAL(wp) ::  direct_exposure             =   0.0_wp  !< calculated exposure by direct solar beam
    REAL(wp) ::  orientation_angle           =   0.0_wp  !< orientation of front/face of the human model
    REAL(wp) ::  projection_area_direct_beam =   0.0_wp  !< projection area for direct solar beam
    REAL(wp) ::  saa                         = 180.0_wp  !< solar azimuth angle
    REAL(wp) ::  startpos_human              =   0.0_wp  !< start value for azimuth interpolation of human geometry array
    REAL(wp) ::  startpos_saa_float          =   0.0_wp  !< start value for azimuth interpolation of radiance array
    REAL(wp) ::  sza                         =  20.0_wp  !< solar zenith angle
    REAL(wp) ::  xfactor                     =   0.0_wp  !< relative x-position used for interpolation
    REAL(wp) ::  yfactor                     =   0.0_wp  !< relative y-position used for interpolation

    REAL(wp), DIMENSION(0:2)  ::  irradiance =   0.0_wp  !< iradiance values extracted from irradiance lookup table

    REAL(wp), DIMENSION(0:2,0:90) ::  irradiance_lookup_table      = 0.0_wp  !< irradiance lookup table
    REAL(wp), DIMENSION(0:35,0:9) ::  integration_array            = 0.0_wp  !< solid angle factors for hemispherical integration
    REAL(wp), DIMENSION(0:35,0:9) ::  projection_area              = 0.0_wp  !< projection areas of a human (all directions)
    REAL(wp), DIMENSION(0:35,0:9) ::  projection_area_lookup_table = 0.0_wp  !< human geometry lookup table (projection areas)
    REAL(wp), DIMENSION(0:71,0:9) ::  projection_area_direct_temp  = 0.0_wp  !< temporary projection area for direct solar beam
    REAL(wp), DIMENSION(0:71,0:9) ::  projection_area_temp         = 0.0_wp  !< temporary projection area for all directions
    REAL(wp), DIMENSION(0:35,0:9) ::  radiance_array               = 0.0_wp  !< radiance extracted from radiance_lookup_table
    REAL(wp), DIMENSION(0:71,0:9) ::  radiance_array_temp          = 0.0_wp  !< temporary radiance data

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vitd3_exposure  !< result variable for instantaneous vitamin-D weighted exposures
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  vitd3_dose      !< result variable for summation of vitamin-D weighted exposures

    REAL(wp), DIMENSION(0:35,0:9,0:90) ::  radiance_lookup_table   = 0.0_wp  !< radiance lookup table


    PRIVATE

    PUBLIC bio_3d_data_averaging,                                                                  &
           bio_calculate_mrt_grid,                                                                 &
           bio_calculate_thermal_index_maps,                                                       &
           bio_calculate_uv_exposure,                                                              &
           bio_calc_ipt,                                                                           &
           bio_check_data_output,                                                                  &
           bio_check_parameters,                                                                   &
           bio_data_output_2d,                                                                     &
           bio_data_output_3d,                                                                     &
           bio_define_netcdf_grid,                                                                 &
           bio_get_thermal_index_input_ij,                                                         &
           bio_header,                                                                             &
           bio_init,                                                                               &
           bio_init_checks,                                                                        &
           bio_nmrtbl,                                                                             &
           bio_parin,                                                                              &
           bio_rrd_global,                                                                         &
           bio_rrd_local,                                                                          &
           bio_vm_sampling,                                                                        &
           bio_wrd_global,                                                                         &
           bio_wrd_local,                                                                          &
           thermal_comfort,                                                                        &
           uv_exposure

    INTERFACE bio_3d_data_averaging
       MODULE PROCEDURE bio_3d_data_averaging
    END INTERFACE bio_3d_data_averaging

    INTERFACE bio_calculate_mrt_grid
       MODULE PROCEDURE bio_calculate_mrt_grid
    END INTERFACE bio_calculate_mrt_grid

    INTERFACE bio_calculate_thermal_index_maps
       MODULE PROCEDURE bio_calculate_thermal_index_maps
    END INTERFACE bio_calculate_thermal_index_maps

    INTERFACE bio_calc_ipt
       MODULE PROCEDURE bio_calc_ipt
    END INTERFACE bio_calc_ipt

    INTERFACE bio_check_data_output
       MODULE PROCEDURE bio_check_data_output
    END INTERFACE bio_check_data_output

    INTERFACE bio_check_parameters
       MODULE PROCEDURE bio_check_parameters
    END INTERFACE bio_check_parameters

    INTERFACE bio_data_output_2d
       MODULE PROCEDURE bio_data_output_2d
    END INTERFACE bio_data_output_2d

    INTERFACE bio_data_output_3d
       MODULE PROCEDURE bio_data_output_3d
    END INTERFACE bio_data_output_3d

    INTERFACE bio_define_netcdf_grid
       MODULE PROCEDURE bio_define_netcdf_grid
    END INTERFACE bio_define_netcdf_grid

    INTERFACE bio_get_thermal_index_input_ij
       MODULE PROCEDURE bio_get_thermal_index_input_ij
    END INTERFACE bio_get_thermal_index_input_ij

    INTERFACE bio_header
       MODULE PROCEDURE bio_header
    END INTERFACE bio_header

    INTERFACE bio_init
       MODULE PROCEDURE bio_init
    END INTERFACE bio_init

    INTERFACE bio_init_checks
       MODULE PROCEDURE bio_init_checks
    END INTERFACE bio_init_checks

    INTERFACE bio_parin
       MODULE PROCEDURE bio_parin
    END INTERFACE bio_parin

    INTERFACE bio_rrd_global
       MODULE PROCEDURE bio_rrd_global_ftn
       MODULE PROCEDURE bio_rrd_global_mpi
    END INTERFACE bio_rrd_global

    INTERFACE bio_rrd_local
       MODULE PROCEDURE bio_rrd_local_ftn
       MODULE PROCEDURE bio_rrd_local_mpi
    END INTERFACE bio_rrd_local

    INTERFACE bio_vm_sampling
       MODULE PROCEDURE bio_vm_sampling
    END INTERFACE bio_vm_sampling

    INTERFACE bio_wrd_global
       MODULE PROCEDURE bio_wrd_global
    END INTERFACE bio_wrd_global

    INTERFACE bio_wrd_local
       MODULE PROCEDURE bio_wrd_local
    END INTERFACE bio_wrd_local

    INTERFACE bio_calculate_uv_exposure
       MODULE PROCEDURE bio_calculate_uv_exposure
    END INTERFACE bio_calculate_uv_exposure

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average biom input quantities as well as allocate the array necessary for storing
!> the average.
!> There is a considerable difference to the 3d_data_averaging subroutines used by other modules:
!> For the thermal indices, the module needs to average the input conditions, not the result!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_3d_data_averaging( mode, variable )

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  mode     !< Averaging mode: allocate, sum, or average
    CHARACTER (LEN=*) ::  variable !< The variable in question

    INTEGER(iwp) ::  i        !< Running index, x-direction
    INTEGER(iwp) ::  j        !< Running index, y-direction
    INTEGER(iwp) ::  k        !< Running index, z-direction
    INTEGER(iwp) ::  l        !< index used to link radiation arrays to 3d grid arrays


    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'bio_mrt' )

                IF ( .NOT. ALLOCATED( mrt_av_grid ) )  THEN
                   ALLOCATE( mrt_av_grid(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                mrt_av_grid = 0.0_wp
                do_average_mrt = .FALSE.  !< Overwrite if that was enabled somehow


          CASE ( 'bio_perct*', 'bio_utci*', 'bio_pet*', 'bio_mrt*' )

!
!--          Averaging, as well as the allocation of the required grids must be done only once,
!--          independent of for how many thermal indices averaged output is desired.
!--          Therefore we need to memorize which index is the one that controls the averaging
!--          (what must be the first thermal index called).
!--          Indices are in unknown order as depending on the input file, determine first index to
!--          average und update only once.
!
!--          Only proceed here if this was not done for any index before. This is done only once
!--          during the whole model run.
             IF ( .NOT. average_trigger_perct  .AND.                                               &
                  .NOT. average_trigger_utci   .AND.                                               &
                  .NOT. average_trigger_pet    .AND.                                               &
                  .NOT. average_trigger_mrt )  THEN
!
!--             Memorize the first index called to control averaging
                IF ( TRIM( variable ) == 'bio_perct*' )  THEN
                    average_trigger_perct = .TRUE.
                ENDIF
                IF ( TRIM( variable ) == 'bio_utci*' )  THEN
                    average_trigger_utci = .TRUE.
                ENDIF
                IF ( TRIM( variable ) == 'bio_pet*' )  THEN
                    average_trigger_pet = .TRUE.
                ENDIF
                IF ( TRIM( variable ) == 'bio_mrt*' )  THEN
                    average_trigger_mrt = .TRUE.
                ENDIF
             ENDIF
!
!--          Allocation of the input _av grids was moved to the "sum" section to make sure averaging
!--          is only done once!

          CASE ( 'uvem_vitd3dose*' )

             IF ( .NOT. ALLOCATED( vitd3_dose ) )  THEN
                ALLOCATE( vitd3_dose(nysg:nyng,nxlg:nxrg) )
             ENDIF
             vitd3_dose = 0.0_wp

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'bio_mrt' )
!
!--          Consider the case 'bio_mrt' is called after some thermal index. In that case
!            do_average_mrt will be .TRUE. leading to a double-averaging.
             IF ( .NOT. do_average_mrt  .AND.  ALLOCATED( mrt_av_grid ) )  THEN

                IF ( mrt_include_sw )  THEN
                   DO  l = 1, nmrtbl
                      i = mrtbl(ix,l)
                      j = mrtbl(iy,l)
                      k = mrtbl(iz,l)
                      mrt_av_grid(k,j,i) = mrt_av_grid(k,j,i)  +                                   &
                                    ( ( human_absorb * mrtinsw(l) +                                &
                                    mrtinlw(l) ) / ( human_emiss * sigma_sb ) )**0.25_wp - degc_to_k
                   ENDDO
                ELSE
                   DO  l = 1, nmrtbl
                      i = mrtbl(ix,l)
                      j = mrtbl(iy,l)
                      k = mrtbl(iz,l)
                      mrt_av_grid(k,j,i) = mrt_av_grid(k,j,i)  +                                   &
                                    ( mrtinlw(l) / ( human_emiss * sigma_sb ) )**0.25_wp - degc_to_k
                   ENDDO
                ENDIF
             ENDIF

          CASE ( 'bio_perct*', 'bio_utci*', 'bio_pet*', 'bio_mrt*' )
!
!--          Only continue if the current index is the one to trigger the input averaging, see
!--          above.
             IF ( average_trigger_perct  .AND.  TRIM( variable ) /= 'bio_perct*')    RETURN
             IF ( average_trigger_utci   .AND.  TRIM( variable ) /= 'bio_utci*' )    RETURN
             IF ( average_trigger_pet    .AND.  TRIM( variable ) /= 'bio_pet*'  )    RETURN
             IF ( average_trigger_mrt    .AND.  TRIM( variable ) /= 'bio_mrt*'  )    RETURN
!
!--          Now memorize which of the input grids are not averaged by other modules. Set averaging
!--          switch to .TRUE. and allocate the respective grid in that case.
             IF ( .NOT. ALLOCATED( pt_av ) )  THEN  !< if not averaged by other module
                ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                do_average_theta = .TRUE.  !< memorize, that bio is responsible
                pt_av = 0.0_wp
             ENDIF
             IF ( ALLOCATED( pt_av )  .AND.  do_average_theta )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         pt_av(k,j,i) = pt_av(k,j,i) + pt(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( .NOT. ALLOCATED( q_av ) )  THEN
                ALLOCATE( q_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                do_average_q = .TRUE.
                q_av = 0.0_wp
             ENDIF
             IF ( ALLOCATED( q_av )  .AND.  do_average_q )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         q_av(k,j,i) = q_av(k,j,i) + q(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

!
!--          u_av, v_av and w_av are always allocated
             IF ( .NOT. ALLOCATED( u_av ) )  THEN
                ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                do_average_u = .TRUE.
                u_av = 0.0_wp
             ENDIF
             IF ( ALLOCATED( u_av )  .AND.  do_average_u )  THEN
                DO  i = nxlg, nxrg       !< yes, ghost points are required here!
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         u_av(k,j,i) = u_av(k,j,i) + u(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( .NOT. ALLOCATED( v_av ) )  THEN
                ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                do_average_v = .TRUE.
                v_av = 0.0_wp
             ENDIF
             IF ( ALLOCATED( v_av )  .AND.  do_average_v )  THEN
                DO  i = nxlg, nxrg       !< yes, ghost points are required here!
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         v_av(k,j,i) = v_av(k,j,i) + v(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( .NOT. ALLOCATED( w_av ) )  THEN
                ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                do_average_w = .TRUE.
                w_av = 0.0_wp
             ENDIF
             IF ( ALLOCATED( w_av )  .AND.  do_average_w )  THEN
                DO  i = nxlg, nxrg       !< yes, ghost points are required here!
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         w_av(k,j,i) = w_av(k,j,i) + w(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( .NOT. ALLOCATED( mrt_av_grid ) )  THEN
                ALLOCATE( mrt_av_grid(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                do_average_mrt = .TRUE.
                mrt_av_grid = 0.0_wp
             ENDIF
             IF ( ALLOCATED( mrt_av_grid )  .AND.  do_average_mrt )  THEN

                IF ( mrt_include_sw )  THEN
                   DO  l = 1, nmrtbl
                      i = mrtbl(ix,l)
                      j = mrtbl(iy,l)
                      k = mrtbl(iz,l)
                      mrt_av_grid(k,j,i) = mrt_av_grid(k,j,i)  +                                   &
                                    ( ( human_absorb * mrtinsw(l) +                                &
                                    mrtinlw(l) ) /                                                 &
                                      ( human_emiss  * sigma_sb ) )**0.25_wp - degc_to_k
                   ENDDO
                ELSE
                   DO  l = 1, nmrtbl
                      i = mrtbl(ix,l)
                      j = mrtbl(iy,l)
                      k = mrtbl(iz,l)
                      mrt_av_grid(k,j,i) = mrt_av_grid(k,j,i) +                                    &
                                    ( mrtinlw(l) /                                                 &
                                    ( human_emiss * sigma_sb ) )**0.25_wp - degc_to_k
                   ENDDO
                ENDIF
             ENDIF
!
!--       This is a cumulated dose. No mode == 'average' for this quantity.
          CASE ( 'uvem_vitd3dose*' )
             IF ( ALLOCATED( vitd3_dose ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      vitd3_dose(j,i) = vitd3_dose(j,i) + vitd3_exposure(j,i)
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

          CASE ( 'bio_mrt' )
!
!--          Consider the case 'bio_mrt' is called after some thermal index. In that case
!--          do_average_mrt will be .TRUE. leading to a double-averaging.
             IF ( .NOT. do_average_mrt  .AND.  ALLOCATED( mrt_av_grid ) )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         mrt_av_grid(k,j,i) = mrt_av_grid(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'bio_perct*', 'bio_utci*', 'bio_pet*', 'bio_mrt*' )
!
!--          Only continue if update index, see above
             IF ( average_trigger_perct  .AND.                                                     &
                TRIM( variable ) /= 'bio_perct*' )  RETURN
             IF ( average_trigger_utci  .AND.                                                      &
                TRIM( variable ) /= 'bio_utci*' )  RETURN
             IF ( average_trigger_pet   .AND.                                                      &
                TRIM( variable ) /= 'bio_pet*' )  RETURN
             IF ( average_trigger_mrt   .AND.                                                      &
                TRIM( variable ) /= 'bio_mrt*' )  RETURN

             IF ( ALLOCATED( pt_av )  .AND.  do_average_theta )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         pt_av(k,j,i) = pt_av(k,j,i) /  REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( ALLOCATED( q_av )  .AND.  do_average_q )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         q_av(k,j,i) = q_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( ALLOCATED( u_av )  .AND.  do_average_u )  THEN
                DO  i = nxlg, nxrg       !< yes, ghost points are required here!
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         u_av(k,j,i) = u_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( ALLOCATED( v_av )  .AND.  do_average_v )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         v_av(k,j,i) = v_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( ALLOCATED( w_av )  .AND.  do_average_w )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         w_av(k,j,i) = w_av(k,j,i) / REAL( average_count_3d, KIND = wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

             IF ( ALLOCATED( mrt_av_grid )  .AND.  do_average_mrt )  THEN
                DO  i = nxl, nxr
                   DO  j = nys, nyn
                      DO  k = nzb, nzt+1
                         mrt_av_grid(k,j,i) = mrt_av_grid(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

!
!--     No averaging for UVEM since we are calculating a dose (only sum is calculated and saved to
!--     av.nc file)
        END SELECT

    ENDIF


 END SUBROUTINE bio_3d_data_averaging



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check data output for biometeorology model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_check_data_output( var, unit, i, j, ilen, k )

    USE control_parameters,                                                                        &
        ONLY: data_output,                                                                         &
              message_string

    IMPLICIT NONE

    CHARACTER (LEN=*) ::  unit     !< The unit for the variable var
    CHARACTER (LEN=*) ::  var      !< The variable in question

    INTEGER(iwp), INTENT(IN) ::  i     !< Current element of data_output
    INTEGER(iwp), INTENT(IN) ::  j     !< Average quantity? 0 = no, 1 = yes
    INTEGER(iwp), INTENT(IN) ::  ilen  !< Length of current entry in data_output
    INTEGER(iwp), INTENT(IN) ::  k     !< Output is xy mode? 0 = no, 1 = yes


    SELECT CASE ( TRIM( var ) )
!
!--    Allocate a temporary array with the desired output dimensions.
!--    Arrays for time-averaged thermal indices are also allocated here because they are not running
!--    through the standard averaging procedure in bio_3d_data_averaging as the values of the
!--    averaged thermal indices are derived in a single step based on priorly averaged arrays (see
!--    bio_calculate_thermal_index_maps).
       CASE ( 'bio_mrt', 'bio_mrt*' )
          unit = 'degree_C'
          thermal_comfort = .TRUE.  !< enable thermal_comfort if user forgot to do so
          IF ( .NOT. ALLOCATED( tmrt_grid ) )  THEN
             ALLOCATE( tmrt_grid (nys:nyn,nxl:nxr) )
             tmrt_grid = output_fill_value
          ENDIF
          IF ( TRIM( var ) == 'bio_mrt*' )  THEN
             do_calculate_mrt2d = .TRUE.
          END IF

       CASE ( 'bio_perct*' )
          unit = 'degree_C'
          thermal_comfort = .TRUE.
          IF ( j == 0 )  THEN                !< if instantaneous input
             do_calculate_perct = .TRUE.
             IF ( .NOT. ALLOCATED( perct ) )  THEN
                ALLOCATE( perct (nys:nyn,nxl:nxr) )
                perct = output_fill_value
             ENDIF
          ELSE                               !< if averaged input
             do_calculate_perct_av = .TRUE.
             IF ( .NOT. ALLOCATED( perct_av ) )  THEN
                ALLOCATE( perct_av (nys:nyn,nxl:nxr) )
                perct_av = 0.0_wp
             ENDIF
          ENDIF

       CASE ( 'bio_utci*' )
          unit = 'degree_C'
          thermal_comfort = .TRUE.
          IF ( j == 0 )  THEN
             do_calculate_utci = .TRUE.
             IF ( .NOT. ALLOCATED( utci ) )  THEN
                ALLOCATE( utci (nys:nyn,nxl:nxr) )
                utci = output_fill_value
             ENDIF
          ELSE
             do_calculate_utci_av = .TRUE.
             IF ( .NOT. ALLOCATED( utci_av ) )  THEN
                ALLOCATE( utci_av (nys:nyn,nxl:nxr) )
                utci_av = 0.0_wp
             ENDIF
          ENDIF

       CASE ( 'bio_pet*' )
          unit = 'degree_C'
          thermal_comfort = .TRUE.
          IF ( j == 0 )  THEN
             do_calculate_pet = .TRUE.
             IF ( .NOT. ALLOCATED( pet ) )  THEN
                ALLOCATE( pet (nys:nyn,nxl:nxr) )
                pet = output_fill_value
             ENDIF
          ELSE
             do_calculate_pet_av = .TRUE.
             IF ( .NOT. ALLOCATED( pet_av ) )  THEN
                ALLOCATE( pet_av (nys:nyn,nxl:nxr) )
                pet_av = 0.0_wp
             ENDIF
          ENDIF


       CASE ( 'uvem_vitd3*' )
          IF ( .NOT. uv_exposure )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requires uv_exposure = .TRUE.' // &
                              '&in namelist "biometeorology_parameters"'
             CALL message( 'uvem_check_data_output', 'BIO0001', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
             message_string = 'illegal value for data_output: "' //                                &
                              TRIM( var ) // '" & only 2d-horizontal ' //                          &
                              'cross sections are allowed for this value'
             CALL message( 'check_parameters', 'BIO0002', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'IU/s'
          IF ( .NOT. ALLOCATED( vitd3_exposure ) )  THEN
             ALLOCATE( vitd3_exposure(nysg:nyng,nxlg:nxrg) )
          ENDIF
          vitd3_exposure = 0.0_wp

       CASE ( 'uvem_vitd3dose*' )
          IF (  .NOT.  uv_exposure )  THEN
             message_string = 'output of "' // TRIM( var ) // '" requires uv_exposure = .TRUE.' // &
                              '&in namelist "biometeorology_parameters"'
             CALL message( 'uvem_check_data_output', 'BIO0001', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
             message_string = 'illegal value for data_output: "' //                                &
                              TRIM( var ) // '" & only 2d-horizontal ' //                          &
                              'cross sections are allowed for this value'
             CALL message( 'check_parameters', 'BIO0002', 1, 2, 0, 6, 0 )
          ENDIF
          unit = 'IU/av-h'
          IF ( .NOT. ALLOCATED( vitd3_dose ) )  THEN
             ALLOCATE( vitd3_dose(nysg:nyng,nxlg:nxrg) )
          ENDIF
          vitd3_dose = 0.0_wp

       CASE DEFAULT
          unit = 'illegal'

    END SELECT

!
!-- Further checks if thermal comfort output is desired.
    IF ( thermal_comfort  .AND.  unit(1:6) == 'degree' )  THEN
!
!--    Break if required modules "radiation" is not available.
       IF ( .NOT.  radiation )  THEN
          message_string = 'output of "' // TRIM( var ) // '" requires radiation = .TRUE.'
          CALL message( 'check_parameters', 'BIO0003', 1, 2, 0, 6, 0 )
          unit = 'illegal'
       ENDIF
!
!--    All "thermal_comfort" outputs except from 'bio_mrt' will also need  humidity input. Check
!--    also for that.
       IF ( TRIM( var ) /= 'bio_mrt' )  THEN
          IF ( .NOT.  humidity )  THEN
             message_string = 'thermal comfort requires humidity = .TRUE.'
             CALL message( 'check_parameters', 'BIO0004', 1, 2, 0, 6, 0 )
             unit = 'illegal'
          ENDIF
       ENDIF


    ENDIF

 END SUBROUTINE bio_check_data_output

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for biom module
!> Currently unused but might come in handy for future checks?
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_check_parameters


    IMPLICIT NONE

!
!-- Check settings for UV exposure part
    IF ( uv_exposure )  THEN

!
!--    Input file not present
       IF ( .NOT. input_pids_uvem )  THEN
          WRITE( message_string, * ) 'uv_exposure = .TRUE. but input file "' //                    &
                                     TRIM( input_file_uvem ) // '" is not present'
          CALL message( 'bio_check_parameters', 'BIO0005', 1, 2, 0, 6, 0 )
       ELSE

!
!--       Required variables not given in input file
          IF ( .NOT. uvem_integration_f%from_file  .OR.  .NOT. uvem_irradiance_f%from_file  .OR.   &
               .NOT. uvem_projarea_f%from_file  .OR.  .NOT. uvem_radiance_f%from_file )  THEN
             WRITE( message_string, * ) 'uv_exposure = .TRUE. but one or more required input ' //  &
                                        'varaibles are not present in file "' //                   &
                                        TRIM( input_file_uvem ) // '"'
             CALL message( 'bio_check_parameters', 'BIO0006', 1, 2, 0, 6, 0 )
          ENDIF

!
!--       Obstruction requested but not given
          IF ( consider_obstructions  .AND.  .NOT. building_obstruction_f%from_file )  THEN
             WRITE( message_string, * ) 'consider_obstructions = .TRUE. but varaible ' //          &
                                        '"obstruction" is not present in file "' //                &
                                        TRIM( input_file_uvem ) // '"'
             CALL message( 'bio_check_parameters', 'BIO0007', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
    ENDIF

 END SUBROUTINE bio_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 2D output variables
!> data_output_2d 1188ff
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_data_output_2d( av, variable, found, grid, local_pf, two_d, nzb_do, nzt_do)


    USE kinds


    IMPLICIT NONE
!
!-- Input variables
    CHARACTER (LEN=*), INTENT(IN) ::  variable    !< Char identifier to select var for output
    INTEGER(iwp), INTENT(IN)      ::  av          !< Use averaged data? 0 = no, 1 = yes?
    INTEGER(iwp), INTENT(IN)      ::  nzb_do      !< Unused. 2D. nz bottom to nz top
    INTEGER(iwp), INTENT(IN)      ::  nzt_do      !< Unused.
!
!-- Output variables
    CHARACTER (LEN=*), INTENT(OUT) ::  grid   !< Grid type (always "zu1" for biom)

    LOGICAL, INTENT(OUT)           ::  found  !< Output found?
    LOGICAL, INTENT(OUT)           ::  two_d  !< Flag parameter that indicates 2D variables,
                                              !< horizontal cross sections, must be .TRUE. for thermal indices and uv
    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !< Temp. result grid to return
!
!-- Internal variables
    INTEGER(iwp) ::  i        !< Running index, x-dir
    INTEGER(iwp) ::  j        !< Running index, y-dir
    INTEGER(iwp) ::  k        !< Running index, z-dir
    INTEGER(iwp) ::  l        !< Running index, radiation grid


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )


        CASE ( 'bio_mrt_xy' )
           grid = 'zu1'
           two_d = .FALSE.  !< can be calculated for several levels
           DO  l = 1, nmrtbl
              i = mrtbl(ix,l)
              j = mrtbl(iy,l)
              k = mrtbl(iz,l)
              IF ( k < nzb_do  .OR.  k > nzt_do  .OR.  j < nys  .OR.                               &
                 j > nyn  .OR.  i < nxl  .OR.  i > nxr )  CYCLE
              IF ( av == 0 )  THEN
                 IF ( mrt_include_sw )  THEN
                    local_pf(i,j,k) = ( ( human_absorb * mrtinsw(l) +                              &
                                      mrtinlw(l) ) /                                               &
                                      ( human_emiss * sigma_sb ) )**0.25_wp - degc_to_k
                 ELSE
                    local_pf(i,j,k) = ( mrtinlw(l) /                                               &
                                      ( human_emiss * sigma_sb ) )**0.25_wp - degc_to_k
                 ENDIF
              ELSE
                 local_pf(i,j,k) = mrt_av_grid(k,j,i)
              ENDIF
           ENDDO

        CASE ( 'bio_mrt*_xy' )        ! 2d-array
           grid = 'zu1'
           two_d = .TRUE.
           IF ( av == 0 )  THEN
              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    local_pf(i,j,nzb+1) = tmrt_grid(j,i)
                 ENDDO
              ENDDO
           ELSE
              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    local_pf(i,j,nzb+1) = tmrt_av_grid(j,i)
                 ENDDO
              ENDDO
           ENDIF


        CASE ( 'bio_perct*_xy' )        ! 2d-array
           grid = 'zu1'
           two_d = .TRUE.
           IF ( av == 0 )  THEN
              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    local_pf(i,j,nzb+1) = perct(j,i)
                 ENDDO
              ENDDO
           ELSE
              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    local_pf(i,j,nzb+1) = perct_av(j,i)
                 ENDDO
              ENDDO
           ENDIF


        CASE ( 'bio_utci*_xy' )        ! 2d-array
           grid = 'zu1'
           two_d = .TRUE.
           IF ( av == 0 )  THEN
              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    local_pf(i,j,nzb+1) = utci(j,i)
                 ENDDO
              ENDDO
           ELSE
              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    local_pf(i,j,nzb+1) = utci_av(j,i)
                 ENDDO
              ENDDO
           ENDIF


        CASE ( 'bio_pet*_xy' )        ! 2d-array
           grid = 'zu1'
           two_d = .TRUE.
           IF ( av == 0 )  THEN
              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    local_pf(i,j,nzb+1) = pet(j,i)
                 ENDDO
              ENDDO
           ELSE
              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    local_pf(i,j,nzb+1) = pet_av(j,i)
                 ENDDO
              ENDDO
           ENDIF

!
!--    Before data is transfered to local_pf, transfer is in 2D dummy variable and exchange ghost
!--    points therein. However, at this point this is only required for instantaneous arrays,
!--    time-averaged quantities are already exchanged.
       CASE ( 'uvem_vitd3*_xy' )        ! 2d-array
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = vitd3_exposure(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'

       CASE ( 'uvem_vitd3dose*_xy' )        ! 2d-array
          IF ( av == 1 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   local_pf(i,j,nzb+1) = vitd3_dose(j,i)
                ENDDO
             ENDDO
          ENDIF

          two_d = .TRUE.
          grid = 'zu1'


       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT


 END SUBROUTINE bio_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 3D output variables (dummy, always 2d!)
!> data_output_3d 709ff
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )

    USE indices

    USE kinds


    IMPLICIT NONE
!
!-- Input variables
    CHARACTER (LEN=*), INTENT(IN) ::  variable   !< Char identifier to select var for output

    INTEGER(iwp), INTENT(IN) ::  av       !< Use averaged data? 0 = no, 1 = yes?
    INTEGER(iwp), INTENT(IN) ::  nzb_do   !< Unused. 2D. nz bottom to nz top
    INTEGER(iwp), INTENT(IN) ::  nzt_do   !< Unused.
!
!-- Output variables
    LOGICAL, INTENT(OUT) ::  found   !< Output found?

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf   !< Temp. result grid to return
!
!-- Internal variables
    INTEGER(iwp) ::  l    !< Running index, radiation grid
    INTEGER(iwp) ::  i    !< Running index, x-dir
    INTEGER(iwp) ::  j    !< Running index, y-dir
    INTEGER(iwp) ::  k    !< Running index, z-dir

!    REAL(wp) ::  mrt  !< Buffer for mean radiant temperature

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

        CASE ( 'bio_mrt' )
            DO  l = 1, nmrtbl
               i = mrtbl(ix,l)
               j = mrtbl(iy,l)
               k = mrtbl(iz,l)
               IF ( k < nzb_do  .OR.  k > nzt_do  .OR.  j < nys  .OR.                              &
                  j > nyn  .OR.  i < nxl  .OR.  i > nxr )  CYCLE
               IF ( av == 0 )  THEN
                  IF ( mrt_include_sw )  THEN
                     local_pf(i,j,k) = ( ( human_absorb * mrtinsw(l) + mrtinlw(l) ) /              &
                                         ( human_emiss * sigma_sb )                                &
                                       )**0.25_wp - degc_to_k
                  ELSE
                     local_pf(i,j,k) = ( mrtinlw(l) / ( human_emiss * sigma_sb ) )**0.25_wp -      &
                                       degc_to_k
                  ENDIF
               ELSE
                  local_pf(i,j,k) = mrt_av_grid(k,j,i)
               ENDIF
            ENDDO

       CASE DEFAULT
          found = .FALSE.

    END SELECT

 END SUBROUTINE bio_data_output_3d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf_interface_mod.
!> netcdf_interface_mod 918ff
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE
!
!-- Input variables
    CHARACTER (LEN=*), INTENT(IN)  ::  var      !< Name of output variable
!
!-- Output variables
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x   !< x grid of output variable
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y   !< y grid of output variable
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z   !< z grid of output variable

    LOGICAL, INTENT(OUT)           ::  found    !< Flag if output var is found
!
!-- Local variables
    INTEGER(iwp) :: l     !< Length of the var array

    LOGICAL      :: is2d  !< Var is 2d?

    found  = .FALSE.
    grid_x = 'none'
    grid_y = 'none'
    grid_z = 'none'

    l = MAX( 2, LEN_TRIM( var ) )
    is2d = ( var(l-1:l) == 'xy' )

    IF ( var(1:4) == 'bio_' )  THEN
       found  = .TRUE.
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'zu'
       IF ( is2d  .AND.  var(1:7) /= 'bio_mrt' )  grid_z = 'zu1'
    ENDIF

    IF ( is2d  .AND.  var(1:4) == 'uvem' )  THEN
       grid_x = 'x'
       grid_y = 'y'
       grid_z = 'zu1'
    ENDIF

 END SUBROUTINE bio_define_netcdf_grid

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Header output for biom module
!> header 982
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_header( io )

    IMPLICIT NONE
!
!-- Input variables
    INTEGER(iwp), INTENT(IN) ::  io           !< Unit of the output file
!
!-- Internal variables
    CHARACTER (LEN=86) ::  output_height_chr  !< String for output height

    WRITE( output_height_chr, '(F8.1,7X)' )  bio_output_height
!
!-- Write biom header
    WRITE( io, 1 )
    WRITE( io, 2 )  TRIM( output_height_chr )
    WRITE( io, 3 )  TRIM( ACHAR( bio_cell_level ) )

1   FORMAT (//' Human thermal comfort module information:'/                                        &
              ' ------------------------------'/)
2   FORMAT ('    --> All indices calculated for a height of (m): ', A )
3   FORMAT ('    --> This corresponds to cell level : ', A )

 END SUBROUTINE bio_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the HTCM
!> init_3d_model 1987ff
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_init

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  netcdf_data_input_uvem

    IMPLICIT NONE
!
!-- Internal vriables
    REAL ( wp )  :: height  !< current height in meters

    IF ( debug_output )  CALL debug_message( 'bio_init', 'start' )
!
!-- Determine cell level corresponding to 1.1 m above ground level (gravimetric center of sample
!-- human)

    bio_cell_level = 0
    bio_output_height = 0.5_wp * dz(1)
    height = 0.0_wp

    bio_cell_level = INT( 1.099_wp / dz(1) )
    bio_output_height = bio_output_height + bio_cell_level * dz(1)
!
!-- Set radiation level if not done by user
    IF ( mrt_minlevel == 0  .AND.  mrt_nlevels == 0 )  THEN
       mrt_minlevel = bio_cell_level
       mrt_nlevels = 1
    ELSE
       IF ( bio_cell_level < mrt_minlevel  .OR.  bio_cell_level > mrt_minlevel + mrt_nlevels - 1 )  &
       THEN
          WRITE( message_string, * ) 'combination of mrt_minlevel = ', mrt_minlevel,               &
                                     ' and mrt_nlevels = ', mrt_nlevels, ' does not include',      &
                                     ' bio_cell_level (gravimetric center of the sample human',    &
                                     ' body 1,1 m)'
          CALL message( 'bio_init', 'BIO0008', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF
!
!-- Init UVEM and load lookup tables
    IF ( uv_exposure )  CALL netcdf_data_input_uvem

!
!-- Check parameters
!-- WARNING This is a WORKAROUND! Due to the design of the module, checks are called at this point
!-- rather than within module_interface_check_parameters.
    CALL bio_check_parameters

    IF ( debug_output )  CALL debug_message( 'bio_init', 'end' )

 END SUBROUTINE bio_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Checks done after the Initialization
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_init_checks

    USE control_parameters,                                                                        &
        ONLY: message_string

    IF ( (.NOT. radiation_interactions) .AND. ( thermal_comfort ) )  THEN
       message_string = 'thermal comfort requires radiation_interaction = .TRUE.'
       CALL message( 'bio_init_checks', 'BIO0009', 1, 2, 0, 6, 0 )
    ENDIF


 END SUBROUTINE bio_init_checks


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &biometeorology_parameters for reading biomet parameters
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_parin

    IMPLICIT NONE

!
!-- Internal variables
    CHARACTER (LEN=100) ::  line  !< Dummy string for current line in parameter file

    INTEGER(iwp) ::  io_status   !< Status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /biometeorology_parameters/  clothing,                                                &
                                          consider_obstructions,                                   &
                                          orientation_angle,                                       &
                                          sun_in_south,                                            &
                                          switch_off_module,                                       &
                                          thermal_comfort,                                         &
                                          turn_to_sun,                                             &
                                          uv_exposure

!
!-- Move to the beginning of the namelist file and try to find and read the namelist named
!-- biometeorology_parameters.
    REWIND( 11 )
    READ( 11, biometeorology_parameters, IOSTAT=io_status )
!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    biometeorology_parameters namelist was found and read correctly. Set flag that
!--    biometeorology_mod is switched on.
       IF ( .NOT. switch_off_module )  biometeorology = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    biometeorology_parameters namelist was found, but contained errors. Print an error message
!--    containing the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'biometeorology_parameters', line )

    ENDIF

 END SUBROUTINE bio_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_rrd_global_ftn( found )

    USE control_parameters,                                                                        &
        ONLY:  length,                                                                             &
               restart_string


    IMPLICIT NONE

    LOGICAL, INTENT(OUT) ::  found      !< variable found? yes = .T., no = .F.

    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

!
!--    Read control flags to determine if input grids need to be averaged.
       CASE ( 'do_average_theta' )
          READ ( 13 )  do_average_theta

       CASE ( 'do_average_q' )
          READ ( 13 )  do_average_q

       CASE ( 'do_average_u' )
          READ ( 13 )  do_average_u

       CASE ( 'do_average_v' )
          READ ( 13 )  do_average_v

       CASE ( 'do_average_w' )
          READ ( 13 )  do_average_w

       CASE ( 'do_average_mrt' )
          READ ( 13 )  do_average_mrt

!
!--    Read control flags to determine which thermal index needs to trigger averaging.
       CASE ( 'average_trigger_perct' )
          READ ( 13 )  average_trigger_perct

       CASE ( 'average_trigger_utci' )
          READ ( 13 )  average_trigger_utci

       CASE ( 'average_trigger_pet' )
          READ ( 13 )  average_trigger_pet

       CASE ( 'average_trigger_mrt' )
          READ ( 13 )  average_trigger_mrt


       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE bio_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_rrd_global_mpi


!
!--  Read control flags to determine if input grids need to be averaged
     CALL rrd_mpi_io( 'do_average_theta', do_average_theta )
     CALL rrd_mpi_io( 'do_average_q', do_average_q )
     CALL rrd_mpi_io( 'do_average_u', do_average_u )
     CALL rrd_mpi_io( 'do_average_v', do_average_v )
     CALL rrd_mpi_io( 'do_average_w', do_average_w )
     CALL rrd_mpi_io( 'do_average_mrt', do_average_mrt )
!
!--  Rad control flags to determine which thermal index needs to trigger averaging
     CALL rrd_mpi_io( 'average_trigger_perct', average_trigger_perct )
     CALL rrd_mpi_io( 'average_trigger_utci', average_trigger_utci )
     CALL rrd_mpi_io( 'average_trigger_pet', average_trigger_pet )
     CALL rrd_mpi_io( 'average_trigger_mrt', average_trigger_mrt )

 END SUBROUTINE bio_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_rrd_local_ftn( found )


    USE control_parameters,                                                                        &
        ONLY:  length,                                                                             &
               restart_string


    IMPLICIT NONE


    LOGICAL, INTENT(OUT) ::  found      !< variable found? yes = .TRUE., no = .FALSE.

    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'mrt_av_grid' )
          IF ( .NOT. ALLOCATED( mrt_av_grid ) )  THEN
             ALLOCATE( mrt_av_grid(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
             mrt_av_grid = 0.0_wp
          ENDIF
          READ ( 13 )  mrt_av_grid

       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE bio_rrd_local_ftn


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!------------------------------------------------------------------------------!
 SUBROUTINE bio_rrd_local_mpi

    USE control_parameters

    USE indices

    USE kinds


    IMPLICIT NONE

    LOGICAL      ::  array_found  !<

!
!-- Note, restart input of time-averaged quantities is skipped in case of cyclic-fill
!-- initialization. This case, input of time-averaged data is useless and can lead to faulty
!-- averaging.
    IF ( .NOT. cyclic_fill_initialization )  THEN
       CALL rd_mpi_io_check_array( 'mrt_av_grid' , found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( mrt_av_grid ) )                                                    &
             ALLOCATE( mrt_av_grid(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          CALL rrd_mpi_io( 'mrt_av_grid', mrt_av_grid )
       ENDIF
    ENDIF

 END SUBROUTINE bio_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write global restart data for the biometeorology module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_wrd_global

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL wrd_write_string( 'do_average_theta' )
       WRITE ( 14 )  do_average_theta
       CALL wrd_write_string( 'do_average_q' )
       WRITE ( 14 )  do_average_q
       CALL wrd_write_string( 'do_average_u' )
       WRITE ( 14 )  do_average_u
       CALL wrd_write_string( 'do_average_v' )
       WRITE ( 14 )  do_average_v
       CALL wrd_write_string( 'do_average_w' )
       WRITE ( 14 )  do_average_w
       CALL wrd_write_string( 'do_average_mrt' )
       WRITE ( 14 )  do_average_mrt
       CALL wrd_write_string( 'average_trigger_perct' )
       WRITE ( 14 )  average_trigger_perct
       CALL wrd_write_string( 'average_trigger_utci' )
       WRITE ( 14 )  average_trigger_utci
       CALL wrd_write_string( 'average_trigger_pet' )
       WRITE ( 14 )  average_trigger_pet
       CALL wrd_write_string( 'average_trigger_mrt' )
       WRITE ( 14 )  average_trigger_mrt

    ELSEIF ( TRIM( restart_data_format_output(1:3) ) == 'mpi' )  THEN

       CALL wrd_mpi_io( 'do_average_theta', do_average_theta )
       CALL wrd_mpi_io( 'do_average_q', do_average_q )
       CALL wrd_mpi_io( 'do_average_u', do_average_u )
       CALL wrd_mpi_io( 'do_average_v', do_average_v )
       CALL wrd_mpi_io( 'do_average_w', do_average_w )
       CALL wrd_mpi_io( 'do_average_mrt', do_average_mrt )
       CALL wrd_mpi_io( 'average_trigger_perct', average_trigger_perct )
       CALL wrd_mpi_io( 'average_trigger_utci', average_trigger_utci )
       CALL wrd_mpi_io( 'average_trigger_pet', average_trigger_pet )
       CALL wrd_mpi_io( 'average_trigger_mrt', average_trigger_mrt )

    ENDIF

 END SUBROUTINE bio_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write local restart data for the biometeorology module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_wrd_local

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       IF ( ALLOCATED( mrt_av_grid ) )  THEN
          CALL wrd_write_string( 'mrt_av_grid' )
          WRITE ( 14 )  mrt_av_grid
       ENDIF

    ELSEIF ( TRIM( restart_data_format_output(1:3) ) == 'mpi' )  THEN

       IF ( ALLOCATED( mrt_av_grid ) )  CALL wrd_mpi_io( 'mrt_av_grid', mrt_av_grid )

    ENDIF

 END SUBROUTINE bio_wrd_local

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate biometeorology MRT for all 2D grid
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_calculate_mrt_grid ( av )

    IMPLICIT NONE

    LOGICAL, INTENT(IN)         ::  av    !< use averaged input?
!
!-- Internal variables
    INTEGER(iwp)                ::  i     !< Running index, x-dir, radiation coordinates
    INTEGER(iwp)                ::  j     !< Running index, y-dir, radiation coordinates
    INTEGER(iwp)                ::  k     !< Running index, y-dir, radiation coordinates
    INTEGER(iwp)                ::  l     !< Running index, radiation coordinates


!
!-- We need to differentiate if averaged input is desired (av == .TRUE.) or not.
    IF ( av )  THEN
!
!--    Make sure tmrt_av_grid is present and set initial value
       IF ( .NOT. ALLOCATED( tmrt_av_grid ) )  THEN
          ALLOCATE( tmrt_av_grid(nys:nyn,nxl:nxr) )
       ENDIF
       tmrt_av_grid = bio_initial_value

!
!--    mrt_av_grid should always be allcoated here, but better make sure ist actually is.
       IF ( ALLOCATED( mrt_av_grid ) )  THEN
!
!--       Iterate over the radiation grid (radiation coordinates) and fill the tmrt_av_grid
!--       (x, y coordinates) where appropriate: tmrt_av_grid is written for all i / j if level (k)
!--       matches output height.
          DO  l = 1, nmrtbl
             i = mrtbl(ix,l)
             j = mrtbl(iy,l)
             k = mrtbl(iz,l)
             IF ( k - topo_top_ind(j,i,0) == bio_cell_level + 1 )  THEN
!
!--             Averaging was done before, so we can just copy the result here.
                tmrt_av_grid(j,i) = mrt_av_grid(k,j,i)

             ENDIF
          ENDDO
       ENDIF

!
!-- In case instantaneous input is desired, mrt values will be re-calculated.
    ELSE
!
!--    Calculate biometeorology MRT from local radiation fluxes calculated by RTM and assign into 2D
!--    grid. Depending on selected output quantities, tmrt_grid might not have been allocated in
!--    bio_check_data_output yet.
       IF ( .NOT. ALLOCATED( tmrt_grid ) )  THEN
          ALLOCATE( tmrt_grid (nys:nyn,nxl:nxr) )
       ENDIF
       tmrt_grid = bio_initial_value

       DO  l = 1, nmrtbl
          i = mrtbl(ix,l)
          j = mrtbl(iy,l)
          k = mrtbl(iz,l)
          IF ( k - topo_top_ind(j,i,0) == bio_cell_level + 1 )  THEN
             IF ( mrt_include_sw )  THEN
                tmrt_grid(j,i) = ( ( human_absorb * mrtinsw(l) +                                   &
                                 mrtinlw(l) )  /                                                   &
                                 ( human_emiss * sigma_sb ) )**0.25_wp -                           &
                                 degc_to_k
             ELSE
                tmrt_grid(j,i) = ( mrtinlw(l)  /                                                   &
                                 ( human_emiss * sigma_sb ) )**0.25_wp -                           &
                                 degc_to_k
             ENDIF
          ENDIF
       ENDDO
    ENDIF

 END SUBROUTINE bio_calculate_mrt_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate static thermal indices for 2D grid point i, j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_get_thermal_index_input_ij( average_input, i, j, ta, vp, ws, pair, tmrt )

    IMPLICIT NONE
!
!-- Input variables
    LOGICAL,      INTENT ( IN ) ::  average_input  !< Determine averaged input conditions?
    INTEGER(iwp), INTENT ( IN ) ::  i     !< Running index, x-dir
    INTEGER(iwp), INTENT ( IN ) ::  j     !< Running index, y-dir
!
!-- Output parameters
    REAL(wp), INTENT ( OUT )    ::  pair  !< Air pressure                    (hPa)
    REAL(wp), INTENT ( OUT )    ::  ta    !< Air temperature                 (degree_C)
    REAL(wp), INTENT ( OUT )    ::  tmrt  !< Mean radiant temperature        (degree_C)
    REAL(wp), INTENT ( OUT )    ::  vp    !< Vapour pressure                 (hPa)
    REAL(wp), INTENT ( OUT )    ::  ws    !< Wind speed    (local level)     (m/s)
!
!-- Internal variables
    INTEGER(iwp)                ::  k     !< Running index, z-dir
    INTEGER(iwp)                ::  k_wind  !< Running index, z-dir, wind speed only

    REAL(wp)                    ::  vp_sat  !< Saturation vapor pressure     (hPa)

!
!-- Determine cell level closest to 1.1m above ground by making use of truncation due to int cast.
    k = INT( topo_top_ind(j,i,0) + bio_cell_level )  !< Vertical cell center closest to 1.1m

!
!-- Avoid non-representative horizontal u and v of 0.0 m/s too close to ground
    k_wind = k
    IF ( bio_cell_level < 1 )  THEN
       k_wind = k + 1
    ENDIF
!
!-- Determine local values:
    IF ( average_input )  THEN
!
!--    Calculate ta from Tp assuming dry adiabatic laps rate
       ta = bio_initial_value
       IF ( ALLOCATED( pt_av ) )  THEN
          ta = pt_av(k,j,i) - ( 0.0098_wp * dz(1) * ( k + 0.5_wp ) ) - degc_to_k
       ENDIF

       vp = bio_initial_value
       IF ( humidity  .AND.  ALLOCATED( q_av ) )  THEN
          vp = q_av(k,j,i)
       ENDIF

       ws = bio_initial_value
       IF ( ALLOCATED( u_av )  .AND.  ALLOCATED( v_av )  .AND.                                     &
          ALLOCATED( w_av ) )  THEN
             ws = ( 0.5_wp * ABS( u_av(k_wind,j,i) + u_av(k_wind,j,i+1) ) +                        &
                    0.5_wp * ABS( v_av(k_wind,j,i) + v_av(k_wind,j+1,i) ) +                        &
                    0.5_wp * ABS( w_av(k_wind,j,i) + w_av(k_wind+1,j,i) ) )
       ENDIF
    ELSE
!
!--    Calculate ta from Tp assuming dry adiabatic laps rate
       ta = pt(k,j,i) - ( 0.0098_wp * dz(1) * (  k + 0.5_wp ) ) - degc_to_k

       vp = bio_initial_value
       IF ( humidity )  THEN
          vp = q(k,j,i)
       ENDIF

       ws = ( 0.5_wp * ABS( u(k_wind,j,i) + u(k_wind,j,i+1) )  +                                   &
              0.5_wp * ABS( v(k_wind,j,i) + v(k_wind,j+1,i) )  +                                   &
              0.5_wp * ABS( w(k_wind,j,i) + w(k_wind+1,j,i) ) )

    ENDIF
!
!-- Local air pressure
    pair = surface_pressure
!
!-- Calculate water vapour pressure at saturation and convert to hPa.
!-- The magnus formula is limited to temperatures up to 333.15 K to avoid negative values of vp_sat.
    IF ( vp > -998.0_wp )  THEN
       vp_sat = 0.01_wp * magnus( MIN( ta + degc_to_k, 333.15_wp ) )
       vp  = vp * pair / ( vp + 0.622_wp )
       IF ( vp > vp_sat )  vp = vp_sat
    ENDIF
!
!-- Local mtr value at [i,j]. Inital values (fill values) will be kept inside terrain/buildings.
    tmrt = bio_initial_value
    IF ( .NOT. average_input )  THEN
!
!--    Use MRT from RTM precalculated in tmrt_grid
       tmrt = tmrt_grid(j,i)
    ELSE
       tmrt = tmrt_av_grid(j,i)
    ENDIF

 END SUBROUTINE bio_get_thermal_index_input_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate static thermal indices for any point within a 2D grid time_integration.f90: 1065ff
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_calculate_thermal_index_maps( av )

    IMPLICIT NONE
!
!-- Input attributes
    LOGICAL, INTENT ( IN ) ::  av  !< Calculate based on averaged input conditions?
!
!-- Internal variables
    INTEGER(iwp) ::  i, j     !< Running index

    REAL(wp) ::  clo          !< Clothing index                (no dimension)
    REAL(wp) ::  pair         !< Air pressure                     (hPa)
    REAL(wp) ::  perct_ij     !< Perceived temperature            (degree_C)
    REAL(wp) ::  pet_ij       !< Physiologically equivalent temperature  (degree_C)
    REAL(wp) ::  ta           !< Air temperature                  (degree_C)
    REAL(wp) ::  tmrt_ij      !< Mean radiant temperature         (degree_C)
    REAL(wp) ::  utci_ij      !< Universal thermal climate index  (degree_C)
    REAL(wp) ::  vp           !< Vapour pressure                  (hPa)
    REAL(wp) ::  ws           !< Wind speed    (local level)      (m/s)

!
!-- Check if some thermal index is desired. Don't do anything if, e.g. only bio_mrt is desired.
    IF ( do_calculate_perct    .OR.  do_calculate_perct_av  .OR.  do_calculate_utci    .OR.        &
         do_calculate_utci_av  .OR.  do_calculate_pet       .OR.  do_calculate_pet_av  .OR.        &
         do_calculate_mrt2d )  THEN

!
!--    fill out the MRT 2D grid from appropriate source (RTM, RRTMG,...)
       CALL bio_calculate_mrt_grid ( av )

       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Determine local input conditions
             tmrt_ij = bio_initial_value
             vp      = bio_initial_value
!
!--          Determine local meteorological conditions
             CALL bio_get_thermal_index_input_ij ( av, i, j, ta, vp, ws, pair, tmrt_ij )
!
!--          Only proceed if input is available
             pet_ij   = bio_initial_value
             perct_ij = bio_initial_value
             utci_ij  = bio_initial_value
             IF ( .NOT. ( tmrt_ij <= -998.0_wp  .OR.  vp <= -998.0_wp  .OR.                        &
                          ws <= -998.0_wp  .OR.  ta <= -998.0_wp ) )  THEN
!
!--             Calculate static thermal indices based on local tmrt
                clo = bio_initial_value

                IF ( do_calculate_perct  .OR.  do_calculate_perct_av )  THEN
!
!--                Estimate local perceived temperature
                   CALL calculate_perct_static( ta, vp, ws, tmrt_ij, pair, clo, perct_ij )
                ENDIF

                IF ( do_calculate_utci  .OR.  do_calculate_utci_av )  THEN
!
!--                Estimate local universal thermal climate index
                   CALL calculate_utci_static( ta, vp, ws, tmrt_ij, bio_output_height, utci_ij )
                ENDIF

                IF ( do_calculate_pet  .OR.  do_calculate_pet_av )  THEN
!
!--                Estimate local physiologically equivalent temperature
                   CALL calculate_pet_static( ta, vp, ws, tmrt_ij, pair, pet_ij )
                ENDIF
             ENDIF


             IF ( av )  THEN
!
!--             Write results for selected averaged indices
                IF ( do_calculate_perct_av )  THEN
                   perct_av(j, i) = perct_ij
                ENDIF
                IF ( do_calculate_utci_av )  THEN
                   utci_av(j, i) = utci_ij
                ENDIF
                IF ( do_calculate_pet_av )  THEN
                   pet_av(j, i)  = pet_ij
                ENDIF
             ELSE
!
!--             Write result for selected indices
                IF ( do_calculate_perct )  THEN
                   perct(j, i) = perct_ij
                ENDIF
                IF ( do_calculate_utci )  THEN
                   utci(j, i) = utci_ij
                ENDIF
                IF ( do_calculate_pet )  THEN
                   pet(j, i)  = pet_ij
                ENDIF
             ENDIF

          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE bio_calculate_thermal_index_maps

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate dynamic thermal indices (currently only iPT, but expandable)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_calc_ipt( ta, vp, ws, pair, tmrt, dt, energy_storage, t_clo, clo, actlev, age,     &
                          weight, height, work, sex, ipt )

    IMPLICIT NONE
!
!-- Input parameters
    INTEGER(iwp), INTENT ( IN ) ::  sex  !< Sex of agent (1 = male, 2 = female)

    REAL(wp), INTENT ( IN )  ::  age     !< Age of agent                     (y)
    REAL(wp), INTENT ( IN )  ::  dt      !< Time past since last calculation (s)
    REAL(wp), INTENT ( IN )  ::  height  !< Height of agent                  (m)
    REAL(wp), INTENT ( IN )  ::  pair    !< Air pressure                     (hPa)
    REAL(wp), INTENT ( IN )  ::  ta      !< Air temperature                  (degree_C)
    REAL(wp), INTENT ( IN )  ::  tmrt    !< Mean radiant temperature         (degree_C)
    REAL(wp), INTENT ( IN )  ::  vp      !< Vapour pressure                  (hPa)
    REAL(wp), INTENT ( IN )  ::  weight  !< Weight of agent                  (Kg)
    REAL(wp), INTENT ( IN )  ::  work    !< Mechanical workload of agent  (without metabolism!) (W)
    REAL(wp), INTENT ( IN )  ::  ws      !< Wind speed  (local level)        (m/s)

!
!-- Both, input and output parameters
    Real(wp), INTENT ( INOUT )  ::  actlev            !< Individuals activity level
                                                      !< per unit surface area      (W/m²)
    Real(wp), INTENT ( INOUT )  ::  clo               !< Current clothing in sulation
    Real(wp), INTENT ( INOUT )  ::  energy_storage    !< Energy storage   (W/m²)
    Real(wp), INTENT ( INOUT )  ::  t_clo             !< Clothing temperature       (degree_C)
!
!-- Output parameters
    REAL(wp), INTENT ( OUT ) ::  ipt    !< Instationary perceived temp.   (degree_C)
!
!-- Return immediatelly if nothing to do!
    IF ( .NOT. thermal_comfort )  THEN
        RETURN
    ENDIF
!
!-- If clo equals the initial value, this is the initial call
    IF ( clo <= -998.0_wp )  THEN
!
!--    Initialize nonstationary perceived temperature with personalized PT as an initial guess, set
!--    actlev and clo
       CALL ipt_init( age, weight, height, sex, work, actlev, clo, ta, vp, ws, tmrt, pair, dt,     &
                      energy_storage, t_clo, ipt )
    ELSE
!
!--    Estimate local nonstatinoary perceived temperature
       CALL ipt_cycle ( ta, vp, ws, tmrt, pair, dt, energy_storage, t_clo, clo, actlev, work, ipt )
    ENDIF

 END SUBROUTINE bio_calc_ipt



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> SUBROUTINE for calculating UTCI Temperature (UTCI)
!> computed by a 6th order approximation
!>
!> UTCI regression equation according to
!> Bröde P, Fiala D, Blazejczyk K, Holmér I, Jendritzky G, Kampmann B, Tinz B, Havenith G (2012)
!> Deriving the operational procedure for the Universal Thermal Climate Index (UTCI). International
!> Journal of Biometeorology 56 (3):481-494. doi:10.1007/s00484-011-0454-1
!>
!> original source available at:
!> www.utci.org
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calculate_utci_static( ta_in, vp, ws_hag, tmrt, hag, utci_ij )

    IMPLICIT NONE
!
!-- Type of input of the argument list
    REAL(WP), INTENT ( IN )  ::  hag      !< Height of wind speed input (m)
    REAL(WP), INTENT ( IN )  ::  ta_in    !< Local air temperature (degree_C)
    REAL(WP), INTENT ( IN )  ::  tmrt     !< Local mean radiant temperature (degree_C)
    REAL(WP), INTENT ( IN )  ::  vp       !< Loacl vapour pressure (hPa)
    REAL(WP), INTENT ( IN )  ::  ws_hag   !< Incident wind speed (m/s)
!
!-- Type of output of the argument list
    REAL(WP) ::  d_tmrt       !< delta-tmrt               (degree_C)
    REAL(WP) ::  d_tmrt2      !< 2 times d_tmrt
    REAL(WP) ::  d_tmrt3      !< 3 times d_tmrt
    REAL(WP) ::  d_tmrt4      !< 4 times d_tmrt
    REAL(WP) ::  d_tmrt5      !< 5 times d_tmrt
    REAL(WP) ::  d_tmrt6      !< 6 times d_tmrt
    REAL(WP) ::  offset       !< utci deviation by ta cond. exceeded      (degree_C)
    REAL(WP) ::  pa           !< air pressure in kPa      (kPa)
    REAL(WP) ::  pa2          !< 2 times pa
    REAL(WP) ::  pa3          !< 3 times pa
    REAL(WP) ::  pa4          !< 4 times pa
    REAL(WP) ::  pa5          !< 5 times pa
    REAL(WP) ::  pa6          !< 6 times pa
    REAL(WP) ::  part_d_tmrt  !< Mean radiant temp. related part of the reg.
    REAL(WP) ::  part_pa      !< Air pressure related part of the regression
    REAL(WP) ::  part_pa2     !< Air pressure^2 related part of the regression
    REAL(WP) ::  part_pa3     !< Air pressure^3 related part of the regression
    REAL(WP) ::  part_pa46    !< Air pressure^4-6 related part of the regression
    REAL(WP) ::  part_ta      !< Air temperature related part of the regression
    REAL(WP) ::  part_va      !< Vapour pressure related part of the regression
    REAL(WP) ::  ta           !< air temperature modified by offset (degree_C)
    REAL(WP) ::  ta2          !< 2 times ta
    REAL(WP) ::  ta3          !< 3 times ta
    REAL(WP) ::  ta4          !< 4 times ta
    REAL(WP) ::  ta5          !< 5 times ta
    REAL(WP) ::  ta6          !< 6 times ta
    REAL(WP) ::  va           !< wind speed at 10 m above ground level    (m/s)
    REAL(WP) ::  va2          !< 2 times va
    REAL(WP) ::  va3          !< 3 times va
    REAL(WP) ::  va4          !< 4 times va
    REAL(WP) ::  va5          !< 5 times va
    REAL(WP) ::  va6          !< 6 times va


    REAL(wp), INTENT ( OUT ) ::  utci_ij  !< Universal Thermal Climate Index (degree_C)

!
!-- Initialize
    offset = 0.0_wp
    ta = ta_in
    d_tmrt = tmrt - ta_in
!
!-- Use vapour pressure in kpa
    pa = vp / 10.0_wp
!
!-- Wind altitude correction from hag to 10m after Broede et al. (2012), eq.3
!-- z(0) is set to 0.01 according to UTCI profile definition
    va = ws_hag *  log ( 10.0_wp / 0.01_wp ) / log ( hag / 0.01_wp )
!
!-- Check if input values in range after Broede et al. (2012)
    IF ( ( d_tmrt > 70.0_wp )  .OR.  ( d_tmrt < -30.0_wp )  .OR.  ( vp >= 50.0_wp ) )  THEN
       utci_ij = bio_initial_value
       RETURN
    ENDIF
!
!-- Apply eq. 2 in Broede et al. (2012) for ta out of bounds
    IF ( ta > 50.0_wp )  THEN
       offset = ta - 50.0_wp
       ta = 50.0_wp
    ENDIF
    IF ( ta < -50.0_wp )  THEN
       offset = ta + 50.0_wp
       ta = -50.0_wp
    ENDIF
!
!-- For routine application. For wind speeds and relative humidity values below 0.5 m/s or 5%,
!-- respectively, the user is advised to use the lower bounds for the calculations.
    IF ( va < 0.5_wp )  va = 0.5_wp
    IF ( va > 17.0_wp )  va = 17.0_wp

!
!-- Pre-calculate multiples of input parameters to save time later
    ta2 = ta  * ta
    ta3 = ta2 * ta
    ta4 = ta3 * ta
    ta5 = ta4 * ta
    ta6 = ta5 * ta

    va2 = va  * va
    va3 = va2 * va
    va4 = va3 * va
    va5 = va4 * va
    va6 = va5 * va

    d_tmrt2 = d_tmrt  * d_tmrt
    d_tmrt3 = d_tmrt2 * d_tmrt
    d_tmrt4 = d_tmrt3 * d_tmrt
    d_tmrt5 = d_tmrt4 * d_tmrt
    d_tmrt6 = d_tmrt5 * d_tmrt

    pa2 = pa  * pa
    pa3 = pa2 * pa
    pa4 = pa3 * pa
    pa5 = pa4 * pa
    pa6 = pa5 * pa

!
!-- Pre-calculate parts of the regression equation
    part_ta = (  6.07562052e-01_wp )       +                                                       &
              ( -2.27712343e-02_wp ) * ta  +                                                       &
              (  8.06470249e-04_wp ) * ta2 +                                                       &
              ( -1.54271372e-04_wp ) * ta3 +                                                       &
              ( -3.24651735e-06_wp ) * ta4 +                                                       &
              (  7.32602852e-08_wp ) * ta5 +                                                       &
              (  1.35959073e-09_wp ) * ta6

    part_va = ( -2.25836520e+00_wp ) *       va  +                                                 &
              (  8.80326035e-02_wp ) * ta  * va  +                                                 &
              (  2.16844454e-03_wp ) * ta2 * va  +                                                 &
              ( -1.53347087e-05_wp ) * ta3 * va  +                                                 &
              ( -5.72983704e-07_wp ) * ta4 * va  +                                                 &
              ( -2.55090145e-09_wp ) * ta5 * va  +                                                 &
              ( -7.51269505e-01_wp ) *       va2 +                                                 &
              ( -4.08350271e-03_wp ) * ta  * va2 +                                                 &
              ( -5.21670675e-05_wp ) * ta2 * va2 +                                                 &
              (  1.94544667e-06_wp ) * ta3 * va2 +                                                 &
              (  1.14099531e-08_wp ) * ta4 * va2 +                                                 &
              (  1.58137256e-01_wp ) *       va3 +                                                 &
              ( -6.57263143e-05_wp ) * ta  * va3 +                                                 &
              (  2.22697524e-07_wp ) * ta2 * va3 +                                                 &
              ( -4.16117031e-08_wp ) * ta3 * va3 +                                                 &
              ( -1.27762753e-02_wp ) *       va4 +                                                 &
              (  9.66891875e-06_wp ) * ta  * va4 +                                                 &
              (  2.52785852e-09_wp ) * ta2 * va4 +                                                 &
              (  4.56306672e-04_wp ) *       va5 +                                                 &
              ( -1.74202546e-07_wp ) * ta  * va5 +                                                 &
              ( -5.91491269e-06_wp ) * va6

    part_d_tmrt = (  3.98374029e-01_wp ) *             d_tmrt  +                                   &
                  (  1.83945314e-04_wp ) * ta  *       d_tmrt  +                                   &
                  ( -1.73754510e-04_wp ) * ta2 *       d_tmrt  +                                   &
                  ( -7.60781159e-07_wp ) * ta3 *       d_tmrt  +                                   &
                  (  3.77830287e-08_wp ) * ta4 *       d_tmrt  +                                   &
                  (  5.43079673e-10_wp ) * ta5 *       d_tmrt  +                                   &
                  ( -2.00518269e-02_wp ) *       va  * d_tmrt  +                                   &
                  (  8.92859837e-04_wp ) * ta  * va  * d_tmrt  +                                   &
                  (  3.45433048e-06_wp ) * ta2 * va  * d_tmrt  +                                   &
                  ( -3.77925774e-07_wp ) * ta3 * va  * d_tmrt  +                                   &
                  ( -1.69699377e-09_wp ) * ta4 * va  * d_tmrt  +                                   &
                  (  1.69992415e-04_wp ) *       va2 * d_tmrt  +                                   &
                  ( -4.99204314e-05_wp ) * ta  * va2 * d_tmrt  +                                   &
                  (  2.47417178e-07_wp ) * ta2 * va2 * d_tmrt  +                                   &
                  (  1.07596466e-08_wp ) * ta3 * va2 * d_tmrt  +                                   &
                  (  8.49242932e-05_wp ) *       va3 * d_tmrt  +                                   &
                  (  1.35191328e-06_wp ) * ta  * va3 * d_tmrt  +                                   &
                  ( -6.21531254e-09_wp ) * ta2 * va3 * d_tmrt  +                                   &
                  ( -4.99410301e-06_wp ) * va4 *       d_tmrt  +                                   &
                  ( -1.89489258e-08_wp ) * ta  * va4 * d_tmrt  +                                   &
                  (  8.15300114e-08_wp ) *       va5 * d_tmrt  +                                   &
                  (  7.55043090e-04_wp ) *             d_tmrt2 +                                   &
                  ( -5.65095215e-05_wp ) * ta  *       d_tmrt2 +                                   &
                  ( -4.52166564e-07_wp ) * ta2 *       d_tmrt2 +                                   &
                  (  2.46688878e-08_wp ) * ta3 *       d_tmrt2 +                                   &
                  (  2.42674348e-10_wp ) * ta4 *       d_tmrt2 +                                   &
                  (  1.54547250e-04_wp ) *       va  * d_tmrt2 +                                   &
                  (  5.24110970e-06_wp ) * ta  * va  * d_tmrt2 +                                   &
                  ( -8.75874982e-08_wp ) * ta2 * va  * d_tmrt2 +                                   &
                  ( -1.50743064e-09_wp ) * ta3 * va  * d_tmrt2 +                                   &
                  ( -1.56236307e-05_wp ) *       va2 * d_tmrt2 +                                   &
                  ( -1.33895614e-07_wp ) * ta  * va2 * d_tmrt2 +                                   &
                  (  2.49709824e-09_wp ) * ta2 * va2 * d_tmrt2 +                                   &
                  (  6.51711721e-07_wp ) *       va3 * d_tmrt2 +                                   &
                  (  1.94960053e-09_wp ) * ta  * va3 * d_tmrt2 +                                   &
                  ( -1.00361113e-08_wp ) *       va4 * d_tmrt2 +                                   &
                  ( -1.21206673e-05_wp ) *             d_tmrt3 +                                   &
                  ( -2.18203660e-07_wp ) * ta  *       d_tmrt3 +                                   &
                  (  7.51269482e-09_wp ) * ta2 *       d_tmrt3 +                                   &
                  (  9.79063848e-11_wp ) * ta3 *       d_tmrt3 +                                   &
                  (  1.25006734e-06_wp ) *       va  * d_tmrt3 +                                   &
                  ( -1.81584736e-09_wp ) * ta  * va  * d_tmrt3 +                                   &
                  ( -3.52197671e-10_wp ) * ta2 * va  * d_tmrt3 +                                   &
                  ( -3.36514630e-08_wp ) *       va2 * d_tmrt3 +                                   &
                  (  1.35908359e-10_wp ) * ta  * va2 * d_tmrt3 +                                   &
                  (  4.17032620e-10_wp ) *       va3 * d_tmrt3 +                                   &
                  ( -1.30369025e-09_wp ) *             d_tmrt4 +                                   &
                  (  4.13908461e-10_wp ) * ta  *       d_tmrt4 +                                   &
                  (  9.22652254e-12_wp ) * ta2 *       d_tmrt4 +                                   &
                  ( -5.08220384e-09_wp ) *       va  * d_tmrt4 +                                   &
                  ( -2.24730961e-11_wp ) * ta  * va  * d_tmrt4 +                                   &
                  (  1.17139133e-10_wp ) *       va2 * d_tmrt4 +                                   &
                  (  6.62154879e-10_wp ) *             d_tmrt5 +                                   &
                  (  4.03863260e-13_wp ) * ta  *       d_tmrt5 +                                   &
                  (  1.95087203e-12_wp ) *       va  * d_tmrt5 +                                   &
                  ( -4.73602469e-12_wp ) *             d_tmrt6

    part_pa = (  5.12733497e+00_wp ) *                       pa +                                  &
              ( -3.12788561e-01_wp ) * ta  *                 pa +                                  &
              ( -1.96701861e-02_wp ) * ta2 *                 pa +                                  &
              (  9.99690870e-04_wp ) * ta3 *                 pa +                                  &
              (  9.51738512e-06_wp ) * ta4 *                 pa +                                  &
              ( -4.66426341e-07_wp ) * ta5 *                 pa +                                  &
              (  5.48050612e-01_wp ) *       va  *           pa +                                  &
              ( -3.30552823e-03_wp ) * ta  * va  *           pa +                                  &
              ( -1.64119440e-03_wp ) * ta2 * va  *           pa +                                  &
              ( -5.16670694e-06_wp ) * ta3 * va  *           pa +                                  &
              (  9.52692432e-07_wp ) * ta4 * va  *           pa +                                  &
              ( -4.29223622e-02_wp ) *       va2 *           pa +                                  &
              (  5.00845667e-03_wp ) * ta  * va2 *           pa +                                  &
              (  1.00601257e-06_wp ) * ta2 * va2 *           pa +                                  &
              ( -1.81748644e-06_wp ) * ta3 * va2 *           pa +                                  &
              ( -1.25813502e-03_wp ) *       va3 *           pa +                                  &
              ( -1.79330391e-04_wp ) * ta  * va3 *           pa +                                  &
              (  2.34994441e-06_wp ) * ta2 * va3 *           pa +                                  &
              (  1.29735808e-04_wp ) *       va4 *           pa +                                  &
              (  1.29064870e-06_wp ) * ta  * va4 *           pa +                                  &
              ( -2.28558686e-06_wp ) *       va5 *           pa +                                  &
              ( -3.69476348e-02_wp ) *             d_tmrt  * pa +                                  &
              (  1.62325322e-03_wp ) * ta  *       d_tmrt  * pa +                                  &
              ( -3.14279680e-05_wp ) * ta2 *       d_tmrt  * pa +                                  &
              (  2.59835559e-06_wp ) * ta3 *       d_tmrt  * pa +                                  &
              ( -4.77136523e-08_wp ) * ta4 *       d_tmrt  * pa +                                  &
              (  8.64203390e-03_wp ) *       va  * d_tmrt  * pa +                                  &
              ( -6.87405181e-04_wp ) * ta  * va  * d_tmrt  * pa +                                  &
              ( -9.13863872e-06_wp ) * ta2 * va  * d_tmrt  * pa +                                  &
              (  5.15916806e-07_wp ) * ta3 * va  * d_tmrt  * pa +                                  &
              ( -3.59217476e-05_wp ) *       va2 * d_tmrt  * pa +                                  &
              (  3.28696511e-05_wp ) * ta  * va2 * d_tmrt  * pa +                                  &
              ( -7.10542454e-07_wp ) * ta2 * va2 * d_tmrt  * pa +                                  &
              ( -1.24382300e-05_wp ) *       va3 * d_tmrt  * pa +                                  &
              ( -7.38584400e-09_wp ) * ta  * va3 * d_tmrt  * pa +                                  &
              (  2.20609296e-07_wp ) *       va4 * d_tmrt  * pa +                                  &
              ( -7.32469180e-04_wp ) *             d_tmrt2 * pa +                                  &
              ( -1.87381964e-05_wp ) * ta  *       d_tmrt2 * pa +                                  &
              (  4.80925239e-06_wp ) * ta2 *       d_tmrt2 * pa +                                  &
              ( -8.75492040e-08_wp ) * ta3 *       d_tmrt2 * pa +                                  &
              (  2.77862930e-05_wp ) *       va  * d_tmrt2 * pa +                                  &
              ( -5.06004592e-06_wp ) * ta  * va  * d_tmrt2 * pa +                                  &
              (  1.14325367e-07_wp ) * ta2 * va  * d_tmrt2 * pa +                                  &
              (  2.53016723e-06_wp ) *       va2 * d_tmrt2 * pa +                                  &
              ( -1.72857035e-08_wp ) * ta  * va2 * d_tmrt2 * pa +                                  &
              ( -3.95079398e-08_wp ) *       va3 * d_tmrt2 * pa +                                  &
              ( -3.59413173e-07_wp ) *             d_tmrt3 * pa +                                  &
              (  7.04388046e-07_wp ) * ta  *       d_tmrt3 * pa +                                  &
              ( -1.89309167e-08_wp ) * ta2 *       d_tmrt3 * pa +                                  &
              ( -4.79768731e-07_wp ) *       va  * d_tmrt3 * pa +                                  &
              (  7.96079978e-09_wp ) * ta  * va  * d_tmrt3 * pa +                                  &
              (  1.62897058e-09_wp ) *       va2 * d_tmrt3 * pa +                                  &
              (  3.94367674e-08_wp ) *             d_tmrt4 * pa +                                  &
              ( -1.18566247e-09_wp ) * ta *        d_tmrt4 * pa +                                  &
              (  3.34678041e-10_wp ) *       va  * d_tmrt4 * pa +                                  &
              ( -1.15606447e-10_wp ) *             d_tmrt5 * pa

    part_pa2 = ( -2.80626406e+00_wp ) *                       pa2 +                                &
               (  5.48712484e-01_wp ) * ta  *                 pa2 +                                &
               ( -3.99428410e-03_wp ) * ta2 *                 pa2 +                                &
               ( -9.54009191e-04_wp ) * ta3 *                 pa2 +                                &
               (  1.93090978e-05_wp ) * ta4 *                 pa2 +                                &
               ( -3.08806365e-01_wp ) *       va *            pa2 +                                &
               (  1.16952364e-02_wp ) * ta  * va *            pa2 +                                &
               (  4.95271903e-04_wp ) * ta2 * va *            pa2 +                                &
               ( -1.90710882e-05_wp ) * ta3 * va *            pa2 +                                &
               (  2.10787756e-03_wp ) *       va2 *           pa2 +                                &
               ( -6.98445738e-04_wp ) * ta  * va2 *           pa2 +                                &
               (  2.30109073e-05_wp ) * ta2 * va2 *           pa2 +                                &
               (  4.17856590e-04_wp ) *       va3 *           pa2 +                                &
               ( -1.27043871e-05_wp ) * ta  * va3 *           pa2 +                                &
               ( -3.04620472e-06_wp ) *       va4 *           pa2 +                                &
               (  5.14507424e-02_wp ) *             d_tmrt  * pa2 +                                &
               ( -4.32510997e-03_wp ) * ta  *       d_tmrt  * pa2 +                                &
               (  8.99281156e-05_wp ) * ta2 *       d_tmrt  * pa2 +                                &
               ( -7.14663943e-07_wp ) * ta3 *       d_tmrt  * pa2 +                                &
               ( -2.66016305e-04_wp ) *       va  * d_tmrt  * pa2 +                                &
               (  2.63789586e-04_wp ) * ta  * va  * d_tmrt  * pa2 +                                &
               ( -7.01199003e-06_wp ) * ta2 * va  * d_tmrt  * pa2 +                                &
               ( -1.06823306e-04_wp ) *       va2 * d_tmrt  * pa2 +                                &
               (  3.61341136e-06_wp ) * ta  * va2 * d_tmrt  * pa2 +                                &
               (  2.29748967e-07_wp ) *       va3 * d_tmrt  * pa2 +                                &
               (  3.04788893e-04_wp ) *             d_tmrt2 * pa2 +                                &
               ( -6.42070836e-05_wp ) * ta  *       d_tmrt2 * pa2 +                                &
               (  1.16257971e-06_wp ) * ta2 *       d_tmrt2 * pa2 +                                &
               (  7.68023384e-06_wp ) *       va  * d_tmrt2 * pa2 +                                &
               ( -5.47446896e-07_wp ) * ta  * va  * d_tmrt2 * pa2 +                                &
               ( -3.59937910e-08_wp ) *       va2 * d_tmrt2 * pa2 +                                &
               ( -4.36497725e-06_wp ) *             d_tmrt3 * pa2 +                                &
               (  1.68737969e-07_wp ) * ta  *       d_tmrt3 * pa2 +                                &
               (  2.67489271e-08_wp ) *       va  * d_tmrt3 * pa2 +                                &
               (  3.23926897e-09_wp ) *             d_tmrt4 * pa2

    part_pa3 = ( -3.53874123e-02_wp ) *                       pa3 +                                &
               ( -2.21201190e-01_wp ) * ta  *                 pa3 +                                &
               (  1.55126038e-02_wp ) * ta2 *                 pa3 +                                &
               ( -2.63917279e-04_wp ) * ta3 *                 pa3 +                                &
               (  4.53433455e-02_wp ) *       va  *           pa3 +                                &
               ( -4.32943862e-03_wp ) * ta  * va  *           pa3 +                                &
               (  1.45389826e-04_wp ) * ta2 * va  *           pa3 +                                &
               (  2.17508610e-04_wp ) *       va2 *           pa3 +                                &
               ( -6.66724702e-05_wp ) * ta  * va2 *           pa3 +                                &
               (  3.33217140e-05_wp ) *       va3 *           pa3 +                                &
               ( -2.26921615e-03_wp ) *             d_tmrt  * pa3 +                                &
               (  3.80261982e-04_wp ) * ta  *       d_tmrt  * pa3 +                                &
               ( -5.45314314e-09_wp ) * ta2 *       d_tmrt  * pa3 +                                &
               ( -7.96355448e-04_wp ) *       va  * d_tmrt  * pa3 +                                &
               (  2.53458034e-05_wp ) * ta  * va  * d_tmrt  * pa3 +                                &
               ( -6.31223658e-06_wp ) *       va2 * d_tmrt  * pa3 +                                &
               (  3.02122035e-04_wp ) *             d_tmrt2 * pa3 +                                &
               ( -4.77403547e-06_wp ) * ta  *       d_tmrt2 * pa3 +                                &
               (  1.73825715e-06_wp ) *       va  * d_tmrt2 * pa3 +                                &
               ( -4.09087898e-07_wp ) *             d_tmrt3 * pa3

    part_pa46 = (  6.14155345e-01_wp ) *                       pa4 +                               &
                ( -6.16755931e-02_wp ) * ta  *                 pa4 +                               &
                (  1.33374846e-03_wp ) * ta2 *                 pa4 +                               &
                (  3.55375387e-03_wp ) *       va  *           pa4 +                               &
                ( -5.13027851e-04_wp ) * ta  * va  *           pa4 +                               &
                (  1.02449757e-04_wp ) *       va2 *           pa4 +                               &
                ( -1.48526421e-03_wp ) *             d_tmrt  * pa4 +                               &
                ( -4.11469183e-05_wp ) * ta  *       d_tmrt  * pa4 +                               &
                ( -6.80434415e-06_wp ) *       va  * d_tmrt  * pa4 +                               &
                ( -9.77675906e-06_wp ) *             d_tmrt2 * pa4 +                               &
                (  8.82773108e-02_wp ) *                       pa5 +                               &
                ( -3.01859306e-03_wp ) * ta  *                 pa5 +                               &
                (  1.04452989e-03_wp ) *       va  *           pa5 +                               &
                (  2.47090539e-04_wp ) *             d_tmrt  * pa5 +                               &
                (  1.48348065e-03_wp ) *                       pa6
!
!-- Calculate 6th order polynomial as approximation
    utci_ij = ta + part_ta + part_va + part_d_tmrt + part_pa + part_pa2 + part_pa3 + part_pa46
!
!-- Consider offset in result
    utci_ij = utci_ij + offset

 END SUBROUTINE calculate_utci_static




!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate_perct_static: Estimation of perceived temperature (PT, degree_C)
!> Value of perct is the Perceived Temperature, degree centigrade
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calculate_perct_static( ta, vp, ws, tmrt, pair, clo, perct_ij )

    IMPLICIT NONE
!
!-- Type of input of the argument list
    REAL(wp), INTENT ( IN )  :: pair !< Local barometric air pressure (hPa)
    REAL(wp), INTENT ( IN )  :: ta   !< Local air temperature (degC)
    REAL(wp), INTENT ( IN )  :: tmrt !< Local mean radiant temperature (degC)
    REAL(wp), INTENT ( IN )  :: vp   !< Local vapour pressure (hPa)
    REAL(wp), INTENT ( IN )  :: ws   !< Local wind velocitry (m/s)
!
!-- Type of output of the argument list
    REAL(wp), INTENT ( OUT ) :: clo       !< Clothing index (dimensionless)
    REAL(wp), INTENT ( OUT ) :: perct_ij  !< Perceived temperature (degC)
!
!-- Parameters for standard "Klima-Michel"
    REAL(wp), PARAMETER :: actlev = 134.6862_wp  !< Workload by activity per standardized surface (A_Du)
    REAL(wp), PARAMETER :: eta = 0.0_wp          !< Mechanical work efficiency for walking on flat
                                                 !< ground (compare to Fanger (1972) pp 24f)
!
!-- Type of program variables
    REAL(wp), PARAMETER :: eps = 0.0005  !< Accuracy in clothing insulation (clo) for evaluation the root of Fanger's PMV (pmva=0)

    INTEGER(iwp) :: ncount      !< running index
    INTEGER(iwp) :: nerr_cold   !< error number (cold conditions)
    INTEGER(iwp) :: nerr        !< error number

    LOGICAL :: sultrieness

    REAL(wp) ::  clon           !< clo for neutral conditions   (clo)
    REAL(wp) ::  d_pmv          !< PMV deviation (dimensionless --> PMV)
    REAL(wp) ::  dgtcm          !< Mean deviation dependent on perct
    REAL(wp) ::  dgtcstd        !< Mean deviation plus its standard deviation
    REAL(wp) ::  d_std          !< factor to threshold for sultriness
    REAL(wp) ::  ireq_minimal   !< Minimal required clothing insulation (clo)
    REAL(wp) ::  pmv_s          !< Fangers predicted mean vote for summer clothing
    REAL(wp) ::  pmv_w          !< Fangers predicted mean vote for winter clothing
    REAL(wp) ::  pmva           !< adjusted predicted mean vote
    REAL(wp) ::  pmvs           !< pred. mean vote considering sultrieness
    REAL(wp) ::  ptc            !< perceived temp. for cold conditions (degree_C)
    REAL(wp) ::  sclo           !< summer clothing insulation
    REAL(wp) ::  svp_ta         !< saturation vapor pressure    (hPa)
    REAL(wp) ::  sult_lim       !< threshold for sultrieness    (hPa)
    REAL(wp) ::  wclo           !< winter clothing insulation

!
!-- Initialise
    perct_ij = bio_initial_value

    nerr     = 0
    ncount   = 0
    sultrieness  = .FALSE.
!
!-- Tresholds: clothing insulation (account for model inaccuracies)
!-- Summer clothing
    sclo     = 0.44453_wp
!
!-- Winter clothing
    wclo     = 1.76267_wp
!
!-- Eecision: first calculate for winter or summer clothing
    IF ( ta <= 10.0_wp )  THEN
!
!--    First guess: winter clothing insulation: cold stress
       clo = wclo
       CALL fanger ( ta, tmrt, vp, ws, pair, clo, actlev, eta, pmva )
       pmv_w = pmva

       IF ( pmva > 0.0_wp )  THEN
!
!--       Case summer clothing insulation: heat load ?
          clo = sclo
          CALL fanger ( ta, tmrt, vp, ws, pair, clo, actlev, eta, pmva )
          pmv_s = pmva
          IF ( pmva <= 0.0_wp )  THEN
!
!--          Case: comfort achievable by varying clothing insulation between winter and summer set
!--                values
             CALL iso_ridder ( ta, tmrt, vp, ws, pair, actlev, eta, sclo, pmv_s, wclo, pmv_w, eps, &
                               pmva, ncount, clo )
             IF ( ncount < 0 )  THEN
                nerr = -1
                RETURN
             ENDIF
          ELSE IF ( pmva > 0.06_wp )  THEN
             clo = 0.5_wp
             CALL fanger ( ta, tmrt, vp, ws, pair, clo, actlev, eta,  pmva )
          ENDIF
       ELSE IF ( pmva < - 0.11_wp )  THEN
          clo = 1.75_wp
          CALL fanger ( ta, tmrt, vp, ws, pair, clo, actlev, eta, pmva )
       ENDIF
    ELSE
!
!--    First guess: summer clothing insulation: heat load
       clo = sclo
       CALL fanger ( ta, tmrt, vp, ws, pair, clo, actlev, eta, pmva )
       pmv_s = pmva

       IF ( pmva < 0.0_wp )  THEN
!
!--       Case winter clothing insulation: cold stress ?
          clo = wclo
          CALL fanger ( ta, tmrt, vp, ws, pair, clo, actlev, eta, pmva )
          pmv_w = pmva

          IF ( pmva >= 0.0_wp )  THEN
!
!--          Case: comfort achievable by varying clothing insulation between winter and summer set
!--                values
             CALL iso_ridder ( ta, tmrt, vp, ws, pair, actlev, eta, sclo, pmv_s, wclo, pmv_w, eps, &
                               pmva, ncount, clo )
             IF ( ncount < 0 )  THEN
                nerr = -1
                RETURN
             ENDIF
          ELSE IF ( pmva < - 0.11_wp )  THEN
             clo = 1.75_wp
             CALL fanger ( ta, tmrt, vp, ws, pair, clo, actlev, eta, pmva )
          ENDIF
       ELSE IF ( pmva > 0.06_wp )  THEN
          clo = 0.5_wp
          CALL fanger ( ta, tmrt, vp, ws, pair, clo, actlev, eta, pmva )
       ENDIF

    ENDIF
!
!-- Determine perceived temperature by regression equation + adjustments
    pmvs = pmva
    CALL perct_regression( pmva, clo, perct_ij )
    ptc = perct_ij
    IF ( clo >= 1.75_wp  .AND.  pmva <= - 0.11_wp )  THEN
!
!--    Adjust for cold conditions according to Gagge 1986
       CALL dpmv_cold ( pmva, ta, ws, tmrt, nerr_cold, d_pmv )
       IF ( nerr_cold > 0 )  nerr = -5
       pmvs = pmva - d_pmv
       IF ( pmvs > - 0.11_wp )  THEN
          d_pmv  = 0.0_wp
          pmvs   = - 0.11_wp
       ENDIF
       CALL perct_regression( pmvs, clo, perct_ij )
    ENDIF
!     clo_fanger = clo
    clon = clo
    IF ( clo > 0.5_wp  .AND.  perct_ij <= 8.73_wp )  THEN
!
!--    Required clothing insulation (ireq) is exclusively defined for perceived temperatures (perct)
!--    less 10 (C) for a reference wind of 0.2 m/s according to 8.73 (C) for 0.1 m/s.
       clon = ireq_neutral ( perct_ij, ireq_minimal, nerr )
       clo = clon
    ENDIF
    CALL calc_sultr( ptc, dgtcm, dgtcstd, sult_lim )
    sultrieness    = .FALSE.
    d_std = -99.0_wp
    IF ( pmva > 0.06_wp  .AND.  clo <= 0.5_wp )  THEN
!
!--    Adjust for warm/humid conditions according to Gagge 1986
       CALL saturation_vapor_pressure ( ta, svp_ta )
       d_pmv  = deltapmv ( pmva, ta, vp, svp_ta, tmrt, ws, nerr )
       pmvs   = pmva + d_pmv
       CALL perct_regression( pmvs, clo, perct_ij )
       IF ( sult_lim < 99.0_wp )  THEN
          IF ( (perct_ij - ptc) > sult_lim )  sultrieness = .TRUE.
!
!--       Set factor to threshold for sultriness
          IF ( ABS( dgtcstd ) > 0.00001_wp )  THEN
             d_std = ( ( perct_ij - ptc ) - dgtcm ) / dgtcstd
          ENDIF
       ENDIF
    ENDIF

 END SUBROUTINE calculate_perct_static

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> The SUBROUTINE calculates the (saturation) water vapour pressure (hPa = hecto Pascal) for a given
!> temperature ta (degC).
!>'ta' can be the air temperature or the dew point temperature. The first will result in the current
!> vapor pressure (hPa), the latter will calulate the saturation vapor pressure (hPa).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE saturation_vapor_pressure( ta, svp_ta )

    IMPLICIT NONE

    REAL(wp), INTENT ( IN )  ::  ta     !< ambient air temperature (degC)
    REAL(wp), INTENT ( OUT ) ::  svp_ta !< water vapour pressure (hPa)

    REAL(wp)      ::  b
    REAL(wp)      ::  c


    IF ( ta < 0.0_wp )  THEN
!
!--    ta  < 0 (degC): water vapour pressure over ice
       b = 17.84362_wp
       c = 245.425_wp
    ELSE
!
!--    ta >= 0 (degC): water vapour pressure over water
       b = 17.08085_wp
       c = 234.175_wp
    ENDIF
!
!-- Saturation water vapour pressure
    svp_ta = 6.1078_wp * EXP( b * ta / ( c + ta ) )

 END SUBROUTINE saturation_vapor_pressure

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Find the clothing insulation value clo_res (clo) to make Fanger's Predicted Mean Vote (PMV) equal
!> comfort (pmva=0) for actual meteorological conditions (ta,tmrt, vp, ws, pair) and values of
!> individual's activity level.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE iso_ridder( ta, tmrt, vp, ws, pair, actlev, eta, sclo, pmv_s, wclo, pmv_w, eps, pmva,  &
                        nerr, clo_res )

    IMPLICIT NONE
!
!-- Input variables of argument list:
    REAL(wp), INTENT ( IN )  :: actlev   !< Individuals activity level per unit surface area (W/m2)
    REAL(wp), INTENT ( IN )  :: eps      !< (0.05) accuracy in clothing insulation (clo) for evaluation the root of Fanger's PMV (pmva=0)
    REAL(wp), INTENT ( IN )  :: eta      !< Individuals work efficiency (dimensionless)
    REAL(wp), INTENT ( IN )  :: pair     !< Barometric air pressure (hPa)
    REAL(wp), INTENT ( IN )  :: pmv_s    !< Fanger's PMV corresponding to sclo
    REAL(wp), INTENT ( IN )  :: pmv_w    !< Fanger's PMV corresponding to wclo
    REAL(wp), INTENT ( IN )  :: sclo     !< Lower threshold of bracketing clothing insulation (clo)
    REAL(wp), INTENT ( IN )  :: ta       !< Ambient temperature (degC)
    REAL(wp), INTENT ( IN )  :: tmrt     !< Mean radiant temperature (degC)
    REAL(wp), INTENT ( IN )  :: vp       !< Water vapour pressure (hPa)
    REAL(wp), INTENT ( IN )  :: wclo     !< Upper threshold of bracketing clothing insulation (clo)
    REAL(wp), INTENT ( IN )  :: ws       !< Wind speed (m/s) 1 m above ground
!
!-- Output variables of argument list:
    INTEGER(iwp), INTENT ( OUT ) :: nerr !< Error status / quality flag
                                         !< nerr >= 0, o.k., and nerr is the number of iterations for convergence
                                         !< nerr = -1: error = malfunction of Ridder's convergence method
                                         !< nerr = -2: error = maximum iterations (max_iteration) exceeded
                                         !< nerr = -3: error = root not bracketed between sclo and wclo

    REAL(wp), INTENT ( OUT ) :: clo_res  !< Resulting clothing insulation value (clo)
    REAL(wp), INTENT ( OUT ) :: pmva     !< 0 (set to zero, because clo is evaluated for comfort)
!
!-- Type of program variables
    INTEGER(iwp), PARAMETER  ::  max_iteration = 15  !< max number of iterations

    REAL(wp),     PARAMETER  ::  guess_0       = -1.11e30_wp  !< initial guess

    INTEGER(iwp) ::  j       !< running index

    REAL(wp) ::  clo_lower   !< lower limit of clothing insulation      (clo)
    REAL(wp) ::  clo_upper   !< upper limit of clothing insulation      (clo)
    REAL(wp) ::  sroot       !< sqrt of PMV-guess
    REAL(wp) ::  x_average   !< average of x_lower and x_upper          (clo)
    REAL(wp) ::  x_lower     !< lower guess for clothing insulation     (clo)
    REAL(wp) ::  x_new       !< preliminary result for clothing insulation (clo)
    REAL(wp) ::  x_ridder    !< current guess for clothing insulation   (clo)
    REAL(wp) ::  x_upper     !< upper guess for clothing insulation     (clo)
    REAL(wp) ::  y_average   !< average of y_lower and y_upper
    REAL(wp) ::  y_new       !< preliminary result for pred. mean vote
    REAL(wp) ::  y_lower     !< predicted mean vote for summer clothing
    REAL(wp) ::  y_upper     !< predicted mean vote for winter clothing
!
!-- Initialise
    nerr    = 0
!
!-- Set pmva = 0 (comfort): Root of PMV depending on clothing insulation
    x_ridder    = bio_initial_value
    pmva        = 0.0_wp
    clo_lower   = sclo
    y_lower     = pmv_s
    clo_upper   = wclo
    y_upper     = pmv_w
    IF ( ( y_lower > 0.0_wp .AND. y_upper < 0.0_wp )  .OR.                                         &
         ( y_lower < 0.0_wp .AND. y_upper > 0.0_wp ) )  THEN
       x_lower  = clo_lower
       x_upper  = clo_upper
       x_ridder = guess_0

       DO  j = 1, max_iteration
          x_average = 0.5_wp * ( x_lower + x_upper )
          CALL fanger ( ta, tmrt, vp, ws, pair, x_average, actlev, eta, y_average )
          sroot = SQRT( y_average**2 - y_lower * y_upper )
          IF ( ABS( sroot ) < 0.00001_wp )  THEN
             clo_res = x_average
             nerr = j
             RETURN
          ENDIF
          x_new = x_average + ( x_average - x_lower ) *                                            &
                  ( SIGN ( 1.0_wp, y_lower - y_upper ) * y_average / sroot )
          IF ( ABS( x_new - x_ridder ) <= eps )  THEN
             clo_res = x_ridder
             nerr       = j
             RETURN
          ENDIF
          x_ridder = x_new
          CALL fanger ( ta, tmrt, vp, ws, pair, x_ridder, actlev, eta, y_new )
          IF ( ABS( y_new ) < 0.00001_wp )  THEN
             clo_res = x_ridder
             nerr       = j
             RETURN
          ENDIF
          IF ( ABS( SIGN( y_average, y_new ) - y_average ) > 0.00001_wp )  THEN
             x_lower = x_average
             y_lower = y_average
             x_upper  = x_ridder
             y_upper  = y_new
          ELSE IF ( ABS( SIGN( y_lower, y_new ) - y_lower ) > 0.00001_wp )  THEN
             x_upper  = x_ridder
             y_upper  = y_new
          ELSE IF ( ABS( SIGN( y_upper, y_new ) - y_upper ) > 0.00001_wp )  THEN
             x_lower = x_ridder
             y_lower = y_new
          ELSE
!
!--          Never get here in x_ridder: singularity in y
             nerr    = -1
             clo_res = x_ridder
             RETURN
          ENDIF
          IF ( ABS( x_upper - x_lower ) <= eps )  THEN
             clo_res = x_ridder
             nerr    = j
             RETURN
          ENDIF
       ENDDO
!
!--    x_ridder exceed maximum iterations
       nerr       = -2
       clo_res = y_new
       RETURN
    ELSE IF ( ABS( y_lower ) < 0.00001_wp )  THEN
       x_ridder = clo_lower
    ELSE IF ( ABS( y_upper ) < 0.00001_wp )  THEN
       x_ridder = clo_upper
    ELSE
!
!--    x_ridder not bracketed by u_clo and o_clo
       nerr = -3
       clo_res = x_ridder
       RETURN
    ENDIF

 END SUBROUTINE iso_ridder

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Regression relations between perceived temperature (perct) and (adjusted) PMV. The regression
!> presumes the Klima-Michel settings for reference individual and reference environment.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE perct_regression( pmv, clo, perct_ij )

    IMPLICIT NONE

    REAL(wp), INTENT ( IN ) ::  clo   !< clothing insulation index (clo)
    REAL(wp), INTENT ( IN ) ::  pmv   !< Fangers predicted mean vote (dimensionless)

    REAL(wp), INTENT ( OUT ) ::  perct_ij   !< perct (degC) corresponding to given PMV / clo

    IF ( pmv <= - 0.11_wp )  THEN
       perct_ij = 5.805_wp + 12.6784_wp * pmv
    ELSE
       IF ( pmv >= + 0.01_wp )  THEN
          perct_ij = 16.826_wp + 6.163_wp * pmv
       ELSE
          perct_ij = 21.258_wp - 9.558_wp * clo
       ENDIF
    ENDIF

 END SUBROUTINE perct_regression

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> FANGER.F90
!>
!> SI-VERSION: ACTLEV W m-2, VAPOUR PRESSURE hPa
!> Calculates the current Predicted Mean Vote according to Fanger.
!> The case of free convection (ws < 0.1 m/s) is dealt with ws = 0.1 m/s
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fanger( ta, tmrt, pa, in_ws, pair, in_clo, actlev, eta, pmva )

    IMPLICIT NONE
!
!-- Input variables of argument list:
    REAL(wp), INTENT ( IN ) ::  actlev   !< Individuals activity level per unit surface area (W/m2)
    REAL(wp), INTENT ( IN ) ::  eta      !< Individuals mechanical work efficiency (dimensionless)
    REAL(wp), INTENT ( IN ) ::  in_clo   !< Clothing insulation (clo)
    REAL(wp), INTENT ( IN ) ::  in_ws    !< Wind speed (m/s) 1 m above ground
    REAL(wp), INTENT ( IN ) ::  pa       !< Water vapour pressure (hPa)
    REAL(wp), INTENT ( IN ) ::  pair     !< Barometric pressure (hPa) at site
    REAL(wp), INTENT ( IN ) ::  ta       !< Ambient air temperature (degC)
    REAL(wp), INTENT ( IN ) ::  tmrt     !< Mean radiant temperature (degC)

!
!-- Output variables of argument list:
    REAL(wp), INTENT ( OUT ) ::  pmva    !< Actual Predicted Mean Vote (PMV,
                                         !< dimensionless) according to Fanger corresponding to meteorological
                                         !< (ta,tmrt,pa,ws,pair) and individual variables (clo, actlev, eta)
!
!-- Internal variables
    INTEGER(iwp) :: i         !< running index

    REAL(wp) ::  activity     !< persons activity  (must stay == actlev, W)
    REAL(wp) ::  bc           !< preliminary result storage
    REAL(wp) ::  cc           !< preliminary result storage
    REAL(wp) ::  clo          !< clothing insulation index              (clo)
    REAL(wp) ::  dc           !< preliminary result storage
    REAL(wp) ::  ec           !< preliminary result storage
    REAL(wp) ::  f_cl         !< Increase in surface due to clothing    (factor)
    REAL(wp) ::  gc           !< preliminary result storage
    REAL(wp) ::  heat_convection  !< energy loss by autocnvection       (W)
    REAL(wp) ::  hr           !< radiational heat resistence
    REAL(wp) ::  t_clothing   !< clothing temperature                   (degree_C)
    REAL(wp) ::  t_skin_aver  !< average skin temperature               (degree_C)
    REAL(wp) ::  ws           !< wind speed                             (m/s)
    REAL(wp) ::  z1           !< Empiric factor for the adaption of the heat
                              !< ballance equation to the psycho-physical scale (Equ. 40 in FANGER)
    REAL(wp) ::  z2           !< Water vapour diffution through the skin
    REAL(wp) ::  z3           !< Sweat evaporation from the skin surface
    REAL(wp) ::  z4           !< Loss of latent heat through respiration
    REAL(wp) ::  z5           !< Loss of radiational heat
    REAL(wp) ::  z6           !< Heat loss through forced convection

!
!-- Clo must be > 0. to avoid div. by 0!
    clo = in_clo
    IF ( clo <= 0.0_wp )  clo = .001_wp
!
!-- f_cl = increase in surface due to clothing
    f_cl = 1.0_wp + 0.15_wp * clo
!
!-- Case of free convection (ws < 0.1 m/s ) not considered
    ws = in_ws
    IF ( ws < 0.1_wp )  THEN
       ws = 0.1_wp
    ENDIF
!
!-- Heat_convection = forced convection
    heat_convection = 12.1_wp * SQRT( ws * pair / 1013.25_wp )
!
!-- Activity = inner heat production per standardized surface
    activity = actlev * ( 1.0_wp - eta )
!
!-- t_skin_aver = average skin temperature
    t_skin_aver = 35.7_wp - 0.0275_wp * activity
!
!-- Calculation of constants for evaluation below
    bc = 0.155_wp * clo * 3.96_wp * 10.0_wp**( -8 ) * f_cl
    cc = f_cl * heat_convection
    ec = 0.155_wp * clo
    dc = ( 1.0_wp + ec * cc ) / bc
    gc = ( t_skin_aver + bc * ( tmrt + degc_to_k )**4 + ec * cc * ta ) / bc
!
!-- Calculation of clothing surface temperature (t_clothing) based on Newton-approximation with air
!-- temperature as initial guess.
    t_clothing = ta
    DO  i = 1, 3
       t_clothing = t_clothing - ( ( t_clothing + degc_to_k )**4 + t_clothing * dc - gc ) /        &
                    ( 4.0_wp * ( t_clothing + degc_to_k )**3 + dc )
    ENDDO
!
!-- Empiric factor for the adaption of the heat ballance equation to the psycho-physical scale (Equ.
!-- 40 in FANGER)
    z1 = ( 0.303_wp * EXP( - 0.036_wp * actlev ) + 0.0275_wp )
!
!-- Water vapour diffution through the skin
    z2 = 0.31_wp * ( 57.3_wp - 0.07_wp * activity-pa )
!
!-- Sweat evaporation from the skin surface
    z3 = 0.42_wp * ( activity - 58.0_wp )
!
!-- Loss of latent heat through respiration
    z4 = 0.0017_wp * actlev * ( 58.7_wp - pa ) + 0.0014_wp * actlev *                              &
         ( 34.0_wp - ta )
!
!-- Loss of radiational heat
    z5 = 3.96e-8_wp * f_cl * ( ( t_clothing + degc_to_k )**4 - ( tmrt + degc_to_k )**4 )
    IF ( ABS( t_clothing - tmrt ) > 0.0_wp )  THEN
       hr = z5 / f_cl / ( t_clothing - tmrt )
    ELSE
       hr = 0.0_wp
    ENDIF
!
!-- Heat loss through forced convection cc*(t_clothing-TT)
    z6 = cc * ( t_clothing - ta )
!
!-- Predicted Mean Vote
    pmva = z1 * ( activity - z2 - z3 - z4 - z5 - z6 )

 END SUBROUTINE fanger

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> For pmva > 0 and clo =0.5 the increment (deltapmv) is calculated that converts pmva into Gagge's
!> et al. (1986) PMV*.
!--------------------------------------------------------------------------------------------------!
 REAL(wp) FUNCTION deltapmv( pmva, ta, vp, svp_ta, tmrt, ws, nerr )

    IMPLICIT NONE

!
!-- Input variables of argument list:
    REAL(wp),     INTENT ( IN )  :: pmva     !< Actual Predicted Mean Vote (PMV) according to Fanger
    REAL(wp),     INTENT ( IN )  :: svp_ta   !< Saturation water vapour pressure (hPa) at ta
    REAL(wp),     INTENT ( IN )  :: ta       !< Ambient temperature (degC) at screen level
    REAL(wp),     INTENT ( IN )  :: tmrt     !< Mean radiant temperature (degC) at screen level
    REAL(wp),     INTENT ( IN )  :: vp       !< Water vapour pressure (hPa) at screen level
    REAL(wp),     INTENT ( IN )  :: ws       !< Wind speed (m/s) 1 m above ground

!
!-- Output variables of argument list:
    INTEGER(iwp), INTENT ( OUT ) :: nerr     !< Error status / quality flag
                                             !<  0 = o.k.
                                             !< -2 = pmva outside valid regression range
                                             !< -3 = rel. humidity set to 5 % or 95 %, respectively
                                             !< -4 = deltapmv set to avoid pmvs < 0

!
!-- Internal variables:
    INTEGER(iwp) :: nreg      !<

    REAL(wp) ::  apa          !< natural logarithm of pa (with hard lower border)
    REAL(wp) ::  dapa         !< difference of apa and pa_p50
    REAL(wp) ::  dpmv_1       !<
    REAL(wp) ::  dpmv_2       !<
    REAL(wp) ::  dtmrt        !< difference mean radiation to air temperature
    REAL(wp) ::  pa           !< vapor pressure (hPa) with hard bounds
    REAL(wp) ::  pa_p50       !< ratio actual water vapour pressure to that of relative humidity of
                              !< 50 %
    REAL(wp) ::  pmv          !< temp storage og predicted mean vote
    REAL(wp) ::  pmvs         !<
    REAL(wp) ::  p10          !< lower bound for pa
    REAL(wp) ::  p95          !< upper bound for pa
    REAL(wp) ::  sqvel        !< square root of local wind velocity
    REAL(wp) ::  weight       !<
    REAL(wp) ::  weight2      !<

!
!-- Regression coefficients:
    REAL(wp), DIMENSION(0:7), PARAMETER ::  bpmv = (/                                              &
     - 0.0556602_wp, - 0.1528680_wp, - 0.2336104_wp, - 0.2789387_wp,                               &
     - 0.3551048_wp, - 0.4304076_wp, - 0.4884961_wp, - 0.4897495_wp /)

    REAL(wp), DIMENSION(0:7), PARAMETER ::  bpa_p50 = (/                                           &
     - 0.1607154_wp, - 0.4177296_wp, - 0.4120541_wp, - 0.0886564_wp,                               &
       0.4285938_wp,   0.6281256_wp,   0.5067361_wp,   0.3965169_wp /)

    REAL(wp), DIMENSION(0:7), PARAMETER ::  bpa = (/                                               &
       0.0580284_wp,   0.0836264_wp,   0.1009919_wp,   0.1020777_wp,                               &
       0.0898681_wp,   0.0839116_wp,   0.0853258_wp,   0.0866589_wp /)

    REAL(wp), DIMENSION(0:7), PARAMETER ::  bapa = (/                                              &
     - 1.7838788_wp, - 2.9306231_wp, - 1.6350334_wp,    0.6211547_wp,                              &
       3.3918083_wp,   5.5521025_wp,   8.4897418_wp,   16.6265851_wp /)

    REAL(wp), DIMENSION(0:7), PARAMETER ::  bdapa = (/                                             &
       1.6752720_wp,   2.7379504_wp,   1.2940526_wp, -  1.0985759_wp,                              &
     - 3.9054732_wp, - 6.0403012_wp, - 8.9437119_wp, - 17.0671201_wp /)

    REAL(wp), DIMENSION(0:7), PARAMETER ::  bsqvel = (/                                            &
     - 0.0315598_wp, - 0.0286272_wp, - 0.0009228_wp,   0.0483344_wp,                               &
       0.0992366_wp,   0.1491379_wp,   0.1951452_wp,   0.2133949_wp /)

    REAL(wp), DIMENSION(0:7), PARAMETER ::  bta = (/                                               &
       0.0953986_wp,   0.1524760_wp,   0.0564241_wp, - 0.0893253_wp,                               &
     - 0.2398868_wp, - 0.3515237_wp, - 0.5095144_wp, - 0.9469258_wp /)

    REAL(wp), DIMENSION(0:7), PARAMETER ::  bdtmrt = (/                                            &
     - 0.0004672_wp, - 0.0000514_wp, - 0.0018037_wp, - 0.0049440_wp,                               &
     - 0.0069036_wp, - 0.0075844_wp, - 0.0079602_wp, - 0.0089439_wp /)

    REAL(wp), DIMENSION(0:7), PARAMETER ::  aconst = (/                                            &
       1.8686215_wp,   3.4260713_wp,    2.0116185_wp, -  0.7777552_wp,                             &
     - 4.6715853_wp, - 7.7314281_wp, - 11.7602578_wp, - 23.5934198_wp /)


!
!-- Test for compliance with regression range
    IF ( pmva < -1.0_wp  .OR.  pmva > 7.0_wp )  THEN
       nerr = -2
    ELSE
       nerr = 0
    ENDIF
!
!-- Initialise classic PMV
    pmv  = pmva
!
!-- Water vapour pressure of air
    p10  = 0.05_wp * svp_ta
    p95  = 1.00_wp * svp_ta
    IF ( vp >= p10  .AND.  vp <= p95 )  THEN
       pa = vp
    ELSE
       nerr = -3
       IF ( vp < p10 )  THEN
!
!--       Due to conditions of regression: r.H. >= 5 %
          pa = p10
       ELSE
!
!--       Due to conditions of regression: r.H. <= 95 %
          pa = p95
       ENDIF
    ENDIF
    IF ( pa > 0.0_wp )  THEN
!
!--    Natural logarithm of pa
       apa = LOG( pa )
    ELSE
       apa = -5.0_wp
    ENDIF
!
!-- Ratio actual water vapour pressure to that of a r.H. of 50 %
    pa_p50   = 0.5_wp * svp_ta
    IF ( pa_p50 > 0.0_wp  .AND.  pa > 0.0_wp )  THEN
       dapa   = apa - LOG( pa_p50 )
       pa_p50 = pa / pa_p50
    ELSE
       dapa   = -5.0_wp
       pa_p50 = 0.0_wp
    ENDIF
!
!-- Square root of wind velocity
    IF ( ws >= 0.0_wp )  THEN
       sqvel = SQRT( ws )
    ELSE
       sqvel = 0.0_wp
    ENDIF
!
!-- Difference mean radiation to air temperature
    dtmrt = tmrt - ta
!
!-- Select the valid regression coefficients
    nreg = INT( pmv )
    IF ( nreg < 0 )  THEN
!
!--    Value of the FUNCTION in the case pmv <= -1
       deltapmv = 0.0_wp
       RETURN
    ENDIF
    weight = MOD ( pmv, 1.0_wp )
    IF ( weight < 0.0_wp )  weight = 0.0_wp
    IF ( nreg > 5 )  THEN
       nreg  = 5
       weight   = pmv - 5.0_wp
       weight2  = pmv - 6.0_wp
       IF ( weight2 > 0 )  THEN
          weight = ( weight - weight2 ) / weight
       ENDIF
    ENDIF
!
!-- Regression valid for 0. <= pmv <= 6., bounds are checked above
    dpmv_1 =                                                                                       &
             + bpa(nreg)     * pa                                                                  &
             + bpmv(nreg)    * pmv                                                                 &
             + bapa(nreg)    * apa                                                                 &
             + bta(nreg)     * ta                                                                  &
             + bdtmrt(nreg)  * dtmrt                                                               &
             + bdapa(nreg)   * dapa                                                                &
             + bsqvel(nreg)  * sqvel                                                               &
             + bpa_p50(nreg) * pa_p50                                                              &
             + aconst(nreg)

!    dpmv_2 = 0.0_wp
!    IF ( nreg < 6 )  THEN  !< nreg is always <= 5, see above
    dpmv_2 =                                                                                       &
             + bpa(nreg+1)     * pa                                                                &
             + bpmv(nreg+1)    * pmv                                                               &
             + bapa(nreg+1)    * apa                                                               &
             + bta(nreg+1)     * ta                                                                &
             + bdtmrt(nreg+1)  * dtmrt                                                             &
             + bdapa(nreg+1)   * dapa                                                              &
             + bsqvel(nreg+1)  * sqvel                                                             &
             + bpa_p50(nreg+1) * pa_p50                                                            &
             + aconst(nreg+1)
!    ENDIF
!
!-- Calculate pmv modification
    deltapmv = ( 1.0_wp - weight ) * dpmv_1 + weight * dpmv_2
    pmvs = pmva + deltapmv
    IF ( ( pmvs ) < 0.0_wp )  THEN
!
!--    Prevent negative pmv* due to problems with clothing insulation
       nerr = -4
       IF ( pmvs > -0.11_wp )  THEN
!
!--       Threshold from perct_regression for winter clothing insulation
          deltapmv = deltapmv + 0.11_wp
       ELSE
!
!--       Set pmvs to "0" for compliance with summer clothing insulation
          deltapmv = -1.0_wp * pmva
       ENDIF
    ENDIF

 END FUNCTION deltapmv

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> The subroutine "calc_sultr" returns a threshold value to perceived temperature allowing to decide
!> whether the actual perceived temperature is linked to perecption of sultriness. The threshold
!> values depends on the Fanger's classical PMV, expressed here as perceived temperature perct.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_sultr( perct_ij, dperctm, dperctstd, sultr_res )

    IMPLICIT NONE
!
!-- Input of the argument list:
    REAL(wp), INTENT ( IN )  ::  perct_ij   !< Classical perceived temperature: Base is Fanger's PMV
!
!-- Additional output variables of argument list:
    REAL(wp), INTENT ( OUT ) ::  dperctm    !< Mean deviation perct (classical gt) to gt* (rational
                                            !< gt calculated based on Gagge's rational PMV*)
    REAL(wp), INTENT ( OUT ) ::  dperctstd  !< dperctm plus its standard deviation times a factor
                                            !< determining the significance to perceive sultriness
    REAL(wp), INTENT ( OUT ) ::  sultr_res
!
!-- Types of coefficients mean deviation: third order polynomial
    REAL(wp), PARAMETER ::  dperctka =   7.5776086_wp
    REAL(wp), PARAMETER ::  dperctkb = - 0.740603_wp
    REAL(wp), PARAMETER ::  dperctkc =   0.0213324_wp
    REAL(wp), PARAMETER ::  dperctkd = - 0.00027797237_wp
!
!-- Types of coefficients mean deviation plus standard deviation
!-- regression coefficients: third order polynomial
    REAL(wp), PARAMETER ::  dperctsa =   0.0268918_wp
    REAL(wp), PARAMETER ::  dperctsb =   0.0465957_wp
    REAL(wp), PARAMETER ::  dperctsc = - 0.00054709752_wp
    REAL(wp), PARAMETER ::  dperctsd =   0.0000063714823_wp
!
!-- Factor to mean standard deviation defining SIGNificance for
!-- sultriness
    REAL(wp), PARAMETER :: faktor = 1.0_wp
!
!-- Initialise
    sultr_res = 99.0_wp
    dperctm   = 0.0_wp
    dperctstd = 999999.0_wp

    IF ( perct_ij < 16.826_wp  .OR.  perct_ij > 56.0_wp )  THEN
!
!--    Unallowed value of classical perct!
       RETURN
    ENDIF
!
!-- Mean deviation dependent on perct
    dperctm = dperctka + dperctkb * perct_ij + dperctkc * perct_ij**2.0_wp + dperctkd *            &
              perct_ij**3.0_wp
!
!-- Mean deviation plus its standard deviation
    dperctstd = dperctsa + dperctsb * perct_ij + dperctsc * perct_ij**2.0_wp + dperctsd *          &
                perct_ij**3.0_wp
!
!-- Value of the FUNCTION
    sultr_res = dperctm + faktor * dperctstd
    IF ( ABS( sultr_res ) > 99.0_wp )  sultr_res = +99.0_wp

 END SUBROUTINE calc_sultr

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Multiple linear regression to calculate an increment delta_cold, to adjust Fanger's classical PMV
!> (pmva) by Gagge's 2 node model, applying Fanger's convective heat transfer coefficient, hcf.
!> Wind velocitiy of the reference environment is 0.10 m/s
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE dpmv_cold( pmva, ta, ws, tmrt, nerr, dpmv_cold_res )

    IMPLICIT NONE
!
!-- Type of input arguments
    REAL(wp), INTENT ( IN ) ::  pmva   !< Fanger's classical predicted mean vote
    REAL(wp), INTENT ( IN ) ::  ta     !< Air temperature 2 m above ground (degC)
    REAL(wp), INTENT ( IN ) ::  tmrt   !< Mean radiant temperature (degC)
    REAL(wp), INTENT ( IN ) ::  ws     !< Relative wind velocity 1 m above ground (m/s)
!
!-- Type of output argument
    INTEGER(iwp), INTENT ( OUT ) ::  nerr !< Error indicator: 0 = o.k., +1 = denominator for
                                          !< intersection = 0

    REAL(wp),     INTENT ( OUT ) ::  dpmv_cold_res    !< Increment to adjust pmva according to the
                                                      !< results of Gagge's 2 node model depending on the input
!
!-- Type of program variables
    INTEGER(iwp) ::  i          !< running index
    INTEGER(iwp) ::  i_bin      !< result row number

    REAL(wp) ::  delta_cold(3)
    REAL(wp) ::  dtmrt          !< delta mean radiant temperature
    REAL(wp) ::  pmv_cross(2)
    REAL(wp) ::  reg_a(3)
    REAL(wp) ::  r_denominator  !< the regression equations denominator
    REAL(wp) ::  sqrt_ws        !< sqare root of wind speed

!    REAL(wp) ::  coeff(3,5)  !< unsafe! array is (re-)writable!
!    coeff(1,1:5) =                                                             &
!       (/ +0.161_wp,   +0.130_wp, -1.125E-03_wp, +1.106E-03_wp, -4.570E-04_wp /)
!    coeff(2,1:5) =                                                             &
!       (/  0.795_wp,    0.713_wp, -8.880E-03_wp, -1.803E-03_wp, -2.816E-03_wp /)
!    coeff(3,1:5) =                                                             &
!       (/ +0.05761_wp, +0.458_wp, -1.829E-02_wp, -5.577E-03_wp, -1.970E-03_wp /)

!
!-- Coefficient of the 3 regression lines:
!      1:const      2:*pmva    3:*ta          4:*sqrt_ws     5:*dtmrt
    REAL(wp), DIMENSION(1:3,1:5), PARAMETER ::  coeff = RESHAPE( (/                                &
        0.161_wp,   0.130_wp, -1.125E-03_wp,  1.106E-03_wp, -4.570E-04_wp,                         &
        0.795_wp,   0.713_wp, -8.880E-03_wp, -1.803E-03_wp, -2.816E-03_wp,                         &
        0.05761_wp, 0.458_wp, -1.829E-02_wp, -5.577E-03_wp, -1.970E-03_wp                          &
       /), SHAPE( coeff ), order=(/ 2, 1 /)                    )
!
!-- Initialise
    nerr           = 0
    dpmv_cold_res  = 0.0_wp
    dtmrt          = tmrt - ta
    sqrt_ws        = ws
    IF ( sqrt_ws < 0.1_wp )  THEN
       sqrt_ws = 0.1_wp
    ELSE
       sqrt_ws = SQRT( sqrt_ws )
    ENDIF

    delta_cold = 0.0_wp
    pmv_cross = pmva

!
!-- Determine regression constant for given meteorological conditions
    DO  i = 1, 3
       reg_a(i) = coeff(i,1) + coeff(i,3) * ta + coeff(i,4) * sqrt_ws + coeff(i,5)*dtmrt
       delta_cold(i) = reg_a(i) + coeff(i,2) * pmva
    ENDDO
!
!-- Intersection points of regression lines in terms of Fanger's PMV
    DO  i = 1, 2
       r_denominator = coeff(i,2) - coeff(i+1,2)
       IF ( ABS( r_denominator ) > 0.00001_wp )  THEN
          pmv_cross(i) = ( reg_a(i+1) - reg_a(i) ) / r_denominator
       ELSE
          nerr = 1
          RETURN
       ENDIF
    ENDDO
!
!-- Select result row number
    i_bin = 3
    DO  i = 1, 2
       IF ( pmva > pmv_cross(i) )  THEN
          i_bin = i
          EXIT
       ENDIF
    ENDDO
!
!-- Adjust to operative temperature scaled according to classical PMV (Fanger)
    dpmv_cold_res = delta_cold(i_bin) - dpmv_cold_adj(pmva)

 END SUBROUTINE dpmv_cold

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates the summand dpmv_cold_adj adjusting to the operative temperature scaled according to
!> classical PMV (Fanger) for cold conditions. Valid for reference environment: v (1m) = 0.10 m/s,
!> dTMRT = 0 K, r.h. = 50 %
!--------------------------------------------------------------------------------------------------!
 REAL(wp) FUNCTION dpmv_cold_adj( pmva )

    IMPLICIT NONE

    INTEGER(iwp)  ::  i       !< running index
    INTEGER(iwp)  ::  thr     !< thermal range

    REAL(wp), INTENT ( IN ) ::  pmva        !< (adjusted) Predicted Mean Vote

    REAL(wp)  ::  pmv     !< pmv-part of the regression

!
!-- Provide regression coefficients for three thermal ranges:
!--                                                    slightly cold  cold           very cold
    REAL(wp), DIMENSION(1:3,0:3), PARAMETER ::  coef = RESHAPE( (/                                 &
                                                       0.0941540_wp, -0.1506620_wp, -0.0871439_wp, &
                                                       0.0783162_wp, -1.0612651_wp,  0.1695040_wp, &
                                                       0.1350144_wp, -1.0049144_wp, -0.0167627_wp, &
                                                       0.1104037_wp, -0.2005277_wp, -0.0003230_wp  &
                                                                 /), SHAPE(coef), order=(/ 1, 2 /) )
!
!-- Select thermal range
    IF ( pmva <= -2.1226_wp )  THEN     !< very cold
       thr = 3
    ELSE IF ( pmva <= -1.28_wp )  THEN  !< cold
       thr = 2
    ELSE                                !< slightly cold
       thr = 1
    ENDIF
!
!-- Initialize
    dpmv_cold_adj = coef(thr,0)
    pmv           = 1.0_wp
!
!-- Calculate pmv adjustment (dpmv_cold_adj)
    DO  i = 1, 3
       pmv           = pmv * pmva
       dpmv_cold_adj = dpmv_cold_adj + coef(thr,i) * pmv
    ENDDO

    RETURN
 END FUNCTION dpmv_cold_adj

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Based on perceived temperature (perct) as input, ireq_neutral determines the required clothing
!> insulation (clo) for thermally neutral conditions (neither body cooling nor body heating). It is
!> related to the Klima-Michel activity level (134.682 W/m2). IREQ_neutral is only defined for perct
!> < 10 (degC)
!--------------------------------------------------------------------------------------------------!
 REAL(wp) FUNCTION ireq_neutral( perct_ij, ireq_minimal, nerr )

    IMPLICIT NONE
!
!-- Type declaration of arguments
    INTEGER(iwp), INTENT ( OUT ) ::  nerr

    REAL(wp),     INTENT ( IN )  ::  perct_ij
    REAL(wp),     INTENT ( OUT ) ::  ireq_minimal
!
!-- Type declaration for internal varables
    REAL(wp)                     ::  perct02
!
!-- Initialise
    nerr = 0
!
!-- Convert perceived temperature from basis 0.1 m/s to basis 0.2 m/s
    perct02 = 1.8788_wp + 0.9296_wp * perct_ij
!
!-- IREQ neutral conditions (thermal comfort)
    ireq_neutral = 1.62_wp - 0.0564_wp * perct02
!
!-- Regression only defined for perct <= 10 (degC)
    IF ( ireq_neutral < 0.5_wp )  THEN
       IF ( ireq_neutral < 0.48_wp )  THEN
          nerr = 1
       ENDIF
       ireq_neutral = 0.5_wp
    ENDIF
!
!-- Minimal required clothing insulation: maximal acceptable body cooling
    ireq_minimal = 1.26_wp - 0.0588_wp * perct02
    IF ( nerr > 0 )  THEN
       ireq_minimal = ireq_neutral
    ENDIF

    RETURN
 END FUNCTION ireq_neutral


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> The SUBROUTINE surface area calculates the surface area of the individual according to its height
!> (m), weight (kg), and age (y)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE surface_area( height_cm, weight, age, surf )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in)  ::  age

    REAL(wp)    , INTENT(in)  ::  height_cm
    REAL(wp)    , INTENT(in)  ::  weight

    REAL(wp)    , INTENT(out) ::  surf

    REAL(wp)                  ::  height

    height = height_cm * 100.0_wp
!
!-- According to Gehan-George, for children
    IF ( age < 19 )  THEN
       IF ( age < 5 )  THEN
          surf = 0.02667_wp * height**0.42246_wp * weight**0.51456_wp
          RETURN
       ENDIF
       surf = 0.03050_wp * height**0.35129_wp * weight**0.54375_wp
       RETURN
    ENDIF
!
!-- DuBois D, DuBois EF: A formula to estimate the approximate surface area if height and weight be
!>  known. In: Arch. Int. Med.. 17, 1916, pp. 863:871.
    surf = 0.007184_wp * height**0.725_wp * weight**0.425_wp
    RETURN

 END SUBROUTINE surface_area

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> The SUBROUTINE persdat calculates
!>  - the total internal energy production = metabolic + workload,
!>  - the total internal energy production for a standardized surface (actlev)
!>  - the DuBois - area (a_surf [m2])
!> from
!>  - the persons age (years),
!>  - weight (kg),
!>  - height (m),
!>  - sex (1 = male, 2 = female),
!>  - work load (W)
!> for a sample human.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE persdat( age, weight, height, sex, work, a_surf, actlev )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(in) ::  sex


    REAL(wp), INTENT(in) ::  age
    REAL(wp), INTENT(in) ::  height
    REAL(wp), INTENT(in) ::  weight
    REAL(wp), INTENT(in) ::  work

    REAL(wp), INTENT(out) ::  actlev

    REAL(wp) ::  a_surf
    REAL(wp) ::  basic_heat_prod
    REAL(wp) ::  energy_prod
    REAL(wp) ::  factor
    REAL(wp) ::  s


    CALL surface_area( height, weight, INT( age ), a_surf )
    s = height * 100.0_wp / ( weight**( 1.0_wp / 3.0_wp ) )
    factor = 1.0_wp + .004_wp  * ( 30.0_wp - age )
    basic_heat_prod = 0.0_wp
    IF ( sex == 1 )  THEN
       basic_heat_prod = 3.45_wp * weight**( 3.0_wp / 4.0_wp ) * ( factor + 0.01_wp                &
                         * ( s - 43.4_wp ) )
    ELSE IF ( sex == 2 )  THEN
       basic_heat_prod = 3.19_wp * weight**( 3.0_wp / 4.0_wp ) * ( factor + 0.018_wp               &
                         * ( s - 42.1_wp ) )
    ENDIF

    energy_prod = work + basic_heat_prod
    actlev = energy_prod / a_surf

 END SUBROUTINE persdat


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> SUBROUTINE ipt_init
!> initializes the nonstationary perceived temperature
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE ipt_init( age, weight, height, sex, work, actlev, clo, ta, vp, ws, tmrt, pair, dt,     &
                      storage, t_clothing, ipt )

    IMPLICIT NONE
!
!-- Input parameters

    INTEGER(iwp), INTENT(in)  :: sex    !< Persons sex (1 = male, 2 = female)

    REAL(wp), INTENT(in) ::  age        !< Persons age          (years)
    REAL(wp), INTENT(in) ::  dt         !< Timestep             (s)
    REAL(wp), INTENT(in) ::  height     !< Persons height       (m)7
    REAL(wp), INTENT(in) ::  pair       !< Air pressure         (hPa)
    REAL(wp), INTENT(in) ::  ta         !< Air Temperature      (degree_C)
    REAL(wp), INTENT(in) ::  tmrt       !< Mean radiant temperature   (degree_C)
    REAL(wp), INTENT(in) ::  vp         !< Vapor pressure       (hPa)
    REAL(wp), INTENT(in) ::  weight     !< Persons weight       (kg)
    REAL(wp), INTENT(in) ::  work       !< Current workload     (W)
    REAL(wp), INTENT(in) ::  ws         !< Wind speed in approx. 1.1m (m/s)
!
!-- Output parameters
    REAL(wp), INTENT(out) ::  actlev
    REAL(wp), INTENT(out) ::  clo
    REAL(wp), INTENT(out) ::  ipt
    REAL(wp), INTENT(out) ::  storage
    REAL(wp), INTENT(out) ::  t_clothing
!
!-- Internal variables
    REAL(wp), PARAMETER :: eps = 0.0005_wp
    REAL(wp), PARAMETER :: eta = 0.0_wp

    INTEGER(iwp) ::  ncount
    INTEGER(iwp) ::  nerr_cold
    INTEGER(iwp) ::  nerr

    LOGICAL ::  sultrieness

!    REAL(wp) ::  acti
    REAL(wp) ::  a_surf
!    REAL(wp) ::  clo_fanger
    REAL(wp) ::  clon
    REAL(wp) ::  d_pmv
    REAL(wp) ::  d_std
    REAL(wp) ::  dgtcm
    REAL(wp) ::  dgtcstd
    REAL(wp) ::  ireq_minimal
    REAL(wp) ::  pmv_s
    REAL(wp) ::  pmv_w
    REAL(wp) ::  pmva
    REAL(wp) ::  pmvs
    REAL(wp) ::  ptc
    REAL(wp) ::  sclo
    REAL(wp) ::  sult_lim
    REAL(wp) ::  svp_ta
    REAL(wp) ::  wclo


    storage = 0.0_wp
    CALL persdat( age, weight, height, sex, work, a_surf, actlev )
!
!-- Initialise
    t_clothing = bio_initial_value
    ipt        = bio_initial_value
    nerr       = 0_wp
    ncount     = 0_wp
    sultrieness    = .FALSE.
!
!-- Tresholds: clothing insulation (account for model inaccuracies)
!-- Summer clothing
    sclo     = 0.44453_wp
!-- Winter clothing
    wclo     = 1.76267_wp
!
!-- Decision: firstly calculate for winter or summer clothing
    IF ( ta <= 10.0_wp )  THEN
!
!--    First guess: winter clothing insulation: cold stress
       clo = wclo
       t_clothing = bio_initial_value  ! force initial run
       CALL fanger_s_acti ( ta, tmrt, vp, ws, pair, clo, actlev, work, t_clothing, storage, dt,    &
                            pmva )
       pmv_w = pmva

       IF ( pmva > 0.0_wp )  THEN
!
!--       Case summer clothing insulation: heat load ?
          clo = sclo
          t_clothing = bio_initial_value  ! force initial run
          CALL fanger_s_acti ( ta, tmrt, vp, ws, pair, clo, actlev, work, t_clothing, storage, dt, &
                               pmva )
          pmv_s = pmva
          IF ( pmva <= 0.0_wp )  THEN
!
!--          Case: comfort achievable by varying clothing insulation between winter and summer set
!--                values
             CALL iso_ridder ( ta, tmrt, vp, ws, pair, actlev, eta , sclo, pmv_s, wclo, pmv_w, eps,&
                               pmva, ncount, clo )
             IF ( ncount < 0 )  THEN
                nerr = -1
                RETURN
             ENDIF
          ELSE IF ( pmva > 0.06_wp )  THEN
             clo = 0.5_wp
             t_clothing = bio_initial_value
             CALL fanger_s_acti ( ta, tmrt, vp, ws, pair, clo, actlev, work, t_clothing, storage,  &
                                  dt, pmva )
          ENDIF
       ELSE IF ( pmva < - 0.11_wp )  THEN
          clo = 1.75_wp
          t_clothing = bio_initial_value
          CALL fanger_s_acti ( ta, tmrt, vp, ws, pair, clo, actlev, work, t_clothing, storage, dt, &
                               pmva )
       ENDIF

    ELSE
!
!--    First guess: summer clothing insulation: heat load
       clo = sclo
       t_clothing = bio_initial_value
       CALL fanger_s_acti ( ta, tmrt, vp, ws, pair, clo, actlev, work, t_clothing, storage, dt,    &
                            pmva )
       pmv_s = pmva

       IF ( pmva < 0.0_wp )  THEN
!
!--       Case winter clothing insulation: cold stress ?
          clo = wclo
          t_clothing = bio_initial_value
          CALL fanger_s_acti ( ta, tmrt, vp, ws, pair, clo, actlev, work, t_clothing, storage, dt, &
                               pmva )
          pmv_w = pmva

          IF ( pmva >= 0.0_wp )  THEN
!
!--          Case: comfort achievable by varying clothing insulation between winter and summer set
!--                values
             CALL iso_ridder ( ta, tmrt, vp, ws, pair, actlev, eta, sclo, pmv_s, wclo, pmv_w, eps, &
                               pmva, ncount, clo )
             IF ( ncount < 0_wp )  THEN
                nerr = -1
                RETURN
             ENDIF
          ELSE IF ( pmva < - 0.11_wp )  THEN
             clo = 1.75_wp
             t_clothing = bio_initial_value
             CALL fanger_s_acti ( ta, tmrt, vp, ws, pair, clo, actlev, work, t_clothing, storage,  &
                                  dt, pmva )
          ENDIF
       ELSE IF ( pmva > 0.06_wp )  THEN
          clo = 0.5_wp
          t_clothing = bio_initial_value
          CALL fanger_s_acti ( ta, tmrt, vp, ws, pair, clo, actlev, work, t_clothing, storage, dt, &
                               pmva )
       ENDIF

    ENDIF
!
!-- Determine perceived temperature by regression equation + adjustments
    pmvs = pmva
    CALL perct_regression( pmva, clo, ipt )
    ptc = ipt
    IF ( clo >= 1.75_wp  .AND.  pmva <= - 0.11_wp )  THEN
!
!--    Adjust for cold conditions according to Gagge 1986
       CALL dpmv_cold ( pmva, ta, ws, tmrt, nerr_cold, d_pmv )
       IF ( nerr_cold > 0 )  nerr = -5
       pmvs = pmva - d_pmv
       IF ( pmvs > - 0.11_wp )  THEN
          d_pmv  = 0.0_wp
          pmvs   = - 0.11_wp
       ENDIF
       CALL perct_regression( pmvs, clo, ipt )
    ENDIF
!     clo_fanger = clo
    clon = clo
    IF ( clo > 0.5_wp  .AND.  ipt <= 8.73_wp )  THEN
!
!--    Required clothing insulation (ireq) is exclusively defined for perceived temperatures (ipt)
!--    less 10 (C) for a reference wind of 0.2 m/s according to 8.73 (C) for 0.1 m/s
       clon = ireq_neutral ( ipt, ireq_minimal, nerr )
       clo = clon
    ENDIF
    CALL calc_sultr( ptc, dgtcm, dgtcstd, sult_lim )
    sultrieness    = .FALSE.
    d_std      = - 99.0_wp
    IF ( pmva > 0.06_wp  .AND.  clo <= 0.5_wp )  THEN
!
!--    Adjust for warm/humid conditions according to Gagge 1986
       CALL saturation_vapor_pressure ( ta, svp_ta )
       d_pmv  = deltapmv ( pmva, ta, vp, svp_ta, tmrt, ws, nerr )
       pmvs   = pmva + d_pmv
       CALL perct_regression( pmvs, clo, ipt )
       IF ( sult_lim < 99.0_wp )  THEN
          IF ( (ipt - ptc) > sult_lim )  sultrieness = .TRUE.
       ENDIF
    ENDIF


 END SUBROUTINE ipt_init

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> SUBROUTINE ipt_cycle
!> Calculates one timestep for the nonstationary version of perceived temperature (iPT, degree_C) for
!>  - standard measured/predicted meteorological values and TMRT as input;
!>  - regressions for determination of PT;
!>  - adjustment to Gagge's PMV* (2-node-model, 1986) as base of PT under warm/humid conditions
!>    (Icl= 0.50 clo) and under cold conditions (Icl= 1.75 clo)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ipt_cycle( ta, vp, ws, tmrt, pair, dt, storage, t_clothing, clo, actlev, work, ipt )

    IMPLICIT NONE
!
!-- Type of input of the argument list
    REAL(wp), INTENT ( IN )  ::  actlev  !< Internal heat production    (W)
    REAL(wp), INTENT ( IN )  ::  clo     !< Clothing index              (no dim)
    REAL(wp), INTENT ( IN )  ::  dt      !< Timestep                    (s)
    REAL(wp), INTENT ( IN )  ::  pair    !< Air pressure                (hPa)
    REAL(wp), INTENT ( IN )  ::  ta      !< Air temperature             (degree_C)
    REAL(wp), INTENT ( IN )  ::  tmrt    !< Mean radiant temperature    (degree_C)
    REAL(wp), INTENT ( IN )  ::  vp      !< Vapor pressure              (hPa)
    REAL(wp), INTENT ( IN )  ::  work    !< Mechanical work load        (W)
    REAL(wp), INTENT ( IN )  ::  ws      !< Wind speed                  (m/s)
!
!-- In and output parameters
    REAL(wp), INTENT (INOUT) ::  storage     !< Heat storage            (W)
    REAL(wp), INTENT (INOUT) ::  t_clothing  !< Clothig temperature     (degree_C)
!
!-- Type of output of the argument list
    REAL(wp), INTENT ( OUT ) ::  ipt  !< Instationary perceived temperature (degree_C)
!
!-- Type of internal variables
    INTEGER(iwp) ::  nerr
    INTEGER(iwp) ::  nerr_cold

    LOGICAL ::  sultrieness

    REAL(wp) ::  d_pmv
    REAL(wp) ::  d_std
    REAL(wp) ::  dgtcm
    REAL(wp) ::  dgtcstd
    REAL(wp) ::  pmva
    REAL(wp) ::  pmvs
    REAL(wp) ::  ptc
    REAL(wp) ::  sult_lim
    REAL(wp) ::  svp_ta
!
!-- Initialise
    ipt = bio_initial_value

    nerr     = 0
    sultrieness  = .FALSE.
!
!-- Determine pmv_adjusted for current conditions
    CALL fanger_s_acti ( ta, tmrt, vp, ws, pair, clo, actlev, work, t_clothing, storage, dt, pmva )
!
!-- Determine perceived temperature by regression equation + adjustments
    CALL perct_regression( pmva, clo, ipt )
!
!-- Consider cold conditions
    IF ( clo >= 1.75_wp  .AND.  pmva <= -0.11_wp )  THEN
!
!--    Adjust for cold conditions according to Gagge 1986
       CALL dpmv_cold ( pmva, ta, ws, tmrt, nerr_cold, d_pmv )
       IF ( nerr_cold > 0 )  nerr = -5
       pmvs = pmva - d_pmv
       IF ( pmvs > - 0.11_wp )  THEN
          d_pmv  = 0.0_wp
          pmvs   = - 0.11_wp
       ENDIF
       CALL perct_regression( pmvs, clo, ipt )
    ENDIF
!
!-- Consider sultriness if appropriate
    ptc = ipt
    CALL calc_sultr( ptc, dgtcm, dgtcstd, sult_lim )
    sultrieness = .FALSE.
    d_std       = - 99.0_wp
    IF ( pmva > 0.06_wp  .AND.  clo <= 0.5_wp )  THEN
!
!--    Adjust for warm/humid conditions according to Gagge 1986
       CALL saturation_vapor_pressure ( ta, svp_ta )
       d_pmv = deltapmv ( pmva, ta, vp, svp_ta, tmrt, ws, nerr )
       pmvs  = pmva + d_pmv
       CALL perct_regression( pmvs, clo, ipt )
       IF ( sult_lim < 99.0_wp )  THEN
          IF ( (ipt - ptc) > sult_lim )  sultrieness = .TRUE.
       ENDIF
    ENDIF

 END SUBROUTINE ipt_cycle

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> SUBROUTINE fanger_s calculates the actual Predicted Mean Vote (dimensionless) according to Fanger
!> corresponding to meteorological (ta,tmrt,pa,ws,pair) and individual variables (clo, actlev, eta)
!> considering a storage and clothing temperature for a given timestep.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE fanger_s_acti( ta, tmrt, pa, in_ws, pair, in_clo, actlev, activity, t_cloth, s, dt,    &
                           pmva )

    IMPLICIT NONE
!
!--  Input argument types
    REAL(wp), INTENT ( IN )  ::  activity !< Work load                (W/m²)
    REAL(wp), INTENT ( IN )  ::  actlev   !< Metabolic + work energy  (W/m²)
    REAL(wp), INTENT ( IN )  ::  dt       !< Timestep                 (s)
    REAL(wp), INTENT ( IN )  ::  in_clo   !< Clothing index (clo)     (no dim)
    REAL(wp), INTENT ( IN )  ::  in_ws    !< Wind speed               (m/s)
    REAL(wp), INTENT ( IN )  ::  pa       !< Vapour pressure          (hPa)
    REAL(wp), INTENT ( IN )  ::  pair     !< Air pressure             (hPa)
    REAL(wp), INTENT ( IN )  ::  ta       !< Air temperature          (degree_C)
    REAL(wp), INTENT ( IN )  ::  tmrt     !< Mean radiant temperature (degree_C)
!
!-- Output argument types
    REAL(wp), INTENT ( OUT ) ::  pmva  !< actual Predicted Mean Vote (no dim)

    REAL(wp), INTENT (INOUT) ::  s  !< storage var. of energy balance (W/m2)
    REAL(wp), INTENT (INOUT) ::  t_cloth  !< clothing temperature (degree_C)
!
!-- Internal variables
    REAL(wp), PARAMETER  ::  time_equil = 7200.0_wp

    INTEGER(iwp) :: i         !< running index
    INTEGER(iwp) ::  niter    !< Running index

    REAL(wp) ::  adjustrate        !< Max storage adjustment rate
    REAL(wp) ::  adjustrate_cloth  !< max clothing temp. adjustment rate
    REAL(wp) ::  bc                !< preliminary result storage
    REAL(wp) ::  cc                !< preliminary result storage
    REAL(wp) ::  clo               !< clothing insulation index              (clo)
    REAL(wp) ::  d_s               !< Storage delta                          (W)
    REAL(wp) ::  dc                !< preliminary result storage
    REAL(wp) ::  en                !< Energy ballance                        (W)
    REAL(wp) ::  ec                !< preliminary result storage
    REAL(wp) ::  f_cl              !< Increase in surface due to clothing    (factor)
    REAL(wp) ::  gc                !< preliminary result storage
    REAL(wp) ::  heat_convection   !< energy loss by autocnvection       (W)
!    REAL(wp) ::  hr                !< radiational heat resistence
    REAL(wp) ::  t_clothing        !< clothing temperature                   (degree_C)
    REAL(wp) ::  t_skin_aver       !< average skin temperature               (degree_C)
    REAL(wp) ::  ws                !< wind speed                             (m/s)
    REAL(wp) ::  z1                !< Empiric factor for the adaption of the heat
                                   !< ballance equation to the psycho-physical scale
                                   !< (Equ. 40 in FANGER)
    REAL(wp) ::  z2                !< Water vapour diffution through the skin
    REAL(wp) ::  z3                !< Sweat evaporation from the skin surface
    REAL(wp) ::  z4                !< Loss of latent heat through respiration
    REAL(wp) ::  z5                !< Loss of radiational heat
    REAL(wp) ::  z6                !< Heat loss through forced convection




!
!-- Clo must be > 0. to avoid div. by 0!
    clo = in_clo
    IF ( clo < 001.0_wp )  clo = 0.001_wp
!
!-- Increase in surface due to clothing
    f_cl = 1.0_wp + 0.15_wp * clo
!
!-- Case of free convection (ws < 0.1 m/s ) not considered
    ws = in_ws
    IF ( ws < 0.1_wp )  THEN
       ws = 0.1_wp
    ENDIF
!
!-- Heat_convection = forced convection
    heat_convection = 12.1_wp * SQRT( ws * pair / 1013.25_wp )
!
!-- Average skin temperature
    t_skin_aver = 35.7_wp - 0.0275_wp * activity
!
!-- Calculation of constants for evaluation below
    bc = 0.155_wp * clo * 3.96_wp * 10.0_wp**( -8.0_wp ) * f_cl
    cc = f_cl * heat_convection
    ec = 0.155_wp * clo
    dc = ( 1.0_wp + ec * cc ) / bc
    gc = ( t_skin_aver + bc * ( tmrt + 273.2_wp )**4.0_wp + ec * cc * ta ) / bc
!
!-- Calculation of clothing surface temperature (t_clothing) based on Newton-approximation with air
!-- temperature as initial guess
    niter = INT( dt * 10.0_wp, KIND=iwp )
    IF ( niter < 1 )  niter = 1
    adjustrate = 1.0_wp - EXP( -1.0_wp * ( 10.0_wp / time_equil ) * dt )
    IF ( adjustrate >= 1.0_wp )  adjustrate = 1.0_wp
    adjustrate_cloth = adjustrate * 30.0_wp
    t_clothing = t_cloth
!
!-- Set initial values for niter, adjustrates and t_clothing if this is the first call
    IF ( t_cloth <= -998.0_wp )  THEN  ! If initial run
       niter = 3
       adjustrate = 1.0_wp
       adjustrate_cloth = 1.0_wp
       t_clothing = ta
    ENDIF
!
!-- Update clothing temperature
    DO  i = 1, niter
       t_clothing = t_clothing - adjustrate_cloth * ( ( t_clothing + 273.2_wp )**4.0_wp  +         &
                    t_clothing * dc - gc ) / ( 4.0_wp * ( t_clothing + 273.2_wp )**3.0_wp + dc )
    ENDDO
!
!-- Empiric factor for the adaption of the heat ballance equation to the psycho-physical scale
!-- (Equ. 40 in FANGER)
    z1 = ( 0.303_wp * EXP( - 0.036_wp * actlev ) + 0.0275_wp )
!
!-- Water vapour diffution through the skin
    z2 = 0.31_wp * ( 57.3_wp - 0.07_wp * activity-pa )
!
!-- Sweat evaporation from the skin surface
    z3 = 0.42_wp * ( activity - 58.0_wp )
!
!-- Loss of latent heat through respiration
    z4 = 0.0017_wp * actlev * ( 58.7_wp - pa ) + 0.0014_wp * actlev * ( 34.0_wp - ta )
!
!-- Loss of radiational heat
    z5 = 3.96e-8_wp * f_cl * ( ( t_clothing + 273.2_wp )**4 - ( tmrt + 273.2_wp )**4 )
!
!-- Heat loss through forced convection
    z6 = cc * ( t_clothing - ta )
!
!-- Write together as energy ballance
    en = activity - z2 - z3 - z4 - z5 - z6
!
!-- Manage storage
    d_s = adjustrate * en + ( 1.0_wp - adjustrate ) * s
!
!-- Predicted Mean Vote
    pmva = z1 * d_s
!
!-- Update storage
    s = d_s
    t_cloth = t_clothing

 END SUBROUTINE fanger_s_acti



!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Physiologically Equivalent Temperature (PET),
!> stationary (calculated based on MEMI),
!> Subroutine based on PETBER vers. 1.5.1996 by P. Hoeppe
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE calculate_pet_static( ta, vpa, v, tmrt, pair, pet_ij )

    IMPLICIT NONE
!
!-- Input arguments:
    REAL(wp), INTENT( IN ) ::  pair  !< Air pressure                (hPa)
    REAL(wp), INTENT( IN ) ::  ta    !< Air temperature             (degree_C)
    REAL(wp), INTENT( IN ) ::  tmrt  !< Mean radiant temperature    (degree_C)
    REAL(wp), INTENT( IN ) ::  v     !< Wind speed                  (m/s)
    REAL(wp), INTENT( IN ) ::  vpa   !< Vapor pressure              (hPa)
!
!-- Output arguments:
    REAL(wp), INTENT ( OUT ) ::  pet_ij  !< PET                     (degree_C)
!
!-- Internal variables:
    REAL(wp) ::  acl        !< clothing area                        (m²)
    REAL(wp) ::  adu        !< Du Bois area                         (m²)
    REAL(wp) ::  aeff       !< effective area                       (m²)
    REAL(wp) ::  ere        !< energy ballance                      (W)
    REAL(wp) ::  erel       !< latent energy ballance               (W)
    REAL(wp) ::  esw        !< Energy-loss through sweat evap.      (W)
    REAL(wp) ::  facl       !< Surface area extension through clothing (factor)
    REAL(wp) ::  feff       !< Surface modification by posture      (factor)
    REAL(wp) ::  rdcl       !< Diffusion resistence of clothing     (factor)
    REAL(wp) ::  rdsk       !< Diffusion resistence of skin         (factor)
    REAL(wp) ::  rtv
    REAL(wp) ::  vpts       !< Sat. vapor pressure over skin        (hPa)
    REAL(wp) ::  tcl        !< Clothing temperature                 (degree_C)
    REAL(wp) ::  tsk        !< Skin temperature                     (degree_C)
    REAL(wp) ::  wetsk      !< Fraction of wet skin                 (factor)
!
!-- Variables:
    REAL(wp) :: int_heat    !< Internal heat        (W)
!
!-- MEMI configuration
    REAL(wp) :: age         !< Persons age          (a)
    REAL(wp) :: clo         !< Clothing insulation index (clo)
    REAL(wp) :: eta         !< Work efficiency      (dimensionless)
    REAL(wp) :: fcl         !< Surface area modification by clothing (factor)
    REAL(wp) :: ht          !< Persons height       (m)
    REAL(wp) :: mbody       !< Persons body mass    (kg)
    REAL(wp) :: work        !< Work load            (W)
!    INTEGER(iwp) :: pos     !< Posture: 1 = standing, 2 = sitting
!    INTEGER(iwp) :: sex     !< Sex: 1 = male, 2 = female
!
!-- Configuration, keep standard parameters!
    age   = 35.0_wp
    mbody = 75.0_wp
    ht    =  1.75_wp
    work  = 80.0_wp
    eta   =  0.0_wp
    clo   =  0.9_wp
    fcl   =  1.15_wp
!
!-- Call subfunctions
    CALL in_body( age, eta, ere, erel, ht, int_heat, mbody, pair, rtv, ta, vpa, work )

    CALL heat_exch( acl, adu, aeff, clo, ere, erel, esw, facl, fcl, feff, ht, int_heat, mbody,     &
                    pair, rdcl, rdsk, ta, tcl, tmrt, tsk, v, vpa, vpts, wetsk )

    CALL pet_iteration( acl, adu, aeff, esw, facl, feff, int_heat, pair, rdcl, rdsk, rtv, ta, tcl, &
                        tsk, pet_ij, vpts, wetsk )


 END SUBROUTINE calculate_pet_static


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate internal energy ballance
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE in_body( age, eta, ere, erel, ht, int_heat, mbody, pair, rtv, ta, vpa, work )
!
!-- Input arguments:
    REAL(wp), INTENT( IN )  ::  age       !< Persons age              (a)
    REAL(wp), INTENT( IN )  ::  eta       !< Work efficiency     (dimensionless)
    REAL(wp), INTENT( IN )  ::  ht        !< Persons height           (m)
    REAL(wp), INTENT( IN )  ::  mbody     !< Persons body mass        (kg)
    REAL(wp), INTENT( IN )  ::  pair      !< air pressure             (hPa)
    REAL(wp), INTENT( IN )  ::  ta        !< air temperature          (degree_C)
    REAL(wp), INTENT( IN )  ::  vpa       !< vapor pressure           (hPa)
    REAL(wp), INTENT( IN )  ::  work      !< Work load                (W)
!
!-- Output arguments:
    REAL(wp), INTENT( OUT ) ::  ere       !< energy ballance          (W)
    REAL(wp), INTENT( OUT ) ::  erel      !< latent energy ballance   (W)
    REAL(wp), INTENT( OUT ) ::  int_heat  !< internal heat production (W)
    REAL(wp), INTENT( OUT ) ::  rtv       !< respiratory volume
!
!-- Internal variables:
    REAL(wp) ::  eres                     !< Sensible respiratory heat flux (W)
    REAL(wp) ::  met
    REAL(wp) ::  tex
    REAL(wp) ::  vpex

!
!-- Metabolic heat production
    met = 3.45_wp * mbody**( 3.0_wp / 4.0_wp ) * (1.0_wp + 0.004_wp *                              &
          ( 30.0_wp - age) + 0.010_wp * ( ( ht * 100.0_wp /                                        &
          ( mbody**( 1.0_wp / 3.0_wp ) ) ) - 43.4_wp ) )
    met = work + met
    int_heat = met * (1.0_wp - eta)
!
!-- Sensible respiration energy
    tex  = 0.47_wp * ta + 21.0_wp
    rtv  = 1.44_wp * 10.0_wp**(-6.0_wp) * met
    eres = c_p * (ta - tex) * rtv
!
!-- Latent respiration energy
    vpex = 6.11_wp * 10.0_wp**( 7.45_wp * tex / ( 235.0_wp + tex ) )
    erel = 0.623_wp * l_v / pair * ( vpa - vpex ) * rtv
!
!-- Sum of the results
    ere = eres + erel

 END SUBROUTINE in_body


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate heat gain or loss
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE heat_exch( acl, adu, aeff, clo, ere, erel, esw, facl, fcl, feff, ht, int_heat, mbody,  &
                       pair, rdcl, rdsk, ta, tcl, tmrt, tsk, v, vpa, vpts, wetsk )

!
!-- Input arguments:
    REAL(wp), INTENT( IN )  ::  clo    !< clothing insulation      (clo)
    REAL(wp), INTENT( IN )  ::  fcl    !< factor for surface area increase by clothing
    REAL(wp), INTENT( IN )  ::  ere    !< Energy ballance          (W)
    REAL(wp), INTENT( IN )  ::  erel   !< Latent energy ballance   (W)
    REAL(wp), INTENT( IN )  ::  ht     !< height                   (m)
    REAL(wp), INTENT( IN )  ::  int_heat  !< internal heat production (W)
    REAL(wp), INTENT( IN )  ::  mbody  !< body mass                (kg)
    REAL(wp), INTENT( IN )  ::  pair   !< Air pressure             (hPa)
    REAL(wp), INTENT( IN )  ::  ta     !< Air temperature          (degree_C)
    REAL(wp), INTENT( IN )  ::  tmrt   !< Mean radiant temperature (degree_C)
    REAL(wp), INTENT( IN )  ::  v      !< Wind speed               (m/s)
    REAL(wp), INTENT( IN )  ::  vpa    !< Vapor pressure           (hPa)
!
!-- Output arguments:
    REAL(wp), INTENT( OUT ) ::  acl    !< Clothing surface area        (m²)
    REAL(wp), INTENT( OUT ) ::  adu    !< Du-Bois area                 (m²)
    REAL(wp), INTENT( OUT ) ::  aeff   !< Effective surface area       (m²)
    REAL(wp), INTENT( OUT ) ::  esw    !< Energy-loss through sweat evap. (W)
    REAL(wp), INTENT( OUT ) ::  facl   !< Surface area extension through clothing (factor)
    REAL(wp), INTENT( OUT ) ::  feff   !< Surface modification by posture (factor)
    REAL(wp), INTENT( OUT ) ::  rdcl   !< Diffusion resistence of clothing (factor)
    REAL(wp), INTENT( OUT ) ::  rdsk   !< Diffusion resistence of skin (factor)
    REAL(wp), INTENT( OUT ) ::  tcl    !< Clothing temperature         (degree_C)
    REAL(wp), INTENT( OUT ) ::  tsk    !< Skin temperature             (degree_C)
    REAL(wp), INTENT( OUT ) ::  vpts   !< Sat. vapor pressure over skin (hPa)
    REAL(wp), INTENT( OUT ) ::  wetsk  !< Fraction of wet skin (dimensionless)
!
!-- Cconstants:
!     REAL(wp), PARAMETER :: cair = 1010.0_wp      !< replaced by c_p
    REAL(wp), PARAMETER :: cb   = 3640.0_wp        !<
    REAL(wp), PARAMETER :: emcl =    0.95_wp      !< Longwave emission coef. of cloth
    REAL(wp), PARAMETER :: emsk =    0.99_wp      !< Longwave emission coef. of skin
!    REAL(wp), PARAMETER :: evap = 2.42_wp * 10.0_wp **6.0_wp  !< replaced by l_v
    REAL(wp), PARAMETER :: food =    0.0_wp        !< Heat gain by food        (W)
    REAL(wp), PARAMETER :: po   = 1013.25_wp      !< Air pressure at sea level (hPa)
    REAL(wp), PARAMETER :: rob  =    1.06_wp      !<
!
!-- Internal variables
    INTEGER(iwp) ::  count1     !< running index
    INTEGER(iwp) ::  count3     !< running index
    INTEGER(iwp) ::  j          !< running index
    INTEGER(iwp) ::  i          !< running index

    LOGICAL ::  skipincreasecount   !< iteration control flag

    REAL(wp) ::  cbare          !< Convection through bare skin
    REAL(wp) ::  cclo           !< Convection through clothing
    REAL(wp) ::  csum           !< Convection in total
    REAL(wp) ::  di             !< difference between r1 and r2
    REAL(wp) ::  ed             !< energy transfer by diffusion     (W)
    REAL(wp) ::  enbal          !< energy ballance                  (W)
    REAL(wp) ::  enbal2         !< energy ballance (storage, last cycle)
    REAL(wp) ::  eswdif         !< difference between sweat production and evaporation potential
    REAL(wp) ::  eswphy         !< sweat created by physiology
    REAL(wp) ::  eswpot         !< potential sweat evaporation
    REAL(wp) ::  fec            !<
    REAL(wp) ::  hc             !<
    REAL(wp) ::  he             !<
    REAL(wp) ::  htcl           !<
    REAL(wp) ::  r1             !<
    REAL(wp) ::  r2             !<
    REAL(wp) ::  rbare          !< Radiational loss of bare skin    (W/m²)
    REAL(wp) ::  rcl            !<
    REAL(wp) ::  rclo           !< Radiational loss of clothing     (W/m²)
    REAL(wp) ::  rclo2          !< Longwave radiation gain or loss  (W/m²)
    REAL(wp) ::  rsum           !< Radiational loss or gain         (W/m²)
    REAL(wp) ::  sw             !<
!    REAL(wp) ::  swf            !< female factor, currently unused
    REAL(wp) ::  swm            !<
    REAL(wp) ::  tbody          !<
    REAL(wp) ::  vb             !<
    REAL(wp) ::  vb1            !<
    REAL(wp) ::  vb2            !<
    REAL(wp) ::  wd             !<
    REAL(wp) ::  wr             !<
    REAL(wp) ::  ws             !<
    REAL(wp) ::  wsum           !<
    REAL(wp) ::  xx             !< modification step                (K)
    REAL(wp) ::  y              !< fraction of bare skin

    REAL(wp) ::  c(0:10)        !< Core temperature array           (degree_C)
    REAL(wp) ::  tcore(1:7)     !<

!
!-- Initialize
    wetsk = 0.0_wp  !< skin is dry everywhere on init (no non-evaporated sweat)
!
!-- Set Du Bois Area for the sample person
    adu = 0.203_wp * mbody**0.425_wp * ht**0.725_wp
!
!-- Initialize convective heat considering local air preassure
    hc = 2.67_wp + ( 6.5_wp * v**0.67_wp )
    hc = hc * ( pair / po )**0.55_wp
!
!-- Set surface modification by posture (the person will always stand)
    feff = 0.725_wp                     !< Posture: 0.725 for stading
!
!-- Set surface modification by clothing
    facl = ( - 2.36_wp + 173.51_wp * clo - 100.76_wp * clo * clo + 19.28_wp * ( clo**3.0_wp ) )    &
           / 100.0_wp
    IF ( facl > 1.0_wp )  facl = 1.0_wp
!
!-- Initialize heat resistences
    rcl = ( clo / 6.45_wp ) / facl
    IF ( clo >= 2.0_wp )  y  = 1.0_wp
    IF ( ( clo > 0.6_wp )   .AND.  ( clo < 2.0_wp ) )   y = ( ht - 0.2_wp ) / ht
    IF ( ( clo <= 0.6_wp )  .AND.  ( clo > 0.3_wp ) )  y = 0.5_wp
    IF ( ( clo <= 0.3_wp )  .AND.  ( clo > 0.0_wp ) )   y = 0.1_wp
    r2   = adu * ( fcl - 1.0_wp + facl ) / ( 2.0_wp * 3.14_wp * ht * y )
    r1   = facl * adu / ( 2.0_wp * 3.14_wp * ht * y )
    di   = r2 - r1

!
!-- Estimate skin temperatur
    DO  j = 1, 7

       tsk    = 34.0_wp
       count1 = 0
       tcl    = ( ta + tmrt + tsk ) / 3.0_wp
       count3 = 1
       enbal2 = 0.0_wp

       DO  i = 1, 100  ! allow for 100 iterations max
          acl   = adu * facl + adu * ( fcl - 1.0_wp )
          rclo2 = emcl * sigma_sb * ( ( tcl + degc_to_k )**4.0_wp -                                &
                  ( tmrt + degc_to_k )**4.0_wp ) * feff
          htcl  = 6.28_wp * ht * y * di / ( rcl * LOG( r2 / r1 ) * acl )
          tsk   = 1.0_wp / htcl * ( hc * ( tcl - ta ) + rclo2 ) + tcl
!
!--       Radiation saldo
          aeff  = adu * feff
          rbare = aeff * ( 1.0_wp - facl ) * emsk * sigma_sb *                                     &
                  ( ( tmrt + degc_to_k )**4.0_wp - ( tsk + degc_to_k )**4.0_wp )
          rclo  = feff * acl * emcl * sigma_sb *                                                   &
                  ( ( tmrt + degc_to_k )**4.0_wp - ( tcl + degc_to_k )**4.0_wp )
          rsum  = rbare + rclo
!
!--       Convection
          cbare = hc * ( ta - tsk ) * adu * ( 1.0_wp - facl )
          cclo  = hc * ( ta - tcl ) * acl
          csum  = cbare + cclo
!
!--       Core temperature
          c(0)  = int_heat + ere
          c(1)  = adu * rob * cb
          c(2)  = 18.0_wp - 0.5_wp * tsk
          c(3)  = 5.28_wp * adu * c(2)
          c(4)  = 0.0208_wp * c(1)
          c(5)  = 0.76075_wp * c(1)
          c(6)  = c(3) - c(5) - tsk * c(4)
          c(7)  = - c(0) * c(2) - tsk * c(3) + tsk * c(5)
          c(8)  = c(6) * c(6) - 4.0_wp * c(4) * c(7)
          c(9)  = 5.28_wp * adu - c(5) - c(4) * tsk
          c(10) = c(9) * c(9) - 4.0_wp * c(4) * ( c(5) * tsk - c(0) - 5.28_wp * adu * tsk )

          IF ( ABS( tsk - 36.0_wp ) < 0.00001_wp )  tsk = 36.01_wp
          tcore(7) = c(0) / ( 5.28_wp * adu + c(1) * 6.3_wp / 3600.0_wp ) + tsk
          tcore(3) = c(0) / ( 5.28_wp * adu + ( c(1) * 6.3_wp / 3600.0_wp ) /   &
                     ( 1.0_wp + 0.5_wp * ( 34.0_wp - tsk ) ) ) + tsk
          IF ( c(10) >= 0.0_wp )  THEN
             tcore(6) = ( - c(9) - c(10)**0.5_wp ) / ( 2.0_wp * c(4) )
             tcore(1) = ( - c(9) + c(10)**0.5_wp ) / ( 2.0_wp * c(4) )
          ENDIF

          IF ( c(8) >= 0.0_wp )  THEN
             tcore(2) = ( - c(6) + ABS( c(8) )**0.5_wp ) / ( 2.0_wp * c(4) )
             tcore(5) = ( - c(6) - ABS( c(8) )**0.5_wp ) / ( 2.0_wp * c(4) )
             tcore(4) = c(0) / ( 5.28_wp * adu + c(1) * 1.0_wp / 40.0_wp ) + tsk
          ENDIF
!
!--       Transpiration
          tbody = 0.1_wp * tsk + 0.9_wp * tcore(j)
          swm   = 304.94_wp * ( tbody - 36.6_wp ) * adu / 3600000.0_wp
          vpts  = 6.11_wp * 10.0_wp**( 7.45_wp * tsk / ( 235.0_wp + tsk ) )

          IF ( tbody <= 36.6_wp )  swm = 0.0_wp  !< no need for sweating

          sw = swm
          eswphy = - sw * l_v
          he     = 0.633_wp * hc / ( pair * c_p )
          fec    = 1.0_wp / ( 1.0_wp + 0.92_wp * hc * rcl )
          eswpot = he * ( vpa - vpts ) * adu * l_v * fec
          wetsk  = eswphy / eswpot

          IF ( wetsk > 1.0_wp )  wetsk = 1.0_wp
!
!--       Sweat production > evaporation?
          eswdif = eswphy - eswpot

          IF ( eswdif <= 0.0_wp )  esw = eswpot     !< Limit is evaporation
          IF ( eswdif > 0.0_wp )   esw = eswphy     !< Limit is sweat production
          IF ( esw  > 0.0_wp )     esw = 0.0_wp     !< Sweat can't be evaporated, no more cooling
                                                    !< effect
!
!--       Diffusion
          rdsk = 0.79_wp * 10.0_wp**7.0_wp
          rdcl = 0.0_wp
          ed   = l_v / ( rdsk + rdcl ) * adu * ( 1.0_wp - wetsk ) * ( vpa - vpts )
!
!--       Max vb
          vb1 = 34.0_wp - tsk
          vb2 = tcore(j) - 36.6_wp

          IF ( vb2 < 0.0_wp )  vb2 = 0.0_wp
          IF ( vb1 < 0.0_wp )  vb1 = 0.0_wp
          vb = ( 6.3_wp + 75.0_wp * vb2 ) / ( 1.0_wp + 0.5_wp * vb1 )
!
!--       Energy ballence
          enbal = int_heat + ed + ere + esw + csum + rsum + food
!
!--       Clothing temperature
          xx = 0.001_wp
          IF ( count1 == 0 )  xx = 1.0_wp
          IF ( count1 == 1 )  xx = 0.1_wp
          IF ( count1 == 2 )  xx = 0.01_wp
          IF ( count1 == 3 )  xx = 0.001_wp

          IF ( enbal > 0.0_wp )  tcl = tcl + xx
          IF ( enbal < 0.0_wp )  tcl = tcl - xx

          skipincreasecount = .FALSE.
          IF ( ( (enbal <= 0.0_wp )  .AND.  (enbal2 > 0.0_wp ) )  .OR.                             &
             ( ( enbal >= 0.0_wp )   .AND.  ( enbal2 < 0.0_wp ) ) )  THEN
             skipincreasecount = .TRUE.
          ELSE
             enbal2 = enbal
             count3 = count3 + 1
          ENDIF

          IF ( ( count3 > 200 )  .OR.  skipincreasecount )  THEN
             IF ( count1 < 3 )  THEN
                count1 = count1 + 1
                enbal2 = 0.0_wp
             ELSE
                EXIT
             ENDIF
          ENDIF
       ENDDO

       IF ( count1 == 3 )  THEN
          SELECT CASE ( j )
             CASE ( 2, 5)
                IF ( .NOT. ( ( tcore(j) >= 36.6_wp )  .AND.  ( tsk <= 34.050_wp ) ) )  CYCLE
             CASE ( 6, 1 )
                IF ( c(10) < 0.0_wp ) CYCLE
                IF ( .NOT. ( ( tcore(j) >= 36.6_wp )  .AND.  ( tsk > 33.850_wp ) ) )  CYCLE
             CASE ( 3 )
                IF ( .NOT. ( ( tcore(j) < 36.6_wp )  .AND.  ( tsk <= 34.000_wp ) ) )  CYCLE
             CASE ( 7 )
                IF ( .NOT. ( ( tcore(j) < 36.6_wp )  .AND.  ( tsk > 34.000_wp ) ) )  CYCLE
             CASE default
          END SELECT
       ENDIF

       IF ( ( j /= 4 )  .AND.  ( vb >= 91.0_wp ) )  CYCLE
       IF ( ( j == 4 )  .AND.  ( vb < 89.0_wp ) )  CYCLE
       IF ( vb > 90.0_wp ) vb = 90.0_wp
!
!--    Loses by water
       ws = sw * 3600.0_wp * 1000.0_wp
       IF ( ws > 2000.0_wp )  ws = 2000.0_wp
       wd = ed / l_v * 3600.0_wp * ( -1000.0_wp )
       wr = erel / l_v * 3600.0_wp * ( -1000.0_wp )

       wsum = ws + wr + wd

       RETURN
    ENDDO
 END SUBROUTINE heat_exch

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculate PET
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE pet_iteration( acl, adu, aeff, esw, facl, feff, int_heat, pair, rdcl, rdsk, rtv, ta,   &
                           tcl, tsk, pet_ij, vpts, wetsk )
!
!-- Input arguments:
    REAL(wp), INTENT( IN ) ::  acl       !< clothing surface area        (m²)
    REAL(wp), INTENT( IN ) ::  adu       !< Du-Bois area                 (m²)
    REAL(wp), INTENT( IN ) ::  esw       !< energy-loss through sweat evap. (W)
    REAL(wp), INTENT( IN ) ::  facl      !< surface area extension through clothing (factor)
    REAL(wp), INTENT( IN ) ::  feff      !< surface modification by posture (factor)
    REAL(wp), INTENT( IN ) ::  int_heat  !< internal heat production (W)
    REAL(wp), INTENT( IN ) ::  pair      !< air pressure                 (hPa)
    REAL(wp), INTENT( IN ) ::  rdcl      !< diffusion resistence of clothing (factor)
    REAL(wp), INTENT( IN ) ::  rdsk      !< diffusion resistence of skin (factor)
    REAL(wp), INTENT( IN ) ::  rtv       !< respiratory volume
    REAL(wp), INTENT( IN ) ::  ta        !< air temperature              (degree_C)
    REAL(wp), INTENT( IN ) ::  tcl       !< clothing temperature         (degree_C)
    REAL(wp), INTENT( IN ) ::  tsk       !< skin temperature             (degree_C)
    REAL(wp), INTENT( IN ) ::  vpts      !< sat. vapor pressure over skin (hPa)
    REAL(wp), INTENT( IN ) ::  wetsk     !< fraction of wet skin (dimensionless)
!
!-- Output arguments:
    REAL(wp), INTENT( OUT ) ::  aeff     !< effective surface area       (m²)
    REAL(wp), INTENT( OUT ) ::  pet_ij   !< PET                          (degree_C)
!
!-- Cconstants:
    REAL(wp), PARAMETER :: emcl =    0.95_wp      !< Longwave emission coef. of cloth
    REAL(wp), PARAMETER :: emsk =    0.99_wp      !< Longwave emission coef. of skin
    REAL(wp), PARAMETER :: po   = 1013.25_wp      !< Air pressure at sea level (hPa)
!
!-- Internal variables
    INTEGER ( iwp ) ::  count1        !< running index
    INTEGER ( iwp ) ::  i             !< running index

    REAL ( wp ) ::  cbare             !< Convection through bare skin
    REAL ( wp ) ::  cclo              !< Convection through clothing
    REAL ( wp ) ::  csum              !< Convection in total
    REAL ( wp ) ::  ed                !< Diffusion                      (W)
    REAL ( wp ) ::  enbal             !< Energy ballance                (W)
    REAL ( wp ) ::  enbal2            !< Energy ballance (last iteration cycle)
    REAL ( wp ) ::  ere               !< Energy ballance result         (W)
    REAL ( wp ) ::  erel              !< Latent energy ballance         (W)
    REAL ( wp ) ::  eres              !< Sensible respiratory heat flux (W)
    REAL ( wp ) ::  hc                !<
    REAL ( wp ) ::  rbare             !< Radiational loss of bare skin  (W/m²)
    REAL ( wp ) ::  rclo              !< Radiational loss of clothing   (W/m²)
    REAL ( wp ) ::  rsum              !< Radiational loss or gain       (W/m²)
    REAL ( wp ) ::  tex               !< Temperat. of exhaled air       (degree_C)
    REAL ( wp ) ::  vpex              !< Vapor pressure of exhaled air  (hPa)
    REAL ( wp ) ::  xx                !< Delta PET per iteration        (K)


    pet_ij = ta
    enbal2 = 0.0_wp

    DO  count1 = 0, 3
       DO  i = 1, 125  ! 500 / 4
          hc = 2.67_wp + 6.5_wp * 0.1_wp**0.67_wp
          hc = hc * ( pair / po )**0.55_wp
!
!--       Radiation
          aeff  = adu * feff
          rbare = aeff * ( 1.0_wp - facl ) * emsk * sigma_sb *                  &
                  ( ( pet_ij + degc_to_k )**4.0_wp - ( tsk + degc_to_k )**4.0_wp )
          rclo  = feff * acl * emcl * sigma_sb *                               &
                  ( ( pet_ij + degc_to_k )**4.0_wp - ( tcl + degc_to_k )**4.0_wp )
          rsum  = rbare + rclo
!
!--       Covection
          cbare = hc * ( pet_ij - tsk ) * adu * ( 1.0_wp - facl )
          cclo  = hc * ( pet_ij - tcl ) * acl
          csum  = cbare + cclo
!
!--       Diffusion
          ed = l_v / ( rdsk + rdcl ) * adu * ( 1.0_wp - wetsk ) * ( 12.0_wp - vpts )
!
!--       Respiration
          tex  = 0.47_wp * pet_ij + 21.0_wp
          eres = c_p * ( pet_ij - tex ) * rtv
          vpex = 6.11_wp * 10.0_wp**( 7.45_wp * tex / ( 235.0_wp + tex ) )
          erel = 0.623_wp * l_v / pair * ( 12.0_wp - vpex ) * rtv
          ere  = eres + erel
!
!--       Energy ballance
          enbal = int_heat + ed + ere + esw + csum + rsum
!
!--       Iteration concerning ta
          xx = 0.001_wp
          IF ( count1 == 0 )  xx = 1.0_wp
          IF ( count1 == 1 )  xx = 0.1_wp
          IF ( count1 == 2 )  xx = 0.01_wp
!           IF ( count1 == 3 )  xx = 0.001_wp
          IF ( enbal > 0.0_wp )  pet_ij = pet_ij - xx
          IF ( enbal < 0.0_wp )  pet_ij = pet_ij + xx
          IF ( ( enbal <= 0.0_wp )  .AND.  ( enbal2 > 0.0_wp ) )  EXIT
          IF ( ( enbal >= 0.0_wp )  .AND.  ( enbal2 < 0.0_wp ) )  EXIT

          enbal2 = enbal
       ENDDO
    ENDDO
 END SUBROUTINE pet_iteration

!
!-- UVEM specific subroutines

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Module-specific routine for new module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE uvem_solar_position

    USE control_parameters,                                                                        &
       ONLY:  latitude,                                                                            &
              longitude,                                                                           &
              time_since_reference_point

    IMPLICIT NONE

    INTEGER(iwp) ::  day_of_year = 0       !< day of year

    REAL(wp) ::  alpha         = 0.0_wp    !< solar azimuth angle in radiant
    REAL(wp) ::  declination   = 0.0_wp    !< declination
    REAL(wp) ::  dtor          = 0.0_wp    !< factor to convert degree to radiant
    REAL(wp) ::  js            = 0.0_wp    !< parameter for solar position calculation
    REAL(wp) ::  lat           = 52.39_wp  !< latitude
    REAL(wp) ::  lon           = 9.7_wp    !< longitude
    REAL(wp) ::  second_of_day = 0.0_wp    !< current second of the day
    REAL(wp) ::  thetar        = 0.0_wp    !< angle for solar zenith angle calculation
    REAL(wp) ::  thetasr       = 0.0_wp    !< angle for solar azimuth angle calculation
    REAL(wp) ::  zgl           = 0.0_wp    !< calculated exposure by direct beam
    REAL(wp) ::  woz           = 0.0_wp    !< calculated exposure by diffuse radiation
    REAL(wp) ::  wsp           = 0.0_wp    !< calculated exposure by direct beam


    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year,                     &
                        second_of_day = second_of_day )
    dtor = pi / 180.0_wp
    lat = latitude
    lon = longitude
!
!-- Calculation of js, necessary for calculation of equation of time (zgl) :
    js=  72.0_wp * ( REAL( day_of_year, KIND = wp ) + ( second_of_day / 86400.0_wp ) ) / 73.0_wp
!
!-- Calculation of equation of time (zgl):
    zgl = 0.0066_wp + 7.3525_wp * COS( ( js + 85.9_wp ) * dtor ) + 9.9359_wp *                     &
    COS( ( 2.0_wp * js + 108.9_wp ) * dtor ) + 0.3387_wp * COS( ( 3 * js + 105.2_wp ) * dtor )
!
!-- Calculation of apparent solar time woz:
    woz = ( ( second_of_day / 3600.0_wp ) - ( 4.0_wp * ( 15.0_wp - lon ) ) / 60.0_wp ) +           &
          ( zgl / 60.0_wp )
!
!-- Calculation of hour angle (wsp):
    wsp = ( woz - 12.0_wp ) * 15.0_wp
!
!-- Calculation of declination:
    declination = 0.3948_wp - 23.2559_wp * COS( ( js + 9.1_wp ) * dtor ) -                         &
                  0.3915_wp * COS( ( 2.0_wp * js + 5.4_wp ) * dtor ) - 0.1764_wp *                 &
                  COS( ( 3.0_wp * js + 26.0_wp ) * dtor )
!
!-- Calculation of solar zenith angle
    thetar  = ACOS( SIN( lat * dtor) * SIN( declination * dtor ) + COS( wsp * dtor ) *             &
              COS( lat * dtor ) * COS( declination * dtor ) )
    thetasr = ASIN( SIN( lat * dtor) * SIN( declination * dtor ) + COS( wsp * dtor ) *             &
              COS( lat * dtor ) * COS( declination * dtor ) )
    sza = thetar / dtor
!
!-- calculation of solar azimuth angle
    IF (woz <= 12.0_wp) alpha = pi - ACOS( ( SIN(thetasr) * SIN( lat * dtor ) -                    &
                                SIN( declination * dtor ) ) / ( COS(thetasr) * COS( lat * dtor ) ) )
    IF (woz > 12.0_wp)  alpha = pi + ACOS( ( SIN(thetasr) * SIN( lat * dtor ) -                    &
                                SIN( declination * dtor ) ) / ( COS(thetasr) * COS( lat * dtor ) ) )
    saa = alpha / dtor

 END SUBROUTINE uvem_solar_position


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Module-specific routine for new module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_calculate_uv_exposure

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !< loop index in x direction
    INTEGER(iwp) ::  j     !< loop index in y direction
    INTEGER(iwp) ::  szai  !< loop index for different sza values

    CALL uvem_solar_position

    IF (sza  >=  90)  THEN
       vitd3_exposure(:,:) = 0.0_wp
    ELSE

       DO  ai = 0, 35
          DO  zi = 0, 9
             projection_area_lookup_table(ai,zi) = uvem_projarea_f%var(clothing,zi,ai)
          ENDDO
       ENDDO
       DO  ai = 0, 35
          DO  zi = 0, 9
             integration_array(ai,zi) = uvem_integration_f%var(zi,ai)
          ENDDO
       ENDDO
       DO  ai = 0, 2
          DO  zi = 0, 90
             irradiance_lookup_table(ai,zi) = uvem_irradiance_f%var(zi,ai)
          ENDDO
       ENDDO
       DO  ai = 0, 35
          DO  zi = 0, 9
             DO  szai = 0, 90
                radiance_lookup_table(ai,zi,szai) = uvem_radiance_f%var(szai,zi,ai)
             ENDDO
          ENDDO
       ENDDO


!--    Rotate 3D-Model human to desired direction
       projection_area_temp( 0:35,:) = projection_area_lookup_table
       projection_area_temp(36:71,:) = projection_area_lookup_table
       IF (  .NOT.  turn_to_sun ) startpos_human = orientation_angle / 10.0_wp
       IF (       turn_to_sun ) startpos_human = saa / 10.0_wp
       DO  ai = 0, 35
          xfactor = ( startpos_human ) - INT( startpos_human )
          DO  zi = 0, 9
             projection_area(ai,zi) = ( projection_area_temp( 36 -                                 &
                                                               INT( startpos_human ) - 1 + ai , zi)&
                                        * ( xfactor ) )                                            &
                                      + ( projection_area_temp( 36 -                               &
                                                                INT( startpos_human ) + ai , zi)   &
                                        * ( 1.0_wp - xfactor ) )
          ENDDO
       ENDDO
!
!
!--    Interpolate to accurate Solar Zenith Angle
       DO  ai = 0, 35
          xfactor = ( sza )-INT( sza )
          DO  zi = 0, 9
             radiance_array(ai,zi) = ( radiance_lookup_table(ai, zi, INT( sza ) ) *                &
                                     ( 1.0_wp - xfactor) ) +                                       &
                                     ( radiance_lookup_table(ai,zi,INT( sza ) + 1) * xfactor )
          ENDDO
       ENDDO
       DO  iq = 0, 2
          irradiance(iq) = ( irradiance_lookup_table(iq, INT( sza ) ) * ( 1.0_wp - xfactor)) +     &
                           ( irradiance_lookup_table(iq, INT( sza ) + 1) * xfactor )
       ENDDO
!
!--    Interpolate to accurate Solar Azimuth Angle
       IF ( sun_in_south )  THEN
          startpos_saa_float = 180.0_wp / 10.0_wp
       ELSE
          startpos_saa_float = saa / 10.0_wp
       ENDIF
       radiance_array_temp( 0:35,:) = radiance_array
       radiance_array_temp(36:71,:) = radiance_array
       xfactor = (startpos_saa_float) - INT( startpos_saa_float )
       DO  ai = 0, 35
          DO  zi = 0, 9
             radiance_array(ai,zi) = ( radiance_array_temp(36 -                                    &
                                                           INT( startpos_saa_float ) - 1 + ai, zi) &
                                       * ( xfactor ) )                                             &
                                     + ( radiance_array_temp(36 -                                  &
                                                             INT( startpos_saa_float ) + ai, zi)   &
                                       * ( 1.0_wp - xfactor ) )
          ENDDO
       ENDDO

!
!--    Calculate Projectionarea for direct beam
       projection_area_direct_temp( 0:35,:) = projection_area
       projection_area_direct_temp(36:71,:) = projection_area
       yfactor = ( sza / 10.0_wp ) - INT( sza / 10.0_wp )
       xfactor = ( startpos_saa_float ) - INT( startpos_saa_float )
       projection_area_direct_beam = ( projection_area_direct_temp( INT(startpos_saa_float)    ,INT(sza/10.0_wp)  ) * &
                                     ( 1.0_wp - xfactor ) * ( 1.0_wp - yfactor ) ) +                                  &
                                     ( projection_area_direct_temp( INT(startpos_saa_float) + 1,INT(sza/10.0_wp)  ) * &
                                     (          xfactor ) * ( 1.0_wp - yfactor ) ) +                                  &
                                     ( projection_area_direct_temp( INT(startpos_saa_float)    ,INT(sza/10.0_wp)+1) * &
                                     ( 1.0_wp - xfactor ) * (          yfactor ) ) +                                  &
                                     ( projection_area_direct_temp( INT(startpos_saa_float) + 1,INT(sza/10.0_wp)+1) * &
                                     (          xfactor ) * (          yfactor ) )


       DO  i = nxl, nxr
          DO  j = nys, nyn
!
!--          Extract obstruction from IBSET-Integer_Array
             IF (consider_obstructions )  THEN
                obstruction_temp1 = building_obstruction_f%var_3d(:,j,i)
                IF ( obstruction_temp1(0)  /=  9 )  THEN
                   DO  pobi = 0, 44
                      DO  bi = 0, 7
                         IF ( BTEST( obstruction_temp1(pobi), bi )  .EQV.  .TRUE.)  THEN
                            obstruction_temp2( ( pobi * 8 ) + bi ) = 1
                         ELSE
                            obstruction_temp2( ( pobi * 8 ) + bi ) = 0
                         ENDIF
                      ENDDO
                   ENDDO
                   DO  zi = 0, 9
                      obstruction(:,zi) = obstruction_temp2( zi * 36 :( zi * 36) + 35 )
                   ENDDO
                ELSE
                   obstruction(:,:) = 0
                ENDIF
             ENDIF
!
!--          Calculated human exposure
             diffuse_exposure = SUM( radiance_array * projection_area * integration_array *        &
                                     obstruction )

             obstruction_direct_beam = obstruction( NINT( startpos_saa_float),                     &
                                                    NINT( sza / 10.0_wp ) )
             IF (sza  >=  89.99_wp)  THEN
                sza = 89.99999_wp
             ENDIF
!
!--          Calculate direct normal irradiance (direct beam)
             direct_exposure = ( irradiance(1) / COS( pi * sza / 180.0_wp ) ) *                    &
                               projection_area_direct_beam * obstruction_direct_beam

             vitd3_exposure(j,i) = ( diffuse_exposure + direct_exposure ) / 1000.0_wp * 70.97_wp
!            unit = international units vitamin D per second
          ENDDO
       ENDDO
    ENDIF

 END SUBROUTINE bio_calculate_uv_exposure


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sampling of biometeorological variables along customized measurement coordinates.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE bio_vm_sampling( variable, var_atmos, i_atmos, j_atmos, k_atmos, ns_atmos,             &
                             var_soil, i_soil, j_soil, k_soil, ns_soil, sampled )

    CHARACTER(LEN=*) ::  variable  !< treated variable

    INTEGER(iwp) ::  i         !< grid index in x-direction
    INTEGER(iwp) ::  j         !< grid index in y-direction
    INTEGER(iwp) ::  m         !< running index over all virtual observation coordinates
    INTEGER(iwp) ::  ns_atmos  !< number of sampling points for atmosphere and surface variables
    INTEGER(iwp) ::  ns_soil   !< number of sampling points for soil variables

    INTEGER(iwp), DIMENSION(1:ns_atmos) ::  i_atmos  !< sampling index in x-direction for atmosphere variables
    INTEGER(iwp), DIMENSION(1:ns_atmos) ::  j_atmos  !< sampling index in y-direction for atmosphere variables
    INTEGER(iwp), DIMENSION(1:ns_atmos) ::  k_atmos  !< sampling index in z-direction for atmosphere variables

    INTEGER(iwp), DIMENSION(1:ns_soil) ::   i_soil   !< sampling index in x-direction for soil variables
    INTEGER(iwp), DIMENSION(1:ns_soil) ::   j_soil   !< sampling index in y-direction for soil variables
    INTEGER(iwp), DIMENSION(1:ns_soil) ::   k_soil   !< sampling index in z-direction for soil variables

    LOGICAL ::  sampled !< flag indicating whether a variable has been sampled

    REAL(wp), DIMENSION(1:ns_atmos) ::  var_atmos  !< array to store atmosphere variables

    REAL(wp), DIMENSION(1:ns_soil ) ::  var_soil   !< array to store soil variables


    SELECT CASE ( TRIM( variable ) )
!
!--    Mean radiant temperature.
       CASE ( 't_mrt' )
          IF ( ALLOCATED( tmrt_grid ) )  THEN
             DO  m = 1, ns_atmos
                j = j_atmos(m)
                i = i_atmos(m)
                var_atmos(m) = tmrt_grid(j,i)
             ENDDO
             sampled = .TRUE.
          ENDIF
!
!--    Perceived temperature.
       CASE ( 't_perceived' )
          IF ( ALLOCATED( perct ) )  THEN
             DO  m = 1, ns_atmos
                j = j_atmos(m)
                i = i_atmos(m)
                var_atmos(m) = perct(j,i)
             ENDDO
             sampled = .TRUE.
          ENDIF
!
!--    Physiological equivalent temperature.
       CASE ( 't_pet' )
          IF ( ALLOCATED( pet ) )  THEN
             DO  m = 1, ns_atmos
                j = j_atmos(m)
                i = i_atmos(m)
                var_atmos(m) = pet(j,i)
             ENDDO
             sampled = .TRUE.
          ENDIF
!
!--    UTCI.
       CASE ( 't_utci' )
          IF ( ALLOCATED( utci ) )  THEN
             DO  m = 1, ns_atmos
                j = j_atmos(m)
                i = i_atmos(m)
                var_atmos(m) = utci(j,i)
             ENDDO
             sampled = .TRUE.
          ENDIF

       CASE DEFAULT

    END SELECT
!
!-- Avoid compiler warning for unused variables by constructing an if condition which is never
!-- fulfilled.
    IF ( .FALSE.  .AND.  ns_soil < 0 )  THEN
       i_soil = i_soil
       j_soil = j_soil
       k_soil = k_soil
       k_atmos = k_atmos
       var_soil = var_soil
    ENDIF

 END SUBROUTINE bio_vm_sampling

 END MODULE biometeorology_mod
