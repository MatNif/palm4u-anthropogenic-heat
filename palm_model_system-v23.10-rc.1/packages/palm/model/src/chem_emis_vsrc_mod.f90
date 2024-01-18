!> @file chem_emis_vsrc_mod.f90
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
! Copyright 2018-2021 Leibniz Universitaet Hannover
! Copyright 2018-2021 Karlsruhe Institute of Technology
! Copyright 2018-2021 Freie Universitaet Berlin
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Edward C. Chan
! @author Ilona Jäkel
!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of emission volume sources data structures and API
!--------------------------------------------------------------------------------------------------!

 MODULE chem_emis_vsrc_mod

    USE arrays_3d,                                                                                 &
        ONLY:  hyp,                                                                                &
               pt,                                                                                 &
               exner

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  p_0,                                                                                &
               rd_d_cp,                                                                            &
               rgas_univ

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar,                                                                               &
               spc_names

    USE kinds

    IMPLICIT NONE
    SAVE 
    PRIVATE

!
!-- data structure relating volume source locations for an individual emission mode
    TYPE, PUBLIC ::  chem_emis_vsrc_pos
       INTEGER(iwp)      ::  i
       INTEGER(iwp)      ::  j
       INTEGER(iwp)      ::  k
       INTEGER(KIND=idp) ::  key
    END TYPE chem_emis_vsrc_pos

!
!-- default values for initialization
    INTEGER(iwp),      PARAMETER ::  default_index =    -1         !< ijk index
    INTEGER(KIND=idp), PARAMETER ::  default_key   = -9999         !< hash key
    REAL(KIND=dp),     PARAMETER ::  default_value =     0.0       !< volume source value

!
!-- netCDF file attributes
    CHARACTER(LEN=*), PARAMETER ::  nc_dim_nvsrc       = 'nvsrc'
    CHARACTER(LEN=*), PARAMETER ::  nc_var_vsrc_i      = 'vsrc_i'
    CHARACTER(LEN=*), PARAMETER ::  nc_var_vsrc_j      = 'vsrc_j'
    CHARACTER(LEN=*), PARAMETER ::  nc_var_vsrc_k      = 'vsrc_k'
    CHARACTER(LEN=*), PARAMETER ::  nc_var_vsrc_prefix = 'vsrc_'

!
!-- volume source locations across all activated emission modes
    INTEGER(iwp)                                   ::  num_vsrc    !< # unique volume source locations
    INTEGER(KIND=idp), ALLOCATABLE, DIMENSION(:)   ::  vsrc_key    !< volume source hash key (kji based)
    REAL(KIND=dp),     ALLOCATABLE, DIMENSION(:,:) ::  vsrc_value  !< vol. source values for all KPP variable species

!
!-- interface to object construction and destruction
    INTERFACE chem_emis_vsrc_init_pos_vector
       MODULE PROCEDURE chem_emis_vsrc_init_pos_vector
    END INTERFACE chem_emis_vsrc_init_pos_vector

    INTERFACE chem_emis_vsrc_cleanup
       MODULE PROCEDURE chem_emis_vsrc_cleanup
    END INTERFACE chem_emis_vsrc_cleanup

!
!-- interface to index-key conversion
    INTERFACE chem_emis_vsrc_key_from_indices
       MODULE PROCEDURE chem_emis_vsrc_key_from_indices
    END INTERFACE chem_emis_vsrc_key_from_indices

    INTERFACE chem_emis_vsrc_indices_from_key
       MODULE PROCEDURE chem_emis_vsrc_indices_from_key
    END INTERFACE chem_emis_vsrc_indices_from_key

!
!-- interace to volume source manipulation
    INTERFACE chem_emis_vsrc_append_positions
       MODULE PROCEDURE chem_emis_vsrc_append_positions
    END INTERFACE chem_emis_vsrc_append_positions

    INTERFACE chem_emis_vsrc_reset_source
       MODULE PROCEDURE chem_emis_vsrc_reset_source
    END INTERFACE chem_emis_vsrc_reset_source

    INTERFACE chem_emis_vsrc_update_source
       MODULE PROCEDURE chem_emis_vsrc_update_source_by_key
       MODULE PROCEDURE chem_emis_vsrc_update_source_by_species
    END INTERFACE chem_emis_vsrc_update_source

    INTERFACE chem_emis_vsrc_assign_source
       MODULE PROCEDURE chem_emis_vsrc_assign_source
    END INTERFACE chem_emis_vsrc_assign_source

!
!-- public variables
    PUBLIC ::  default_index, default_value,                                                       &
               nc_dim_nvsrc, nc_var_vsrc_i, nc_var_vsrc_j, nc_var_vsrc_k,                          &
               nc_var_vsrc_prefix

!
!-- public methods
    PUBLIC ::  chem_emis_vsrc_append_positions, chem_emis_vsrc_assign_source,                      &
               chem_emis_vsrc_cleanup, chem_emis_vsrc_indices_from_key,                            &
               chem_emis_vsrc_init_pos_vector, chem_emis_vsrc_key_from_indices,                    &
               chem_emis_vsrc_reset_source, chem_emis_vsrc_update_source
               

 CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  PUBLIC METHODS
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> allocate memory for and initialize position vector
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE chem_emis_vsrc_init_pos_vector( pos_vector, npos )

    IMPLICIT NONE

!
!-- arguments in order of appearance
    TYPE(chem_emis_vsrc_pos), INTENT(OUT), ALLOCATABLE, DIMENSION(:) ::  pos_vector  !< volume source position vector
    INTEGER(iwp),             INTENT(IN)                             ::  npos        !< element count

!
!-- local variables
    INTEGER(iwp) ::  k   !< counter

    IF ( ALLOCATED( pos_vector )  .OR.  ( npos<=0 ) )  RETURN   ! sanity check
    ALLOCATE ( pos_vector(npos) )

    DO  k = 1, npos
       pos_vector(k)%i   = default_index
       pos_vector(k)%j   = default_index
       pos_vector(k)%k   = default_index
       pos_vector(k)%key = default_key     
    END DO

 END SUBROUTINE chem_emis_vsrc_init_pos_vector


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> deallocate memory for overall volume source positions and source term values
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE chem_emis_vsrc_cleanup( )

    IMPLICIT NONE

!
!-- local variables
    LOGICAL, PARAMETER ::  dummy = .FALSE.

    IF ( dummy ) CALL diag_print_vsrc ( )
    IF ( ALLOCATED( vsrc_key)   )  DEALLOCATE( vsrc_key   )
    IF ( ALLOCATED( vsrc_value) )  DEALLOCATE( vsrc_value )

 END SUBROUTINE chem_emis_vsrc_cleanup


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> convert ijk indices to hash key and vice versa
!> NOTE: to conform with loop conventions in prognostic equations, key will be generated 
!>       in the order kji, k being the fasting changing index
!--------------------------------------------------------------------------------------------------!

!
!-- get key from ijk indices
 FUNCTION chem_emis_vsrc_key_from_indices ( i, j, k ) RESULT  ( key )

    USE kinds

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               ny

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(KIND=idp)            ::  key  !< hash key
    INTEGER(iwp),     INTENT(IN) ::  i    !< i index
    INTEGER(iwp),     INTENT(IN) ::  j    !< j index
    INTEGER(iwp),     INTENT(IN) ::  k    !< k index

    key = i + ( nx + 1 ) * ( j + ( ny + 1 ) * k )

 END FUNCTION chem_emis_vsrc_key_from_indices

!
!-- get indices from key

 SUBROUTINE chem_emis_vsrc_indices_from_key ( i, j, k, key )

    USE kinds

    USE indices,                                                                                   &
        ONLY:  nx,                                                                                 &
               ny

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp),      INTENT(OUT) ::  i      !< i index
    INTEGER(iwp),      INTENT(OUT) ::  j      !< j index
    INTEGER(iwp),      INTENT(OUT) ::  k      !< k index
    INTEGER(KIND=idp), INTENT(IN)  ::  key    !< hash key

!
!-- local variables
    INTEGER(iwp) ::  mod_k               !< temporary variable
    INTEGER(iwp) ::  slice_size          !< temporary variable
    INTEGER(iwp) ::  nx1                 !< nx + 1

    nx1        = nx + 1
    slice_size = nx1 * ( ny + 1 )

    k          = key / slice_size
    mod_k      = key - ( slice_size * k )
    j          = mod_k / nx1
    i          = mod_k - ( nx1 * j )

 END SUBROUTINE chem_emis_vsrc_indices_from_key


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> append volume source positions from position vector specified by emission mode
!> volume source values will be initialized accordingly
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE chem_emis_vsrc_append_positions ( pos_vector )

    USE kinds

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    IMPLICIT NONE

!
!-- arguments in order of appearance
    TYPE(chem_emis_vsrc_pos), INTENT(IN), DIMENSION(:) ::  pos_vector

!
!-- local variables
    INTEGER(iwp)     ::  buf_size            !< temporary array size
    INTEGER(iwp)     ::  num_vsrc_new        !< temporary array size
    INTEGER(iwp)     ::  k                   !< generic counter

    INTEGER(KIND=idp), ALLOCATABLE, DIMENSION(:) ::  buf_key_vector      !< temporary hash key vector

    buf_size = SIZE( pos_vector )

!
!-- for first time assignment (i.e. vsrc_key / _value not allocated)
!-- simply transfer everything directly from position vector sorted
    IF ( .NOT. vsrc_is_allocated (  ) )  THEN
        num_vsrc = buf_size
        CALL init_vsrc ( num_vsrc, nvar )
        DO  k = 1, num_vsrc
           vsrc_key(k) = pos_vector(k)%key
        END DO
        CALL qsort ( vsrc_key )
        RETURN
    END IF

!
!-- if vsrc_key / _value already present
    buf_size = buf_size + num_vsrc          ! create buffer and transfer existing hash keys
    ALLOCATE( buf_key_vector(buf_size) )
    buf_key_vector(1:num_vsrc) = vsrc_key
    
    num_vsrc_new = num_vsrc                 ! append unique keys from position vector

    DO  k = 1, SIZE( pos_vector )
       IF ( locate_vsrc_key_index ( buf_key_vector, pos_vector(k)%key ) < 0 )  THEN
          num_vsrc_new = num_vsrc_new + 1
          buf_key_vector(num_vsrc_new) = pos_vector(k)%key
       END IF
    END DO

    num_vsrc = num_vsrc_new                 ! reassign buffer to vsrc_key / _value and resort
    CALL chem_emis_vsrc_cleanup ( )
    CALL init_vsrc ( num_vsrc, nvar )
    vsrc_key = buf_key_vector(1:num_vsrc)
    CALL qsort ( vsrc_key )

    DEALLOCATE( buf_key_vector )           ! cleaning up
    
 END SUBROUTINE chem_emis_vsrc_append_positions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> resets all source values to zero at the beginning of time integration
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE chem_emis_vsrc_reset_source ( )

    IMPLICIT NONE

    IF ( num_vsrc > 0 )  vsrc_value = default_value

 END SUBROUTINE chem_emis_vsrc_reset_source


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates a single volume source value from each emission mode for a given species
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE chem_emis_vsrc_update_source_by_key ( ispecies, key, vsrc )

    USE kinds

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp),      INTENT(IN) ::  ispecies  !< species index
    INTEGER(KIND=idp), INTENT(IN) ::  key       !< hash key
    REAL(KIND=dp),     INTENT(IN) ::  vsrc      !< volume source value

!
!-- local variables
    INTEGER(iwp) ::  k   !< position index

    IF ( num_vsrc == 0 )                      RETURN  ! sanity checks
    IF ( (ispecies<1)  .OR.  (ispecies>nvar) )  RETURN

    k = locate_vsrc_key_index ( vsrc_key, key )
    IF ( k /= default_index )  vsrc_value(ispecies,k) = vsrc_value(ispecies,k) + vsrc

 END SUBROUTINE chem_emis_vsrc_update_source_by_key


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> updates all volume source values from each emission mode for a given species
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE chem_emis_vsrc_update_source_by_species ( ispecies, pos_vector, vsrc_vector )

    USE kinds

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp),             INTENT(IN)               ::  ispecies     !< species index
    TYPE(chem_emis_vsrc_pos), INTENT(IN), DIMENSION(:) ::  pos_vector   !< emission volume source positions
    REAL(KIND=dp),            INTENT(IN), DIMENSION(:) ::  vsrc_vector  !< volume sources

!
!-- local variables
    INTEGER(iwp) ::  k  !< counter

    IF ( num_vsrc == 0 )                          RETURN   ! sanity checks
    IF ( SIZE(pos_vector) /= SIZE(vsrc_vector) )  RETURN
    IF ( (ispecies<1)  .OR.  (ispecies>nvar)   )  RETURN

    DO  k = 1, SIZE(pos_vector)
       CALL chem_emis_vsrc_update_source_by_key ( ispecies, pos_vector(k)%key, vsrc_vector(k) )
    END DO

 END SUBROUTINE chem_emis_vsrc_update_source_by_species


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> assigns source terms directly to prognostic equation variables
!> interface to chemistry_model_mod > chem_non_advective_processes*
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE chem_emis_vsrc_assign_source ( i, j )

    USE kinds

    USE statistics,                                                                                &
        ONLY:  weight_pres

    USE indices,                                                                                   &
        ONLY:  nzb,                                                                                &
               nzt

    USE control_parameters,                                                                        &
        ONLY:  dt_3d,                                                                              &
               intermediate_timestep_count

    USE chem_modules,                                                                              &
        ONLY:  chem_species,                                                                       &
               call_chem_at_all_substeps

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp), INTENT(IN) ::  i                        !< horizontal cell index
    INTEGER(iwp), INTENT(IN) ::  j                        !< horizontal cell index

!
!-- local variables
    INTEGER(iwp) ::  ispecies                             !< counter
    INTEGER(iwp) ::  k                                    !< counter
    INTEGER(iwp) ::  nzb1                                 !< first unstaggerd vertical grid index

    LOGICAL      ::  is_gas_phase                         !< whether species is particulate

    REAL(KIND=dp) ::  dt_chem                             !< chemical time step
    REAL(KIND=dp) ::  source_term                         !< source term value
    REAL(wp), DIMENSION(nzb:nzt+1) ::  temperature        !< columnar thermodynamic temperature
    REAL(wp), DIMENSION(nzb:nzt+1) ::  conv_ratio         !< conversion factor mol/m3 to ppmv
                                                          !< for gas phase species
    REAL(wp), PARAMETER            ::  fr2ppm = 1.0E6_wp  !< conversion molar fraction to ppmv
                                                          !< for gas phase species

    IF ( num_vsrc == 0 )  RETURN       ! get out right away if there are no volume sources

    nzb1 = nzb + 1

    dt_chem = dt_3d
    IF ( call_chem_at_all_substeps )  dt_chem = dt_chem * weight_pres(intermediate_timestep_count)

!
!-- calculate columnar conversion factor from mol/m3 to ppmv for gas phase species
!-- (PM and pollen species are exempt)

    temperature(nzb:nzt+1) = pt(nzb:nzt+1,j,i) * exner(nzb:nzt+1)
    conv_ratio(nzb:nzt+1)  = fr2ppm * rgas_univ * temperature(nzb:nzt+1) / hyp(nzb:nzt+1)

    DO  ispecies = 1, nvar             ! source terms for fixed species are zero

       is_gas_phase = ( .NOT. (                                                                     &
                         ( TRIM( spc_names(ispecies)(1:2) ) == "PM"   ) .OR.                        &
                         ( TRIM( spc_names(ispecies)(1:4) ) == "POL_" )                             &
                      ) )

       DO  k = nzb1,nzt
          source_term = dt_chem * fetch_species_source_by_indices( ispecies, i, j, k )
          IF  ( is_gas_phase )  source_term = source_term * conv_ratio(k)
          chem_species(ispecies)%conc(k,j,i) = chem_species(ispecies)%conc(k,j,i) + source_term
       ENDDO
    ENDDO

 END SUBROUTINE chem_emis_vsrc_assign_source


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  PRIVATE METHODS
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> check if all vsrc arrays have been allocated
!--------------------------------------------------------------------------------------------------!

 FUNCTION vsrc_is_allocated ( )  RESULT  ( is_allocated )

    IMPLICIT NONE

!
!-- arguments in order of appearance
    LOGICAL ::  is_allocated

    is_allocated = ALLOCATED(vsrc_key)  .AND.  ALLOCATED(vsrc_value)

 END FUNCTION vsrc_is_allocated


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> allocate memory for vsrc_ variables in the module and set default values
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE init_vsrc ( num_vsrc_pos, num_species )

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp), INTENT(IN) ::  num_vsrc_pos, num_species   !< # vector sizes

    IF ( (num_vsrc_pos<=0)  .OR.  (num_species<=0) )  RETURN   ! sanity check

    IF ( .NOT. ALLOCATED( vsrc_key ) )  THEN
       ALLOCATE( vsrc_key(num_vsrc_pos) )
       vsrc_key = default_key
    END IF

    IF ( .NOT. ALLOCATED( vsrc_value ) )  THEN
       ALLOCATE( vsrc_value(num_species,num_vsrc_pos) )
       vsrc_value = default_value
    END IF

 END SUBROUTINE init_vsrc


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> returns source associated with ijk indicies
!> defaults to zero in case of failure or error
!--------------------------------------------------------------------------------------------------!

 FUNCTION fetch_species_source_by_indices ( ispecies, i, j, k )  RESULT  ( src )

    USE kinds

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nys,                                                                                &
               nyn,                                                                                &
               nzb,                                                                                &
               nzt

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar

    IMPLICIT NONE

!
!-- arguments in order of appearance
    REAL(KIND=dp)            ::  src       !< volume source
    INTEGER(iwp), INTENT(IN) ::  ispecies  !< species index
    INTEGER(iwp), INTENT(IN) ::  i         !< cell index
    INTEGER(iwp), INTENT(IN) ::  j         !< cell index
    INTEGER(iwp), INTENT(IN) ::  k         !< cell index

!
!-- local variables
    INTEGER(iwp)      ::  isrc      !< volume source location
    INTEGER(KIND=idp) ::  key       !< hash key


    src = default_value                                   ! default return volume source value to zero

    IF ( ispecies > nvar )  RETURN                        ! get out if species not in mechanism

    IF ( ((i<nxl)  .OR.  (i>nxr))  .OR.       &           ! get out if cell indices out of range
         ((j<nys)  .OR.  (j>nyn))  .OR.       &
         ((k<nzb)  .OR.  (k>nzt)) )  RETURN

    key = chem_emis_vsrc_key_from_indices (i,j,k)     ! get out if key out of key location set
    IF ( (key<vsrc_key(1))  .OR.  (key>vsrc_key(num_vsrc)) )  RETURN

    isrc = locate_vsrc_key_index ( vsrc_key, key )        ! and only then fetch volume source
    IF ( isrc /= default_index )  src = vsrc_value(ispecies,isrc)

 END FUNCTION fetch_species_source_by_indices


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> in-place quick sort for large arrays using median of three pivoting method
!> reverts to insertion sort at sufficiently small partition size (hardcoded to 4 elements)
!> this is the same as the C built-in qsort function, without the void* and comp function
!--------------------------------------------------------------------------------------------------!

 RECURSIVE SUBROUTINE qsort ( x, i0, i1 )

    USE kinds
    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(KIND=idp), DIMENSION(:)      ::  x
    INTEGER(iwp), INTENT(IN), OPTIONAL   ::  i0
    INTEGER(iwp), INTENT(IN), OPTIONAL   ::  i1

!
!-- local variables
    INTEGER(iwp) ::  j0        !< partition bound index
    INTEGER(iwp) ::  j1        !< partition bound index

    INTEGER(iwp) ::  k0        !< median of three pivot indices
    INTEGER(iwp) ::  k1        !< median of three pivot indices
    INTEGER(iwp) ::  kp        !< median of three pivot indices

    INTEGER(iwp), PARAMETER ::  kswitch = 4   !< partition size threshold for switching to insert sort

    INTEGER(KIND=idp)  ::  x0        !< median of three values
    INTEGER(KIND=idp)  ::  x1        !< median of three values
    INTEGER(KIND=idp)  ::  xp        !< median of three values
    INTEGER(KIND=idp)  ::  xtemp     !< swap variable

!
!-- define sort region
    k0 = 1
    k1 = SIZE(x)

    IF ( PRESENT(i0) )  k0 = i0
    IF ( PRESENT(i1) )  k1 = i1

!
!-- check if insertion sort should be used
    IF ( (k1-k0) <= kswitch )  THEN
       CALL isort ( x, k0, k1 )
       RETURN
    END IF

!
!-- median of three pivot construction
    kp = ( k0 + k1 ) / 2
    x0 = MIN( x(k0), x(kp), x(k1) )
    x1 = MAX( x(k0), x(kp), x(k1) )
    xp = MAX( MIN(x(k0),x(k1)), MIN(x(kp),MAX(x(k0),x(k1))) )

    x(k0) = x0   ! rearrange in sorted order
    x(k1) = x1
    x(kp) = xp

!
!-- partition sort bound around pivot
    j0 = k0
    j1 = k1

    DO

       DO WHILE ( x(j0) < xp )
          j0 = j0 + 1
       END DO

       DO WHILE ( x(j1) > xp )
          j1 = j1 - 1
       END DO

       IF ( j0 >= j1 )  EXIT

       xtemp = x(j0)
       x(j0) = x(j1)
       x(j1) = xtemp

       j0 = j0 + 1
       j1 = j1 - 1      

    END DO

    IF ( j0 == j1 ) j0 = j0 + 1

!
!-- further sort subpartitions
    CALL qsort ( x, k0, j0-1 )
    CALL qsort ( x, j0, k1   )

 END SUBROUTINE qsort


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> in-place insertion sort for small arrays
!--------------------------------------------------------------------------------------------------!

 RECURSIVE SUBROUTINE isort ( x, i0, i1 )

    USE kinds

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(KIND=idp), DIMENSION(:) ::  x
    INTEGER(iwp), INTENT(IN), OPTIONAL   ::  i0
    INTEGER(iwp), INTENT(IN), OPTIONAL   ::  i1

!
!-- local variables
    INTEGER(iwp)      ::  k       !< search index
    INTEGER(iwp)      ::  k0      !< search index
    INTEGER(iwp)      ::  k1      !< search index
    INTEGER(KIND=idp) ::  xtemp   !< swap variable

!
!- define sort region
    k0 = 1
    k1 = SIZE(x)

    IF ( PRESENT(i0) )  k0 = i0
    IF ( PRESENT(i1) )  k1 = i1

    IF ( k1 == k0 )  RETURN       ! termination condition

    xtemp = x(k1)                 ! successively reduce search bound
    k     = k1 - 1
    CALL isort ( x, k0, k )

    DO WHILE ( k >= k0 )          ! element by element comparison
       IF ( x(k) <= xtemp ) EXIT
       x(k+1) = x(k)
       k = k - 1
    END DO

    x(k+1) = xtemp

 END SUBROUTINE isort


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> locate index of volume source key given hash key
!> returns -1 if no match is found
!--------------------------------------------------------------------------------------------------!

 RECURSIVE FUNCTION locate_vsrc_key_index ( x, xk, i0, i1 ) RESULT ( k )

    USE kinds

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp)                                ::  k
    INTEGER(KIND=idp), INTENT(IN), DIMENSION(:) ::  x
    INTEGER(KIND=idp), INTENT(IN)               ::  xk
    INTEGER(iwp),      INTENT(IN), OPTIONAL     ::  i0
    INTEGER(iwp),      INTENT(IN), OPTIONAL     ::  i1

!
!-- local variables
    INTEGER(iwp) ::  k0     !< start, end and center index of search region
    INTEGER(iwp) ::  k1     !< start, end and center index of search region

!
!-- define search region
    k0 = 1
    k1 = SIZE(x)

    IF ( PRESENT(i0) )  k0 = i0
    IF ( PRESENT(i1) )  k1 = i1

!
!-- termination conditions
    k = default_index  ! no match
    IF ( k0 > k1 )     RETURN

    k = ( k0 + k1 ) / 2
    IF ( x(k) == xk )  RETURN

!
!-- bisection search
    IF ( x(k) > xk )  THEN       ! key in left half of search bound
       k = locate_vsrc_key_index ( x, xk, k0,  k-1 )
    ELSE                         ! key in right half of search bound
       k = locate_vsrc_key_index ( x, xk, k+1, k1  )
    END IF

 END FUNCTION locate_vsrc_key_index


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> diagnostic function for printing out source term values
!--------------------------------------------------------------------------------------------------!

 SUBROUTINE diag_print_vsrc ( ispecies, ipos )

    USE chem_gasphase_mod,                                                                         &
        ONLY:  nvar,                                                                               &
               spc_names

    IMPLICIT NONE

!
!-- arguments in order of appearance
    INTEGER(iwp), INTENT(IN), OPTIONAL ::  ispecies
    INTEGER(iwp), INTENT(IN), OPTIONAL ::  ipos

!
!-- local variables
    INTEGER(iwp)      ::  i        !< cell index
    INTEGER(iwp)      ::  j        !< cell index
    INTEGER(iwp)      ::  k        !< cell index
    INTEGER(iwp)      ::  kpos
    INTEGER(iwp)      ::  ksp
    INTEGER(iwp)      ::  pos0     !< display bound
    INTEGER(iwp)      ::  pos1     !< display bound
    INTEGER(iwp)      ::  sp0      !< display bound
    INTEGER(iwp)      ::  sp1      !< display bound

    INTEGER(KIND=idp) ::  key      !< hash key

!
!-- assign species index
    sp0  = 1
    sp1  = nvar 
    pos0 = 1
    pos1 = num_vsrc

    IF ( PRESENT(ispecies) )  THEN
       sp0 = ispecies
       sp1 = ispecies
    END IF 

    IF ( PRESENT(ipos) )  THEN
       pos0 = ipos
       pos1 = ipos
    END IF 

!
!-- print
    PRINT *, ' *** volume source positions diagnostics:'
    DO  kpos = pos0, pos1
       CALL chem_emis_vsrc_indices_from_key( i, j, k, vsrc_key(kpos) )  ! get ijk from key
       key = chem_emis_vsrc_key_from_indices( i, j, k )                 ! see if key matches
       WRITE (*,*) kpos, '  (', i, j, k, ' )  ', vsrc_key(kpos), vsrc_key(kpos) - key
    END DO

    WRITE (*, fmt="(4x)", advance="no") 
    DO  ksp = sp0, sp1
       WRITE (*, fmt="(A10)", advance="no") spc_names(ksp)
    END DO
    WRITE (*,*)

    DO  kpos = pos0, pos1
       WRITE (*, fmt="(i2,2x)", advance="no") kpos
 
       DO  ksp = sp0, sp1
          WRITE (*, fmt="(E9.2,2x)", advance="no") vsrc_value(ksp,kpos)
       END DO
       WRITE (*,*)
    END DO

 END SUBROUTINE diag_print_vsrc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  END OF MODULE
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 END MODULE chem_emis_vsrc_mod
