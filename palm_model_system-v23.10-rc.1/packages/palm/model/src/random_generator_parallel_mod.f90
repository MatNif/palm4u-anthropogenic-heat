!> @file random_generator_parallel_mod.f90
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
!> This module contains and supports the random number generating routine ran_parallel.
!> ran_parallel returns a uniform random deviate between 0.0 and 1.0 (exclusive of the end point
!> values). Moreover, it contains a routine returning a normally distributed random number with mean
!> zero and unity standard deviation. Additionally, it provides the generator with five integers for
!> use as initial state space. The first tree integers (iran, jran, kran) are maintained as non
!> negative values, while the last two (mran, nran) have 32-bit nonzero values. Also provided by
!> this module is support for initializing or reinitializing the state space to a desired standard
!> sequence number, hashing the initial values to random values and allocating and deallocating the
!> internal workspace.
!>
!> This routine is taken from the "numerical recipies vol. 2"
!--------------------------------------------------------------------------------------------------!
 MODULE random_generator_parallel

    USE kinds

    IMPLICIT NONE

    INTEGER(isp), SAVE ::  lenran = 0  !<
    INTEGER(isp), SAVE ::  seq = 0     !<
    INTEGER(isp), SAVE ::  iran0       !<
    INTEGER(isp), SAVE ::  jran0       !<
    INTEGER(isp), SAVE ::  kran0       !<
    INTEGER(isp), SAVE ::  mran0       !<
    INTEGER(isp), SAVE ::  nran0       !<
    INTEGER(isp), SAVE ::  rans        !<

    INTEGER(isp), DIMENSION(:, :), POINTER, SAVE ::  ranseeds  !<

    INTEGER(isp), DIMENSION(:), POINTER, SAVE ::  iran  !<
    INTEGER(isp), DIMENSION(:), POINTER, SAVE ::  jran  !<
    INTEGER(isp), DIMENSION(:), POINTER, SAVE ::  kran  !<
    INTEGER(isp), DIMENSION(:), POINTER, SAVE ::  mran  !<
    INTEGER(isp), DIMENSION(:), POINTER, SAVE ::  nran  !<
    INTEGER(isp), DIMENSION(:), POINTER, SAVE ::  ranv  !<

    INTEGER(isp), DIMENSION(:,:), ALLOCATABLE   ::  id_random_array   !<
    INTEGER(isp), DIMENSION(:,:,:), ALLOCATABLE ::  seq_random_array  !<

    REAL(wp), SAVE ::  amm  !<

    REAL(wp) ::  random_dummy = 0.0  !<

    INTERFACE init_parallel_random_generator
       MODULE PROCEDURE init_parallel_random_generator_1d
       MODULE PROCEDURE init_parallel_random_generator_2d
    END INTERFACE

    INTERFACE random_number_parallel
       MODULE PROCEDURE ran0_s
    END INTERFACE

    INTERFACE random_number_parallel_gauss
       MODULE PROCEDURE gasdev_s
    END INTERFACE

    INTERFACE random_seed_parallel
       MODULE PROCEDURE random_seed_parallel
    END INTERFACE

    INTERFACE ran_hash
       MODULE PROCEDURE ran_hash_v
    END INTERFACE

    INTERFACE reallocate
       MODULE PROCEDURE reallocate_iv,reallocate_im
    END INTERFACE

    INTERFACE arth
       MODULE PROCEDURE arth_i
    END INTERFACE

    PRIVATE

    PUBLIC random_number_parallel,                                                                 &
           random_seed_parallel,                                                                   &
           random_dummy,                                                                           &
           id_random_array,                                                                        &
           seq_random_array,                                                                       &
           init_parallel_random_generator,                                                         &
           random_number_parallel_gauss

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize the parallel random number generator for a 1-dimensional array.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_parallel_random_generator_1d( nxy, ns, ne, id_rand, seq_rand )

    USE control_parameters,                                                                        &
        ONLY:  ensemble_member_nr

    INTEGER(isp), INTENT(IN) ::  nxy  !< constant scaling with grid dimensions
    INTEGER(isp), INTENT(IN) ::  ns   !< start index on subdomain
    INTEGER(isp), INTENT(IN) ::  ne   !< end index on subdomain

    INTEGER(iwp) ::  j  !< loop index

    INTEGER(isp), DIMENSION(ns:ne)   ::  id_rand   !< initial IDs
    INTEGER(isp), DIMENSION(5,ns:ne) ::  seq_rand  !< initial random seeds


!
!-- Asigning an ID to every vertical gridpoint column dependig on the ensemble run number.
    DO  j = ns, ne
       id_rand(j) = j * ( nxy + 1.0_wp ) + 1.0_wp + 1E6 * ( ensemble_member_nr - 1000 )
    ENDDO
!
!-- Initializing with random_seed_parallel for every vertical gridpoint column.
    random_dummy = 0
    DO  j = ns, ne
       CALL random_seed_parallel( random_sequence=id_rand(j) )
       CALL random_number_parallel( random_dummy )
       CALL random_seed_parallel( get=seq_rand(:,j) )
    ENDDO

 END SUBROUTINE init_parallel_random_generator_1d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize the parallel random number generator for a specific subdomain.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE init_parallel_random_generator_2d( nx_l, nys_l, nyn_l, nxl_l, nxr_l )

    USE kinds

    USE control_parameters,                                                                        &
        ONLY:  ensemble_member_nr

    IMPLICIT NONE

    INTEGER(isp), INTENT(IN) ::  nx_l   !< constant
    INTEGER(isp), INTENT(IN) ::  nys_l  !< local lower subdomain bound index in y-direction
    INTEGER(isp), INTENT(IN) ::  nyn_l  !< local upper subdomain bound index in y-direction
    INTEGER(isp), INTENT(IN) ::  nxl_l  !< local lower subdomain bound index in x-direction
    INTEGER(isp), INTENT(IN) ::  nxr_l  !< local upper subdomain bound index in x-direction

    INTEGER(iwp) ::  i  !< grid index x-direction
    INTEGER(iwp) ::  j  !< grid index y-direction


!
!-- Allocate ID-array and state-space-array
    ALLOCATE ( seq_random_array(5,nys_l:nyn_l,nxl_l:nxr_l) )
    ALLOCATE ( id_random_array(nys_l:nyn_l,nxl_l:nxr_l) )
    seq_random_array = 0
    id_random_array  = 0

!-- Asigning an ID to every vertical gridpoint column dependig on the ensemble run number.
    DO  i = nxl_l, nxr_l
       DO  j = nys_l, nyn_l
          id_random_array(j,i) = j * ( nx_l + 1.0_wp ) + i + 1.0_wp + 1E6 *                        &
                                 ( ensemble_member_nr - 1000 )
       ENDDO
    ENDDO
!-- Initializing with random_seed_parallel for every vertical gridpoint column.
    random_dummy = 0
    DO  i = nxl_l, nxr_l
       DO  j = nys_l, nyn_l
          CALL random_seed_parallel( random_sequence = id_random_array(j,i) )
          CALL random_number_parallel( random_dummy )
          CALL random_seed_parallel( get = seq_random_array(:,j,i) )
       ENDDO
    ENDDO

 END SUBROUTINE init_parallel_random_generator_2d

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Lagged Fibonacci generator combined with a Marsaglia shift sequence. Returns as harvest a uniform
!> random deviate between 0.0 and 1.0 (exclusive of the end point values). This generator has the
!> same calling and initialization conventions as Fortran 90's random_number routine.
!> Use random_seed_parallel to initialize or reinitialize to a particular sequence.
!> The period of this generator is about 2.0 x 10^28, and it fully vectorizes.
!> Validity of the integer model assumed by this generator is tested at initialization.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ran0_s( harvest )

    USE kinds

    IMPLICIT NONE

    REAL(wp), INTENT(OUT) ::  harvest  !<
!
!-- Initialization routine in ran_state.
    IF ( lenran < 1 )  CALL ran_init( 1 )
!
!-- Update Fibonacci generator, which has period p^2 + p + 1, p = 2^31 - 69.
    rans = iran0 - kran0

    IF ( rans < 0 )  rans = rans + 2147483579_isp

    iran0 = jran0
    jran0 = kran0
    kran0 = rans
!
!-- Update Marsaglia shift sequence with period 2^32 - 1.
    nran0 = IEOR( nran0, ISHFT( nran0, 13 ) )
    nran0 = IEOR( nran0, ISHFT( nran0, -17 ) )
    nran0 = IEOR( nran0, ISHFT( nran0, 5 ) )
!
!-- Combine the generators.
    rans  = IEOR( nran0, rans )
!
!-- Make the result positive definite (note that amm is negative).
    harvest = amm * MERGE( rans, NOT( rans ), rans < 0 )

 END SUBROUTINE ran0_s

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Returns in harvest a normally distributed deviate with zero mean and unit variance, using ran0_s
!> as the source of uniform deviates. Following Numerical Recipes in Fortran90 (Press et al., 2nd
!> edition, 1996, pp 1152ff). Note, instead of ran1_s, ran0_s is used.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE gasdev_s( harvest )

    REAL(wp), INTENT(OUT) ::  harvest  !<

    REAL(wp) ::  rsq  !<
    REAL(wp) ::  v1   !<
    REAL(wp) ::  v2   !<
    REAL(wp), SAVE ::  g  !<

    LOGICAL, SAVE ::  gaus_stored = .FALSE.  !<
!
!-- We have an extra deviate handy, so return it, and unset the flag.
    IF ( gaus_stored )  THEN
       harvest = g
       gaus_stored = .FALSE.
!
!-- We don’t have an extra deviate handy, so pick two uniform numbers in the square extending
!-- from -1 to +1 in each direction
    ELSE
       DO
          CALL ran0_s( v1 )
          CALL ran0_s( v2 )
          v1 = 2.0_wp * v1 - 1.0_wp
          v2 = 2.0_wp * v2 - 1.0_wp
!
!--       See if they are in the unit circle
          rsq = v1**2 + v2**2
!
!--       Otherwise try again.
          IF ( rsq > 0.0  .AND.  rsq < 1.0 )  EXIT
       ENDDO
!
!--    Now make the Box-Muller transformation to get two normal deviates.
!--    Return one and save the other for next time. Set flag.
       rsq = SQRT( - 2.0_sp * LOG( rsq ) / rsq )
       harvest = v1 * rsq
       g = v2 * rsq
       gaus_stored = .TRUE.
    ENDIF

 END SUBROUTINE gasdev_s

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize or reinitialize the random generator state space to vectors of size length.
!> The saved variable seq is hashed (via calls to the module routine ran_hash) to create unique
!> starting seeds, different for each vector component.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ran_init( length )

    USE kinds

    IMPLICIT NONE

    INTEGER(isp), INTENT(IN) ::  length  !<

    INTEGER(isp) ::  hg    !<
    INTEGER(isp) ::  hgm   !<
    INTEGER(isp) ::  hgng  !<

    INTEGER(isp) ::  new   !<
    INTEGER(isp) ::  j     !<
    INTEGER(isp) ::  hgt   !<
!
!-- Simply return if enough space is already allocated.
    IF ( length < lenran )  RETURN

    hg = HUGE( 1_isp )
    hgm = - hg
    hgng = hgm - 1
    hgt = hg

!-- The following lines check that kind value isp is in fact a 32-bit integer with the usual
!-- properties that we expect it to have (under negation and wrap-around addition).
!-- If all of these tests are satisfied, then the routines that use this module are portable, even
!-- though they go beyond Fortran 90's integer model.

    IF ( hg /= 2147483647   )  CALL ran_error( 'arithmetic assumption 1 failed' )
    IF ( hgng >= 0          )  CALL ran_error( 'arithmetic assumption 2 failed' )
    IF ( hgt + 1 /= hgng    )  CALL ran_error( 'arithmetic assumption 3 failed' )
    IF ( NOT( hg ) >= 0     )  CALL ran_error( 'arithmetic assumption 4 failed' )
    IF ( NOT( hgng ) < 0    )  CALL ran_error( 'arithmetic assumption 5 failed' )
    IF ( hg + hgng >= 0     )  CALL ran_error( 'arithmetic assumption 6 failed' )
    IF ( NOT( - 1_isp ) < 0 )  CALL ran_error( 'arithmetic assumption 7 failed' )
    IF ( NOT( 0_isp ) >= 0  )  CALL ran_error( 'arithmetic assumption 8 failed' )
    IF ( NOT( 1_isp ) >= 0  )  CALL ran_error( 'arithmetic assumption 9 failed' )

    IF ( lenran > 0 )  THEN                          ! Reallocate space, or ...

       ranseeds => reallocate( ranseeds, length, 5 )
       ranv => reallocate( ranv, length - 1 )
       new = lenran + 1

    ELSE                                            ! Allocate space.

       ALLOCATE( ranseeds(length,5) )
       ALLOCATE( ranv(length-1) )
       new = 1                               !- Index of first location not yet initialized.
       amm = NEAREST( 1.0_wp, - 1.0_wp ) / hgng
!
!--    Use of nearest is to ensure that returned random deviates are strictly less than 1.0.
       IF  ( amm * hgng >= 1.0  .OR.  amm * hgng <= 0.0 )                                          &
          CALL ran_error( 'arithmetic assumption 10 failed' )

    END IF
!
!-- Set starting values, unique by seq and vector component.
    ranseeds(new:,1) = seq
    ranseeds(new:,2:5) = SPREAD( arth( new, 1, SIZE( ranseeds(new:,1) ) ), 2, 4 )

    DO  j = 1, 4   ! Hash them.
       CALL ran_hash( ranseeds(new:,j), ranseeds(new:,j+1) )
    END DO

    WHERE ( ranseeds(new:,1:3) < 0 )                                                              &
       ranseeds(new:,1:3) = NOT( ranseeds(new:,1:3) )  ! Enforce nonnegativity.

    WHERE ( ranseeds(new:,4:5) == 0 ) ranseeds(new:,4:5) = 1 ! Enforce nonzero.

    IF ( new == 1 )  THEN ! Set scalar seeds.

       iran0 = ranseeds(1,1)
       jran0 = ranseeds(1,2)
       kran0 = ranseeds(1,3)
       mran0 = ranseeds(1,4)
       nran0 = ranseeds(1,5)
       rans  = nran0

    END IF

    IF ( length > 1 )  THEN   ! Point to vector seeds.

       iran => ranseeds(2:,1)
       jran => ranseeds(2:,2)
       kran => ranseeds(2:,3)
       mran => ranseeds(2:,4)
       nran => ranseeds(2:,5)
       ranv = nran

    END IF

    lenran = length

 END SUBROUTINE ran_init

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> User interface to release the workspace used by the random number routines.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ran_deallocate

    IF ( lenran > 0 )  THEN

       DEALLOCATE( ranseeds, ranv )
       NULLIFY( ranseeds, ranv, iran, jran, kran, mran, nran )
       lenran = 0

    END IF

 END SUBROUTINE ran_deallocate

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> User interface for seeding the random number routines.
!> Syntax is exactly like Fortran 90's random_seed routine, with one additional argument keyword:
!> random_sequence, set to any integer value, causes an immediate new initialization, seeded by that
!> integer.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE random_seed_parallel( random_sequence, state_size, put, get )

    IMPLICIT NONE

    INTEGER(isp), OPTIONAL, INTENT(IN)  ::  random_sequence  !<
    INTEGER(isp), OPTIONAL, INTENT(OUT) ::  state_size       !<

    INTEGER(isp), DIMENSION(:), OPTIONAL, INTENT(IN)  ::  put  !<
    INTEGER(isp), DIMENSION(:), OPTIONAL, INTENT(OUT) ::  get  !<

    IF ( PRESENT( state_size ) )  THEN

       state_size = 5 * lenran

    ELSE IF ( PRESENT( put ) )  THEN

       IF ( lenran == 0 )  RETURN

       ranseeds = RESHAPE( put, SHAPE( ranseeds ) )
!
!--    Enforce nonnegativity and nonzero conditions on any user-supplied seeds.
       WHERE ( ranseeds(:,1:3) < 0 ) ranseeds(:,1:3) = NOT( ranseeds(:,1:3) )

       WHERE ( ranseeds(:,4:5) == 0 ) ranseeds(:,4:5) = 1

       iran0 = ranseeds(1,1)
       jran0 = ranseeds(1,2)
       kran0 = ranseeds(1,3)
       mran0 = ranseeds(1,4)
       nran0 = ranseeds(1,5)

    ELSE IF ( PRESENT( get ) )  THEN

       IF ( lenran == 0 )  RETURN

       ranseeds(1,1:5) = (/ iran0, jran0, kran0, mran0, nran0 /)
       get = RESHAPE( ranseeds, SHAPE( get ) )

    ELSE IF ( PRESENT( random_sequence ) ) THEN

       CALL ran_deallocate
       seq = random_sequence

    END IF

 END SUBROUTINE random_seed_parallel

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> DES-like hashing of two 32-bit integers, using shifts, xor's, and adds to make the internal
!> nonlinear function.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ran_hash_v( il, ir )

    IMPLICIT NONE

    INTEGER(isp), DIMENSION(:), INTENT(INOUT) ::  il  !<
    INTEGER(isp), DIMENSION(:), INTENT(INOUT) ::  ir  !<

    INTEGER(isp), DIMENSION(size(il)) ::  is  !<

    INTEGER(isp) :: j  !<

    DO  j = 1, 4

       is = ir
       ir = IEOR( ir, ISHFT( ir, 5 ) ) + 1422217823
       ir = IEOR( ir, ISHFT( ir, - 16 ) ) + 1842055030
       ir = IEOR( ir, ISHFT( ir, 9 ) ) + 80567781
       ir = IEOR( il, ir )
       il = is
    END DO

 END SUBROUTINE ran_hash_v


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> User interface to process error-messages produced by the random_number_parallel module.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE ran_error( string )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    CHARACTER(LEN=*), INTENT(IN) ::  string  !< Error message string

    message_string = 'incompatible integer arithmetic: "' // TRIM( string ) // '"'
    CALL message( 'ran_init', 'PAC0266', 1, 2, 0, 6, 0 )

 END SUBROUTINE ran_error


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reallocates the generators state space "ranseeds" to vectors of size length.
!--------------------------------------------------------------------------------------------------!
 FUNCTION reallocate_iv( p, n )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    INTEGER(isp), DIMENSION(:), POINTER ::  p              !<
    INTEGER(isp), DIMENSION(:), POINTER ::  reallocate_iv  !<

    INTEGER(isp), INTENT(IN) ::  n  !<

    INTEGER(isp) ::  nold  !<
    INTEGER(isp) ::  ierr  !<


    ALLOCATE( reallocate_iv(n), stat = ierr )

    IF ( ierr /= 0 )  THEN
       message_string = 'problem in attempt to allocate memory'
       CALL message( 'reallocate_iv', 'PAC0267', 1, 2, 0, 6, 0 )
    END IF

    IF ( .NOT. ASSOCIATED( p ) )  RETURN

    nold = SIZE( p )

    reallocate_iv(1:MIN( nold, n )) = p(1:MIN( nold, n ))

    DEALLOCATE( p )

 END FUNCTION reallocate_iv


 FUNCTION reallocate_im( p, n, m )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    INTEGER(isp), DIMENSION(:,:), POINTER ::  p              !<
    INTEGER(isp), DIMENSION(:,:), POINTER ::  reallocate_im  !<

    INTEGER(isp), INTENT(IN) ::  m  !<
    INTEGER(isp), INTENT(IN) ::  n  !<

    INTEGER(isp) ::  mold  !<
    INTEGER(isp) ::  nold  !<
    INTEGER(isp) ::  ierr  !<


    ALLOCATE( reallocate_im(n,m), stat = ierr )

    IF ( ierr /= 0 )  THEN
       message_string = 'problem in attempt to allocate memory'
       CALL message( 'reallocate_im', 'PAC0267', 1, 2, 0, 6, 0 )
    END IF

    IF ( .NOT. ASSOCIATED( p ) )  RETURN

    nold = SIZE( p, 1 )
    mold = SIZE( p, 2 )

    reallocate_im(1:MIN( nold, n ),1:MIN( mold, m )) = p(1:MIN( nold, n ),1:MIN( mold, m ))

    DEALLOCATE(p)

 END FUNCTION reallocate_im

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Reallocates the generators state space "ranseeds" to vectors of size length.
!--------------------------------------------------------------------------------------------------!
 FUNCTION arth_i( first, increment, n )

    INTEGER(isp), INTENT(IN) ::  first      !<
    INTEGER(isp), INTENT(IN) ::  increment  !<
    INTEGER(isp), INTENT(IN) ::  n          !<

    INTEGER(isp), DIMENSION(n) ::  arth_i  !<

    INTEGER(isp) ::  k     !<
    INTEGER(isp) ::  k2    !<
    INTEGER(isp) ::  temp  !<

    INTEGER(isp), PARAMETER ::  npar_arth = 16  !<
    INTEGER(isp), PARAMETER ::  npar2_arth = 8  !<

    IF ( n > 0 )  arth_i(1) = first

    IF ( n <= npar_arth )  THEN

       DO  k = 2, n
          arth_i(k) = arth_i(k-1) + increment
       END DO

    ELSE

       DO  k = 2, npar2_arth
          arth_i(k) = arth_i(k-1) + increment
       END DO

       temp = increment * npar2_arth
       k = npar2_arth

       DO
          IF ( k >= n )  EXIT
          k2 = k + k
          arth_i(k+1:MIN( k2, n )) = temp + arth_i(1:MIN( k, n - k ))
          temp = temp + temp
          k = k2
       END DO

    END IF

 END FUNCTION arth_i

END MODULE random_generator_parallel
