! NOTE: This routine is a vectorised version of Rosenbrock
!       and can only be used with KP4.
!
!Current revisions:
!------------------
!
!
! Former revisions:
! -----------------------
! $Id$
! commented USE kp4_compress                       (30.10.2018, forkel)
!
! initial version (Rev. 3185)                      (June 2018, ketelsen)
!

SUBROUTINE Rosenbrock(N,Y,Tstart,Tend, &
           AbsTol,RelTol,              &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!    Solves the system y'=F(t,y) using a Rosenbrock method defined by:
!
!     G = 1/(H*gamma(1)) - Jac(t0,Y0)
!     T_i = t0 + Alpha(i)*H
!     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j
!     G * K_i = Fun( T_i, Y_i ) + \sum_{j=1}^S C(i,j)/H * K_j +
!         gamma(i)*dF/dT(t0, Y0)
!     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j
!
!    For details on Rosenbrock methods and their implementation consult:
!      E. Hairer and G. Wanner
!      "Solving ODEs II. Stiff and differential-algebraic problems".
!      Springer series in computational mathematics, Springer-Verlag, 1996.
!    The codes contained in the book inspired this implementation.
!
!    (C)  Adrian Sandu, August 2004
!    Virginia Polytechnic Institute and State University
!    Contact: sandu@cs.vt.edu
!    Revised by Philipp Miehe and Adrian Sandu, May 2006                  
!    This implementation is part of KPP - the Kinetic PreProcessor
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>   INPUT ARGUMENTS:
!
!-     Y(N)    = vector of initial conditions (at T=Tstart)
!-    [Tstart,Tend]  = time range of integration
!     (if Tstart>Tend the integration is performed backwards in time)
!-    RelTol, AbsTol = user precribed accuracy
!- SUBROUTINE  Fun( T, Y, Ydot ) = ODE function,
!                       returns Ydot = Y' = F(T,Y)
!- SUBROUTINE  Jac( T, Y, Jcb ) = Jacobian of the ODE function,
!                       returns Jcb = dFun/dY
!-    ICNTRL(1:20)    = integer inputs parameters
!-    RCNTRL(1:20)    = real inputs parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     OUTPUT ARGUMENTS:
!
!-    Y(N)    -> vector of final states (at T->Tend)
!-    ISTATUS(1:20)   -> integer output parameters
!-    RSTATUS(1:20)   -> real output parameters
!-    IERR            -> job status upon return
!                        success (positive value) or
!                        failure (negative value)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~>     INPUT PARAMETERS:
!
!    Note: For input parameters equal to zero the default values of the
!       corresponding variables are used.
!
!    ICNTRL(1) = 1: F = F(y)   Independent of T (AUTONOMOUS)
!              = 0: F = F(t,y) Depends on T (NON-AUTONOMOUS)
!
!    ICNTRL(2) = 0: AbsTol, RelTol are N-dimensional vectors
!              = 1: AbsTol, RelTol are scalars
!
!    ICNTRL(3)  -> selection of a particular Rosenbrock method
!        = 0 :    Rodas3 (default)
!        = 1 :    Ros2
!        = 2 :    Ros3
!        = 3 :    Ros4F
!        = 4 :    Rodas3
!        = 5 :    Rodas4
!
!    ICNTRL(4)  -> maximum number of integration steps
!        For ICNTRL(4)=0) the default value of 100000 is used
!
!    RCNTRL(1)  -> Hmin, lower bound for the integration step size
!          It is strongly recommended to keep Hmin = ZERO
!    RCNTRL(2)  -> Hmax, upper bound for the integration step size
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!
!    RCNTRL(4)  -> FacMin, lower bound on step decrease factor (default=0.2)
!    RCNTRL(5)  -> FacMax, upper bound on step increase factor (default=6)
!    RCNTRL(6)  -> FacRej, step decrease factor after multiple rejections
!                          (default=0.1)
!    RCNTRL(7)  -> FacSafe, by which the new step is slightly smaller
!         than the predicted value  (default=0.9)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!    OUTPUT ARGUMENTS:
!    -----------------
!
!    T           -> T value for which the solution has been computed
!                     (after successful return T=Tend).
!
!    Y(N)        -> Numerical solution at T
!
!    IDID        -> Reports on successfulness upon return:
!                    = 1 for success
!                    < 0 for error (value equals error code)
!
!    ISTATUS(1)  -> No. of function calls
!    ISTATUS(2)  -> No. of jacobian calls
!    ISTATUS(3)  -> No. of steps
!    ISTATUS(4)  -> No. of accepted steps
!    ISTATUS(5)  -> No. of rejected steps (except at very beginning)
!    ISTATUS(6)  -> No. of LU decompositions
!    ISTATUS(7)  -> No. of forward/backward substitutions
!    ISTATUS(8)  -> No. of singular matrix decompositions
!
!    RSTATUS(1)  -> Texit, the time corresponding to the
!                     computed Y upon return
!    RSTATUS(2)  -> Hexit, last accepted step before exit
!    RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
!                   For multiple restarts, use Hnew as Hstart 
!                     in the subsequent run
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  IMPLICIT NONE

!~~~>  Arguments
   INTEGER,       INTENT(IN)    :: N
   REAL(kind=dp), INTENT(INOUT) :: Y(:,:)
   REAL(kind=dp), INTENT(IN)    :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)    :: AbsTol(N),RelTol(N)
   INTEGER,       INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)    :: RCNTRL(20)
   INTEGER,       INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR
!~~~>  Parameters of the Rosenbrock method, up to 6 stages
   INTEGER ::  ros_S, rosMethod
   INTEGER, PARAMETER :: RS2=1, RS3=2, RS4=3, RD3=4, RD4=5, RG3=6
   REAL(kind=dp) :: ros_A(15), ros_C(15), ros_M(6), ros_E(6), &
                    ros_Alpha(6), ros_Gamma(6), ros_ELO
   LOGICAL :: ros_NewF(6)
   CHARACTER(LEN=12) :: ros_Name
!~~~>  Local variables
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart
   REAL(kind=dp) :: Texit
   INTEGER       :: i, UplimTol, Max_no_steps
   LOGICAL       :: Autonomous, VectorTol
!~~~>   Parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
   INTEGER        :: solver_method


!~~~>  Autonomous or time dependent ODE. Default is time dependent.
   Autonomous = .NOT.(ICNTRL(1) == 0)

!~~~>  For Scalar tolerances (ICNTRL(2).NE.0)  the code uses AbsTol(1) and RelTol(1)
!   For Vector tolerances (ICNTRL(2) == 0) the code uses AbsTol(1:N) and RelTol(1:N)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
      UplimTol  = N
   ELSE
      VectorTol = .FALSE.
      UplimTol  = 1
   END IF

   solver_method = ICNTRL(3)

!~~~>   Initialize the particular Rosenbrock method selected
   SELECT CASE (solver_method)
     CASE (1)
       CALL Ros2
     CASE (2)
       CALL Ros3
     CASE (3)
       CALL Ros4
     CASE (0,4)
       CALL Rodas3
     CASE (5)
       CALL Rodas4
     CASE (6)
       CALL Rang3
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(3)=',ICNTRL(3) 
       CALL ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       IERRV(:) = -2
       RETURN
   END SELECT

!~~~>   The maximum number of steps admitted
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 200000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      IERRV(:) = -1
      RETURN
   END IF

!~~~>  Unit roundoff (1+Roundoff>1)
   Roundoff = WLAMCH('E')

!~~~>  Lower bound on the step size: (positive value)
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      IERRV(:) = -3
      RETURN
   END IF
!~~~>  Upper bound on the step size: (positive value)
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      IERRV(:) = -3
      RETURN
   END IF
!~~~>  Starting step size: (positive value)
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
!~~~>  Step size can be changed s.t.  FacMin < Hnew/Hold < FacMax
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      IERRV(:) = -4
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      IERRV(:) = -4
      RETURN
   END IF
!~~~>   FacRej: Factor to decrease step after 2 succesive rejections
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      IERRV(:) = -4
      RETURN
   END IF
!~~~>   FacSafe: Safety Factor in the computation of new step size
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      IERRV(:) = -4
      RETURN
   END IF
!~~~>  Check if tolerances are reasonable
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO


!~~~>  CALL Rosenbrock method
   CALL ros_Integrator(Y, Tstart, Tend, Texit,   &
        AbsTol, RelTol,                          &
!  Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!  Error indicator
        IERR)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS !  SUBROUTINES internal to Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 SUBROUTINE ros_ErrorMsg(Code,T,H,IERR)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
!    Handles all error messages
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
  
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR
   
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:' 
     
   SELECT CASE (Code)
    CASE (-1)    
      PRINT * , '--> Improper value for maximal no of steps'
    CASE (-2)    
      PRINT * , '--> Selected Rosenbrock method not implemented'
    CASE (-3)    
      PRINT * , '--> Hmin/Hmax/Hstart must be positive'
    CASE (-4)    
      PRINT * , '--> FacMin/FacMax/FacRej must be positive'
    CASE (-5) 
      PRINT * , '--> Improper tolerance values'
    CASE (-6) 
      PRINT * , '--> No of steps exceeds maximum bound'
    CASE (-7) 
      PRINT * , '--> Step size too small: T + 10*H = T', &
            ' or H < Roundoff'
      PRINT *, 'T=', T, 'and H=', H
    CASE (-8)    
      PRINT * , '--> Matrix is repeatedly singular'
    CASE DEFAULT
      PRINT *, 'Unknown Error code: ', Code
   END SELECT
   
     
 END SUBROUTINE ros_ErrorMsg

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE ros_Integrator (Y, Tstart, Tend, Tout,  &
        AbsTol, RelTol,                          &
!~~~> Integration parameters
        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart,            &
        FacMin, FacMax, FacRej, FacSafe,         &
!~~~> Error indicator
        IERR_out )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Template for the implementation of a generic Rosenbrock method
!      defined by ros_S (no of stages)
!      and its coefficients ros_{A,C,M,E,Alpha,Gamma}
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! USE kp4_compress,    ONLY: kco_initialize, kco_compress, &
!                            kco_finalize, cell_done
  IMPLICIT NONE

!~~~> Input: the initial condition at Tstart; Output: the solution at T
   REAL(kind=dp), INTENT(INOUT) :: Y(:,:)
!~~~> Input: integration interval
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
!~~~> Output: time at which the solution is returned (T=Tend if success)
   REAL(kind=dp), INTENT(OUT) ::  Tout
!~~~> Input: tolerances
   REAL(kind=dp), INTENT(IN) ::  AbsTol(N), RelTol(N)
!~~~> Input: integration parameters
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
!~~~> Output: Error indicator
   INTEGER, INTENT(OUT) :: IERR_out
! ~~~~ Local variables
   REAL(kind=dp) :: Ynew(VL_DIM,N), Fcn0(VL_DIM,N), Fcn(VL_DIM,N)
   REAL(kind=dp) :: K(VL_DIM,N*ros_S), dFdT(VL_DIM,N)
   REAL(kind=dp) :: Jac0(VL_DIM,LU_NONZERO), Ghimj(VL_DIM,LU_NONZERO)
!kk   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Yerr(VL_DIM,N)
   REAL(kind=dp) :: Y_save(VL_DIM,N),Ynew_save(VL_DIM,N)
   REAL(kind=dp) :: Hnew_save(VL_DIM)
   INTEGER :: Pivot(N), Direction, ioffset, j, istage, i
   LOGICAL :: Singular
!~~~>  Local parameters
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp

   LOGICAL,SAVE             :: lfirst = .TRUE.
!~~~>  Locally called functions
!    REAL(kind=dp) WLAMCH
!    EXTERNAL WLAMCH
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  vector compression algorithm

   REAL(kind=dp), dimension(VL_dim)      :: T
   REAL(kind=dp), dimension(VL_dim)      :: H, Hnew, HC, HG, Fac, Tau, Err,ros_v
   LOGICAL, dimension(VL_dim)            :: RejectLastH, RejectMoreH
   LOGICAL, dimension(VL_dim)            :: Accept_step,y_copied
   integer                               :: vl_save
   integer                               :: loop_count

   if (lfirst)  then
      write(9,'(a,i4,/)') 'kpp4palm vector version   ==> Vector length = ',vl_dim
   end if

!~~~>  Initial preparations
   T(1:VL) = Tstart
   RSTATUS(Nhexit) = ZERO
   Tout    = Tend
   H(1:VL) = MIN( MAX(ABS(Hmin),ABS(Hstart)) , ABS(Hmax) )
   where (ABS(H(1:VL)) <= 10.0_dp*Roundoff) H(1:VL) = DeltaMin


   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   H = Direction*H

   RejectLastH(1:VL)=.FALSE.
   RejectMoreH(1:VL)=.FALSE.

   loop_count = 0
!~~~> Time loop begins below

   IERRV(1:VL) = 1  ! successful integration
   vl_save          = VL
   Kacc       = 0
   Krej       = 0

!kk   IF (.not. l_fixed_step) then
      call kco_initialize (VL, y, IERRV, Kacc, Krej)
!kk   end if

TimeLoop: DO
!kk Eventuell als Where einbauen
   where ( (Direction > 0).AND.((T(1:VL)-Tend)+Roundoff <= ZERO) &
        .OR. (Direction < 0).AND.((Tend-T(1:VL))+Roundoff <= ZERO) )
      cell_done(1:VL) = .FALSE.
   elsewhere
      cell_done(1:VL) = .TRUE.
   end where
!   do i=1,vl
!      if ( (Direction > 0).AND.((T(i)-Tend)+Roundoff <= ZERO) &
!        .OR. (Direction < 0).AND.((Tend-T(i))+Roundoff <= ZERO) )  then
!         cell_done(i) = .FALSE.
!      else
!         cell_done(i) = .TRUE.
!       end if
!   end do

   if(count(cell_done(1:VL)) == vl_save)   then                      ! EXIT 1, all cells converge at the same time
      EXIT                                                           ! Iteration terminates without without compress
   end if

   loop_count = loop_count+1

   call kco_compress (VL, T, H, Hnew, IERRV, y, RCONST, FIX, RejectLastH, RejectMoreH, Kacc, Krej)

   if(all(cell_done(1:VL)))  then                                    ! EXIT 2, all cells converge have been converged
     EXIT
   end if

   IF ( ISTATUS(Nstp) > Max_no_steps ) THEN  ! Too many steps
      CALL ros_ErrorMsg(-6,T(1),H(1),IERR)
      IERRV(1:VL) = -6
      RETURN
   END IF

   IERRV = 0
   WHERE ( ((T(1:VL)+0.1_dp*H(1:VL)) == T(1)).OR.(H(1:VL) <= Roundoff) )   ! Step size too small
      IERRV(1:VL) = -7
   END WHERE

   IF(MINVAL(IERRV) == -7)   then
      CALL ros_ErrorMsg(-7,T(1),H(1),IERR)
      RETURN
    END IF

!~~~>  Limit H if necessary to avoid going beyond Tend
   H(1:VL) = MIN(H(1:VL),ABS(Tend-T(1:VL)))

!~~~>   Compute the function at current time
   CALL FunTemplate(T,Y,Fcn0)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1

!~~~>  Compute the function derivative with respect to T
   IF (.NOT.Autonomous) THEN
      CALL ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
   END IF

!~~~>   Compute the Jacobian at current time
   CALL JacTemplate(T,Y,Jac0)
   ISTATUS(Njac) = ISTATUS(Njac) + 1                            ! This is not correct in vector mode, please use xnacc and xnrej

!~~~>  Repeat step calculation until current step accepted
    Accept_step(1:VL) = .FALSE.
    y_copied(1:VL)    = .FALSE.

UntilAccepted: DO

   CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular)
   IF (Singular) THEN ! More than 5 consecutive failed decompositions
       CALL ros_ErrorMsg(-8,T(1),H(1),IERR)
       RETURN
   END IF

!~~~>   Compute the stages
Stage: DO istage = 1, ros_S

!      Current istage offset. Current istage vector is K(ioffset+1:ioffset+N)
       ioffset = N*(istage-1)

!      For the 1st istage the function has been computed previously
       IF ( istage == 1 ) THEN
         Fcn(1:VL,1:N) = Fcn0(1:VL,1:N) 
      ! istage>1 and a new function evaluation is needed at the current istage
       ELSEIF ( ros_NewF(istage) ) THEN
	   Ynew(1:VL,1:N) = Y(1:VL,1:N)
        DO j = 1, istage-1
!           ros_v = ros_A((istage-1)*(istage-2)/2+j)
!           CALL WAXPY(N,ros_v,K(:,N*(j-1)+1:N*(j-1)+N),1,Ynew,1)
            Ynew(1:vl,1:N) = Ynew(1:vl,1:N)+ros_A((istage-1)*(istage-2)/2+j)*K(1:VL,N*(j-1)+1:N*(j-1)+N)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL FunTemplate(Tau,Ynew,Fcn)
         ISTATUS(Nfun) = ISTATUS(Nfun) + 1
       END IF ! if istage == 1 elseif ros_NewF(istage)
       K(1:VL,ioffset+1:ioffset+N) = Fcn(1:VL,1:N)
       DO j = 1, istage-1
         HC(1:vl) = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H(1:vl))
!         CALL WAXPY(N,HC,K(:,N*(j-1)+1:N*(j-1)+N),1,K(:,ioffset+1:ioffset+n),1)
         DO i=1,N
            K(1:vl,ioffset+i) = K(1:vl,ioffset+i) + HC(1:vl)*K(1:vl,N*(j-1)+i)
         END DO

       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG(1:vl) = Direction*H(1:vl)*ros_Gamma(istage)
         CALL WAXPY(N,HG,dFdT,1,K(:,ioffset+1:ioffset+n),1)
       END IF

       CALL ros_Solve(Ghimj, K(:,ioffset+1:ioffset+n))

   END DO Stage


!~~~>  Compute the new solution
   Ynew(1:vl,1:N) = Y(1:vl,1:N)
   DO j=1,ros_S
         ros_v = ros_M(j)
         CALL WAXPY(N,ros_v,K(:,N*(j-1)+1:),1,Ynew,1)
   END DO

!~~~>  Compute the error estimation
   Yerr(1:vl,1:N) = ZERO
   DO j=1,ros_S
        ros_v = ros_E(j)
        CALL WAXPY(N,ros_v,K(:,N*(j-1)+1:N*(j-1)+N),1,Yerr,1)
   END DO
   Err(1:VL) = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

!~~~> New step size is bounded by FacMin <= Hnew/H <= FacMax
   Fac(1:vl)  = MIN(FacMax,MAX(FacMin,FacSafe/Err(1:vl)**(ONE/ros_ELO)))
   Hnew(1:vl) = H(1:vl)*Fac(1:vl)

!~~~>  Check the error magnitude and adjust step size
   ISTATUS(Nstp) = ISTATUS(Nstp) + 1
!   write(9,*) 'Hmin ',Hmin,Hmax
 do i=1,VL
    IF(Accept_step(i))   then
       Ynew(i,:) = Ynew_save(i,:)
       Hnew(i)   = Hnew_save(i)
    ELSE IF ( (Err(i) <= ONE).OR.(H(i) <= Hmin) ) THEN  !~~~> Accept step
      ISTATUS(Nacc) = ISTATUS(Nacc) + 1
      Y(i,1:N) = Ynew(i,1:N)
      T(i) = T(i) + Direction*H(i)
      Hnew(i) = MAX(Hmin,MIN(Hnew(i),Hmax))
      IF (RejectLastH(i)) THEN  ! No step size increase after a rejected step
         Hnew(i) = MIN(Hnew(i),H(i))
      END IF
      RSTATUS(Nhexit) = H(i)                    ! ueberdenken
      RSTATUS(Nhnew)  = Hnew(i)
      RSTATUS(Ntexit) = T(i)                    ! ueberdenken ende
      RejectLastH(i) = .FALSE.
      RejectMoreH(i) = .FALSE.
      H(i) = Hnew(i)

      Accept_step(i) = .TRUE.                   ! EXIT THE LOOP: WHILE STEP NOT ACCEPTED
      Kacc(i) = Kacc(i) +1
   ELSE           !~~~> Reject step
      IF (RejectMoreH(i)) THEN
         Hnew (i)= H(i)*FacRej
      END IF
      RejectMoreH(i) = RejectLastH(i)
      RejectLastH(i) = .TRUE.
      H(i) = Hnew(i)

      IF (ISTATUS(Nacc) >= 1)  ISTATUS(Nrej) = ISTATUS(Nrej) + 1
      Krej(i) = Krej(i)+1

   END IF ! Err <= 1
  END DO
   if(all(Accept_step(1:VL)))   EXIT

   do i=1,VL
      IF(Accept_step(i))   then
         if(.not. y_copied(i))  then
            y_save(i,:)    = y(i,:)
            Ynew_save(i,:) = Ynew(i,:)
            Hnew_save(i)   = Hnew(i)
            y_copied(i) = .TRUE.
         end if
      end if
   end do


   END DO UntilAccepted

    do i=1,VL
       if(y_copied(i))  then
           y(i,:)  = y_save(i,:)
       end if
    end do

   END DO TimeLoop

   call kco_finalize (y, IERRV, Kacc, Krej)

   VL = VL_save

!~~~> Succesful exit
   IERR_out = 1  !~~~> The integration was successful

   lfirst = .FALSE.

  END SUBROUTINE ros_Integrator


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  FUNCTION ros_ErrorNorm ( Y, Ynew, Yerr, &
               AbsTol, RelTol, VectorTol )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> Computes the "scaled norm" of the error vector Yerr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

! Input arguments
   REAL(kind=dp), INTENT(IN) :: Y(:,:), Ynew(:,:), &
          Yerr(:,:), AbsTol(N), RelTol(N)
   LOGICAL, INTENT(IN) ::  VectorTol
! Local variables
   REAL(kind=dp), dimension(VL)       :: Err, Scale, Ymax
   INTEGER  :: i,j
   REAL(kind=dp), PARAMETER           :: ZERO = 0.0_dp
   REAL(kind=dp), dimension(VL)       :: ros_ErrorNorm

   Err = ZERO
   DO i=1,NVAR
       do j=1,VL
           Ymax(j) = MAX( ABS(Y(j,i)),ABS(Ynew(j,i)))
           IF (VectorTol) THEN
               Scale(j) = AbsTol(i)+ RelTol(i)*Ymax(j)
           ELSE
               Scale(j) = AbsTol(1)+ RelTol(1)*Ymax(1)
           END IF
           Err(j) = Err(j)+ (Yerr(j,i)/Scale(j))**2
       end do
   END DO
   Err  = SQRT(Err/N)

   ros_ErrorNorm = max(err, 1.0d-10)

  END FUNCTION ros_ErrorNorm


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_FunTimeDerivative ( T, Roundoff, Y, &
                Fcn0, dFdT )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~> The time partial derivative of the function by finite differences
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) :: T(:), Roundoff, Y(:,:), Fcn0(:,:)
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: dFdT(:,:)
!~~~> Local variables
   REAL(kind=dp) :: Delta(1:vl_dim),one_v(1:vl_dim)
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp
   INTEGER   :: i

   Delta(1:vl) = SQRT(Roundoff)*MAX(DeltaMin,ABS(T(1:vl)))
   CALL FunTemplate(T(1:vl)+Delta(1:vl),Y,dFdT)
   ISTATUS(Nfun) = ISTATUS(Nfun) + 1
   one_v = -ONE
   CALL WAXPY(N,one_v,Fcn0,1,dFdT,1)
!   CALL WSCAL(N,(ONE/Delta),dFdT,1)
   do i=1,size(dFdT,2)
      dFdT(1:vl,i) = 1.0_dp/Delta(1:vl)*dFdT(1:vl,i)
   end do

  END SUBROUTINE ros_FunTimeDerivative


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular )
! --- --- --- --- --- --- --- --- --- --- --- --- ---
!  Prepares the LHS matrix for stage calculations
!  1.  Construct Ghimj = 1/(H*ham) - Jac0
!      "(Gamma H) Inverse Minus Jacobian"
!  2.  Repeat LU decomposition of Ghimj until successful.
!       -half the step size if LU decomposition fails and retry
!       -exit after 5 consecutive fails
! --- --- --- --- --- --- --- --- --- --- --- --- ---
   IMPLICIT NONE

!~~~> Input arguments
   REAL(kind=dp), INTENT(IN) ::  Jac0(:,:)
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction
!~~~> Output arguments
   REAL(kind=dp), INTENT(OUT) :: Ghimj(:,:)

   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(N)
!~~~> Inout arguments
   REAL(kind=dp), INTENT(INOUT) :: H(:)   ! step size is decreased when LU fails
!~~~> Local variables
   INTEGER  :: i, ISING, Nconsecutive
   REAL(kind=dp) :: ghinv(1:VL)

   Nconsecutive = 0
   Singular = .FALSE.

!kk   DO WHILE (Singular)                        !test of Singular currently disabled

!~~~>    Construct Ghimj = 1/(H*gam) - Jac0
   Ghimj(1:VL,1:LU_NONZERO) = -JAC0(1:VL,1:LU_NONZERO)
   ghinv(1:VL) = 1.0_dp/(H(1:VL)*gam)
   DO i=1,N
     Ghimj(1:VL,LU_DIAG(i)) = Ghimj(1:VL,LU_DIAG(i))+ghinv(1:VL)
   END DO

!~~~>    Compute LU decomposition
   CALL KppDecomp( Ghimj, ising )

!kk Original code, implement singular check if necessary
!!~~~>    Compute LU decomposition
!     CALL ros_Decomp( Ghimj, Pivot, ISING )
!     IF (ISING == 0) THEN
!!~~~>    If successful done
!        Singular = .FALSE.
!     ELSE ! ISING .ne. 0
!!~~~>    If unsuccessful half the step size; if 5 consecutive fails then return
!        ISTATUS(Nsng) = ISTATUS(Nsng) + 1
!        Nconsecutive = Nconsecutive+1
!        Singular = .TRUE.
!        PRINT*,'Warning: LU Decomposition returned ISING = ',ISING
!        IF (Nconsecutive <= 5) THEN ! Less than 5 consecutive failed decompositions
!           H = H*HALF
!        ELSE  ! More than 5 consecutive failed decompositions
!           RETURN
!        END IF  ! Nconsecutive
!      END IF    ! ISING

!    END DO ! WHILE Singular

     Pivot = 0                                !Set to avoid compiler warnings

  END SUBROUTINE ros_PrepareMatrix


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  SUBROUTINE ros_Decomp( A, Pivot, ISING )
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!  Template for the LU decomposition
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   IMPLICIT NONE
!!~~~> Inout variables
!#ifdef FULL_ALGEBRA
!   REAL(kind=dp), INTENT(INOUT) :: A(N,N)
!#else
!   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)
!#endif
!!~~~> Output variables
!   INTEGER, INTENT(OUT) :: Pivot(N), ISING
!
!#ifdef FULL_ALGEBRA
!   CALL  DGETRF( N, N, A, N, Pivot, ISING )
!#else
!   CALL KppDecomp ( A, ISING )
!   Pivot(1) = 1
!#endif
!   ISTATUS(Ndec) = ISTATUS(Ndec) + 1
!
!  END SUBROUTINE ros_Decomp


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE ros_Solve( A, b )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Template for the forward/backward substitution (using pre-computed LU decomposition)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   IMPLICIT NONE
!~~~> Input variables
   REAL(kind=dp), INTENT(IN)         :: A(:,:)
!~~~> InOut variables
   REAL(kind=dp), INTENT(INOUT)      :: b(:,:)

   CALL KppSolve( A, b )

   ISTATUS(Nsol) = ISTATUS(Nsol) + VL

  END SUBROUTINE ros_Solve



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 2 stages, order 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   DOUBLE PRECISION g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    rosMethod = RS2
!~~~> Name of the method
    ros_Name = 'ROS-2'
!~~~> Number of stages
    ros_S = 2

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
!~~~> M_i = Coefficients for new step solution
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
! E_i = Coefficients for error estimator
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus one
    ros_ELO = 2.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE Ros2


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- AN L-STABLE METHOD, 3 stages, order 3, 2 function evaluations
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE
   rosMethod = RS3
!~~~> Name of the method
   ros_Name = 'ROS-3'
!~~~> Number of stages
   ros_S = 3

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514_dp
! E_i = Coefficients for error estimator
   ros_E(1) =  0.5_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199_dp
!~~~> ros_ELO = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1)= 0.0_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1)= 0.43586652150845899941601945119356_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE Ros3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Ros4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     L-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 4 STAGES
!     L-STABLE EMBEDDED ROSENBROCK METHOD OF ORDER 3
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1990)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RS4
!~~~> Name of the method
   ros_Name = 'ROS-4'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626_dp
   ros_C(6) =-0.6949742501781779_dp
!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792_dp
   ros_M(3) = 0.4353179431840180_dp
   ros_M(4) = 0.1093502252409163E+01_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) =-0.2815431932141155_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311_dp
   ros_E(4) =-0.1093502252409163E+01_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 4.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900_dp
   ros_Alpha(4) = ros_Alpha(3)
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5728200000000000_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482_dp
   ros_Gamma(4) =-0.1049021087100450_dp

  END SUBROUTINE Ros4

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! --- A STIFFLY-STABLE METHOD, 4 stages, order 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

   rosMethod = RD3
!~~~> Name of the method
   ros_Name = 'RODAS-3'
!~~~> Number of stages
   ros_S = 4

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:
!       A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!       C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

   ros_A(1) = 0.0_dp
   ros_A(2) = 2.0_dp
   ros_A(3) = 0.0_dp
   ros_A(4) = 2.0_dp
   ros_A(5) = 0.0_dp
   ros_A(6) = 1.0_dp

   ros_C(1) = 4.0_dp
   ros_C(2) = 1.0_dp
   ros_C(3) =-1.0_dp
   ros_C(4) = 1.0_dp
   ros_C(5) =-1.0_dp
   ros_C(6) =-(8.0_dp/3.0_dp)

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.
!~~~> M_i = Coefficients for new step solution
   ros_M(1) = 2.0_dp
   ros_M(2) = 0.0_dp
   ros_M(3) = 1.0_dp
   ros_M(4) = 1.0_dp
!~~~> E_i  = Coefficients for error estimator
   ros_E(1) = 0.0_dp
   ros_E(2) = 0.0_dp
   ros_E(3) = 0.0_dp
   ros_E(4) = 1.0_dp
!~~~> ros_ELO  = estimator of local order - the minimum between the
!    main and the embedded scheme orders plus 1
   ros_ELO  = 3.0_dp
!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.0_dp
   ros_Alpha(3) = 1.0_dp
   ros_Alpha(4) = 1.0_dp
!~~~> Gamma_i = \sum_j  gamma_{i,j}
   ros_Gamma(1) = 0.5_dp
   ros_Gamma(2) = 1.5_dp
   ros_Gamma(3) = 0.0_dp
   ros_Gamma(4) = 0.0_dp

  END SUBROUTINE Rodas3

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     STIFFLY-STABLE ROSENBROCK METHOD OF ORDER 4, WITH 6 STAGES
!
!      E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL
!      EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS.
!      SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS,
!      SPRINGER-VERLAG (1996)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RD4
!~~~> Name of the method
    ros_Name = 'RODAS-4'
!~~~> Number of stages
    ros_S = 6

!~~~> Y_stage_i ~ Y( T + H*Alpha_i )
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp

!~~~> Gamma_i = \sum_j  gamma_{i,j}
    ros_Gamma(1) = 0.2500000000000000_dp
    ros_Gamma(2) =-0.1043000000000000_dp
    ros_Gamma(3) = 0.1035000000000000_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp

!~~~> The coefficient matrices A and C are strictly lower triangular.
!   The lower triangular (subdiagonal) elements are stored in row-wise order:
!   A(2,1) = ros_A(1), A(3,1)=ros_A(2), A(3,2)=ros_A(3), etc.
!   The general mapping formula is:  A(i,j) = ros_A( (i-1)*(i-2)/2 + j )
!                  C(i,j) = ros_C( (i-1)*(i-2)/2 + j )

    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826_dp
    ros_A(3) = 0.2557011698983284_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915_dp
    ros_C(4) =-0.1073529058151375_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp

!~~~> M_i = Coefficients for new step solution
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0_dp
    ros_M(6) = 1.0_dp

!~~~> E_i  = Coefficients for error estimator
    ros_E(1) = 0.0_dp
    ros_E(2) = 0.0_dp
    ros_E(3) = 0.0_dp
    ros_E(4) = 0.0_dp
    ros_E(5) = 0.0_dp
    ros_E(6) = 1.0_dp

!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 4.0_dp

  END SUBROUTINE Rodas4
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! STIFFLY-STABLE W METHOD OF ORDER 3, WITH 4 STAGES
!
! J. RANG and L. ANGERMANN
! NEW ROSENBROCK W-METHODS OF ORDER 3
! FOR PARTIAL DIFFERENTIAL ALGEBRAIC
!        EQUATIONS OF INDEX 1                                             
! BIT Numerical Mathematics (2005) 45: 761-787
!  DOI: 10.1007/s10543-005-0035-y
! Table 4.1-4.2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IMPLICIT NONE

    rosMethod = RG3
!~~~> Name of the method
    ros_Name = 'RANG-3'
!~~~> Number of stages
    ros_S = 4

    ros_A(1) = 5.09052051067020d+00;
    ros_A(2) = 5.09052051067020d+00;
    ros_A(3) = 0.0d0;
    ros_A(4) = 4.97628111010787d+00;
    ros_A(5) = 2.77268164715849d-02;
    ros_A(6) = 2.29428036027904d-01;

    ros_C(1) = -1.16790812312283d+01;
    ros_C(2) = -1.64057326467367d+01;
    ros_C(3) = -2.77268164715850d-01;
    ros_C(4) = -8.38103960500476d+00;
    ros_C(5) = -8.48328409199343d-01;
    ros_C(6) =  2.87009860433106d-01;

    ros_M(1) =  5.22582761233094d+00;
    ros_M(2) = -5.56971148154165d-01;
    ros_M(3) =  3.57979469353645d-01;
    ros_M(4) =  1.72337398521064d+00;

    ros_E(1) = -5.16845212784040d+00;
    ros_E(2) = -1.26351942603842d+00;
    ros_E(3) = -1.11022302462516d-16;
    ros_E(4) =  2.22044604925031d-16;

    ros_Alpha(1) = 0.0d00;
    ros_Alpha(2) = 2.21878746765329d+00;
    ros_Alpha(3) = 2.21878746765329d+00;
    ros_Alpha(4) = 1.55392337535788d+00;

    ros_Gamma(1) =  4.35866521508459d-01;
    ros_Gamma(2) = -1.78292094614483d+00;
    ros_Gamma(3) = -2.46541900496934d+00;
    ros_Gamma(4) = -8.05529997906370d-01;


!~~~> Does the stage i require a new function evaluation (ros_NewF(i)=TRUE)
!   or does it re-use the function evaluation from stage i-1 (ros_NewF(i)=FALSE)
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.

!~~~> ros_ELO  = estimator of local order - the minimum between the
!        main and the embedded scheme orders plus 1
    ros_ELO = 3.0_dp

  END SUBROUTINE Rang3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   End of the set of internal Rosenbrock subroutines
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
END SUBROUTINE Rosenbrock
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





