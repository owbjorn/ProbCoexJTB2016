program main
!*****************************************************************************80
!! MAIN is the main program for REPLI3D, i.e.
!    the 3-species replicator system function for RKF45 library.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!    25 February 2015
!
!  Author:
!    Björn Vessman (Based on examples for rkf45 by John Burkardt)
!
  implicit none
  interface
    subroutine repRSolv( t_stop, n_runs, n_step, mdelflag )
      double precision     :: t_stop
      integer ( kind = 4 ) :: n_runs, n_step, mdelflag
    end subroutine repRSolv
    subroutine rSolv( t_stop, n_step, mdelflag )
      double precision     :: t_stop
      integer ( kind = 4 ) :: n_step, mdelflag
    end subroutine rSolv
    subroutine rFun ( t, y, yp )
      double precision :: t
      double precision :: y(:), yp(:)
    end subroutine rFun
    subroutine initRep ( mdlflag, eta, y, E )
      integer          :: mdlflag
      double precision :: eta
      double precision :: y(:)
      double precision :: E(:,:)
    end subroutine initRep
  end interface

  integer          :: neqs
  integer          :: n_step, n_reps, mdlflag
  double precision :: t_stop

  write ( *, '(a)' ) ' 3-SPECIES REPLICATOR SYSTEM'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  '

  ! Make test call
  t_stop  = 1.0d6
  n_step  = 10000
  n_reps  = 100000
  mdlflag = 3

  ! Call repeated solver
  call repRSolv( t_stop, n_reps, n_step, mdlflag )
  
  ! Call single-system solver for a phase portrait
  call rSolv( t_stop, n_step, mdlflag )

  ! Print output info
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Start and end states Y0, YT are found in file ''repli3Y ... ''.'
  write ( *, '(a)' ) '  Coordinates Gamma found in file ''repli3Gam ... ''.'

  stop
end

subroutine repRSolv( t_stop, n_runs, n_step, mdelflag )
!*****************************************************************************80
!! REPRSOLV is Repeated REPlicator SOLVer.
  implicit none
  integer ( kind = 4 ), parameter :: neqn = 3
  integer ( kind = 4 ) :: i, k, flag, i_step, n_runs, n_step, mdelflag
  external rFun 
  double precision :: y(neqn), y0(neqn), yp(neqn) ! Solver params
  double precision :: abserr, relerr, t, t_out, t_start, t_stop
  double precision :: E(neqn, 0:neqn)   ! Model params
  double precision :: gam, kap, eta
  COMMON /metabs/ eta, gam, E ! Shared with RHS function f(t,y,yp)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8_RKF45 repeated solution of 3D replicator system:'
  write ( *, '(a)' ) '  DY = Y( f(Y) - dot(Y,f(Y)) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Repeats        Tstop     Nsteps     dt '
  write ( *, '(4x,I6,2x,F10.0,2x,I8,2x,F10.5)' ), n_runs, t_stop, n_step, t_stop/dble(n_step)

  ! Check that the chosen model is available
  if ((mdelflag .GT. 4) .OR. (mdelflag .LT. 1)) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)', advance='no' ) '  Non-supported model flag,'
    write ( *, '(a)' ) ' Flag = 1,2,3,4 supported. '
    write ( *, '(a)' ) '  Aborting ... '
    STOP
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Parameter model: '
  SELECT CASE (mdelflag) ! Print model choice
    CASE (1)
      write( *, '(a)') '    Independent, uniformly distributed'
    CASE (2)
      write( *, '(a)') '    Tree model, uniformly distributed'
    CASE (3)
      write( *, '(a)') '    Independent, exponentially distributed'
    CASE (4)
      write( *, '(a)') '    Tree model, exponentially distributed'
    CASE default
  END SELECT

  ! Fixed initialisations
  gam = 0.03d0
  kap = 0.25d0
  eta = kap/(kap +gam)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  gamma   Kappa     eta'
  write ( *, '(2x,F4.2,4x,F4.2,4x,F8.6)' ), gam, kap, eta

  ! Open data files. File units used:
  !    1    System states y0(1:neqn) yT(neqn)
  !    2a)  Xi-coordinates Xi_a Xi_b        if neqn == 2
  !         as defined in Lundh&Gerlee (2013) Eqs (18), (19)
  !    2b)  Gamma-coordinates G12, G23, G31 if neqn == 3
  !         as defined JUST BELOW Lundh&Gerlee (2013) Eq (24)
  open( unit = 1,       file = './repli3Y_N1e5_expE'   )
  write ( 1, '(a)' )   &
   ' FLAG     Y0(1)     Y0(2)     Y0(3)        Y(1)              Y(2)              Y(3)   '

  open( unit = 2,   file = './repli3Gam_N1e5_expE' )
  write ( 2, '(a)' ) '             G12             G23             G31  '

  open( unit = 3,   file = './repli3LD_N1e5_expE' )
  write ( 3, '(a)' ) '    l21-l22     l32-l33     l13-l11     l11-l12     l22-l23     l33-l31  '

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Running solver ... '
  ! Repeat solver ... 
  do i=1,n_runs
    abserr = sqrt ( epsilon ( abserr ) )
    relerr = sqrt ( epsilon ( relerr ) )
    flag = 1
    t_start = 0.0d0
    t       = 0.0d0
    t_out   = 0.0d0

    !!! Initialisation
    call initRep( mdelflag, eta, y, E )

    ! Call initial value of the replicator RHS function
    y0 = y
    call rFun ( t, y, yp )
    do i_step = 1, n_step ! Call solver for t = (1, 2, ... , Nt)dt
      t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start &
          + real (          i_step - 1, kind = 8 ) * t_stop ) &
          / real ( n_step,              kind = 8 )

      t_out = ( real ( n_step - i_step, kind = 8 ) * t_start &
              + real (          i_step, kind = 8 ) * t_stop ) &
              / real ( n_step,          kind = 8 )

      call r8_rkf45( rFun, neqn, y, yp, t, t_out, relerr, abserr, flag )
    end do ! i_step=1, n_step
    write ( 1, '(2x,I2,2x,3F10.6,3g18.8E3)' ) flag, y0(1), y0(2), y0(3), y(1), y(2), y(3)
  end do ! i=1, n_runs

  ! Close data files
  close( 1 )       ! Species frequencies
  close( 2 )       ! Xi / gamma coordinates
  close( 3 )       ! Lambda-difference quantities

  return
end subroutine repRSolv

subroutine rSolv( t_stop, n_step, mdelflag )
!*****************************************************************************80
!! RSOLV solves the replicator ODE system only once and prints trajectory.
  implicit none
  integer ( kind = 4 ), parameter :: neqn = 3
  integer ( kind = 4 )            :: k, flag, i_step, n_step, mdelflag
  external rFun 
  double precision :: y(neqn), yp(neqn) ! Solver params
  double precision :: abserr, relerr, t, t_out, t_start, t_stop
  double precision :: E(neqn, 0:neqn)   ! Model params
  double precision :: gam, kap, eta
  COMMON /metabs/ eta, gam, E ! Shared with RHS function f(t,y,yp)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8_RKF45 solution of 3D replicator system:'
  write ( *, '(a)' ) '  DY = Y * ( f(Y) - dot( Y, f(Y) ) )'

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  !!! Initialisation
  call initRep( mdelflag, eta, y, E )

  ! Print E1, E2 to sdout
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  E_i' 
  write ( *, '(4x,3F10.6)' ) E(1,0), E(2,0), E(3,0)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  E_ij' 
  write ( *, '(4x,F9.6,2x,F9.6,2x,F9.6)' ) E(1,1), E(1,2), E(1,3)
  write ( *, '(4x,F9.6,2x,F9.6,2x,F9.6)' ) E(2,1), E(2,2), E(2,3)
  write ( *, '(4x,F9.6,2x,F9.6,2x,F9.6)' ) E(3,1), E(3,2), E(3,3)

  ! Initialise system constants
  gam = 0.03d0
  kap = 0.25d0
  eta = kap/(kap +gam)

  ! Initialisations for solver
  flag    = 1
  t_start = 0.0d0
  t       = 0.0d0
  t_out   = 0.0d0

  ! Call initial value of the replicator RHS function
  call rFun ( t, y, yp )

  ! Print initial state Y0
  open( unit = 4, file = './repli3YT_single_expE' )
  write ( 4, '(a)' ) '  FLAG     T       Y(1)       Y(2)       Y(3)'
  write ( 4, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2), y(3)

  do i_step = 1, n_step
    t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start &
        + real (          i_step - 1, kind = 8 ) * t_stop ) &
        / real ( n_step,              kind = 8 )

    t_out = ( real ( n_step - i_step, kind = 8 ) * t_start &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )

    call r8_rkf45( rFun, neqn, y, yp, t, t_out, relerr, abserr, flag )
    write ( 4, '(i4,2x,4g16.6)' ) flag, t, y(1), y(2), y(3)
  end do

  close( 4 )
  return
end subroutine rSolv

!*****************************************************************************80

subroutine rFun ( t, y, yp )
!*********************************************************************
!! RFUN evaluates the derivative for the Replicator ODE:
!        Y' = Y * ( f(y) - <f(y)> )
!
!  Parameters:
!    Input, double T,  the value of the independent variable.
!    Input, double Y,  the value of the dependent variable.
!
!    Output, double YP, the value of the derivative
!                      dY(1:NEQN)/dT.
!
  implicit none
  double precision :: t
  double precision :: y(3), phi(3), yp(3)
  double precision :: eta, gam, avePhi
  double precision :: E(3,0:3)
  COMMON /metabs/ eta, gam, E

  ! Compute fitness function values
  phi    = eta*gam*( E(:,0) +eta*MATMUL(TRANSPOSE(E(:,1:3)),y) )
  avePhi = DOT_PRODUCT( y, phi ) ! Compute average fitness
  yp     = y * ( phi - avePhi )  ! Compute rate of change
  return
end subroutine rFun

!*****************************************************************************80

subroutine initRep ( mdlflag, eta, y, E )
  implicit none
  integer          :: mdlflag, i, j
  double precision :: mu, eta
  double precision :: y(3), ytemp(2), unirand(3)
  double precision :: GamVec(3), lD22(3), ld23(3)
  double precision :: E(3,0:3)

!  ! Initialise seed for RNG
!   CALL INIT_RANDOM_SEED()

  ! Initialise (uniformly) random state Y0
  CALL RANDOM_NUMBER( ytemp )          ! NOTE! Rejection sampling needs 
  DO WHILE ( SUM( ytemp ) .GT. 1.0d0 ) !       resampling 50% of the times, 
    CALL RANDOM_NUMBER( ytemp )        !       on average. Find better.
  END DO
  y( 1:2 ) =             ytemp
  y(   3 ) = 1.0d0 -SUM( ytemp )

  ! E contains metabolism parameters as [E0 E1 E2 ... ]
  ! Flag 'mdlflag' can take the following numbers
  ! representing the different models for the parameters
  !      1       Independent uniformly distributed
  !      2       Tree model, uniformly distributed
  !      3       Independent exponentially distributed
  !      4       Tree model, exponentially distributed
  CALL RANDOM_NUMBER(E) ! Uniform parameter model
  SELECT CASE (mdlflag)
    CASE (2)            ! Uniform tree model
      DO j = 1, 3
        ! E(:,i) < E(:,0) enforced by E(:,0)<1
        E(:,j) = E(:,j)*E(:,0)
      END DO
    CASE (3)            ! Independent exponential model
      mu = 2.0d0        ! Exp(mu)-distributed E_i, E_ij. mu=2 gives mean=1/2
      E  = -LOG(E)/mu   ! T = -log(u)/mu ~ Exp(mu) when u~Uni(0,1)
    CASE (4)            ! Exponential tree model no1
      CALL RANDOM_NUMBER(unirand)
      DO j = 1, 3
        ! E(:,i) < E(:,0) enforced by unirand<1
        E(:,j) = E(:,j)*unirand
      END DO
    CASE default
  END SELECT

  ! Collect quantities 
  !      lD22(i) = e[i] -e[i+1] +eta( e[i+1,i] -e[i+1,i+1] ) 
  !      lD23(i) = e[i] -e[i+1] +eta( e[  i,i] -e[  i,i+1]   ) 
  ! per Eq (22) and Eq (23) of Lundh & Gerlee 2013
  DO i = 1,2
    lD22(i) = E(i,0) -E(i+1,0) +eta*( E(i+1,i) -E(i+1,i+1) )
    lD23(i) = E(i,0) -E(i+1,0) +eta*( E(i,i)   -E(i,i+1)   )
  END DO
  lD22(3)   = E(3,0) -E(1,0)   +eta*( E(1,3) -E(1,1) ) ! Special case for
  lD23(3)   = E(3,0) -E(1,0)   +eta*( E(3,3) -E(3,1) ) !  pair (3,1)
  write ( 3, '(6F12.6)' ) lD22(1), lD22(2), lD22(3), lD23(1), lD23(2), lD23(3)

  ! Gamma coordinates of Lundh & Gerlee (2013)
  GamVec(1:3) = lD22(1:3) / lD23(1:3) ! G12, G23, G31 coordinate
  write ( 2, '(2x,3F16.6)' ) GamVec(1), GamVec(2), GamVec(3)
  
end subroutine initRep
