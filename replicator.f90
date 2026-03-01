program main
!*****************************************************************************80
!! MAIN is the main program for REPLICATOR, i.e.
!    the replicator equation function for RKF45 library.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!    05 February 2015
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

  integer          :: n_step, n_reps, mdlflag
  double precision :: t_stop

  write ( *, '(a)' ) 'REPLICATOR'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  '

  ! Make test call
  t_stop  = 1.0d6
  n_step  = 10000
  n_reps  = 10
  mdlflag = 3

  call repRSolv( t_stop, n_reps, n_step, mdlflag )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Start and end states Y0, YT are found in file ''repliY''.'
  write ( *, '(a)' ) '  Coordinates Xi found in file ''repliXi''.'

  stop
end

subroutine repRSolv( t_stop, n_runs, n_step, mdelflag )
!*****************************************************************************80
!! REPRSOLV is Repeated REPlicator SOLVer.
  implicit none
  integer ( kind = 4 ), parameter :: neqn = 2
  integer ( kind = 4 )            :: i, k, flag, i_step, n_runs, n_step, mdelflag
  external rFun 
  double precision :: y(neqn), y0(neqn), yp(neqn), Xi(2) ! Solver params
  double precision :: abserr, relerr, t, t_out, t_start, t_stop
  double precision :: E(neqn, neqn+1)   ! Model params
  double precision :: gam, kap, eta, lambda
  COMMON /metabs/ eta, gam, E ! Shared with RHS function f(t,y,yp)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8_RKF45 repeated solution of 2D replicator system:'
  write ( *, '(a)' ) '  DY = Y( f(Y) - dot(Y,f(Y)) )'
  write ( *, '(a)' ) '  Repeats        Tstop     Nsteps     dt '
  write ( *, '(4x,I6,2x,F10.0,2x,I8,2x,F10.5)' ), n_runs, t_stop, n_step, t_stop/dble(n_step)

  ! Check that the chosen model is available
  if ((mdelflag .GT. 4) .OR. (mdelflag .LT. 1)) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)', advance='no' ) '  Non-supported model flag, '
    write ( *, '(a)' ) ' Flag = 1,2,3,4 supported. '
    write ( *, '(a)' ) ' Aborting ... '
    STOP
  end if

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
  lambda = 2.0d0 ! For Exp(lambda)-distributed E_i, E_ij
  gam = 0.03d0
  kap = 0.25d0
  eta = kap/(kap +gam)
  write ( *, '(a)' ) '  Gamma     Kappa     Eta     Lambda'
  write ( *, '(4x,F5.2,2x,F5.2,2x,F8.6,2x,F5.2)' ), gam, kap, eta, lambda

  ! Open data files. File units used:
  !    1    System states y0(1:neqn) yT(neqn)
  !    2    Xi-coordinates Xi_a Xi_b        if neqn == 2
  !         as defined in Lundh&Gerlee (2013) Eqs (18), (19)
  open( unit = 1,  file = './repliY_N1e5_unifE'   )
  write ( 1, '(a)' )  '  Y0(1)     Y0(2)        Y(1)              Y(2)   '
  open( unit = 2, file = './repliXi_N1e5_expE' )
  write ( 2, '(a)' ) '   xi(1)      xi(2)'

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
    write ( 1, '(2x,I2,2x,2F10.6,2x,2g18.8E3)' ) flag, y0(1), y0(2), y(1), y(2)
  end do ! i=1, n_runs
  close( 1 )
  close( 2 )

  return
end subroutine repRSolv

subroutine rSolv( t_stop, n_step, mdelflag )
!*****************************************************************************80
!! RSOLV solves the replicator ODE system only once and prints trajectory.
  implicit none
  integer ( kind = 4 ), parameter :: neqn = 2
  integer ( kind = 4 )            :: k, flag, i_step, n_step, mdelflag
  external rFun 
  double precision :: y(neqn), yp(neqn) ! Solver params
  double precision :: abserr, relerr, t, t_out, t_start, t_stop
  double precision :: E(neqn, neqn+1)   ! Model params
  double precision :: gam, kap, eta
  COMMON /metabs/ eta, gam, E ! Shared with RHS function f(t,y,yp)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R8_RKF45 solution of 2D replicator system:'
  write ( *, '(a)' ) '  DY(1) = Y(1)( f1(Y) - dot( Y, f(Y) ) )'
  write ( *, '(a)' ) '  DY(2) = Y(2)( f2(Y) - dot( Y, f(Y) ) )'

  abserr = sqrt ( epsilon ( abserr ) )
  relerr = sqrt ( epsilon ( relerr ) )

  !!! Initialisation
! E contains metabolism parameters as [E1 E2 E3 ... ]
!                 !  E(:,1) = 0.5d0
!                 !  E(1,2) = 0.242687824361421d0
!                 !  E(1,3) = 0.478583474121473d0
!                 !  E(2,2) = 0.400140234444400d0
!                 !  E(2,3) = 0.070943169313608d0

  CALL RANDOM_NUMBER(E) ! Uniform parameter model
  if( mdelflag .EQ. 1 ) then
    E(:,2) = E(:,2)*E(:,1) ! Tree hierarchy model
    E(:,3) = E(:,3)*E(:,1)
  end if

  ! Print E1, E2 to sdout
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  E1' 
  write ( *, '(4x,F8.6)' ) E(:,1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  E2' 
  write ( *, '(4x,F8.6,2x,F8.6)' ) E(1,2), E(1,3)
  write ( *, '(4x,F8.6,2x,F8.6)' ) E(2,2), E(2,3)

  gam = 0.03d0
  kap = 0.25d0
  eta = kap/(kap +gam)
  do k = 1,neqn
    y(k) = 1.0d0/neqn
  end do

  flag = 1

  t_start = 0.0d0
  t       = 0.0d0
  t_out   = 0.0d0

  ! Call initial value of the replicator RHS function
  call rFun ( t, y, yp )

  open( unit = 1, file = './testrun' )
  write ( 1, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
  write ( 1, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

  do i_step = 1, n_step
    t = ( real ( n_step - i_step + 1, kind = 8 ) * t_start &
        + real (          i_step - 1, kind = 8 ) * t_stop ) &
        / real ( n_step,              kind = 8 )

    t_out = ( real ( n_step - i_step, kind = 8 ) * t_start &
            + real (          i_step, kind = 8 ) * t_stop ) &
            / real ( n_step,          kind = 8 )

    call r8_rkf45( rFun, neqn, y, yp, t, t_out, relerr, abserr, flag )
    write ( 1, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)
  end do
  close( 1 )

  return
end subroutine rSolv

subroutine rFun ( t, y, yp )
!*****************************************************************************80
!! RFUN evaluates the derivative for the Replicator ODE:
!        Y' = Y * ( f(y) - <f(y)> )
!           = Y * ( 1 - Y  )*( phi_a(Y) - phi_b(y) )         [2D case]
!  Parameters:
!    Input, double T, the value of the independent variable.
!    Input, double Y, the value of the dependent variable.
!
!    Output, double Y, the value of the derivative
!    dY(1:NEQN)/dT.
!
  implicit none
  double precision :: t
  double precision :: y(2), phi(2), yp(2)
  double precision :: eta, gam, avePhi
  double precision :: E(2,3)
  COMMON /metabs/ eta, gam, E

  ! Compute fitness function values
  phi    = eta*gam*( E(:,1) +eta*MATMUL(TRANSPOSE(E(:,2:3)),Y) )
  avePhi = DOT_PRODUCT( y, phi ) ! Compute average fitness
  yp     = y * ( phi - avePhi )  ! Compute rate of change
  return
end subroutine rFun

subroutine initRep ( mdlflag, eta, y, E )
  implicit none
  integer            :: mdlflag, i
  integer, parameter :: neqs = 2
  double precision :: mu, eta
  double precision :: y(neqs), ytemp(neqs-1), unirand(neqs)
  double precision :: Xi(2)
  double precision :: E(neqs,0:neqs)

  ! Initialise (uniformly) random state Y0
  CALL RANDOM_NUMBER( ytemp )      ! NOTE! Rejection sampling needs 
  DO WHILE ( SUM( ytemp ) .GT. 1 ) !       resampling 50% of the times, 
    CALL RANDOM_NUMBER( ytemp )    !       on average
  END DO
  y(1:(neqs-1)) = ytemp
  y(neqs)       = 1 -SUM(ytemp)

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
      DO i = 1, neqs
        ! E(:,i) < E(:,0) enforced by E(:,0)<1
        E(:,i) = E(:,i)*E(:,0)
      END DO
    CASE (3)            ! Independent exponential model
      mu = 2.0d0        ! Exp(mu)-distributed E_i, E_ij. mu=2 gives mean=1/2
      E  = -LOG(E)/mu   ! T = -log(u)/mul ~ Exp(mu) when u~Uni(0,1)
    CASE (4)            ! Exponential tree model
      CALL RANDOM_NUMBER(unirand)
      DO i = 1, neqs
        ! E(:,i) < E(:,0) enforced by E(:,0)<1
        E(:,i) = E(:,i)*E(:,0)
      END DO
    CASE default
  END SELECT

  ! Xi coordinates of Lundh & Gerlee (2013)
  Xi(1) = eta*(E(1,2)-E(1,1))-(E(1,0)-E(2,0))
  Xi(2) = eta*(E(2,1)-E(2,2))-(E(2,0)-E(1,0)) 
  write ( 2, '(2x,2F12.6)' ) Xi(1), Xi(2)
  
end subroutine initRep
