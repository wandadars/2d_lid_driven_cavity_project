!************************************************************************
!*
!* Conjugate Gradient Method (CG Method)
!*
!************************************************************************
subroutine cg_method (    &     ! Conjugate Gradient Method
                       n, &     ! Size of the linear system
                       a, &     ! System matrix
                       y, &     ! right hand side
                       x, &     ! solution vector
                       fehler & ! error code
                     )
!integer, parameter :: SIZE=24
real*8, parameter :: ZERO=0.d0, MACH_EPS=2.d-16
real*8  a(1:n,1:n),x(1:n),y(1:n)
integer fehler
!************************************************************************
!* Input parameters:
!*
!* =================
!*
!* n  Size of the linear system
!*
!* a  [1..n,1..n] system matrix A. Only the upper triangle of A is
!*
!*    used.
!*
!* y  [1..n] vector of the right hand side
!*
!*
!*
!* Output parameters:
!*
!* ==================
!*
!* x  [1..n] vector giving the solution
!*
!*
!*
!* Return value:
!*
!* =============
!*
!* = 0: all is ok
!*
!* = 1: n < 2 or other disallowed input parameters
!*
!* = 2: memory exceeded
!*
!*
!*
!************************************************************************
  real*8 d(1:n), &   ! (1..n) auxiliary vectors d and g
         g(1:n), &
         AmalD(1:n)  ! (1..n) auxiliary vector A * d
  real*8 alpha,    &   ! coefficient
         beta,     &   ! coefficient
         dividend, &   ! numerator and denominator of a fraction
         divisor,  &   ! respectively, used to compute alpha, beta
         hilf,     &   ! auxiliary variables
         hilf2,    &
         abstand,  &   ! distance of two successive approximations
                       ! for the solution vector x (taken in the
                       ! euclidean norm)
         xnorm         ! euklidean norm of x
  integer k, i, j      ! loop variables

  if (n < 2) then      ! invalid parameter?
    fehler=1
        return
  end if

  !------------------------------------------------------------------
  ! start with x at the origin
  !------------------------------------------------------------------
  do i = n, 1, -1
    x(i) = ZERO
  end do

  !------------------------------------------------------------------
  ! initialize  d and g :
  ! d = -g = -(a*x - y) = y (since x = 0)
  !------------------------------------------------------------------
  do i = n, 1, -1
    hilf = y(i)
    d(i) = hilf
    g(i) = -hilf
  end do


  !------------------------------------------------------------------
  ! perform at most n steps of the CG Method
  !------------------------------------------------------------------
  do k = n+1, 1, -1

    !----------------------------------------------------------------
    ! compute new alpha:
    ! alpha = -(d(transp) * g) / (d(transp) * (a * d))
    !----------------------------------------------------------------

    dividend = ZERO
    divisor  = ZERO

    do i = n, 1, -1
      dividend = dividend + d(i) * g(i)
          hilf = ZERO
      do j = 1, i
        hilf = hilf + a(j,i) * d(j)
      end do
      do j = i+1, n
        hilf = hilf + a(i,j) * d(j)
      end do
      AmalD(i) = hilf
      divisor = divisor + d(i) * hilf
    end do

    if (divisor.eq.ZERO) then
          fehler=0
      return
    end if

    alpha = -dividend / divisor

    !----------------------------------------------------------------
    ! compute the norm of x und  alpha * d  and find a new x:
    ! x  =  x + alpha * d, then check whether x is close enough,
    ! in order to stop the process before n complete steps
    !----------------------------------------------------------------
    xnorm   = ZERO
    abstand = ZERO

    do i = n, 1, -1
      hilf =  x(i)
      xnorm   = xnorm + hilf*hilf
      hilf2   =  alpha * d(i)
      abstand = abstand + hilf2*hilf2
      x(i)    =  hilf + hilf2
    end do

    if (abstand < MACH_EPS * xnorm) then
          fehler=0
      return
    end if


    !----------------------------------------------------------------
    ! compute new g:   g  =  g + alpha * (a * d)
    !----------------------------------------------------------------
    do i = n, 1, -1
      g(i) = g(i) + alpha * AmalD(i)
    end do

    !----------------------------------------------------------------
    ! compute new beta :
    ! beta = (g(transp) * (a * d)) / (d(transp) * (a * d))
    !----------------------------------------------------------------
    dividend = ZERO
    do i = n, 1, -1
      dividend = dividend + g(i) * AmalD(i)
    end do

    beta = dividend / divisor

    !----------------------------------------------------------------
    ! compute new d :   d  =  - g + beta * d
    !----------------------------------------------------------------
    do i = n, 1, -1
      d(i) = -g(i) + beta * d(i)
    end do

  end do  !k loop

  fehler=0
  return
end subroutine cg_method

