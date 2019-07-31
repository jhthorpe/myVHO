!//////////////////////////////////////////////////////////////////
!///
!///            Module that contains subroutines for 
!///            fitting and evaluating polynomial functions
!///            to an arbitrary order
!///
!//////////////////////////////////////////////////////////////////

MODULE fit
  USE linal
  USE val
  USE nints

CONTAINS

!---------------------------------------------------------------------
! fit
!       -fits surface with desired fitting functions
!---------------------------------------------------------------------
! Variables
! func          : int, function type 
! Fx            : 1D real*8, function to fit
! x             : 1D real*8, list of x values
! np            : int, number of points
! tol           : real*8, RMS target tollerance
! fit           : int, type of basis functions
! ord           : int, on exit, order of the polynomial used 
! coef          : 2D real*8, list of coefficients used
! a             : real*8, alpha
! m             : real*8, mass
! nb            : int, number of HO basis functions
! qeq           : real*8, equilibrium position
! error         : bool, true on exit if there was a problem

SUBROUTINE fit_surf(func,Fx,x,np,tol,ord,coef,a,m,nb,qeq,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: coef
  REAL(KIND=8), DIMENSION(0:np-1), INTENT(INOUT) :: Fx,x
  REAL(KIND=8), INTENT(IN) :: tol,a,m,qeq
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(INOUT) :: ord
  INTEGER, INTENT(IN) :: np,func,nb

  WRITE(*,*) &
"---------------------------------------------------------------------"
  WRITE(*,*) "Fitting the surface"
  WRITE(*,*) 
  WRITE(*,*) "Fitting options are:"
  WRITE(*,*) " -1             : Display abscissa needed"
  WRITE(*,*) "  0             : Precalculated abscissa"
  WRITE(*,*) "  1             : a*exp(-b*r)"
  WRITE(*,*) "  2             : a*x^i"
  WRITE(*,*) "  3             : cubic spline"
  WRITE(*,*) "  4             : no fitting" 
  WRITE(*,*)
  IF (func .EQ. -1) THEN
    WRITE(*,*) "Option (-1) selected"
    CALL absc_calc(x,Fx,a,m,nb,qeq,error)
    RETURN

  ELSE IF (func .EQ. 1) THEN
    WRITE(*,*) "Option (1) selected"
    WRITE(*,*) "Sorry, that method not yet implimented"
    WRITE(*,*) 
    error = .TRUE.
    RETURN
    !CALL exp_fit(Fx(0:N-1),x(0:N-1),N,tol,ord,coef,error)

  ELSE IF (func .EQ. 2) THEN
    WRITE(*,*) "Option (2) selected"
    CALL poly_fit(Fx(0:np-1),x(0:np-1),np,tol,ord,coef,error)

  ELSE IF (func .EQ. 0) THEN
    WRITE(*,*) "Option (0) selected"
    WRITE(*,*) "Sorry, that option not implemented yet"
    WRITE(*,*) 
    error = .TRUE.
    RETURN

  ELSE IF (func .EQ. 3) THEN
    WRITE(*,*) "Option (3) selected"
    CALL spline_fit(Fx(0:np-1),x(0:np-1),np,a,m,nb,qeq,error)

  ELSE IF (func .EQ. 4) THEN
    WRITE(*,*) "Option (4) selected"
    RETURN

  ELSE
    WRITE(*,*) 
    WRITE(*,*) "That fitting function not supported"
    WRITE(*,*) 
    error = .TRUE. 
    RETURN
  END IF

END SUBROUTINE fit_surf

!---------------------------------------------------------------------
! absc_calc 
!       -precalcuate how many and which values of abscissa will be
!        needed to fit nb HO basis functions
!---------------------------------------------------------------------
! Variables
! absc          : 1D real*8, x values of abscissa
! W             : 1D real*8, weights of abscissa
! a             : real*8, alpha
! m             : real*8, mass
! nb            : int, number of basis functions
! error         : bool, true on exit if problem

SUBROUTINE absc_calc(absc,W,a,m,nb,qeq,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: absc,W
  REAL(KIND=8), INTENT(IN) :: a,m,qeq
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nb

  REAL(KIND=8) :: A2B
  INTEGER :: na,i

  error = .FALSE.
  A2B = 1.88973

  na = nb+10 !this is, I think, probably overkill?"
  !na = 2
  !na = 1 

  WRITE(*,*)
  WRITE(*,*) "Calculating abscissa and weights needed",&
              " to evaluate HO integrals"
  WRITE(*,*) "Number of basis           :", nb
  WRITE(*,*) "Number of abscissa        :", na

  CALL gauher(absc(0:na-1),W(0:na-1),na,error)
  IF (error) THEN
    WRITE(*,*) "HO1D_abscissa -- error from gauher"
    error = .TRUE.
    RETURN
  END IF

  absc = absc/a

  WRITE(*,*) "Writting abscissa and weights to abscissa.dat"
  WRITE(*,*) "These values are shifted from a qeq of", qeq/A2B
  WRITE(*,*)

  OPEN(file='abscissa.dat',unit=101,status='replace')
  WRITE(101,*) "point     x(Å)      W"
  DO i=0,na-1
    WRITE(101,'(2x,I4,2x,F20.16,2x,ES24.15)') i,(absc(i)+qeq)/A2B,W(i)
    WRITE(*,'(2x,I4,2x,F20.16,2x,ES24.15)') i,(absc(i)+qeq)/A2B,W(i)
  END DO
  CLOSE(unit=101)

END SUBROUTINE absc_calc 

!---------------------------------------------------------------------
! spline_fit
!       - uses cubic spline interpolation to obtain the needed values
!         of Vx and x for Gaussian Quadrature
!       - note that we are working in y = a*x space, hence the weird
!         factors of a
!---------------------------------------------------------------------
! Variables
! Fx            : 1D real*8, function to fit
! x             : 1D real*8, list of x values
! np            : int, number of points
! a             : real*8, alpha
! m             : real*8, mass
! nb            : int, nubmer of HO basis functions
! error         : bool, true on exit if there was a problem

SUBROUTINE spline_fit(Fx,x,np,a,m,nb,qeq,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:np-1), INTENT(IN) :: Fx,x
  REAL(KIND=8), INTENT(IN) :: a,m,qeq
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: np,nb

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ab,Vab,W,y2
  REAL(KIND=8) :: qmin,qmax,yp1,ypn,h,yf,yb
  INTEGER :: na,i

  error = .FALSE.

  WRITE(*,*) 

  na = nb+10 !probably an overestimation?
  na = 13
  WRITE(*,*) "FIXED ABSCISSA AT ",na
  !na = 50 !probably an overestimation?
  !na = 10
  !na = 2
  ALLOCATE(Vab(0:na-1))
  ALLOCATE(ab(0:na-1))
  ALLOCATE(W(0:na-1))
  ALLOCATE(y2(0:np-1))
  qmin = x(1)
  qmax = x(np-2)
 
  WRITE(*,*) "Determining location of abscissa"
  CALL gauher(ab(0:na-1),W(0:na-1),na,error) !this is in y=a*x space
  IF (error) THEN
    WRITE(*,*) "fit:spline_fit -- error out of gauher"
    RETURN
  END IF

  WRITE(*,*) "Performing cubic spline interpolation of surface"
  !derivatives for spline
  yp1 = (Fx(2) - Fx(0))/(x(2)-x(0))
  ypn = (Fx(np-1) - Fx(np-3))/(x(np-1)-x(np-3)) 
  
  CALL spline(x(1:np-2),Fx(1:np-2),np-2,yp1,ypn,y2(1:np-2),np-2)
  DO i=0,na-1
    IF (ab(i)/a .LT. qmin .OR. ab(i)/a .GT. qmax) THEN
      Vab(i) = 0.0D0 
      !WRITE(*,*) "WARNING, outside of qmin or qmax", ab(i),ab(i)/a
    ELSE
      CALL splint(x(1:np-2),Fx(1:np-2),y2(1:np-2),np,ab(i)/a,Vab(i),error)
    END IF
    IF (error) THEN
      WRITE(*,*) "fit:spline_fit -- error out of splint" 
      DEALLOCATE(Vab)
      DEALLOCATE(ab)
      DEALLOCATE(W)
      DEALLOCATE(y2)
      RETURN
    END IF
  END DO 

  !write out abscissa and values to intermediate file
  WRITE(*,*) "Cublic spline interpolation stored in gauss.dat"
  OPEN(file='gauss.dat',unit=110,status='replace')
  WRITE(110,*) na
  DO i=0,na-1
    WRITE(110,*) ab(i), W(i), Vab(i)
  END DO  
  CLOSE(unit=110)

  DEALLOCATE(Vab)
  DEALLOCATE(ab)
  DEALLOCATE(W)
  DEALLOCATE(y2)

END SUBROUTINE spline_fit

!---------------------------------------------------------------------
! splint
!       -splilnt subroutine from NRF77
!Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC
!COMPUTING (ISBN 0-521-43064-X)
!Copyright (C) 1986-1992 by Cambridge University Press.
!Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
!Permission is granted for internet users to make one paper copy for their own
!personal use. Further reproduction, or any copying of machinereadable
!files (including this one) to any server
!computer, is strictly prohibited. To order Numerical Recipes books,
!diskettes, or CDROMs
!visit website http://www.nr.com or call 1-800-872-7423 (North America only),
!or send email to trade@cup.cam.ac.uk (outside North America).
!---------------------------------------------------------------------
SUBROUTINE splint(xa,ya,y2a,n,x,y,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(1:n), INTENT(IN) :: xa,y2a,ya
  REAL(KIND=8), INTENT(INOUT) :: y
  REAL(KIND=8), INTENT(IN) :: x
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: n

  REAL(KIND=8) ::  a,b,h
  INTEGER :: k,khi,klo

  error = .FALSE.

  klo=1
  khi=n
! 1 if (khi-klo.gt.1) then
!    k=(khi+klo)/2
!    if(xa(k).gt.x)then
!      khi=k
!    else
!      klo=k
!    endif
!  goto 1
!  endif

!my own, really really stupid search for now
  do k=1,n
    if (xa(k) .GT. x) then
      klo = k-1
      khi = k
      EXIT
    end if
  end do
 
  IF (khi .EQ. 1) THEN
    WRITE(*,*) "The search for ",x
    WRITE(*,*) "ended with khi = 1"
  END IF

  h=xa(khi)-xa(klo)
  if (h.eq.0.) then
    write(*,*) "fit:splint -- bad xa input"
    error = .TRUE.
    return
  end if
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+&
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  return

END SUBROUTINE splint

!---------------------------------------------------------------------
! spline
!       -spline subroutine from NRF77
!Sample page from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC
!COMPUTING (ISBN 0-521-43064-X)
!Copyright (C) 1986-1992 by Cambridge University Press.
!Programs Copyright (C) 1986-1992 by Numerical Recipes Software.
!Permission is granted for internet users to make one paper copy for their own
!personal use. Further reproduction, or any copying of machinereadable
!files (including this one) to any server
!computer, is strictly prohibited. To order Numerical Recipes books,
!diskettes, or CDROMs
!visit website http://www.nr.com or call 1-800-872-7423 (North America only),
!or send email to trade@cup.cam.ac.uk (outside North America).
!---------------------------------------------------------------------
SUBROUTINE spline(x,y,n,yp1,ypn,y2,NMAX)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(1:n), INTENT(INOUT) :: y2
  REAL(KIND=8), DIMENSION(1:n), INTENT(IN) :: x,y
  REAL(KIND=8), INTENT(IN) :: yp1,ypn
  INTEGER, INTENT(IN) :: n,NMAX

  REAL(KIND=8), DIMENSION(1:NMAX) :: u
  REAL(KIND=8) :: p,qn,sig,un
  INTEGER :: i,k

  if (yp1.gt..99e30) then
    y2(1)=0.0
    u(1)=0.0
  else
    y2(1)=-0.5
    u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  end if
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.
    y2(i)=(sig-1.)/p
    u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))&
          /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo
  if (ypn.gt..99e30) then
    qn=0.
    un=0.
  else
    qn=0.5
    un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo
  return

END SUBROUTINE spline
!---------------------------------------------------------------------
! exp_fit
!       - subroutine that fits a function to a sum of exponentials
!       - employs the SVD algorithm of Zeiger-McEwan and Kung
!---------------------------------------------------------------------
! Variables
! Fx            : 1D real*8, function to fit
! x             : 1D real*8, list of x values
! N             : int, number of points
! tol           : real*8, RMS target tollerance
! fit           : int, type of basis functions
! ord           : int, on exit, order of the polynomial used 
! coef          : 2D real*8, list of coefficients used
! error         : bool, true on exit if there was a problem

SUBROUTINE exp_fit(Fx,x,N,tol,ord,coef,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: coef
  REAL(KIND=8), DIMENSION(0:N-1), INTENT(IN) :: Fx,x
  REAL(KIND=8), INTENT(IN) :: tol
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(INOUT) :: ord
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: H,U,S,Vt
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: dd
  REAL(KIND=8) :: rms,rms_old,au2cm
  INTEGER :: i,j,nc

  error = .FALSE.
  au2cm = 2.1947E5

  WRITE(*,*) "Generating exponential fit of surface"
  WRITE(*,*) "Using Σa*exp(-b*r)"
  WRITE(*,*) "RMS Convergence :", tol

  ord = 2
  ALLOCATE(H(0:N-1,0:ord-1)) 

  !loop until convergence is reached
  !DO i=0,N/2-1
    
  !END DO

  CALL exp_clean(H,U,S,Vt,dd)
  

END SUBROUTINE exp_fit
!---------------------------------------------------------------------
! exp_clean
!       -cleans memory for exp subroutine
!---------------------------------------------------------------------
SUBROUTINE exp_clean(H,U,S,Vt,dd)
  IMPLICIT NONE
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: H,U,S,Vt
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: dd

  IF (ALLOCATED(H)) DEALLOCATE(H)
  IF (ALLOCATED(U)) DEALLOCATE(S)
  IF (ALLOCATED(Vt)) DEALLOCATE(Vt)
  IF (ALLOCATED(dd)) DEALLOCATE(dd)
  
END SUBROUTINE exp_clean

!---------------------------------------------------------------------
! poly_fit
!       - subrotuine that fits a function to an increasing number of
!         polynomials, until a threshold RMS is met
!---------------------------------------------------------------------
! Variables
! Fx            : 1D real*8, function to fit
! x             : 1D real*8, list of x values
! N             : int, number of points
! tol           : real*8, RMS target tollerance
! fit           : int, type of basis functions
! ord           : int, on exit, order of the polynomial used 
! coef          : 2D real*8, list of coefficients used
! error         : bool, true on exit if there was a problem

SUBROUTINE poly_fit(Fx,x,N,tol,ord,coef,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: coef
  REAL(KIND=8), DIMENSION(0:N-1), INTENT(IN) :: Fx,x
  REAL(KIND=8), INTENT(IN) :: tol
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(INOUT) :: ord
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A,b
  REAL(KIND=8) :: rms,rms_old,au2cm
  INTEGER :: i,j,nc,dof

  error = .FALSE.
  au2cm = 2.1947E5

  WRITE(*,*) 
  WRITE(*,*) "Generating polynomial fit of surface"
  WRITE(*,*) "Using Σa*x^i"
  WRITE(*,*) "RMS Convergence :", tol

  ALLOCATE(A(0:N-1,0:N-1))
  dof = 1
  
  ALLOCATE(b(0:N-1,0:0))
  A = 0
  b = 0

  !loop until we get below tolerance
  !DO i=0,1!N-1 
  DO i=0,N/2-1
    ord = i
    IF (ord .EQ. 0) THEN
      nc = 1
      A(0:N-1,ord) = (/ (1.0D0,j=0,N-1) /) 
      CALL lsqr(A(0:N-1,0:nc-1),Fx(0:N-1),N,nc,dof,b(0:nc-1,0:dof-1),error) 
      IF (error) THEN
        WRITE(*,*) "poly:poly_fit -- lsqr exited with error" 
        error = .TRUE. 
        DEALLOCATE(A)
        DEALLOCATE(b)
        RETURN
      END IF 
      CALL poly_rms(Fx,x,N,dof,b(0:nc-1,0:dof-1),nc,ord,rms,error)
      IF (error) THEN
        WRITE(*,*) "poly:poly_fit -- poly_rms exited with error" 
        error = .TRUE. 
        DEALLOCATE(A)
        DEALLOCATE(b)
        RETURN
      END IF 
    ELSE
      nc = ord + 1
      A(0:N-1,ord) = (/ (A(j,ord-1)*x(j) ,j=0,N-1) /) 
      CALL lsqr(A(0:N-1,0:nc-1),Fx(0:N-1),N,nc,dof,b(0:nc-1,0:dof-1),error) 
      IF (error) THEN
        WRITE(*,*) "poly:poly_fit -- lsqr exited with error" 
        error = .TRUE. 
        DEALLOCATE(A)
        DEALLOCATE(b)
        RETURN
      END IF 
      CALL poly_rms(Fx,x,N,dof,b(0:nc-1,0:dof-1),nc,ord,rms,error)
      IF (error) THEN
        WRITE(*,*) "poly:poly_fit -- poly_rms exited with error" 
        error = .TRUE. 
        DEALLOCATE(A)
        DEALLOCATE(b)
        RETURN
      END IF 
    END IF

!    WRITE(*,*) "At iteration ", ord," rms is", rms
    IF (rms .LT. tol) EXIT
    IF (rms .GT. rms_old .AND. i .GT. 1) EXIT
    rms_old = rms
    
  END DO

  IF (rms .GT. tol) THEN
    error = .TRUE.
    WRITE(*,*) "poly:poly_fit -- failed to get below RMS target"
    WRITE(*,*) 
    WRITE(*,'(2x,A24,I4,A11)') "RMS did not converge in ", ord, " iterations"
    WRITE(*,'(2x,A4,2x,F11.8,2x,A2,2x,F8.4,2x,A4)') &
                "RMS=",rms,"au",rms*au2cm,"cm-1"
    ALLOCATE(coef(0:nc-1,0:dof-1))
    coef(0:nc-1,0:dof-1) = b(0:nc-1,0:dof-1)
  ELSE
    WRITE(*,*) 
    WRITE(*,'(2x,A17,I4,A11)') "RMS converged in ", ord, " iterations"
    WRITE(*,'(2x,A4,2x,F11.8,2x,A2,2x,F8.4,2x,A4)') &
                "RMS=",rms,"au",rms*au2cm,"cm-1"
    ALLOCATE(coef(0:nc-1,0:dof-1))
    coef(0:nc-1,0:dof-1) = b(0:nc-1,0:dof-1)
  END IF
  OPEN(file='fit.dat',unit=105,status='replace')
  WRITE(105,*) ord
  DO i=0,ord
    WRITE(105,*) i, b(i,0)
  END DO
  CLOSE(unit=105)
  WRITE(*,*) "Coefficients written to fit.dat"

  CALL poly_plot(Fx,x,N,ord,coef,error)

  DEALLOCATE(A)
  DEALLOCATE(b)

END SUBROUTINE poly_fit

!---------------------------------------------------------------------
! poly_plot
!       - plots the polynomial fit of Vq
!---------------------------------------------------------------------
SUBROUTINE poly_plot(Fx,x,N,ord,coef,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: coef
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Fx,x
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,ord

  INTEGER :: i

  error = .FALSE.
  WRITE(*,*) "Saving fit runfile to fitplot"

  OPEN(unit=101,file='plot.dat',status='replace')
  DO i=0,N-1
    WRITE(101,*) x(i),Fx(i)
  END DO
  CLOSE(unit=101)

  OPEN(unit=100,file='fitplot',status='replace')
  WRITE(100,*) "set terminal png"
  WRITE(100,*) "set output 'fit.png'"
  WRITE(100,*) "set xrange [", x(0), ":", x(N-1),"]"
  WRITE(100,*) "set yrange [", MINVAL(Fx),":",MAXVAL(Fx),"]"
  WRITE(100,*) "plot 'plot.dat' u 1:2 t 'Vq',\"
  DO i=0,ord-1
    WRITE(100,*) coef(i,0),"*x**",i," + \"
  END DO
  WRITE(100,*) coef(ord,0),"*x**",ord
  CLOSE(unit=100)

END SUBROUTINE poly_plot
!---------------------------------------------------------------------
! poly_rms
!       - returns rms of poynomial fit
!---------------------------------------------------------------------
SUBROUTINE poly_rms(Fx,x,N,dof,coef,nc,ord,rms,error)
  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: coef
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Fx,x
  REAL(KIND=8), INTENT(INOUT) :: rms
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,ord,dof,nc

  REAL(KIND=8) :: val
  INTEGER :: i

  error = .FALSE.

  rms = 0.0D0
  OPEN(file='rms.dat',unit=107,status='replace')
  DO i=0,N-1
    CALL poly_eval(ord,coef(0:nc-1,0:dof-1),x(i),val,error) 
    IF (error) RETURN
    WRITE(107,*) x(i), SQRT((Fx(i) - val)**2.0/N)
    rms = rms + (Fx(i) - val)**2.0
  END DO
  rms = SQRT(rms/N)
  CLOSE(unit=107)

END SUBROUTINE poly_rms

!---------------------------------------------------------------------
! poly_eval
!       - given order and coefficients, returns the value of 
!         a fitted polynomial at x
!---------------------------------------------------------------------
SUBROUTINE poly_eval(ord,coef,x,val,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: coef
  REAL(KIND=8), INTENT(INOUT) :: val
  REAL(KIND=8), INTENT(IN) :: x
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ord

  INTEGER :: i
  
  error = .FALSE.

  val = 0
  DO i=0,ord !remember that ord = 1 is a 1st order polynomial, a*x + b
    val = val + coef(i,0)*x**i
  END DO

  !check for inf/nan
  CALL checkval(val,error)
  IF (error) THEN
    WRITE(*,*) "poly:poly_eval -- bad value in evalulation" 
    RETURN
  END IF

END SUBROUTINE poly_eval

!---------------------------------------------------------------------

END MODULE fit
