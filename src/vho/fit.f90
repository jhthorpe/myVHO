!---------------------------------------------------------------------
! fit
!       - module containing subroutines for fitting and interpolation
!---------------------------------------------------------------------
MODULE fit

CONTAINS

!---------------------------------------------------------------------
! fit_splint
!       - spline subroutine from NRF77
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
SUBROUTINE fit_splint(xa,ya,y2a,n,x,y,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(1:n), INTENT(IN) :: xa,y2a,ya
  REAL(KIND=8), INTENT(INOUT) :: y
  REAL(KIND=8), INTENT(IN) :: x
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: n

  REAL(KIND=8) ::  a,b,h
  INTEGER :: k,khi,klo

  error = 0 

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
    write(*,*) "ERROR"
    write(*,*) "fit_splint  : bad xa input"
    error = 1
    return
  end if
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+&
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  return

END SUBROUTINE fit_splint

!---------------------------------------------------------------------
! fit_spline
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
SUBROUTINE fit_spline(x,y,n,yp1,ypn,y2,NMAX)
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

END SUBROUTINE fit_spline

!---------------------------------------------------------------------

END MODULE fit
!---------------------------------------------------------------------
