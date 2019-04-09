!//////////////////////////////////////////////////////////////////
!///
!///            Module containing numerical integration 
!///            subroutines
!///
!//////////////////////////////////////////////////////////////////

MODULE nints

CONTAINS

!---------------------------------------------------------------------
! gauher
!       - calculates Gauss-Herminet abscissas and weights for ints of
!         type exp(-a^2 x^2)f(x)
!
!       - subroutine largely taken from NR77
!
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
!
!---------------------------------------------------------------------
! Variables
! x             : 1D real*8, list of abscissa
! w             : 1D real*8, list of weights
! n             : int, number of points in integration

SUBROUTINE gauher(x,w,n,error)
  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(1:n), INTENT(INOUT) :: x,w
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: n

  REAL(KIND=8) :: EPS, PIM4, p1,p2,p3,pp,z,z1
  INTEGER :: i,its,j,m,MAXIT
  error = .FALSE.
  EPS = 3.0D-14
  PIM4 = 0.7511255444649425D0
  MAXIT = 10

  m = (n+1)/2

  DO i=1,m
    IF (i .EQ. 1) THEN
      z = SQRT(FLOAT(2*n+1)) - 1.85575*(2*n+1)**(-0.16667)
    ELSE IF (i .EQ. 2) THEN
      z=z-1.14*n**.426/z
    ELSE IF (i .EQ. 3) THEN
      z=1.86*z-.86*x(1)
    ELSE IF (i .EQ. 4) THEN
      z=1.91*z-.91*x(2)
    ELSE 
      z=2.*z-x(i-2)
    END IF

    DO its=1,MAXIT
      p1 = PIM4
      p2 = 0.0D0
      DO j=1,n
        p3 = p2
        p2 = p1
        p1 = z*sqrt(2.d0/j)*p2-sqrt(dble(j-1)/dble(j))*p3
      END DO
      pp=sqrt(2.d0*n)*p2
      z1=z
      z=z1-p1/pp
      if(abs(z-z1).le.EPS) EXIT
      IF (its .EQ. MAXIT) THEN
        WRITE(*,*) "ERROR"
        WRITE(*,*) "Too many iterations in nint:gauher"
        error = .TRUE.
        RETURN
      END IF
    END DO
    x(i)=z
    x(n+1-i)=-z
    w(i)=2.d0/(pp*pp)
    w(n+1-i)=w(i)
  END DO

END SUBROUTINE gauher 
!---------------------------------------------------------------------
END MODULE nints
