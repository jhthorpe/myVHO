!---------------------------------------------------------------------
! gauss
!       - module containing subroutines for gaussian quadrature
!---------------------------------------------------------------------
MODULE gauss

CONTAINS
!---------------------------------------------------------------------
! guass_generate
!       - generates or reads gaussian quadrature abscissa and weights
!---------------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB
! error         : int, error code
! nabs          : 1D int, number of abscissa for each dimension
! q             : 1D real*8, list of abscissa [abscissa,dimension]
! W             : 1D real*8, list of weights  [abscissa,dimension]

SUBROUTINE gauss_generate(job,ndim,nbas,mem,nabs,q,W,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: q,W
  INTEGER, DIMENSION(:), ALLOCATABLE :: nabs
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim

  REAL(KINd=8), DIMENSION(0:ndim-1) :: line
  INTEGER ::i,j,cnt

  error = 0
  WRITE(*,*) 
  WRITE(*,*) "gauss_generate : called"
  WRITE(*,*) "Generating Gauss-Hermite abscissa and weights"
  
  ALLOCATE(nabs(0:ndim-1))
  DO i=0,ndim-1
    nabs(i) = nbas(i) + 10
  END DO

  WRITE(*,*) "Number of abscissa to be used", nabs
  ALLOCATE(q(0:MAXVAL(nabs)-1,0:ndim-1))
  ALLOCATE(W(0:MAXVAL(nabs)-1,0:ndim-1))
  q = 0
  W = 0

  DO j=0,ndim-1
    CALL gauss_hermite(nabs(j),q(0:nabs(j)-1,j),W(0:nabs(j)-1,j),error)  
  END DO

  WRITE(*,*) "Writting abscissa and weights to abscissa.dat"
  WRITE(*,*) 
  OPEN(file='abscissa.dat',unit=101,status='replace')
  cnt = 0
  DO j=0,ndim-1
    line = 0.0D0
    DO i=0,nabs(j)-1
      cnt = cnt + 1
      line(j) = q(i,j)
      WRITE(101,*) cnt,line,W(i,j)
    END DO
  END DO 
  CLOSE(unit=101)

END SUBROUTINE gauss_generate

!---------------------------------------------------------------------
! gauss_hermite
!       - calculates abscissa and weights for gauss-hermite 
!         integration
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
! n             : int, number of abscissa
! x             : 1D real*8, list of abscissa
! w             : 1D real*8, list of weights
! error         : int, error code

SUBROUTINE gauss_hermite(n,x,w,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(1:), INTENT(INOUT) :: x,w
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: n

  REAL(KIND=8) :: EPS, PIM4, p1,p2,p3,pp,z,z1
  INTEGER :: i,its,j,m,MAXIT
  error = 0 
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
        WRITE(*,*) "gauss_hermite  : failed on this abscissa",i
        error = i 
        RETURN 
      END IF
    END DO
    x(i)=z
    x(n+1-i)=-z
    w(i)=2.d0/(pp*pp)
    w(n+1-i)=w(i)
  END DO

END SUBROUTINE gauss_hermite
!---------------------------------------------------------------------


END MODULE gauss
