!---------------------------------------------------------------------
! gauss
!       - module containing subroutines for gaussian quadrature
!---------------------------------------------------------------------
MODULE gauss
  USE input
  USE fname

CONTAINS
!---------------------------------------------------------------------
! guass_generate
!       - generates or reads gaussian quadrature abscissa and weights
!---------------------------------------------------------------------
! job           : int, job type
! bas           : int, basis type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB
! error         : int, error code
! nabs          : int, number of abscissa 
! q             : 1D real*8, list of abscissa
! W             : 1D real*8, list of weights 

SUBROUTINE gauss_generate(job,bas,ndim,nbas,mem,nabs,q,W,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: q,W
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error,nabs
  INTEGER, INTENT(IN) :: job,bas,ndim

  REAL(KINd=8), DIMENSION(0:ndim-1) :: line
  INTEGER ::i

  error = 0
  IF (bas .EQ. 1) THEN
    WRITE(*,*) "Generating Gauss-Hermite abscissa and weights"
    WRITE(*,*) 
    nabs = MAXVAL(nbas)+10

    WRITE(*,*) "Number of abscissa to be used", nabs
    ALLOCATE(q(0:nabs-1))
    ALLOCATE(W(0:nabs-1))
    q = 0.0D0
    W = 0.0D0

    CALL gauss_hermite(nabs,q,W,error)  
    IF (error .NE. 0) RETURN
    IF (job .EQ. -2) THEN
      CALL gauss_2_xyz(ndim,nabs,q,error)
    END IF

  ELSE
    WRITE(*,*) "gauss_generate  : ERROR"
    WRITE(*,*) "Sorry, only HO basis is coded now" 
    error = 1
    RETURN
  END IF

  WRITE(*,*) "Writting abscissa and weights to abscissa.dat"
  WRITE(*,*) 
  OPEN(file='abscissa.dat',unit=101,status='replace')
  DO i=0,nabs-1
    WRITE(101,*) i+1,q(i),W(i)
  END DO
  CLOSE(unit=101)
  WRITE(*,*) "----------------------------------------------"
  WRITE(*,*)

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
! gauss_2_xyz
!       - converts gaussian quadrature to xyz coordinates
!       - reads QUADRATURE file from CFOUR
!---------------------------------------------------------------------
! ndim          : int, number of dimensions
! nabs          : int, number of abscissa
! q             : 1D real*8, abscsissa
! error         : int, exit code

SUBROUTINE gauss_2_xyz(ndim,nabs,q,error)
  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nabs

  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: qmat
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: q0,xyz
  REAL(KIND=8), DIMENSION(0:ndim-1) :: basK
  CHARACTER(LEN=1024) :: fname
  LOGICAL :: ex
  INTEGER :: natoms,fid,i,j,k
  
  error = 0
  natoms = -1
  
  INQUIRE(file='QUADRATURE',EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "gauss_2_xyz  : ERROR"
    WRITE(*,*) "Could not find QUADRATURE output from cfour"
    error = 1
    RETURN
  END IF

  CALL input_QUAD_natoms(natoms,error)

  ALLOCATE(q0(0:natoms-1,0:2))
  ALLOCATE(qmat(0:natoms-1,0:2,0:ndim-1))
  ALLOCATE(xyz(0:natoms-1,0:2))
  
  CALL input_QUAD_q0(natoms,q0,error)
  IF (error .NE. 0) RETURN
  CALL input_QUAD_bask(ndim,bask,error)
  IF (error .NE. 0) RETURN
  CALL input_QUAD_qmat(ndim,natoms,qmat,error)
  IF (error .NE. 0) RETURN

  WRITE(*,*) "The following normal coordinates were read"
  DO j=0,ndim-1
    WRITE(*,*)
    WRITE(*,'(1x,A11,1x,F16.10)') "Frequency :", bask(j)
    DO i=0,natoms-1
      WRITE(*,'(1x,3(f14.10,2x))') qmat(i,0:2,j)
    END DO  
  END DO
  WRITE(*,*) 

  WRITE(*,*) "Writing points to 'pointsX.txt'"
  WRITE(*,*) "Format is atom1_x,atom1_y,atom1_z,atom2_x,...."
  DO j=0,ndim-1
    fid = 500 + j
    CALL fname_pointstxt(j+1,fname,error)
    OPEN(file=TRIM(fname),unit=fid,status='replace')
    DO i=0,nabs-1 
      xyz = q0 + q(i)*qmat(0:natoms-1,0:2,j)
      !This is just hilariously bad code, but who cares
      WRITE(fid,*) TRANSPOSE(xyz)
    END DO
    CLOSE(unit=fid)
  END DO

END SUBROUTINE gauss_2_xyz

!---------------------------------------------------------------------
END MODULE gauss
