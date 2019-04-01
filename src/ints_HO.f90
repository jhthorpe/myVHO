!//////////////////////////////////////////////////////////////////
!///
!///		Module for calculating harmonic oscillator 1D 
!///		integrals
!///
!//////////////////////////////////////////////////////////////////

MODULE ints_HO
  IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------
!       HO1D_integrals 
!               -calculates 1D harmonic oscillator integrals 
!---------------------------------------------------------------------
! Variables
! N             : int, number of harmonic oscillator basis functions
! Hij		: 2D real*8, Hamiltonian Matrix
! Vq            : 1D real*8, 1D potential energy surface
! qmin          : real*8, minimum r
! qmax          : real*8, max r
! qeq           : real*8, equilibriu q
! Ncon          : 1D real*8, list of normalization constants
! error         : bool, true if error

SUBROUTINE HO1D_integrals(N,Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Vq,q
  REAL(KIND=8), INTENT(IN) :: qmin,qmax,qeq,k,m,V_off,a
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N, npoints

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Ncon
  INTEGER :: i
  error = .FALSE.

  WRITE(*,*) 
  WRITE(*,*) "Calculating HO integrals"

  ALLOCATE(Hij(0:N-1,0:N-1))
  ALLOCATE(Ncon(0:N-1))
  Hij = 0.0D0
  Ncon = 0.0D0

  CALL HO1D_potential(Hij,N,Vq,q,npoints,k,m,V_off,a,error)
  CALL HO1D_harmonic(Hij,N,k,m,error)

END SUBROUTINE HO1D_integrals

!---------------------------------------------------------------------
!	HO1D_potential
!		-calculates 1D HO potential energy integrals
!		-stored upper triangular
!		-using x form of harmonic oscillator
!---------------------------------------------------------------------
! Variables
! N 	       : int, number of harmonic oscillator basis functions
! Hij		: 2D real*8, Hamiltonian Matrix
! Vq            : 1D real*8, 1D potential energy surface
! q		: 1D real*8, list of q values
! npoints	: int, number of points in Vq and q
! error		: bool, true if error

SUBROUTINE HO1D_potential(Hij,N,Vq,q,npoints,k,m,V_off,a,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Vq,q
  REAL(KIND=8), INTENT(IN) :: k,m,a,V_off
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N, npoints

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Sij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Htab 
  REAL(KIND=8) :: dq,dx,xlim,x,tol
  INTEGER :: i,j,u,xpoints

  error = .FALSE.
  WRITE(*,*) "Calculating potential energy numerical integrals"
  WRITE(*,*) "Number of points :", npoints
 
  ALLOCATE(Sij(0:N-1,0:N-1))
  ALLOCATE(Htab(0:N-1))


  !DO u=1,npoints-2
  DO u=0,npoints-2
  !DO u=0,npoints-1
   
    !Construct Hermitian table
    CALL build_Htab(N,a*q(u),Htab(0:N-1))
    dq = ABS(q(u+1) - q(u))

    DO j=0,N-1
      DO i=j,N-1
        Hij(i,j) = Hij(i,j) + dq*Htab(i)*Htab(j)*EXP(-a**2.0*q(u)**2.0)* (&
                  Vq(u) - 0.5*k*q(u)**2.0D0)
      END DO
    END DO

  END DO 

  Sij = 0
  xlim = 3
  xpoints = 500000
  dx = (2*xlim)/xpoints

  WRITE(*,*) 
  WRITE(*,*) "Calculating overlap integrals"
  WRITE(*,*) "Number of points :", xpoints 
  WRITE(*,*) "between          :", -xlim,",", xlim 
  WRITE(*,*) "with dx          :", dx

  x = -xlim

  DO u=0,xpoints-2

    CALL build_Htab(N,a*x,Htab(0:N-1))
    DO j=0,N-1
        Sij(j,j) = Sij(j,j) + Htab(j)*Htab(j)*EXP(-a**2.0*x**2.0)*dx
    END DO

  !  IF (u .EQ. xpoints/2) THEN
  !    WRITE(*,*) "at halfway point"
  !    WRITE(*,*) "x = ", x
  !    WRITE(*,*) "S00", Sij(0,0)
  !    WRITE(*,*) "S10", Sij(1,0)
  !    WRITE(*,*) "S11", Sij(1,1)
  !    WRITE(*,*) 
  !  END IF

    x = x + dx
  END DO

  DO j=0,N-1
    DO i=j,N-1
      Hij(i,j) = Hij(i,j) / SQRT((Sij(i,i)*Sij(j,j)))
    END DO
  END DO

  !WRITE(*,*) "S00", Sij(0,0)
  !WRITE(*,*) "S10", Sij(1,0)
  !WRITE(*,*) "S11", Sij(1,1)
  !WRITE(*,*) 

  DEALLOCATE(Htab)
  DEALLOCATE(Sij)

END SUBROUTINE HO1D_potential
!---------------------------------------------------------------------
!	build_Htab
!		-constructs list of Hermite polynomials at
!		a particular q
!---------------------------------------------------------------------
! Variables
! Htab		: 1D real*8, list of hermite polynomials eval. at q
! q		: real*8, value to evaluate Htab at
! N		: int, max Hermite polynomial to calc

SUBROUTINE build_Htab(N,q,Htab)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Htab 
  REAL(KIND=8), INTENT(IN) :: q
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8) :: H0, H1
  INTEGER :: i

  IF (N .EQ. 1) THEN
    Htab(0) = 1.0D0
    RETURN
  ELSE IF (N .EQ. 2) THEN
    Htab(0) = 1.0D0
    Htab(1) = 2.0*q
    RETURN
  ELSE 
    Htab(0) = 1.0D0
    Htab(1) = 2.0*q
    DO i=1,N-2
      Htab(i+1) = 2*q*Htab(i) - 2.*i*Htab(i-1)
    END DO
  END IF 

END SUBROUTINE build_Htab

!---------------------------------------------------------------------
!       HO1D_normalize
!               -normalizes the matrix elements in Hij
!               -creates list of normalization constants
!---------------------------------------------------------------------
! Variables
! Hij           : 2D real*8, matrix to be normalized
! N             : int, number of basis functions
! a             : real*8, alpha of basis functions 
! Ncon          : 1D real*8, list of normalization constants
! error         : bool, true on exit if problem

SUBROUTINE  HO1D_normalize(Hij,N,a,Ncon,error)
  IMPLICIT NONE
 
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Ncon
  REAL(KIND=8), INTENT(IN) :: a
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8) :: pi, const
  INTEGER(KIND=8) :: pwr, fct
  INTEGER :: i,j

  error = .FALSE.
  pi = 3.1415926535897932
  const = a/SQRT(pi)

  Ncon(0) = SQRT(1*const) 
  
  pwr = 1
  fct = 1
  DO i=1,N-1
    pwr = pwr * 2
    fct = i*fct 
    Ncon(i) = SQRT(1.0/(pwr*fct)*const)
  END DO

  DO j=0,N-1
    DO i=j,N-1
      Hij(i,j) = Hij(i,j)*Ncon(i)*Ncon(j)
    END DO
  END DO

  WRITE(*,*) "N00", Ncon(0)*Ncon(0)
  WRITE(*,*) "N10", Ncon(1)*Ncon(0)
  WRITE(*,*) "N11", Ncon(1)*Ncon(1)

END SUBROUTINE HO1D_normalize

!---------------------------------------------------------------------
!       HO1D_kinetic
!               -adds in (normalized) kinetic energy to Hij
!               - this abuses the viral theorum, and the fact
!               - that our basis functions are orthogonal HO
!---------------------------------------------------------------------
! Variables
! Hij           : 2D real*8, matrix of integrals
! N             : int, number of basis functions
! k             : real*8, k of basis functions
! m             : real*8, m of basis functions
! a             : real*8, alpha of basis functions
! Ncon          : 1D real*8, list of normalization constants
! error         : bool, true on exit if problem

SUBROUTINE HO1D_kinetic(Hij,N,k,m,a,Ncon,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Ncon
  REAL(KIND=8), INTENT(IN) :: k,m,a
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N
  REAL(KIND=8) :: w ,pi,C
  INTEGER :: i

  WRITE(*,*) "Calculating Kinetic Integrals"
  error = .FALSE.
  pi = 3.1415926535897932
  w = SQRT(k/m)
  C = SQRT(1/(2*w*m))

  DO i=0,N-1
    Hij(i,i) = Hij(i,i) + w*(1.0*i+0.5)/2.0 !using viral theorum
    IF (i+2 .LE. N-1) THEN
      Hij(i+2,i) = Hij(i+2,i) + (C**2.0*a**4.0*SQRT(n+2.)*SQRT(n+1.)&
                   - Ncon(i+2)*2*a**3.0*(i+2.)*C/Ncon(i+1))/(2.*m)
    END IF

    !+2 off diagonal terms... 
    ! I have been lazy and do not have analytical derivatives...
    !H(i+2,i) = H(i+2,j) + nint_fdif_2drv_HO(i,i+2,a,m,qmin,qmax,qeq)
  END DO

END SUBROUTINE HO1D_kinetic

!---------------------------------------------------------------------
!       HO1D_harmonic
!               -adds in harmonic energy terms to diagaonals 
!---------------------------------------------------------------------
! Variables
! Hij           : 2D real*8, matrix of integrals
! N             : int, number of basis functions
! k             : real*8, k of basis functions
! m             : real*8, m of basis functions
! a             : real*8, alpha of basis functions
! Ncon          : 1D real*8, list of normalization constants
! error         : bool, true on exit if problem

SUBROUTINE HO1D_harmonic(Hij,N,k,m,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), INTENT(IN) :: k,m
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N
  REAL(KIND=8) :: w
  INTEGER :: i

  error = .FALSE.
  w = SQRT(k/m)

  DO i=0,N-1
    Hij(i,i) = Hij(i,i) + w*(1.0*i+0.5) !using viral theorum
  END DO

END SUBROUTINE HO1D_harmonic

!---------------------------------------------------------------------

END MODULE ints_HO
