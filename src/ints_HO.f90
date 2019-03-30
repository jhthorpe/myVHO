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
! N        : int, number of harmonic oscillator basis functions
! Hij		: 2D real*8, Hamiltonian Matrix
! Vq            : 1D real*8, 1D potential energy surface
! qmin          : real*8, minimum r
! qmax          : real*8, max r
! qeq           : real*8, equilibriu q
! error         : bool, true if error

SUBROUTINE HO1D_integrals(N,Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Vq,q
  REAL(KIND=8), INTENT(IN) :: qmin,qmax,qeq,k,m,V_off,a
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N, npoints

  error = .FALSE.

  WRITE(*,*) 
  WRITE(*,*) "Calculating HO integrals"

  ALLOCATE(Hij(0:N-1,0:N-1))
  Hij = 0.0D0

  !CALL HO1D_kinetic()
  CALL HO1D_potential(Hij,N,Vq,q,npoints,k,m,V_off,a,error)
  !CALL HO1D_normalize()

  !TESTING TESTING TESTING
  WRITE(*,*) "Testing potential energy"
  !TESTING TESTING TESTING


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

  !REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Htemp
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Htab 
  REAL(KIND=8) :: dq
  INTEGER :: i,j,u

  error = .FALSE.
  WRITE(*,*) "Calculating potential energy numerical integrals"
 
  !ALLOCATE(Htemp(0:N-1,0:N-1))
  ALLOCATE(Htab(0:N-1))
  
  DO u=0,npoints-1
   
    !Construct Hermitian table
    CALL build_Htab(N,a*q(u),Htab(0:N-1))
    dq = q(u+1) - q(u)

    DO j=0,N-1
      DO i=j,N-1
        Hij(i,j) = Hij(i,j) + Htab(i)*Htab(j)*EXP(-a**2.0*q(u)**2.0)*Vq(u)*dq
      END DO
    END DO

  END DO 
 
  DEALLOCATE(Htab)
  !DEALLOCATE(Htemp) 

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

  IF (N .EQ. 0) THEN
    Htab(0) = 1.0D0
    RETURN
  ELSE IF (N .EQ. 1) THEN
    Htab(0) = 1.0D0
    Htab(1) = 2.0*q
    RETURN
  ELSE 
    Htab(0) = 1.0D0
    Htab(1) = 2.0*q
    DO i=1,N-2
      Htab(i+1) = 2*q*Htab(i) - 2*(i)*Htab(i-1)
    END DO
  END IF 

END SUBROUTINE build_Htab
!---------------------------------------------------------------------
END MODULE ints_HO
