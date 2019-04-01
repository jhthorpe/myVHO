!//////////////////////////////////////////////////////////////////
!///	
!///		Module for linear algebra subroutines 
!///			James H. Thorpe
!///
!//////////////////////////////////////////////////////////////////

MODULE linal

CONTAINS


!---------------------------------------------------------------------
!	diag
!		-diagonalizes the vibrational hamiltonian
!---------------------------------------------------------------------
! Variables
! N             : int, number of harmonic oscillator basis functions
! vmax          : int, max vibrational quantum number
! Hij           : 2D real*8, Hamiltonian Matrix
! Ei		: 1D real*8, vibrational energy levels
! Cij		: 2D real*8, eigenvectors (coefficients)
! Voff          : real*8, potential energy offset
! error         : bool, true if error

SUBROUTINE diag(N,vmax,Hij,Ei,Cij,Voff,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Cij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Ei
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8) , INTENT(IN) :: Voff 
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,vmax

  INTEGER :: i,j

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: WORK
  REAL(KIND=8) , DIMENSION(0:1) :: dummy
  CHARACTER(LEN=1) :: JOBZ, UPLO,RNGE
  INTEGER :: LDA,LWORK,INFO,IL,UL

  error = .FALSE.
  
  JOBZ = 'V'
  RNGE='I'
  UPLO = 'L'
  LWORK = -1
  LDA = N
  IL = 0
  UL = vmax
  ALLOCATE(WORK(0:1))
  ALLOCATE(Ei(0:N-1))
  ALLOCATE(Cij(0:N-1,0:N-1))

  CALL DSYEV(JOBZ,UPLO,N,Hij(0:N-1,0:N-1),LDA,Ei(0:N-1),WORK(0:1),LWORK,INFO) 
  !CALL DSYEVX(JOBZ,RNGE,UPLO,N,Hij(0:N-1,0:N-1),LDA,dummy,dummy, &
  !            IL,UL,Ei(0:N-1), WORK(0:1),LWORK,INFO) 
  LWORK = WORK(0)
  DEALLOCATE(WORK)
  ALLOCATE(WORK(0:LWORK-1))
  CALL DSYEV(JOBZ,UPLO,N,Hij(0:N-1,0:N-1),LDA,Ei(0:N-1),WORK(0:1),LWORK,INFO) 

  WRITE(*,*) 
  WRITE(*,*) "Vibrational Eigenvalues"
  DO j=0,vmax-1
    WRITE(*,*) Ei(j) + Voff
  END DO 

  !WRITE(*,*) 
  !WRITE(*,*) "Eigenvectors"
  !DO i=0,N-1
  !  WRITE(*,*) Hij(i,0:vmax-1)
  !END DO

  DEALLOCATE(WORK)
  
END SUBROUTINE diag
!---------------------------------------------------------------------
END MODULE linal
