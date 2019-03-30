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
! Hij           : 2D real*8, Hamiltonian Matrix
! Ei		: 1D real*8, vibrational energy levels
! Cij		: 2D real*8, eigenvectors (coefficients)
! error         : bool, true if error

SUBROUTINE diag(N,Hij,Ei,Cij,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Cij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Ei
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: WORK
  CHARACTER(LEN=1) :: JOBZ, UPLO
  INTEGER :: LDA,LWORK,INFO

  error = .FALSE.
  
  JOBZ = 'V'
  UPLO = 'L'
  LWORK = -1
  LDA = N
  ALLOCATE(WORK(0:1))
  ALLOCATE(Ei(0:N-1))

  CALL DSYEV(JOBZ,UPLO,N,Hij(0:N-1,0:N-1),LDA,Ei(0:N-1),WORK(0:1),LWORK,INFO) 
  WRITE(*,*) WORK(0) 
  WRITE(*,*) WORK(1) 

  DEALLOCATE(WORK)
  
END SUBROUTINE diag
!---------------------------------------------------------------------
END MODULE linal
