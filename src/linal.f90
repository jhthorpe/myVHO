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

SUBROUTINE diag(N,vmax,Hij,Ei,Cij,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Cij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Ei
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
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

  !WRITE(*,*) "The Hamiltonian"
  !DO i=0,N-1
  !  WRITE(*,*) Hij(i,0:N-1)
  !END DO
  !WRITE(*,*) 

  WRITE(*,*) 
  WRITE(*,*) "Diagonalizing the Hamiltonian"
  CALL DSYEV(JOBZ,UPLO,N,Hij(0:N-1,0:N-1),LDA,Ei(0:N-1),WORK(0:1),LWORK,INFO) 
  !CALL DSYEVX(JOBZ,RNGE,UPLO,N,Hij(0:N-1,0:N-1),LDA,dummy,dummy, &
  !            IL,UL,Ei(0:N-1), WORK(0:1),LWORK,INFO) 

  LWORK = CEILING(WORK(0))
  DEALLOCATE(WORK)
  ALLOCATE(WORK(0:LWORK-1))

  CALL DSYEV(JOBZ,UPLO,N,Hij(0:N-1,0:N-1),LDA,Ei(0:N-1),WORK(0:LWORK-1),LWORK,INFO) 
  WRITE(*,*) "Info is:", INFO 
  IF (INFO .NE. 0) THEN
    error = .TRUE.
    DEALLOCATE(WORK)
    DEALLOCATE(Ei)
    DEALLOCATE(Cij)
    RETURN
  END IF

  WRITE(*,*) 
  WRITE(*,*) "Vibrational Eigenvalues"
  WRITE(*,*) "Eigenvalues written to eigs.dat"
  OPEN(unit=103,file='eigs.dat',status='replace')
  DO j=0,MIN(vmax,N)-1
    WRITE(*,*) j, Ei(j)
    WRITE(103,*) j, Ei(j)
  END DO 
  DO j=MIN(vmax,N)+1,N-1
    WRITE(103,*) j, Ei(j)
  END DO
  CLOSE(unit=103)


  WRITE(*,*) 
  WRITE(*,*) "Eigenvectors written to evec.dat"
  OPEN(unit=102,file='evec.dat',status='replace')
  DO i=0,N-1
    WRITE(102,*) "Vibrational Eigenstate, eigenvalue : ", i, Ei(i)
    DO j=0,N-1
      WRITE(102,'(2x,I4,4x,F11.8)') j, Hij(i,j)
    END DO
    WRITE(102,*) "--------------------------"
  END DO
  CLOSE(unit=102)

  DEALLOCATE(WORK)
  
END SUBROUTINE diag
!---------------------------------------------------------------------
END MODULE linal
