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
! Hij           : 2D real*8, Hamiltonian Matrix, on exit, coeffs
! Ei		: 1D real*8, vibrational energy levels
! Voff          : real*8, potential energy offset
! error         : bool, true if error

SUBROUTINE diag(N,vmax,Hij,Ei,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Ei
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,vmax

  INTEGER :: i,j

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: WORK
  REAL(KIND=8) , DIMENSION(0:1) :: dummy
  REAL(KIND=8) :: au2cm
  CHARACTER(LEN=1) :: JOBZ, UPLO,RNGE
  INTEGER :: LDA,LWORK,INFO,IL,UL

  error = .FALSE.
  au2cm = 2.1947E5
  
  JOBZ = 'V'
  RNGE='I'
  UPLO = 'L'
  LWORK = -1
  LDA = N
  IL = 0
  UL = vmax
  ALLOCATE(WORK(0:1))
  ALLOCATE(Ei(0:N-1))

  !WRITE(*,*) "The Hamiltonian"
  !DO i=0,N-1
  !  WRITE(*,*) Hij(i,0:N-1)
  !END DO
  !WRITE(*,*) 

  WRITE(*,*) 
  WRITE(*,*) "Diagonalizing the Hamiltonian"
  CALL DSYEV(JOBZ,UPLO,N,Hij(0:N-1,0:N-1),LDA,Ei(0:N-1),WORK(0:1),LWORK,INFO) 

  LWORK = CEILING(WORK(0))
  DEALLOCATE(WORK)
  ALLOCATE(WORK(0:LWORK-1))

  CALL DSYEV(JOBZ,UPLO,N,Hij(0:N-1,0:N-1),LDA,Ei(0:N-1),WORK(0:LWORK-1),LWORK,INFO) 
  WRITE(*,*) "Info is:", INFO 
  IF (INFO .NE. 0) THEN
    error = .TRUE.
    DEALLOCATE(WORK)
    DEALLOCATE(Ei)
    RETURN
  END IF

  !Standardize the eigenvectors
  DO j=0,N-1
    i = MAXLOC(ABS(Hij(0:N-1,j)),1)-1
    IF (Hij(i,j) .LT. 0) Hij(0:N-1,j) = -1.0D0*Hij(0:N-1,j)
  END DO

  WRITE(*,*) 
  WRITE(*,*) "Vibrational Eigenvalues (a.u, cm-1)"
  WRITE(*,*) "Eigenvalues written to eigs.dat"
  OPEN(unit=103,file='eigs.dat',status='replace')
  DO j=0,MIN(vmax+1,N)-1
    WRITE(*,*) j, Ei(j), Ei(j)*au2cm
    WRITE(103,*) j, Ei(j),Ei(j)*au2cm
  END DO 
  DO j=MIN(vmax+1,N)+1,N-1
    WRITE(103,*) j, Ei(j), Ei(j)*au2cm
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
! lsqr
!       - calculates least squares coefficients via LAPACK
!       - note, N > M
!---------------------------------------------------------------------
! Variables
! A             : 2d real*8, M basis functions evaluated at N abscissa
! b             : 1d real*8, function to fit evaluated at N abscissa
! np            : int, number of abscissa (rows of A)
! nc            : int, nubmer of basis functions/coefs (cols of A)
! coef          : 1d real*8, list of M coefficients
! error         : bool, true on exit if problem

SUBROUTINE lsqr(A,b,np,nc,coef,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: A
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: coef
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: b
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: np,nc

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Ax,bv
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: S,WORK
  INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK
  REAL(KIND=8) :: RCOND
  INTEGER :: NRHS,LDA,LDB,RANK,LWORK,LIWORK,INFO,N,M

  error = .FALSE.

  IF (np .LT. nc) THEN
    error = .TRUE.
    WRITE(*,*) "ERROR"
    WRITE(*,*) "linal:lsqr -- tried to solve least squares with cols > rows"
    RETURN
  END IF

  M = np
  N = nc
  NRHS = 1
  LDA = M
  LDB = MAX(N,M)
  RCOND = -1
  LWORK = -1
  LIWORK = 2
  INFO = 0

  ALLOCATE(Ax(0:LDA-1,0:N-1))
  ALLOCATE(bv(0:LDB-1,0))
  ALLOCATE(S(0:MIN(M,N)-1))
  ALLOCATE(WORK(0:1))
  ALLOCATE(IWORK(0:LIWORK-1))

  Ax = A
  bv(0:M-1,0) = b(0:np-1)

  !call stuff
  CALL DGELSD(M,N,NRHS,Ax(0:LDA-1,0:N-1),LDA,bv(0:LDB-1,0:NRHS-1),LDB,S,&
              RCOND,RANK,WORK(0:1),LWORK,IWORK(0:LIWORK-1),INFO)

  IF (INFO .NE. 0) THEN
    error = .TRUE.
    WRITE(*,*) "ERROR"
    WRITE(*,*) "linal:lsqr -- dgelsd exited with value ", INFO
    GOTO 11
  END IF

  LWORK = CEILING(WORK(0))
  LIWORK = IWORK(0)
  DEALLOCATE(WORK)
  DEALLOCATE(IWORK)
  ALLOCATE(WORK(0:LWORK-1))
  ALLOCATE(IWORK(0:LIWORK-1))

  !call stuff again
  CALL DGELSD(M,N,NRHS,Ax(0:LDA-1,0:N-1),LDA,bv(0:LDB-1,0:NRHS-1),LDB,S,&
              RCOND,RANK,WORK(0:LWORK-1),LWORK,IWORK(0:LIWORK-1),INFO)

  IF (INFO .NE. 0) THEN
    error = .TRUE.
    WRITE(*,*) "ERROR"
    WRITE(*,*) "linal:lsqr -- dgelsd exited with value ", INFO
    GOTO 11 
  END IF

  !process coefficients
  coef(0:nc-1) = bv(0:N-1,0) 

11  DEALLOCATE(Ax)
  DEALLOCATE(bv)
  DEALLOCATE(S)
  DEALLOCATE(WORK)
  DEALLOCATE(IWORK)

END SUBROUTINE lsqr

!---------------------------------------------------------------------

END MODULE linal
