!//////////////////////////////////////////////////////////////////
!///
!///            Program to calculate transition dipole moments from
!///            output of myVHO
!///
!///            James H. Thorpe
!///
!//////////////////////////////////////////////////////////////////

! Variables
! FC_A          : 2D real*8, frank-condon data for A vibrator
! FC_B          : 2D real*8, frank-condon data for B vibrator
! dq            : real*8, stepsize
! nsteps        : int, number of steps
! vmax_A        : int, maximum quantum number of A vibrator
! vmax_B        : int, maximum quantum number of B vibrator

PROGRAM tdp
  IMPLICIT NONE

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: FC_A, FC_B
  REAL(KIND=8) :: dq,temp,q,qmin,qeq
  INTEGER :: vmax_A,vmax_B, nsteps
  INTEGER :: i,j,u,vA,vB
  LOGICAL :: error

  WRITE(*,*) 
  WRITE(*,*) "Calculating Transition Dipole Moments"
  WRITE(*,*) 

  !check files exist
  CALL check_files(error) 
  IF (error) RETURN

  !read files 
  OPEN(file='FC_param.dat',unit=100,status='old')
  READ(100,*) nsteps
  READ(100,*) vmax_A
  READ(100,*) dq
  READ(100,*) qmin
  READ(100,*) qeq
  CLOSE(unit=100)

  ALLOCATE(FC_A(0:nsteps-1,0:vmax_A))

  OPEN(file='FC_ints.dat',unit=102,form='unformatted',status='old')
  READ(102) FC_A(0:nsteps-1,0:vmax_A)
  CLOSE(unit=102)

  WRITE(*,*) "   A(vi)   A(vj)      |<vi|μ|vj>|^2"
  DO i=0,vmax_A
    ! i -> i-1
    !j = i-1 
    DO j=0,i-1
      temp = 0.0D0
      q = qmin - qeq
      DO u=0,nsteps-1
        temp = temp + FC_A(u,i)*(q - qeq)*FC_A(u,j)*dq
        q = q + dq
      END DO
      WRITE(*,'(2x,I4,4x,I4,4x,F11.8)') i, j, temp**2.0
    END DO
  
    ! i -> i forbidden
    WRITE(*,'(2x,I4,4x,I4,4x,F11.8)') i, i, 0.0

    ! i -> i+1
    !j = i+1
    DO j=i+1,vmax_A
      temp = 0.0D0
      q = qmin - qeq
      DO u=0,nsteps-1
        temp = temp + FC_A(u,i)*(q - qeq)*FC_A(u,j)*dq
        q = q + dq
      END DO
      WRITE(*,'(2x,I4,4x,I4,4x,F11.8)') i, j, temp**2.0
    END DO
    WRITE(*,*) "------------------------------"
  END DO

  DEALLOCATE(FC_A)

CONTAINS

!---------------------------------------------------------------------
!       check_files
!               -checks that we have all the files we need
!---------------------------------------------------------------------
! Variables
! error         : bool, true on exit if we're missing a file
SUBROUTINE check_files(error)
  IMPLICIT NONE
  LOGICAL, INTENT(INOUT) :: error
  LOGICAL :: exists
  error = .FALSE.
  INQUIRE(file='FC_ints.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    error = .TRUE.
    WRITE(*,*) "You are missing 'FC_ints.dat'"
  END IF
  INQUIRE(file='FC_param.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    error = .TRUE.
    WRITE(*,*) "You are missing 'FC_param.dat'"
  END IF
END SUBROUTINE check_files
!---------------------------------------------------------------------

END PROGRAM tdp
