!//////////////////////////////////////////////////////////////////
!///
!///            Program to calculate FC spectra from 
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

PROGRAM fcs
  IMPLICIT NONE

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: FC_A, FC_B
  REAL(KIND=8) :: dq,temp,q,qmin,qeq
  INTEGER :: vmax_A,vmax_B, nsteps
  INTEGER :: i,j,u,vA,vB
  LOGICAL :: error

  WRITE(*,*) 
  WRITE(*,*) "Calculating Frank-Condon Factors"
  WRITE(*,*) 

  !check files exist
  CALL check_files(error) 
  IF (error) RETURN

  !read files 
  OPEN(file='A_param.dat',unit=100,status='old')
  READ(100,*) nsteps
  READ(100,*) vmax_A
  READ(100,*) dq
  READ(100,*) qmin
  READ(100,*) qeq
  CLOSE(unit=100)
  OPEN(file='B_param.dat',unit=101,status='old')
  READ(101,*) 
  READ(101,*) vmax_B
  CLOSE(unit=101)

  ALLOCATE(FC_A(0:nsteps-1,0:vmax_A))
  ALLOCATE(FC_B(0:nsteps-1,0:vmax_B))

  OPEN(file='A_ints.dat',unit=102,form='unformatted',status='old')
  READ(102) FC_A(0:nsteps-1,0:vmax_A)
  CLOSE(unit=102)
  OPEN(file='B_ints.dat',unit=103,form='unformatted',status='old')
  READ(103) FC_B(0:nsteps-1,0:vmax_B)
  CLOSE(unit=103)

  !the actual integrals
  WRITE(*,*) "   vA      vB      FC(vA,vB)"
  DO vA=0,vmax_A
    DO vB=0,vmax_B
      temp = 0.0D0
      DO i=0,nsteps-1
        temp = temp + FC_A(i,vA)*FC_B(i,vB)*dq
      END DO
      WRITE(*,'(2x,I4,4x,I4,4x,F11.8)') vA, vB, temp**2.0
    END DO 
    WRITE(*,*) "------------------------------"
  END DO

  DEALLOCATE(FC_A)
  DEALLOCATE(FC_B)

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
  INQUIRE(file='A_ints.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    error = .TRUE.
    WRITE(*,*) "You are missing 'A_ints.dat'"
  END IF
  INQUIRE(file='B_ints.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    error = .TRUE.
    WRITE(*,*) "You are missing 'B_ints.dat'"
  END IF
  INQUIRE(file='A_param.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    error = .TRUE.
    WRITE(*,*) "You are missing 'A_param.dat'"
  END IF
  INQUIRE(file='B_param.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    error = .TRUE.
    WRITE(*,*) "You are missing 'B_param.dat'"
  END IF
END SUBROUTINE check_files
!---------------------------------------------------------------------

END PROGRAM fcs
