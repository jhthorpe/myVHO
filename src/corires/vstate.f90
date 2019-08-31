!------------------------------------------------------------
! vstate        
!       - module containing subroutines for dealing with 
!         vibrational states
!------------------------------------------------------------
MODULE vstate

CONTAINS

!------------------------------------------------------------
! vstate_set
!       - sets the nubmer and ids of states to be 
!         calculated
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! h2l           : 1D int, converts labels -> harmonics
! nstates       : int, number of states to calculate
! states        : 2D int, list of vib quantum number states
!                         (quantum numbers, state)
SUBROUTINE vstate_set(nvib,h2l,nstates,states)
  IMPLICIT NONE
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: states
  INTEGER, DIMENSION(0:), INTENT(IN) :: h2l
  INTEGER, INTENT(INOUT) :: nstates
  INTEGER, INTENT(IN) :: nvib
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: temp
  INTEGER :: i,j,k,nread
  LOGICAL :: match

  IF (ALLOCATED(states)) DEALLOCATE(states)
  IF (ALLOCATED(temp)) DEALLOCATE(temp)

  READ(*,*) nread
  DO WHILE(nread .LE. 0) 
    WRITE(*,*) "nstates must be positive"
    READ(*,*) nread
  END DO
  nstates = nread

  ALLOCATE(temp(0:nvib-1,0:nread-1))

  !read in states
  DO i=0,nread-1
    READ(*,*) temp(0:nvib-1,i)
    DO k=0,i-1
      IF (ALL(temp(0:nvib-1,k) .EQ. temp(0:nvib-1,i))) THEN
        nstates = nstates - 1 
        EXIT
      END IF
    END DO
  END DO
  ALLOCATE(states(0:nvib-1,0:nstates-1))

  !save unique states
  j = 0
  DO i=0,nread-1
    match = .FALSE.
    DO k=0,i-1
      IF (ALL(temp(0:nvib-1,k) .EQ. temp(0:nvib-1,i))) THEN
        match = .TRUE. 
        EXIT
      END IF
    END DO
    IF (match) CYCLE
    states(0:nvib-1,j) = temp(0:nvib-1,i)
    j = j + 1 
  END DO

  WRITE(*,*) "Using States..."
  DO i=0,nstates-1
    WRITE(*,'(1x,999(I2x,1x))') states(0:nvib-1,i)
  END DO
  DEALLOCATE(temp)
END SUBROUTINE vstate_set

!------------------------------------------------------------
! vstate_level
!       - sets up calcluations of all states of a certain
!         level
!------------------------------------------------------------
SUBROUTINE vstates_level(nvib,h2l,nstates,states)
  IMPLICIT NONE
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: states
  INTEGER, DIMENSION(0:), INTENT(IN) :: h2l
  INTEGER, INTENT(INOUT) :: nstates
  INTEGER, INTENT(IN) :: nvib
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: temp
  INTEGER :: i,j,level

  IF (ALLOCATED(states)) DEALLOCATE(states)
  IF (ALLOCATED(temp)) DEALLOCATE(temp)

  READ(*,*) level
  DO WHILE(level .LE. 0) 
    WRITE(*,*) "level must be positive"
    READ(*,*) level
  END DO
  nstates = nvib

  ALLOCATE(states(0:nvib-1,0:nstates-1))
  DO j=0,nstates-1
    states(0:nvib-1,j) = 0
    states(j,j) = level 
  END DO

  WRITE(*,*) "Using States..."
  DO i=0,nstates-1
    WRITE(*,'(1x,999(I2x,1x))') states(0:nvib-1,i)
  END DO

END SUBROUTINE vstates_level
!------------------------------------------------------------

END MODULE vstate
!------------------------------------------------------------
