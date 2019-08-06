!------------------------------------------------------------
! rota
!       - module containing subroutines concering rotational
!         terms
!------------------------------------------------------------
MODULE rota
  USE input

CONTAINS

!------------------------------------------------------------
! rota_get
!------------------------------------------------------------
! nrota         : int, number of rotational constants
! rota          : 1D real*8, rotational constants
! error         : int, error code 

SUBROUTINE rota_get(nrota,rota,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: rota
  INTEGER, INTENT(INOUT) :: nrota,error
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  INTEGER :: voff,i,j
  LOGICAL :: ex
  error = 0
  fname = 'rota'

  CALL input_nline(nrota,fname)
  ALLOCATE(rota(0:nrota-1))
  IF (nrota .GT. 0) THEN
    OPEN(file=TRIM(fname),unit=400,status='old')
    DO i=0,nrota-1
      READ(400,*) j,val
      IF (j .GT. 0 .AND. j .LE. 3) THEN
        rota(j-1) = val
      ELSE
        WRITE(*,*) "rota_get  : ERROR"
        WRITE(*,*) "At line",i,"in file rota"
        WRITE(*,*) "Poorly formatted rotational constant"
      END IF
    END DO
    CLOSE(unit=400)
  END IF
  
  IF (nrota .GT. 0) THEN
     WRITE(*,*) "Rotational Constants"
     DO i=0,nrota-1
       WRITE(*,'(2x,I1,4x,F24.15)') i+1,rota(i)
     END DO
     WRITE(*,*)
  END IF
END SUBROUTINE rota_get

!------------------------------------------------------------
! rota_eval
!       - evaluates the contribution of the rotational terms
!         to the Watson hamiltonian
!
!       - currently this is just at second order 
!------------------------------------------------------------
! nrota         : int, number of rotational constants
! rota          : 1D real*8, rotational constants
! val           : real*8, the value to add

SUBROUTINE rota_eval(nrota,rota,val)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota
  REAL(KIND=8), INTENT(INOUT) :: val
  INTEGER, INTENT(IN) :: nrota
  val = SUM(rota(0:nrota-1))
  val = -0.25D0*val
END SUBROUTINE rota_eval

!------------------------------------------------------------

END MODULE rota
!------------------------------------------------------------
