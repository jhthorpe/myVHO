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
  CHARACTER(LEN=1024) :: fname,line
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
      !READ(400,*) j,val
      READ(400,*) j,line
      IF (j .GT. 0 .AND. j .LE. 3) THEN
        IF (TRIM(line) .EQ. "Infinity") THEN
          val = 0.0D0
        ELSE
          READ(line,'(F24.15)') val
        END IF
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
END MODULE rota
!------------------------------------------------------------
