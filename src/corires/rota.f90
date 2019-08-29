!------------------------------------------------------------
! rota
!       - module containing subroutines concering rotational
!         terms
!------------------------------------------------------------
MODULE rota
  USE input
  USE conv

CONTAINS

!------------------------------------------------------------
! rota_get
!------------------------------------------------------------
! nrota         : int, number of rotational constants
! rota          : 1D real*8, rotational constants
! error         : int, error code 

SUBROUTINE rota_get(nrota,rota,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: rota
  INTEGER, INTENT(INOUT) :: nrota,error
  CHARACTER(LEN=1024) :: fname,line
  REAL(KIND=8) :: val
  INTEGER :: voff,i,j
  LOGICAL :: ex
  error = 0
  fname = 'rota'
  nrota = 3 
  INQUIRE(FILE=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "rota_get  : ERROR"
    WRITE(*,*) "You are missing the 'rota' file"
    error = error + 1
    RETURN
  END IF
 
  OPEN(file=TRIM(fname),unit=400,status='old')
  DO i=0,2
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
      error = error + 1
    END IF
  END DO
  CLOSE(unit=400)
  IF (error .NE. 0) RETURN

   WRITE(*,*) "Rotational Constants (cm-1)"
   WRITE(*,'(2x,A1,4x,F24.15)') "X",rota(0)
   WRITE(*,'(2x,A1,4x,F24.15)') "Y",rota(1)
   WRITE(*,'(2x,A1,4x,F24.15)') "Z",rota(2)
   WRITE(*,*)
   WRITE(*,*) "Rotational Constants (MHz)"
   WRITE(*,'(2x,A1,4x,F24.8)') "X",conv_cm2MHz(rota(0))
   WRITE(*,'(2x,A1,4x,F24.8)') "Y",conv_cm2MHz(rota(1))
   WRITE(*,'(2x,A1,4x,F24.8)') "Z",conv_cm2MHz(rota(2))
END SUBROUTINE rota_get

!------------------------------------------------------------
END MODULE rota
!------------------------------------------------------------

