!------------------------------------------------------------
! input
!       - module for parsing input
!------------------------------------------------------------
MODULE input

CONTAINS

!------------------------------------------------------------
! input_fline
!------------------------------------------------------------
SUBROUTINE input_fline(fline,fname,error)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(IN) :: fname
  INTEGER, INTENT(INOUT) :: fline,error
  !Internal
  INTEGER :: io
  LOGICAL :: ex
  error = 0
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You need to create the input file : ", TRIM(fname)
    error = 1
    fline = -1
    RETURN
  END IF
  fline = 0
  io = 0
  OPEN(unit=999,file=TRIM(fname),status='old',access='sequential')
  DO WHILE (io .EQ. 0)
    READ(999,*,iostat=io)
    IF (io .EQ. 0) fline = fline + 1
  END DO
  CLOSE(unit=999)
END SUBROUTINE input_fline

!------------------------------------------------------------
! input_nline
!       - checks how many lines are in a file 
!       - if file does not exist, returns 0
!------------------------------------------------------------

SUBROUTINE input_nline(nline,fname)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(IN) :: fname
  INTEGER, INTENT(INOUT) :: nline
  !Internal
  INTEGER :: io
  LOGICAL :: ex
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    nline = 0
    RETURN
  END IF
  nline = 0
  io = 0
  OPEN(unit=999,file=TRIM(fname),status='old',access='sequential')
  DO WHILE (io .EQ. 0)
    READ(999,*,iostat=io)
    IF (io .EQ. 0) nline = nline + 1
  END DO
  CLOSE(unit=999)
END SUBROUTINE input_nline
!------------------------------------------------------------

END MODULE input
!------------------------------------------------------------
