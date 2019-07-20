!program that conceverts input from Å and a.u to borh and cm-1

PROGRAM conv
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Q,En
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: A2B,au2cm,m2mp,mass,m1,m2
  LOGICAL :: flag
  INTEGER :: nlines,i,dummy

  WRITE(*,*) "THIS PROGRAM SHOULD NOT BE TRUSTED"
  STOP
  
  A2B = 3.7794519772  
  au2cm = 2.1947E5 
  m2mp = 1836

  WRITE(*,*) "Transforming Vq from Å and a.u. to Bohr*mp^1/2 and cm-1"

  OPEN(file='mass',unit=100)
  READ(100,*) m1
  READ(100,*) m2
  CLOSE(unit=100)

  mass = (m1*m2)*m2mp/(m1+m2)

  fname = "Vq"
  call getfline(nlines,fname,flag) 
  IF (flag) RETURN

  ALLOCATE(Q(0:nlines-1))
  ALLOCATE(En(0:nlines-1))
 
  OPEN(file='Vq',unit=101,status='old')
  DO i=0,nlines-1
    READ(101,*) dummy, Q(i), En (i)
  END DO
  CLOSE(unit=101)

  Q = Q*A2B!i*SQRT(mass)
  !En = En*au2cm

  OPEN(file='Vq_new',unit=102,status='replace')
  DO i=0,nlines-1
    WRITE(102,*) i, Q(i), En(i)
  END DO 
  CLOSE(unit=102)

  DEALLOCATE(Q)
  DEALLOCATE(En)

  WRITE(*,*) "Vq has been written to Vq_new in the standard coordinates"
  WRITE(*,*) "Remember to set mass to 1. Eigenvalues will be in cm-1"

CONTAINS

SUBROUTINE getfline(fline,fname,flag)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(IN) :: fname
  INTEGER, INTENT(INOUT) :: fline
  LOGICAL, INTENT(INOUT) :: flag
  !Internal
  INTEGER :: io
  LOGICAL :: ex
  flag = .FALSE.
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You need to create the input file : ", TRIM(fname)
    flag = .TRUE.
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

END SUBROUTINE getfline

END PROGRAM conv
