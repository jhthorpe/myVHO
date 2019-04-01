!//////////////////////////////////////////////////////////////////
!///
!///		Module Containing subroutines for reading
!///		input for myVSCF
!///
!///		James H. Thorpe
!///
!//////////////////////////////////////////////////////////////////

MODULE input
  IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------
!	read_input
!		-reads input
!---------------------------------------------------------------------
! Variables
! N		: int, number of harmonic oscillator basis functions
! vmax          : int, maximum vibrational quantum number
! Vq		: 1D real*8, 1D potential energy surface
! qmin		: real*8, minimum r
! qmax		: real*8, max r
! qeq		: real*8, equilibrium q
! k	: real*8, k value of basis functions
! m	: real*8, m balue of basis functions
! Voff		: real*8, basis potential offset below qeq
! a	: real*8, alpha value of basis functions 
! error		: bool, true if error

SUBROUTINE read_input(N,vmax,Vq,q,qmin,qmax,qeq,npoints,k,m,Voff,a,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: Vq,q
  REAL(KIND=8), INTENT(INOUT) :: qmin,qmax,qeq,k,m,Voff,a
  INTEGER, INTENT(INOUT) :: N, npoints,vmax
  LOGICAL, INTENT(INOUT)  :: error
  
  CHARACTER(LEN=1024) :: fname,word 
  REAL(KIND=8) :: temp
  INTEGER :: dummy,i,ueq
  
  error = .FALSE.
  fname = "Vq"

  OPEN(file='input',unit=100,status='old')
  READ(100,*) word, N
  READ(100,*) word, vmax
  READ(100,*) word, qeq
  READ(100,*) word, m
  CLOSE(unit=100)
  CALL getfline(npoints,fname,error)
  

  ALLOCATE(Vq(0:npoints-1)) 
  ALLOCATE(q(0:npoints-1))

  OPEN(file="Vq",unit=101,status='old')
  !read first line
  DO i=0,npoints-1
    READ(101,*) dummy, q(i), Vq(i)
  END DO
  qmin = q(0)
  qmax = q(npoints-1)
  CLOSE(unit=101)
   
  Voff = -1.0*MINVAL(Vq)
  Vq = Vq + Voff

  DO i=0,npoints-1
    IF (q(i) .GT. qeq) THEN
      ueq = i-1
      EXIT
    END IF
  END DO 
  k = (Vq(ueq+1) - 2*Vq(ueq) + Vq(ueq - 1))/(q(ueq+1)-q(ueq))**2.0
  IF (ABS(q(ueq+1)-q(ueq) - (q(ueq)-q(ueq-1))) .GT. 1.0D-15) THEN
    WRITE(*,*) "WARNING : potentially bad k value"
    WRITE(*,*) "forwards and backwards differences not equal"
    WRITE(*,*) q(ueq+1)-q(ueq), q(ueq)-q(ueq-1)
    WRITE(*,*) 
  END IF
  

  q = q - qeq

  a = SQRT(m*SQRT(k/m))
  WRITE(*,*) "Number of basis functions :", N
  WRITE(*,*) "Equilibrium position      :", qeq
  WRITE(*,*) "Basis function k          :", k
  WRITE(*,*) "Basis function m          :", m
  WRITE(*,*) "Basis function offset     :", Voff
  WRITE(*,*) "Basis function alpha      :", a

  !DO i=0,npoints-1
  !WRITE(*,*) q(i), Vq(i)
  !END DO


END SUBROUTINE read_input 

!---------------------------------------------------------------------
! getfline
!       James H. Thorpe
!       - checks if a file exists, and if so, how many lines are in it
!---------------------------------------------------------------------

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

!---------------------------------------------------------------------


END MODULE input
