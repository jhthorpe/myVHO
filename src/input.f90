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
! Vq		: 1D real*8, 1D potential energy surface
! qmin		: real*8, minimum r
! qmax		: real*8, max r
! qeq		: real*8, equilibrium q
! k	: real*8, k value of basis functions
! m	: real*8, m balue of basis functions
! V_off		: real*8, basis potential offset below qeq
! a	: real*8, alpha value of basis functions 
! error		: bool, true if error

SUBROUTINE read_input(N,Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: Vq,q
  REAL(KIND=8), INTENT(INOUT) :: qmin,qmax,qeq,k,m,V_off,a
  INTEGER, INTENT(INOUT) :: N, npoints
  LOGICAL, INTENT(INOUT)  :: error
  
  CHARACTER(LEN=1024) :: fname,word 
  REAL(KIND=8) :: temp
  INTEGER :: dummy,i
  
  error = .FALSE.
  fname = "Vq"

  OPEN(file='input',unit=100,status='old')
  READ(100,*) word, N
  READ(100,*) word, qeq
  READ(100,*) word, k
  READ(100,*) word, m
  READ(100,*) word, V_off
  CLOSE(unit=100)
  a = m*SQRT(k/m)
  WRITE(*,*) "Number of basis functions : ", N
  WRITE(*,*) "Equilibrium position : ", qeq
  WRITE(*,*) "Basis function k :", k
  WRITE(*,*) "Basis function m :", m
  WRITE(*,*) "Basis function offset :", V_off

  CALL getfline(npoints,fname,error)

  ALLOCATE(Vq(0:npoints-1)) 
  ALLOCATE(q(0:npoints-1))

  OPEN(file="Vq",unit=101,status='old')
  !read first line
  READ(101,*) dummy, q(0), temp
  Vq(0) = temp - V_off
  qmin = q(0) 

  DO i=1,npoints-2
    READ(101,*) dummy, q(i), temp
    Vq(i) = temp - V_off 
  END DO

  !read last line
  READ(101,*) dummy, q(npoints-1), temp
  Vq(npoints-1) = temp - V_off
  qmax = q(npoints-1)
   
  CLOSE(unit=101)

  q = q - qeq

  DO i=0,npoints-1
  WRITE(*,*) q(i), Vq(i)
  END DO

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
