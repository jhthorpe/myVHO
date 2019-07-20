!//////////////////////////////////////////////////////////////////
!///
!///		Module Containing subroutines for reading
!///		input for myVHO
!///
!///		James H. Thorpe
!///
!//////////////////////////////////////////////////////////////////

MODULE input
  USE fit
  USE val
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
! k     	: real*8, k value of basis functions
! m     	: real*8, m balue of basis functions
! Voff		: real*8, basis potential offset below qeq
! a     	: real*8, alpha value of basis functions 
! func           : int, type of functing function
! conv          : real*8, convergence for func
! units         : int, units (and coordinate system)
! error		: bool, true if error

SUBROUTINE read_input(N,vmax,Vq,q,qmin,qmax,qeq,np,k,m,Voff,a,func,conv,units,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: Vq,q
  REAL(KIND=8), INTENT(INOUT) :: qmin,qmax,qeq,k,m,Voff,a,conv
  INTEGER, INTENT(INOUT) :: N, np,vmax,func,units
  LOGICAL, INTENT(INOUT)  :: error
  
  CHARACTER(LEN=1024) :: fname,word 
  REAL(KIND=8) :: temp,infty,A2B,amu2me
  INTEGER :: dummy,i,ueq,exitval
  LOGICAL :: exists
  
  error = .FALSE.
  fname = "Vq"
  infty = HUGE(qmin)
  A2B = 1.88973
  !amu2me = 1836 
  amu2me = 1822.89
  
  !check input file exists
  INQUIRE(file='vho.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    WRITE(*,*) "You need to create the input file, 'vho.dat'"
    error = .TRUE.
    RETURN
  END IF

  !read in basic data
  OPEN(file='vho.dat',unit=100,status='old')
  READ(100,*) word, N
  READ(100,*) word, vmax
  READ(100,*) word, qeq
  READ(100,*) word, m
  READ(100,*) word, units
  READ(100,*) word, func
  READ(100,*) word, conv
  CLOSE(unit=100)
  CALL getfline(np,fname,error)
  
  ALLOCATE(Vq(0:np-1)) 
  ALLOCATE(q(0:np-1))

  !read in q and Vq
  OPEN(file="Vq",unit=101,status='old')
  DO i=0,np-1
    READ(101,*) dummy, q(i), Vq(i)
  END DO
  CLOSE(unit=101)

  !Convert units
  IF (units .EQ. 1) THEN
    WRITE(*,*) "Converting Å -> Bohrs"
    WRITE(*,*) "Converting μ(gfm) -> μ(me)" 
    q = q*A2B
    qeq = qeq*A2B
    m = m*amu2me
  ELSE IF (units .EQ. 0) THEN
    WRITE(*,*) "No units will be converted... be very careful"
  ELSE IF (units .EQ. 2) THEN
    WRITE(*,*) "Using Dimensionless Normal Coordinates"
    WRITE(*,*) "Converting Hartrees -> cm-1"
  ELSE
    WRITE(*,*) "That unit conversion not supported!"
    WRITE(*,*) "Possible unit options are:"
    WRITE(*,*) " 0  :  no conversions at all"
    WRITE(*,*) " 1  :  Å -> Bohr, gfm -> me "
    WRITE(*,*) " 2  :  DNC, Hartree -> cm-1 "
    error = .TRUE.
    RETURN
  END IF
  WRITE(*,*) 

  qmin = q(0)
  qmax = q(np-1)
   
  Voff = -1.0*MINVAL(Vq)
  Vq = Vq + Voff

  IF (units .NE. 2) THEN !anything other than normal coordinates
    !Get k
    CALL get_k(Vq(0:np-1),q(0:np-1),qeq,np,k,error)
    a = SQRT(m*SQRT(k/m))
    q = q - qeq
  ELSE !we are in normal coordinates
    m = 1.0D0
    a = 1.0D0
    qeq = 0.0D0
    Vq = Vq * 219474.63 !convert to cm-1
    Voff = Voff * 219474.63
    INQUIRE(file='k',EXIST=exists)
    IF (.NOT. exists) THEN
      WRITE(*,*) "You need to create the input file 'k'"
      error = .TRUE.
    ELSE
      OPEN(file='k',unit=100,status='old')
        READ(100,*) k
      CLOSE(unit=100) 
    END IF 
    IF (error) RETURN
  END IF

  !print output
  WRITE(*,*) "Number of basis functions :", N
  WRITE(*,*) "Equilibrium position      :", qeq
  WRITE(*,*) "Basis function k          :", k
  WRITE(*,*) "Basis function m          :", m
  WRITE(*,*) "Basis function offset     :", Voff
  WRITE(*,*) "Basis function alpha      :", a
  WRITE(*,*) "Fit function type         :", func
  WRITE(*,*) "Fit RMS convergence       :", conv

END SUBROUTINE read_input 

!---------------------------------------------------------------------
! get_k
!       - calculates k via cubic spline interpolation, or reads from
!         file if provided
!---------------------------------------------------------------------
! Variables
! Vq            : 1D real*8, points on surface
! q             : 1D real*8, distances on surface
! qeq           : real*8, value of equilibrium position
! np            : int, number of points
! k             : real*8, value of k 
! error         : bool, true on exit if error

SUBROUTINE get_k(Vq,q,qeq,np,k,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Vq,q
  REAL(KIND=8), INTENT(INOUT) :: k
  REAL(KIND=8), INTENT(IN) :: qeq
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: np

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y2
  REAL(KIND=8) :: yp1,ypn,yf,yb,ym,h
  LOGICAL :: exists

  INQUIRE(FILE='k',EXIST=exists)
  IF (exists) THEN
    OPEN(file='k',unit=105,status='old')
    READ(105,*) k
    CLOSE(unit=105)

  ELSE 

    !if k file not present, calculate from potential at qeq
    ALLOCATE(y2(0:np-1))
    yp1 = (Vq(2) - Vq(0))/(q(2)-q(0))
    ypn = (Vq(np-1) - Vq(np-3))/(q(np-1) - q(np-3))
    h = 1.0D-8

    CALL spline(q(1:np-2),Vq(1:np-2),np-2,yp1,ypn,y2(1:np-2),np-2)
    CALL splint(q(1:np-2),Vq(1:np-2),y2(1:np-2),np-2,qeq+h,yf,error)
    CALL splint(q(1:np-2),Vq(1:np-2),y2(1:np-2),np-2,qeq,ym,error)
    CALL splint(q(1:np-2),Vq(1:np-2),y2(1:np-2),np-2,qeq-h,yb,error)

    !from central differences, 2nd derivative
    k = (yf - 2*ym + yb)/(h)**2.0

    DEALLOCATE(y2)
    !check k values
    IF (ABS((yf-ym) - (ym-yb)) .GT. 1.0D-15) THEN
      WRITE(*,*) "WARNING : potentially bad k value"
      WRITE(*,*) "forwards and backwards differences not equal:"
      WRITE(*,*) yf-ym, ym-yb
      WRITE(*,*) 
    END IF
 
    CALL checkval(k,error)
    IF (error) THEN
      WRITE(*,*) "input:get_k -- error out of checkval"
      RETURN
    END IF
 
  END IF

END SUBROUTINE get_k
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

END MODULE input
