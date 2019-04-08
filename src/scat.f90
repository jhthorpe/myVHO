!//////////////////////////////////////////////////////////////////
!///
!///	        Program for numerically calculating scattering
!///            wavefunctions above a potential	
!///
!///		James H. Thorpe
!///
!//////////////////////////////////////////////////////////////////

PROGRAM scat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: Vq,q,U
  REAL(KIND=8) :: qmin,qmax,qeq,m,Voff,En
  INTEGER :: npoints,vmax
  LOGICAL :: error

  WRITE(*,*) 
  WRITE(*,*) "Numerically Solving Scattering Problem"
  WRITE(*,*) "James H. Thorpe"
  WRITE(*,*) 
  CALL read_input(Vq,q,qmin,qmax,qeq,npoints,Voff,m,En,error) 
  IF (error) THEN
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    RETURN
  END IF

  CALL scat_solve(U,Vq,q,En,m,qmin,qmax,npoints,error)
  IF (error) THEN
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF (ALLOCATED(U)) DEALLOCATE(U)
    RETURN
  END IF
  
  CALL scat_plot(U,Vq,q,qeq,npoints) 
  
  IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
  IF (ALLOCATED(q)) DEALLOCATE(q)
  IF (ALLOCATED(U)) DEALLOCATE(U)

CONTAINS

!---------------------------------------------------------------------
!	read_input
!		-reads input
!---------------------------------------------------------------------
! Variables
! Vq		: 1D real*8, 1D potential energy surface
! q             : 1D real*8, 1D list of x values for Vq
! qmin		: real*8, minimum r
! qmax		: real*8, max r
! qeq		: real*8, equilibrium q
! Voff		: real*8, basis potential offset below qeq
! m             : real*8, mass
! En            : real*8, scattering energy
! error		: bool, true if error

SUBROUTINE read_input(Vq,q,qmin,qmax,qeq,npoints,Voff,m,En,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)  :: Vq,q
  REAL(KIND=8), INTENT(INOUT) :: qmin,qmax,qeq,Voff,En,m
  INTEGER, INTENT(INOUT) :: npoints
  LOGICAL, INTENT(INOUT)  :: error
  
  CHARACTER(LEN=1024) :: fname,word 
  REAL(KIND=8) :: temp,infty,A2B,me2mp
  INTEGER :: dummy,i,ueq,exitval,units
  LOGICAL :: exists
  
  error = .FALSE.
  fname = "Vq"
  infty = HUGE(qmin)
  A2B = 3.7794519772
  me2mp = 1836 
  
  !check input file exists
  INQUIRE(file='scat.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    WRITE(*,*) "You need to create the input file, 'scat.dat'"
    error = .TRUE.
    RETURN
  END IF

  !read in basic data
  OPEN(file='scat.dat',unit=100,status='old')
  READ(100,*) word, En
  READ(100,*) word, qeq
  READ(100,*) word, m
  READ(100,*) word, units
  CLOSE(unit=100)
  CALL getfline(npoints,fname,error)
  
  ALLOCATE(Vq(0:npoints-1)) 
  ALLOCATE(q(0:npoints-1))

  !read in q and Vq
  OPEN(file="Vq",unit=101,status='old')
  DO i=0,npoints-1
    READ(101,*) dummy, q(i), Vq(i)
  END DO

  !Convert units
  IF (units .EQ. 1) THEN
    WRITE(*,*) "Converting Å -> Bohrs"
    WRITE(*,*) "Converting μ(gfm) -> μ(me)" 
    q = q*A2B
    qeq = qeq*A2B
    m = m*me2mp
  ELSE IF (units .EQ. 0) THEN
    WRITE(*,*) "No units will be converted... be very careful"
  ELSE
    WRITE(*,*) "That unit conversion not supported!"
    WRITE(*,*) "Possible unit options are:"
    WRITE(*,*) " 0  :  no conversions at all"
    WRITE(*,*) " 1  :  Å -> Bohr, gfm -> me"
    error = .TRUE.
    RETURN
  END IF
  WRITE(*,*) 

  qmin = q(0)
  qmax = q(npoints-1)
  CLOSE(unit=101)

  Voff = -1.0*MINVAL(Vq)
  Vq = Vq + Voff

  !center around qeq
  q = q - qeq

  !print output
  WRITE(*,*) "Energy above minima       :", En
  WRITE(*,*) "Equilibrium position      :", qeq
  WRITE(*,*) "Reduced mass              :", m
  WRITE(*,*) "Basis function offset     :", Voff

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
!       scat_solve
!               -numerically solves for scattering wavefunction
!---------------------------------------------------------------------
! Variables
! U             : 1D real*8, normalized wavefunction
! Vq		: 1D real*8, 1D potential energy surface
! q             : 1D real*8, 1D list of x values for Vq
! En            : real*8, scattering energy
! m             : real*8, mass
! qmin		: real*8, minimum r
! qmax		: real*8, max r
! error		: bool, true if error

SUBROUTINE scat_solve(U,Vq,q,En,m,qmin,qmax,npoints,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: U
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Vq,q
  REAL(KIND=8), INTENT(IN) :: qmin,qmax,m,En
  LOGICAL,INTENT(INOUT) :: error
  INTEGER,INTENT(IN) :: npoints

  REAL(KIND=8) :: temp,h,infty
  INTEGER :: i

  error = .FALSE.
  infty = HUGE(h)
 
  WRITE(*,*) 
  WRITE(*,*) "Numerically solving for scattering wavefunction"
  WRITE(*,*) "Remember, this algorithm assumes that..."
  WRITE(*,*) "  1. qmin and V(qmin) ensure |U(E,qmin)> = 0"
  WRITE(*,*) "  2. dq is sufficiently small" 
  WRITE(*,*) 

  ALLOCATE(U(0:npoints-1))

  U(0) = 0.0
  !U(1) = 0.001 !arbitrarily chosen
  U(1) = 1.0 !arbitrarily chosen
  temp = U(1)**2.0
  !WRITE(*,*) "i+1,h, En, Vq, U(i), U(i-1),U(i+1)"
  DO i=1,npoints-2
    h = q(i+1) - q(i)

    U(i+1) = (2 + h**2.0*(-2*m*(En - Vq(i))) )*U(i) - U(i-1) 
!    WRITE(*,*) i+1,h,h**2.0,En, Vq(i), U(i), U(i-1), U(i+1)

    !check for NaN and inf
    IF (U(i+1) .GT. infty) THEN
      WRITE(*,'(2x,A1,2x,I4,2x,A12)') "U",i+1," is infinite"
      error = .TRUE.
      EXIT
    ELSE IF (U(i+1) .NE. U(i+1)) THEN
      WRITE(*,'(2x,A1,2x,I4,2x,A7)') "U",i+1," is NaN" 
      error = .TRUE.
      EXIT
    END IF

   ! temp = temp + U(i+1)**2.0

   ! IF (temp .GT. infty) THEN
   !   WRITE(*,'(2x,A1,2x,I4,2x,A12)') "N",i+1," is infinite"
   !   error = .TRUE.
   !   EXIT
   ! ELSE IF (temp .NE. temp) THEN
   !   WRITE(*,'(2x,A1,2x,I4,2x,A7)') "N",i+1," is NaN" 
   !   error = .TRUE.
   !   EXIT
   ! END IF

  END DO

  IF (error) RETURN

  !WRITE(*,*) "Normalizing the wavefunction"
  !U = U/SQRT(temp) !I don't know how to feel about this normalization...


END SUBROUTINE scat_solve
!---------------------------------------------------------------------
!       scat_plot
!               -plot scattering wavefunction
!---------------------------------------------------------------------
! Variables
! U             : 1D real*8, normalized wavefunction
! Vq            : 1D real*8, potential 
! q             : 1D real*8, list of points
! qeq           : real*8, equilibrium position
! npoints       : real*8, number of points
SUBROUTINE scat_plot(U,Vq,q,qeq,npoints)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: U,Vq,q
  REAL(KIND=8), INTENT(IN) :: qeq
  INTEGER, INTENT(IN) :: npoints

  INTEGER :: i

  WRITE(*,*) 
  WRITE(*,*) "Gennerating scattering plot data: scat_plot.dat"
  WRITE(*,*) 

  OPEN(file='scat_plot.dat',unit=103,status='replace')
  DO i=0,npoints-1
    WRITE(103,*) q(i) + qeq, U(i), Vq(i)
  END DO
  CLOSE(unit=103)

END SUBROUTINE scat_plot
!---------------------------------------------------------------------

END PROGRAM scat
