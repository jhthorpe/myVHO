!//////////////////////////////////////////////////////////////////
!///		
!///	Module for plotting output of model system	
!///		
!//////////////////////////////////////////////////////////////////

MODULE plots
  IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------
!	make_gnuplot
!		-makes gnuplot file for the model system
!---------------------------------------------------------------------
! Variables
! N             : int, number of harmonic oscillator basis functions
! vmax          : int, max vibrational quantum number
! Vq            : 1D real*8, 1D potential energy surface
! qmin          : real*8, minimum r
! qmax          : real*8, max r
! qeq           : real*8, equilibrium q
! k	      	: real*8, k value of basis functions
! m    		: real*8, m balue of basis functions
! V_off         : real*8, basis potential offset below qeq
! a		: real*8, alpha value of basis functions 
! Ei		: 1D real*8, eigenvalues
! Cij           : 2D real*8, coefficients of basis functions
! Ni            : 1D real*8, normalization constants 
! error         : bool, true if error
 
SUBROUTINE make_gnuplot(N,vmax,Vq,q,qmin,qmax,qeq,npoints, &
           k,m,V_off,a,Ei,Cij,Ni,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Cij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Vq,q,Ei,Ni
  REAL(KIND=8), INTENT(IN) :: qmin,qmax,qeq,k,m,V_off,a
  INTEGER, INTENT(IN) :: N,npoints,vmax
  LOGICAL, INTENT(INOUT) :: error

  INTEGER :: i 

  error = .FALSE.
  WRITE(*,*)
  WRITE(*,*) "Generating plot.dat and gnurun" 
  WRITE(*,*) "Visualize with:"
  WRITE(*,*) "        gnuplot < gnurun"
  WRITE(*,*) "        display plot.png"

  OPEN(unit=101,file='plot.dat',status='replace')
    DO i=0,npoints-1
      WRITE(101,*) q(i)+qeq,Vq(i)
    END DO 
  CLOSE(unit=101) 
 
  OPEN(unit=100,file="gnurun",status='replace')
  WRITE(100,*) "set terminal png"
  WRITE(100,*) "set output 'plot.png'"
  WRITE(100,*) "set xrange [", qmin, ":", qmax,"]"
  WRITE(100,*) "set yrange [", MINVAL(Vq),":",MAXVAL(Vq),"]"
  WRITE(100,*) "k = ", k
  WRITE(100,*) "xeq = ", qeq
  WRITE(100,*) "plot 0.5 * k * (x - xeq)**2 t 'basis', \"
  !WRITE(100,*) "plot 0.5 * k * (x - xeq)**2 t 'basis', \"
  !WRITE(100,*) "'Vq' u 2:3 t 'Vq' ,\"
  WRITE(100,*) "'plot.dat' u 1:2 t 'Vq',\"
  DO i=0,MIN(vmax+1,N)-2
    WRITE(100,*) Ei(i),"t ","'v=",i,"'",",\"
  END DO
  WRITE(100,*) Ei(MIN(vmax+1,N)-1), "t ","'v=", MIN(vmax+1,N)-1,"'"
  CLOSE(unit=100) 

  CALL plot_wave(N,vmax,qmin,qmax,qeq,a,Cij,Ni,error)

END SUBROUTINE make_gnuplot

!---------------------------------------------------------------------
!       plot_wave
!               -plots wavefunctions up to vmax
!---------------------------------------------------------------------
SUBROUTINE plot_wave(N,vmax,qmin,qmax,qeq,a,Cij,Ni,error)
! Variables
! N             : int, number of harmonic oscillator basis functions
! vmax          : int, max vibrational quantum number
! qmin          : real*8, minimum r
! qmax          : real*8, max r
! qeq           : real*8, equilibrium q
! a		: real*8, alpha value of basis functions 
! Cij           : 2D real*8, coefficients of basis functions
! Ni            : 1D real*8, normalization constants 
! FC            : 2D real*8, integrals for Frank-Condon specta
! error         : bool, true if error
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Cij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Ni
  REAL(KIND=8), INTENT(IN) :: qmin,qmax,qeq,a
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,vmax
  
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: FC
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Htab
  REAL(KIND=8) :: q,dq,temp
  INTEGER :: nsteps,u,foff
  INTEGER :: i,j

  error = .FALSE.
  WRITE(*,*)
  WRITE(*,*) "Wavefunction plots saved as 'wave_X.dat'"
  WRITE(*,*) "Frank-Condon data saved as 'FC_ints.dat'"
  WRITE(*,*) "Frank-Condon parameters save as 'FC_param.dat'"
  CALL EXECUTE_COMMAND_LINE('rm wave_[0-9]*.dat')
  
  nsteps = 1000
  dq = (qmax - qmin)/nsteps  
  q = qmin - qeq
  ALLOCATE(Htab(0:N-1))
  ALLOCATE(FC(0:nsteps-1,0:vmax))
  FC = 0

  foff = 200
  CALL open_wave(vmax,foff,error)
  DO u=0,nsteps-1
    CALL build_Htab(N,a*q,Htab(0:N-1)) 
    DO i=0,vmax
      temp = 0.0
      DO j=0,N-1
        temp = temp + 1.0/Ni(j)*Cij(i,j)*Htab(j)*EXP(-a**2.0*q**2.0/2) 
      END DO
      WRITE(foff+i,*) q+qeq,temp
      FC(u,i) = temp
    END DO    

    q = q + dq  
  END DO
  CALL close_wave(vmax,foff,error)

  OPEN(file='FC_ints.dat',unit=106,form='unformatted',status='replace')
  WRITE(106) FC(0:nsteps-1,0:vmax)
  CLOSE(unit=106)

  OPEN(file='FC_param.dat',unit=107,status='replace')
  WRITE(107,*) nsteps
  WRITE(107,*) vmax
  WRITE(107,*) dq
  CLOSE(unit=107)
 
  
  DEALLOCATE(Htab) 
  DEALLOCATE(FC)

END SUBROUTINE plot_wave

!---------------------------------------------------------------------
!       open_wave
!               -opens files for writing wavefunction plots
!---------------------------------------------------------------------
SUBROUTINE open_wave(vmax,foff,error)
  IMPLICIT NONE

  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: vmax,foff

  CHARACTER(LEN=1024) :: fname
  INTEGER :: i,fid

  error = .FALSE.
  
  DO i=0,vmax
    fid = foff + i
    CALL get_fname(i,fname,error)
    OPEN(file=TRIM(fname),unit=fid,status='replace')
  END DO 

END SUBROUTINE open_wave

!---------------------------------------------------------------------
!       close_wave
!               -closes files for writing wavefunction plots
!---------------------------------------------------------------------
SUBROUTINE close_wave(vmax,foff,error)
  IMPLICIT NONE

  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: vmax,foff

  INTEGER :: i,fid

  error = .FALSE.
  
  DO i=0,vmax
    fid = foff + i
    CLOSE(unit=fid)
  END DO 

END SUBROUTINE close_wave
!---------------------------------------------------------------------
! get_fname
!       James H. Thorpe
!       
!       - formats file name
!---------------------------------------------------------------------
! Variables
! id            : int, id
! fname         : int, filename
! error         : bool, true on exit if error

SUBROUTINE get_fname(id,fname,error)
  IMPLICIT NONE
  !Inout
  CHARACTER(LEN=1024), INTENT(INOUT) :: fname
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: id
  !Internal
  CHARACTER(LEN=1024) :: str_fmt

  error = .FALSE.

  IF (id .LT. 10) THEN
    str_fmt = "(A5,I1,A4)"
  ELSE IF (id .GE. 10 .AND. id .LT. 100) THEN
    str_fmt = "(A5,I2,A4)"
  ELSE IF (id .GE. 100 .AND. id .LT. 1000) THEN
    str_fmt = "(A5,I3,A4)"
  ELSE IF (id .GE. 1000 .AND. id .LT. 10000) THEN
    str_fmt = "(A5,I4,A4)"
  ELSE
    WRITE(*,*) "More cases needed in get_fname"
    error = .TRUE.
    RETURN
  END IF

  WRITE(fname,str_fmt) "wave_",id,".dat"

END SUBROUTINE get_fname
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!       build_Htab
!               -constructs list of Hermite polynomials at
!               a particular q
!---------------------------------------------------------------------
! Variables
! Htab          : 1D real*8, list of hermite polynomials eval. at q
! q             : real*8, value to evaluate Htab at
! N             : int, max Hermite polynomial to calc

SUBROUTINE build_Htab(N,q,Htab)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Htab
  REAL(KIND=8), INTENT(IN) :: q
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8) :: H0, H1
  INTEGER :: i

  IF (N .EQ. 1) THEN
    Htab(0) = 1.0D0
    RETURN
  ELSE IF (N .EQ. 2) THEN
    Htab(0) = 1.0D0
    Htab(1) = 2.0*q
    RETURN
  ELSE
    Htab(0) = 1.0D0
    Htab(1) = 2.0*q
    DO i=1,N-2
      Htab(i+1) = 2*q*Htab(i) - 2.*i*Htab(i-1)
    END DO
  END IF

END SUBROUTINE build_Htab

!---------------------------------------------------------------------

END MODULE plots
