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
! error         : bool, true if error
 
SUBROUTINE make_gnuplot(N,vmax,Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,Ei,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Vq,q,Ei
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
      WRITE(101,*) q(i),Vq(i)
    END DO 
  CLOSE(unit=101) 
 
  OPEN(unit=100,file="gnurun",status='replace')
  WRITE(100,*) "set terminal png"
  WRITE(100,*) "set output 'plot.png'"
  WRITE(100,*) "set xrange [", qmin, ":", qmax,"]"
  WRITE(100,*) "k = ", k
  WRITE(100,*) "plot 0.5 * k * x**2 t 'basis', \"
  WRITE(100,*) "'plot.dat' u 1:2 t 'Vq',\"
  DO i=0,vmax-2
    WRITE(100,*) Ei(i),"t ","'v=",i,"'",",\"
  END DO
  WRITE(100,*) Ei(vmax-1), "t ","'v=", vmax-1,"'"
  CLOSE(unit=100) 

END SUBROUTINE make_gnuplot

!---------------------------------------------------------------------

END MODULE plots
