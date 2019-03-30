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
! nbasis        : int, number of harmonic oscillator basis functions
! Vq            : 1D real*8, 1D potential energy surface
! qmin          : real*8, minimum r
! qmax          : real*8, max r
! qeq           : real*8, equilibrium q
! k	      	: real*8, k value of basis functions
! m    		: real*8, m balue of basis functions
! V_off         : real*8, basis potential offset below qeq
! a		:real*8, alpha value of basis functions 
! error         : bool, true if error
 
SUBROUTINE make_gnuplot(Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: Vq,q
  REAL(KIND=8) :: qmin,qmax,qeq,k,m,V_off,a
  INTEGER :: N, npoints
  LOGICAL :: error

  WRITE(*,*) "Generating datafile and gnuplot runner" 
 
  OPEN(unit=100,file="gnurun",status='replace')

  CLOSE(unit=100) 


END SUBROUTINE make_gnuplot

!---------------------------------------------------------------------

END MODULE plots
