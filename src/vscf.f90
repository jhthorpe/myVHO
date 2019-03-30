PROGRAM vscf
  USE input
  USE ints_HO 
  USE plots
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Hij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: Vq,q
  REAL(KIND=8) :: qmin,qmax,qeq,k,m,V_off,a
  INTEGER :: N, npoints
  LOGICAL :: error

  WRITE(*,*) 
  WRITE(*,*) "Starting VSCF calcuation"
  WRITE(*,*) "James H. Thorpe"
  WRITE(*,*)

  CALL read_input(N,Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,error)
  CALL HO1D_integrals(N,Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,Hij,error)
  CALL make_gnuplot(Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,error)

END PROGRAM vscf
