PROGRAM vscf
  USE input
  USE ints_HO 
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Hij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: Vq,q
  REAL(KIND=8) :: qmin,qmax,qeq
  INTEGER :: N, npoints
  LOGICAL :: error

  WRITE(*,*) 
  WRITE(*,*) "Starting VSCF calcuation"
  WRITE(*,*) "James H. Thorpe"
  WRITE(*,*)

  CALL read_input(N,Vq,q,qmin,qmax,qeq,npoints,error)
  CALL HO1D_integrals(N,Vq,q,qmin,qmax,qeq,npoints,Hij,error)

END PROGRAM vscf
