PROGRAM vscf
  USE input
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: Vq,q
  REAL(KIND=8) :: qmin,qmax
  INTEGER :: nbasis, npoints
  LOGICAL :: error

  WRITE(*,*) 
  WRITE(*,*) "Starting VSCF calcuation"
  WRITE(*,*) "James H. Thorpe"
  WRITE(*,*)

  CALL read_input(nbasis,Vq,q,qmin,qmax,npoints,error)

END PROGRAM vscf
