PROGRAM vho
  USE input
  USE ints_HO 
  USE linal
  USE proc
  USE nints
  USE poly
  USE val
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Hij,Cij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: Vq,q,Ei,Ni,coef
  REAL(KIND=8) :: qmin,qmax,qeq,k,m,V_off,a
  INTEGER :: N, npoints,vmax,ord
  LOGICAL :: error

  WRITE(*,*) 
  WRITE(*,*) "Starting VSCF calcuation"
  WRITE(*,*) "James H. Thorpe"
  WRITE(*,*)

  CALL read_input(N,vmax,Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,error)
  IF (error) THEN
    WRITE(*,*) 
    WRITE(*,*) "ERROR ERROR ERROR"
    WRITE(*,*) "There was an error reading input"
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    STOP 1 
  END IF

  CALL poly_fit(Vq,q,npoints,1.0D-6,ord,coef,error)
  IF (error) THEN
    WRITE(*,*) "ERROR ERROR ERROR"
    WRITE(*,*) "There was an error fitting the surface"
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF (ALLOCATED(coef)) DEALLOCATE(coef)
  END IF

  CALL HO1D_integrals(N,Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,Hij,Ni,error)
  IF (error) THEN 
    WRITE(*,*) 
    WRITE(*,*) "ERROR ERROR ERROR"
    WRITE(*,*) "There was an error processing integrals"
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF(ALLOCATED(Hij)) DEALLOCATE(Hij)
    IF(ALLOCATED(Ni)) DEALLOCATE(Ni)
    STOP 2 
  END IF

  CALL diag(N,vmax,Hij,Ei,error)
  IF (error) THEN
    WRITE(*,*) 
    WRITE(*,*) "ERROR ERROR ERROR"
    WRITE(*,*) "There was an error diagonalizing matrix"
    IF(ALLOCATED(Hij)) DEALLOCATE(Hij)
    IF (ALLOCATED(Cij)) DEALLOCATE(Cij)
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF (ALLOCATED(Ei)) DEALLOCATE(Ei)
    STOP 3 
  END IF

  CALL make_gnuplot(N,vmax,Vq,q,qmin,qmax,qeq,npoints,&
                    k,m,V_off,a,Ei,Hij,Ni,error)
  IF(ALLOCATED(Hij)) DEALLOCATE(Hij)
  IF (ALLOCATED(Cij)) DEALLOCATE(Cij)
  IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
  IF (ALLOCATED(q)) DEALLOCATE(q)
  IF (ALLOCATED(Ei)) DEALLOCATE(Ei)
  IF (error) STOP 4

  WRITE(*,*) 

END PROGRAM vho