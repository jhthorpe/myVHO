PROGRAM vho
  USE input
  USE ints_HO 
  USE linal
  USE proc
  USE nints
  USE fit
  USE val
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Hij,Cij,coef
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: Vq,q,Ei,Ni,W
  REAL(KIND=8) :: qmin,qmax,qeq,k,m,V_off,a,conv
  INTEGER :: N, npoints,vmax,ord,func
  LOGICAL :: error

  WRITE(*,*) 
  WRITE(*,*) "Starting VHO calcuation"
  WRITE(*,*) "James H. Thorpe"
  WRITE(*,*)

  CALL read_input(N,vmax,Vq,q,qmin,qmax,qeq,npoints,k,m,V_off,a,&
                 func,conv,error)
  IF (error) THEN
    WRITE(*,*) 
    WRITE(*,*) "ERROR ERROR ERROR"
    WRITE(*,*) "There was an error reading input"
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    STOP 1 
  END IF

  CALL fit_surf(func,Vq,q,npoints,conv,ord,coef,a,m,N,qeq,error)
  IF (error) THEN
    WRITE(*,*) "ERROR ERROR ERROR"
    WRITE(*,*) "There was an error fitting the surface"
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF (ALLOCATED(coef)) DEALLOCATE(coef)
    STOP 2
  END IF

  CALL HO1D_integrals(func,N,Vq,q,qmin,qmax,qeq,&
                      npoints,k,m,V_off,a,Hij,Ni,error)
  IF (error) THEN 
    WRITE(*,*) 
    WRITE(*,*) "ERROR ERROR ERROR"
    WRITE(*,*) "There was an error processing integrals"
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF(ALLOCATED(Hij)) DEALLOCATE(Hij)
    IF(ALLOCATED(Ni)) DEALLOCATE(Ni)
    IF (ALLOCATED(coef)) DEALLOCATE(coef)
    STOP 3 
  END IF

  IF (func .EQ. -1) THEN
    IF(ALLOCATED(Hij)) DEALLOCATE(Hij)
    IF (ALLOCATED(Cij)) DEALLOCATE(Cij)
    IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
    IF (ALLOCATED(q)) DEALLOCATE(q)
    IF (ALLOCATED(Ei)) DEALLOCATE(Ei)
    IF (ALLOCATED(coef)) DEALLOCATE(coef)
    STOP 0
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
    IF (ALLOCATED(coef)) DEALLOCATE(coef)
    STOP 4 
  END IF

  CALL make_gnuplot(N,vmax,Vq,q,qmin,qmax,qeq,npoints,&
                    k,m,V_off,a,Ei,Hij,Ni,error)
  IF(ALLOCATED(Hij)) DEALLOCATE(Hij)
  IF (ALLOCATED(Cij)) DEALLOCATE(Cij)
  IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
  IF (ALLOCATED(q)) DEALLOCATE(q)
  IF (ALLOCATED(Ei)) DEALLOCATE(Ei)
  IF (ALLOCATED(coef)) DEALLOCATE(coef)
  IF (error) STOP 5

  WRITE(*,*) 

END PROGRAM vho
