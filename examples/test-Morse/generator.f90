PROGRAM generator 
  IMPLICIT NONE
  
  REAL(KIND=8) :: r,req,k,rmin,rmax,dr,a,ke,De
  INTEGER :: i,nsteps

  req = 0.5
  rmin = 0.05
  rmax = 3
  nsteps = 50000
  De = 10.0
  ke = 50.0

  a = SQRT(ke/(2*De))
  dr = (rmax - rmin)/nsteps

  r = rmin
  DO i=0,nsteps
    WRITE(*,*) i, r, De*(EXP(-2*a*(r-req))-2*EXP(-1.0*a*(r-req)))
    r = r + dr
  END DO 
  

END PROGRAM generator
