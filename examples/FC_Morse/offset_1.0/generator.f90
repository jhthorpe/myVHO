PROGRAM generator 
  IMPLICIT NONE
  
  REAL(KIND=8) :: r,req,k,rmin,rmax,dr,a,ke,De
  INTEGER :: i,nsteps

  req = 1.0
  rmin = 0.1
  rmax = 4
  nsteps = 50000
  De = 100.0
  ke = 1200.0

  a = SQRT(ke/(2*De))
  dr = (rmax - rmin)/nsteps

  r = rmin
  DO i=0,nsteps
    WRITE(*,*) i, r, De*(EXP(-2*a*(r-req))-2*EXP(-1.0*a*(r-req)))
    r = r + dr
  END DO 
  

END PROGRAM generator
