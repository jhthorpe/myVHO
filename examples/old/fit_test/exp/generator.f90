PROGRAM ho_pot
  IMPLICIT NONE
  
  REAL(KIND=8) :: r,req,k,rmin,rmax,dr,a,b
  INTEGER :: i,nsteps

  k = 1.0
  req = 0.0
  rmin = 1.0
  rmax = 3.0
  nsteps = 10
  dr = (rmax - rmin)/nsteps
  a = 1.0
  b = 1.0

  r = rmin
  DO i=0,nsteps
    WRITE(*,*) i, r, 1.0 + a*EXP(-b*r) 
    r = r + dr
  END DO 
  

END PROGRAM ho_pot
