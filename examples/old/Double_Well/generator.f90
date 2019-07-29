PROGRAM generator 
  IMPLICIT NONE
  
  REAL(KIND=8) :: rr,rl,a,b,c
  REAL(KIND=8) :: r, rmin, rmax, dr
  INTEGER :: i,nsteps

  rmin = -10
  rmax = 10
  nsteps = 50000
  rr = 3.5
  rl = -3.5
  a = 1.0
  b = 5.0
  c = 1.0/(rr**2.0)

  dr = (rmax - rmin)/nsteps

  r = rmin
  DO i=0,nsteps
    WRITE(*,*) i, r, (1.0/2.0*c*(r-rr)**2.0 + 1.0/2.0*c*(r-rl)**2.0 &
                      + b*Exp(-a*r**2))
    r = r + dr
  END DO 
  

END PROGRAM generator
