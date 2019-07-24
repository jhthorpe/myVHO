program generate
  real(kind=8) :: qmin,qmax,qstep,k,q
  integer :: i

  write(*,*) "k = "
  read(*,*) k
  write(*,*) "qmin = "
  read(*,*) qmin
  write(*,*) "qmax = "
  read(*,*) qmax
  write(*,*) "qstep = "
  read(*,*) qstep

  q = qmin
  i = 0
  do while (q .LE. qmax)
    i = i+1
    write(100,*) i,q, 0.5*k*q**2.0D0
    q = q + qstep
  END DO 

  write(*,*) "output saved in fort.100"
  

end program generate
