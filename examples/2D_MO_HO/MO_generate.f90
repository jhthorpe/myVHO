program generate
  real(kind=8) :: qmin,qmax,qstep,k,q,De,al
  integer :: i

  write(*,*) "k = "
  read(*,*) k
  write(*,*) "De = "
  read(*,*) De
  write(*,*) "qmin = "
  read(*,*) qmin
  write(*,*) "qmax = "
  read(*,*) qmax
  write(*,*) "qstep = "
  read(*,*) qstep

  al = SQRT(k/(2.0D0*De))
  q = qmin
  do while (q .LE. qmax)
    write(100,*) q, De + De*(EXP(-2.0D0*al*q) - 2.0D0*EXP(-1.0D0*al*q)) 
    q = q + qstep
  END DO 

  write(*,*) "output saved in fort.100"
  

end program generate
