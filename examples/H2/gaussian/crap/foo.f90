program foo
  real(kind=8) :: a

  a = 14.0D-124

  write(*,*) a
  write(*,'(F20.16)') a
  WRITE(*,'(ES22.15)') a
end program foo
