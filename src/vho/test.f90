program test
  USE cori
  implicit none

  integer, dimension(0:5) :: ids
  integer, dimension(0:2) :: psi
  real(kind=8) :: val
  integer :: ndim

  ndim = 3
  val = 0.0D0
 
  write(*,*) "test called"

  write(*,*) "enter psi"
  read(*,*) psi
  write(*,*) "enter ints"
  read(*,*) ids
  ids = ids - 1
  write(*,*) "psi  :", psi
  write(*,*) "ints :", ids+1

  CALL cori_HO_O4_ints(ndim,psi,psi,ids,val) 

end program test
