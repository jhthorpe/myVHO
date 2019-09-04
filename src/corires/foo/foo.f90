program foo
  character(len=1), dimension(0:3) :: a
  integer :: n
 
  n = 4
  a = ['i','j','k','l']

contains

recursive function combo(b,n) result (c)
  character(len=1), dimension(0:n-1) :: b

end program foo
