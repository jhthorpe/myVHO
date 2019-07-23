!------------------------------------------------------------
! ints
!       - module supporting integral calculations
!------------------------------------------------------------
MODULE ints

CONTAINS
!------------------------------------------------------------
! ints_key
!       - generate the proper key for determining quantum
!         numbers
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, nubmer of basis functions per d
! key           : 1D int, key
! error         : int, exit code

SUBROUTINE ints_key(ndim,nbas,key,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: key
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i
  error = 0
  key(ndim-1) = 1
  DO i=ndim-2,0,-1
    key(i) = key(i+1)*nbas(i+1)
  END DO
END SUBROUTINE ints_key

!------------------------------------------------------------
! ints_qnum
!       - returns quantum number array given the proper key 
!       - this assumes the product array has the final index
!         as the inner most loop
!------------------------------------------------------------
! ndim          : int, number of dimensions
! it            : int, location in 1D product array
! nbas          : 1D int, basis functons per dimension
! key           : 1D int, key vector
! qnum          : 1D int, quantum number vector
! error         : int, exti code

SUBROUTINE ints_qnum(ndim,it,nbas,key,qnum,error)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: qnum
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,key
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,it
  INTEGER :: i
  error = 0
  DO i=0,ndim-1
    qnum(i) = MOD((it)/key(i),nbas(i))
  END DO 
END SUBROUTINE ints_qnum
!------------------------------------------------------------

END MODULE ints
!------------------------------------------------------------
