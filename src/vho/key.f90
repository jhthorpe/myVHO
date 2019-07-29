!------------------------------------------------------------
! key
!       - module containing subroutines for unwrapping
!         nested loops into one large loop
!
!       - "ids" are the indices of the nested loops
!
!       - "idx" is the index of the unwrapped loop
!------------------------------------------------------------
MODULE key

CONTAINS 
!------------------------------------------------------------
! key_generate
!       - generate the proper key for determining quantum
!         numbers
!------------------------------------------------------------
! ndim          : int, number of dimensions (loops)
! nelem         : 1D int, nubmer of elemens in each loop 
! key           : 1D int, key

SUBROUTINE key_generate(ndim,nelem,key)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: key
  INTEGER, DIMENSION(0:), INTENT(IN) :: nelem
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i
  key(ndim-1) = 1
  DO i=ndim-2,0,-1
    key(i) = key(i+1)*nelem(i+1)
  END DO
END SUBROUTINE key_generate 

!------------------------------------------------------------
! key_idx2ids 
!       - returns loop ids array given an idx in the loop 
!       - this assumes the nested loops have the final index
!         as the inner most loop
!------------------------------------------------------------
! ndim          : int, number of dimensions
! idx           : int, location in 1D product array
! nelem         : 1D int, basis functons per dimension
! key           : 1D int, key vector
! ids           : 1D int, quantum number vector

SUBROUTINE key_idx2ids(ndim,idx,nelem,key,ids)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: ids
  INTEGER, DIMENSION(0:), INTENT(IN) :: nelem,key
  INTEGER, INTENT(IN) :: ndim,idx
  INTEGER :: i
  DO i=0,ndim-1
    ids(i) = MOD((idx)/key(i),nelem(i))
  END DO
END SUBROUTINE key_idx2ids 

!------------------------------------------------------------
! key_ids2idx
!       - given loop ids array, returns the unwrapped index 
!       - this assumes nested loops have the final index
!         as the innermost loop
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nelem          : 1D int, number of basis functions per dimension
! key           : 1D int, key vector
! ids           : 1D int, quantum number vector
! idx           : int, index to return

SUBROUTINE key_ids2idx(ndim,nelem,key,ids,idx)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: nelem,key,ids
  INTEGER, INTENT(INOUT) :: idx
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i
  idx = 0
  DO i=0,ndim-1
    idx = idx + key(i)*ids(i)
  END DO 
END SUBROUTINE key_ids2idx

!------------------------------------------------------------

END MODULE key
!------------------------------------------------------------
