!------------------------------------------------------------
! sort
!       - module containing sorting routines
!------------------------------------------------------------
MODULE sort

CONTAINS
!------------------------------------------------------------
! sort_int_ij
!       - sorts length 2 array in ascending order
!------------------------------------------------------------
SUBROUTINE sort_int_ij(A)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: A
  INTEGER :: k
  IF (A(0) .GT. A(1)) THEN
    k = A(1)
    A(1) = A(0)
    A(0) = k
  END IF
END SUBROUTINE sort_int_ij 

!------------------------------------------------------------
! sort_int_ijk
!       - sorts length 3 array in ascending order
!------------------------------------------------------------
SUBROUTINE sort_int_ijk(A)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: A
  INTEGER, DIMENSION(0:2) :: B
  IF (A(0) .GE. A(1) .AND. A(0) .GE. A(2)) THEN
    IF (A(1) .GE. A(2)) THEN
      B = [A(2), A(1),A(0)] 
      A = B
    ELSE
      B = [A(1),A(2),A(0)]
      A = B
    END IF
  ELSE IF (A(1) .GE. A(0) .AND. A(1) .GE. A(2)) THEN
    IF (A(0) .GE. A(2)) THEN
      B = [A(2),A(0),A(1)]
      A = B
    ELSE
      B = [A(0),A(2),A(1)]
      A = B
    END IF
  ELSE IF (A(0) .GE. A(1)) THEN
    B = [A(1),A(0),A(2)]
    A = B
  END IF
END SUBROUTINE sort_int_ijk
!------------------------------------------------------------
! sort_int_ijkl
!       - sorts length 4 array in ascending order
!------------------------------------------------------------
SUBROUTINE sort_int_ijkl(A)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: A
  INTEGER, DIMENSION(0:2) :: B
  IF (A(0) .GE. A(1) .AND. A(0) .GE. A(2) .AND. A(0) .GE. A(3)) THEN
    B = [A(1),A(2),A(3)]
    CALL sort_int_ijk(B)
    A(3) = A(0)
    A(0:2) = B
  ELSE IF (A(1) .GE. A(0) .AND. A(1) .GE. A(2) .AND. A(1) .GE. A(3)) THEN
    B = [A(0),A(2),A(3)]
    CALL sort_int_ijk(B)
    A(3) = A(1)
    A(0:2) = B
  ELSE IF (A(2) .GE. A(0) .AND. A(2) .GE. A(1) .AND. A(2) .GE. A(3)) THEN
    B = [A(0),A(1),A(3)]
    CALL sort_int_ijk(B)
    A(3) = A(2)
    A(0:2) = B
  ELSE 
    B = [A(0),A(1),A(2)]
    CALL sort_int_ijk(B)
    A(0:2) = B
  END IF
END SUBROUTINE sort_int_ijkl
!------------------------------------------------------------

END MODULE sort
!------------------------------------------------------------
