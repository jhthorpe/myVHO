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
! sort_int_ijklm
!	- sorts length 5 array in ascending order
!------------------------------------------------------------
SUBROUTINE sort_int_ijklm(A)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: A
  INTEGER, DIMENSION(0:3) :: B
  IF (A(0) .GE. A(1) .AND. A(0) .GE. A(2) .AND. &
      A(0) .GE. A(3) .AND. A(0) .GE. A(4)) THEN
    B = [A(1),A(2),A(3),A(4)]
    CALL sort_int_ijkl(B)
    A(4) = A(0) 
    A(0:3) = B
  ELSE IF (A(1) .GE. A(0) .AND. A(1) .GE. A(2) .AND. &
           A(1) .GE. A(3) .AND. A(1) .GE. A(4)) THEN
    B = [A(0),A(2),A(3),A(4)]
    CALL sort_int_ijkl(B)
    A(4) = A(1)
    A(0:3) = B
  ELSE IF (A(2) .GE. A(0) .AND. A(2) .GE. A(1) .AND. &
           A(2) .GE. A(3) .AND. A(2) .GE. A(4)) THEN
    B = [A(0),A(1),A(3),A(4)]
    CALL sort_int_ijkl(B)
    A(4) = A(2)
    A(0:3) = B
  ELSE IF (A(3) .GE. A(0) .AND. A(3) .GE. A(1) .AND. &
           A(3) .GE. A(2) .AND. A(3) .GE. A(4)) THEN
    B = [A(0),A(1),A(2),A(4)]
    CALL sort_int_ijkl(B)
    A(4) = A(3) 
    A(0:3) = B
  ELSE
    B = [A(0),A(1),A(2),A(3)]
    CALL sort_int_ijkl(B)
    A(0:3) = B
  END IF
END SUBROUTINE sort_int_ijklm

!------------------------------------------------------------
! sort_int_ijklmn
!	- sorts length 5 array in ascending order
!------------------------------------------------------------
SUBROUTINE sort_int_ijklmn(A)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: A
  INTEGER, DIMENSION(0:4) :: B
  IF (A(0) .GE. A(1) .AND. A(0) .GE. A(2) .AND. &
      A(0) .GE. A(3) .AND. A(0) .GE. A(4) .AND. &
      A(0) .GE. A(5)) THEN
    B = [A(1),A(2),A(3),A(4),A(5)]
    CALL sort_int_ijklm(B)
    A(5) = A(0) 
    A(0:4) = B
  ELSE IF (A(1) .GE. A(0) .AND. A(1) .GE. A(2) .AND. &
           A(1) .GE. A(3) .AND. A(1) .GE. A(4) .AND. &
           A(1) .GE. A(5)) THEN
    B = [A(0),A(2),A(3),A(4),A(5)]
    CALL sort_int_ijklm(B)
    A(5) = A(1)
    A(0:4) = B
  ELSE IF (A(2) .GE. A(0) .AND. A(2) .GE. A(1) .AND. &
           A(2) .GE. A(3) .AND. A(2) .GE. A(4) .AND. &
           A(2) .GE. A(5)) THEN
    B = [A(0),A(1),A(3),A(4),A(5)]
    CALL sort_int_ijklm(B)
    A(5) = A(2)
    A(0:4) = B
  ELSE IF (A(3) .GE. A(0) .AND. A(3) .GE. A(1) .AND. &
           A(3) .GE. A(2) .AND. A(3) .GE. A(4) .AND. &
           A(3) .GE. A(5)) THEN
    B = [A(0),A(1),A(2),A(4),A(5)]
    CALL sort_int_ijklm(B)
    A(5) = A(3) 
    A(0:4) = B
  ELSE IF (A(4) .GE. A(0) .AND. A(4) .GE. A(1) .AND. &
           A(4) .GE. A(2) .AND. A(4) .GE. A(3) .AND. &
           A(4) .GE. A(5)) THEN
    B = [A(0),A(1),A(2),A(3),A(5)]
    CALL sort_int_ijklm(B)
    A(5) = A(4) 
    A(0:4) = B
  ELSE
    B = [A(0),A(1),A(2),A(3),A(4)]
    CALL sort_int_ijklm(B)
    A(0:4) = B
  END IF
END SUBROUTINE sort_int_ijklmn

!------------------------------------------------------------
! sort_dirty_1Dint_1Dreal8
!       - sorts a 1D real8 array into order based on a
!         1D int array, and checks that it is the 
!         appropriate length
!       - we assume that this array is quite small, so that
!         this algorithm doesn't impede performance too much
!------------------------------------------------------------
! n             : int, number of elements that should be found
! m             : int, length of arrays passed in 
! Ivec          : 1D int, integer vector
! Rvec          : 1D real*8, real vector
! Svec          : 1D real*8, real sorted vector
! error         : int, error code
SUBROUTINE sort_dirty_1Dint_1Dreal8(n,m,Ivec,Rvec,Svec,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Svec
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Rvec
  INTEGER, DIMENSION(0:), INTENT(IN) :: Ivec
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: n,m
  INTEGER :: off,loc
  INTEGER :: i,j
  error = 0
  IF (n .NE. m) THEN
    WRITE(*,*) "sort_dirty_1Dint_1Dreal8  : ERROR"
    WRITE(*,*) "Requested length,",n,"and given length",m,"do not match"
    error = 1
    RETURN
  END IF
  Svec = 0.0D0
  off = MINVAL(Ivec) 
  DO i=0,n-1
    loc = -1
    DO j=0,m-1
      IF (Ivec(j) .EQ. i) THEN
        loc = j
        EXIT
      END IF 
    END DO
    IF (loc .LT. 0) THEN
      WRITE(*,*) "sort_dirty_1Dint_1Dreal8  : ERROR"
      WRITE(*,*) "An element in Ivec,Rvec was missing"
      error = 1
      RETURN
    ELSE
      Svec(i) = Rvec(loc)
    END IF
  END DO 
END SUBROUTINE sort_dirty_1Dint_1Dreal8

!------------------------------------------------------------
END MODULE sort
!------------------------------------------------------------
