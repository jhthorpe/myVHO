!------------------------------------------------------------
! fcon  
!       - module for dealing with force constants
!------------------------------------------------------------
MODULE fcon
  USE input

CONTAINS
!------------------------------------------------------------
! fcon_get
!       - read in force constants
!       - naming scheme is as follows: aXY
!       - a = n,q,nothing. n -> number, q -> list of quant.
!             number, nothing -> actual values
!       - X = Q,P,QP,PQ
!       - Y = 1,2,3,4, order. 
!       - example: nQ4 -> number of Q^4 type force constants 
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nXY           : int, number of constant x order y 
! qXY           : 1D int, quantum numbers of constant x order y 
! XY            : 1D real*8, force constants of x order y
! error         : int, exit code

SUBROUTINE fcon_get(ndim,nQ1,qQ1,Q1,nQ2,qQ2,Q2,nQ3,qQ3,Q3,&
                 nQ4,qQ4,Q4,nP1,qP1,P1,nP2,qP2,P2,nQP,qQP,QP,&
                 nPQ,qPQ,PQ,error)
  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Q1,Q2,Q3,Q4,&
                                                            P1,P2,QP,PQ
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qQ1,qQ2,qQ3,qQ4,&
                                                       qP1,qP2,qQP,qPQ
  INTEGER, INTENT(INOUT) :: nQ1,nQ2,nQ3,nQ4,nP1,nP2,nQP,nPQ
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim

  error = 0

  CALL fcon_read_Q1(ndim,nQ1,qQ1,Q1,error)
  IF (error .NE. 0) RETURN
  CALL fcon_read_Q2(ndim,nQ2,qQ2,Q2,error)
  IF (error .NE. 0) RETURN
  CALL fcon_read_Q3(ndim,nQ3,qQ3,Q3,error)
  IF (error .NE. 0) RETURN
  CALL fcon_read_Q4(ndim,nQ4,qQ4,Q4,error)
  IF (error .NE. 0) RETURN
  CALL fcon_read_P1(ndim,nP1,qP1,P1,error)
  IF (error .NE. 0) RETURN
  CALL fcon_read_P2(ndim,nP2,qP2,P2,error)
  IF (error .NE. 0) RETURN
  CALL fcon_read_QP(ndim,nQP,qQP,QP,error)
  IF (error .NE. 0) RETURN
  CALL fcon_read_PQ(ndim,nPQ,qPQ,PQ,error)
  IF (error .NE. 0) RETURN

END SUBROUTINE fcon_get

!------------------------------------------------------------
! fcon_read_Q1
!       - read info for force constants of type q^1
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nQ1           : int, number of q^1 force con
! qQ1           : 1D int, quantum numbers
! Q1            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_Q1(ndim,nQ1,qQ1,Q1,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Q1
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qQ1
  INTEGER, INTENT(INOUT) :: nQ1,error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  INTEGER :: i,j,fid
  LOGICAL :: ex
  error = 0
  fname = 'Q1.in'
  fid = 401 
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(nQ1,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(qQ1(0:nQ1-1))
    ALLOCATE(Q1(0:nQ1-1))
    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,nQ1-1
      READ(fid,*) j,Q1(i)
      IF (j .LT. 1 .OR. j .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In Q1.in, input",i,", is outside range [1:ndim]" 
        error = 1
      ELSE
        qQ1(i) = j-1
      END IF
    END DO
    CLOSE(unit=fid)
  ELSE
    nQ1 = 0
  END IF

END SUBROUTINE fcon_read_Q1

!------------------------------------------------------------
! fcon_read_Q2
!       - read info for force constants of type q^2
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nQ2           : int, number of q^2 force con
! qQ2           : 1D int, quantum numbers
! Q2            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_Q2(ndim,nQ2,qQ2,Q2,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Q2
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qQ2
  INTEGER, INTENT(INOUT) :: nQ2,error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  INTEGER :: i,j,k,fid
  LOGICAL :: ex
  error = 0
  fname = 'Q2.in'
  fid = 402 
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(nQ2,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(qQ2(0:2*nQ2-1))
    ALLOCATE(Q2(0:nQ2-1))
    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,nQ2-1
      READ(fid,*) j,k,Q2(i)
      IF (j .LT. 1 .OR. j .GT. ndim .OR. k .LT. 1 .OR. k .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In Q2.in, input",i,", is outside range [1:ndim]" 
        error = 1
      ELSE
        qQ2(2*i) = j-1
        qQ2(2*i+1) = k-1
      END IF
    END DO
    CLOSE(unit=fid)
  ELSE
    nQ2 = 0
  END IF

  !WRITE(*,*) "TESTING TESTING"
  !DO i=0,nQ2-1
  !  WRITE(*,*) qQ2(2*i),qQ2(2*i+1),Q2(i)
  !END DO

END SUBROUTINE fcon_read_Q2

!------------------------------------------------------------
! fcon_read_Q3
!       - read info for force constants of type q^3
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nQ3           : int, number of q^3 force con
! qQ3           : 1D int, quantum numbers
! Q3            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_Q3(ndim,nQ3,qQ3,Q3,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Q3
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qQ3
  INTEGER, INTENT(INOUT) :: nQ3,error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  INTEGER :: i,j,k,l,fid
  LOGICAL :: ex
  error = 0
  fname = 'Q3.in'
  fid = 403 
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(nQ3,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(qQ3(0:3*nQ3-1))
    ALLOCATE(Q3(0:nQ3-1))
    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,nQ3-1
      READ(fid,*) j,k,l,Q3(i)
      IF (j .LT. 1 .OR. j .GT. ndim .OR. k .LT. 1 .OR. k .GT. ndim &
          .OR. l .LT. 1 .OR. l .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In Q3.in, input",i,", is outside range [1:ndim]" 
        error = 1
      ELSE
        qQ3(3*i) = j-1
        qQ3(3*i+1) = k-1
        qQ3(3*i+2) = l-1
      END IF
    END DO
    CLOSE(unit=fid)
  ELSE
    nQ3 = 0
  END IF

  !WRITE(*,*) "TESTING TESTING"
  !DO i=0,nQ3-1
  !  WRITE(*,*) qQ3(3*i:3*i+2),Q3(i)
  !END DO

END SUBROUTINE fcon_read_Q3

!------------------------------------------------------------
! fcon_read_Q4
!       - read info for force constants of type q^4
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nQ4           : int, number of q^4 force con
! qQ4           : 1D int, quantum numbers
! Q4            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_Q4(ndim,nQ4,qQ4,Q4,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Q4
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qQ4
  INTEGER, INTENT(INOUT) :: nQ4,error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  INTEGER :: i,j,k,l,m,fid
  LOGICAL :: ex
  error = 0
  fname = 'Q4.in'
  fid = 404 
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(nQ4,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(qQ4(0:4*nQ4-1))
    ALLOCATE(Q4(0:nQ4-1))
    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,nQ4-1
      READ(fid,*) j,k,l,m,Q4(i)
      IF (j .LT. 1 .OR. j .GT. ndim .OR. k .LT. 1 .OR. k .GT. ndim &
          .OR. l .LT. 1 .OR. l .GT. ndim .OR. &
          m .LT. 1 .OR. m .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In Q4.in, input",i,", is outside range [1:ndim]" 
        error = 1
      ELSE
        qQ4(4*i) = j-1
        qQ4(4*i+1) = k-1
        qQ4(4*i+2) = l-1
        qQ4(4*i+3) = m-1
      END IF
    END DO
    CLOSE(unit=fid)
  ELSE
    nQ4 = 0
  END IF

  !WRITE(*,*) "TESTING TESTING"
  !DO i=0,nQ4-1
  !  WRITE(*,*) qQ4(4*i:4*i+3),Q4(i)
  !END DO

END SUBROUTINE fcon_read_Q4

!------------------------------------------------------------
! fcon_read_P1
!       - read info for force constants of type p^1
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nP1           : int, number of p^1 force con
! qP1           : 1D int, quantum numbers
! P1            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_P1(ndim,nP1,qP1,P1,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: P1
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qP1
  INTEGER, INTENT(INOUT) :: nP1,error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  INTEGER :: i,j,fid
  LOGICAL :: ex
  error = 0
  fname = 'P1.in'
  fid = 405
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(nP1,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(qP1(0:nP1-1))
    ALLOCATE(P1(0:nP1-1))
    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,nP1-1
      READ(fid,*) j,P1(i)
      IF (j .LT. 1 .OR. j .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In P1.in, input",i,", is outside range [1:ndim]" 
        error = 1
      ELSE
        qP1(i) = j-1
      END IF
    END DO
    CLOSE(unit=fid)
  ELSE
    nP1 = 0
  END IF

  !WRITE(*,*) "TESTING TESTING TESTING"
  !DO i=0,nP1-1
  !  WRITE(*,*) qP1(i),P1(i)
  !END DO

END SUBROUTINE fcon_read_P1

!------------------------------------------------------------
! fcon_read_P2
!       - read info for force constants of type p^2
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nP2           : int, number of p^2 force con
! qP2           : 1D int, quantum numbers
! P2            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_P2(ndim,nP2,qP2,P2,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: P2
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qP2
  INTEGER, INTENT(INOUT) :: nP2,error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  INTEGER :: i,j,k,fid
  LOGICAL :: ex
  error = 0
  fname = 'P2.in'
  fid = 406 
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(nP2,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(qP2(0:2*nP2-1))
    ALLOCATE(P2(0:nP2-1))
    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,nP2-1
      READ(fid,*) j,k,P2(i)
      IF (j .LT. 1 .OR. j .GT. ndim .OR. k .LT. 1 .OR. k .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In P2.in, input",i,", is outside range [1:ndim]" 
        error = 1
      ELSE
        qP2(2*i) = j-1
        qP2(2*i+1) = k-1
      END IF
    END DO
    CLOSE(unit=fid)
  ELSE
    nP2 = 0
  END IF

  !WRITE(*,*) "TESTING TESTING"
  !DO i=0,nP2-1
  !  WRITE(*,*) qP2(2*i),qP2(2*i+1),P2(i)
  !END DO

END SUBROUTINE fcon_read_P2

!------------------------------------------------------------
! fcon_read_QP
!       - read info for force constants of type QP
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nQP           : int, number of QP force con
! qQP           : 1D int, quantum numbers
! QP            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_QP(ndim,nQP,qQP,QP,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: QP
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qQP
  INTEGER, INTENT(INOUT) :: nQP,error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  INTEGER :: i,j,k,fid
  LOGICAL :: ex
  error = 0
  fname = 'QP.in'
  fid = 407 
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(nQP,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(qQP(0:2*nQP-1))
    ALLOCATE(QP(0:nQP-1))
    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,nQP-1
      READ(fid,*) j,k,QP(i)
      IF (j .LT. 1 .OR. j .GT. ndim .OR. k .LT. 1 .OR. k .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In QP.in, input",i,", is outside range [1:ndim]" 
        error = 1
      ELSE
        qQP(2*i) = j-1
        qQP(2*i+1) = k-1
      END IF
    END DO
    CLOSE(unit=fid)
  ELSE
    nQP = 0
  END IF

 ! WRITE(*,*) "TESTING TESTING"
 ! DO i=0,nQP-1
 !   WRITE(*,*) qQP(2*i),qQP(2*i+1),QP(i)
 ! END DO

END SUBROUTINE fcon_read_QP

!------------------------------------------------------------
! fcon_read_PQ
!       - read info for force constants of type PQ
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nPQ           : int, number of PQ force con
! qPQ           : 1D int, quantum numbers
! PQ            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_PQ(ndim,nPQ,qPQ,PQ,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: PQ
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qPQ
  INTEGER, INTENT(INOUT) :: nPQ,error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  INTEGER :: i,j,k,fid
  LOGICAL :: ex
  error = 0
  fname = 'PQ.in'
  fid = 408 
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(nPQ,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(qPQ(0:2*nPQ-1))
    ALLOCATE(PQ(0:nPQ-1))
    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,nPQ-1
      READ(fid,*) j,k,PQ(i)
      IF (j .LT. 1 .OR. j .GT. ndim .OR. k .LT. 1 .OR. k .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In QP.in, input",i,", is outside range [1:ndim]" 
        error = 1
      ELSE
        qPQ(2*i) = j-1
        qPQ(2*i+1) = k-1
      END IF
    END DO
    CLOSE(unit=fid)
  ELSE
    nPQ = 0
  END IF

  !WRITE(*,*) "TESTING TESTING"
  !DO i=0,nPQ-1
  !  WRITE(*,*) qPQ(2*i),qPQ(2*i+1),PQ(i)
  !END DO

END SUBROUTINE fcon_read_PQ
!------------------------------------------------------------
END MODULE fcon
!------------------------------------------------------------
