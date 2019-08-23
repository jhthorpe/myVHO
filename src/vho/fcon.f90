!------------------------------------------------------------
! fcon  
!       - module for dealing with force constants
!------------------------------------------------------------
MODULE fcon
  USE input
  USE sort

CONTAINS
!------------------------------------------------------------
! fcon_get
!       - read in force constants
!       - naming scheme is as follows: aXY
!       - a = n,q,nothing. n -> number, q -> list of quant.
!             number, nothing -> actual values
!       - X = Q,P,QP,PQ
!       - Y = 1,2,3,4, order. 
!       - example: nquar -> number of Q^4 type force constants 
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nquad           : int, number of quadratic force constants
! qquad           : 1D int, QN of quadratic force constants
! quad            : 1D real*8, quadratic force constants
! ncubi           : int, number of cubic force constants
! qcubi           : 1D int, QN of cubic force constants
! cubi            : 1D real*8, cubic force constants
! nquar           : int, number of quartic force constants
! qquar           : 1D int, QN of quartic force constants
! quar            : 1D real*8, quartic force constants
! error         : int, exit code
!

SUBROUTINE fcon_get(job,ndim,nquad,qquad,quad,ncubi,qcubi,cubi,&
                 nquar,qquar,quar,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: quad,cubi,quar
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qquad,qcubi,qquar
  INTEGER, INTENT(INOUT) :: nquad,ncubi,nquar
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,job
  INTEGER :: voff
  LOGICAL :: ex
  error = 0
  IF (job .NE. 1 .AND. job .NE. 2 .AND. job .NE. 3) RETURN
  WRITE(*,*) "Reading force constants"
  INQUIRE(file='voff.in',EXIST=ex)
  IF (ex) THEN
    OPEN(file='voff.in',unit=100,status='old')
    READ(100,*) voff
    CLOSE(unit=100)
  ELSE
    voff = 0
  END IF
  WRITE(*,*) "Numbering offset is:", voff
  WRITE(*,*)
 
  CALL fcon_read_quad(ndim,voff,nquad,qquad,quad,error)
  IF (error .NE. 0) RETURN
  CALL fcon_read_cubi(job,ndim,voff,ncubi,qcubi,cubi,error)
  IF (error .NE. 0) RETURN
  CALL fcon_read_quar(job,ndim,voff,nquar,qquar,quar,error)
  IF (error .NE. 0) RETURN
END SUBROUTINE fcon_get

!------------------------------------------------------------
! fcon_read_Q1
!       - read info for force constants of type q^1
!------------------------------------------------------------
! ndim          : int, number of dimensions
! voff          : int, numbering offset
! nQ1           : int, number of q^1 force con
! qQ1           : 1D int, quantum numbers
! Q1            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_Q1(ndim,voff,nQ1,qQ1,Q1,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Q1
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qQ1
  INTEGER, INTENT(INOUT) :: nQ1,error
  INTEGER, INTENT(IN) :: ndim,voff
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
      IF (j-voff .LT. 1 .OR. j-voff .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In Q1.in, input",i,", is outside range [1:ndim]" 
        WRITE(*,*) "Are you sure 'voff.in' is correct?"
        error = 1
      ELSE
        qQ1(i) = j-voff-1
      END IF
    END DO
    CLOSE(unit=fid)
  ELSE
    nQ1 = 0
  END IF

END SUBROUTINE fcon_read_Q1

!------------------------------------------------------------
! fcon_read_quad
!       - read quadratic force constants from cfour 
!         'quadratic' file
!       - note that quadratic terms are diagonal in normal
!         coordinates
!------------------------------------------------------------
! ndim          : int, number of dimensions
! voff          : int, numbering offset
! nquad           : int, number of q^2 force con
! qquad           : 1D int, quantum numbers
! quad            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_quad(ndim,voff,nquad,qquad,quad,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: quad
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qquad
  INTEGER, INTENT(INOUT) :: nquad,error
  INTEGER, INTENT(IN) :: ndim,voff
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Rtemp
  INTEGER, DIMENSION(:), ALLOCATABLE :: Itemp
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  INTEGER :: i,j,k,fid,fline
  LOGICAL :: ex
  error = 0
  fname = 'quadratic'
  fid = 402 
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN

    CALL input_fline(fline,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read_quad  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(Itemp(0:fline-1))
    ALLOCATE(Rtemp(0:fline-1))
    Itemp = -1
    Rtemp = 0.0D0 
    nquad = 0

    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,fline-1
     
      !read and sort
      READ(fid,*) j,val
      j = j - voff
      IF (j .LT. 1 .OR. j .GT. ndim) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In quadratic, input",i,", is outside range [1:ndim]" 
        WRITE(*,*) "Are you sure 'voff.in' is correct?"
        error = 1
      ELSE IF (ALL(j-1 .NE. Itemp(0:nquad-1))) THEN
        nquad = nquad + 1
        Itemp(nquad-1) = j-1
        Rtemp(nquad-1) = val
      END IF
    END DO
    CLOSE(unit=fid)
    
    !now, put them all in efficient order
    ALLOCATE(qquad(0:nquad-1))
    ALLOCATE(quad(0:nquad-1))
    DO i=0,nquad-1
      qquad(i) = Itemp(i)
      quad(i) = Rtemp(i)
    END DO

  ELSE
    nquad = 0
  END IF

  WRITE(*,*) "Quadratic Force Constants"
  DO i=0,nquad-1
    WRITE(*,'(1x,I3,4x,F24.15)') qquad(i)+voff+1,quad(i)
  END DO
  WRITE(*,*)

  IF (ALLOCATED(Itemp)) DEALLOCATE(Itemp)
  IF (ALLOCATED(Rtemp)) DEALLOCATE(Rtemp)

  IF (ndim .NE. nquad) THEN
    WRITE(*,*) "fcon_read_quad  : ERROR"
    WRITE(*,*) "You seem to be missing some quadratic force constants"
    WRITE(*,*) "Number of dimensions :", ndim
    WRITE(*,*) "Number of quadratic phi :", nquad
  END IF

END SUBROUTINE fcon_read_quad

!------------------------------------------------------------
! fcon_read_cubi
!       - read info for force constants of type q^3
!       - if the jobtype is 2, diagonal terms will be
!         neglected
!------------------------------------------------------------
! job           : int, jobtype
! ndim          : int, number of dimensions
! voff          : int, numbering offset
! ncubi           : int, number of q^3 force con
! qcubi           : 1D int, quantum numbers
! cubi            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_cubi(job,ndim,voff,ncubi,qcubi,cubi,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: cubi
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qcubi
  INTEGER, INTENT(INOUT) :: ncubi,error
  INTEGER, INTENT(IN) :: job,ndim,voff
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Rtemp
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Itemp
  INTEGER, DIMENSION(0:2) :: idx
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  INTEGER :: i,j,k,l,fid,fline
  LOGICAL :: ex,match
  error = 0
  fname = 'cubic'
  fid = 403 
  
  IF (job .EQ. 3) RETURN
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(fline,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read_cubi  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(Itemp(0:2,0:fline-1))
    ALLOCATE(Rtemp(0:fline-1))
    ncubi = 0

    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,fline-1

      !read and sort
      READ(fid,*) idx,val
      idx = idx - voff
      !skip diagonal terms for jobtype 2
      !IF (job .EQ. 2 .AND. idx(0) .EQ. idx(1) .AND.  &
      !    idx(1) .EQ. idx(2)) CYCLE 
      CALL sort_int_ijk(idx)
      IF (ANY(idx .LT. 1) .OR. ANY(idx .GT. ndim)) THEN
        WRITE(*,*) "fcon_read_cubi  : ERROR"
        WRITE(*,*) "In cubic, input",i,", is outside range [1:ndim]" 
        WRITE(*,*) "Are you sure 'voff.in' is correct?"
        error = 1
      END IF

      !check this is something we haven't seen before
      match = .FALSE.
      DO j=0,ncubi-1
        IF (ALL(idx-1 .EQ. Itemp(0:2,j))) match = .TRUE.
      END DO
      IF (.NOT. match) THEN 
        ncubi = ncubi + 1
        Itemp(0:2,ncubi-1) = idx-1
        Rtemp(ncubi-1) = val
      END IF

    END DO
    CLOSE(unit=fid)

    !put into order
    ALLOCATE(qcubi(0:3*ncubi-1))
    ALLOCATE(cubi(0:ncubi-1))
    DO i=0,ncubi-1
      qcubi(3*i) = Itemp(0,i)
      qcubi(3*i+1) = Itemp(1,i)
      qcubi(3*i+2) = Itemp(2,i)
      cubi(i) = Rtemp(i)
    END DO
  ELSE
    ncubi = 0
  END IF

  WRITE(*,*) "Cubic Force Constants" 
  DO i=0,ncubi-1
    WRITE(*,'(1x,3(I3,2x),4x,F24.15)') qcubi(3*i:3*i+2)+voff+1,cubi(i)
  END DO
  WRITE(*,*)

  IF (ALLOCATED(Itemp)) DEALLOCATE(Itemp)
  IF (ALLOCATED(Rtemp)) DEALLOCATE(Rtemp)

END SUBROUTINE fcon_read_cubi

!------------------------------------------------------------
! fcon_read_quar
!       - read info for force constants of type q^4
!       - if jobtype = 2, skip diagonal terms
!------------------------------------------------------------
! job           : int, jobtype
! ndim          : int, number of dimensions
! voff          : int, numbering offset
! nquar           : int, number of q^4 force con
! qquar           : 1D int, quantum numbers
! quar            : 1D real*8, force con
! error         : exit code

SUBROUTINE fcon_read_quar(job,ndim,voff,nquar,qquar,quar,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: quar
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qquar
  INTEGER, INTENT(INOUT) :: nquar,error
  INTEGER, INTENT(IN) :: job,ndim,voff
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Rtemp
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Itemp
  INTEGER, DIMENSION(0:3) :: idx
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  INTEGER :: i,j,k,l,fid,fline
  LOGICAL :: ex,match
  error = 0
  fname = 'quartic'
  fid = 404 
  IF (job .EQ. 3) RETURN
  
  INQUIRE(file=TRIM(fname),EXIST=ex)
  IF (ex) THEN
    CALL input_fline(fline,fname,error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "fcon_read_quar  : ERROR"
      WRITE(*,*) "There is some problem in ",TRIM(fname)
    END IF
    ALLOCATE(Itemp(0:3,0:fline-1))
    ALLOCATE(Rtemp(0:fline-1))
    nquar = 0

    OPEN(file=TRIM(fname),unit=fid,status='old')
    DO i=0,fline-1

      !read and sort
      READ(fid,*) idx,val
      idx = idx - voff
      !skip diagonal terms for jobtype 2
      !IF (job .EQ. 2 .AND. idx(0) .EQ. idx(1) .AND.  &
      !    idx(1) .EQ. idx(2) .AND. idx(2) .EQ. idx(3)) CYCLE 
      CALL sort_int_ijkl(idx)
      IF (ANY(idx .LT. 1) .OR. ANY(idx .GT. ndim)) THEN
        WRITE(*,*) "fcon_read  : ERROR"
        WRITE(*,*) "In quartic, input",i,", is outside range [1:ndim]" 
        WRITE(*,*) "Are you sure 'voff.in' is correct?"
        error = 1
      END IF

      !check this is something we haven't seen before
      match = .FALSE.
      DO j=0,nquar-1
        IF (ALL(idx-1 .EQ. Itemp(0:3,j))) match = .TRUE.
      END DO
      IF (.NOT. match) THEN 
        nquar = nquar + 1
        Itemp(0:3,nquar-1) = idx-1
        Rtemp(nquar-1) = val
      END IF

    END DO
    CLOSE(unit=fid)

    !put into order
    ALLOCATE(qquar(0:4*nquar-1))
    ALLOCATE(quar(0:nquar-1))
    DO i=0,nquar-1
      qquar(4*i) = Itemp(0,i)
      qquar(4*i+1) = Itemp(1,i)
      qquar(4*i+2) = Itemp(2,i)
      qquar(4*i+3) = Itemp(3,i)
      quar(i) = Rtemp(i)
    END DO
    
  ELSE
    nquar = 0
  END IF

  WRITE(*,*) "Quartic Force Constants" 
  DO i=0,nquar-1
    WRITE(*,'(1x,4(I3,2x),4x,F24.15)') qquar(4*i:4*i+3)+voff+1,quar(i)
  END DO
  WRITE(*,*)

  IF (ALLOCATED(Itemp)) DEALLOCATE(Itemp)
  IF (ALLOCATED(Rtemp)) DEALLOCATE(Rtemp)

END SUBROUTINE fcon_read_quar

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
  
  STOP 5
  
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

  STOP 6
  
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

  STOP 7
  
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

  STOP 8
  
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
