!------------------------------------------------------------
! didq
!	- module containing subroutines for dealing with 
!	  derivatives of the intertia tensor wrt 
!	  coordinates
!------------------------------------------------------------
MODULE didq
  USE input
  USE sort

CONTAINS
!------------------------------------------------------------
! didq_get
!	- reads in the derivatives of the intertia tensor wrt 
!	  dimensionless normal coordinates
!
!	- these are stored like:
!         didq(m,3*i+j) -> a^{i,j}_m
!------------------------------------------------------------
! ndim		: int, number of dimensions
! ndidq		: 1D int, number of didq
! qdidq		: 2D int, quantum numbers of didq
! didq		: 2D real*8, values of didq
! error		: int, exit code

SUBROUTINE didq_get(ndim,ndidq,qdidq,didq,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: didq
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: qdidq
  INTEGER, DIMENSION(:), ALLOCATABLE :: ndidq
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: Rtemp
  INTEGER, DIMENSION(:,:,:),ALLOCATABLE :: Itemp
  INTEGER, DIMENSION(0:1) :: idx
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  LOGICAL :: ex,match
  INTEGER :: fline,fid,voff
  INTEGER :: i,j,k,m,a,b

  error = 0
  fname = 'didq'
  fid = 301

  !If there is no didq file 
  INQUIRE(file=TRIM(fname),exist=ex)
  IF (ndim .EQ. 1 .OR. .NOT. ex) THEN
    ALLOCATE(ndidq(0:8))
    ALLOCATE(qdidq(0:0,0:8))
    ALLOCATE(didq(0:0,0:8))
    ndidq = 0
    qdidq = 0
    didq = 0.0D0
    RETURN
  END IF

  !Read in vibrational offset 
  INQUIRE(file='voff.in',EXIST=ex)
  IF (ex) THEN
    OPEN(file='voff.in',unit=100,status='old')
    READ(100,*) voff
    CLOSE(unit=100)
  ELSE
    voff = 0
  END IF
 
  !Initialize
  CALL input_nline(fline,fname)
  ALLOCATE(Itemp(0:ndim-1,0:2,0:2))
  ALLOCATE(Rtemp(0:ndim-1,0:2,0:2))
  ALLOCATE(ndidq(0:8))
  ndidq = 0
  Itemp = -1
  Rtemp = 0.0D0

  !read in data
  OPEN(file=TRIM(fname),unit=fid,status='old')
  DO k=0,fline-1
    READ(fid,*) idx,i,val
    m = i - voff - 1
    idx = idx - 1
    
    !sort and check
    !CALL sort_int_ij(idx)
    a = idx(0)
    b = idx(1)
    IF (ANY(idx .LT. 0) .OR. ANY(idx .GT. 2) .OR. m .LT. 0 .OR. m .GT. ndim-1) THEN
      WRITE(*,*) "didq_get  : ERROR"
      WRITE(*,*) "In didq, input",k,", is outside range [1:ndim] or [1:3]"
      WRITE(*,*) "Are you sure 'voff.in' is correct?"
      error = 1
    END IF

    !check we haven't seen this before
    match = .FALSE.
    DO j=0,ndidq(3*a+b)-1
      IF (m .EQ. Itemp(j,a,b)) match = .TRUE.
    END DO
    IF (.NOT. match) THEN
      ndidq(3*a+b) = ndidq(3*a+b) + 1
      Itemp(ndidq(3*a+b)-1,a,b) = m 
      Rtemp(ndidq(3*a+b)-1,a,b) = val
    END IF

  END DO
  CLOSE(unit=fid)

  !Put into order
  ALLOCATE(qdidq(0:MAXVAL(ndidq)-1,0:8))
  ALLOCATE(didq(0:MAXVAL(ndidq)-1,0:8))
  qdidq = 0
  didq = 0.0D0
  DO a=0,2
    DO b=0,2
      j = 3*a+b
      DO i=0,ndidq(j)-1
        qdidq(i,j) = Itemp(i,a,b)
        didq(i,j) = Rtemp(i,a,b)
      END DO
    END DO
  END DO

  !Print
  IF (MAXVAL(ndidq) .GT. 0) THEN
    WRITE(*,*) "dI/dq Constants" 
    DO a=0,2
      DO b=0,2
        j = 3*a+b
        DO i=0,ndidq(j)-1      
          WRITE(*,'(2x,I1,5x,I1,2x,I4,2x,F24.15)') a+1,b+1,qdidq(i,j)+voff+1,didq(i,j)
        END DO
      END DO
    END DO
    WRITE(*,*)
  END IF

  IF (ALLOCATED(Itemp)) DEALLOCATE(Itemp)
  IF (ALLOCATED(Rtemp)) DEALLOCATE(Rtemp)

END SUBROUTINE didq_get

!------------------------------------------------------------
! didq_get_debug
!	- reads in the derivatives of the intertia tensor wrt 
!	  dimensionless normal coordinates
!
!	- these are stored like:
!         didq(m,a,b) -> a^{a,b}_m
!------------------------------------------------------------
! ndim		: int, number of dimensions
! didq		: 3D real*8, values of didq
! error		: int, exit code

SUBROUTINE didq_get_debug(ndim,didq,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: didq
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  LOGICAL :: ex
  INTEGER :: fline,fid,voff
  INTEGER :: k,m,a,b
  
  error = 0
  fname = 'didq'
  fid = 301

  ALLOCATE(didq(0:ndim-1,0:2,0:2))
  didq = 0.0D0

  !If there is no didq file
  INQUIRE(file=TRIM(fname),exist=ex)
  IF (ndim .EQ. 1 .OR. .NOT. ex) THEN
    RETURN
  END IF

  !Read in vibrational offset 
  INQUIRE(file='voff.in',EXIST=ex)
  IF (ex) THEN
    OPEN(file='voff.in',unit=100,status='old')
    READ(100,*) voff
    CLOSE(unit=100)
  ELSE
    voff = 0
  END IF

  CALL input_nline(fline,fname)
  OPEN(file=TRIM(fname),unit=fid,status='old')

  !Read in data
  DO k=0,fline-1
    READ(fid,*) a,b,m,val

    a = a - 1
    b = b - 1
    m = m - voff - 1

    IF (a .LT. 0 .OR. a .GT. 2 .OR. b .LT. 0 .OR. b .GT. 2 .OR. &
        m .LT. 0 .OR. m .GT. ndim-1) THEN
      WRITE(*,*) "didq_get_debug  : ERROR"
      WRITE(*,*) "In didq, line ",k,", has a bad value"
      error = 1
    END IF

    didq(m,a,b) = val
  END DO
  CLOSE(unit=fid)

  WRITE(*,*) "dI/dq Constants -- DEBUG"
  DO a=0,2
    DO b=0,2
      DO m=0,ndim-1
        IF (ABS(didq(m,a,b)) .GT. 1.0D-15) THEN
          WRITE(*,'(2x,I1,2x,2(I4,2x),F24.15)') a+1,b+1,m+voff+1,didq(m,a,b)
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE didq_get_debug

!------------------------------------------------------------
END MODULE didq
!------------------------------------------------------------
