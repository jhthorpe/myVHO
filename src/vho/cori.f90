!------------------------------------------------------------
! cori
!       - module containing subroutines concerning 
!         coriolis coupling
!------------------------------------------------------------
MODULE cori
  USE input
  USE sort
  
CONTAINS
!------------------------------------------------------------
! cori_get
!       - reads in, sorts, and stores the coriolis coupling
!         constants
!       - the qantum numbers of the constants are stored
!         [(i,j),rotational mode]
!------------------------------------------------------------
! ndim          : int, number of dimensions
! ncori         : 1D int, number of coriolis terms per
!                         rotational mode
! qcori         : 2D real*8, quantum numbers of coriolis 
!                            terms per rotational mode [qn,rotation]
! cori          : 2D real*8, coriolis zetas [zeta,rotation]
! error         : int, exit code

SUBROUTINE cori_get(ndim,ncori,qcori,cori,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: cori
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: qcori
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ncori 
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Rtemp
  INTEGER, DIMENSION(:,:,:),ALLOCATABLE :: Itemp 
  INTEGER,DIMENSION(0:1) :: idx,idx0
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  LOGICAL :: ex,match
  INTEGER :: fline,fid,voff
  INTEGER :: i,j,k,a

  error = 0
  fname = 'coriolis'
  fid = 300

  !If there is only one dimension, or no coriolis file
  INQUIRE(file=TRIM(fname),exist=ex)
  IF (ndim .EQ. 1 .OR. .NOT. ex) THEN
    ALLOCATE(ncori(0:2))
    ALLOCATE(qcori(0:0,0:2))
    ALLOCATE(cori(0:0,0:2))
    qcori = 0
    cori = 0.0D0
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
  ALLOCATE(Itemp(0:1,0:fline-1,0:2))
  ALLOCATE(Rtemp(0:fline-1,0:2))
  ALLOCATE(ncori(0:2))
  ncori = 0
  Itemp = -1
  Rtemp = 0.0D0

  !read in data
  OPEN(file=TRIM(fname),unit=fid,status='old')
  DO k=0,fline-1
    READ(fid,*) i,idx0,val 
    a = i - 1
    idx = idx0
    idx = idx - voff 

    !sort and check
    CALL sort_int_ij(idx)
    IF (ANY(idx .LT. 1) .OR. ANY(idx .GT. ndim)) THEN
      WRITE(*,*) "cori_get  : ERROR"
      WRITE(*,*) "In coriolis, input",k,", is outside range [1:ndim]"
      WRITE(*,*) "Are you sure 'voff.in' is correct?"
      error = 1
    END IF

    !check we haven't seen this before
    match = .FALSE. 
    DO j=0,ncori(a)-1
      IF (ALL(idx-1 .EQ. Itemp(0:1,j,a))) match = .TRUE. 
    END DO
    IF (.NOT. match) THEN
      IF (idx(0) .NE. idx0(0)) val = -1.0D0*val !if there was a swap
      ncori(a) = ncori(a) + 1
      Itemp(0:1,ncori(a)-1,a) = idx - 1
      Rtemp(ncori(a)-1,a) = val
    END IF

  END DO
  CLOSE(unit=fid)

  !put into order
  ALLOCATE(qcori(0:2*MAXVAL(ncori)-1,0:2))
  ALLOCATE(cori(0:MAXVAL(ncori)-1,0:2))
  qcori = 0
  cori = 0.0D0
  DO a=0,2
    DO i=0,ncori(a)-1
      qcori(2*i,a) = Itemp(0,i,a)
      qcori(2*i+1,a) = Itemp(1,i,a)
      cori(i,a) = Rtemp(i,a)
    END DO
  END DO

  !Print coriolis constants
  IF (MAXVAL(ncori) .GT. 0) THEN
    WRITE(*,*) "Coriolis Constants"
    DO a=0,2
      DO j=0,ncori(a)-1
        WRITE(*,'(2x,I1,2x,2(I4,2x),F24.15)') a,qcori(2*j:2*j+1,a),cori(j,a)
      END DO
    END DO
    WRITE(*,*) 
  END IF

  IF (ALLOCATED(Itemp)) DEALLOCATE(Itemp)
  IF (ALLOCATED(Rtemp)) DEALLOCATE(Rtemp)
   
END SUBROUTINE cori_get

!------------------------------------------------------------

END MODULE cori
!------------------------------------------------------------
