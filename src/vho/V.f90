!------------------------------------------------------------
! V
!       - module containing subroutines dealing with the
!         the potential energy of the dimensions
!------------------------------------------------------------
MODULE V
  USE fname
  USE input

CONTAINS
!------------------------------------------------------------
! V_get
!       - gets the potential energy of dimensions
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of dimensions
! nabs          : int, number of abscissa
! q             : 1D real*8, abscissa
! Vij           : 2D real*8, potential energy
! error         : int, exit code          

SUBROUTINE V_get(job,ndim,nabs,q,Vij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Vij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim,nabs

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Vtemp,qtemp
  INTEGER, DIMENSION(0:ndim-1) :: npot
  INTEGER :: i,j

  error = 0
  !read in potential energies
  IF (job .EQ. 1) THEN
    CALL V_read(job,ndim,npot,qtemp,Vtemp,error)
  END IF
  IF (error .NE. 0) RETURN 


END SUBROUTINE V_get

!------------------------------------------------------------
! V_read
!       - read in potential energies from Vx.in files
!------------------------------------------------------------
! job           : int, jobtype
! ndim          : int, number of dimensions
! npot          : 1D int, number of points per dimension
! Vtemp         : 2D real*8, potentials at points [potential, dimension]
! qtemp         : 2D real*8, normco at points [normco, dimension]
! error         : int, error code

SUBROUTINE V_read(job,ndim,npot,qtemp,Vtemp,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Vtemp,qtemp
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: npot
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,job

  CHARACTER(LEN=1024) :: fname
  INTEGER :: foff,fid 
  INTEGER :: i,j

  error = 0
  foff = 200
  IF (job .EQ. 1) THEN  !if we are not already in abscissa form
    !check files exist and get the number of lines
    DO i=0,ndim-1
      CALL fname_Vin(i+1,fname,error)
      IF (error .NE. 0) RETURN
      CALL input_fline(npot(i),fname,error) 
      IF (error .NE. 0) RETURN
    END DO

    ALLOCATE(Vtemp(0:MAXVAL(npot)-1,0:ndim-1))
    ALLOCATE(qtemp(0:MAXVAL(npot)-1,0:ndim-1))

    DO j=0,ndim-1
      fid = foff + j
      CALL fname_Vin(j+1,fname,error)
      IF (error .NE. 0) RETURN
      OPEN(file=TRIM(fname),unit=fid,status='old')
      DO i=0,npot(j)-1
        READ(fid,*) qtemp(i,j),Vtemp(i,j)
      END DO
      CLOSE(unit=fid)
    END DO

  ELSE !we are already in abscissa form
    WRITE(*,*) "ERROR"
    WRITE(*,*) "Vread  : not coded to deal with precalced abscissa yet"
    error = 1
    RETURN
  END IF

END SUBROUTINE V_read
!------------------------------------------------------------

END MODULE V

!------------------------------------------------------------

