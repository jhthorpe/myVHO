!------------------------------------------------------------
! V
!       - module containing subroutines dealing with the
!         the potential energy of the dimensions
!------------------------------------------------------------
MODULE V
  USE fname
  USE input
  USE fit
  USE valu
  USE key

CONTAINS
!------------------------------------------------------------
! V_get
!       - gets the potential energy of dimensions
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of dimensions
! nabs          : 1D int, number of abscissa
! Vij           : 2D real*8, potential energy [abscissa,dimension]
! Vq            : 1D real*8, potential at abscissa [abscissa]
! error         : int, exit code          

SUBROUTINE V_get(job,ndim,nabs,Vij,Vq,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Vij
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Vq
  INTEGER, DIMENSION(0:), INTENT(IN) :: nabs
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim

  INTEGER, DIMENSION(0:ndim-1) :: npot
  INTEGER :: i,j,N

  error = 0
  N = PRODUCT(nabs)
  !read in potential energies
  CALL V_read(job,ndim,nabs,Vij,Vq,error)
  IF (error .NE. 0) RETURN 
  
  !adjust potentials to be zero
  IF (job .EQ. 2) THEN 
    DO j=0,ndim-1
      Vij(0:nabs(j)-1,j) = Vij(0:nabs(j)-1,j) - MINVAL(Vij(0:nabs(j)-1,j))
    END DO
    Vij = Vij*219474.63 !convert hartrees to cm-1
  ELSE
    Vq(0:N-1) = Vq(0:N-1) - MINVAL(Vq(0:N-1))
    Vq(0:N-1) = Vq(0:N-1)*219474.63
  END IF

  !WRITE(*,*) "TESTING TESTING TESTING"
  !DO i=0,ndim-1
  !  WRITE(*,*) "Dimension", i
  !  WRITE(*,*) "------------------------"
  !  DO j=0,nabs(i)-1
  !    WRITE(*,*) j,Vij(j,i)
  !  END DO
  !END DO

END SUBROUTINE V_get

!------------------------------------------------------------
! V_read
!       - read in potential energies from Vx.in files
!------------------------------------------------------------
! job           : int, jobtype
! ndim          : int, number of dimensions
! nabs          : int, number of abscissa
! npot          : 1D int, number of points per dimension
! Vij           : 2D real*8, potentials at points [potential, dimension]
! Vq            : 1D real*8, potentials at abscissa [potential]
! error         : int, error code

SUBROUTINE V_read(job,ndim,nabs,Vij,Vq,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Vij
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Vq
  INTEGER, DIMENSION(0:), INTENT(IN) :: nabs
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,job

  INTEGER, DIMENSION(0:ndim-1) :: key,ids
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val1
  INTEGER :: foff,fid,dummy 
  LOGICAL :: ex
  INTEGER :: i,j,N

  error = 0
  foff = 200
  IF (job .EQ. 2) THEN

    DO j=0,ndim-1
      fid = foff + j
      CALL fname_Vin(j+1,fname,error)
      IF (error .NE. 0) RETURN
      OPEN(file=TRIM(fname),unit=fid,status='old')
      DO i=0,nabs(j)-1
        READ(fid,*) dummy, val1
        IF (dummy .LT. 1 .OR. dummy .GT. nabs(j)) THEN
          WRITE(*,*) "V_read  : error"
          WRITE(*,*) "Line :",i,"of file",TRIM(fname)
          WRITE(*,*) "Input is outside of range [1:nabs(j)]"
          error = 1
          EXIT 
        END IF
        dummy = dummy - 1
        Vij(dummy,j) = val1
      END DO
      CLOSE(unit=fid)
    END DO

  ELSE IF (job .EQ. 3) THEN
    N = PRODUCT(nabs)
    INQUIRE(file='V.in',exist=ex)
    IF (.NOT. ex) THEN
      WRITE(*,*) "V_read  : ERROR"
      WRITE(*,*) "You need to create the file V.in" 
      error = 1
      RETURN
    END IF
    OPEN(file='V.in',unit=300,status='old')
    DO i=0,N-1
      READ(300,*) dummy,val1
      IF (dummy .LT. 1 .OR. dummy .GT. N) THEN
        WRITE(*,*) "V_read  : error"
        WRITE(*,*) "Line :",i,"of file 'V.in'"
        WRITE(*,*) "Input is outside of range [1:N]"
        error = 1
        EXIT 
      END IF
      dummy = dummy-1
      Vq(dummy) = val1 
    END DO 
    CLOSE(unit=300) 

  ELSE !we are in some other form 
    WRITE(*,*) "V_read  : ERROR"
    WRITE(*,*) "Bad jobtype" 
    error = 1
    RETURN
  END IF

END SUBROUTINE V_read

!------------------------------------------------------------
! V_spline
!       - calculates potential at abscissa via cubic
!         spline interpolation
!       - the code looks strange because I am mixing the 
!         order of dimensions and abscissa
!       - the code is not terribly efficient, but these
!         arrays are small and this won't cause problems
!------------------------------------------------------------
! ndim          : int, nubmer of dimensions
! nabs          : int, nubmer of abscissa
! npot          : 1D int, points per potential
! q             : 1D real*8, abscissa
! qtemp         : 2D real*8, input q's              [npot,ndim]
! Vtemp         : 2D real*8, input V's              [npot,ndim]
! Vij           : 2D real*8, potential at abscissa  [ndim,nabs]
! error         : int, exit code

SUBROUTINE V_spline(ndim,nabs,npot,q,qtemp,Vtemp,Vij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Vij 
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: qtemp,Vtemp
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q
  INTEGER, DIMENSION(0:), INTENT(IN) :: npot
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nabs
 
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y2
  INTEGER, DIMENSION(:), ALLOCATABLE :: Vtype
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: yp1,ypn,beta,alpha,val,v0,v1,v2,vn0,vn1,vn2
  INTEGER :: np,fid,foff
  INTEGER :: i,j

  error = 0
  WRITE(*,*) "Performing cubic spline interpolation of surface"
  WRITE(*,*) "Sorry, this seems to be bugged and is disabled"
  STOP

  !Determine Vtype for each dimension
  ! type 1      : V -> +inf as q -> +/- 
  !               model +/- as a*x^b
  ! type 2      : V -> De as q -> +, V -> +inf as q-> -
  !               model + as alpha*exp(-beta*q),  
  !               model - as a*x^b
  ! type 3      : V -> De as q -> -, V -> +inf as q-> +
  !               model - as alpha*exp(-beta*q),  
  !               model + as a*x^b
  ALLOCATE(Vtype(0:ndim-1))
  Vtype = -1
  DO j=0,ndim-1
    np = npot(j)
    v0 = Vtemp(0,j) 
    v1 = Vtemp(1,j)
    v2 = Vtemp(2,j)
    vn0 = Vtemp(np-1,j)
    vn1 = Vtemp(np-2,j)
    vn2 = Vtemp(np-3,j)
    IF (v0 - v1 .GT. v1 - v2 ) THEN
      IF( vn0 - vn1 .GT. vn1 - vn2) THEN
        Vtype(j) = 1
      ELSE
        Vtype(j) = 2
      END IF
    ELSE 
      IF( vn0 - vn1 .GT. vn1 - vn2) THEN
        Vtype(j) = 3
      ELSE
        WRITE(*,*) "There is a case I haven't predicted" 
        error = 1
      END IF
    END IF
  END DO
  WRITE(*,*) "Vtype is:", Vtype
  IF (error .NE. 0) RETURN

  !WRITE(*,*) "TESTING TESTING TESTING" 
  !WRITE(*,*) "qtemp0 is", qtemp(:,0)
  !WRITE(*,*) "qtemp1 is", qtemp(:,1)
 
  ALLOCATE(y2(0:MAXVAL(npot)-1))
  !generate y2 for each dimension 
  DO j=0,ndim-1
    np = npot(j)
    yp1 = 0.5D0*( Vtemp(1,j) - Vtemp(0,j) )/( qtemp(1,j) - qtemp(0,j) )
    yp1 = yp1 + 0.5D0*( Vtemp(2,j) - Vtemp(1,j) )/( qtemp(2,j) - qtemp(1,j) )
    ypn = 0.5D0*( Vtemp(np-2,j) - Vtemp(np-3,j) )/( qtemp(np-2,j) - qtemp(np-3,j) )
    ypn = ypn + 0.5D0*( Vtemp(np-1,j) -Vtemp(np-2,j) )/( qtemp(np-1,j) - qtemp(np-2,j) )
    CALL fit_spline(qtemp(1:np-2,j),Vtemp(1:np-2,j),np-2,yp1,ypn,y2(1:npot(j)-2),np-2)

    !outside of trusted potential
    
    !interpolate for each abscissa
    DO i=0,nabs-1

      !catch edge cases
      IF (q(i) .LT. qtemp(1,j)) THEN

        IF (Vtype(j) .EQ. 1) THEN
          beta = LOG(Vtemp(0,j) / Vtemp(1,j))/LOG(qtemp(0,j) / qtemp(1,j))
          alpha = Vtemp(0,j)/ABS(qtemp(0,j))**beta
          val = alpha*ABS(q(i))**beta

        ELSE IF (Vtype(j) .EQ. 2) THEN
          beta = LOG(Vtemp(0,j) / Vtemp(1,j))/LOG(qtemp(0,j) / qtemp(1,j))
          alpha = Vtemp(0,j)/ABS(qtemp(0,j))**beta
          val = alpha*ABS(q(i))**beta

        ELSE IF (Vtype(j) .EQ. 3) THEN
          beta = -1.0D0*LOG(Vtemp(1,j)/Vtemp(0,j))/(qtemp(1,j) - qtemp(0,j))
          alpha = Vtemp(0,j)/EXP(-1.0D0*beta*qtemp(0,j))
          val = alpha*EXP(-1.0D0*beta*q(i))

        ELSE
          WRITE(*,*) "ERROR"
          WRITE(*,*) "V_spline  : only coded to handle potential types 1,2,3"
          error = 1
          RETURN
        END IF

      ELSE IF (q(i) .GT. qtemp(np-2,j)) THEN
        IF (Vtype(j) .EQ. 1) THEN
          beta = LOG(Vtemp(np-1,j) / Vtemp(np-2,j))/LOG(qtemp(np-1,j) / qtemp(np-2,j))
          alpha = Vtemp(np-1,j)/ABS(qtemp(np-1,j))**beta
          val = alpha*ABS(q(i))**beta

        ELSE IF (Vtype(j) .EQ. 2) THEN
          beta = -1.0D0*LOG(Vtemp(np-1,j)/Vtemp(np-2,j))/(qtemp(np-1,j) - qtemp(np-2,j))
          alpha = Vtemp(np-1,j)/EXP(-1.0D0*beta*qtemp(np-1,j))
          val = alpha*EXP(-1.0D0*beta*q(i))
          
        ELSE IF (Vtype(j) .EQ. 3) THEN
          beta = LOG(Vtemp(np-1,j) / Vtemp(np-2,j))/LOG(qtemp(np-1,j) / qtemp(np-2,j))
          alpha = Vtemp(np-1,j)/ABS(qtemp(np-1,j))**beta
          val = alpha*ABS(q(i))**beta
    
        ELSE
          WRITE(*,*) "ERROR"
          WRITE(*,*) "V_spline  : only coded to handle potential types 1,2,3"
          error = 1
          RETURN
        END IF
 
      ELSE
        CALL fit_splint(qtemp(1:np-2,j),Vtemp(1:np-2,j),y2(1:np-2),np,q(i),val,error)
        IF (error .NE. 0) THEN
          WRITE(*,*) "ERROR"
          WRITE(*,*) "V_spline  : error out of fit_splint"
          RETURN
         END IF
       END IF
       CALL valu_check(val,error)
       IF (error .NE. 0) THEN
         WRITE(*,*) "V_spline  : Abscissa", i," of dimension", j," had a bad value"
         RETURN
        ELSE
          Vij(i,j) = val
        END IF

    END DO
    
  END DO 

  WRITE(*,*) "Writing spline output to spline.dat"
  foff = 300
  DO j=0,ndim-1
    fid = foff + j 
    CALL fname_splinedat(j+1,fname,error)
    IF (error .NE. 0) RETURN
    OPEN(file=TRIM(fname),unit=fid,status='replace')
    DO i=0,nabs-1
      WRITE(fid,*) q(i),Vij(i,j)    
    END DO
    CLOSE(unit=fid)
  END DO
  WRITE(*,*)

  DEALLOCATE(y2)
  DEALLOCATE(Vtype)

END SUBROUTINE V_spline
!------------------------------------------------------------

END MODULE V

!------------------------------------------------------------

