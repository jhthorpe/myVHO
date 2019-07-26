!------------------------------------------------------------
! V
!       - module containing subroutines dealing with the
!         the potential energy of the dimensions
!------------------------------------------------------------
MODULE V
  USE fname
  USE input
  USE fit
  USE val

CONTAINS
!------------------------------------------------------------
! V_get
!       - gets the potential energy of dimensions
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of dimensions
! nabs          : int, number of abscissa
! q             : 1D real*8, abscissa        
! Vij           : 2D real*8, potential energy [abscissa,dimension]
! error         : int, exit code          

SUBROUTINE V_get(job,ndim,nabs,q,Vij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Vij
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: q
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim,nabs

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Vtemp,qtemp
  INTEGER, DIMENSION(0:ndim-1) :: npot
  INTEGER :: i,j

  error = 0
  !read in potential energies
  CALL V_read(job,ndim,nabs,npot,qtemp,Vtemp,error)
  IF (error .NE. 0) RETURN 
  IF (job .EQ. 0) THEN
    Vij = Vtemp
    q = qtemp(:,0)
  END IF
  
  IF (job .EQ. 1) THEN
    CALL V_spline(ndim,nabs,npot,q,qtemp,Vtemp,Vij,error)
    IF (error .NE. 0) THEN
      RETURN
    END IF 
  END IF
 
  !adjust potentials to be zero
  DO j=0,ndim-1
    Vij(0:nabs-1,j) = Vij(0:nabs-1,j) - MINVAL(Vij(0:nabs-1,j))
  END DO

  Vij = Vij*219474.63 !convert hartrees to cm-1

  IF(ALLOCATED(Vtemp)) DEALLOCATE(Vtemp)
  IF(ALLOCATED(qtemp)) DEALLOCATE(qtemp)

END SUBROUTINE V_get

!------------------------------------------------------------
! V_read
!       - read in potential energies from Vx.in files
!------------------------------------------------------------
! job           : int, jobtype
! ndim          : int, number of dimensions
! nabs          : int, number of abscissa
! npot          : 1D int, number of points per dimension
! Vtemp         : 2D real*8, potentials at points [potential, dimension]
! qtemp         : 2D real*8, normco at points [normco, dimension]
! error         : int, error code

SUBROUTINE V_read(job,ndim,nabs,npot,qtemp,Vtemp,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Vtemp,qtemp
  INTEGER, DIMENSION(0:), INTENT(INOUT) :: npot
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,job,nabs

  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val1,val2
  INTEGER :: foff,fid,dummy 
  INTEGER :: i,j

  error = 0
  foff = 200
  IF (job .EQ. 1 .OR. job .EQ. 0) THEN  !if we are not already in abscissa form
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
        READ(fid,*) dummy, val1, val2
        dummy = dummy - 1
        qtemp(dummy,j) = val1
        Vtemp(dummy,j) = val2
        !READ(fid,*) dummy, qtemp(i,j),Vtemp(i,j)
      END DO
      CLOSE(unit=fid)
    END DO

  ELSE !we are in some other form 
    WRITE(*,*) "V_read  : ERROR"
    WRITE(*,*) "Bad jobtype" 
    error = 1
    RETURN
  END IF

  IF (job .EQ. 0) THEN
    IF (ANY(npot .NE. nabs)) THEN
      WRITE(*,*) "V_read  : ERROR"
      WRITE(*,*) "Incorrect number of abscissa provided"
    END IF
  END IF

  !WRITE(*,*) "TESTING TESTING TESTING"
  !DO i=0,ndim-1
  !  WRITE(*,*) "------------------------"
  !  DO j=0,npot(i)-1
  !    WRITE(*,*) j,qtemp(j,i),Vtemp(j,i)
  !  END DO
  !END DO

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
       CALL val_check(val,error)
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

