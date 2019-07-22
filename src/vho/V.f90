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
! Vij           : 2D real*8, potential energy [ndim,nabs]
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
  CALL V_read(job,ndim,npot,qtemp,Vtemp,error)
  IF (error .NE. 0) RETURN 

  !Generate the V at abscissa
  IF (job .EQ. 1) THEN
    CALL V_spline(ndim,nabs,npot,q,qtemp,Vtemp,Vij,error)
    IF (error .NE. 0) THEN
      RETURN
    END IF 
  END IF

  IF(ALLOCATED(Vtemp)) DEALLOCATE(Vtemp)
  IF(ALLOCATED(qtemp)) DEALLOCATE(qtemp)

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
! npot          : 1D int, number of points per potential
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
  REAL(KIND=8) :: yp1,ypn,beta,alpha,val
  INTEGER :: np
  INTEGER :: i,j

  error = 0
  WRITE(*,*) 
  WRITE(*,*) "Performing cubic spline interpolation of surface"

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
  !DO i=0,ndim-1
  !  Vtype = 1
  !END DO
  Vtype(0) = 1
  Vtype(1) = 2

  !WRITE(*,*) "TESTING TESTING TESTING" 
  !WRITE(*,*) "qtemp0 is", qtemp(:,0)
  !WRITE(*,*) "qtemp1 is", qtemp(:,1)
 
  ALLOCATE(y2(0:MAXVAL(npot)-1))
  !generate y2 for each dimension 
  DO i=0,ndim-1
    np = npot(i)
    yp1 = 0.5D0*( Vtemp(1,i) - Vtemp(0,i) )/( qtemp(1,i) - qtemp(0,i) )
    yp1 = yp1 + 0.5D0*( Vtemp(2,i) - Vtemp(1,i) )/( qtemp(2,i) - qtemp(1,i) )
    ypn = 0.5D0*( Vtemp(np-2,i) - Vtemp(np-3,i) )/( qtemp(np-2,i) - qtemp(np-3,i) )
    ypn = ypn + 0.5D0*( Vtemp(np-1,i) -Vtemp(np-2,i) )/( qtemp(np-1,i) - qtemp(np-2,i) )
    CALL fit_spline(qtemp(1:np-2,i),Vtemp(1:np-2,i),np-2,yp1,ypn,y2(1:npot(i)-2),np-2)
    
    !interpolate for each abscissa
    DO j=0,nabs-1

      !catch edge cases
      IF (q(j) .LT. qtemp(1,i)) THEN

        IF (Vtype(i) .EQ. 1) THEN
          beta = LOG(Vtemp(0,i) / Vtemp(1,i))/LOG(qtemp(0,i) / qtemp(1,i))
          alpha = Vtemp(0,i)/ABS(qtemp(0,i))**beta
          val = alpha*ABS(q(j))**beta

        ELSE IF (Vtype(i) .EQ. 2) THEN
          beta = LOG(Vtemp(0,i) / Vtemp(1,i))/LOG(qtemp(0,i) / qtemp(1,i))
          alpha = Vtemp(0,i)/ABS(qtemp(0,i))**beta
          val = alpha*ABS(q(j))**beta
          WRITE(*,*) "alpha is", alpha
          WRITE(*,*) "beta is", beta

        ELSE IF (Vtype(i) .EQ. 3) THEN
          beta = -1.0D0*LOG(Vtemp(1,i)/Vtemp(0,i))/(qtemp(1,i) - qtemp(0,i))
          alpha = Vtemp(0,i)/EXP(-1.0D0*beta*qtemp(0,i))
          val = alpha*EXP(-1.0D0*beta*q(j))

        ELSE
          WRITE(*,*) "ERROR"
          WRITE(*,*) "V_spline  : only coded to handle potential types 1,2,3"
          error = 1
          RETURN
        END IF

      ELSE IF (q(j) .GT. qtemp(np-2,i)) THEN
        IF (Vtype(i) .EQ. 1) THEN
          beta = LOG(Vtemp(np-1,i) / Vtemp(np-2,i))/LOG(qtemp(np-1,i) / qtemp(np-2,i))
          alpha = Vtemp(np-1,i)/ABS(qtemp(np-1,i))**beta
          val = alpha*ABS(q(j))**beta

        ELSE IF (Vtype(i) .EQ. 2) THEN
          beta = -1.0D0*LOG(Vtemp(np-1,i)/Vtemp(np-2,i))/(qtemp(np-1,i) - qtemp(np-2,i))
          alpha = Vtemp(np-1,i)/EXP(-1.0D0*beta*qtemp(np-1,i))
          val = alpha*EXP(-1.0D0*beta*q(j))
          
        ELSE IF (Vtype(i) .EQ. 3) THEN
          beta = LOG(Vtemp(np-1,i) / Vtemp(np-2,i))/LOG(qtemp(np-1,i) / qtemp(np-2,i))
          alpha = Vtemp(np-1,i)/ABS(qtemp(np-1,i))**beta
          val = alpha*ABS(q(j))**beta
    
        ELSE
          WRITE(*,*) "ERROR"
          WRITE(*,*) "V_spline  : only coded to handle potential types 1,2,3"
          error = 1
          RETURN
        END IF
 
      ELSE
        CALL fit_splint(qtemp(1:np-2,i),Vtemp(1:np-2,i),y2(1:np-2),np,q(j),val,error)
        IF (error .NE. 0) THEN
          WRITE(*,*) "ERROR"
          WRITE(*,*) "V_spline  : error out of fit_splint"
          !RETURN
         END IF
       END IF
       CALL val_check(val,error)
       IF (error .NE. 0) THEN
         WRITE(*,*) "V_spline  : Abscissa", j," of dimension", i," had a bad value"
         !RETURN
        ELSE
          Vij(i,j) = val
        END IF

    END DO
    
    !adjust potentials to be zero
    Vij(i,0:nabs-1) = Vij(i,0:nabs-1) - MINVAL(Vij(i,0:nabs-1))

  END DO 

  WRITE(*,*) "Writing spline output to spline.dat"
  OPEN(file="spline.dat",unit=150,status='replace')
  DO j=0,nabs-1
    WRITE(150,*) q(j),Vij(0:ndim-1,j)
  END DO
  CLOSE(unit=150)

  DEALLOCATE(y2)
  DEALLOCATE(Vtype)

  STOP

END SUBROUTINE V_spline
!------------------------------------------------------------

END MODULE V

!------------------------------------------------------------

