!------------------------------------------------------------
! evec
!       - module containing subroutines for dealing with 
!         eigenvectors
!------------------------------------------------------------
MODULE evec
  USE ints_HO

CONTAINS
!------------------------------------------------------------
! evec_print    
!       - prints egenvectors
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! N             : int, length of eigenvectors
! enum          : int, number of eigenvectors
! eval          : 1D real*8, list of eigenvalues
! Cij           : 2D real*8, eigenvectors [N,enum]
! error         : int, exit code

SUBROUTINE evec_print(ndim,nbas,N,enum,eval,Cij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Cij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: eval
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,enum,ndim
  INTEGER, DIMENSION(0:ndim-1) :: Psi,key
  CHARACTER(LEN=1024) :: str_fmt
  INTEGER :: i,j
  error = 0

  CALL ints_HO_key(ndim,nbas,key,error)
  IF (error .NE. 0) RETURN

  str_fmt = "(2x,"
  DO i=0,ndim-1
    str_fmt = TRIM(str_fmt) // "I4,2x,"
  END DO
  str_fmt = TRIM(str_fmt) // "2x,F8.5)"
  str_fmt = TRIM(str_fmt)

  OPEN(file='evec.dat',unit=100,status='replace')
  DO i=0,enum-1

    !write to stdout
    WRITE(*,'(2x,A17,2x,I5)') "Vibrational level", i
    WRITE(*,'(2x,A6,2x,F16.2,1x,A4)') "Energy",eval(i),"cm-1"
    WRITE(*,'(2x,A6,2x,F16.2,1x,A4)') "v0 -> ",eval(i)-eval(0),"cm-1"
    WRITE(*,*) "--------------------------------"
    CALL evec_order(ndim,nbas,N,Cij(0:N-1,i),str_fmt,key,error) 
    WRITE(*,*) "------------------------------------------------"
    WRITE(*,*)

    !write to evec.dat
    WRITE(100,'(2x,A17,2x,I5)') "Vibrational level", i
    WRITE(100,'(2x,A6,2x,F16.2,1x,A4)') "Energy",eval(i),"cm-1"
    WRITE(100,*) "--------------------------------"
    DO j=0,N-1
      CALL ints_HO_qnum(ndim,j,nbas,key,Psi,error)
      IF (error .NE. 0) RETURN
      WRITE(100,str_fmt) Psi,Cij(j,i)
    END DO
    WRITE(100,*) "------------------------------------------------"
    WRITE(100,*)
  END DO
  CLOSE(unit=100)
  
END SUBROUTINE evec_print

!------------------------------------------------------------
! evec_order
!       - orders and prints eigenvector
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, number of basis functions
! N             : int, length of eigenvector
! evec          : 1D real*8, eigenvector
! str3_fmt       : char*1024, string format
! key           : 1D int, keys used for getting Psi
! error         : int, exit code

SUBROUTINE evec_order(ndim,nbas,N,evec,str3_fmt,key,error) 
  IMPLICIT NONE
  CHARACTER(LEN=1024), INTENT(IN) :: str3_fmt
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: evec
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,key
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,N

  REAL(KIND=8), DIMENSION(0:9) :: sigval,tsigval 
  INTEGER, DIMENSION(0:ndim-1) :: Psi
  INTEGER, DIMENSION(0:9) :: sigloc,sgn,tsigloc,tsgn
  CHARACTER(LEN=1024) :: str1_fmt,str1,nstr
  INTEGER :: loc,mnlc
  INTEGER :: i
  
  error = 0

  !str1_fmt = "(2x,"
  !str1 = ""
  !DO i=0,ndim-1
  !  str1_fmt = TRIM(str1_fmt) // "A2,4x,"
  !  WRITE(nstr,'(A1,I1)') "v",i+1
  !  !str1 = TRIM(str1) // TRIM(nstr)
  !  str1 = TRIM(str1) // nstr
  !END DO
  !str1_fmt = TRIM(str1_fmt) // "4x,A4)"
  !str1 = TRIM(str1)//"Coef"
  !!str1 = str1//"Coef"
  !str1_fmt = TRIM(str1_fmt)
  !!str1 = TRIM(str1)
  !str1 = str1

  !search though and find most 10 most significant contributions
  sigval = 0.0D0
  sigloc = -1
  sgn = 1
  mnlc = 0
  DO i=0,N-1
    IF (ABS(evec(i)) .GT. sigval(mnlc)) THEN
      sigval(mnlc) = ABS(evec(i))
      sigloc(mnlc) = i 
      IF (evec(i) .LT. 0.0D0) THEN
        sgn(mnlc) = -1
      ELSE
        sgn(mnlc) = 1
      END IF
      mnlc = MINLOC(sigval,1)-1
    END IF
  END DO

  !dirty sort
  DO i=0,8
    loc = MAXLOC(sigval(i:9),1)-1+i
    IF (loc .NE. i) THEN
      tsigval(i) = sigval(loc)
      tsgn(i) = sgn(loc)
      tsigloc(i) = sigloc(loc)

      sigval(loc) = sigval(i)
      sgn(loc) = sgn(i)
      sigloc(loc) = sigloc(i)

      sigval(i) = tsigval(i)
      sgn(i) = tsgn(i)
      sigloc(i) = tsigloc(i)
      
    END IF
  END DO

  !resign
  DO i=0,9
    sigval(i) = sigval(i)*sgn(i)
  END DO 

  !WRITE(*,str1_fmt) TRIM(str1) 
  !WRITE(*,*) str1
  WRITE(*,*) "v's, Coef"
  WRITE(*,*) "------------------------------------------------"
  DO i=0,MIN(N,9)
    CALL ints_HO_qnum(ndim,sigloc(i),nbas,key,Psi,error)
    IF (error .NE. 0) RETURN
    WRITE(*,str3_fmt) Psi,sigval(i) 
  END DO

END SUBROUTINE evec_order
!------------------------------------------------------------

END MODULE evec
!------------------------------------------------------------
