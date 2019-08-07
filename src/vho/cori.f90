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
    ncori = 0
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
! cori_eval
!       - evaluates coriolis contribution to a particular
!         PsiL and PsiR
!       - we assume that the orthogonality checks for 
!         uninvolved dimensions have already been done
!
!       Integrals are stored like:
!       Q1(i,k) -> <i+1|Qk|i>
!       Q2(2*i,k) -> <i|Qk^2|i>    , Q2(2*i+1,k) -> <i+2|Qk^2|i>
!       P1(i,k) -> <i+1|Pk|i>
!       P2(2*i,k) -> <i|Pk^2|i>    , P2(2*i+1,k)   -> <i+2|Pk^2|i>
!       QP(2*i,k) -> <i|QPk|i>    , QP(2*i+1,k)   -> <i+2|QPk|i>
!       PQ(2*i,k) -> <i|PQk|i>    , PQ(2*i+1,k)   -> <i+2|PQk|i>
!
!       Note that care must be taken to account for factors 
!         of -1 and 1/i
!------------------------------------------------------------
! ndim          : int, number of dimensions
! i,j,k,l       : int, the dimensions invovled in the i,j and k,l
!                      coriolis zeta
! omega         : 1D real*8, harmonic frequencies
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! zetaL         : real*8, LHS zeta 
! zetaR         : real*8, RHS zeta
! Q1            : 2D real*8, <i|q|i'> integrals
! Q2            : 2D real*8, <i|q^2|i'> integrals
! P1            : 2D real*8, <i|p|i'> integrals
! P2            : 2D real*8, <i|p^2|i'> integrals
! QP            : 2D real*8, <i|qp|i'> integrals
! PQ            : 2D real*8, <i|pq|i'> integrals
! val           : real*8, value to iterate
! error         : int, error code
SUBROUTINE cori_eval(ndim,i,j,k,l,omega,PsiL,PsiR,&
                     ZetaL,ZetaR,Q1,Q2,&
                     P1,P2,QP,PQ,val,error)
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Q1,Q2,P1,&
                                              P2,QP,PQ
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: omega
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: val
  REAL(KIND=8), INTENT(IN) :: ZetaL,ZetaR
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,i,j,k,l
  REAL(KIND=8) :: sgni,sgnj,sgnk,sgnl
  INTEGER :: m,n,o,p
  error = 0
   
  !type 1 : i j i j 
  IF (i .NE. j .AND. k .EQ. i .AND. l .EQ. j) THEN
    !i, i
    IF (PsiL(i) .EQ. PsiR(i)) THEN 
      !j, j
      IF (PsiL(j) .EQ. PsiR(j)) THEN
        sgni = 1.0D0*SIGN(1,PsiR(i)-PsiL(i))
        sngj = 1.0D0*SIGN(1,PsiR(j)-PsiL(j))
        m = 2*PsiR(i)
        n = 2*PsiR(j)
        val = val + ZetaL*ZetaR*(omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
              + omega(i)/omega(j)*P2(m,i)*Q2(n,j) + &
              sgni*sngj*QP(m,i)*PQ(n,j) + sgni*sgnj*PQ(m,i)*QP(n,j))
      !j+2 j and j  j+2
      ELSE IF (ABS(PsiL(j) - PsiR(j)) .EQ. 2) THEN
        sgni = 1.0D0*SIGN(1,PsiR(i)-PsiL(i))
        sngj = 1.0D0*SIGN(1,PsiR(j)-PsiL(j))
        m = 2*PsiR(i)
        n = 2*PsiR(j)+1
        val = val + ZetaL*ZetaR*(omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
              + omega(i)/omega(j)*P2(m,i)*Q2(n,j) + &
              sgni*sgnj*QP(m,i)*PQ(n,j) + sgni*sgnj*PQ(m,i)*QP(n,j))
      END IF

    !i+2 i and i i+2
    ELSE IF (ABS(PsiL(i) - PsiR(i)) .EQ. 2) THEN
      !j, j
      IF (PsiL(j) .EQ. PsiR(j)) THEN
        sgni = 1.0D0*SIGN(1,PsiR(i)-PsiL(i))
        sngj = 1.0D0*SIGN(1,PsiR(j)-PsiL(j))
        m = 2*PsiR(i)+1
        n = 2*PsiR(j)
        val = val + ZetaL*ZetaR*(omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
              + omega(i)/omega(j)*P2(m,i)*Q2(n,j) + &
              sgni*sngj*QP(m,i)*PQ(n,j) + sgni*sgnj*PQ(m,i)*QP(n,j))
      !j+2 j and j  j+2
      ELSE IF (ABS(PsiL(j) - PsiR(j)) .EQ. 2) THEN
        sgni = 1.0D0*SIGN(1,PsiR(i)-PsiL(i))
        sngj = 1.0D0*SIGN(1,PsiR(j)-PsiL(j))
        m = 2*PsiR(i)+1
        n = 2*PsiR(j)+1
        val = val + ZetaL*ZetaR*(omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
              + omega(i)/omega(j)*P2(m,i)*Q2(n,j) + &
              sgni*sgnj*QP(m,i)*PQ(n,j) + sgni*sgnj*PQ(m,i)*QP(n,j))
      END IF

    END IF

  !type 2 i j i l 
  ELSE IF (i .NE. j .AND. i .EQ. k .AND. i .NE. l .AND. &
           j .NE. l) THEN
    !i i
    IF (PsiL(i) .EQ. PsiR(i)) THEN
      !j j+1 and l l+1 
      IF (ABS(PsiL(j) - PsiR(j)) .EQ. 1 .AND. &
          ABS(PsiL(l) - PsiL(l)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnl = 1.0D0*SIGN(1,PsiL(l)-PsiR(l))
        m = 2*PsiR(i)
        n = PsiR(j)
        p = PsiR(l)
        val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
         - sgni*sgnj*SQRT(omega(j)/omega(l))*QP(m,i)*P1(n,j)*Q1(p,l) &
         + sgnj*sgnl*SQRT(omega(j)*omega(l))/omega(i)*Q2(m,i)*P1(n,j)*P1(p,l) &
         - omega(i)/SQRT(omega(j)*omega(l))*P2(m,i)*Q1(n,j)*Q1(p,l) &
         - sgni*sgnl*SQRT(omega(l)/omega(j))*PQ(m,i)*Q1(n,j)*P1(p,l) )
      END IF

    !i i+2 and i+2 i
    ELSE IF (ABS(PsiL(i) - PsiR(i)) .EQ. 2) THEN
      !j j+1 and l l+1 
      IF (ABS(PsiL(j) - PsiR(j)) .EQ. 1 .AND. &
          ABS(PsiL(l) - PsiL(l)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnl = 1.0D0*SIGN(1,PsiL(l)-PsiR(l))
        m = 2*PsiR(i)+1
        n = PsiR(j)
        p = PsiR(l)
        val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
         - sgni*sgnj*SQRT(omega(j)/omega(l))*QP(m,i)*P1(n,j)*Q1(p,l) &
         + sgnj*sgnl*SQRT(omega(j)*omega(l))/omega(i)*Q2(m,i)*P1(n,j)*P1(p,l) &
         - omega(i)/SQRT(omega(j)*omega(l))*P2(m,i)*Q1(n,j)*Q1(p,l) &
         - sgni*sgnl*SQRT(omega(l)/omega(j))*PQ(m,i)*Q1(n,j)*P1(p,l) )
      END IF

    END IF

  !type 3 i j k i
  ELSE IF (i .NE. j .AND. i .NE. k .AND. i .EQ. l &
           .AND. j .NE. k) THEN
    !i  i
    IF (PsiL(i) .EQ. PsiR(i)) THEN
      !j j+1 and k k+1
      IF (ABS(PsiL(j) - PsiR(j)) .EQ. 1 .AND. &
          ABS(PsiL(k) - PsiL(k)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnk = 1.0D0*SIGN(1,PsiL(k)-PsiR(k))
        m = 2*PsiR(i)
        n = PsiR(j)
        o = PsiR(k)
        val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
         - sgni*sgnj*SQRT(omega(j)/omega(k))*QP(m,i)*P1(n,j)*Q1(o,k) &
         + sgnj*sgnk*SQRT(omega(j)*omega(k))/omega(i)*Q2(m,i)*P1(n,j)*P1(o,k) &
         - omega(i)/SQRT(omega(j)*omega(k))*P2(m,i)*Q1(n,j)*Q1(o,k) &
         - sgni*sgnk*SQRT(omega(k)/omega(j))*PQ(m,i)*Q1(n,j)*P1(o,k) )
      END IF

    !i i+2 and i+2 i
    ELSE IF (ABS(PsiL(i) - PsiR(i)) .EQ. 2) THEN
      !j j+1 and k k+1
      IF (ABS(PsiL(j) - PsiR(j)) .EQ. 1 .AND. &
          ABS(PsiL(k) - PsiL(k)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnk = 1.0D0*SIGN(1,PsiL(k)-PsiR(k))
        m = 2*PsiR(i)+1
        n = PsiR(j)
        o = PsiR(k)
        val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
         - sgni*sgnj*SQRT(omega(j)/omega(k))*QP(m,i)*P1(n,j)*Q1(o,k) &
         + sgnj*sgnk*SQRT(omega(j)*omega(k))/omega(i)*Q2(m,i)*P1(n,j)*P1(o,k) &
         - omega(i)/SQRT(omega(j)*omega(k))*P2(m,i)*Q1(n,j)*Q1(o,k) &
         - sgni*sgnk*SQRT(omega(k)/omega(j))*PQ(m,i)*Q1(n,j)*P1(o,k) )
      END IF

    END IF

  !type 4 i j j l 
  ELSE IF (i .NE. j .AND. i .NE. l .AND. j .EQ. k .AND. j .NE. l) THEN
    !j j
    IF (PsiL(j) .EQ. PsiR(j)) THEN
      !i i+1  l l+1
      IF (ABS(PsiL(i) - PsiR(i)) .EQ. 1 .AND. &
          ABS(PsiL(l) - PsiL(l)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnl = 1.0D0*SIGN(1,PsiL(l)-PsiR(l))
        m = PsiR(i)
        n = 2*PsiR(j)
        p = PsiR(l)
        val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
         + omega(j)/SQRT(omega(i)*omega(l))*Q1(m,i)*P2(n,j)*Q1(p,l) &
         + sgnj*sgnl*SQRT(omega(l)/omega(i))*Q1(m,i)*PQ(n,j)*Q1(p,l) &
         + sgni*sgnj*SQRT(omega(i)/omega(l))*P1(m,i)*QP(n,j)*Q1(p,l) &
         - sgni*sgnl*SQRT(omega(i)*omega(l))/omega(j)*P1(m,i)*Q2(n,j)*P1(p,l) )
      END IF

    !j j+2
    ELSE IF (ABS(PsiL(j) - PsiR(j)) .EQ. 2) THEN
      !i i+1  l l+1
      IF (ABS(PsiL(i) - PsiR(i)) .EQ. 1 .AND. &
          ABS(PsiL(l) - PsiL(l)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnl = 1.0D0*SIGN(1,PsiL(l)-PsiR(l))
        m = PsiR(i)
        n = 2*PsiR(j)+1
        p = PsiR(l)
        val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
         + omega(j)/SQRT(omega(i)*omega(l))*Q1(m,i)*P2(n,j)*Q1(p,l) &
         + sgnj*sgnl*SQRT(omega(l)/omega(i))*Q1(m,i)*PQ(n,j)*Q1(p,l) &
         + sgni*sgnj*SQRT(omega(i)/omega(l))*P1(m,i)*QP(n,j)*Q1(p,l) &
         - sgni*sgnl*SQRT(omega(i)*omega(l))/omega(j)*P1(m,i)*Q2(n,j)*P1(p,l) )
      END IF

    END IF

  !type 5 i j k j
  ELSE IF (i .NE. j .AND. i .NE. k .AND. j .NE. k .AND. j .EQ. l) THEN
    !j j 
    IF (PsiL(j) .EQ. PsiR(j)) THEN
      !i i+1  k k+1
      IF (ABS(PsiL(i) - PsiR(i)) .EQ. 1 .AND. &
          ABS(PsiL(k) - PsiL(k)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnk = 1.0D0*SIGN(1,PsiL(k)-PsiR(k))
        m = PsiR(i)
        n = 2*PsiR(j)
        o = PsiR(k)
        val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
         + omega(j)/SQRT(omega(i)*omega(k))*Q1(m,i)*P2(n,j)*Q1(o,k) &
         + sgnj*sgnk*SQRT(omega(k)/omega(i))*Q1(m,i)*PQ(n,j)*Q1(o,k) &
         + sgni*sgnj*SQRT(omega(i)/omega(k))*P1(m,i)*QP(n,j)*Q1(o,k) &
         - sgni*sgnk*SQRT(omega(i)*omega(k))/omega(j)*P1(m,i)*Q2(n,j)*P1(o,k) )
      END IF

    !j j+2
    ELSE IF (ABS(PsiL(j) - PsiR(j)) .EQ. 2) THEN
      !i i+1  k k+1
      IF (ABS(PsiL(i) - PsiR(i)) .EQ. 1 .AND. &
          ABS(PsiL(k) - PsiL(k)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnk = 1.0D0*SIGN(1,PsiL(k)-PsiR(k))
        m = PsiR(i)
        n = 2*PsiR(j)+1
        o = PsiR(k)
        val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
         + omega(j)/SQRT(omega(i)*omega(k))*Q1(m,i)*P2(n,j)*Q1(o,k) &
         + sgnj*sgnk*SQRT(omega(k)/omega(i))*Q1(m,i)*PQ(n,j)*Q1(o,k) &
         + sgni*sgnj*SQRT(omega(i)/omega(k))*P1(m,i)*QP(n,j)*Q1(o,k) &
         - sgni*sgnk*SQRT(omega(i)*omega(k))/omega(j)*P1(m,i)*Q2(n,j)*P1(o,k) )
      END IF

    END IF

  !type 6 i j k l
  ELSE IF (i .NE. j .AND. i .NE. k .AND. j .NE. l .AND.&
           j .NE. k .AND. j .NE. l .AND. k .NE. l) THEN
    !all must be within +/- i
    IF (ABS(PsiL(i) - PsiR(i)) .EQ. 1 .AND.&
        ABS(PsiL(j) - PsiR(j)) .EQ. 1 .AND. &
        ABS(PsiL(k) - PsiR(k)) .EQ. 1 .AND. &
        ABS(PsiL(l) - PsiR(l)) .EQ. 1) THEN
      sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
      sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
      sgnk = 1.0D0*SIGN(1,PsiL(k)-PsiR(k))
      sgnl = 1.0D0*SIGN(1,PsiL(l)-PsiR(l))
      m = PsiR(i)
      n = PsiR(j)
      o = PsiR(k)
      p = PsiR(l)
      val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
        - sgnj*sgnl*SQRT(omega(j)*omega(l)/omega(i)/omega(k))*&
                              Q1(m,i)*P1(n,j)*Q1(o,k)*P1(p,l) &
        + sgnj*sgnk*SQRT(omega(j)*omega(k)/omega(i)/omega(l))*&
                              Q1(m,i)*P1(n,j)*P1(o,k)*Q1(p,l) &
        + sgni*sgnl*SQRT(omega(i)*omega(l)/omega(j)/omega(k))*&
                              P1(m,i)*Q1(n,j)*Q1(o,k)*P1(p,l) &
        - sgni*sgnk*SQRT(omega(i)*omega(k)/omega(j)/omega(l))*&
                              P1(m,i)*Q1(n,j)*P1(o,k)*Q1(p,l) )

    END IF

  ELSE
    WRITE(*,*) "cori_eval  : ERROR"
    WRITE(*,*) "James, you missed a case!"
    WRITE(*,*) "PsiL", PsiL
    WRITE(*,*) "i,j",i,j
    WRITE(*,*) "PsiR", PsiR
    WRITE(*,*) "k,l",k,l
    error = 1
    RETURN
  END IF 

END SUBROUTINE cori_eval

!------------------------------------------------------------
END MODULE cori
!------------------------------------------------------------
