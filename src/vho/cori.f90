!------------------------------------------------------------
! cori
!       - module containing subroutines concerning 
!         coriolis coupling
!------------------------------------------------------------
MODULE cori
  USE input
  USE sort
  USE ints_HO

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
      IF (idx(0) .NE. idx0(0)-voff) val = -1.0D0*val !if there was a swap
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
        WRITE(*,'(2x,I1,2x,2(I4,2x),F24.15)') a+1,qcori(2*j:2*j+1,a)+voff+1,cori(j,a)
      END DO
    END DO
    WRITE(*,*) 
  END IF

  IF (ALLOCATED(Itemp)) DEALLOCATE(Itemp)
  IF (ALLOCATED(Rtemp)) DEALLOCATE(Rtemp)
   
END SUBROUTINE cori_get

!------------------------------------------------------------
! cori_get
!       - reads in, sorts, and stores the coriolis coupling
!         constants as a 3 index matrix
!       - the qantum numbers of the constants are stored
!         [i,j,alpha], all indexed from zero
!------------------------------------------------------------
! ndim          : int, number of dimensions
! cori          : 2D real*8, coriolis zetas [zeta,rotation]
! error         : int, exit code

SUBROUTINE cori_get_debug(ndim,cori,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: cori
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  CHARACTER(LEN=1024) :: fname
  REAL(KIND=8) :: val
  LOGICAL :: ex
  INTEGER :: fline,fid,voff
  INTEGER :: i,j,a,k

  error = 0
  fname = 'coriolis'
  fid = 300

  ALLOCATE(cori(0:ndim-1,0:ndim-1,0:2))
  cori = 0.0D0
  
  !If only one dimension
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

  !Initialize
  CALL input_nline(fline,fname)

  !read in data
  OPEN(file=TRIM(fname),unit=fid,status='old')
  DO k=0,fline-1 
    READ(fid,*) a,i,j,val

    i = i - voff - 1
    j = j - voff - 1
    a = a - 1

    IF (a .LT. 0 .OR. a .GT. 2 .OR. i .LT. 0 .OR. i .GT. ndim-1 .OR. &
      j .LT. 0 .OR. j .GT. ndim-1) THEN
      WRITE(*,*) "cori_get_debug  : ERROR"
      WRITE(*,*) "In coriolis, line ",k,", has a bad value"
      error = 1
    END IF

    cori(i,j,a) = val
  END DO
  CLOSE(unit=fid)

  WRITE(*,*) "Coriolis Constants -- DEBUG"
  DO a=0,2
    DO i=0,ndim-1
      DO j=i+1,ndim-1
        IF (ABS(cori(i,j,a)) .GT. 1.0D-15) THEN
         WRITE(*,'(2x,I1,2x,2(I4,2x),F24.15)') a+1,i+voff+1,j+voff+1,cori(i,j,a)
         WRITE(*,'(2x,I1,2x,2(I4,2x),F24.15)') a+1,j+voff+1,i+voff+1,cori(j,i,a)
        END IF
      END DO
    END DO
  END DO 
  WRITE(*,*)


END SUBROUTINE cori_get_debug

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
    !!WRITE(*,*) "case 1"
    !!WRITE(*,*) "type: ", i,j,k,l
    !!WRITE(*,*) "PsiL(i),PsiR(i)",PsiL(i),PsiR(i)
    !!WRITE(*,*) "PsiL(j),PsiR(j)",PsiL(j),PsiR(j)
    !i, i
    IF (PsiL(i) .EQ. PsiR(i)) THEN 
      !j, j
      IF (PsiL(j) .EQ. PsiR(j)) THEN
        sgni = 1.0D0*SIGN(1,PsiR(i)-PsiL(i))
        sgnj = 1.0D0*SIGN(1,PsiR(j)-PsiL(j))
      !  m = 2*PsiR(i)
      !  n = 2*PsiR(j)
        m = 2*MIN(PsiL(i),PsiR(i))
        n = 2*MIN(PsiL(j),PsiR(j))
        val = val + ZetaL*ZetaR*(0.0D0 &
              + omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
              + omega(i)/omega(j)*P2(m,i)*Q2(n,j) &
              + sgni*sgnj*QP(m,i)*PQ(n,j) &
              + sgni*sgnj*PQ(m,i)*QP(n,j))
!      !WRITE(*,*) "Contribution to val:",&
! ZetaL*ZetaR*(omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
!+ omega(i)/omega(j)*P2(m,i)*Q2(n,j) + &
!sgni*sgnj*QP(m,i)*PQ(n,j) + sgni*sgnj*PQ(m,i)*QP(n,j))
! ZetaL*ZetaR*(omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
!+ omega(i)/omega(j)*P2(m,i)*Q2(n,j) + &
!QP(m,i)*PQ(n,j) + PQ(m,i)*QP(n,j))
      !j+2 j and j  j+2
      ELSE IF (ABS(PsiL(j) - PsiR(j)) .EQ. 2) THEN
        sgni = 1.0D0*SIGN(1,PsiR(i)-PsiL(i))
        sgnj = 1.0D0*SIGN(1,PsiR(j)-PsiL(j))
       ! m = 2*PsiR(i)
       ! n = 2*PsiR(j)+1
        m = 2*MIN(PsiL(i),PsiR(i))+1
        n = 2*MIN(PsiL(j),PsiR(j))
        val = val + ZetaL*ZetaR*(omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
              + omega(i)/omega(j)*P2(m,i)*Q2(n,j) + &
              sgni*sgnj*QP(m,i)*PQ(n,j) + sgni*sgnj*PQ(m,i)*QP(n,j))
      END IF

    !i+2 i and i i+2
    ELSE IF (ABS(PsiL(i) - PsiR(i)) .EQ. 2) THEN
      !j, j
      IF (PsiL(j) .EQ. PsiR(j)) THEN
        sgni = 1.0D0*SIGN(1,PsiR(i)-PsiL(i))
        sgnj = 1.0D0*SIGN(1,PsiR(j)-PsiL(j))
        !m = 2*PsiR(i)+1
        !n = 2*PsiR(j)
        m = 2*MIN(PsiL(i),PsiR(i))+1
        n = 2*MIN(PsiL(j),PsiR(j))
        val = val + ZetaL*ZetaR*(omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
              + omega(i)/omega(j)*P2(m,i)*Q2(n,j) + &
              sgni*sgnj*QP(m,i)*PQ(n,j) + sgni*sgnj*PQ(m,i)*QP(n,j))
      !j+2 j and j  j+2
      ELSE IF (ABS(PsiL(j) - PsiR(j)) .EQ. 2) THEN
        sgni = 1.0D0*SIGN(1,PsiR(i)-PsiL(i))
        sgnj = 1.0D0*SIGN(1,PsiR(j)-PsiL(j))
       ! m = 2*PsiR(i)+1
       ! n = 2*PsiR(j)+1
        m = 2*MIN(PsiL(i),PsiR(i))+1
        n = 2*MIN(PsiL(j),PsiR(j))+1
        val = val + ZetaL*ZetaR*(omega(j)/omega(i)*Q2(m,i)*P2(n,j) &
              + omega(i)/omega(j)*P2(m,i)*Q2(n,j) + &
              sgni*sgnj*QP(m,i)*PQ(n,j) + sgni*sgnj*PQ(m,i)*QP(n,j))
      END IF

    END IF

  !type 2 i j i l 
  ELSE IF (i .NE. j .AND. i .EQ. k .AND. i .NE. l .AND. &
           j .NE. l) THEN
    !!WRITE(*,*) "case 2"
    !!WRITE(*,*) "type: ", i,j,k,l
    !!WRITE(*,*) "PsiL(i),PsiR(i)",PsiL(i),PsiR(i)
    !!WRITE(*,*) "PsiL(j),PsiR(j)",PsiL(j),PsiR(j)
    !!WRITE(*,*) "PsiL(l),PsiR(l)",PsiL(l),PsiR(l)
    !i i
    IF (PsiL(i) .EQ. PsiR(i)) THEN
      !j j+1 and l l+1 
      IF (ABS(PsiL(j) - PsiR(j)) .EQ. 1 .AND. &
          ABS(PsiL(l) - PsiL(l)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnl = 1.0D0*SIGN(1,PsiL(l)-PsiR(l))
       ! m = 2*PsiR(i)
       ! n = PsiR(j)
       ! p = PsiR(l)
        m = 2*MIN(PsiL(i),PsiR(i))
        n = MIN(PsiL(j),PsiR(j))
        p = MIN(PsiL(l),PsiR(l))
       ! val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
        val = val + ZetaL*ZetaR*(0.0D0 &
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
      !  m = 2*PsiR(i)+1
      !  n = PsiR(j)
      !  p = PsiR(l)
        m = 2*MIN(PsiL(i),PsiR(i))+1
        n = MIN(PsiL(j),PsiR(j))
        p = MIN(PsiL(l),PsiR(l))
       ! val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
        val = val + ZetaL*ZetaR*(0.0D0 &
         - sgni*sgnj*SQRT(omega(j)/omega(l))*QP(m,i)*P1(n,j)*Q1(p,l) &
         + sgnj*sgnl*SQRT(omega(j)*omega(l))/omega(i)*Q2(m,i)*P1(n,j)*P1(p,l) &
         - omega(i)/SQRT(omega(j)*omega(l))*P2(m,i)*Q1(n,j)*Q1(p,l) &
         - sgni*sgnl*SQRT(omega(l)/omega(j))*PQ(m,i)*Q1(n,j)*P1(p,l) )
      END IF

    END IF

  !type 3 i j k i
  ELSE IF (i .NE. j .AND. i .NE. k .AND. i .EQ. l &
           .AND. j .NE. k) THEN
    !WRITE(*,*) "case 3"
    !WRITE(*,*) "type: ", i,j,k,l
    !WRITE(*,*) "PsiL(i),PsiR(i)",PsiL(i),PsiR(i)
    !WRITE(*,*) "PsiL(j),PsiR(j)",PsiL(j),PsiR(j)
    !WRITE(*,*) "PsiL(k),PsiR(k)",PsiL(k),PsiR(k)
    !i  i
    IF (PsiL(i) .EQ. PsiR(i)) THEN
      !j j+1 and k k+1
      IF (ABS(PsiL(j) - PsiR(j)) .EQ. 1 .AND. &
          ABS(PsiL(k) - PsiL(k)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnk = 1.0D0*SIGN(1,PsiL(k)-PsiR(k))
       ! m = 2*PsiR(i)
       ! n = PsiR(j)
       ! o = PsiR(k)
        m = 2*MIN(PsiL(i),PsiR(i))
        n = MIN(PsiL(j),PsiR(j))
        o = MIN(PsiL(k),PsiR(k))
       ! val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
        val = val + ZetaL*ZetaR*(0.0D0 &
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
       ! m = 2*PsiR(i)+1
       ! n = PsiR(j)
       ! o = PsiR(k)
        m = 2*MIN(PsiL(i),PsiR(i))+1
        n = MIN(PsiL(j),PsiR(j))
        o = MIN(PsiL(k),PsiR(k))
       ! val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
        val = val + ZetaL*ZetaR*(0.0D0 &
         - sgni*sgnj*SQRT(omega(j)/omega(k))*QP(m,i)*P1(n,j)*Q1(o,k) &
         + sgnj*sgnk*SQRT(omega(j)*omega(k))/omega(i)*Q2(m,i)*P1(n,j)*P1(o,k) &
         - omega(i)/SQRT(omega(j)*omega(k))*P2(m,i)*Q1(n,j)*Q1(o,k) &
         - sgni*sgnk*SQRT(omega(k)/omega(j))*PQ(m,i)*Q1(n,j)*P1(o,k) )
      END IF

    END IF

  !type 4 i j j l 
  ELSE IF (i .NE. j .AND. i .NE. l .AND. j .EQ. k .AND. j .NE. l) THEN
    !WRITE(*,*) "case 4"
    !WRITE(*,*) "type: ", i,j,k,l
    !WRITE(*,*) "PsiL(i),PsiR(i)",PsiL(i),PsiR(i)
    !WRITE(*,*) "PsiL(j),PsiR(j)",PsiL(j),PsiR(j)
    !WRITE(*,*) "PsiL(l),PsiR(l)",PsiL(l),PsiR(l)
    !j j
    IF (PsiL(j) .EQ. PsiR(j)) THEN
      !i i+1  l l+1
      IF (ABS(PsiL(i) - PsiR(i)) .EQ. 1 .AND. &
          ABS(PsiL(l) - PsiL(l)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnl = 1.0D0*SIGN(1,PsiL(l)-PsiR(l))
       ! m = PsiR(i)
       ! n = 2*PsiR(j)
       ! p = PsiR(l)
        m = MIN(PsiL(i),PsiR(i))
        n = 2*MIN(PsiL(j),PsiR(j))
        p = MIN(PsiL(l),PsiR(l))
       ! val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
        val = val + ZetaL*ZetaR*(0.0D0 &
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
       ! m = PsiR(i)
       ! n = 2*PsiR(j)+1
       ! p = PsiR(l)
        m = MIN(PsiL(i),PsiR(i))
        n = 2*MIN(PsiL(j),PsiR(j))+1
        p = MIN(PsiL(l),PsiR(l))
       ! val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
        val = val + ZetaL*ZetaR*(0.0D0 &
         + omega(j)/SQRT(omega(i)*omega(l))*Q1(m,i)*P2(n,j)*Q1(p,l) &
         + sgnj*sgnl*SQRT(omega(l)/omega(i))*Q1(m,i)*PQ(n,j)*Q1(p,l) &
         + sgni*sgnj*SQRT(omega(i)/omega(l))*P1(m,i)*QP(n,j)*Q1(p,l) &
         - sgni*sgnl*SQRT(omega(i)*omega(l))/omega(j)*P1(m,i)*Q2(n,j)*P1(p,l) )
      END IF

    END IF

  !type 5 i j k j
  ELSE IF (i .NE. j .AND. i .NE. k .AND. j .NE. k .AND. j .EQ. l) THEN
    !WRITE(*,*) "case 5"
    !WRITE(*,*) "type: ", i,j,k,l
    !WRITE(*,*) "PsiL(i),PsiR(i)",PsiL(i),PsiR(i)
    !WRITE(*,*) "PsiL(j),PsiR(j)",PsiL(j),PsiR(j)
    !WRITE(*,*) "PsiL(k),PsiR(k)",PsiL(k),PsiR(k)
    !j j 
    IF (PsiL(j) .EQ. PsiR(j)) THEN
      !i i+1  k k+1
      IF (ABS(PsiL(i) - PsiR(i)) .EQ. 1 .AND. &
          ABS(PsiL(k) - PsiL(k)) .EQ. 1) THEN 
        sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
        sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
        sgnk = 1.0D0*SIGN(1,PsiL(k)-PsiR(k))
       ! m = PsiR(i)
       ! n = 2*PsiR(j)
       ! o = PsiR(k)
        m = MIN(PsiL(i),PsiR(i))
        n = 2*MIN(PsiL(j),PsiR(j))
        o = MIN(PsiL(k),PsiR(k))
       ! val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
        val = val + ZetaL*ZetaR*(0.0D0 &
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
      !  m = PsiR(i)
      !  n = 2*PsiR(j)+1
      !  o = PsiR(k)
        m = MIN(PsiL(i),PsiR(i))
        n = 2*MIN(PsiL(j),PsiR(j))+1
        o = MIN(PsiL(k),PsiR(k))
       ! val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
        val = val + ZetaL*ZetaR*(0.0D0 &
         + omega(j)/SQRT(omega(i)*omega(k))*Q1(m,i)*P2(n,j)*Q1(o,k) &
         + sgnj*sgnk*SQRT(omega(k)/omega(i))*Q1(m,i)*PQ(n,j)*Q1(o,k) &
         + sgni*sgnj*SQRT(omega(i)/omega(k))*P1(m,i)*QP(n,j)*Q1(o,k) &
         - sgni*sgnk*SQRT(omega(i)*omega(k))/omega(j)*P1(m,i)*Q2(n,j)*P1(o,k) )
      END IF

    END IF

  !type 6 i j k l
  ELSE IF (i .NE. j .AND. i .NE. k .AND. j .NE. l .AND.&
           j .NE. k .AND. j .NE. l .AND. k .NE. l) THEN
    !WRITE(*,*) "case 6"
    !WRITE(*,*) "type: ", i,j,k,l
    !WRITE(*,*) "PsiL(i),PsiR(i)",PsiL(i),PsiR(i)
    !WRITE(*,*) "PsiL(j),PsiR(j)",PsiL(j),PsiR(j)
    !WRITE(*,*) "PsiL(k),PsiR(k)",PsiL(k),PsiR(k)
    !WRITE(*,*) "PsiL(l),PsiR(l)",PsiL(l),PsiR(l)
    !all must be within +/- i
    IF (ABS(PsiL(i) - PsiR(i)) .EQ. 1 .AND.&
        ABS(PsiL(j) - PsiR(j)) .EQ. 1 .AND. &
        ABS(PsiL(k) - PsiR(k)) .EQ. 1 .AND. &
        ABS(PsiL(l) - PsiR(l)) .EQ. 1) THEN
      sgni = 1.0D0*SIGN(1,PsiL(i)-PsiR(i))
      sgnj = 1.0D0*SIGN(1,PsiL(j)-PsiR(j))
      sgnk = 1.0D0*SIGN(1,PsiL(k)-PsiR(k))
      sgnl = 1.0D0*SIGN(1,PsiL(l)-PsiR(l))
      m = MIN(PsiL(i),PsiR(i))
      n = MIN(PsiL(j),PsiR(j))
      o = MIN(PsiL(k),PsiR(k))
      p = MIN(PsiL(l),PsiR(l))
      !m = PsiR(i)
      !n = PsiR(j)
      !o = PsiR(k)
      !p = PsiR(l)
      !val = val + 2.0D0*ZetaL*ZetaR*(0.0D0 &
      val = val + ZetaL*ZetaR*(0.0D0 &
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
! cori_HO_O2_old
!	- evalute 2nd order coriolis contributions in the
!	  harmonic oscillator basis
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nrota         : int, number of rotational terms
! rota          : 1D real*8, rotational terms
! omega         : 1D real*8, harmonic frequencies stored in
!                            order omega[0] -> omega for dim 0 
! ncori         : 1D int, number of coriolis interactions per
!                         rotational mode
! qcori         : 2D int, coriolis quantum numbers for 
!                         each rotational mode
! cori          : 2D real*8, coriolis zetas
! Q1int         : 2D real*8, <i|q|i'> integrals
! Q2int         : 2D real*8, <i|q^2|i'> integrals
! P1int         : 2D real*8, <i|p|i'> integrals
! P2int         : 2D real*8, <i|p^2|i'> integrals
! QPint         : 2D real*8, <i|qp|i'> integrals
! PQint         : 2D real*8, <i|pq|i'> integrals
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! cval          : real*8, value to output
! error         : int, error code

SUBROUTINE cori_HO_O2_old(ndim,nrota,rota,omega,ncori,qcori,cori,Q1int,&
                            Q2int,P1int,P2int,QPint,PQint,PsiL,PsiR,cval,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: cori,Q1int,Q2int,P1int,&
                                                P2int,QPint,PQint
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota,omega
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: qcori
  INTEGER, DIMENSION(0:), INTENT(IN) :: ncori,PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: cval
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nrota
  INTEGER, DIMENSION(0:3) :: v
  REAL(KIND=8) :: val
  INTEGER :: b,c,a,i,j,k,l
  !WRITE(*,*) "ncori is",ncori
  error = 0
  cval = 0.0D0

  !WRITE(*,*) 
  !WRITE(*,*) 
  !WRITE(*,*) 
  !WRITE(*,*) 
  !WRITE(*,*) "TESTING TESTING TESTING"
  !WRITE(*,*) "PsiL is", PsiL
  !WRITE(*,*) "PsiR is", PsiR
 
  DO a=0,nrota-1
  !  WRITE(*,*) "Rotational mode:", a
  !  WRITE(*,*) "Be^a",rota(a)
    val = 0.0D0
    DO b=0,ncori(a)-1
  !    WRITE(*,*) "RHS Coriolis term"
  !    WRITE(*,*) qcori(2*b,a),qcori(2*b+1,a),cori(b,a)
      k = qcori(2*b,a)
      l = qcori(2*b+1,a)
      !DO c=b,ncori(a)-1
      DO c=0,ncori(a)-1
  !      WRITE(*,*)
  !      WRITE(*,*) "LHS Coriolis term"
  !      WRITE(*,*) qcori(2*c,a),qcori(2*c+1,a),cori(c,a) 
        i = qcori(2*c,a)
        j = qcori(2*c+1,a)
        v = [i,j,k,l]
        CALL sort_int_ijkl(v)
        !check orthognality
        IF (ANY(PsiL(0:v(0)-1) .NE. PsiR(0:v(0)-1)) .OR. &
            ANY(PsiL(v(0)+1:v(1)-1) .NE. PsiR(v(0)+1:v(1)-1)) .OR. &
            ANY(PsiL(v(1)+1:v(2)-1) .NE. PsiR(v(1)+1:v(2)-1)) .OR. &
            ANY(PsiL(v(2)+1:v(3)-1) .NE. PsiR(v(2)+1:v(3)-1)) .OR. &
            ANY(PsiL(v(3)+1:ndim-1) .NE. PsiR(v(3)+1:ndim-1))) THEN
  !       WRITE(*,*) "Orthogonal, skipping"
         CYCLE
        ELSE
 !         WRITE(*,*) "val before is",val 
          CALL cori_eval(ndim,i,j,k,l,omega,PsiL,PsiR,&
                         cori(c,a),cori(b,a),Q1int,Q2int,&
                         P1int,P2int,QPint,PQint,val,error)
 !       WRITE(*,*) "val after is",val 
          IF (error .NE. 0) RETURN
        END IF
      END DO
 !     WRITE(*,*) "--------------------------------------------------"
    END DO
 !   WRITE(*,*) "val added is",rota(a)*val 
    cval = cval + rota(a)*val
 !   WRITE(*,*) "===================================================="
  END DO

  !WRITE(*,*) "coriolis is  =",cval
END SUBROUTINE cori_HO_O2_old

!------------------------------------------------------------
! cori_HO_O2_ints
!	- evaluate 2nd order coriolis integrals in the 
!	  harmonic oscillator basis
!	- this assumes the orthogonality checks have 
!    	  been performed
!------------------------------------------------------------
! ndim		: int, number of dimensions
! PsiL		: 1D int, LHS quantum numbers
! PsiR		: 1D int, RHF quantum numbers
! ids		: 1D int, [i,j,k,l] ids
! val		: real*8, value of integrals

SUBROUTINE cori_HO_O2_ints(ndim,PsiL,PsiR,ids,val)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR,ids
  REAL(KIND=8), INTENT(INOUT) :: val
  INTEGER, INTENT(IN) :: ndim
  LOGICAL, DIMENSION(0:4-1,0:4-1) :: isq,isp
  INTEGER, DIMENSION(0:4-1,0:4-1) :: positions
  INTEGER, DIMENSION(0:4-1) :: values,counts
  REAL(KIND=8) :: temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 4
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0
 
  !WRITE(*,*) "int ids",ids+7

  !Gather information on positions
  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2) THEN
          isq(0,j) = .TRUE. 
          isp(0,j) = .FALSE.
        ELSE
          isq(0,j) = .FALSE. 
          isp(0,j) = .TRUE.
        END IF
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2) THEN
          isq(counts(j),j) = .TRUE. 
          isp(counts(j),j) = .FALSE.
        ELSE
          isq(counts(j),j) = .FALSE. 
          isp(counts(j),j) = .TRUE.
        END IF
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO

  !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> P
  ! 2 -> Q
  ! 3 -> P
  DO i=0,nids-1
    IF (counts(i) .EQ. 0) THEN 
      EXIT
    
    ELSE IF (counts(i) .EQ. 1) THEN
      !q
      IF (isq(0,i)) THEN
        temp = ints_HO_q(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "q",values(i)
      !p
      ELSE IF (isp(0,i)) THEN
        temp = ints_HO_p(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "p",values(i)
      ELSE
        WRITE(*,*) "cori_HO_O2_ints  : ERROR"
        WRITE(*,*) "There is a 1 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 2) THEN
      !qq
      IF (isq(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_qq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qq",values(i)
      !qp
      ELSE IF (isq(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_qp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qp",values(i)
      !pq
      ELSE IF (isp(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_pq(PsiL(values(i)),PsiR(values(i))) 
!        WRITE(*,*) "pq",values(i)
      !pp
      ELSE IF (isp(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_pp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pp",values(i)
      ELSE
        WRITE(*,*) "cori_HO_O2_ints  : ERROR"
        WRITE(*,*) "There is a 2 index case that has been missed"
        STOP
      END IF

    ELSE
      WRITE(*,*) "cori_HO_O2_ints  : ERROR"
      WRITE(*,*) "There is a bad case here"
      WRITE(*,*) "PsiL",PsiL
      WRITE(*,*) "PsiR",PsiR
      WRITE(*,*) "ids",ids
      STOP
    END IF

    val = val*temp
  END DO

  !accont for 1/i factor of the two p's
  val = -1.0D0*val

END SUBROUTINE cori_HO_O2_ints

!------------------------------------------------------------
! cori_HO_O3_pup_ints
!       - evalutes integrals of the kind 
!         <qi pj qm qk pl>
!	- assumes that the orthogonality checks 
!	  have been performed
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers 
! ids           : 1D int, [i,j,m,k,l] ids 
! val           : real*8, value of the integral

SUBROUTINE cori_HO_O3_pup_ints(ndim,PsiL,PsiR,ids,val)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR,ids
  REAL(KIND=8), INTENT(INOUT) :: val
  INTEGER, INTENT(IN) :: ndim
  LOGICAL, DIMENSION(0:5-1,0:5-1) :: isq,isp
  INTEGER, DIMENSION(0:5-1,0:5-1) :: positions
  INTEGER, DIMENSION(0:5-1) :: values,counts
  REAL(KIND=8) :: temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 5
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0

!  WRITE(*,*) "int ids:",ids+7
 
  !Gather information on positions
  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2 .OR. i .EQ. 3) THEN
          isq(0,j) = .TRUE. 
          isp(0,j) = .FALSE.
        ELSE
          isq(0,j) = .FALSE. 
          isp(0,j) = .TRUE.
        END IF
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2 .OR. i .EQ. 3) THEN
          isq(counts(j),j) = .TRUE. 
          isp(counts(j),j) = .FALSE.
        ELSE
          isq(counts(j),j) = .FALSE. 
          isp(counts(j),j) = .TRUE.
        END IF
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO

  !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> P
  ! 2 -> Q
  ! 3 -> Q
  ! 4 -> P
  DO i=0,nids-1
    
    IF (counts(i) .EQ. 0) THEN
      EXIT
    
    ELSE IF (counts(i) .EQ. 1) THEN
      !q
      IF (isq(0,i)) THEN
        temp = ints_HO_q(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "q",values(i)
      !p
      ELSE IF (isp(0,i)) THEN
        temp = ints_HO_p(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "p",values(i)
      ELSE
        WRITE(*,*) "cori_HO_O3_pup_ints  : ERROR"
        WRITE(*,*) "There is a 1 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 2) THEN
      !qq
      IF (isq(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_qq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qq",values(i)
      !qp
      ELSE IF (isq(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_qp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qp",values(i)
      !pq
      ELSE IF (isp(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_pq(PsiL(values(i)),PsiR(values(i))) 
!        WRITE(*,*) "pq",values(i)
      !pp
      ELSE IF (isp(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_pp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pp",values(i)
      ELSE
        WRITE(*,*) "cori_HO_O3_pup_ints  : ERROR"
        WRITE(*,*) "There is a 2 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 3) THEN
      !qqq
      IF (isq(0,i) .AND. isq(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_qqq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqq",values(i)
      !qqp
      ELSE IF (isq(0,i) .AND. isq(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_qqp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqp",values(i)
      !qpq
      ELSE IF (isq(0,i) .AND. isp(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_qpq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qpq",values(i)
      !qpp
      ELSE IF (isq(0,i) .AND. isp(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_qpp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qpp",values(i)
      !pqq
      ELSE IF (isp(0,i) .AND. isq(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_pqq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pqq",values(i)
      !pqp
      ELSE IF (isp(0,i) .AND. isq(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_pqp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pqp",values(i)
      ELSE 
        WRITE(*,*) "cori_HO_O3_pup_ints  : ERROR"
        WRITE(*,*) "There is a 3 index case that has been missed"
        STOP
      END IF

    ELSE
      WRITE(*,*) "cori_HO_O3_pup_ints  : ERROR"
      WRITE(*,*) "There is a bad case here"
      STOP
    END IF
    val = val*temp
  END DO

  !accont for 1/i factor of the two p's
  val = -1.0D0*val

END SUBROUTINE cori_HO_O3_pup_ints 

!------------------------------------------------------------
! cori_HO_O3_upp_ints
!       - evalutes integrals of the kind 
!         <qi pj qm qk pl>
!	- assumes that the orthogonality checks 
!	  have been performed
!       - this is the μ.π.π case of the Watson Hamiltonian
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers 
! ids           : 1D int, [i,j,m,k,l] ids 
! val           : real*8, value of the integral

SUBROUTINE cori_HO_O3_upp_ints(ndim,PsiL,PsiR,ids,val)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR,ids
  REAL(KIND=8), INTENT(INOUT) :: val
  INTEGER, INTENT(IN) :: ndim
  LOGICAL, DIMENSION(0:5-1,0:5-1) :: isq,isp
  INTEGER, DIMENSION(0:5-1,0:5-1) :: positions
  INTEGER, DIMENSION(0:5-1) :: values,counts
  REAL(KIND=8) :: temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 5
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0

!  WRITE(*,*) "int ids:",ids+7
 
  !Gather information on positions
  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        IF (i .EQ. 0 .OR. i .EQ. 1 .OR. i .EQ. 3) THEN
          isq(0,j) = .TRUE. 
          isp(0,j) = .FALSE.
        ELSE
          isq(0,j) = .FALSE. 
          isp(0,j) = .TRUE.
        END IF
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        IF (i .EQ. 0 .OR. i .EQ. 1 .OR. i .EQ. 3) THEN
          isq(counts(j),j) = .TRUE. 
          isp(counts(j),j) = .FALSE.
        ELSE
          isq(counts(j),j) = .FALSE. 
          isp(counts(j),j) = .TRUE.
        END IF
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO

  !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> Q
  ! 2 -> P
  ! 3 -> Q
  ! 4 -> P
  DO i=0,nids-1
    
    IF (counts(i) .EQ. 0) THEN
      EXIT
    
    ELSE IF (counts(i) .EQ. 1) THEN
      !q
      IF (isq(0,i)) THEN
        temp = ints_HO_q(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "q",values(i)
      !p
      ELSE IF (isp(0,i)) THEN
        temp = ints_HO_p(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "p",values(i)
      ELSE
        WRITE(*,*) "cori_HO_O3_upp_ints  : ERROR"
        WRITE(*,*) "There is a 1 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 2) THEN
      !qq
      IF (isq(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_qq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qq",values(i)
      !qp
      ELSE IF (isq(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_qp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qp",values(i)
      !pq
      ELSE IF (isp(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_pq(PsiL(values(i)),PsiR(values(i))) 
!        WRITE(*,*) "pq",values(i)
      !pp
      ELSE IF (isp(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_pp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pp",values(i)
      ELSE
        WRITE(*,*) "cori_HO_O3_upp_ints  : ERROR"
        WRITE(*,*) "There is a 2 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 3) THEN
      !qqq
      IF (isq(0,i) .AND. isq(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_qqq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqq",values(i)
      !qqp
      ELSE IF (isq(0,i) .AND. isq(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_qqp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqp",values(i)
      !qpq
      ELSE IF (isq(0,i) .AND. isp(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_qpq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qpq",values(i)
      !qpp
      ELSE IF (isq(0,i) .AND. isp(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_qpp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qpp",values(i)
      ELSE 
        WRITE(*,*) "cori_HO_O3_upp_ints  : ERROR"
        WRITE(*,*) "There is a 3 index case that has been missed"
        STOP
      END IF

    ELSE
      WRITE(*,*) "cori_HO_O3_upp_ints  : ERROR"
      WRITE(*,*) "There is a bad case here"
      STOP
    END IF
    val = val*temp
  END DO

  !accont for 1/i factor of the two p's
  val = -1.0D0*val

END SUBROUTINE cori_HO_O3_upp_ints 

!------------------------------------------------------------
! cori_HO_O4_pup_ints
!       - evalutes integrals of the kind 
!         <qi pj qm qn qk pl>
!	- assumes that the orthogonality checks 
!	  have been performed
!       - this code is for the case of π.μ.π in the Watson
!         hamiltonian
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers 
! ids           : 1D int, [i,j,m,k,l] ids 
! val           : real*8, value of the integral

SUBROUTINE cori_HO_O4_pup_ints(ndim,PsiL,PsiR,ids,val)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR,ids
  REAL(KIND=8), INTENT(INOUT) :: val
  INTEGER, INTENT(IN) :: ndim
  LOGICAL, DIMENSION(0:6-1,0:6-1) :: isq,isp
  INTEGER, DIMENSION(0:6-1,0:6-1) :: positions
  INTEGER, DIMENSION(0:6-1) :: values,counts
  REAL(KIND=8) :: temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 6
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0

!  WRITE(*,*) "PsiL", PsiL
!  WRITE(*,*) "PsiR", PsiR
!  WRITE(*,*) "int ids",ids
 
  !Gather information on positions
  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2 .OR. i .EQ. 3 .OR. i .EQ. 4) THEN
          isq(0,j) = .TRUE. 
          isp(0,j) = .FALSE.
        ELSE
          isq(0,j) = .FALSE. 
          isp(0,j) = .TRUE.
        END IF
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        IF (i .EQ. 0 .OR. i .EQ. 2 .OR. i .EQ. 3 .OR. i .EQ. 4) THEN
          isq(counts(j),j) = .TRUE. 
          isp(counts(j),j) = .FALSE.
        ELSE
          isq(counts(j),j) = .FALSE. 
          isp(counts(j),j) = .TRUE.
        END IF
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO

  !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> P
  ! 2 -> Q
  ! 3 -> Q
  ! 4 -> Q
  ! 5 -> P
  DO i=0,nids-1

    IF (counts(i) .EQ. 0) THEN
      temp = 1.0D0
      !EXIT
    
    ELSE IF (counts(i) .EQ. 1) THEN
      !q
      IF (isq(0,i)) THEN
        temp = ints_HO_q(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "q",values(i)+1,PsiL(values(i)),temp
      !p
      ELSE IF (isp(0,i)) THEN
        temp = ints_HO_p(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "p",values(i)+1,PsiL(values(i)),temp
      ELSE
        WRITE(*,*) "cori_HO_O4_pup_ints  : ERROR"
        WRITE(*,*) "There is a 1 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 2) THEN
      !qq
      IF (isq(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_qq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qq",values(i)+1,PsiL(values(i)),temp
      !qp
      ELSE IF (isq(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_qp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qp",values(i)+1,PsiL(values(i)),temp
      !pq
      ELSE IF (isp(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_pq(PsiL(values(i)),PsiR(values(i))) 
!        WRITE(*,*) "pq",values(i)+1,PsiL(values(i)),temp
      !pp
      ELSE IF (isp(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_pp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pp",values(i)+1,PsiL(values(i)),temp
      ELSE
        WRITE(*,*) "cori_HO_O4_pup_ints  : ERROR"
        WRITE(*,*) "There is a 2 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 3) THEN
      !qqq
      IF (isq(0,i) .AND. isq(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_qqq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqq",values(i)+1,PsiL(values(i)),temp
      !qqp
      ELSE IF (isq(0,i) .AND. isq(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_qqp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqp",values(i)+1,PsiL(values(i)),temp
      !qpq
      ELSE IF (isq(0,i) .AND. isp(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_qpq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qpq",values(i)+1,PsiL(values(i)),temp
      !qpp
      ELSE IF (isq(0,i) .AND. isp(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_qpp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qpp",values(i)+1,PsiL(values(i)),temp
      !pqq
      ELSE IF (isp(0,i) .AND. isq(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_pqq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pqq",values(i)+1,PsiL(values(i)),temp
      !pqp
      ELSE IF (isp(0,i) .AND. isq(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_pqp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pqp",values(i)+1,PsiL(values(i)),temp
      ELSE 
        WRITE(*,*) "cori_HO_O4_pup_ints  : ERROR"
        WRITE(*,*) "There is a 3 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 4) THEN
      !qqqq
      IF (isq(0,i) .AND. isq(1,i) .AND. isq(2,i) .AND. isq(3,i)) THEN
        temp = ints_HO_qqqq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqqq",values(i)+1,PsiL(values(i)),temp
      !qqqp
      ELSE IF (isq(0,i) .AND. isq(1,i) .AND. isq(2,i) .AND. isp(3,i)) THEN
        temp = ints_HO_qqqp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqqp",values(i)+1,PsiL(values(i)),temp
      !pqqq
      ELSE IF (isp(0,i) .AND. isq(1,i) .AND. isq(2,i) .AND. isq(3,i)) THEN
        temp = ints_HO_pqqq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pqqq",values(i)+1,PsiL(values(i)),temp
      !pqqp
      ELSE IF (isp(0,i) .AND. isq(1,i) .AND. isq(2,i) .AND. isp(3,i)) THEN
        temp = ints_HO_pqqp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pqqp",values(i)+1,PsiL(values(i)),temp
      ELSE
        WRITE(*,*) "cori_HO_O4_pup_ints  : ERROR"
        WRITE(*,*) "There is a 4 index case that has been missed"
        STOP
      END IF

    ELSE
      WRITE(*,*) "cori_HO_O4_pup_ints  : ERROR"
      WRITE(*,*) "There is a bad case here"
      STOP
    END IF
    val = val*temp
  END DO

  !accont for 1/i factor of the two p's
  val = -1.0D0*val
!  write(*,*) "final int :", val

END SUBROUTINE cori_HO_O4_pup_ints 

!------------------------------------------------------------
! cori_HO_O4_upp_ints
!       - evalutes integrals of the kind 
!         <qi pj qm qn qk pl>
!	- assumes that the orthogonality checks 
!	  have been performed
!       - this code is for the case of μ.π.π in the Watson
!         hamiltonian
!------------------------------------------------------------
! ndim          : int, number of dimensions
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers 
! ids           : 1D int, [i,j,m,k,l] ids 
! val           : real*8, value of the integral

SUBROUTINE cori_HO_O4_upp_ints(ndim,PsiL,PsiR,ids,val)
  IMPLICIT NONE
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR,ids
  REAL(KIND=8), INTENT(INOUT) :: val
  INTEGER, INTENT(IN) :: ndim
  LOGICAL, DIMENSION(0:6-1,0:6-1) :: isq,isp
  INTEGER, DIMENSION(0:6-1,0:6-1) :: positions
  INTEGER, DIMENSION(0:6-1) :: values,counts
  REAL(KIND=8) :: temp
  INTEGER :: nids
  INTEGER :: i,j
  nids = 6
  counts = 0
  values = -1
  positions = -1
  val = 1.0D0

!  WRITE(*,*) "PsiL", PsiL
!  WRITE(*,*) "PsiR", PsiR
!  WRITE(*,*) "int ids",ids
 
  !Gather information on positions
  DO i=0,nids-1
    DO j=0,nids-1
      IF (counts(j) .EQ. 0) THEN !a new element
        values(j) = ids(i)
        positions(0,j) = i
        IF (i .EQ. 0 .OR. i .EQ. 1 .OR. i .EQ. 2 .OR. i .EQ. 4) THEN
          isq(0,j) = .TRUE. 
          isp(0,j) = .FALSE.
        ELSE
          isq(0,j) = .FALSE. 
          isp(0,j) = .TRUE.
        END IF
        counts(j) = 1
        EXIT
      ELSE IF (ids(i) .EQ. values(j)) THEN !an old element
        positions(counts(j),j) = i
        IF (i .EQ. 0 .OR. i .EQ. 1 .OR. i .EQ. 2 .OR. i .EQ. 4) THEN
          isq(counts(j),j) = .TRUE. 
          isp(counts(j),j) = .FALSE.
        ELSE
          isq(counts(j),j) = .FALSE. 
          isp(counts(j),j) = .TRUE.
        END IF
        counts(j) = counts(j) + 1
        EXIT
      END IF
    END DO
  END DO

  !Send off to the right subroutines
  ! positions (ordered list)
  ! 0 -> Q
  ! 1 -> Q
  ! 2 -> Q
  ! 3 -> P
  ! 4 -> Q
  ! 5 -> P
  DO i=0,nids-1

    IF (counts(i) .EQ. 0) THEN
      temp = 1.0D0
      !EXIT
    
    ELSE IF (counts(i) .EQ. 1) THEN
      !q
      IF (isq(0,i)) THEN
        temp = ints_HO_q(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "q",values(i)+1,PsiL(values(i)),temp
      !p
      ELSE IF (isp(0,i)) THEN
        temp = ints_HO_p(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "p",values(i)+1,PsiL(values(i)),temp
      ELSE
        WRITE(*,*) "cori_HO_O4_upp_ints  : ERROR"
        WRITE(*,*) "There is a 1 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 2) THEN
      !qq
      IF (isq(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_qq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qq",values(i)+1,PsiL(values(i)),temp
      !qp
      ELSE IF (isq(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_qp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qp",values(i)+1,PsiL(values(i)),temp
      !pq
      ELSE IF (isp(0,i) .AND. isq(1,i)) THEN
        temp = ints_HO_pq(PsiL(values(i)),PsiR(values(i))) 
!        WRITE(*,*) "pq",values(i)+1,PsiL(values(i)),temp
      !pp
      ELSE IF (isp(0,i) .AND. isp(1,i)) THEN
        temp = ints_HO_pp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "pp",values(i)+1,PsiL(values(i)),temp
      ELSE
        WRITE(*,*) "cori_HO_O4_upp_ints  : ERROR"
        WRITE(*,*) "There is a 2 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 3) THEN
      !qqq
      IF (isq(0,i) .AND. isq(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_qqq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqq",values(i)+1,PsiL(values(i)),temp
      !qqp
      ELSE IF (isq(0,i) .AND. isq(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_qqp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqp",values(i)+1,PsiL(values(i)),temp
      !qpq
      ELSE IF (isq(0,i) .AND. isp(1,i) .AND. isq(2,i)) THEN
        temp = ints_HO_qpq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qpq",values(i)+1,PsiL(values(i)),temp
      !qpp
      ELSE IF (isq(0,i) .AND. isp(1,i) .AND. isp(2,i)) THEN
        temp = ints_HO_qpp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qpp",values(i)+1,PsiL(values(i)),temp
      ELSE 
        WRITE(*,*) "cori_HO_O4_upp_ints  : ERROR"
        WRITE(*,*) "There is a 3 index case that has been missed"
        STOP
      END IF

    ELSE IF (counts(i) .EQ. 4) THEN
      !qqqq
      IF (isq(0,i) .AND. isq(1,i) .AND. isq(2,i) .AND. isq(3,i)) THEN
        temp = ints_HO_qqqq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqqq",values(i)+1,PsiL(values(i)),temp
      !qqqp
      ELSE IF (isq(0,i) .AND. isq(1,i) .AND. isq(2,i) .AND. isp(3,i)) THEN
        temp = ints_HO_qqqp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqqp",values(i)+1,PsiL(values(i)),temp
      !qqpq
      ELSE IF (isq(0,i) .AND. isq(1,i) .AND. isp(2,i) .AND. isq(3,i)) THEN
        temp = ints_HO_qqpq(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqpq",values(i)+1,PsiL(values(i)),temp
      !qqpp
      ELSE IF (isq(0,i) .AND. isq(1,i) .AND. isp(2,i) .AND. isp(3,i)) THEN
        temp = ints_HO_qqpp(PsiL(values(i)),PsiR(values(i)))
!        WRITE(*,*) "qqpp",values(i)+1,PsiL(values(i)),temp
      ELSE
        WRITE(*,*) "cori_HO_O4_upp_ints  : ERROR"
        WRITE(*,*) "There is a 4 index case that has been missed"
        STOP
      END IF

    ELSE
      WRITE(*,*) "cori_HO_O4_upp_ints  : ERROR"
      WRITE(*,*) "There is a bad case here"
      STOP
    END IF
    val = val*temp
  END DO

  !accont for 1/i factor of the two p's
  val = -1.0D0*val
!  write(*,*) "final int :", val

END SUBROUTINE cori_HO_O4_upp_ints 
!------------------------------------------------------------
! cori_HO_O2_eval
!       - evalute 2nd order coriolis contributions in the
!         harmonic oscillator basis
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nrota         : int, number of rotational terms
! rota          : 1D real*8, rotational terms
! omega         : 1D real*8, harmonic frequencies stored in
!                            order omega[0] -> omega for dim 0
! ncori         : 1D int, number of coriolis interactions per
!                         rotational mode
! qcori         : 2D int, coriolis quantum numbers for
!                         each rotational mode
! cori          : 2D real*8, coriolis zetas
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! cval          : real*8, value to output
! error         : int, error code

SUBROUTINE cori_HO_O2_eval(ndim,nrota,rota,omega,ncori,qcori,cori,&
                           PsiL,PsiR,cval,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: cori
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota,omega
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: qcori
  INTEGER, DIMENSION(0:), INTENT(IN) :: ncori,PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: cval
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nrota
  INTEGER, DIMENSION(0:3) :: v
  REAL(KIND=8) :: val,temp
  INTEGER :: b,c,a,i,j,k,l
  error = 0
  cval = 0.0D0

  DO a=0,nrota-1
    val = 0.0D0
    DO b=0,ncori(a)-1
      k = qcori(2*b,a)
      l = qcori(2*b+1,a)
      DO c=0,ncori(a)-1
        i = qcori(2*c,a)
        j = qcori(2*c+1,a)
        v = [i,j,k,l]
        CALL sort_int_ijkl(v)
        !check orthognality
        IF (ANY(PsiL(0:v(0)-1) .NE. PsiR(0:v(0)-1)) .OR. &
            ANY(PsiL(v(0)+1:v(1)-1) .NE. PsiR(v(0)+1:v(1)-1)) .OR. &
            ANY(PsiL(v(1)+1:v(2)-1) .NE. PsiR(v(1)+1:v(2)-1)) .OR. &
            ANY(PsiL(v(2)+1:v(3)-1) .NE. PsiR(v(2)+1:v(3)-1)) .OR. &
            ANY(PsiL(v(3)+1:ndim-1) .NE. PsiR(v(3)+1:ndim-1))) THEN
         CYCLE
        ELSE
          CALL cori_HO_O2_ints(ndim,PsiL,PsiR,[i,j,k,l],temp)
          !val = val + cori(c,a)*cori(b,a)*&
          !      SQRT(omega(j)*omega(l)/(omega(i)*omega(k)))*&
          !      temp 
          cval = cval + rota(a)*cori(c,a)*cori(b,a)*&
                SQRT(omega(j)*omega(l)/(omega(i)*omega(k)))*&
                temp 
          CALL cori_HO_O2_ints(ndim,PsiL,PsiR,[j,i,k,l],temp)
          !val = val - cori(c,a)*cori(b,a)*&
          !      SQRT(omega(i)*omega(l)/(omega(j)*omega(k)))*&
          !      temp 
          cval = cval - rota(a)*cori(c,a)*cori(b,a)*&
                SQRT(omega(i)*omega(l)/(omega(j)*omega(k)))*&
                temp 
          CALL cori_HO_O2_ints(ndim,PsiL,PsiR,[i,j,l,k],temp)
          !val = val - cori(c,a)*cori(b,a)*&
          !      SQRT(omega(j)*omega(k)/(omega(i)*omega(l)))*&
          !      temp 
          cval = cval - rota(a)*cori(c,a)*cori(b,a)*&
                SQRT(omega(j)*omega(k)/(omega(i)*omega(l)))*&
                temp 
          CALL cori_HO_O2_ints(ndim,PsiL,PsiR,[j,i,l,k],temp)
          !val = val + cori(c,a)*cori(b,a)*&
          !      SQRT(omega(i)*omega(k)/(omega(j)*omega(l)))*&
          !      temp 
          cval = cval + rota(a)*cori(c,a)*cori(b,a)*&
                SQRT(omega(i)*omega(k)/(omega(j)*omega(l)))*&
                temp 
        END IF
      END DO
    END DO
    !cval = cval + rota(a)*val
  END DO

END SUBROUTINE cori_HO_O2_eval


!------------------------------------------------------------
! cori_HO_O2_eval_debug
!       - evalutes second order coriolis contributions to the
!         Watson hamiltonian expressed in harmonic oscillator
!         basis
!       
!       - uses debug style arrays
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nrota         : int, number of non-inf rotational axis
! rota          : 1D int, non-inf rotational constants
! omega         : 1D real*8, HO frequencies
! cori_d        : 2D real*8, coriolis constants in debug form
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! val           : real*8, value
! error         : int, error code

SUBROUTINE cori_HO_O2_eval_debug(ndim,nrota,rota,omega,&
                                 cori_d,PsiL,PsiR,val,error)
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: cori_d 
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota,omega
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: val
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nrota
  INTEGER, DIMENSION(0:4-1) :: v
  REAL(KIND=8) :: temp,foo
  INTEGER :: i,j,k,l,a

  error = 0
  val = 0.0D0

  DO a=0,nrota-1
    DO i=0,ndim-1
      DO j=0,ndim-1
        IF (ABS(cori_d(i,j,a)) .LT. 1.0D-15) CYCLE
        !IF (i .EQ. j) CYCLE
        DO k=0,ndim-1
          DO l=0,ndim-1
        !    IF (k .EQ. l) CYCLE
            IF (ABS(cori_d(k,l,a)) .LT. 1.0D-15) CYCLE
            v = [i,j,k,l]
            CALL sort_int_ijkl(v)
            !check orthognality
            IF (ANY(PsiL(0:v(0)-1) .NE. PsiR(0:v(0)-1)) .OR. &
                ANY(PsiL(v(0)+1:v(1)-1) .NE. PsiR(v(0)+1:v(1)-1)) .OR. &
                ANY(PsiL(v(1)+1:v(2)-1) .NE. PsiR(v(1)+1:v(2)-1)) .OR. &
                ANY(PsiL(v(2)+1:v(3)-1) .NE. PsiR(v(2)+1:v(3)-1)) .OR. &
                ANY(PsiL(v(3)+1:ndim-1) .NE. PsiR(v(3)+1:ndim-1))) THEN
             CYCLE
            ELSE
              CALL cori_HO_O2_ints(ndim,PsiL,PsiR,[i,j,k,l],temp)
              foo = rota(a)*cori_d(i,j,a)*cori_d(k,l,a)*&
                    SQRT(omega(j)*omega(l)/(omega(i)*omega(k)))*&
                    temp
              val = val + foo
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE cori_HO_O2_eval_debug

!------------------------------------------------------------
! cori_HO_O3_eval
!	- evalutes third order coriolis contributions to the 
!	  Watson Hamiltonian expressed in the Harmonic 
! 	  Oscillator basis
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nrota         : int, number of rotational terms
! rota          : 1D real*8, rotational terms
! omega         : 1D real*8, harmonic frequencies stored in
!                            order omega[0] -> omega for dim 0
! ndidq		: 1D int, number of didq terms
! qdidq		: 2D int, QN of didq terms
! didq		: 2D real*8, didq terms
! ncori         : 1D int, number of coriolis interactions per
!                         rotational mode
! qcori         : 2D int, coriolis quantum numbers for
!                         each rotational mode
! cori          : 2D real*8, coriolis zetas
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! cval          : real*8, value to output
! error         : int, error code

SUBROUTINE cori_HO_O3_eval(ndim,nrota,rota,omega,ndidq,qdidq,didq,&
                           ncori,qcori,cori,PsiL,PsiR,cval,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: cori,didq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota,omega
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: qcori,qdidq
  INTEGER, DIMENSION(0:), INTENT(IN) :: ncori,PsiL,PsiR,ndidq
  REAL(KIND=8), INTENT(INOUT) :: cval
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nrota
  INTEGER, DIMENSION(0:4) :: v
  REAL(KIND=8) :: val,temp
  INTEGER :: b,c,a,i,j,k,l,m,ii,jj,kk,ll,mm
  error = 0
  cval = 0.0D0
  val = 0.0D0

  DO a=0,nrota-1
    DO b=0,nrota-1
      DO mm=0,ndidq(3*a+b)-1
        m = qdidq(mm,3*a+b)
        DO jj=0,ncori(b)-1
          k = qcori(2*jj,b)
          l = qcori(2*jj+1,b)
          DO ii=0,ncori(a)-1
            i = qcori(2*ii,a)
            j = qcori(2*ii+1,a)
            !orthogonality checks
            v = [i,j,m,k,l]
            CALL sort_int_ijklm(v)
            IF (ANY(PsiL(0:v(0)-1) .NE. PsiR(0:v(0)-1)) .OR.  &
                ANY(PsiL(v(0)+1:v(1)-1) .NE. PsiR(v(0)+1:v(1)-1)) .OR.  &
                ANY(PsiL(v(1)+1:v(2)-1) .NE. PsiR(v(1)+1:v(2)-1)) .OR.  &
                ANY(PsiL(v(2)+1:v(3)-1) .NE. PsiR(v(2)+1:v(3)-1)) .OR.  &
                ANY(PsiL(v(3)+1:v(4)-1) .NE. PsiR(v(3)+1:v(4)-1)) .OR.  &
                ANY(PsiL(v(4)+1:ndim-1) .NE. PsiR(v(4)+1:ndim-1))) THEN
              CYCLE
            ELSE
              CALL cori_HO_O3_pup_ints(ndim,PsiL,PsiR,[i,j,m,k,l],temp)
!              CALL cori_HO_O3_upp_ints(ndim,PsiL,PsiR,[m,i,j,k,l],temp)
              val = val - rota(a)*rota(b)*didq(mm,3*a+b)*&
                    SQRT(omega(j)*omega(l)/(omega(i)*omega(k)))*&
                    cori(ii,a)*cori(jj,b)*temp

              CALL cori_HO_O3_pup_ints(ndim,PsiL,PsiR,[j,i,m,k,l],temp)
!              CALL cori_HO_O3_upp_ints(ndim,PsiL,PsiR,[m,j,i,k,l],temp)
              val = val + rota(a)*rota(b)*didq(mm,3*a+b)*&
                    SQRT(omega(i)*omega(l)/(omega(j)*omega(k)))*&
                    cori(ii,a)*cori(jj,b)*temp

              CALL cori_HO_O3_pup_ints(ndim,PsiL,PsiR,[i,j,m,l,k],temp)
!              CALL cori_HO_O3_upp_ints(ndim,PsiL,PsiR,[m,i,j,l,k],temp)
              val = val + rota(a)*rota(b)*didq(mm,3*a+b)*&
                    SQRT(omega(j)*omega(k)/(omega(i)*omega(l)))*&
                    cori(ii,a)*cori(jj,b)*temp

              CALL cori_HO_O3_pup_ints(ndim,PsiL,PsiR,[j,i,m,l,k],temp)
!              CALL cori_HO_O3_upp_ints(ndim,PsiL,PsiR,[m,j,i,l,k],temp)
              val = val - rota(a)*rota(b)*didq(mm,3*a+b)*&
                    SQRT(omega(i)*omega(k)/(omega(j)*omega(l)))*&
                    cori(ii,a)*cori(jj,b)*temp
            END IF
          END DO
        END DO 
      END DO
    END DO
  END DO
  cval = val

END SUBROUTINE cori_HO_O3_eval

!------------------------------------------------------------
! cori_HO_O4_eval
!	- evalutes fourth order coriolis contributions to the 
!	  Watson Hamiltonian expressed in the Harmonic 
! 	  Oscillator basis
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nrota         : int, number of rotational terms
! rota          : 1D real*8, rotational terms
! omega         : 1D real*8, harmonic frequencies stored in
!                            order omega[0] -> omega for dim 0
! ndidq		: 1D int, number of didq terms
! qdidq		: 2D int, QN of didq terms
! didq		: 2D real*8, didq terms
! ncori         : 1D int, number of coriolis interactions per
!                         rotational mode
! qcori         : 2D int, coriolis quantum numbers for
!                         each rotational mode
! cori          : 2D real*8, coriolis zetas
! PsiL          : 1D int, LHS quantum numbers
! PsiR          : 1D int, RHS quantum numbers
! cval          : real*8, value to output
! error         : int, error code

SUBROUTINE cori_HO_O4_eval(ndim,nrota,rota,omega,ndidq,qdidq,didq,&
                           ncori,qcori,cori,PsiL,PsiR,cval,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: cori,didq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota,omega
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: qcori,qdidq
  INTEGER, DIMENSION(0:), INTENT(IN) :: ncori,PsiL,PsiR,ndidq
  REAL(KIND=8), INTENT(INOUT) :: cval
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nrota
  INTEGER, DIMENSION(0:5) :: v
  REAL(KIND=8) :: val,temp,foo
  INTEGER :: b,a,g,i,j,k,l,m,n,nn,ii,jj,mm
  error = 0
  cval = 0.0D0
  val = 0.0D0

  DO a=0,nrota-1
    DO b=0,nrota-1
      DO g=0,nrota-1
        DO mm=0,ndidq(3*a+g)-1
          m = qdidq(mm,3*a+g)
          DO nn=0,ndidq(3*g+b)-1
            n = qdidq(nn,3*g+b)
            DO jj=0,ncori(b)-1
              k = qcori(2*jj,b)
              l = qcori(2*jj+1,b)
              DO ii=0,ncori(a)-1
                i = qcori(2*ii,a)
                j = qcori(2*ii+1,a)
                !orthogonality checks
                v = [i,j,m,n,k,l]
                CALL sort_int_ijklmn(v)
                IF (ANY(PsiL(0:v(0)-1) .NE. PsiR(0:v(0)-1)) .OR.  &
                    ANY(PsiL(v(0)+1:v(1)-1) .NE. PsiR(v(0)+1:v(1)-1)) .OR.  &
                    ANY(PsiL(v(1)+1:v(2)-1) .NE. PsiR(v(1)+1:v(2)-1)) .OR.  &
                    ANY(PsiL(v(2)+1:v(3)-1) .NE. PsiR(v(2)+1:v(3)-1)) .OR.  &
                    ANY(PsiL(v(3)+1:v(4)-1) .NE. PsiR(v(3)+1:v(4)-1)) .OR.  &
                    ANY(PsiL(v(4)+1:v(5)-1) .NE. PsiR(v(4)+1:v(5)-1)) .OR.  &
                    ANY(PsiL(v(5)+1:ndim-1) .NE. PsiR(v(5)+1:ndim-1))) THEN
                  CYCLE
                ELSE
                  CALL cori_HO_O4_pup_ints(ndim,PsiL,PsiR,[i,j,m,n,k,l],temp)
!                  CALL cori_HO_O4_upp_ints(ndim,PsiL,PsiR,[m,n,i,j,k,l],temp)

                  val = val + rota(a)*rota(b)*rota(g)*&
                        didq(mm,3*a+g)*didq(nn,3*g+b)*&
                        SQRT(omega(j)*omega(l)/(omega(i)*omega(k)))*&
                        cori(ii,a)*cori(jj,b)*temp
!             IF (ABS(temp) .GT. 1.0D-15) THEN
!WRITE(*,'(1x,3(I1,1x),4x,3(I1,1x),4x,6(I1,1x),4x,F15.10,4x,F15.10,4x,F15.10)') PsiL(0:ndim-1),a+1,b+1,g+1,i+1,j+1,m+1,n+1,k+1,l+1,temp,foo,0.75D0*val
!             END IF

                  CALL cori_HO_O4_pup_ints(ndim,PsiL,PsiR,[j,i,m,n,k,l],temp)
!                  CALL cori_HO_O4_upp_ints(ndim,PsiL,PsiR,[m,n,j,i,k,l],temp)
                  val = val - rota(a)*rota(b)*rota(g)*&
                        didq(mm,3*a+g)*didq(nn,3*g+b)*&
                        SQRT(omega(i)*omega(l)/(omega(j)*omega(k)))*&
                        cori(ii,a)*cori(jj,b)*temp
!             IF (ABS(temp) .GT. 1.0D-15) THEN
!WRITE(*,'(1x,3(I1,1x),4x,3(I1,1x),4x,6(I1,1x),4x,F15.10,4x,F15.10,4x,F15.10)') PsiL(0:ndim-1),a+1,b+1,g+1,j+1,i+1,m+1,n+1,k+1,l+1,temp,foo,0.75D0*val
!             END IF

                  CALL cori_HO_O4_pup_ints(ndim,PsiL,PsiR,[i,j,m,n,l,k],temp)
!                  CALL cori_HO_O4_upp_ints(ndim,PsiL,PsiR,[m,n,i,j,l,k],temp)
                  val = val - rota(a)*rota(b)*rota(g)*&
                        didq(mm,3*a+g)*didq(nn,3*g+b)*&
                        SQRT(omega(j)*omega(k)/(omega(i)*omega(l)))*&
                        cori(ii,a)*cori(jj,b)*temp
!             IF (ABS(temp) .GT. 1.0D-15) THEN
!WRITE(*,'(1x,3(I1,1x),4x,3(I1,1x),4x,6(I1,1x),4x,F15.10,4x,F15.10,4x,F15.10)') PsiL(0:ndim-1),a+1,b+1,g+1,i+1,j+1,m+1,n+1,l+1,k+1,temp,foo,0.75D0*val
!             END IF

                  CALL cori_HO_O4_pup_ints(ndim,PsiL,PsiR,[j,i,m,n,l,k],temp)
!                  CALL cori_HO_O4_upp_ints(ndim,PsiL,PsiR,[m,n,j,i,l,k],temp)
                  val = val + rota(a)*rota(b)*rota(g)*&
                        didq(mm,3*a+g)*didq(nn,3*g+b)*&
                        SQRT(omega(i)*omega(k)/(omega(j)*omega(l)))*&
                        cori(ii,a)*cori(jj,b)*temp
!             IF (ABS(temp) .GT. 1.0D-15) THEN
!WRITE(*,'(1x,3(I1,1x),4x,3(I1,1x),4x,6(I1,1x),4x,F15.10,4x,F15.10,4x,F15.10)') PsiL(0:ndim-1),a+1,b+1,g+1,j+1,i+1,m+1,n+1,l+1,k+1,temp,foo,0.75D0*val
!             END IF
                END IF
              END DO
            END DO 
          END DO
        END DO
      END DO
    END DO
  END DO
  cval = 0.75D0*val
END SUBROUTINE cori_HO_O4_eval

!------------------------------------------------------------
! cori_HO_O4_eval_debug
!       - calculates 4th order coriolis contributions
!         to the hamiltonian in the harmonic oscillator 
!         basis using debugging arrays
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nrota         : int, number of non-inf rotational constants
! rota          : 1D real*8, non-info rotaitonal constants
! omega         : 1D real*8, frequencies
! didq_d        : 3D real*8, debugging didq
! cori_d        : 3D real*8, debugging cori
! PsiL          : 1D int, LHS vib quantum numbers
! PsiR          : 1D int, RHS vib quantum numbers
! val           : real*8, value to return
! error         : int, exit code

SUBROUTINE cori_HO_O4_eval_debug(ndim,nrota,rota,omega,didq_d,cori_d,&
                                 PsiL,PsiR,val,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: didq_d,cori_d 
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota,omega
  INTEGER, DIMENSION(0:), INTENT(IN) :: PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: val
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nrota
  INTEGER, DIMENSION(0:6-1) :: v
  REAL(KIND=8) :: temp,foo
  INTEGER :: i,j,k,l,m,n,a,b,g
  
  error = 0
  val = 0.0D0

  DO a=0,nrota-1
    DO b=0,nrota-1
      DO g=0,nrota-1
        DO i=0,ndim-1
          DO j=0,ndim-1
            IF (ABS(cori_d(i,j,a)) .LT. 1.0D-15) CYCLE 
            DO k=0,ndim-1
              DO l=0,ndim-1
                IF (ABS(cori_d(k,l,b)) .LT. 1.0D-15) CYCLE 
                DO m=0,ndim-1
                  IF (ABS(didq_d(m,a,g)) .LT. 1.0D-15) CYCLE
                  DO n=0,ndim-1
                    IF (ABS(didq_d(n,g,b)) .LT. 1.0D-15) CYCLE
                    v = [i,j,m,n,k,l]
                    CALL sort_int_ijklmn(v)
                    IF (ANY(PsiL(0:v(0)-1) .NE. PsiR(0:v(0)-1)) .OR.  &
                        ANY(PsiL(v(0)+1:v(1)-1) .NE. PsiR(v(0)+1:v(1)-1)) .OR.  &
                        ANY(PsiL(v(1)+1:v(2)-1) .NE. PsiR(v(1)+1:v(2)-1)) .OR.  &
                        ANY(PsiL(v(2)+1:v(3)-1) .NE. PsiR(v(2)+1:v(3)-1)) .OR.  &
                        ANY(PsiL(v(3)+1:v(4)-1) .NE. PsiR(v(3)+1:v(4)-1)) .OR.  &
                        ANY(PsiL(v(4)+1:v(5)-1) .NE. PsiR(v(4)+1:v(5)-1)) .OR.  &
                        ANY(PsiL(v(5)+1:ndim-1) .NE. PsiR(v(5)+1:ndim-1))) THEN
                      CYCLE
                    ELSE
                      !\pi \mu \pi order
                      CALL cori_HO_O4_pup_ints(ndim,PsiL,PsiR,[i,j,m,n,k,l],temp)
                      ! \mu \pi \pi order
!                      CALL cori_HO_O4_upp_ints(ndim,PsiL,PsiR,[m,n,i,j,k,l],temp)
                      foo = 0.75D0*rota(a)*rota(b)*rota(g)*&
                            didq_d(m,a,g)*didq_d(n,g,b)*&
                            cori_d(i,j,a)*cori_d(k,l,b)*&
                            SQRT(omega(j)*omega(l)/(omega(i)*omega(k)))*&
                            temp
                      val = val + foo
!                      IF (ABS(temp) .GT. 1.0D-15) THEN
WRITE(*,'(1x,3(I1,1x),4x,3(I1,1x),4x,6(I1,1x),4x,F15.10,4x,F15.10,4x,F15.10)') PsiL(0:ndim-1),a+1,b+1,g+1,i+1,j+1,m+1,n+1,k+1,l+1,temp,foo,val
!                      END IF

                    END IF
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE cori_HO_O4_eval_debug

END MODULE cori
!------------------------------------------------------------

