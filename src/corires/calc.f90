!------------------------------------------------------------
! calc
!       - module containing subroutines for calculating
!         Beff
!------------------------------------------------------------
MODULE calc
  USE conv
  USE ints
  USE term
  USE linal

CONTAINS

!------------------------------------------------------------
! calc_states
!       - calculates Beff for given states
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! voff          : int, vibrational mode numbering offset
! nstates       : int, number of states
! l2h           : 1D int, labeling -> harmonic numbering
! states        : 2D int, states to calculated (vQns, state)
! phi2          : 1D real*8, quadratic force constants
! phi3          : 3D real*8, cubic force constants
! Be            : 1D real*8, rotational constants
! zeta          : 3D real*8, coriolis zetas
! mu1           : 3D real*8, order 1 mu terms (vib,rot,rot)
! mu2           : 4D real*8, order 2 mu terms (vib,vib,rot,rot)
SUBROUTINE calc_states(nvib,voff,nstates,l2h,states,phi2,&
                                    phi3,Be,zeta,mu1,mu2,didq)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: mu2
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: phi3,zeta,mu1,didq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2,Be
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: states
  INTEGER, DIMENSION(0:), INTENT(IN) :: l2h
  INTEGER, INTENT(IN) :: nvib,voff,nstates
  REAL(KIND=8), DIMENSION(0:2,0:nstates-1) :: Beff
  INTEGER, DIMENSION(0:nvib-1) :: bra,ket,zero
  REAL(KIND=8), DIMENSION(0:2) :: Beff_0
  REAL(KIND=8) :: val
  INTEGER :: i,j,a,b,n

  !Be
  bra = 0
  ket = 0
  zero = 0
  Beff_0 = Be
  WRITE(*,*) "                                                 X        Y        Z"
  WRITE(*,*) "---------------------------------------------------------------------"
  WRITE(*,'(2x,A2)',ADVANCE='no') "Be"
  WRITE(*,'(8x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(0))
  WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(1))
  WRITE(*,'(2x,F18.7)') conv_cm2MHz(Beff_0(2))

!  TESTING TESTING TESTING
!  WRITE(*,*) "axis    mode   quadratic"
!  DO a=0,2
!    DO i=0,nvib-1
!      ket = 0
!      ket(i) = 1
!      WRITE(*,'(1x,2(I1,1x),2x,F12.7)') a+1,i+voff,term_1(nvib,a,ket,mu2)-term_1(nvib,a,zero,mu2)
!    END DO
!  END DO
!
!  WRITE(*,*) 
!  WRITE(*,*)
!  WRITE(*,*) "axis     mode        coriolis"
!  DO a=0,2
!    DO i=0,nvib-1
!      ket = 0
!      ket(i) = 1
!      WRITE(*,'(1x,2(I1,1x),2x,F12.7)') a+1,i+voff,term_2(nvib,a,ket,Be,phi2,zeta)&
!                            -term_2(nvib,a,zero,Be,phi2,zeta)
!    END DO
!  END DO
!
!  WRITE(*,*)
!  WRITE(*,*)
!  WRITE(*,*) "axis    mode      anharmonic"
!  DO a=0,2
!    DO i=0,nvib-1
!      ket = 0
!      ket(i) = 1
!      WRITE(*,'(1x,2(I1,1x),2x,F12.7)') a+1,i+voff,term_3(nvib,a,ket,phi2,phi3,mu1)-&
!                                        term_3(nvib,a,zero,phi2,phi3,mu1)
!    END DO
!  END DO
!  STOP

  !B0
  !term1
  ket = 0
  DO a=0,2
    Beff_0(a) = Beff_0(a) + term_1(nvib,a,ket,mu2)
  END DO

  !H(1)H(1) term
  DO a=0,2
    Beff_0(a) = Beff_0(a) + term_2(nvib,a,ket,Be,phi2,zeta)
  END DO

  !H(1)H(3) term
  DO a=0,2
    Beff_0(a) = Beff_0(a) + term_3(nvib,a,ket,phi2,phi3,mu1)
  END DO

  WRITE(*,'(2x,A2)',ADVANCE='no') "B0"
  WRITE(*,'(8x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(0))
  WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff_0(1))
  WRITE(*,'(2x,F18.7)') conv_cm2MHz(Beff_0(2))

  !Be - B0
  WRITE(*,'(2x,A7)',ADVANCE='no') "Be - B0"
  WRITE(*,'(3x,F18.7)',ADVANCE='no') conv_cm2MHz(Be(0) - Beff_0(0))
  WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Be(1) - Beff_0(1))
  WRITE(*,'(2x,F18.7)') conv_cm2MHz(Be(2) - Beff_0(2))

  WRITE(*,*)
  WRITE(*,*)

  IF (nstates .EQ. 0) RETURN

  !calculate Be's per state
  WRITE(*,*) " State                                           X        Y        Z"
  WRITE(*,*) "---------------------------------------------------------------------"
  DO n=0,nstates-1
 
    !get harmonic numbering
    DO j=0,nvib-1
      ket(l2h(j)-1) = states(j,n)
    END DO

    !order 0 term
    Beff(0:2,n) = Be

    !term1
    DO a=0,2
      Beff(a,n) = Beff(a,n) + term_1(nvib,a,ket,mu2)
    END DO

    !term 2
    DO a=0,2
      Beff(a,n) = Beff(a,n) + term_2(nvib,a,ket,Be,phi2,zeta)
    END DO

    !term 3
    DO a=0,2
      Beff(a,n) = Beff(a,n) + term_3(nvib,a,ket,phi2,phi3,mu1)
    END DO
    
    !print
    WRITE(*,'(1x,999(5(I2,1x),2x))',ADVANCE='no') states(0:nvib-1,n)
    WRITE(*,'(4x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(0,n))
    WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(1,n))
    WRITE(*,'(2x,F18.7)') conv_cm2MHz(Beff(2,n))
    
  END DO
  
  WRITE(*,*)
  WRITE(*,*)
  
  !calculate Beff transitions
  WRITE(*,*) " B_eff(v0) --> B_eff(vi)                         X        Y        Z"
  WRITE(*,*) "---------------------------------------------------------------------"
  DO n=0,nstates-1

    !get harmonic numbering
    DO j=0,nvib-1
      ket(l2h(j)-1) = states(j,n)
    END DO

    !Print
    WRITE(*,'(1x,999(5(I2,1x),2x))',ADVANCE='no') states(0:nvib-1,n)
    WRITE(*,'(4x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(0,n) - Beff_0(0))
    WRITE(*,'(2x,F18.7)',ADVANCE='no') conv_cm2MHz(Beff(1,n) - Beff_0(1))
    WRITE(*,'(2x,F18.7)') conv_cm2MHz(Beff(2,n) - Beff_0(2))

  END DO

  WRITE(*,*)
  WRITE(*,*)

END SUBROUTINE calc_states

!------------------------------------------------------------
! calc_diag
!       - diagonalizes effective hamiltonian for the 
!         specified states 
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! voff          : int, vibrational mode numbering offset
! nstates       : int, number of states
! l2h           : 1D int, labeling -> harmonic numbering
! states        : 2D int, states to calculated (vQns, state)
! phi2          : 1D real*8, quadratic force constants
! phi3          : 3D real*8, cubic force constants
! Be            : 1D real*8, rotational constants
! zeta          : 3D real*8, coriolis zetas
! mu1           : 3D real*8, order 1 mu terms (vib,rot,rot)
! mu2           : 4D real*8, order 2 mu terms (vib,vib,rot,rot)
! error         : int, error status
SUBROUTINE calc_diag(nvib,voff,nstates,l2h,states,phi2,&
                     phi3,Be,zeta,mu1,mu2,didq,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: mu2
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: phi3,zeta,mu1,didq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2,Be
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: states
  INTEGER, DIMENSION(0:), INTENT(IN) :: l2h
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nvib,voff,nstates
  INTEGER, DIMENSION(0:nvib-1,0:nstates-1) :: bra,ket
  REAL(KIND=8), DIMENSION(0:nstates-1,0:nstates-1,0:2) :: Beff
  REAL(KIND=8), DIMENSION(0:nstates-1,0:2) :: Eval
  REAL(KIND=8), DIMENSION(0:nstates-1) :: Ei
  REAL(KIND=8), DIMENSION(0:2) :: Beff_0,mval
  REAL(KIND=8) :: val1,val2
  CHARACTER(LEN=1), DIMENSION(0:2) :: xyz
  CHARACTER(LEN=100) :: label
  INTEGER, DIMENSION(0:2) :: mloc
  INTEGER :: i,j,k,a,b,n,m

  Beff = 0.0D0
  Beff_0 = Be
  xyz = ["X","Y","Z"]

  !B0
  WRITE(*,*) "                                                 X        Y        Z"
  WRITE(*,*) "---------------------------------------------------------------------"
  !term1
  ket = 0
  DO a=0,2
    Beff_0(a) = Beff_0(a) + term_1(nvib,a,ket(0:nvib-1,0),mu2)
  END DO

  !H(1)H(1) term
  DO a=0,2
    Beff_0(a) = Beff_0(a) + term_2(nvib,a,ket(0:nvib-1,0),Be,phi2,zeta)
  END DO

  !H(1)H(3) term
  DO a=0,2
    Beff_0(a) = Beff_0(a) + term_3(nvib,a,ket(0:nvib-1,0),phi2,phi3,mu1)
  END DO

  !Convert units
  DO a=0,2
    Beff_0(a) = conv_cm2MHz(Beff_0(a))
  END DO  

  WRITE(*,'(2x,A2)',ADVANCE='no') "B0"
  WRITE(*,'(8x,F18.7)',ADVANCE='no') Beff_0(0)
  WRITE(*,'(2x,F18.7)',ADVANCE='no') Beff_0(1)
  WRITE(*,'(2x,F18.7)') Beff_0(2)
  WRITE(*,*)
  WRITE(*,*)

  !Get kets
  DO n=0,nstates-1
    DO j=0,nvib-1
      ket(l2h(j)-1,n) = states(j,n)
      bra(l2h(j)-1,n) = states(j,n)
    END DO
  END DO

  !Get E0 energies
  DO j=0,nstates-1
    Ei(j) = 0.0D0
    DO i=0,nvib-1
      Ei(j) = Ei(j) + phi2(i)*(1.0D0*ket(i,j)+0.5D0)
    END DO
  END DO

  !Get initial Beffs
  DO n=0,nstates-1
    DO a=0,2
      Beff(n,n,a) = Be(a) 
      Beff(n,n,a) = Beff(n,n,a) + term_1(nvib,a,ket(0:nvib-1,n),mu2)  
      Beff(n,n,a) = Beff(n,n,a) + term_2(nvib,a,ket(0:nvib-1,n),Be,phi2,zeta)
      Beff(n,n,a) = Beff(n,n,a) + term_3(nvib,a,ket(0:nvib-1,n),phi2,phi3,mu1) 
    END DO 
  END DO

  !Convert units
  DO a=0,2
    DO n=0,nstates-1
      DO m=0,nstates-1
        Beff(n,m,a) = conv_cm2MHz(Beff(n,m,a))
      END DO
    END DO
  END DO

  WRITE(*,*) "Initial Beff Matrix"
  DO a=0,2
    WRITE(*,'(2x,A1,2x,A5)') xyz(a),"(MHz)"
    WRITE(*,*) "---------------------------------------------------------------------"
    DO n=0,nstates-1
      IF (n .LT. 10) THEN
        WRITE(label,'(A1,I1)') "s",n+1
      ELSE
        WRITE(label,'(A1,I2)') "s",n+1
      END IF
      WRITE(*,'(1x,A3)',ADVANCE='no') TRIM(label)
      WRITE(*,'(1x,999(F11.2,2x))') Beff(n,0:nstates-1,a)
    END DO
    WRITE(*,*) 
    WRITE(*,*)
  END DO

  WRITE(*,*)

  !Go through states
  DO a=0,2
    DO m=0,nstates-2
      DO n=m+1,nstates-1
        !find where they differ
        i = -1
        j = -1
        DO k=0,nvib-1
          IF (bra(k,n) - ket(k,m) .NE. 0) THEN
            IF (i .EQ. -1) THEN
              i = k
            ELSE IF (i .NE. -1 .AND. j .EQ. -1) THEN
              j = k
            ELSE
              i = -2
            END IF
          END IF
        END DO
!        WRITE(*,*) "<",bra(0:nvib-1,n),"|",ket(0:nvib-1,m),">"
        !if they do not differ by two indicies
        IF (i .EQ. -1 .OR. i .EQ. -2 .OR. j .EQ. -1) THEN
          val1 = 0.0D0
          val2 = 0.0D0
!          WRITE(*,*) "case 1"
        ELSE
!          WRITE(*,*) "case 2"
          val1 = 4.0D0*Be(a)*Be(a)*&
                 term_2_aux(i,j,a,bra(0:nvib-1,n),ket(0:nvib-1,m),phi2,zeta)/&
                (Ei(m) - Ei(n))
          val2 = 4.0D0*Be(a)*Be(a)*&
                 term_2_aux(i,j,a,bra(0:nvib-1,m),ket(0:nvib-1,n),phi2,zeta)/&
                 (Ei(n) - Ei(m))
          
        END IF
!        WRITE(*,*) i,j,&
!        term_2_aux(i,j,a,bra(0:nvib-1,n),ket(0:nvib-1,m),phi2,zeta),&
!        term_2_aux(i,j,a,bra(0:nvib-1,m),ket(0:nvib-1,n),phi2,zeta),&
!        (Ei(m)-Ei(n))

        !Put differences on off diagonal
        Beff(m,m,a) = Beff(m,m,a) - conv_cm2MHz(val1)
        Beff(n,n,a) = Beff(n,n,a) - conv_cm2MHz(val2)
        Beff(n,m,a) = 0.25D0*(mu2(i,j,a,a)*ints_q(bra(i,n),ket(i,m))*&
                                           ints_q(bra(j,n),ket(j,m))&
                            + mu2(j,i,a,a)*ints_q(bra(j,n),ket(j,m))*&
                                           ints_q(bra(i,n),ket(i,m)))
        Beff(n,m,a) = conv_cm2MHz(Beff(n,m,a))
        Beff(m,n,a) = Beff(n,m,a)
      END DO
    END DO
  END DO

  WRITE(*,*) "Effective Rotational Constant Matrix" 
  WRITE(*,*)
  DO a=0,2
    WRITE(*,'(2x,A1,2x,A5)') xyz(a),"(MHz)"
    WRITE(*,*) "---------------------------------------------------------------------"
    DO n=0,nstates-1
      IF (n .LT. 10) THEN
        WRITE(label,'(A1,I1)') "s",n+1
      ELSE
        WRITE(label,'(A1,I2)') "s",n+1
      END IF
      WRITE(*,'(1x,A3)',ADVANCE='no') TRIM(label)
      WRITE(*,'(1x,999(F11.2,2x))') Beff(n,0:nstates-1,a)
    END DO
    WRITE(*,*) 
    WRITE(*,*)
  END DO

  !Diagonalize
  DO a=0,2
    CALL linal_dsyev(nstates,Beff(0:nstates-1,0:nstates-1,a),&
                     3*nstates-1,Eval(0:nstates-1,a),error)   
    IF (error .NE. 0) RETURN
  END DO

  WRITE(*,*) 
  WRITE(*,*) "Diagonalized Rotational Constant Matrix"
  WRITE(*,*)
  DO a=0,2
    WRITE(*,'(2x,A1,1x,A5)',ADVANCE="no") xyz(a),"(MHz)"
    DO n=0,nstates-2
      WRITE(*,'(2x,F11.2)',ADVANCE="no") Eval(n,a)
    END DO
    WRITE(*,'(2x,F11.2)') Eval(nstates-1,a)
    WRITE(*,*) "---------------------------------------------------------------------"
    DO n=0,nstates-1
      IF (n .LT. 10) THEN
        WRITE(label,'(A1,I1)') "s",n+1
      ELSE
        WRITE(label,'(A1,I2)') "s",n+1
      END IF
      WRITE(*,'(4x,A3,2x)',ADVANCE='no') TRIM(label)
      WRITE(*,'(1x,999(F11.2,2x))') Beff(n,0:nstates-1,a)
    END DO
    WRITE(*,*)
    WRITE(*,*)
  END DO

  WRITE(*,*) "Dominant State Transitions"
  
  !calculate Beff transitions
  WRITE(*,*) " B_eff(v0) --> B_eff(vi)                         X        Y        Z"
  WRITE(*,*) "---------------------------------------------------------------------"
  DO n=0,nstates-1

    !Find the states we're looking for
    mloc = -1 
    mval = 0.0
    DO a=0,2
      DO m=0,nstates-1
        IF (ABS(Beff(m,n,a)) .GT. mval(a)) THEN
          mloc(a) = m
          mval(a) = ABS(Beff(m,n,a))
        END IF
      END DO
    END DO
    !WRITE(*,*) "mval:",mval
    !WRITE(*,*) "mloc:",mloc

    !Print
    WRITE(*,'(1x,999(5(I2,1x),2x))',ADVANCE='no') states(0:nvib-1,n)
    WRITE(*,'(4x,F18.7)',ADVANCE='no') Eval(mloc(0),0) - Beff_0(0)
    WRITE(*,'(2x,F18.7)',ADVANCE='no') Eval(mloc(1),1) - Beff_0(1)
    WRITE(*,'(2x,F18.7)') Eval(mloc(2),2) - Beff_0(2)

  END DO

  WRITE(*,*)
  WRITE(*,*)



END SUBROUTINE calc_diag

!------------------------------------------------------------
! calc_Heff
!       - calculates Heff between states
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! voff          : int, vibrational mode numbering offset
! nstates       : int, number of states
! l2h           : 1D int, labeling -> harmonic numbering
! states        : 2D int, states to calculated (vQns, state)
! phi2          : 1D real*8, quadratic force constants
! phi3          : 3D real*8, cubic force constants
! phi4          : 3D real*8, quartic force constants
! Be            : 1D real*8, rotational constants
! zeta          : 3D real*8, coriolis zetas
! mu0           : 2D real*8, order 0 mu terms (rot,rot)
! mu1           : 3D real*8, order 1 mu terms (vib,rot,rot)
! mu2           : 4D real*8, order 2 mu terms (vib,vib,rot,rot)
! error         : int, error status
 
SUBROUTINE calc_Heff(nvib,voff,nstates,l2h,states,phi2,&
                     phi3,phi4,Be,zeta,mu0,mu1,mu2,didq,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: mu2,phi4
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: phi3,zeta,mu1,didq
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: mu0
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2,Be
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: states
  INTEGER, DIMENSION(0:), INTENT(IN) :: l2h
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nvib,voff,nstates
  INTEGER, DIMENSION(0:nvib-1,0:nstates-1) :: bra,ket
  REAL(KIND=8), DIMENSION(0:nstates-1,0:nstates-1,0:2) :: Heff,H0,H1,H2
  REAL(KIND=8), DIMENSION(0:nstates-1,0:2) :: Eval
  REAL(KIND=8), DIMENSION(0:nstates-1) :: Ei
  REAL(KIND=8), DIMENSION(0:2) :: mval
  REAL(KIND=8) :: val1,val2,Ja,Jb
  CHARACTER(LEN=1), DIMENSION(0:2) :: xyz
  CHARACTER(LEN=100) :: label
  INTEGER, DIMENSION(0:2) :: mloc,Rstates
  INTEGER :: i,j,k,a,b,n,m

  error = 0
  xyz = ["X","Y","Z"]
  Heff = 0.0D0
  H0 = 0.0D0
  H1 = 0.0D0
  H2 = 0.0D0

  !Get kets
  DO n=0,nstates-1
    DO j=0,nvib-1
      ket(l2h(j)-1,n) = states(j,n)
      bra(l2h(j)-1,n) = states(j,n)
    END DO
  END DO

  !Read in Rotational Quantum Numbers
  WRITE(*,*) "Enter Rotational Quantum numbers (x,y,z)"
  READ(*,*) Rstates(0:2)

  !Construct Effective Hamiltonian between states
  DO a=0,2
    Ja = 1.0D0*Rstates(a)
    DO b=0,2
      Jb = 1.0D0*Rstates(b)
      DO m=0,nstates-1
        DO n=m,nstates-1
          !Order 0
          ! mu0
          H0(n,m,a) = H0(n,m,a) + 0.5D0*Ja*&
                      ints_mu0(nvib,a,b,mu0,bra(0:nvib-1,n),ket(0:nvib-1,m))&
                      *Jb 

          ! qr^2 + pr^2 
          H0(n,m,a) = H0(n,m,a) + &
                      (ints_V0(nvib,phi2,bra(0:nvib-1,n),ket(0:nvib-1,m)) +&
                       ints_T0(nvib,phi2,bra(0:nvib-1,n),ket(0:nvib-1,m)))

          !Order 1
          ! J u1 J
          H1(n,m,a) = H1(n,m,a) + 0.5D0*Ja*&
                      ints_mu1(nvib,a,b,mu1,bra(0:nvib-1,n),ket(0:nvib-1,m))*&
                      Jb

          ! J mu0 pi
          H1(n,m,a) = H1(n,m,a) - 0.5D0*Ja*&
                      ints_mu0pi(nvib,a,b,mu0,phi2,zeta,&
                                 bra(0:nvib-1,n),ket(0:nvib-1,m)) 

          ! J pi mu0
          H1(n,m,a) = H1(n,m,a) - 0.5D0*Jb*&
                      ints_mu0pi(nvib,b,a,mu0,phi2,zeta,&
                                 bra(0:nvib-1,n),ket(0:nvib-1,m)) 
          ! phi3
          H1(n,m,a) = H1(n,m,a) + ints_V1(nvib,phi3,&
                      bra(0:nvib-1,n),ket(0:nvib-1,m))

          !Order 2
          ! Ja mu2 Jb
          H2(n,m,a) = H2(n,m,a) + 0.5D0*Ja*&
                      ints_mu2(nvib,a,b,mu2,bra(0:nvib-1,n),ket(0:nvib-1,m))*&
                      Jb

          ! Ja mu1 pi
          H2(n,m,a) = H2(n,m,a) - 0.5D0*Ja*ints_mu1pi(nvib,a,b,mu1,phi2,zeta,&
                                  bra(0:nvib-1,n),ket(0:nvib-1,m))

          ! pi mu1 Jb
          H2(n,m,a) = H2(n,m,a) - 0.5D0*Jb*ints_pimu1(nvib,a,b,mu1,phi2,zeta,&
                                  bra(0:nvib-1,n),ket(0:nvib-1,m))

          ! pi mu0 pi
          H2(n,m,a) = H2(n,m,a) + 0.5D0*ints_pimu0pi(nvib,a,b,mu0,phi2,zeta,&
                                        bra(0:nvib-1,n),ket(0:nvib-1,m))

          ! phi4
          H2(n,m,a) = H2(n,m,a) + ints_V2(nvib,phi4,&
                                  bra(0:nvib-1,n),ket(0:nvib-1,m))

        END DO
      END DO
    END DO
  END DO

  Heff = H0 + H1 + H2

  !Write out Matrix
  WRITE(*,*) "Effective Hamiltonian (cm-1)"
  WRITE(*,*)
  DO a=0,2
    WRITE(*,'(2x,A1,2x,A5)') xyz(a),"(cm-1)"
    WRITE(*,*) "---------------------------------------------------------------------"
    DO n=0,nstates-1
      IF (n .LT. 10) THEN
        WRITE(label,'(A1,I1)') "s",n+1
      ELSE
        WRITE(label,'(A1,I2)') "s",n+1
      END IF
      WRITE(*,'(1x,A3)',ADVANCE='no') TRIM(label)
      WRITE(*,'(1x,999(F11.2,2x))') Heff(n,0:nstates-1,a)
    END DO
    WRITE(*,*) 
    WRITE(*,*)
  END DO
  
  
END SUBROUTINE calc_Heff
!------------------------------------------------------------

END MODULE calc
!------------------------------------------------------------
