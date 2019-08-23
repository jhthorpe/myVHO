!------------------------------------------------------------
! pseu
!	- module containing subroutines concerning the 
!	  Watson Pseudopotential
!------------------------------------------------------------
MODULE pseu
  USE sort
  USE ints_HO

CONTAINS

!------------------------------------------------------------
! pseu_HO_O2_eval 
!	- evalutes the contribution of the second order Watson 
!	  pseudopotential terms to the hamiltonian in the
!	  Harmonic Oscillator basis
!------------------------------------------------------------
! nrota         : int, number of rotational constants
! rota          : 1D real*8, rotational constants
! val           : real*8, the value to add

SUBROUTINE pseu_HO_O2_eval(nrota,rota,val)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota
  REAL(KIND=8), INTENT(INOUT) :: val
  INTEGER, INTENT(IN) :: nrota
  val = SUM(rota(0:nrota-1))
  val = -0.25D0*val
END SUBROUTINE pseu_HO_O2_eval

!------------------------------------------------------------
! pseu_HO_O3_eval
!	- evalutes the third order contributions of the 
!	  Watson Pseudopotential to the vibrational 
!	  Hamiltonian expressed in the Harmonic Oscillator 
! 	  basis
!
!	- The dI/dq constants are stored like:
! 	  didq(m,3*i+j) -> a^{i,j}_m
!------------------------------------------------------------
! ndim		: int, number of dimensions
! nrota		: int, number of rotational constants
! rota		: 1D real*8, rotational constants
! ndidq		: 1D int, number of didq constants
! qdidq		: 2D int, QN of didq constants
! didq		: 2D real*8, didq constants
! PsiL		: 1D int, QN of LHS wavefunction
! PsiR		: 1D int, QN of RHS wavefunction
! pval		: real*8, value to add

SUBROUTINE pseu_HO_O3_eval(ndim,nrota,rota,ndidq,qdidq,&
                           didq,PsiL,PsiR,pval)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: didq
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: qdidq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota
  INTEGER, DIMENSION(0:), INTENT(IN) :: ndidq,PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: pval
  INTEGER, INTENT(IN) :: ndim,nrota
  REAL(KIND=8) :: temp,val
  INTEGER :: i,j,a,v
  pval = 0.0D0
  val = 0.0D0
  DO a=0,nrota-1
!    val = 0.0D0
    j = 3*a+a
    DO i=0,ndidq(j)-1
      v = qdidq(i,j)
      !check orthogonality
      IF (ANY(PsiL(0:v-1) .NE. PsiR(0:v-1)) .OR. &
          ANY(PsiL(v+1:ndim-1) .NE. PsiR(v+1:ndim-1))) THEN
        CYCLE
      ELSE
        temp = ints_HO_q(PsiL(v),PsiR(v)) 
!        val = val + temp*didq(i,j)
        pval = pval + rota(a)*rota(a)*temp*didq(i,j)
      END IF
    END DO
!    pval = val + rota(a)*rota(a)*val
  END DO
  pval = 0.25D0*pval
END SUBROUTINE pseu_HO_O3_eval

!------------------------------------------------------------
! pseu_HO_O4_eval
!	- evalutes the fourth order contributions of the 
!	  Watson Pseudopotential to the vibrational 
!	  Hamiltonian expressed in the Harmonic Oscillator 
! 	  basis
!
!	- The dI/dq constants are stored like:
! 	  didq(m,3*i+j) -> a^{i,j}_m
!------------------------------------------------------------
! ndim		: int, number of dimensions
! nrota		: int, number of rotational constants
! rota		: 1D real*8, rotational constants
! ndidq		: 1D int, number of didq constants
! qdidq		: 2D int, QN of didq constants
! didq		: 2D real*8, didq constants
! PsiL		: 1D int, QN of LHS wavefunction
! PsiR		: 1D int, QN of RHS wavefunction
! pval		: real*8, value to add

SUBROUTINE pseu_HO_O4_eval(ndim,nrota,rota,ndidq,qdidq,didq,&
                           PsiL,PsiR,pval)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: didq
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: qdidq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: rota
  INTEGER, DIMENSION(0:), INTENT(IN) :: ndidq,PsiL,PsiR
  REAL(KIND=8), INTENT(INOUT) :: pval
  INTEGER, INTENT(IN) :: ndim,nrota
  INTEGER, DIMENSION(0:1) :: v
  REAL(KIND=8) :: temp
  INTEGER :: a,b,i,j,m,n
  pval = 0.0D0
  DO a=0,nrota-1
!    WRITE(*,*) "AAAA IS",a
    DO b=0,nrota-1
!      WRITE(*,*) "BBBB IS",b
      DO m=0,ndidq(3*a+b)-1
        i = qdidq(m,3*a+b)
        DO n=0,ndidq(3*b+a)-1
          j = qdidq(n,3*b+a)
    !      WRITE(*,*) "a,b :", a,b
    !      WRITE(*,*) "qi,qj  :", i,j
   !       WRITE(*,*) "Be^a,Be^b",rota(a),rota(b)
          !orthogonality checks
          v = [i,j]
          CALL sort_int_ij(v)
          IF (ANY(PsiL(0:v(0)-1) .NE. PsiR(0:v(0)-1)) .OR.&
              ANY(PsiL(v(0)+1:v(1)-1) .NE. PsiR(v(0)+1:v(1)-1)) .OR.&
              ANY(PsiL(v(1)+1:ndim-1) .NE. PsiR(v(1)+1:ndim-1))) THEN 
!            WRITE(*,*) "orthogonal"
            CYCLE
          ELSE
            IF (i .EQ. j) THEN 
              temp = ints_HO_qq(PsiL(i),PsiR(i))
            ELSE
              temp = ints_HO_q(PsiL(i),PsiR(i))*ints_HO_q(PsiL(j),PsiR(j))
            END IF
  !          WRITE(*,*) "ints  :",temp 
 !           WRITE(*,*) "ζ_i^ab,ζ_j^ba",didq(m,3*a+b),didq(n,3*b+a)
            pval = pval + rota(a)**2.0D0*rota(b)*didq(m,3*a+b)*didq(n,3*b+a)*temp
            
          END IF
!          WRITE(*,*) "---------------"
          
        END DO
      END DO
    END DO
  END DO
  pval = -3.0D0/16.0D0*pval
  
END SUBROUTINE pseu_HO_O4_eval

!------------------------------------------------------------


END MODULE pseu
!------------------------------------------------------------
