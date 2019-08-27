!------------------------------------------------------------
! H_HO
!       - module containing subroutines dealing with the 
!         Hamiltonian in harmonic oscillator basis
!------------------------------------------------------------
MODULE H_HO
  USE V
  USE fcon
  USE basis
  USE ints_HO
  USE key
  USE memory
  USE evec
  USE rota
  USE cori
  USE sort
  USE pseu
  USE didq

CONTAINS
!------------------------------------------------------------
! H_HO_build
!       - builds the Hamiltonian in the harmonic osciallator 
!         basis
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB
! nabs          : 1D int, number of abscissa
! q             : 2D real*8, list of abscissa
! W             : 2D real*8, list of weights
! Hij           : 2D real*8, hamiltonian
! error         : int, error code

SUBROUTINE H_HO_build(job,ndim,nbas,nabs,mem,q,W,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: q,W
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(IN) :: job,ndim
  INTEGER, INTENT(INOUT) :: error

  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: cori_d,didq_d
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Vij,Herm,cori,didq
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Vq,quad,cubi,quar,basK,rota
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: qcori,qdidq
  INTEGER, DIMENSION(:), ALLOCATABLE :: qquad,qcubi,qquar,ncori,ndidq
  REAL(KIND=8) :: qmem
  REAL(KIND=8) :: ti,tf
  INTEGER :: nQ1,nquad,ncubi,nquar,nrota
  INTEGER :: memstat
  INTEGER :: i,j,N,M

  INTEGER, DIMENSION(0:ndim-1) :: PsiL,PsiR
  INTEGER, DIMENSION(0:4) :: foo
  LOGICAL :: ex

  CALL CPU_TIME(ti)
  error = 0
  N = PRODUCT(nbas(0:ndim-1))
  WRITE(*,*) "Beginning Construction of HO Hamiltonian"
  WRITE(*,*) "Dimension of Hamiltonian",N
  WRITE(*,*)

  !analyze the memory situation 
  CALL memory_HObuild(job,mem,N,ndim,nbas,nabs,memstat,error) 
  IF (error .NE. 0) RETURN
  IF (memstat .EQ. 2) THEN
    ALLOCATE(Hij(0:N-1,0:N-1))
  ELSE 
    WRITE(*,*) "Only incore memory coded"
    error = 1
    RETURN
  END IF

  !If we're calculating V with gauss quad
  IF (job .EQ. 2 .OR. job .EQ. 3) THEN
    !Read in basis set information
    CALL basis_get(ndim,basK,error)
    IF (error .NE. 0) RETURN

    !Get potential energies at abscissa 
    IF (job .EQ. 2) THEN
      M = PRODUCT(nabs(0:ndim-1))
      ALLOCATE(Vij(0:M-1,0:ndim-1))
    ELSE IF (job .EQ. 3) THEN
      M = PRODUCT(nabs(0:ndim-1))
      ALLOCATE(Vq(0:M-1)) 
    END IF
    CALL V_get(job,ndim,nabs,Vij,Vq,error)
    IF (error .NE. 0) RETURN
  END IF
  

  !get constants
  CALL fcon_get(job,ndim,nquad,qquad,quad,ncubi,qcubi,cubi,&
                nquar,qquar,quar,error)
  IF (error .NE. 0) RETURN
  CALL rota_get(nrota,rota,error)
  IF (error .NE. 0) RETURN
  CALL cori_get(ndim,ncori,qcori,cori,error)  
  IF (error .NE. 0) RETURN
  CALL didq_get(ndim,ndidq,qdidq,didq,error)
  IF (error .NE. 0) RETURN
  INQUIRE(file='test',exist=ex)
  IF (ex) THEN
    CALL cori_get_debug(ndim,cori_d,error)
    IF (error .NE. 0) RETURN
    CALL didq_get_debug(ndim,didq_d,error)
    IF (error .NE. 0) RETURN
  END IF

  !build the hamiltonian
  IF (job .EQ. 1) THEN
    IF (memstat .EQ. 2) THEN
      CALL H_HO_build_poly_incore(ndim,nbas,nquad,qquad,quad,&
                                  ncubi,qcubi,cubi,&
                                  nquar,qquar,quar,nrota,rota,&
                                  ncori,qcori,cori,ndidq,qdidq,didq,&
                                  cori_d,didq_d,&
                                  Hij,error)
    END IF
  ELSE IF (job .EQ. 2) THEN
    IF (memstat .EQ. 2) THEN
      CALL H_HO_build_diag_incore(ndim,nbas,nabs,q,W,basK,&
                                  nquad,qquad,quad,ncubi,qcubi,cubi,&
                                  nquar,qquar,quar,nrota,rota,&
                                  ncori,qcori,cori,ndidq,qdidq,didq,&
                                  Vij,Hij,error)
    END IF
  ELSE IF (job .EQ. 3) THEN
    IF (memstat .EQ. 2) THEN
      CALL H_HO_build_quad_incore(ndim,nbas,nabs,q,W,basK,&
                                  nquad,qquad,quad,nrota,rota,&
                                  ncori,qcori,cori,ndidq,qdidq,didq,&
                                  Vq,Hij,error)
    END IF
  END IF


  IF (ALLOCATED(Vij)) DEALLOCATE(Vij)
  IF (ALLOCATED(Vq)) DEALLOCATE(Vq)
  IF (ALLOCATED(Herm)) DEALLOCATE(Herm)
  IF (ALLOCATED(rota)) DEALLOCATE(rota)

  CALL CPU_TIME(tf)
  WRITE(*,'(1x,A35,I2,A4,F6.1,1x,A1)') "H_HO_build  : finished with status ", error," in ", tf-ti, "s"
  WRITE(*,*)
  WRITE(*,*) "-----------------------------------------------------" 

END SUBROUTINE H_HO_build

!------------------------------------------------------------
! H_HO_Hermcalc
!       - calculates needed hermite polynomials
!       - stored [abscissa,basis qn, dimension]
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, number of basis functions
! nabs          : 1D int, number of abscissa
! q             : 2D real*8, abscissa [abscissa, dimension]
! Herm          : 3D real*8, hermite poly [abcs,basis qn,dimension]
! error         : int, exit code

SUBROUTINE H_HO_Hermcalc(ndim,nbas,nabs,q,Herm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: Herm
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: q
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=8) :: qval
  INTEGER :: i,j,k
  error = 0
  Herm = 0.0D0
  DO k=0,ndim-1
    DO j=0,nbas(k)-1
      IF (j .EQ. 0) THEN
        Herm(0:nabs(k)-1,0,k) = 1.0D0
      ELSE IF (j .EQ. 1) THEN
        Herm(0:nabs(k)-1,1,k) = 2.0D0*q(0:nabs(k)-1,k)
      ELSE
        Herm(0:nabs(k)-1,j,k) = 2.0D0*q(0:nabs(k)-1,k)*Herm(0:nabs(k)-1,j-1,k)&
                                -  2.0D0*(j-1)*Herm(0:nabs(k)-1,j-2,k)
      END IF
    END DO
  END DO
END SUBROUTINE H_HO_Hermcalc

!------------------------------------------------------------
! H_HO_Hermite_incore
!       - preconstruct all the Hermite polynomials needed
!       - outdated
!------------------------------------------------------------
! nabs          : int, nubmer of abscissa
! nbas          : 1D int, nubmer of basis functions
! q             : 1D real*8, abscissa
! Herm          : 2D real*8, hermite poly [order,abscissa]
! error         : int, exit code

SUBROUTINE H_HO_Hermite_incore(nabs,nbas,q,Herm,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Herm
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nabs
  
  REAL(KIND=8) :: qval
  INTEGER :: i,j,imax

  imax = MAXVAL(nbas)
  
  Herm = 0.0D0
  DO j=0,nabs-1
    qval = q(j)
    IF (imax .EQ. 1) THEN
      Herm(0,j) = 1.0D0
      CYCLE
    ELSE IF (imax .EQ. 2) THEN
      Herm(0,j) = 1.0D0
      Herm(1,j) = 2.0*qval
    ELSE
      Herm(0,j) = 1.0D0
      Herm(1,j) = 2.0*qval
      DO i=1,imax-2
        Herm(i+1,j) = 2.0D0*qval*Herm(i,j) - 2.0D0*i*Herm(i-1,j)
      END DO
    END IF
  END DO

END SUBROUTINE H_HO_Hermite_incore

!------------------------------------------------------------
! H_HO_build_poly_incore
!       - build hamiltonian where the potential is 
!         approximated as a 4th order polynomial in the 
!         harmonic oscillator basis set
!       - all terms are held incore
!       
!       PsiL and PsiR are the quantum nubmers of the left/right
!            wavefunctions
!
!       Integrals are stored like:
!       Q1(i,k) -> <i+1|Qk|i>
!       Q2(2*i,k) -> <i|Qk^2|i>    , Q2(2*i+1,k) -> <i+2|Qk^2|i>
!       Q3(2*i,k) -> <i+1|Qk^3|i>  , Q3(2*i+1,k) -> <i+3|Qk^3|i>
!       Q4(3*i,k) -> <i|Qk^4|i>    , Q4(3*i+1,k) -> <i+2|Qk^4|i>
!                                  , Q4(3*i+2,k) -> <i+4|Qk^4|i>
!       P2(2*i,k) -> <i|Pk^2|i>    , P2(2*i,k)   -> <i+2|Pk^2|i>

!       
!       Where i indicates the i'th quantum number of the 
!         k'th dimension
!
!       The integral
!          < i0 j0 k0 | Φ_122 | i1 j0 k0 >
!       Is then
!           1/6 * Φ_122 * <i1|Qi|i0> * <j0|Qj^2|j0> * <k0|k0>
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! nquad         : int, number of quadratic force constants
! qquad         : 1D int, QN of quadratic force constants
! quad          : 1D real*8, quadratic force constants
! ncubi         : int, number of cubic force constants
! qcubi         : 1D int, QN of cubic force constants
! cubi          : 1D real*8, cubic force constants
! nquar         : int, number of quartic force constants
! qquar         : 1D int, QN of quartic force constants
! quar          : 1D real*8, quartic force constants
! nrota         : int, number of rotational constants
! rota          : 1D real*8, rotational constants
! ncori         : 1D int, number of coriolis zeta
! qcori         : 2D int, coriolis zeta quantum numbers
! cori          : 2D real*8, coriolis zetas
! ndim          : int, number of dimensions
! ndidq         : 1D int, number of didq
! qdidq         : 2D int, quantum numbers of didq
! cori          : 3D real*8, coriolis constants in debugging
! didq          : 3D real*8, didq constants in debugging
! Hij		: 2D real*8, Hamiltonian
! error         : int, exit code

SUBROUTINE H_HO_build_poly_incore(ndim,nbas,nquad,qquad,quad,&
                                  ncubi,qcubi,cubi,&
                                  nquar,qquar,quar,nrota,rota,&
                                  ncori,qcori,cori,ndidq,qdidq,didq,&
                                  cori_d,didq_d,&
                                  Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: cori_d,didq_d
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: cori,didq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: quad,cubi,quar,rota 
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: qcori,qdidq
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER, DIMENSION(0:), INTENT(IN) :: qquad,qcubi,qquar,ncori,ndidq
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nquad,ncubi,nquar,nrota
  INTEGER, INTENT(IN) :: ndim

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Q1int,Q2int,Q3int,Q4int,&
                                               P1int,P2int,QPint,PQint
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: eps
  REAL(KIND=8), DIMENSION(0:ndim-1) :: omega
  INTEGER, DIMENSION(0:ndim-1) :: PsiL,PsiR,keyR,keyL,numL
  REAL(KIND=8) :: val,p2val,p3val,p4val,c2val,c3val,c4val,temp
  LOGICAL :: ex
  INTEGER :: N,M,mbas,il,ir
  INTEGER :: i,j,k,l,a

  error = 0
  mbas = MAXVAL(nbas)
  Hij = 0.0D0

  !Sort the quadratic force constants
  CALL sort_dirty_1Dint_1Dreal8(ndim,nquad,qquad,quad,omega,error)

  ALLOCATE(Q1int(0:mbas-1,0:ndim-1))
  ALLOCATE(Q2int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(Q3int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(Q4int(0:3*mbas-1,0:ndim-1))
  ALLOCATE(P1int(0:mbas-1,0:ndim-1))
  ALLOCATE(P2int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(QPint(0:2*mbas-1,0:ndim-1))
  ALLOCATE(PQint(0:2*mbas-1,0:ndim-1))

  !precalculate the needed integrals
  !we could do this on the fly, but this
  !is probably faster
  WRITE(*,*) "Calculating Integrals..."
  CALL ints_HO_Q1calc(ndim,nbas,Q1int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q2calc(ndim,nbas,Q2int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q3calc(ndim,nbas,Q3int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q4calc(ndim,nbas,Q4int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_P1calc(ndim,nbas,P1int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_P2calc(ndim,nbas,P2int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_QPcalc(ndim,nbas,QPint,error)
  IF (error .NE. 0) RETURN
  CALL ints_HO_PQcalc(ndim,nbas,PQint,error)
  IF (error .NE. 0) RETURN 

  ! the combination of key_generate and key_idx2ids
  ! allows us to transform what could be coded as a 
  ! recursive call into two unrolled loops
  ! The ket loop needs to go over all states
  !   to be lower triangular
  CALL key_generate(ndim,nbas,keyR) 
  N = PRODUCT(nbas(0:ndim-1))
  DO j=0,N-1
    CALL key_idx2ids(ndim,j,nbas,keyR,PsiR)

    !This is a fun little bit of code. It ensures that all
    ! indices on LHS are >= their partners on the RHS. 
    ! Additionally, for this polynomial case, only 
    ! elements within +/- 4 of the RHS can contribute. 
    numL = MIN(nbas - PsiR,5) 
    M = PRODUCT(numL(0:ndim-1)) ! number of elements in bra loop 
    CALL key_generate(ndim,numL,keyL)

    DO k=0,M-1
      CALL key_idx2ids(ndim,k,numL,keyL,PsiL)
      PsiL = PsiR + PsiL

      !evaluate force constants
      CALL ints_HO_polyput(ndim,PsiL,PsiR,nquad,qquad,quad,&
                           ncubi,qcubi,cubi,nquar,qquar,quar,Q1int,Q2int,&
                           Q3int,Q4int,P2int,val,error)       

      !fill in all permutations other than first index
      !we ignore the first index because lower triangular
      !this is doing way more work than we need to to, actually. We 
      !only need to calculate one of these permutations to get the 
      !rest, so this is a terrible way of doing things
      CALL key_ids2idx(ndim,nbas,keyR,PsiL,i) 
      DO l=0,2**(ndim-1)-1
        il = i
        ir = j
        DO a=1,ndim-1
          IF (MOD(l/2**(a-1),2) .EQ. 1) THEN
            il = il - (PsiL(a)-PsiR(a))*keyR(a)
            ir = ir + (PsiL(a)-PsiR(a))*keyR(a)
          END IF    
        END DO
        IF (il .GE. ir) Hij(il,ir) = val
!        Hij(il,ir) = val
      END DO
    END DO 
  END DO 
  WRITE(*,*)


!=============================================
! TESTING ZONE
INQUIRE(file='test',exist=ex)
IF (ex) THEN
  ALLOCATE(eps(0:N-1))
  DO i=0,N-1
    eps(i) = Hij(i,i)
  END DO
  Hij = 0.0D0
  DO i=0,N-1
    Hij(i,i) = eps(i)
  END DO

  c2val = 0.0D0
  c3val = 0.0D0
  c4val = 0.0D0
  p2val = 0.0D0
  p3val = 0.0D0
  p4val = 0.0D0

  WRITE(*,*) "TESTING TESTING TESTING"
  WRITE(*,*) "Enter PsiR:"
  READ(*,*) PsiR
  CALL key_ids2idx(ndim,nbas,keyR,PsiR,k)
  CALL key_idx2ids(ndim,k,nbas,keyR,PsiR)
  WRITE(*,*) "PsiR is", PsiR
  temp = Hij(k,k)
  WRITE(*,*) "Before ", temp

  WRITE(*,*) "Omega is", omega

  !Before diagonal
  DO i=0,k-1
    CALL key_idx2ids(ndim,i,nbas,keyR,PsiL)
    IF (ABS(eps(k) - eps(i)) .LT. 1.0D-15) CYCLE
!    CALL cori_HO_O2_eval(ndim,nrota,rota,omega,ncori,qcori,cori,&
!                         PsiR,PsiL,val,error)
!    CALL cori_HO_O2_eval(ndim,nrota,rota,omega,ncori,qcori,cori,&
!                         PsiL,PsiR,c2val,error)
    CALL cori_HO_O2_eval_debug(ndim,nrota,rota,omega,cori_d,&
                         PsiR,PsiL,val,error)
    CALL cori_HO_O2_eval_debug(ndim,nrota,rota,omega,cori_d,&
                         PsiL,PsiR,c2val,error)
    temp = temp + c2val*val/(eps(k)-eps(i))
    Hij(i,k) = Hij(i,k) + c2val*val/(eps(k)-eps(i))
  END DO

  !diagonal
  PsiL = PsiR
  i = k

  c2val = 0.0D0
  c3val = 0.0D0
  c4val = 0.0D0
  p2val = 0.0D0
  p3val = 0.0D0
  p4val = 0.0D0

  CALL pseu_HO_O2_eval(nrota,rota,p2val)
!  CALL pseu_HO_O3_eval(ndim,nrota,rota,ndidq,qdidq,didq,&
!                       PsiL,PsiR,p3val)
!  CALL pseu_HO_O4_eval(ndim,nrota,rota,ndidq,qdidq,didq,&
!                       PsiL,PsiR,p4val)
  CALL pseu_HO_O4_eval_debug(ndim,nrota,rota,didq_d,&
                             PsiL,PsiR,p4val,error)

!  CALL cori_HO_O2_eval(ndim,nrota,rota,omega,ncori,qcori,cori,&
!                       PsiL,PsiR,c2val,error)
  CALL cori_HO_O2_eval_debug(ndim,nrota,rota,omega,cori_d,&
                       PsiL,PsiR,c2val,error)
  !CALL cori_HO_O3_eval(ndim,nrota,rota,omega,ndidq,qdidq,didq,&
  !                     ncori,qcori,cori,PsiL,PsiR,c3val,error)
!  CALL cori_HO_O4_eval(ndim,nrota,rota,omega,ndidq,qdidq,didq,&
!                       ncori,qcori,cori,PsiL,PsiR,c4val,error)
  CALL cori_HO_O4_eval_debug(ndim,nrota,rota,omega,didq_d,cori_d,&
                             PsiL,PsiR,c4val,error)
  temp = temp + p2val + p3val + p4val + c2val + c3val + c4val
  Hij(i,k) = Hij(i,k) + p2val + p3val + p4val + c2val + c3val + c4val
  
  !After diagonal
  DO i=k+1,N-1
    IF (ABS(eps(k) - eps(i)) .LT. 1.0D-15) CYCLE
    CALL key_idx2ids(ndim,i,nbas,keyR,PsiL)
!    CALL cori_HO_O2_eval(ndim,nrota,rota,omega,ncori,qcori,cori,&
!                         PsiR,PsiL,val,error)
!    CALL cori_HO_O2_eval(ndim,nrota,rota,omega,ncori,qcori,cori,&
!                         PsiL,PsiR,c2val,error)
    CALL cori_HO_O2_eval_debug(ndim,nrota,rota,omega,cori_d,&
                         PsiR,PsiL,val,error)
    CALL cori_HO_O2_eval_debug(ndim,nrota,rota,omega,cori_d,&
                         PsiL,PsiR,c2val,error)
    temp = temp + c2val*val/(eps(k)-eps(i))
    Hij(i,k) = Hij(i,k) + c2val*val/(eps(k)-eps(i))
  END DO

  WRITE(*,*) "PsiR",PsiR
  WRITE(*,*) "VPT4 value :", temp
  WRITE(*,*) 
  WRITE(*,*)
  WRITE(*,*)

  INQUIRE(file='print',EXIST=ex) 
  IF (ex) THEN
    WRITE(*,*)
    WRITE(*,*) "PsiR",PsiR
    DO j=0,N-1
      CALL key_idx2ids(ndim,j,nbas,keyR,PsiL)
      WRITE(*,'(2x,F20.4,2x,A1,3(I4,1x),2x,A1,F20.4))') Hij(j,k),"|",PsiL,"|",eps(j)
    END DO
  END IF
  STOP
END IF
!=============================================
! END TESTING ZONE

  !rotational and coriolis parts
  CALL pseu_HO_O2_eval(nrota,rota,p2val)
  CALL key_generate(ndim,nbas,keyR) 
  DO j=0,N-1
    CALL key_idx2ids(ndim,j,nbas,keyR,PsiR)

    DO i=j,N-1
      CALL key_idx2ids(ndim,i,nbas,keyR,PsiL)

      !evaluate Pseduopotential terms
      ! O(2)
      IF (i .EQ. j) THEN
        Hij(i,j) = Hij(i,j) + p2val
      END IF
    
      ! O(3)
      IF (ALL(PsiL .LE. PsiR+1) .AND. ALL(PsiL .GE. PsiR-1)) THEN  
        CALL pseu_HO_O3_eval(ndim,nrota,rota,ndidq,qdidq,didq,&
                             PsiL,PsiR,p3val) 
        Hij(i,j) = Hij(i,j) + p3val
      END IF
      
      ! O(4) 
      IF (ALL(PsiL .LE. PsiR+2) .AND. ALL(PsiL .GE. PsiR-2)) THEN  
        CALL pseu_HO_O4_eval(ndim,nrota,rota,ndidq,qdidq,didq,&
                             PsiL,PsiR,p4val) 
        Hij(i,j) = Hij(i,j) + p4val
      END IF

      !evalutate Coriolis terms
      ! O(2)
      IF (ALL(PsiL .GE. PsiR-4) .AND. ALL(PsiL .LE. PsiR+4)) THEN
        CALL cori_HO_O2_eval(ndim,nrota,rota,omega,ncori,qcori,cori,&
                             PsiL,PsiR,c2val,error)
        IF (error .NE. 0) RETURN
        Hij(i,j) = Hij(i,j) + c2val
      END IF

      ! O(3)
      IF (ALL(PsiL .GE. PsiR-5) .AND. ALL(PsiL .LE. PsiR+5)) THEN
        CALL cori_HO_O3_eval(ndim,nrota,rota,omega,ndidq,qdidq,didq,&
                             ncori,qcori,cori,PsiL,PsiR,c3val,error)
        IF (error .NE. 0) RETURN
        Hij(i,j) = Hij(i,j) + c3val
      END IF

      ! O(4)
      IF (ALL(PsiL .GE. PsiR-6) .AND. ALL(PsiL .LE. PsiR+6)) THEN
        CALL cori_HO_O4_eval(ndim,nrota,rota,omega,ndidq,qdidq,didq,&
                             ncori,qcori,cori,PsiL,PsiR,c4val,error)
        IF (error .NE. 0) RETURN
        Hij(i,j) = Hij(i,j) + c4val
      END IF
    END DO
  END DO

  INQUIRE(file='print',EXIST=ex) 
  IF (ex) THEN
    WRITE(*,*)
    DO j=0,N-1
      WRITE(*,'(2x,999(F20.4,2x))') Hij(j,0:N-1)
    END DO
  END IF

  DEALLOCATE(Q1int)
  DEALLOCATE(Q2int)
  DEALLOCATE(Q3int)
  DEALLOCATE(Q4int)
  DEALLOCATE(P1int)
  DEALLOCATE(P2int)
  DEALLOCATE(QPint)
  DEALLOCATE(PQint)

END SUBROUTINE H_HO_build_poly_incore

!------------------------------------------------------------
! H_HO_build_diag_incore
!       - constructs hamiltonian where diagonal potential 
!         terms are calcualted via gaussian quadrature, 
!         and cubic/quartic force constants are calcualted 
!         via Harmonic oscillator integrals
!       - all memory is contained incore
!
!       PsiL and PsiR are the quantum nubmers of the left/right
!            wavefunctions
!
!       Integrals are stored like:
!       VT((j-i)+(keyI(i,k),k) -> <j|Vk + Tk|i>
!       Q1(i,k) -> <i+1|Qk|i>
!       Q2(2*i,k) -> <i|Qk^2|i>    , Q2(2*i+1,k) -> <i+2|Qk^2|i>
!       Q3(2*i,k) -> <i+1|Qk^3|i>  , Q3(2*i+1,k) -> <i+3|Qk^3|i>
!       Q4(3*i,k) -> <i|Qk^4|i>    , Q4(3*i+1,k) -> <i+2|Qk^4|i>
!                                  , Q4(3*i+2,k) -> <i+4|Qk^4|i>
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, number of basis functions per dimension
! nabs          : 1D int, number of abscissa per dimensions
! q             : 2D real*8, abscissa   [abscissa,dimension]
! W             : 2D real*8, weights    [abscissa,dimension]
! basK          : 1D real*8, basis function frequencies
! nquad         : int, number of quadratic force constants
! qquad         : 1D int, QNs of quadratic force constants
! quad          : 1D real*8, quadratic force constants
! ncubi         : int, number of cubic force constants
! qcubi         : 1D int, QNs of cubic force constants
! cubi          : 1D real*8, cubic force constants
! nquar         : int, number of quartic force constants
! qquar         : 1D int, QNs of quartic force constants
! quar          : 1D real*8, quartic force constants 
! nrota         : int, number of rotational constants
! rota          : 1D real*8, rotational constnants 
! ncori         : 1D int, number of coriolis zeta
! qcori         : 2D int, coriolis zeta quantum numbers
! cori          : 2D real*8, coriolis zetas
! ndim          : int, number of dimensions
! ndidq         : 1D int, number of didq
! qdidq         : 2D int, quantum numbers of didq
! Vij           : 2D real*8, potential energy   [abscissa,dimension]
! Hij           : 2D real*8, hamiltonian
! error         : int, exit code

SUBROUTINE H_HO_build_diag_incore(ndim,nbas,nabs,q,W,basK,&
                                  nquad,qquad,quad,ncubi,qcubi,cubi,&
                                  nquar,qquar,quar,nrota,rota,&
                                  ncori,qcori,cori,ndidq,qdidq,didq,&
                                  Vij,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: Vij,q,W,cori,didq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: basK,quad,cubi,quar,rota
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: qcori,qdidq
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs,qquad,qcubi,qquar,ncori,ndidq
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nquad,ncubi,nquar,nrota

  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: Herm
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: norm,VTint,Q1int,&
                                               Q2int,Q3int,Q4int,&
                                               P1int,P2int,QPint,PQint
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: keyI
  REAL(KIND=8), DIMENSION(0:ndim-1) :: omega
  INTEGER, DIMENSION(0:ndim-1) :: keyL,keyR,PsiL,PsiR,numL
  INTEGER, DIMENSION(0:ndim-2) :: keyX,PsiX,nelm
  REAL(KIND=8) :: val,p2val,c2val
  LOGICAL :: ex
  INTEGER :: N,M,mbas,mabs
  INTEGER :: i,j,k,l,a,il,ir

  error = 0
  Hij = 0.0D0
  mbas = MAXVAL(nbas)
  mabs = MAXVAL(nabs)

  !Sort the quadratic force constants
  CALL sort_dirty_1Dint_1Dreal8(ndim,nquad,qquad,quad,omega,error)

  !force constant integrals
  WRITE(*,*) "Calculating coupling integrals"
  ALLOCATE(Q1int(0:mbas-1,0:ndim-1))
  ALLOCATE(Q2int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(Q3int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(Q4int(0:3*mbas-1,0:ndim-1))
  ALLOCATE(P1int(0:mbas-1,0:ndim-1))
  ALLOCATE(P2int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(QPint(0:2*mbas-1,0:ndim-1))
  ALLOCATE(PQint(0:2*mbas-1,0:ndim-1))

  CALL ints_HO_Q1calc(ndim,nbas,Q1int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q2calc(ndim,nbas,Q2int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q3calc(ndim,nbas,Q3int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q4calc(ndim,nbas,Q4int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_P1calc(ndim,nbas,P1int,error)
  IF (error .NE. 0) RETURN
  CALL ints_HO_P2calc(ndim,nbas,P2int,error)
  IF (error .NE. 0) RETURN
  CALL ints_HO_QPcalc(ndim,nbas,QPint,error)
  IF (error .NE. 0) RETURN
  CALL ints_HO_PQcalc(ndim,nbas,PQint,error)
  IF (error .NE. 0) RETURN

  !Fill in the force constant integrals
  ! the combination of key_generate and key_idx2ids
  ! allows us to transform what could be coded as a 
  ! recursive call into two unrolled loops
  ! The ket loop needs to go over all states
  !   to be lower triangular
  !CALL key_generate(ndim,nbas,keyR) 
  CALL key_generate(ndim,nbas,keyR)
  N = PRODUCT(nbas(0:ndim-1))
  !DO j=0,-1
  DO j=0,N-1
    CALL key_idx2ids(ndim,j,nbas,keyR,PsiR)

    !This is a fun little bit of code. It ensures that all
    ! indices on LHS are >= their partners on the RHS. 
    numL = MIN(nbas - PsiR,5)
    M = PRODUCT(numL(0:ndim-1)) ! number of elements in bra loop 
    CALL key_generate(ndim,numL,keyL)

    DO k=0,M-1
      CALL key_idx2ids(ndim,k,numL,keyL,PsiL)
      PsiL = PsiR + PsiL
   
     !IF (ALL(PsiL .EQ. PsiR)) CYCLE

      !evaluate force constants
      CALL ints_HO_diagput(ndim,PsiL,PsiR,nquad,qquad,quad,&
                           ncubi,qcubi,cubi,nquar,qquar,quar,Q1int,Q2int,&
                           Q3int,Q4int,val,error)       
     ! WRITE(*,*) "val = ", val

      !fill in all permutations other than first index
      !we ignore the first index because lower triangular
      CALL key_ids2idx(ndim,nbas,keyR,PsiL,i) 
      CALL key_idx2ids(ndim,i,nbas,keyR,PsiL)
      DO l=0,2**(ndim-1)-1
        il = i
        ir = j
        DO a=1,ndim-1
          IF (MOD(l/2**(a-1),2) .EQ. 1) THEN
            il = il - (PsiL(a)-PsiR(a))*keyR(a)
            ir = ir + (PsiL(a)-PsiR(a))*keyR(a)
            !WRITE(*,*) "il,ir", il,ir
          END IF    
        END DO
        IF (il .GE. ir) Hij(il,ir) = val
      END DO
    END DO 
  END DO 

  !Rotational and Coriolis terms
  CALL pseu_HO_O2_eval(nrota,rota,val)
  CALL key_generate(ndim,nbas,keyR)
  DO j=0,N-1
    CALL key_idx2ids(ndim,j,nbas,keyR,PsiR)
    DO i=j,N-1
      CALL key_idx2ids(ndim,i,nbas,keyR,PsiL)

      !evalutate Coriolis terms
      IF (ANY(PsiL .LT. PsiR-5) .OR. ANY(PsiL .GT. PsiR+5)) THEN
        c2val = 0.0D0
      ELSE
        CALL cori_HO_O2_eval(ndim,nrota,rota,omega,ncori,qcori,cori,&
                              PsiL,PsiR,c2val,error)
        IF (error .NE. 0) RETURN
        Hij(i,j) = Hij(i,j) + c2val
      END IF
      !evaluate rotational terms
      IF (i .EQ. j) Hij(i,j) = Hij(i,j) + p2val

    END DO
  END DO

  DEALLOCATE(Q1int)
  DEALLOCATE(Q2int)
  DEALLOCATE(Q3int)
  DEALLOCATE(Q4int)
  DEALLOCATE(P1int)
  DEALLOCATE(P2int)
  DEALLOCATE(QPint)
  DEALLOCATE(PQint)

  INQUIRE(file='print',EXIST=ex) 
  IF (ex) THEN
    WRITE(*,*)
    DO j=0,N-1
      WRITE(*,'(2x,999(F20.4,2x))') Hij(j,0:N-1)
    END DO
  END IF

  !generate key for intermediates
  ALLOCATE(keyI(0:mbas-1,0:ndim-1))
  DO k=0,ndim-1
    keyI(0,k) = 0
    DO j=1,nbas(k)-1
      keyI(j,k) = keyI(j-1,k) + nbas(k) - (j-1)
    END DO
  END DO 

  !precalculate integrals
  ! hermite polynomials, normalization, potential and kinetic terms
  WRITE(*,*) "Calculating V+T integrals"
  ALLOCATE(Herm(0:mabs-1,0:mbas-1,0:ndim-1))
  ALLOCATE(VTint(0:(mbas*(mbas+1))/2-1,0:ndim-1))
  ALLOCATE(norm(0:mbas-1,0:ndim-1))
  CALL H_HO_Hermcalc(ndim,nbas,nabs,q,Herm,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_normcalc(ndim,nbas,nabs,W,Herm,norm,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_VTcalc(ndim,nbas,nabs,q,W,basK,norm,Herm,keyI,Vij,VTint,error)
  IF (error .NE. 0) RETURN 
  DEALLOCATE(norm)
  DEALLOCATE(Herm)

  !Place the VT integrals 
  !Special case for 1 dimension
  IF (ndim .EQ. 1) THEN
    N = PRODUCT(nbas(0:ndim-1))
    DO j=0,N-1
      DO i=j,N-1
        Hij(i,j) = Hij(i,j) + VTint(i-j+keyI(j,0),0)
      END DO
    END DO
  ELSE 
    !The only nonzero terms are <x...k..z|Vk|x..k'...z>
    !where k and k' are the only ones allowed to differ
    CALL key_generate(ndim,nbas,keyR)
    N = PRODUCT(nbas(0:ndim-1))
    ! loop over all dimesions
    DO k=0,ndim-1
    !DO k=0,-1
      !WRITE(*,*) "Dimension :", k
      nelm(0:k-1) = nbas(0:k-1)
      nelm(k:ndim-2) = nbas(k+1:ndim-1) 
      CALL key_generate(ndim-1,nelm,keyX)
      M = PRODUCT(nelm(0:ndim-2))

      !Loop over all <x...z|x...z>
      !x>z ???
      DO a=0,M-1
        CALL key_idx2ids(ndim-1,a,nelm,keyX,PsiX)
        !WRITE(*,*) "PsiX",PsiX

        !loop over |k'>
        DO j=0,nbas(k)-1
          PsiR(0:k-1) = PsiX(0:k-1)
          PsiR(k) = j
          PsiR(k+1:ndim-1) = PsiX(k:ndim-2) 
          CALL key_ids2idx(ndim,nbas,keyR,PsiR,ir)
         ! WRITE(*,*) "PsiR: ", PsiR

          !loop over <k|
          DO i=j,nbas(k)-1
            PsiL(0:k-1) = PsiR(0:k-1)
            PsiL(k) = i
            PsiL(k+1:ndim-1) = PsiR(k+1:ndim-1)
            CALL key_ids2idx(ndim,nbas,keyR,PsiL,il)
          !  WRITE(*,*) "PsiL: ", PsiL
            
            val = VTint(i-j+keyI(j,k),k)
            Hij(il,ir) = Hij(il,ir) + val
            !WRITE(*,*) "il,ir",il,ir
            IF (i .NE. j .AND. k .NE. 0) THEN
              !WRITE(*,*) "ir",ir
              !WRITE(*,*) "il",il
            !  il = il - (i-j)*keyR(k)
            !  ir = ir + (i-j)*keyR(k)
              Hij(il-(i-j)*keyR(k),ir+(i-j)*keyR(k)) = Hij(il-(i-j)*keyR(k),ir+(i-j)*keyR(k)) + val
            END IF
            !fill in other permuations of x<->z??
          END DO
          !WRITE(*,*) "=========="
        END DO
      END DO
    END DO
  END IF

  DEALLOCATE(VTint)

  
  INQUIRE(file='print',EXIST=ex) 
  IF (ex) THEN
    WRITE(*,*)
    DO j=0,N-1
      WRITE(*,'(2x,999(F20.4,2x))') Hij(j,0:N-1)
    END DO
  END IF

  WRITE(*,*)

  IF(ALLOCATED(keyI)) DEALLOCATE(keyI)
  IF(ALLOCATED(VTint)) DEALLOCATE(VTint)
  IF(ALLOCATED(Q1int)) DEALLOCATE(Q1int)
  IF(ALLOCATED(Q2int)) DEALLOCATE(Q2int)
  IF(ALLOCATED(Q3int)) DEALLOCATE(Q3int)
  IF(ALLOCATED(Q4int)) DEALLOCATE(Q4int)

  !WRITE(*,*)
  !DO j=0,N-1
  !  WRITE(*,'(2x,999(F20.4,2x))') Hij(j,0:N-1)
  !END DO

END SUBROUTINE H_HO_build_diag_incore

!------------------------------------------------------------
! H_HO_build_quad_incore
!       - constructs the hamiltonian with full dimensional 
!         integrals evaluated via Gaussian quadrature, 
!         and everything held in memory
!
!       - note that the potential should be already stored 
!         in the order in which it needs to be accessed
!
!       - the potential is stored assuming the last dimension
!         is the innermost loop
!
!       Improvements: 
!       1) Parallelization with OpenMP and/or MPI
!       2) The coriolis terms have symmetry realtionships
!            that could be used, and this has been ignored
!            for now. 
!       3) the full integrals (with v,v-2  v,v  v,v+2) 
!          are currently stored. Again, this could be 
!          made more efficient
!------------------------------------------------------------
! ndim          : int, nubmer of dimensions
! nbas          : 1D int, number of basis functions
! nabs          : 1D int, number of abscissa
! q             : 2D real*8, abscissa           [abscissa,dimension]
! W             : 2D real*8, weights            [weight,dimensions]
! basK          : 1D real*8, basis function FC
! nquad         : int, number of quadratic force constants
! qquad         : 1D int, QN's of quadratic force constants
! quad          : 1D real*8, quadratic force constants
! nrota         : int, number of rotational constants
! rota          : 1D real*8, rotational constants
! ncori         : 1D int, number of coriolis zeta
! qcori         : 2D int, coriolis zeta quantum numbers
! cori          : 2D real*8, coriolis zetas
! ndim          : int, number of dimensions
! ndidq         : 1D int, number of didq
! qdidq         : 2D int, quantum numbers of didq
! Vq            : 1D real*8, potential at abscissa
! Hij           : 2D real*8, the Hamiltonian
! error         : int, exit code 

SUBROUTINE H_HO_build_quad_incore(ndim,nbas,nabs,q,W,basK,&
                                  nquad,qquad,quad,nrota,rota,&
                                  ncori,qcori,cori,ndidq,qdidq,didq,&
                                  Vq,Hij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: q,W,cori,didq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: basK,Vq,rota,quad
  INTEGER, DIMENSION(0:,0:), INTENT(IN) :: qcori,qdidq
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas,nabs,ncori,qquad,ndidq
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,nrota,nquad
  
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: Herm
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Norm,Heff,Q1int,Q2int,&
                                               P1int,P2int,QPint,PQint
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Neff
  REAL(KIND=8), DIMENSION(0:ndim-1) :: omega
  INTEGER, DIMENSION(0:ndim-1) :: key,PsiR,PsiL
  REAL(KIND=8) :: c2val,c3val,c4val,vval,p2val,p3val,p4val
  LOGICAL :: ex
  INTEGER :: i,j,k,N,M,mbas,mabs
  
  error = 0
  WRITE(*,*) "Constructing the Hamiltonian"
  WRITE(*,*) 
  N = PRODUCT(nbas(0:ndim-1))
  M = PRODUCT(nabs(0:ndim-1))
  mbas = MAXVAL(nbas(0:ndim-1))
  mabs = MAXVAL(nabs(0:ndim-1))
  Hij = 0.0D0

  WRITE(*,*) "WARNING WARNING WARNING"
  WRITE(*,*) "Orthogonality checks are wrong"
  WRITE(*,*) "This is probably fixed for coriolis?"
!  STOP
  WRITE(*,*) "WARNING WARNING WARNING"

  !Sort the quadratic force constants
  CALL sort_dirty_1Dint_1Dreal8(ndim,nquad,qquad,quad,omega,error)

  !Allocate Space for needed integrals
  ALLOCATE(Q1int(0:mbas-1,0:ndim-1))
  ALLOCATE(Q2int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(P1int(0:mbas-1,0:ndim-1))
  ALLOCATE(P2int(0:2*mbas-1,0:ndim-1))
  ALLOCATE(QPint(0:2*mbas-1,0:ndim-1))
  ALLOCATE(PQint(0:2*mbas-1,0:ndim-1))
  ALLOCATE(Herm(0:mabs-1,0:mbas-1,0:ndim-1))
  ALLOCATE(Norm(0:mbas-1,0:ndim-1))
  ALLOCATE(Heff(0:mabs-1,0:2*ndim-1))
  ALLOCATE(Neff(0:2*ndim-1))
  Norm = 0.0D0

  !Evaluate the integrals
  CALL ints_HO_Q1calc(ndim,nbas,Q1int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_Q2calc(ndim,nbas,Q2int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_P1calc(ndim,nbas,P1int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_P2calc(ndim,nbas,P2int,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_QPcalc(ndim,nbas,QPint,error)
  IF (error .NE. 0) RETURN 
  CALL ints_HO_PQcalc(ndim,nbas,PQint,error)
  IF (error .NE. 0) RETURN 
  CALL H_HO_Hermcalc(ndim,nbas,nabs,q,Herm,error)
  IF (error .NE. 0) RETURN
  CALL ints_HO_normcalc(ndim,nbas,nabs,W,Herm,norm,error)
  IF (error .NE. 0) RETURN 

  !Fill in potential, kinetic, rotational, and coriolis
  CALL key_generate(ndim,nbas,key)
  CALL pseu_HO_O2_eval(nrota,rota,p2val)

  DO j=0,N-1
    CALL key_idx2ids(ndim,j,nbas,key,PsiR)
    DO k=0,ndim-1
      Heff(0:nabs(k)-1,2*k+1) = Herm(0:nabs(k)-1,PsiR(k),k)
      Neff(2*k+1) = Norm(PsiR(k),k)
    END DO

    DO i=j,N-1
      CALL key_idx2ids(ndim,i,nbas,key,PsiL)
      DO k=0,ndim-1
        Heff(0:nabs(k)-1,2*k) = Herm(0:nabs(k)-1,PsiL(k),k)
        Neff(2*k) = Norm(PsiL(k),k)
      END DO

      !evaluate potential and kinetic terms
      CALL ints_HO_quadput(ndim,nabs,q,W,basK,Neff,Heff,PsiL,PsiR,&
                           Vq,vval,error)
      IF (error .NE. 0) RETURN
      Hij(i,j) = Hij(i,j) + vval

      !evaluate Pseduopotential terms
      ! O(2)
      IF (i .EQ. j) THEN
        Hij(i,j) = Hij(i,j) + p2val
      END IF
    
      ! O(3)
      IF (ALL(PsiL .LE. PsiR+1) .AND. ALL(PsiL .GE. PsiR-1)) THEN  
        CALL pseu_HO_O3_eval(ndim,nrota,rota,ndidq,qdidq,didq,&
                             PsiL,PsiR,p3val) 
        Hij(i,j) = Hij(i,j) + p3val
      END IF
      
      ! O(4) 
      IF (ALL(PsiL .LE. PsiR+2) .AND. ALL(PsiL .GE. PsiR-2)) THEN  
        CALL pseu_HO_O4_eval(ndim,nrota,rota,ndidq,qdidq,didq,&
                             PsiL,PsiR,p4val) 
        Hij(i,j) = Hij(i,j) + p4val
      END IF
      

      !evalutate Coriolis terms
      ! O(2)
      IF (ALL(PsiL .GE. PsiR-4) .AND. ALL(PsiL .LE. PsiR+4)) THEN
        CALL cori_HO_O2_eval(ndim,nrota,rota,omega,ncori,qcori,cori,&
                             PsiL,PsiR,c2val,error)
        IF (error .NE. 0) RETURN
        Hij(i,j) = Hij(i,j) + c2val
      END IF

      ! O(3)
      IF (ALL(PsiL .GE. PsiR-5) .AND. ALL(PsiL .LE. PsiR+5)) THEN
        CALL cori_HO_O3_eval(ndim,nrota,rota,omega,ndidq,qdidq,didq,&
                             ncori,qcori,cori,PsiL,PsiR,c3val,error)
        IF (error .NE. 0) RETURN
        Hij(i,j) = Hij(i,j) + c3val
      END IF

      ! O(4)
      IF (ALL(PsiL .GE. PsiR-6) .AND. ALL(PsiL .LE. PsiR+6)) THEN
        CALL cori_HO_O4_eval(ndim,nrota,rota,omega,ndidq,qdidq,didq,&
                             ncori,qcori,cori,PsiL,PsiR,c4val,error)
        IF (error .NE. 0) RETURN
        Hij(i,j) = Hij(i,j) + c4val
      END IF

    END DO
  END DO

  IF (ALLOCATED(Norm)) DEALLOCATE(Norm) 
  IF (ALLOCATED(Herm)) DEALLOCATE(Herm) 
  IF (ALLOCATED(Heff)) DEALLOCATE(Heff) 
  IF (ALLOCATED(Neff)) DEALLOCATE(Neff) 
  IF (ALLOCATED(Herm)) DEALLOCATE(Herm) 
  IF (ALLOCATED(Q1int)) DEALLOCATE(Q1int)
  IF (ALLOCATED(Q2int)) DEALLOCATE(Q2int)
  IF (ALLOCATED(P1int)) DEALLOCATE(P1int)
  IF (ALLOCATED(P2int)) DEALLOCATE(P2int)
  IF (ALLOCATED(QPint)) DEALLOCATE(QPint)
  IF (ALLOCATED(PQint)) DEALLOCATE(PQint)

  INQUIRE(file='print',EXIST=ex) 
  IF (ex) THEN
    WRITE(*,*)
    DO j=0,N-1
      WRITE(*,'(2x,999(F20.4,2x))') Hij(j,0:N-1)
    END DO
  END IF

END SUBROUTINE H_HO_build_quad_incore
!------------------------------------------------------------
! H_HO_diag
!       - diagonalizes the hamiltonian in the HO basis
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nbas          : 1D int, basis functions per dimension
! enum          : int, number of eigenvalues to calc
! mem           : int*8, memory in MB
! Hij           : 2D real*8, hamiltonian
! eval          : 1D real*8, array of eigenvalues
! Cij           : 2D real*8, coefficients 
! error         : int, error 

SUBROUTINE H_HO_diag(ndim,nbas,enum,mem,Hij,eval,Cij,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Cij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: eval
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,enum

  REAL(KIND=8) :: ti,tf
  INTEGER :: memstat,lwork,N
  INTEGER :: i
  
  CALL CPU_TIME(ti)
  error = 0
  N = PRODUCT(nbas(0:ndim-1))

  CALL memory_Hdiag(ndim,nbas,enum,mem,memstat,lwork,error)
  IF (error .NE. 0) RETURN
  IF (memstat .NE. 2) THEN
    WRITE(*,*) "H_diag  : ERROR"
    WRITE(*,*) "Sorry, only incore diagonalization has been coded"
    error = 1
    RETURN
  END IF
  lwork = MAX(10,lwork)

  ALLOCATE(eval(0:enum-1))
  ALLOCATE(Cij(0:N-1,0:enum-1))
 
  WRITE(*,*) "Diagonalizing the Hamiltonian"
  WRITE(*,*)
  CALL linal_dsyevx(N,enum,lwork,Hij,eval,Cij,error)
  CALL CPU_TIME(tf)
  WRITE(*,'(1x,A27,1x,F7.1,A2)') "Hamiltonian diagonalized in",tf-ti," s"
  WRITE(*,*)
  CALL evec_print(ndim,nbas,N,enum,eval,Cij,error)
  
END SUBROUTINE H_HO_diag

!------------------------------------------------------------
END MODULE H_HO
!------------------------------------------------------------
