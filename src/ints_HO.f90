!//////////////////////////////////////////////////////////////////
!///
!///		Module for calculating harmonic oscillator 1D 
!///		integrals
!///
!//////////////////////////////////////////////////////////////////

MODULE ints_HO
  USE fit
  USE nints
  USE val
  IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------
!       HO1D_integrals 
!               -calculates 1D harmonic oscillator integrals 
!---------------------------------------------------------------------
! Variables
! func          : int, type of fitting for surface
! N             : int, number of harmonic oscillator basis functions
! Hij		: 2D real*8, Hamiltonian Matrix
! Ni		: 1D real*8, Normalization constants
! Vq            : 1D real*8, 1D potential energy surface
! qmin          : real*8, minimum r
! qmax          : real*8, max r
! qeq           : real*8, equilibrium q
! Ncon          : 1D real*8, list of normalization constants
! error         : bool, true if error

SUBROUTINE HO1D_integrals(func,N,Vq,q,qmin,qmax,qeq,np,&
                          k,m,Voff,a,Hij,Ni,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Ni
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Vq,q
  REAL(KIND=8), INTENT(IN) :: qmin,qmax,qeq,k,m,Voff,a
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N, np, func

  INTEGER :: i
  error = .FALSE.
  
  WRITE(*,*) &
"---------------------------------------------------------------------"

  WRITE(*,*) 
  WRITE(*,*) "Calculating HO integrals"

  ALLOCATE(Hij(0:N-1,0:N-1))
  ALLOCATE(Ni(0:N-1))
  Hij = 0.0D0
  Ni = 0.0D0

  IF (func .EQ. -1) THEN
    WRITE(*,*) "Why am I here??" 
    error = .TRUE.
    RETURN 
  ELSE IF (func .EQ. 0) THEN
    WRITE(*,*) "Sorry, that integral type not supported"
    error = .TRUE.
    RETURN 
  ELSE IF (func .EQ. 1) THEN
    WRITE(*,*) "Sorry, that integral type not supported"
    error = .TRUE.  
    RETURN
  ELSE IF (func .EQ. 2) THEN
    WRITE(*,*) "Sorry, that integral type not supported"
    error = .TRUE.  
    RETURN
  ELSE IF (func .EQ. 3) THEN
    CALL gauss_potential(Hij,Ni,N,np,q(1),q(np-2),Vq(np-1),a,k,error) 
    IF (error) RETURN
  ELSE IF (func .EQ. 4) THEN
    CALL RHS_potential(Hij,Ni,N,Vq,q,np,k,m,Voff,a,error)
    IF (error) RETURN 
  ELSE
    error = .TRUE.
    WRITE(*,*) "HO1D_integrals -- this option not implemented"
  END IF

  CALL HO_harmonic(Hij,N,k,m,error) !this is the same no matter what
 
END SUBROUTINE HO1D_integrals

!---------------------------------------------------------------------
!	gauss_potential
!		-calculates 1D HO potential energy integrals
!		-stored upper triangular
!		-using x form of harmonic oscillator
!---------------------------------------------------------------------
! Variables
! Hij		: 2D real*8, Hamiltonian Matrix
! Ni            : 1D real*8, list of normalization consts
! nb 	        : int, number of harmonic oscillator basis functions
! np    	: int, number of points in Vq and q
! qmin		: real*8, minimum trusted position value
! qmax          : real*8, maximum trusted position 
! Vnp           : real*8, value of Vq at the last position
! a             : real*8, alpha
! k             : real*8, k constant
! error		: bool, true if error

SUBROUTINE gauss_potential(Hij,Ni,nb,np,qmin,qmax,Vnp,a,k,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Ni
  REAL(KIND=8), INTENT(IN) :: qmin,qmax,Vnp,a,k
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: nb,np

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Vq,q,Wq,Ht
  INTEGER :: u,na,i,j
  
  error = .FALSE.

  WRITE(*,*) 
  WRITE(*,*) "Starting Gaussian Quadrature for potential energy integrals"

  ALLOCATE(Ht(0:nb-1)) !this is the hermite polynomial table
  
  CALL read_gauss(Vq,q,Wq,na,error)!read in abscissa, weights, and potentials 
  !q was generated as y = a*x
  !we are going to very carefully work in the x, not y integrals
  !where I am implicitly assume that the locations of the abscissa
  !are just scaled

  !WRITE(*,*) "q in y=a*x, x"
  !DO i=0,na-1
  !  WRITE(*,*) q(i), q(i)/a   
  !END DO
  !WRITE(*,*) "qmin is", qmin
  !WRITE(*,*) "qmax is", qmax

  DO u=0,na-1
    CALL build_Htab(nb,q(u),Ht(0:nb-1))    

    !the cases below should probably have been handled in the fitting
    !subroutine...
   
    !we are below qmin
    IF (q(u)/a .LE. qmin) THEN
      DO j=0,nb-1
        Ni(j) = Ni(j) + Wq(u)*Ht(j)*Ht(j)
      END DO 

    !we are above qmax 
    ELSE IF (q(u)/a .GE. qmax) THEN
      DO j=0,nb-1
        DO i=j,nb-1
          Hij(i,j) = Hij(i,j) + Wq(u)*Ht(i)*Ht(j)*&
                   (Vnp - 0.5*k*(q(u)/a)**2.0D0) 
        END DO
        Ni(j) = Ni(j) + Wq(u)*Ht(j)*Ht(j)
      END DO

    !we are at a trusted value 
    ELSE
      DO j=0,nb-1
        DO i=j,nb-1
          Hij(i,j) = Hij(i,j) + Wq(u)*Ht(i)*Ht(j)*&
                   (Vq(u) - 0.5*k*(q(u)/a)**2.0D0) 
        END DO
        Ni(j) = Ni(j) + Wq(u)*Ht(j)*Ht(j)
      END DO
    END IF
  END DO 

  ! account for the change of variables in the integral
  Hij = Hij/a
  Ni = Ni/a
  Ni = 1.0/SQRT(Ni) 

  !check values and normalize
  DO j=0,nb-1
    CALL checkval(Ni(j),error)
    IF (error) THEN
      WRITE(*,*) "gauss_potential -- bad integral at N",j
      DEALLOCATE(Vq)
      DEALLOCATE(q)
      DEALLOCATE(Wq)
      DEALLOCATE(Ht)
      RETURN
    END IF
    
    DO i=j,nb-1
      CALL checkval(Hij(i,j),error)
      IF (error) THEN
        WRITE(*,*) "gauss_potential -- bad integral at H", i,j
        DEALLOCATE(Vq)
        DEALLOCATE(q)
        DEALLOCATE(Wq)
        DEALLOCATE(Ht)
        RETURN
      END IF
      Hij(i,j) = Hij(i,j)*Ni(i)*Ni(j) 
    END DO
  END DO

  
  !WRITE(*,*) "Hamiltonian is..."
  !DO i=0,nb-1
  !  WRITE(*,*) Hij(i,0:i)
  !END DO

END SUBROUTINE gauss_potential

!---------------------------------------------------------------------
! read_gauss
!       - reads in information about abscissa, their weights, and the 
!         potential energies at these points
!---------------------------------------------------------------------
! Variables
! Vq            : 1D real*8, list of potential energies at abscissa
! q             : 1D real*8, list of abscissa
! Wq            : 1D real*8, list of weights at abscissa
! na            : int, number of abscissa
! error         : bool, true on exit if error

SUBROUTINE read_gauss(Vq,q,Wq,na,error)
  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Vq,q,Wq
  INTEGER, INTENT(INOUT) :: na
  LOGICAL, INTENT(INOUT) :: error

  LOGICAL :: exists
  INTEGER :: i

  error = .FALSE.

  INQUIRE(file='gauss.dat',EXIST=exists)
  IF (.NOT. exists) THEN
    WRITE(*,*) "HO1D_ints:read_gauss -- you need the file 'gauss.dat'"
    error = .TRUE. 
    RETURN
  END IF

  OPEN(file='gauss.dat',unit=108,status='old')
  READ(108,*) na
  ALLOCATE(Vq(0:na-1))
  ALLOCATE(q(0:na-1))
  ALLOCATE(Wq(0:na-1))
  DO i=0,na-1
    READ(108,*) q(i), Wq(i), Vq(i)
  END DO
  CLOSE(unit=108)

END SUBROUTINE read_gauss

!---------------------------------------------------------------------
!	RHS_potential
!		-calculates 1D HO potential energy integrals
!		-stored upper triangular
!		-using x form of harmonic oscillator
!---------------------------------------------------------------------
! Variables
! N 	        : int, number of harmonic oscillator basis functions
! Hij		: 2D real*8, Hamiltonian Matrix
! Ni            : 1D real*8, list of normalization consts
! Vq            : 1D real*8, 1D potential energy surface
! q		: 1D real*8, list of q values
! np    	: int, number of points in Vq and q
! error		: bool, true if error

SUBROUTINE RHS_potential(Hij,Ni,N,Vq,q,np,k,m,Voff,a,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:),INTENT(INOUT) :: Ni
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Vq,q
  REAL(KIND=8), INTENT(IN) :: k,m,a,Voff
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N, np

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Htab
  REAL(KIND=8) :: dq,dx,xlim,x,tol
  REAL(KIND=8) :: infty
  INTEGER :: i,j,u,xpoints

  error = .FALSE.
  WRITE(*,*) "Calculating potential energy numerical integrals"
  WRITE(*,*) "Number of points :", np
 
  ALLOCATE(Htab(0:N-1))
  infty = HUGE(tol)

  !DO u=1,np-2
  DO u=0,np-2
  !DO u=0,np-1
   
    !Construct Hermitian table
    CALL build_Htab(N,a*q(u),Htab(0:N-1))
    dq = ABS(q(u+1) - q(u))

    DO j=0,N-1
      DO i=j,N-1
        Hij(i,j) = Hij(i,j) + dq*Htab(i)*Htab(j)*EXP(-a**2.0*q(u)**2.0)* (&
                  Vq(u) - 0.5*k*q(u)**2.0D0)
      END DO
    END DO

  END DO 

  Ni = 0
  xlim = 20
  xpoints = 1000000
  dx = (2*xlim)/xpoints

  WRITE(*,*) 
  WRITE(*,*) "Calculating Normalization Constants" 
  WRITE(*,*) "Number of points :", xpoints 
  WRITE(*,*) "between          :", -xlim,",", xlim 
  WRITE(*,*) "with dx          :", dx

  x = -xlim

  DO u=0,xpoints-2

    CALL build_Htab(N,a*x,Htab(0:N-1))
    DO j=0,N-1
        Ni(j) = Ni(j) + Htab(j)*Htab(j)*EXP(-a**2.0*x**2.0)*dx
    END DO

    x = x + dx
  END DO

  !put into correct form
  Ni = 1./SQRT(Ni)
  
  WRITE(*,*)
  WRITE(*,*) "Normalizing Wavefunction"
  DO j=0,N-1
    DO i=j,N-1
      IF (Hij(i,j) .GT. infty) THEN
        WRITE(*,'(2x,A1,2x,I4,A1,I4,2x,A12)') "H",i,",",j," is infinite"
        error  = .TRUE.
      ELSE IF (Hij(i,j) .NE. Hij(i,j)) THEN 
        WRITE(*,'(2x,A1,2x,I4,A1,I4,2x,A7)') "H",i,",",j," is NaN"
        error = .TRUE.
      END IF
      Hij(i,j) = Hij(i,j)*Ni(i)*Ni(j)
    END DO
  END DO

  DEALLOCATE(Htab)

END SUBROUTINE RHS_potential
!---------------------------------------------------------------------
!	build_Htab
!		-constructs list of Hermite polynomials at
!		a particular q
!---------------------------------------------------------------------
! Variables
! Htab		: 1D real*8, list of hermite polynomials eval. at q
! q		: real*8, value to evaluate Htab at
! N		: int, max Hermite polynomial to calc

SUBROUTINE build_Htab(N,q,Htab)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Htab 
  REAL(KIND=8), INTENT(IN) :: q
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8) :: H0, H1
  INTEGER :: i

  IF (N .EQ. 1) THEN
    Htab(0) = 1.0D0
    RETURN
  ELSE IF (N .EQ. 2) THEN
    Htab(0) = 1.0D0
    Htab(1) = 2.0*q
    RETURN
  ELSE 
    Htab(0) = 1.0D0
    Htab(1) = 2.0*q
    DO i=1,N-2
      Htab(i+1) = 2*q*Htab(i) - 2.*i*Htab(i-1)
    END DO
  END IF 

END SUBROUTINE build_Htab

!---------------------------------------------------------------------
!       RHS_normalize
!               -normalizes the matrix elements in Hij
!               -creates list of normalization constants
!---------------------------------------------------------------------
! Variables
! Hij           : 2D real*8, matrix to be normalized
! N             : int, number of basis functions
! a             : real*8, alpha of basis functions 
! Ncon          : 1D real*8, list of normalization constants
! error         : bool, true on exit if problem

SUBROUTINE  RHS_normalize(Hij,N,a,Ncon,error)
  IMPLICIT NONE
 
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: Ncon
  REAL(KIND=8), INTENT(IN) :: a
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8) :: pi, const
  INTEGER(KIND=8) :: pwr, fct
  INTEGER :: i,j

  error = .FALSE.
  pi = 3.1415926535897932
  const = a/SQRT(pi)

  Ncon(0) = SQRT(1*const) 
  
  pwr = 1
  fct = 1
  DO i=1,N-1
    pwr = pwr * 2
    fct = i*fct 
    Ncon(i) = SQRT(1.0/(pwr*fct)*const)
  END DO

  DO j=0,N-1
    DO i=j,N-1
      Hij(i,j) = Hij(i,j)*Ncon(i)*Ncon(j)
    END DO
  END DO

  WRITE(*,*) "N00", Ncon(0)*Ncon(0)
  WRITE(*,*) "N10", Ncon(1)*Ncon(0)
  WRITE(*,*) "N11", Ncon(1)*Ncon(1)

END SUBROUTINE RHS_normalize

!---------------------------------------------------------------------
!       RHS_kinetic
!               -adds in (normalized) kinetic energy to Hij
!               - this abuses the viral theorum, and the fact
!               - that our basis functions are orthogonal HO
!---------------------------------------------------------------------
! Variables
! Hij           : 2D real*8, matrix of integrals
! N             : int, number of basis functions
! k             : real*8, k of basis functions
! m             : real*8, m of basis functions
! a             : real*8, alpha of basis functions
! Ncon          : 1D real*8, list of normalization constants
! error         : bool, true on exit if problem

SUBROUTINE RHS_kinetic(Hij,N,k,m,a,Ncon,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Ncon
  REAL(KIND=8), INTENT(IN) :: k,m,a
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N
  REAL(KIND=8) :: w ,pi,C
  INTEGER :: i

  WRITE(*,*) "Calculating Kinetic Integrals"
  error = .FALSE.
  pi = 3.1415926535897932
  w = SQRT(k/m)
  C = SQRT(1/(2*w*m))

  DO i=0,N-1
    Hij(i,i) = Hij(i,i) + w*(1.0*i+0.5)/2.0 !using viral theorum
    IF (i+2 .LE. N-1) THEN
      Hij(i+2,i) = Hij(i+2,i) + (C**2.0*a**4.0*SQRT(n+2.)*SQRT(n+1.)&
                   - Ncon(i+2)*2*a**3.0*(i+2.)*C/Ncon(i+1))/(2.*m)
    END IF

    !+2 off diagonal terms... 
    ! I have been lazy and do not have analytical derivatives...
    !H(i+2,i) = H(i+2,j) + nint_fdif_2drv_HO(i,i+2,a,m,qmin,qmax,qeq)
  END DO

END SUBROUTINE RHS_kinetic

!---------------------------------------------------------------------
!       HO_harmonic
!               -adds in harmonic energy terms to diagaonals 
!---------------------------------------------------------------------
! Variables
! Hij           : 2D real*8, matrix of integrals
! N             : int, number of basis functions
! k             : real*8, k of basis functions
! m             : real*8, m of basis functions
! a             : real*8, alpha of basis functions
! Ncon          : 1D real*8, list of normalization constants
! error         : bool, true on exit if problem

SUBROUTINE HO_harmonic(Hij,N,k,m,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hij
  REAL(KIND=8), INTENT(IN) :: k,m
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N
  REAL(KIND=8) :: w
  INTEGER :: i

  error = .FALSE.
  w = SQRT(k/m)

  DO i=0,N-1
    Hij(i,i) = Hij(i,i) + w*(1.0*i+0.5) !using Viral theorum 
!    Hij(i,i) = Hij(i,i) + -1.406912999D-5
  END DO

END SUBROUTINE HO_harmonic

!---------------------------------------------------------------------

END MODULE ints_HO
