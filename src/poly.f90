!//////////////////////////////////////////////////////////////////
!///
!///            Module that contains subroutines for 
!///            fitting and evaluating polynomial functions
!///            to an arbitrary order
!///
!//////////////////////////////////////////////////////////////////

MODULE poly
  USE linal
  USE val

CONTAINS
!---------------------------------------------------------------------
! poly_fit
!       - subrotuine that fits a function to an increasing number of
!         polynomials, until a threshold RMS is met
!---------------------------------------------------------------------
! Variables
! Fx            : 1D real*8, function to fit
! x             : 1D real*8, list of x values
! N             : int, number of points
! tol           : real*8, RMS target tollerance
! fit           : int, type of basis functions
! ord           : int, on exit, order of the polynomial used 
! coef          : 2D real*8, list of coefficients used
! error         : bool, true on exit if there was a problem

SUBROUTINE poly_fit(Fx,x,N,tol,fit,ord,coef,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: coef
  REAL(KIND=8), DIMENSION(0:N-1), INTENT(IN) :: Fx,x
  REAL(KIND=8), INTENT(IN) :: tol
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(INOUT) :: ord
  INTEGER, INTENT(IN) :: N,fit

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A,b
  REAL(KIND=8) :: rms,rms_old,au2cm
  INTEGER :: i,j,nc,dof

  error = .FALSE.
  au2cm = 2.1947E5

  WRITE(*,*) 
  WRITE(*,*) "Generating polynomial fit of surface"
  WRITE(*,*) "Convergence :", tol

  IF (fit .EQ. 0) THEN
    WRITE(*,*) "Using functions of type a*exp(-b*r)"
    WRITE(*,*) "Actually, this is broken, exiting"
    error = .TRUE.
    RETURN
    ALLOCATE(A(0:N-1,0:2*N-1))
  ELSE IF (fit .EQ. 1) THEN
    WRITE(*,*) "Using functions of type a*r" 
    ALLOCATE(A(0:N-1,0:N-1))
  ELSE
    WRITE(*,*) "That kind of fit is not coded, yet"
    error = .TRUE.
    RETURN
  END IF

  dof = 1
  
  ALLOCATE(b(0:N-1,0:0))
  A = 0
  b = 0

  !loop until we get below tolerance
  !DO i=0,1!N-1 
  DO i=0,N/2-1
    ord = i
    IF (ord .EQ. 0) THEN
      IF (fit .EQ. 0) THEN
        nc = 2*(ord+1)
        A(0:N-1,2*ord) = (/ (1.0,j=0,N-1) /)
        A(0:N-1,2*ord+1) = (/ (1.0,j=0,N-1) /)
        CALL lsqr(A(0:N-1,0:nc-1),Fx(0:N-1),N,nc,dof,b(0:nc-1,0:dof-1),error) 
        IF (error) THEN
          WRITE(*,*) "poly:poly_fit -- lsqr exited with error" 
          error = .TRUE. 
          DEALLOCATE(A)
          DEALLOCATE(b)
          RETURN
        END IF 
        CALL poly_rms(Fx,x,N,dof,b(0:nc-1,0:dof-1),nc,ord,fit,rms,error)
        IF (error) THEN
          WRITE(*,*) "poly:poly_fit -- poly_rms exited with error" 
          error = .TRUE. 
          DEALLOCATE(A)
          DEALLOCATE(b)
          RETURN
        END IF 
      ELSE IF (fit .EQ. 1) THEN
        nc = 1
        A(0:N-1,ord) = (/ (1.0D0,j=0,N-1) /) 
        CALL lsqr(A(0:N-1,0:nc-1),Fx(0:N-1),N,nc,dof,b(0:nc-1,0:dof-1),error) 
        IF (error) THEN
          WRITE(*,*) "poly:poly_fit -- lsqr exited with error" 
          error = .TRUE. 
          DEALLOCATE(A)
          DEALLOCATE(b)
          RETURN
        END IF 
        CALL poly_rms(Fx,x,N,dof,b(0:nc-1,0:dof-1),nc,ord,fit,rms,error)
        IF (error) THEN
          WRITE(*,*) "poly:poly_fit -- poly_rms exited with error" 
          error = .TRUE. 
          DEALLOCATE(A)
          DEALLOCATE(b)
          RETURN
        END IF 
      END IF
    ELSE
      IF (fit .EQ. 0) THEN
        nc = 2*(ord+1)
        !A(0:N-1,2*ord) = (/ (1.0D0,j=0,N-1) /)
        A(0:N-1,2*ord) = (/ (EXP(-ord*x(j)),j=0,N-1) /)
        A(0:N-1,2*ord+1) = (/ (EXP(-ord*x(j)),j=0,N-1) /)
        CALL lsqr(A(0:N-1,0:nc-1),Fx(0:N-1),N,nc,dof,b(0:nc-1,0:dof-1),error) 
        IF (error) THEN
          WRITE(*,*) "poly:poly_fit -- lsqr exited with error" 
          error = .TRUE. 
          DEALLOCATE(A)
          DEALLOCATE(b)
          RETURN
        END IF 
        CALL poly_rms(Fx,x,N,dof,b(0:nc-1,0:dof-1),nc,ord,fit,rms,error)
        IF (error) THEN
          WRITE(*,*) "poly:poly_fit -- poly_rms exited with error" 
          error = .TRUE. 
          DEALLOCATE(A)
          DEALLOCATE(b)
          RETURN
        END IF 
      ELSE
        nc = ord + 1
        A(0:N-1,ord) = (/ (A(j,ord-1)*x(j) ,j=0,N-1) /) 
        CALL lsqr(A(0:N-1,0:nc-1),Fx(0:N-1),N,nc,dof,b(0:nc-1,0:dof-1),error) 
        IF (error) THEN
          WRITE(*,*) "poly:poly_fit -- lsqr exited with error" 
          error = .TRUE. 
          DEALLOCATE(A)
          DEALLOCATE(b)
          RETURN
        END IF 
        CALL poly_rms(Fx,x,N,dof,b(0:nc-1,0:dof-1),nc,ord,fit,rms,error)
        IF (error) THEN
          WRITE(*,*) "poly:poly_fit -- poly_rms exited with error" 
          error = .TRUE. 
          DEALLOCATE(A)
          DEALLOCATE(b)
          RETURN
        END IF 
      END IF
    END IF

!    WRITE(*,*) "At iteration ", ord," rms is", rms
    IF (rms .LT. tol) EXIT
    IF (rms .GT. rms_old .AND. i .GT. 1) EXIT
    rms_old = rms
    
  END DO

  IF (rms .GT. tol) THEN
    error = .TRUE.
    WRITE(*,*) 
    WRITE(*,*) "poly:poly_fit -- failed to get below RMS target"
  ELSE
    WRITE(*,*) 
    WRITE(*,'(2x,A17,I4,A11)') "RMS converged in ", ord, " iterations"
    WRITE(*,'(2x,A4,2x,F11.8,2x,A2,2x,F8.4,2x,A4)') &
                "RMS=",rms,"au",rms*au2cm,"cm-1"
    ALLOCATE(coef(0:nc-1,0:dof-1))
    coef(0:nc-1,0:dof-1) = b(0:nc-1,0:dof-1)
  END IF
  OPEN(file='fit.dat',unit=105,status='replace')
  WRITE(105,*) ord
  DO i=0,ord
    !WRITE(105,*) i, b(i,0:dof-1)
    IF (fit .EQ. 0) THEN
      WRITE(105,*) i,b(i,0),b(i+ord+1,0)
    ELSE 
      WRITE(105,*) i, b(i,0)
    END IF
  END DO
  CLOSE(unit=105)
  WRITE(*,*) "Coefficients written to fit.dat"

  CALL poly_plot(Fx,x,N,ord,coef,error)

   DEALLOCATE(A)
   DEALLOCATE(b)

END SUBROUTINE poly_fit

!---------------------------------------------------------------------
! poly_plot
!       - plots the polynomial fit of Vq
!---------------------------------------------------------------------
SUBROUTINE poly_plot(Fx,x,N,ord,coef,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: coef
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Fx,x
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,ord

  INTEGER :: i

  error = .FALSE.
  WRITE(*,*) "Saving fit runfile to fitplot"

  OPEN(unit=101,file='plot.dat',status='replace')
  DO i=0,N-1
    WRITE(101,*) x(i),Fx(i)
  END DO
  CLOSE(unit=101)

  OPEN(unit=100,file='fitplot',status='replace')
  WRITE(100,*) "set terminal png"
  WRITE(100,*) "set output 'fit.png'"
  WRITE(100,*) "set xrange [", x(0), ":", x(N-1),"]"
  WRITE(100,*) "set yrange [", MINVAL(Fx),":",MAXVAL(Fx),"]"
  WRITE(100,*) "plot 'plot.dat' u 1:2 t 'Vq',\"
  DO i=0,ord-1
    WRITE(100,*) coef(i,0),"*x**",i," + \"
  END DO
  WRITE(100,*) coef(ord,0),"*x**",ord
  CLOSE(unit=100)

END SUBROUTINE poly_plot
!---------------------------------------------------------------------
! poly_rms
!       - returns rms of poynomial fit
!---------------------------------------------------------------------
SUBROUTINE poly_rms(Fx,x,N,dof,coef,nc,ord,fit,rms,error)
  IMPLICIT NONE
  
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: coef
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Fx,x
  REAL(KIND=8), INTENT(INOUT) :: rms
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N,ord,fit,dof,nc

  REAL(KIND=8) :: val
  INTEGER :: i

  error = .FALSE.

  rms = 0.0D0
  OPEN(file='rms.dat',unit=107,status='replace')
  DO i=0,N-1
    CALL poly_eval(fit,ord,coef(0:nc-1,0:dof-1),x(i),val,error) 
    IF (error) RETURN
     WRITE(107,*) x(i), (Fx(i) - val)**2.0
    rms = rms + (Fx(i) - val)**2.0
  END DO
  rms = rms/N
  CLOSE(unit=107)

END SUBROUTINE poly_rms

!---------------------------------------------------------------------
! poly_eval
!       - given order and coefficients, returns the value of 
!         a fitted polynomial at x
!---------------------------------------------------------------------
SUBROUTINE poly_eval(fit,ord,coef,x,val,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: coef
  REAL(KIND=8), INTENT(INOUT) :: val
  REAL(KIND=8), INTENT(IN) :: x
  LOGICAL, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ord,fit

  INTEGER :: i
  
  error = .FALSE.

  val = 0
  IF (fit .EQ. 0) THEN
    DO i=0,ord
      val = val + coef(2*i,0)*EXP(-i*coef(2*i+1,0)*x)
      !val = val + coef(i,0)*EXP(-i*coef(i,1)*x**2.0) 
    END DO
  ELSE IF (fit .EQ. 1) THEN
    DO i=0,ord !remember that ord = 1 is a 1st order polynomial, a*x + b
      val = val + coef(i,0)*x**i
    END DO
  ELSE
    WRITE(*,*) "poly:poly_eval -- that fittting function not coded"
    error = .TRUE.
    RETURN
  END IF

  !check for inf/nan
  CALL checkval(val,error)
  IF (error) THEN
    WRITE(*,*) "poly:poly_eval -- bad value in evalulation" 
    RETURN
  END IF

END SUBROUTINE poly_eval

!---------------------------------------------------------------------

END MODULE poly
