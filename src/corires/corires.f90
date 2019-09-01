!------------------------------------------------------------
!  corires
!       - control program for treating coriolis resonances
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! nrot          : int, number of rotational modes
! voff          : int, vibrational numbering offset
! Be            : 1D real*8, rotational constants (x,y,z)
! phi2          : 1D real*8, quadratic force constants (cm-1)
! phi3          : 3D real*8, cubic force constants (cm-1)
! h2l           : 1D int, list of harmonic assignments -> labels
! l2h           : 1D int, list of labels -> harmonic assignments
! mu1           : 3D real*8, μ_i^{a,b} terms (vib,rot,rot) 
! mu2           : 4D real*8, μ_{i,j}^{a,b}, terms (vib,vib,rot,rot)
! zeta          : 3D real*8, coriolis zeta matrix (vib,vib,rot)
! error         : int, exit code

PROGRAM corires
  USE quad
  USE cubi
  USE rota
  USE mu
  USE cori
  USE vstate
  USE calc
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: mu2
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: phi3,mu1,zeta,didq
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: phi2
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: states
  INTEGER, DIMENSION(:), ALLOCATABLE :: h2l,l2h
  REAL(KIND=8), DIMENSION(0:2) :: Be
  CHARACTER(LEN=1024) :: line
  INTEGER :: nvib,nrot,voff,error,nstates
  LOGICAL :: cont,res
  error = 0

  !Print Starting Info 

  !Read data
  !phi2,h2l 
  CALL quad_get(nvib,voff,phi2,h2l,l2h,error)
  IF (error .NE. 0) CALL corires_exit(error)
  CALL rota_get(nrot,Be,error)
  IF (error .NE. 0) CALL corires_exit(error)
  CALL cubi_get(nvib,voff,phi3,error)
  IF (error .NE. 0) CALL corires_exit(error)
  CALL mu_get(nvib,voff,Be,mu1,mu2,didq,error)
  IF (error .NE. 0) CALL corires_exit(error)
  CALL cori_get(nvib,voff,zeta,error)
  IF (error .NE. 0) CALL corires_exit(error)
  WRITE(*,*)

  !Read user input
  cont = .TRUE. 
  nstates = 0
  DO WHILE (cont) 
    IF (error .NE. 0) EXIT
    WRITE(*,'(A2)',ADVANCE='no') "> " 
    READ(*,*) line 
    !Read input option
    IF (TRIM(line) .EQ. "states") THEN
      CALL vstate_set(nvib,h2l,nstates,states)
    ELSE IF (TRIM(line) .EQ. "level") THEN
      CALL vstates_level(nvib,h2l,nstates,states)
    ELSE IF (TRIM(line) .EQ. "diag") THEN 
      res = .TRUE.
    ELSE IF (TRIM(line) .EQ. "quit" .OR. &
             TRIM(line) .EQ. "quit()") THEN
      EXIT
    ELSE IF (TRIM(line) .EQ. "calc") THEN
      IF (res) CALL calc_diag(nvib,voff,nstates,l2h,states,phi2,&
                              phi3,Be,zeta,mu1,mu2,didq,error)
      IF (.NOT. res) CALL calc_states(nvib,voff,nstates,l2h,states,phi2,&
                                      phi3,Be,zeta,mu1,mu2,didq) 
    ELSE IF (TRIM(line) .EQ. "print") THEN
      CALL corires_print(nvib,Be,phi2,phi3,zeta,didq,mu1,mu2,h2l)
    ELSE
      WRITE(*,*) TRIM(line), " was not recognized"
    END IF
  END DO

  !Print Ending Info
  CALL corires_exit(error)

  CONTAINS

!------------------------------------------------------------
!  corires_exit
!       - exits corires
!------------------------------------------------------------
! error         : int, exit code

SUBROUTINE corires_exit(error)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: error
  WRITE(*,'(1x,A29,2x,I2)') "corires completed with status",error
  STOP
END SUBROUTINE corires_exit

!------------------------------------------------------------
! corires_print
!       - prints everything that was read in
!------------------------------------------------------------
! nvib          : int, number of vibrational modes
! Be            : 1D real*8, rotational constants
! phi2          : 1D real*8, quadratic force constants
! phi3          : 3D real*8, cubic force constants
! zeta          : 3D real*8, Coriolis zetas
! didq          : 3D real*8, didq terms
! mu1           : 3D real*8, order 1 μ term
! mu2           : 4D real*8, order 2 μ term
! h2l           : 1D int, harmonic -> label index

SUBROUTINE corires_print(nvib,Be,phi2,phi3,zeta,didq,mu1,mu2,h2l)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:,0:,0:), INTENT(IN) :: mu2
  REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: phi3,mu1,zeta,didq
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: phi2,Be
  INTEGER, DIMENSION(0:), INTENT(IN) :: h2l
  INTEGER, INTENT(IN) :: nvib
  INTEGER :: i,j,k,a,b

  WRITE(*,*)
  WRITE(*,*) "Rotational Constants (cm-1)"
  DO a=0,2
    WRITE(*,'(1x,I2,2x,F24.13)') a+1,Be(a)
  END DO
  WRITE(*,*)
  WRITE(*,*) "Harmonic Frequencies (cm-1)"
  DO i=0,nvib-1
    WRITE(*,'(1x,I2,2x,F24.13)') h2l(i),phi2(i)
  END DO
  WRITE(*,*)
  WRITE(*,*) "Cubic Force Constants (cm-1)"
  DO i=0,nvib-1
    DO j=i,nvib-1
      DO k=j,nvib-1
        WRITE(*,'(1x,3(I2,1x),2x,F24.13)'),h2l(i),h2l(j),h2l(k),phi3(i,j,k)
      END DO
    END DO
  END DO
  WRITE(*,*)
  WRITE(*,*) "Coriolis Zetas ζ_α^{i,j} (cm-1)"
  DO a=0,2
    DO i=0,nvib-1
      DO j=0,i-1
        WRITE(*,'(1x,3(I2,1x),2x,F24.13)') a+1,h2l(i),h2l(j),zeta(i,j,a)  
        WRITE(*,'(1x,3(I2,1x),2x,F24.13)') a+1,h2l(j),h2l(i),zeta(j,i,a)  
      END DO
    END DO
  END DO
  WRITE(*,*)
  WRITE(*,*) "ℏ^2/(2hc) dI_{α,β}/dq_i (cm-1)"
  DO a=0,2
    DO b=a,2
      DO i=0,nvib-1
        WRITE(*,'(1x,3(I2,1x),2x,F24.13)') a+1,b+1,h2l(i),didq(i,a,b)
      END DO
    END DO
  END DO
  WRITE(*,*) 
  WRITE(*,*) "μ_{α,β}^i (cm-1)"
  DO a=0,2
    DO b=a,2
      DO i=0,nvib-1
        WRITE(*,'(1x,3(I2,1x),2x,F24.13)') a+1,b+1,h2l(i),mu1(i,a,b)
      END DO
    END DO
  END DO
  WRITE(*,*) 
  WRITE(*,*) "μ_{α,β}^{i,j} (cm-1)"
  DO a=0,2
    DO b=a,2
      DO i=0,nvib-1
        DO j=i,nvib-1
          WRITE(*,'(1x,4(I2,1x),2x,F24.13)') a+1,b+1,h2l(i),h2l(j),mu2(i,j,a,b)
        END DO
      END DO
    END DO
  END DO
  WRITE(*,*) 
  WRITE(*,*)

END SUBROUTINE corires_print

!------------------------------------------------------------
END PROGRAM corires

!------------------------------------------------------------
