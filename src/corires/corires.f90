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
! eps           : 1D real*8, list of harmonic state energies
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
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: phi3,mu1,zeta
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: phi2,eps
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
  CALL mu_get(nvib,voff,Be,mu1,mu2,error)
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
!    IF (TRIM(line) .EQ. "levels")
    ELSE IF (TRIM(line) .EQ. "diag") THEN 
      res = .TRUE.
    ELSE IF (TRIM(line) .EQ. "quit" .OR. &
             TRIM(line) .EQ. "quit()") THEN
      EXIT
    ELSE IF (TRIM(line) .EQ. "calc") THEN
!      IF (res) CALL calc_diag()
      IF (.NOT. res) CALL calc_states(nvib,voff,nstates,l2h,states,phi2,&
                                      phi3,Be,zeta,mu1,mu2) 
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
END PROGRAM corires
!------------------------------------------------------------
