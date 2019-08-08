!------------------------------------------------------------
!  vho
!       -Control program for variational 
!        harmonic oscillator calculations
!------------------------------------------------------------
! job           : int, job type
! bas           : int, basis type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! enum          : int, number of eigenvalues to compute          
! mem           : int*8, memory in MB
! error         : int, error code
! Hij           : 2D real*8, hamiltonian
! q             : 2D real*8, list of abscissa    [abscissa,dimension]
! W             : 2D real*8, list of weights     [abscissa,dimension]
! nabs          : 1D int, number of abscissa

PROGRAM vho
  USE input
  USE gauss 
  USE H_HO

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Hij,Cij,q,W
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: eval
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbas,nabs
  INTEGER(KIND=8) :: mem
  REAL(KIND=8) :: ti,tf
  INTEGER :: ndim,job,bas,error,enum
  
  CALL CPU_TIME(ti)
  error = 0
  WRITE(*,*) "xvho : called"
  CALL vho_strtmsg()
  
  CALL input_jobinfo(job,bas,ndim,nbas,enum,mem,error)
  IF (error .NE. 0) THEN
    CALL CPU_TIME(tf)
    CALL vho_endmsg(ti,tf,error)
    STOP 1
  END IF
  
  !generate abscissa if needed
  IF (job .EQ. 2 .OR. job .EQ. 3) THEN
    CALL gauss_read(job,ndim,nabs,error)
    IF (error .NE. 0) THEN
      CALL CPU_TIME(tf)
      CALL vho_endmsg(ti,tf,error)
      STOP 1
    END IF
    CALL gauss_generate(job,bas,ndim,mem,nabs,q,W,error)
  ELSE IF (job .LT. 0) THEN
    CALL gauss_generate(job,bas,ndim,mem,nbas,q,W,error)
  END IF
  IF (error .NE. 0) THEN
    CALL CPU_TIME(tf)
    CALL vho_endmsg(ti,tf,error)
    STOP 1
  END IF

  IF (job .LT. 0) THEN
    CALL CPU_TIME(tf)
    CALL vho_endmsg(ti,tf,error)
    STOP 0
  END IF

  IF (bas .EQ. 1) THEN
    CALL H_HO_build(job,ndim,nbas,nabs,mem,q,W,Hij,error)  
    IF (error .NE. 0) THEN
      CALL CPU_TIME(tf)
      CALL vho_endmsg(ti,tf,error)
      STOP 1
    END IF
    CALL H_HO_diag(ndim,nbas,enum,mem,Hij,eval,Cij,error)
  END IF

  CALL CPU_TIME(tf)
  CALL vho_endmsg(ti,tf,error)

CONTAINS

!------------------------------------------------------------
! vho_strtmsg
!       - prints starting message
!------------------------------------------------------------

SUBROUTINE vho_strtmsg()
  IMPLICIT NONE

  WRITE(*,*) 
  WRITE(*,*) "======================================================="
  WRITE(*,*) "Starting xvho"
  WRITE(*,*)  

END SUBROUTINE vho_strtmsg

!------------------------------------------------------------
! vho_endmsg
!       - prints ending vho message
!------------------------------------------------------------
! ti            : real*8, start time
! tf            : real*8, end time
! error         : int, exit code

SUBROUTINE vho_endmsg(ti,tf,error)
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: ti,tf
  INTEGER, INTENT(IN) :: error

  WRITE(*,*)
  WRITE(*,'(1x,A26,I2,A4,F6.1,1x,A1)') "xvho finished with status ", error," in ",tf-ti,"s" 
  WRITE(*,*) "======================================================="
  WRITE(*,*) 

END SUBROUTINE vho_endmsg
!------------------------------------------------------------

END PROGRAM vho
!------------------------------------------------------------
