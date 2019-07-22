!------------------------------------------------------------
!  vho
!       -Control program for variational 
!        harmonic oscillator calculations
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB
! error         : int, error code
! Hij           : 2D real*8, hamiltonian
! Herm          : 2D real*8, Hermite polynomials [abscissa,dimension]
! q             : 2D real*8, list of abscissa    [abscissa,dimension]
! W             : 2D real*8, list of weights     [abscissa,dimension]
! nabs          : 1D int, number of abscissa

PROGRAM vho
  USE input
  USE gauss 
  USE H
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Hij,Herm,q,W
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbas,nabs
  INTEGER(KIND=8) :: mem
  REAL(KIND=8) :: ti,tf
  INTEGER :: ndim,job,error
  
  CALL CPU_TIME(ti)
  error = 0
  WRITE(*,*) "xvho : called"
  CALL vho_strtmsg()
  
  CALL input_jobinfo(job,ndim,nbas,mem,error)
  CALL gauss_generate(job,ndim,nbas,mem,nabs,q,W,error)

  IF (job .NE. -1) THEN
    CALL H_build(job,ndim,nbas,nabs,mem,q,W,Hij,Herm,error)  
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

  WRITE(*,'(1x,A26,I2,A4,F6.1,1x,A1)') "xvho finished with status ", error," in ",tf-ti,"s" 
  WRITE(*,*) "======================================================="
  WRITE(*,*) 

END SUBROUTINE vho_endmsg
!------------------------------------------------------------

END PROGRAM vho
!------------------------------------------------------------
