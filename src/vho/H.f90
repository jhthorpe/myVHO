!------------------------------------------------------------
! H
!       - module containing subroutines dealing with the 
!         Hamiltonian
!------------------------------------------------------------
MODULE H

CONTAINS
!------------------------------------------------------------
! H_build
!       - builds the Hamiltonian
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of normal coords
! nbas          : 1D int, number of basis functions in each dim 
! mem           : int*8, memory in MB
! nabs          : int, number of abscissa
! q             : 1D real*8, list of abscissa
! W             : 1D real*8, list of weights
! H             : 2D real*8, hamiltonian
! error         : int, error code

SUBROUTINE H_build(job,ndim,nbas,nabs,mem,q,W,H,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: H
  REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: q,W
  INTEGER, DIMENSION(0:), INTENT(IN) :: nbas
  INTEGER(KIND=8), INTENT(IN) :: mem
  INTEGER, INTENT(IN) :: job,ndim,nabs
  INTEGER, INTENT(INOUT) :: error

  INTEGER(KIND=8) :: qmem
  REAL(KIND=8) :: ti,tf
  INTEGER :: i,j

  CALL CPU_TIME(ti)
  error = 0
  WRITE(*,*) "H_build  : called"

  !analyze the memory situation 
  qmem = mem*1000000/8 !memory in qwords, 1 real*8 per qword

  !Generate Hermite polynomials

  !Get potential integrals
  IF (job .EQ. 0) THEN  !read the integrals
    WRITE(*,*) "Sorry, this jobtype is not supported yet",0
    error = 1
    RETURN
  ELSE IF (job .EQ. 1) THEN  !generate the potential integrals
    
  END IF

  CALL CPU_TIME(tf)
  WRITE(*,'(1x,A29,I2,A4,F6.1,1x,A1)') "H_build  : finished with status ", error", in ", tf-ti, "s"

END SUBROUTINE H_build



!------------------------------------------------------------
END MODULE H
!------------------------------------------------------------
